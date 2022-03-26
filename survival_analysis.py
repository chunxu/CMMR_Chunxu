#Survival analysis and estimation of hazard ratios in different variants
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import warnings
%matplotlib inline
warnings.filterwarnings('ignore')

#Loading the dataset
df = pd.read_csv('survival_data.csv')
#see what the data look like
df.head() # all three var and time duration are continuous numbers

df.shape
# statistical information on the numbers in the dataset
df.describe()
# datatype information
df.info()
# check for null values
df.isnull().sum() #perfect! No null values found.

# create box plots to check outliers:
fig, ax = plt.subplots(ncols=5, nrows=1, figsize=(10, 5))
index = 0
ax = ax.flatten()

for col, value in df.items():
    sns.boxplot(y=col, data=df, ax=ax[index])
    index += 1
plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=5.0) 

# create dist plot to see distribution: 
fig, ax = plt.subplots(ncols=5, nrows=1, figsize=(10, 5))
index = 0
ax = ax.flatten()

for col, value in df.items():
    sns.distplot(value, ax=ax[index])
    index += 1
plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=5.0)


#Corelation 
corr = df.corr()
plt.figure(figsize=(20,10))
sns.heatmap(corr, annot=True, cmap='coolwarm')
# all three variables are independent

#generate survival curve
from lifelines import KaplanMeierFitter
kmf = KaplanMeierFitter()
kmf.fit(durations = df['T'], event_observed = df['E'])
kmf.plot_survival_function()

# Cox Proportional Hazard Model to estimate risk factors
from lifelines import CoxPHFitter
cph = CoxPHFitter()
cph.fit(df, duration_col = 'T', event_col = 'E')
cph.print_summary()
cph.plot()

