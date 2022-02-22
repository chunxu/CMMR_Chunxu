#this script is written by Chunxu to calculate whether there is relationship between DNA concentrations and Ct values in qPCR

#import data sheet
import pandas as pd
con = pd.read_csv('C:/gao/BCM2021/Conner/Corralation between DNA conc and Ct.csv')

#scatter plot data to visulize data trend
import seaborn as sns
sns.scatterplot(x="Nucleic_Acid(ng/uL)", y="Mean_Ct", data=con)
ax = sns.scatterplot(x="Nucleic_Acid(ng/uL)", y="Mean_Ct", data=con)
ax.set_title("Nucleic Acid Conc. vs. qPCR Ct values")
ax.set_xlabel("Nucleic Acid Conc")

#set best fit line
sns.lmplot(x="Nucleic_Acid(ng/uL)", y="Mean_Ct", data=con);

#test normal distribution to data yes or no
from scipy.stats import shapiro
data1 = con['Nucleic_Acid(ng/uL)']
stat, p = shapiro(data1)
print('stat=%.3f, p=%.3f' % (stat, p))
if p > 0.05:
	print('Probably Gaussian')
else:
	print('Probably not Gaussian') #P =0 indicating not normal data distribution

data2 = con['Mean_Ct']
stat, p = shapiro(data1)
print('stat=%.3f, p=%.3f' % (stat, p))
if p > 0.05:
	print('Probably Gaussian')
else:
	print('Probably not Gaussian') ##P =0 indicating not normal data distribution

#calculate corelation coefficiency
from scipy import stats
#stats.pearsonr(con['Nucleic_Acid(ng/uL)'], con['Mean_Ct']) #for normal data
stats.spearmanr(con['Nucleic_Acid(ng/uL)'], con['Mean_Ct']) #for data not fit normal distribution

#other plot if more groups compared.
#cormat = con.corr()
#round(cormat,2)
#sns.heatmap(cormat)