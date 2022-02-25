#Genetic_test.py
#The purpuse of this script is to calculate the probability of fetal genotype (AA, AB, or BB) based on allele ratio
import numpy as np
import pandas as pd

df = pd.read_excel('C:/Users/gaoch/Downloads/B24_L02_Abias_sm.xlsx',sheet_name='Sheet1')
df.shape # check the data size
df.head() # it seems that re-arrange of the data table will be needed

# get the fetal_fraction data for calculation 
fetal_fraction = df.loc[0].index.tolist()
fetal_fraction.pop(0) #remove the 1st item which is not the value
len(fetal_fraction) #N=25 which is correct
fetal_fraction

# check whether there are abnormal values in the list of fetal fractions 
min(fetal_fraction) #0.0263164348711594
max(fetal_fraction) #0.20479337765572
# it seems that all fetal_fraction data are reasonable

df2 = pd.read_excel('C:/Users/gaoch/Downloads/B24_L02_Abias_sm.xlsx',sheet_name='Sheet1', header=1, index_col='SNP') 
#get the data table that fit commonly used pandas dataframe

df2

# get sample names
df2.columns
sample = df2.columns.to_list()
len(sample) #N=25 which is correct
sample

'''Rationale for the calculation:
    maternal contribution to A allele = 0.5*(1-corresponding fetal_fraction) 
    fetal contribution to A allele = corresponding cell value (if not NA) - maternal contribution to A allele
Then normalize the fetal contribution to A allele by fetal_fraction:
    normalized fetal contribution to A allele = fetal contribution to A allele / corresponding fetal_fraction'''

def cal_normalized_fetal_contribution(corresponding_fetal_fraction, cell_value):
    if cell_value =='NaN':
        return 'NaN'
    else: 
        maternal_contribution = 0.5*(1-corresponding_fetal_fraction)
        fetal_contribution = cell_value - maternal_contribution
        normalized_fetal_contribution = fetal_contribution / corresponding_fetal_fraction
        return normalized_fetal_contribution

#test one cell value using the fucntion
cal_normalized_fetal_contribution(fetal_fraction[0],df2[sample[0]][2]) # 0.17647817498919763 correct

#generate a new dataframe for the result of above calculation
df3=df2.copy()[[]]

#calculate all values in data frame and write the normalized_fetal_contribution to the new dataframe
for i in range(25):
    df3[sample[i]] = df2[sample[i]].apply(lambda x: cal_normalized_fetal_contribution(fetal_fraction[i],x))

df3 # show the data as normalized_fetal_contribution

'''The data above show the contributions of fetal genotypes to A allele. The value should equal to PAA + 0.5*PAB where PAA, PAB, PBB means the probability of fetal genotype (AA, AB, or BB) At the same time: PAA + PAB + PBB =1 Based on these conditions, I am not able to calculate the probability of each genotype as there are three variants and two equations. Considering AA genotype's contribution is 1, AB genotype's contribution is 0.5 while BB genotype's contribution is 0, I will give a probability of ~0 to the genotype whose contribution is farthest from the generated dataframe. And then the other two genotypes could be calculated as below.

use three functions as it would be easier to transfer the data to normal dataframe later instead of returning three items the same time
'''

def genotype_estimate_AA(normalized_fetal_contribution):
    if normalized_fetal_contribution > 1: # out of range in theory: this means even if there is 100% chance of AA genotype, the A allele still cannot be such high. So it might be sequencing or processing errors to get such a high allele A ratio. 
        BB = AB = 0
        AA = 1
    elif normalized_fetal_contribution < 0: # out of range in theory: this means even if there is no fetal contribution to A allele, mother's contribution only still is higher than the measured value.
        AA = AB =0
        BB =1
    elif normalized_fetal_contribution > 0.5:
        BB = 0
        AB = (1-normalized_fetal_contribution)*2
        AA = 2*normalized_fetal_contribution -1
    elif normalized_fetal_contribution <=0.5:
        AA = 0
        AB = 2*normalized_fetal_contribution
        BB = 1-AB
    else:
        AB = BB = AA = 'NaN'
    return AA
        
def genotype_estimate_AB(normalized_fetal_contribution):
    if normalized_fetal_contribution > 1: # out of range in theory
        BB = AB = 0
        AA = 1
    elif normalized_fetal_contribution < 0: # out of range in theory
        AA = AB =0
        BB =1
    elif normalized_fetal_contribution > 0.5:
        BB = 0
        AB = (1-normalized_fetal_contribution)*2
        AA = 2*normalized_fetal_contribution -1
    elif normalized_fetal_contribution <=0.5:
        AA = 0
        AB = 2*normalized_fetal_contribution
        BB = 1-AB
    else:
        AB = BB = AA = 'NaN'
    return AB

def genotype_estimate_BB(normalized_fetal_contribution):
    if normalized_fetal_contribution > 1: # out of range in theory
        BB = AB = 0
        AA = 1
    elif normalized_fetal_contribution < 0: # out of range in theory
        AA = AB =0
        BB =1
    elif normalized_fetal_contribution > 0.5:
        BB = 0
        AB = (1-normalized_fetal_contribution)*2
        AA = 2*normalized_fetal_contribution -1
    elif normalized_fetal_contribution <=0.5:
        AA = 0
        AB = 2*normalized_fetal_contribution
        BB = 1-AB
    else:
        AB = BB = AA = 'NaN'
    return BB
        
#test one normalized_fetal_contribution using the fucntion
genotype_estimate_AA(0.176478)

#generate a new dataframe for the result of above calculation
df4=df3.copy()[[]]
#df4

#calculate all values in data frame and write the fetal genotype data to the new dataframe
for i in range(25):
    df4[str(sample[i])+'_AA'] = df3[sample[i]].apply(lambda x: genotype_estimate_AA(x))
    df4[str(sample[i])+'_AB'] = df3[sample[i]].apply(lambda x: genotype_estimate_AB(x))
    df4[str(sample[i])+'_BB'] = df3[sample[i]].apply(lambda x: genotype_estimate_BB(x))


df4

#write results to excel file
df4.to_excel('C:/Users/gaoch/Downloads/B24_L02_Abias_sm_results-Chunxu.xlsx',sheet_name='Sheet1')

