#!/usr/bin/env python
# coding: utf-8

# In[1]:


#The purpuse of this script is to help data cleaning for MILES V1 sequencing results
import numpy as np
import pandas as pd


# In[2]:


#sample list file
df1 = pd.read_csv('C:/gao/BCM2021/MILES/2. MCF Samples-MCF Submission View -all April 23 2022.csv', index_col=0)


# In[3]:


df1


# In[4]:


df1 = df1.drop(columns=['Visit'])


# In[5]:


df1 = df1.rename(columns={"4. WGS_HGSC_Report": "sample_ID"})


# In[6]:


df1


# In[7]:


df2 =pd.read_csv('C:/gao/BCM2021/MILES/Pool912/p912-Goodarzi.metaphlan2.Bacteria.EstCount.rare.csv')


# In[8]:


df2 = df2.rename(columns={"#OTU_ID":"sample_ID"})


# In[9]:


df2


# In[10]:


df3 = pd.merge(df1,df2,on='sample_ID')


# In[11]:


df3 = df3.drop(columns='MCF_SubjectID')


# In[12]:


df3 = df3.T


# In[13]:


df3


# In[14]:


### write results to excel file
df3.to_csv('C:/gao/BCM2021/MILES/Pool912/p912-Goodarzi.metaphlan2.Bacteria.EstCount.rare.clean.csv')

