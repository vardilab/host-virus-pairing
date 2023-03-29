#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# In[2]:


def load_data(path,file):
    list_columns = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore".split(" ")
    data = pd.read_table(path+file, sep = "\t", names = list_columns)
    #data['cell'] = file.replace(".metapr2.tsv","")
    data['cell'] =  data['qseqid'].str.split("i\d+_", expand = True)[1]
    # in case of pr2 open:
    # data[["sseqid","rrna_type","rrna_loc","clone","Domain","Supergroup","Division","Class","Order","Family","Genus","Species"]] = data['sseqid'].str.split("|",expand = True)
    # in case of metapr2 open:

    data[["sseqid","Supergroup","Division","Class","Order","Family","Genus","Species"]] = data['sseqid'].str.split("|",expand = True)
    return data


# In[6]:


#.drop_duplicates(subset = ['qseqid','family'])
def subset_results(sample_data):
   # sample_data = l_data[num]
    print(sample_data['cell'][0])
    name = sample_data['cell'][0]
 #   print("removing 16s rRNA...")
  #  sample_data = sample_data[(sample_data['rrna_type'] == "18S_rRNA")]
    print("filtering by length and identity...")
    sample_data = sample_data[(sample_data['pident'] >= 99)]
    sample_data = sample_data[(sample_data['length'] >= 100)]
    print("sorting data...")
    sample_data = sample_data.sort_values(by = ['qseqid','bitscore'],ascending = False)

  #  sample_data = sample_data[(sample_data['length'] >= 100)]
    print("Choosing only best hit...")
    sample_data['rank'] = sample_data.groupby(['qseqid'])['bitscore'].rank("min",ascending = False)
    final_run = sample_data[sample_data['rank'] == 1]
    #.drop_duplicates(subset = ["qseqid",'Family']).drop_duplicates('qseqid',keep = False)
    return final_run


# In[7]:


path = '/home/labs/vardi/Collaboration/sc_computational_training/last/230206_marker_genes/230228_assemble_all_cells/sortmerna/'
file = 'all_cells.transcripts.edit.metaPR2.tsv'
df = load_data(path,file)
df_dedup = subset_results(df)


# In[5]:


# annotate cells based on their best contig
idx = df_dedup.groupby(['cell'])['bitscore'].transform(max) == df_dedup['bitscore']
dict_cells = dict(zip(df_dedup[idx].cell, df_dedup[idx].Family))
df_top = df_dedup[idx]
df_top['Annotation'] = df_top['Family'] 
# open if pr2
#df_dedup.to_csv(path+"transcripts.summary.99.pr2.tsv", sep = "\t")
# open if metapr2
df_dedup.to_csv(path+"transcripts.summary.99.metapr2.tsv", sep = "\t")


# In[7]:


def replace(data,level, name, new_name):
    group = set(data[data[level] == name]['Family'])
    value_to_replace = {x : new_name for x in group}
    return value_to_replace


# In[8]:


# list_replacement = []

# list_replacement.append(replace(df_top,'Class', 'Bacillariophyta','Diatoms'))
# list_replacement.append(replace(df_top,'Class', 'Prymnesiophyceae', 'Prymnesiophyceae'))
# list_replacement.append(replace(df_top,'Family', 'Noelaerhabdaceae', 'Prymnesiophyceae;Noelaerhabdaceae'))

# list_replacement.append(replace(df_top,'Class', 'Chrysophyceae','Chrysophyceae'))
# #list_replacement.append(replace(df_small,'Family', 'Chrysophyceae_Clade-C','Chrysophyceae;Clade-C'))
# list_replacement.append(replace(df_top,'Class', "MAST-3","Opalozoa;MAST-3"))
# list_replacement.append(replace(df_top,'Order', 'Strombidiida','Strombidiidae'))
# list_replacement.append(replace(df_top,'Class', 'Dinophyceae','Dinophyceae'))
# list_replacement.append(replace(df_top,'Division', 'Cercozoa','Cercozoa'))
# #list_replacement.append(replace(df_top,'Order', 'Cryomonadida','Cercozoa;Cryomonadida'))
# #list_replacement.append(replace(df_top,'Family', 'Protaspa-lineage','Cercozoa;Cryomonadida;Protaspa'))
# list_replacement.append(replace(df_top,'Order', 'Katablepharidales','Katablepharidales'))
# super_dict = {key:val for d in list_replacement for key,val in d.items()}
# #super_dict.pop("Noelaerhabdaceae", None)
# large_list = set(super_dict.keys())


# # In[9]:


# # find the other families that were not listed before
# top_families = set(df_top.Family.value_counts().head(15).index)
# top_families.remove('Urotrichidae')
# rest_top = {x : x for x in [x for x in top_families if x not in large_list]}
# named_species = set.union(*[large_list,rest_top])
# rest_species = {x : "Other eykaryotes" for x in [x for x in list(df_top.Family) if x not in named_species]}


# # In[11]:


# dicts = [super_dict,
# rest_top,
# rest_species]
# super_dict = {key:val for d in dicts for key,val in d.items()}
# df_top['Annotation'] = df_top['Family'].map(super_dict)


# # In[12]:


# df_top['Annotation'] = df_top['Family'].map(super_dict)


# # In[13]:


# dict_cells = dict(zip(df_top.cell, df_top.Annotation))


# In[15]:


#np.save(path+'dict_cells_new.npy', dict_cells) 

