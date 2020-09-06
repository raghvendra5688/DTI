# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.2'
#       jupytext_version: 1.2.3
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %%
import pubchempy as pcp
import pandas
from pubchempy import Compound, get_compounds
import logging
logging.basicConfig(level=logging.DEBUG)

#Fetch a compound with cid
#c = pcp.Compound.from_cid(5090)
c = Compound.from_cid(1423)
cs1 = get_compounds('Aspirin', 'name')
cs2 = get_compounds('C1=CC2=C(C3=C(C=CC=N3)C=C2)N=C1','smiles')

# %%
#To get 3d information about compounds
cs1 = pcp.get_compounds('Aspirin', 'name', record_type='3d')
cs1

# %%
c.to_dict()

# %%
#Fetch a list of compounds based on cids and print their metadata (schema)
cs_list = []
blocksize = 10
no_blocks = int(250/blocksize)
for i in range(1,no_blocks):
    if (i != no_blocks-1):
        blockrange = list(range((i-1)*(blocksize)+1,i*blocksize))
        cs_list.extend(pcp.get_compounds(blockrange,as_dataframe=False))
    if (i == no_blocks-1):
        blockrange = list(range((i-1)*(blocksize)+1,250))
        cs_list.extend(pcp.get_compounds(blockrange,as_dataframe=False))

# %%
df1 = pcp.compounds_to_frame(cs_list)
df2 = pcp.compounds_to_frame(cs_list, properties=['isomeric_smiles', 'xlogp', 'rotatable_bond_count'])
print(df1.columns)

# %%
fp=open("../Results/metadata_compounds.csv","w")
schema_list = df1.columns.tolist()
for word in schema_list:
    fp.write(word+"\n")
fp.close()

# %%
df1.to_csv("../Results/subset_compounds.csv",index=False)

# %%
