#!/usr/bin/env python
# coding: utf-8

# In[14]:


import pciSeq
import scipy.io
from scipy.sparse import coo_matrix
import pandas as pd
import scanpy as sc

#INPUT: molecules.csv, singlecell.h5ad labels.mat
#OUTPUT: assignments.csv
#From config:
segmentation_method = 'imagej'
molecules = "C:/Users/Habib/Projects/HMGU/tx_project/heart/raw_data/spots_PCW4.5_1.csv"
sc_data = "C:/Users/Habib/Projects/HMGU/tx_project/heart/raw_data/heart_sc.h5ad" 
opts = {'exclude_genes': ['TCIM']}

#Read and format molecules, single cell data, and labels
spots = pd.read_csv(molecules)
spots.columns = ['Gene', 'x', 'y']

adata = sc.read_h5ad(sc_data)
scdata = adata.X
scdata  = pd.DataFrame(scdata.transpose())
scdata.columns = adata.obs['celltype']
scdata.index = adata.var_names

label = scipy.io.loadmat('data/label_{}'.format(segmentation_method))['label']
coo = coo_matrix(label)

#Run through pciSeq
pciSeq.attach_to_log()
cellData, geneData = pciSeq.fit(spots, coo, scdata, opts)

#Save in correct format
assignments = geneData[ ["Gene", "x", "y", "neighbour"] ]
assignments.columns = ["gene", "x", "y", "cell"]
assignments.to_csv("data/assignments_{}_pciseq.csv".format(segmentation_method), index = False)

