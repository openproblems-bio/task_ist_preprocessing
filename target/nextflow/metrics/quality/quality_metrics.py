import numpy as np
import pandas as pd
import spatialdata as sd


def proportion_of_assigned_reads(
    sdata: sd.SpatialData,
) -> [float, pd.Series]:
    """ Calculate the proportion of assigned reads

    Parameters
    ----------
    sdata : sd.SpatialData
        SpatialData object with sdata['transcripts'] including the column 'cell_id'
        
    Returns
    -------
    float
        Proportion of assigned reads
    pd.Series
        Proportion of assigned reads per gene
        
    """
    
    # Proportion of assigned reads
    prop_of_assigned_reads = float(sdata["counts"].layers["counts"].sum() / len(sdata['transcripts']))
    
    # Proportion of assigned reads per gene
    if prop_of_assigned_reads == 1.0:
        prop_of_assigned_reads_per_gene = pd.Series(
            index=sdata['transcripts']['feature_name'].unique().compute().values,
            data=1.0
        )
    elif prop_of_assigned_reads == 0.0:
        prop_of_assigned_reads_per_gene = pd.Series(
            index=sdata['transcripts']['feature_name'].unique().compute().values,
            data=0.0
        )
    else:
        genes, counts = np.unique(sdata['transcripts']['feature_name'], return_counts=True)
        df = pd.DataFrame(index=genes, data = {"fraction":0, "count":counts, "count_assigned":0})
        df.loc[sdata["counts"].var_names, "count_assigned"] = np.array(sdata["counts"].layers["counts"].sum(axis=0))[0,:]
        df["fraction"] = df["count_assigned"] / df["count"]
        prop_of_assigned_reads_per_gene = df["fraction"]
    
    return prop_of_assigned_reads, prop_of_assigned_reads_per_gene

# Previous version only based on transcripts table. 
#
#def proportion_of_assigned_reads(
#    sdata: sd.SpatialData,
#) -> [float, pd.Series]:
#    """ Calculate the proportion of assigned reads
#
#    Parameters
#    ----------
#    sdata : sd.SpatialData
#        SpatialData object with sdata['transcripts'] including the column 'cell_id'
#        
#    Returns
#    -------
#    float
#        Proportion of assigned reads
#    pd.Series
#        Proportion of assigned reads per gene
#        
#    """
#    
#    sdata['transcripts']['assigned'] = sdata['transcripts']['cell_id'] != 0
#    
#    # Proportion of assigned reads
#    prop_of_assigned_reads = float(((sdata['transcripts']['assigned']).sum() / len(sdata['transcripts'])).compute())
#    
#    # Proportion of assigned reads per gene
#    if prop_of_assigned_reads == 1.0:
#        prop_of_assigned_reads_per_gene = pd.Series(
#            index=sdata['transcripts']['feature_name'].unique().compute().values,
#            data=1.0
#        )
#    elif prop_of_assigned_reads == 0.0:
#        prop_of_assigned_reads_per_gene = pd.Series(
#            index=sdata['transcripts']['feature_name'].unique().compute().values,
#            data=0.0
#        )
#    else:
#        df = pd.crosstab(sdata['transcripts']['feature_name'], sdata['transcripts']['assigned'])
#        prop_of_assigned_reads_per_gene = df[True] / (df[False] + df[True])
#    
#    return prop_of_assigned_reads, prop_of_assigned_reads_per_gene
    
    
    
def proportion_of_annotated_cells(
    sdata: sd.SpatialData,
) -> float:
    """ Calculate the proportion of cells that are annotated with cell types

    Parameters
    ----------
    sdata : sd.SpatialData
        SpatialData object with sdata['counts'] including the obs column 'cell_type'
        
    Returns
    -------
    float
        Proportion of cells that are annotated with cell types
    
    """
    
    cts = sdata['counts'].obs["cell_type"]
    
    n_annotated_cells = (cts.notna() & (cts != "None_sp")).sum()
    
    return float(n_annotated_cells / len(cts))
    
    
    

