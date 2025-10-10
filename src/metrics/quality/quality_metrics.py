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
    
    sdata['transcripts']['assigned'] = sdata['transcripts']['cell_id'] != 0
    
    # Proportion of assigned reads
    prop_of_assigned_reads = float(((sdata['transcripts']['assigned']).sum() / len(sdata['transcripts'])).compute())
    
    # Proportion of assigned reads per gene
    df = pd.crosstab(sdata['transcripts']['feature_name'], sdata['transcripts']['assigned'])
    prop_of_assigned_reads_per_gene = df[True] / (df[False] + df[True])
    
    return prop_of_assigned_reads, prop_of_assigned_reads_per_gene
    
    
    
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
    
    
    

