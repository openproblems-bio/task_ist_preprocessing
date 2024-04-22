import argparse
from pathlib import Path
import yaml
import pandas as pd
import anndata as ad
import txsim as tx


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--counts", type=str, help="Path to spatial data h5ad")
    parser.add_argument("-d", "--scd", type=str, help="Path to scRNAseq data reference h5ad")
    parser.add_argument("-o", "--output", type=str, help="Path to output csv with cell type annotations")
    parser.add_argument("-m", "--ct_method", type=str, help="Cell type annotation method")
    parser.add_argument("-p", "--hyper_params", type=str, help="Dictionary of method specific hyper parameters")
    parser.add_argument("-g", "--group_params", type=str, help="Dictionary of method group hyper parameters")
    return parser.parse_args()


def ct_annotation():
    args = parse_args()
    adata_sp = ad.read_h5ad(args.counts)
    adata_sc = ad.read_h5ad(args.scd)
    ct_method = args.ct_method

    # Get hyperparameters
    hparams_defaults_csv = Path(__file__).parent.parent / "configs" / "defaults.yaml"
    with open(hparams_defaults_csv, 'r') as file:
        defaults = yaml.safe_load(file)
        hparams_defaults = defaults[ct_method]
        gparams_defaults = defaults["annotation_params"] # NOTE: currently no params used in this step!
    hyperparams = eval(args.hyperparams)
    hyperparams.update({k:v for k,v in hparams_defaults.items() if k not in hyperparams})
    hyperparams = {k:(v if v != "None" else None) for k,v in hyperparams.items()}
    groupparams = eval(args.groupparams)
    groupparams.update({k:v for k,v in gparams_defaults.items() if k not in groupparams})
    groupparams = {k:(v if v != "None" else None) for k,v in groupparams.items()}
    
    # Annotate cell types
    if ct_method == "majority":
        adata_sp = tx.preprocessing.run_majority_voting(adata_sp, adata_sp.uns["spots"])
    elif ct_method == "ssam":
        adata_sp = tx.preprocessing.run_ssam(adata_sp, adata_sp.uns["spots"], adata_sc, um_p_px=hyperparams["um_p_px"])
    else:
        raise ValueError(f"Invalid ct_method: {ct_method}")
    
    # Save output
    pd.DataFrame(data={
        "cell_id": adata_sp.obs_names,
        "celltype": adata_sp.obs[f"ct_{ct_method}"],
        "score": adata_sp.obs[f"ct_{ct_method}_cert"],
    }).to_csv(args.output, index=False)
    
if __name__ == "__main__":
    ct_annotation()
