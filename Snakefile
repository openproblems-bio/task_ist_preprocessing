from TxsimConfig import *

configfile: 'configs/config.yaml'
defaults = 'configs/defaults.yaml'
parsed = ParsedConfig(config, defaults)
final_files = parsed.gen_file_names()

#for f in final_files: print(f)

#Ensures dataset is a name not a file path, and id_code is an int
wildcard_constraints:
    dataset="[^\/]+",
    replicate="[^\/]+",
    id_code="\d+",
    rep_id="\d+"

ruleorder: aggregate_metrics > metric > aggregate_quality_metrics > quality_metric

# Helper function used to get parameters for snakemake rules
def get_params(
    method: str, 
    id_code: int,
    param_name: str
) -> dict:
    """Get parameters for a given method

    Parameters
    ----------
    method : str
        The name of the method
    id_code : int
        The specific ID code
    param_name : str
        The parameter to get. Leave as 'p' for method-specific hyperparameters    

    Returns
    -------
    dict
        Dictionary of parameters
    """
    if parsed.get_method_params(method, id_code) is not None:
        return parsed.get_method_params(method, id_code).get(param_name)
    return None



#Main rule (will always be run)
rule all:
    input:
        final_files

#Rules corresponding to each method
rule pre_segmented:
    conda:
        "envs/txsim-env.yaml"
    input:
        img = lambda w: parsed.get_replicate_file(w.dataset, 'segmented_image', int(w.rep_id)-1)
    output:
        '{results}/{dataset}/replicate{rep_id}/segments_custom-{id_code}.ome.tif',
        '{results}/{dataset}/replicate{rep_id}/areas_custom-{id_code}.csv'
    params:
        hyp = lambda w: get_params('custom', int(w.id_code), 'p'),
        exp = lambda w: get_params('custom', int(w.id_code), 'expand')
    shell:
        "python3 scripts/segment_image.py "
        "-i {input.img} "
        "-p \"{params.hyp}\" "
        "-o {wildcards.results}/{wildcards.dataset}/replicate{wildcards.rep_id} "
        "-s custom "
        "-id {wildcards.id_code} "
        "-e {params.exp} "
        "-b "

rule watershed:
    conda:
        "envs/txsim-env.yaml"
    input:
        img = lambda w: parsed.get_replicate_file(w.dataset, 'images', int(w.rep_id)-1)
    output:
        '{results}/{dataset}/replicate{rep_id}/segments_watershed-{id_code}.ome.tif',
        '{results}/{dataset}/replicate{rep_id}/areas_watershed-{id_code}.csv'
    params:
        hyp = lambda w: get_params('watershed', int(w.id_code), 'p'),
        exp = lambda w: get_params('watershed', int(w.id_code), 'expand')
    shell:
        "python3 scripts/segment_image.py "
        "-i {input.img} "
        "-p \"{params.hyp}\" "
        "-o {wildcards.results}/{wildcards.dataset}/replicate{wildcards.rep_id} "
        "-s watershed "
        "-id {wildcards.id_code} "
        "-e {params.exp} "

rule binning:
    conda:
        "envs/txsim-env.yaml"
    input:
        img = lambda w: parsed.get_replicate_file(w.dataset, 'images', int(w.rep_id)-1)
    output:
        '{results}/{dataset}/replicate{rep_id}/segments_binning-{id_code}.ome.tif',
        '{results}/{dataset}/replicate{rep_id}/areas_binning-{id_code}.csv'
    params:
        hyp = lambda w: get_params('binning', int(w.id_code), 'p')
    shell:
        "python3 scripts/segment_image.py "
        "-i {input.img} "
        "-p \"{params.hyp}\" "
        "-o {wildcards.results}/{wildcards.dataset}/replicate{wildcards.rep_id} "
        "-s binning "
        "-id {wildcards.id_code} "

rule stardist:
    conda:
        "envs/stardist-env.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: 32000 * attempt
    input:
        img = lambda w: parsed.get_replicate_file(w.dataset, 'images', int(w.rep_id)-1)
    output:
        '{results}/{dataset}/replicate{rep_id}/segments_stardist-{id_code}.ome.tif',
        '{results}/{dataset}/replicate{rep_id}/areas_stardist-{id_code}.csv'
    params:
        hyp = lambda w: get_params('stardist', int(w.id_code), 'p'),
        exp = lambda w: get_params('stardist', int(w.id_code), 'expand')
    shell:
        "python3 scripts/segment_image.py "
        "-i {input.img} "
        "-p \"{params.hyp}\" "
        "-o {wildcards.results}/{wildcards.dataset}/replicate{wildcards.rep_id} "
        "-s stardist "
        "-id {wildcards.id_code} "
        "-e {params.exp} "

rule cellpose:
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: 32000 * attempt
    conda: 
        "envs/cellpose-env.yaml"
    input:
        img = lambda w: parsed.get_replicate_file(w.dataset, 'images', int(w.rep_id)-1)
    output:
        '{results}/{dataset}/replicate{rep_id}/segments_cellpose-{id_code}.ome.tif',
        '{results}/{dataset}/replicate{rep_id}/areas_cellpose-{id_code}.csv'
    params:
        hyp = lambda w: get_params('cellpose', int(w.id_code), 'p'),
        exp = lambda w: get_params('cellpose', int(w.id_code), 'expand')
    shell:
        "python3 scripts/segment_image.py "
        "-i {input.img} "
        "-p \"{params.hyp}\" "
        "-o {wildcards.results}/{wildcards.dataset}/replicate{wildcards.rep_id} "
        "-s cellpose "
        "-id {wildcards.id_code} "
        "-e {params.exp} "

rule clustermap:
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: 32000 * attempt + 32000 * (attempt-1)
    conda:
        "envs/clustermap-env.yaml"
    input:
        img = lambda w: parsed.get_replicate_file(w.dataset, 'images', int(w.rep_id)-1),
        mol = lambda w: parsed.get_replicate_file(w.dataset, 'molecules', int(w.rep_id)-1)
    params:
        hyp = lambda w: get_params('clustermap', int(w.id_code), 'p')
    output:
        '{results}/{dataset}/replicate{rep_id}/assignments_clustermap-{id_code}.csv'
    shell:
        "python3 scripts/run_clustermap.py "
        "-m {input.mol} "
        "-p \"{params.hyp}\" "
        "-i {input.img} "
        "-d {wildcards.results}/{wildcards.dataset}/replicate{wildcards.rep_id} "
        "-id {wildcards.id_code} "

rule pciSeq:
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: 32000 * attempt
    conda:
        "envs/pciSeq-env.yaml"
    input:
        '{results}/{dataset}/replicate{rep_id}/segments_{seg}.ome.tif',
        mol = lambda w: parsed.get_replicate_file(w.dataset, 'molecules', int(w.rep_id)-1),
        scd = lambda w: parsed.get_data_file(w.dataset, 'sc_data')
    params:
        hyp = lambda w: get_params('pciSeq', int(w.id_code), 'p')
    output:
        '{results}/{dataset}/replicate{rep_id}/assignments_{seg}_pciSeq-{id_code}.csv',
        '{results}/{dataset}/replicate{rep_id}/celltypes_{seg}_pciSeq-{id_code}.csv'
    shell:
        "python3 scripts/run_pciseq.py "
        "-m {input.mol} "
        "-p \"{params.hyp}\" "
        "-sc {input.scd} "
        "-d {wildcards.results}/{wildcards.dataset}/replicate{wildcards.rep_id} "
        "-s {wildcards.seg} "
        "-id {wildcards.id_code} "

rule basic_assign:
    conda:
        "envs/txsim-env.yaml"
    input:
        '{results}/{dataset}/replicate{rep_id}/segments_{seg}.ome.tif',
        mol = lambda w: parsed.get_replicate_file(w.dataset, 'molecules', int(w.rep_id)-1)
    params:
        hyp = lambda w: get_params('basic', int(w.id_code), 'p')
    output:
        '{results}/{dataset}/replicate{rep_id}/assignments_{seg}_basic-{id_code}.csv'
    shell:
        "python3 scripts/basic_assignment.py "
        "-m {input.mol} "
        "-p \"{params.hyp}\" "
        "-d {wildcards.results}/{wildcards.dataset}/replicate{wildcards.rep_id} "
        "-s {wildcards.seg} "
        "-id {wildcards.id_code} "

rule baysor_prior:
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: 32000 * attempt    
    #conda:
    #    "envs/base-env.yaml"
    container:
        "singularity_container/txsim_baysor_latest.sif"
        #"docker://louisk92/txsim_baysor:latest"
	#"docker://vpetukhov/baysor:master"
    input: 
        '{results}/{dataset}/replicate{rep_id}/segments_{seg}.ome.tif',
        mol = lambda w: parsed.get_replicate_file(w.dataset, 'molecules', int(w.rep_id)-1)
    params:
        hyp = lambda w: get_params('baysor', int(w.id_code), 'p'),
        tmp = "" if config.get('TEMP') is None else f"--temp {config['TEMP']} "
    output:
        '{results}/{dataset}/replicate{rep_id}/assignments_{seg}_baysor-{id_code}.csv',
        '{results}/{dataset}/replicate{rep_id}/areas_{seg}_baysor-{id_code}.csv'
    shell:
        "python3 scripts/run_baysor.py "
        "-m {input.mol} "
        "-p \"{params.hyp}\" "
        "-d {wildcards.results}/{wildcards.dataset}/replicate{wildcards.rep_id} "
        "-id {wildcards.id_code} "
        "-s {wildcards.seg} "
        "{params.tmp}"

rule baysor_no_prior:
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: 32000 * attempt
    #conda:
    #    "envs/base-env.yaml"
    container:
        "singularity_container/txsim_baysor_latest.sif"
        #"docker://louisk92/txsim_baysor:latest"
        #"docker://vpetukhov/baysor:master"
    input:
        mol = lambda w: parsed.get_replicate_file(w.dataset, 'molecules', int(w.rep_id)-1)
    params:
        hyp = lambda w: get_params('baysor', int(w.id_code), 'p'),
        tmp = "" if config.get('TEMP') is None else f"--temp {config['TEMP']} "
    output:
        '{results}/{dataset}/replicate{rep_id}/assignments_baysor-{id_code}.csv',
        '{results}/{dataset}/replicate{rep_id}/areas_baysor-{id_code}.csv'
    shell:
        "python3 scripts/run_baysor.py "
        "-m {input.mol} "
        "-p \"{params.hyp}\" "
        "-d {wildcards.results}/{wildcards.dataset}/replicate{wildcards.rep_id} "
        "-id {wildcards.id_code} "
        "{params.tmp}"

rule normalize_total:
    conda:
        "envs/txsim-env.yaml"
    input:
        '{results}/{dataset}/replicate{rep_id}/assignments_{assign}.csv',
        scd = '{results}/{dataset}/sc_normalized.h5ad'
    params:
        hyp = lambda w: get_params('total', int(w.id_code), 'p'),
        thr = lambda w: get_params('total', int(w.id_code), 'prior_threshold'),
        ct = lambda w: get_params('total', int(w.id_code), 'ct_method'),
        ctthresh = lambda w: get_params('total', int(w.id_code), 'ct_threshold'),
        pergene = lambda w: get_params('total', int(w.id_code), 'per_gene_correction'),
	    pergene_layer = lambda w: get_params('total', int(w.id_code), 'per_gene_layer')

    output:
        '{results}/{dataset}/replicate{rep_id}/counts_{assign}_total-{id_code}.h5ad'
    shell:
        "python3 scripts/gen_counts.py "
        "-as {wildcards.assign} "
        "--singlecell {input.scd} "
        "-d {wildcards.results}/{wildcards.dataset}/replicate{wildcards.rep_id} "
        "-n total "
        "-id {wildcards.id_code} "
        "-p \"{params.hyp}\" "
        "-t {params.thr} "
        "-c {params.ct} "
        "--ctcertthresh {params.ctthresh} "
        "-g {params.pergene} "
        "-l {params.pergene_layer}"

rule normalize_area:
    threads: 1
    conda:
        "envs/txsim-env.yaml"
    input:
        assign = '{results}/{dataset}/replicate{rep_id}/assignments_{method}.csv',
        scd = '{results}/{dataset}/sc_normalized.h5ad'
    params:
        hyp = lambda w: get_params('area', int(w.id_code), 'p'),
        thr = lambda w: get_params('total', int(w.id_code), 'prior_threshold'),
        ct = lambda w: get_params('area', int(w.id_code), 'ct_method'),
        ctthresh = lambda w: get_params('area', int(w.id_code), 'ct_threshold'),
        pergene = lambda w: get_params('total', int(w.id_code), 'per_gene_correction'),
	    pergene_layer = lambda w: get_params('total', int(w.id_code), 'per_gene_layer')

    output:
        '{results}/{dataset}/replicate{rep_id}/counts_{method}_area-{id_code}.h5ad'
    shell:
        "python3 scripts/gen_counts.py "
        "-as {wildcards.method} "
        "--singlecell {input.scd} "
        "-d {wildcards.results}/{wildcards.dataset}/replicate{wildcards.rep_id} "
        "-n area "
        "-id {wildcards.id_code} "
        "-p \"{params.hyp}\" "
        "-t {params.thr} "
        "-c {params.ct} "
        "--ctcertthresh {params.ctthresh} "
        "-g {params.pergene} "
        "-l {params.pergene_layer}"

rule normalize_sc:
    conda:
        "envs/txsim-env.yaml"
    input:
        lambda w: parsed.get_data_file(w.dataset, 'sc_data')
    output:
        '{results}/{dataset}/sc_normalized.h5ad'
    shell:
        "python3 scripts/normalize_sc.py "
        "-sc {input} "
        "-o {output} "

rule metric:
    conda:
        "envs/txsim-env.yaml"
    input:
        '{results}/{dataset}/{replicate}/counts_{methods}.h5ad',
        scd = '{results}/{dataset}/sc_normalized.h5ad'
    output:
        '{results}/{dataset}/{replicate}/metrics_{methods}.csv'
    shell:
        "python3 scripts/calc_metrics.py "
        "-m {wildcards.methods} "
        "-d {wildcards.results}/{wildcards.dataset}/{wildcards.replicate} "
        "-sc {input.scd} "

rule quality_metric:
    conda:
        "envs/txsim-env.yaml"
    input:
        '{results}/{dataset}/{replicate}/counts_{methods}.h5ad',
    output:
        '{results}/{dataset}/{replicate}/quality_metrics_{methods}.csv'
    shell:
        "python3 scripts/calc_quality_metrics.py "
        "-m {wildcards.methods} "
        "-d {wildcards.results}/{wildcards.dataset}/{wildcards.replicate} "

rule group_metrics:
    conda:
        "envs/txsim-env.yaml"
    input:
        parsed.get_metric_inputs #function that can accept wildcards as input argument
    output:
        '{results}/{dataset}/{replicate}/group_metrics.csv'
    shell:
        "python3 scripts/calc_group_metrics.py "
        "-f all "
        "-d {wildcards.results}/{wildcards.dataset}/{wildcards.replicate} "

rule aggregate_counts:
    conda:
        "envs/txsim-env.yaml"
    input:
        lambda w : parsed.get_aggregate_inputs(w, 'counts') #function can accept wildcards as input argument
    output:
        '{results}/{dataset}/aggregated/counts_{method}.h5ad'
    shell:
        "python3 scripts/aggregate_counts.py "
        "-m {wildcards.method} "
        "-d {wildcards.results}/{wildcards.dataset} "

rule aggregate_metrics:
    conda:
        "envs/txsim-env.yaml"
    input:
        lambda w : parsed.get_aggregate_inputs(w, 'metrics') 
    output:
        '{results}/{dataset}/aggregated/aggregated_metrics_{method}.csv'
    shell:
        "python3 scripts/aggregate_metrics.py "
        "-m {wildcards.method} "
        "-d {wildcards.results}/{wildcards.dataset} "
        "-t metrics"

rule aggregate_quality_metrics:
    conda:
        "envs/txsim-env.yaml"
    input:
        lambda w : parsed.get_aggregate_inputs(w, 'quality_metrics') 
    output:
        '{results}/{dataset}/aggregated/aggregated_quality_metrics_{method}.csv'
    shell:
        "python3 scripts/aggregate_metrics.py "
        "-m {wildcards.method} "
        "-d {wildcards.results}/{wildcards.dataset} "
        "-t quality_metrics "

rule aggregate_group_metrics:
    conda:
        "envs/txsim-env.yaml"
    input:
        lambda w : parsed.get_aggregate_inputs(w, 'group_metrics') 
    output:
        '{results}/{dataset}/aggregated/aggregated_group_metrics.csv'
    shell:
        "python3 scripts/aggregate_group_metrics.py "
        "-d {wildcards.results}/{wildcards.dataset} "


