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
    params_dict = parsed.get_method_params(method, id_code)
    if (params_dict is not None) and (params_dict != {}):
        return parsed.get_method_params(method, id_code).get(param_name)
    return {}



#Main rule (will always be run)
rule all:
    input:
        final_files


#####################
# rules for pre-run #
#####################

rule get_tile_info:
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: 64000 * attempt
    conda:
        "envs/txsim-env.yaml"
    input:
        img = lambda w: parsed.get_replicate_file(w.dataset, 'images', int(w.rep_id)-1),
        mol = lambda w: parsed.get_replicate_file(w.dataset, 'molecules', int(w.rep_id)-1)
    output:
        '{results}/{dataset}/replicate{rep_id}/tile_info.csv',
    params:
        n_spots_baysor = f"{config['MAX_N_SPOTS_PER_TILE_BAYSOR']}",
        n_spots_clustermap = f"{config['MAX_N_SPOTS_PER_TILE_CLUSTERMAP']}",
        extend_n_pixel = f"{config['N_PIXEL_EXTEND_TILE']}"
    shell:
        "python3 scripts/retrieve_tiles_info.py "
        "-i {input.img} "
        "-m {input.mol} "
        "-o {output} "
        "-nb {params.n_spots_baysor} "
        "-nc {params.n_spots_clustermap} "
        "-ne {params.extend_n_pixel}"
            

################
# Segmentation #
################

#Rules corresponding to each method
rule pre_segmented:
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: 64000 * attempt
    conda:
        "envs/txsim-env.yaml"
    input:
        img = lambda w: parsed.get_replicate_file(w.dataset, 'segmented_image', int(w.rep_id)-1)
    output:
        '{results}/{dataset}/replicate{rep_id}/segments_custom-{id_code}.ome.tif',
        '{results}/{dataset}/replicate{rep_id}/areas_custom-{id_code}.csv'
    params:
        hyper_params = lambda w: get_params('custom', int(w.id_code), 'hyper_params'),
        group_params = lambda w: get_params('custom', int(w.id_code), 'group_params')
    shell:
        "python3 scripts/segment_image.py "
        "-i {input.img} "
        "-p \"{params.hyper_params}\" "
        "-g \"{params.group_params}\" "
        "-o {wildcards.results}/{wildcards.dataset}/replicate{wildcards.rep_id} "
        "-s custom "
        "-id {wildcards.id_code} "
        "-b "

rule watershed:
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: 64000 * attempt
    conda:
        "envs/txsim-env.yaml"
    input:
        img = lambda w: parsed.get_replicate_file(w.dataset, 'images', int(w.rep_id)-1)
    output:
        '{results}/{dataset}/replicate{rep_id}/segments_watershed-{id_code}.ome.tif',
        '{results}/{dataset}/replicate{rep_id}/areas_watershed-{id_code}.csv'
    params:
        hyper_params = lambda w: get_params('watershed', int(w.id_code), 'hyper_params'),
        group_params = lambda w: get_params('watershed', int(w.id_code), 'group_params')
    shell:
        "python3 scripts/segment_image.py "
        "-i {input.img} "
        "-p \"{params.hyper_params}\" "
        "-g \"{params.group_params}\" "
        "-o {wildcards.results}/{wildcards.dataset}/replicate{wildcards.rep_id} "
        "-s watershed "
        "-id {wildcards.id_code} "

rule binning:
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: 64000 * attempt
    conda:
        "envs/txsim-env.yaml"
    input:
        img = lambda w: parsed.get_replicate_file(w.dataset, 'images', int(w.rep_id)-1)
    output:
        '{results}/{dataset}/replicate{rep_id}/segments_binning-{id_code}.ome.tif',
        '{results}/{dataset}/replicate{rep_id}/areas_binning-{id_code}.csv'
    params:
        hyper_params = lambda w: get_params('binning', int(w.id_code), 'hyper_params'),
        group_params = lambda w: get_params('binning', int(w.id_code), 'group_params')
    shell:
        "python3 scripts/segment_image.py "
        "-i {input.img} "
        "-p \"{params.hyper_params}\" "
        "-g \"{params.group_params}\" "
        "-o {wildcards.results}/{wildcards.dataset}/replicate{wildcards.rep_id} "
        "-s binning "
        "-id {wildcards.id_code} "

rule stardist:
    conda:
        "envs/stardist-env.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: 64000 * attempt
    input:
        img = lambda w: parsed.get_replicate_file(w.dataset, 'images', int(w.rep_id)-1)
    output:
        '{results}/{dataset}/replicate{rep_id}/segments_stardist-{id_code}.ome.tif',
        '{results}/{dataset}/replicate{rep_id}/areas_stardist-{id_code}.csv'
    params:
        hyper_params = lambda w: get_params('stardist', int(w.id_code), 'hyper_params'),
        group_params = lambda w: get_params('stardist', int(w.id_code), 'group_params')
    shell:
        "python3 scripts/segment_image.py "
        "-i {input.img} "
        "-p \"{params.hyper_params}\" "
        "-g \"{params.group_params}\" "
        "-o {wildcards.results}/{wildcards.dataset}/replicate{wildcards.rep_id} "
        "-s stardist "
        "-id {wildcards.id_code} "

rule cellpose:
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: 64000 * attempt
    conda: 
        "envs/cellpose-env.yaml"
    input:
        img = lambda w: parsed.get_replicate_file(w.dataset, 'images', int(w.rep_id)-1)
    output:
        '{results}/{dataset}/replicate{rep_id}/segments_cellpose-{id_code}.ome.tif',
        '{results}/{dataset}/replicate{rep_id}/areas_cellpose-{id_code}.csv'
    params:
        hyper_params = lambda w: get_params('cellpose', int(w.id_code), 'hyper_params'),
        group_params = lambda w: get_params('cellpose', int(w.id_code), 'group_params')
    shell:
        "python3 scripts/segment_image.py "
        "-i {input.img} "
        "-p \"{params.hyper_params}\" "
        "-g \"{params.group_params}\" "
        "-o {wildcards.results}/{wildcards.dataset}/replicate{wildcards.rep_id} "
        "-s cellpose "
        "-id {wildcards.id_code} "

rule mesmer:
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: 64000 * attempt
    conda: 
        "envs/mesmer-env.yaml"
    input:
        img = lambda w: parsed.get_replicate_file(w.dataset, 'images', int(w.rep_id)-1)
    output:
        '{results}/{dataset}/replicate{rep_id}/segments_mesmer-{id_code}.ome.tif',
        '{results}/{dataset}/replicate{rep_id}/areas_mesmer-{id_code}.csv'
    params:
        hyper_params = lambda w: get_params('cellpose', int(w.id_code), 'hyper_params'),
        group_params = lambda w: get_params('cellpose', int(w.id_code), 'group_params')
    shell:
        "python3 scripts/run_mesmer.py "
        "-i {input.img} "
        "-p \"{params.hyper_params}\" "
        "-g \"{params.group_params}\" "
        "-o {wildcards.results}/{wildcards.dataset}/replicate{wildcards.rep_id} "
        "-id {wildcards.id_code} "

#rule clustermap:
#    threads: 8
#    resources:
#        mem_mb = lambda wildcards, attempt: 200000 + 32000 * (attempt-1)  #64000 * attempt + 32000 * (attempt-1)
#    conda:
#        "envs/clustermap-env.yaml"
#    input:
#        img = lambda w: parsed.get_replicate_file(w.dataset, 'images', int(w.rep_id)-1),
#        mol = lambda w: parsed.get_replicate_file(w.dataset, 'molecules', int(w.rep_id)-1)
#    params:
#        hyper_params = lambda w: get_params('clustermap', int(w.id_code), 'hyper_params')
#    output:
#        '{results}/{dataset}/replicate{rep_id}/assignments_clustermap-{id_code}.csv'
#    shell:
#        "python3 scripts/run_clustermap.py "
#        "-m {input.mol} "
#        "-p \"{params.hyper_params}\" "
#        "-i {input.img} "
#        "-d {wildcards.results}/{wildcards.dataset}/replicate{wildcards.rep_id} "
#        "-id {wildcards.id_code} "

rule pciseq:
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: 64000 * attempt
    conda:
        "envs/pciseq-env.yaml"
    input:
        '{results}/{dataset}/replicate{rep_id}/segments_{seg}.ome.tif',
        mol = lambda w: parsed.get_replicate_file(w.dataset, 'molecules', int(w.rep_id)-1),
        scd = '{results}/{dataset}/sc_normalized.h5ad'
        #scd = lambda w: parsed.get_data_file(w.dataset, 'sc_data')
    params:
        hyper_params = lambda w: get_params('pciseq', int(w.id_code), 'hyper_params')
    output:
        '{results}/{dataset}/replicate{rep_id}/assignments_{seg}_pciseq-{id_code}.csv',
        '{results}/{dataset}/replicate{rep_id}/celltypes_{seg}_pciseq-{id_code}.csv'
    shell:
        "python3 scripts/run_pciseq.py "
        "-m {input.mol} "
        "-p \"{params.hyper_params}\" "
        "-sc {input.scd} "
        "-d {wildcards.results}/{wildcards.dataset}/replicate{wildcards.rep_id} "
        "-s {wildcards.seg} "
        "-id {wildcards.id_code} "


###################
# Spot assignment #
###################

rule basic_assign:
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: 64000 * attempt
    conda:
        "envs/txsim-env.yaml"
    input:
        '{results}/{dataset}/replicate{rep_id}/segments_{seg}.ome.tif',
        mol = lambda w: parsed.get_replicate_file(w.dataset, 'molecules', int(w.rep_id)-1)
    params:
        hyper_params = lambda w: get_params('basic', int(w.id_code), 'hyper_params')
    output:
        '{results}/{dataset}/replicate{rep_id}/assignments_{seg}_basic-{id_code}.csv'
    shell:
        "python3 scripts/basic_assignment.py "
        "-m {input.mol} "
        "-p \"{params.hyper_params}\" "
        "-d {wildcards.results}/{wildcards.dataset}/replicate{wildcards.rep_id} "
        "-s {wildcards.seg} "
        "-id {wildcards.id_code} "


rule generate_tile_of_input_data:
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: 64000 * attempt
    conda:
        "envs/txsim-env.yaml"
    input:
        img = lambda w: parsed.get_replicate_file(w.dataset, 'images', int(w.rep_id)-1),
        mol = lambda w: parsed.get_replicate_file(w.dataset, 'molecules', int(w.rep_id)-1)
    output:
        img_out = '{results}/{dataset}/replicate{rep_id}/dapi_ny{y_tiles}_nx{x_tiles}_{y_id}_{x_id}_px{n_expand_px}.ome.tif',
        mol_out = '{results}/{dataset}/replicate{rep_id}/spots_ny{y_tiles}_nx{x_tiles}_{y_id}_{x_id}_px{n_expand_px}.csv',
    shell:
        "python3 scripts/generate_tile.py "
        "-i {input.img} "
        "-m {input.mol} "
        "-t {wildcards.y_tiles}-{wildcards.x_tiles}-{wildcards.y_id}-{wildcards.x_id}-{wildcards.n_expand_px} "
        "--output_img {output.img_out} "
        "--output_mol {output.mol_out}"


rule generate_tile_of_segmentation:
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: 64000 * attempt
    conda:
        "envs/txsim-env.yaml"
    input:
        '{results}/{dataset}/replicate{rep_id}/segments_{seg}.ome.tif'
    output:
        '{results}/{dataset}/replicate{rep_id}/segments_{seg}_ny{y_tiles}_nx{x_tiles}_{y_id}_{x_id}_px{n_expand_px}.ome.tif'
    shell:
        "python3 scripts/generate_tile.py "
        "-i {input} "
        "-t {wildcards.y_tiles}-{wildcards.x_tiles}-{wildcards.y_id}-{wildcards.x_id}-{wildcards.n_expand_px} "
        "--output_img {output}"


rule clustermap_tile:
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: 200000 + 32000 * (attempt-1)  #64000 * attempt + 32000 * (attempt-1)
    conda:
        "envs/clustermap-env.yaml"
    input: 
        img = '{results}/{dataset}/replicate{rep_id}/dapi_ny{y_tiles}_nx{x_tiles}_{y_id}_{x_id}_px{n_expand_px}.ome.tif',
        mol = '{results}/{dataset}/replicate{rep_id}/spots_ny{y_tiles}_nx{x_tiles}_{y_id}_{x_id}_px{n_expand_px}.csv'
    params:
        hyper_params = lambda w: get_params('clustermap', int(w.id_code), 'hyper_params')
    output:
        '{results}/{dataset}/replicate{rep_id}/assignments_clustermap-{id_code}_ny{y_tiles}_nx{x_tiles}_{y_id}_{x_id}_px{n_expand_px}.csv'
    shell:
        "python3 scripts/run_clustermap.py "
        "-m {input.mol} "
        "-o {output} "
        "-p \"{params.hyper_params}\" "
        "-i {input.img} "
        "-d {wildcards.results}/{wildcards.dataset}/replicate{wildcards.rep_id} "
        "-id {wildcards.id_code}"


rule baysor_prior_tile:
    threads: 22
    resources:
        mem_mb = lambda wildcards, attempt: 200000 + 32000 * attempt
    #conda:
    #    "envs/base-env.yaml"
    container:
        "singularity_container/baysor_v0.6.2bin.sif"
        #"docker://louisk92/txsim_baysor:v0.6.2bin"
    input: 
        img = '{results}/{dataset}/replicate{rep_id}/segments_{seg}_ny{y_tiles}_nx{x_tiles}_{y_id}_{x_id}_px{n_expand_px}.ome.tif',
        mol = '{results}/{dataset}/replicate{rep_id}/spots_ny{y_tiles}_nx{x_tiles}_{y_id}_{x_id}_px{n_expand_px}.csv'
    params:
        hyper_params = lambda w: get_params('baysor', int(w.id_code), 'hyper_params'),
        tmp = f"{config['TEMP']}"
    output:
        assign = '{results}/{dataset}/replicate{rep_id}/assignments_{seg}_baysor-{id_code}_ny{y_tiles}_nx{x_tiles}_{y_id}_{x_id}_px{n_expand_px}.csv',
        areas = '{results}/{dataset}/replicate{rep_id}/areas_{seg}_baysor-{id_code}_ny{y_tiles}_nx{x_tiles}_{y_id}_{x_id}_px{n_expand_px}.csv'
    shell:
        "export JULIA_NUM_THREADS=20; python3 scripts/run_baysor.py "
        "-m {input.mol} "
        "-o {output.assign} "
        "-p \"{params.hyper_params}\" "
        "-d {wildcards.results}/{wildcards.dataset}/replicate{wildcards.rep_id} "
        "-id {wildcards.id_code} "
        "-s {wildcards.seg} "
        "--temp {params.tmp}/{wildcards.dataset}/rep{wildcards.rep_id}/{wildcards.seg}_baysor-{wildcards.id_code}-{wildcards.y_tiles}-{wildcards.x_tiles}-{wildcards.y_id}-{wildcards.x_id}"

rule baysor_no_prior_tile:
    threads: 22
    resources:
        mem_mb = lambda wildcards, attempt: 200000 + 32000 * attempt
    #conda:
    #    "envs/base-env.yaml"
    container:
        "singularity_container/baysor_v0.6.2bin.sif"
        #"docker://louisk92/txsim_baysor:v0.6.2bin"
    input: 
        mol = '{results}/{dataset}/replicate{rep_id}/spots_ny{y_tiles}_nx{x_tiles}_{y_id}_{x_id}_px{n_expand_px}.csv'
    params:
        hyper_params = lambda w: get_params('baysor', int(w.id_code), 'hyper_params'),
        tmp = f"{config['TEMP']}"
    output:
        assign = '{results}/{dataset}/replicate{rep_id}/assignments_baysor-{id_code}_ny{y_tiles}_nx{x_tiles}_{y_id}_{x_id}_px{n_expand_px}.csv',
        areas = '{results}/{dataset}/replicate{rep_id}/areas_baysor-{id_code}_ny{y_tiles}_nx{x_tiles}_{y_id}_{x_id}_px{n_expand_px}.csv'
    shell:
        "export JULIA_NUM_THREADS=20; python3 scripts/run_baysor.py "
        "-m {input.mol} "
        "-o {output.assign} "
        "-p \"{params.hyper_params}\" "
        "-d {wildcards.results}/{wildcards.dataset}/replicate{wildcards.rep_id} "
        "-id {wildcards.id_code} "
        "--temp {params.tmp}/{wildcards.dataset}/rep{wildcards.rep_id}/baysor-{wildcards.id_code}-{wildcards.y_tiles}-{wildcards.x_tiles}-{wildcards.y_id}-{wildcards.x_id}"


rule aggregate_baysor_no_prior_tile_assignments:
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: 64000 * attempt
    conda:
        "envs/txsim-env.yaml"
    input:
        lambda w : parsed.get_tiles_aggregation_input_files(w.dataset, w.rep_id, "baysor", w.id_code, areas=True, seg=None),
        img = lambda w: parsed.get_replicate_file(w.dataset, 'images', int(w.rep_id)-1),
        mol = lambda w: parsed.get_replicate_file(w.dataset, 'molecules', int(w.rep_id)-1)
        #'{results}/{dataset}/replicate{rep_id}/assignments_{method}_ny{y_tiles}_nx{x_tiles}_{y_id}_{x_id}_px{n_expand_px}.csv',
        #'{results}/{dataset}/replicate{rep_id}/areas_{method}_ny{y_tiles}_nx{x_tiles}_{y_id}_{x_id}_px{n_expand_px}.csv'
    output:
        '{results}/{dataset}/replicate{rep_id}/assignments_baysor-{id_code}.csv',
        '{results}/{dataset}/replicate{rep_id}/areas_baysor-{id_code}.csv'
    shell:
        "python3 scripts/aggregate_tile_assignments.py "
        "-p {input} "
        "-o {output} "
        "-m {input.mol} "
        "-i {input.img}"

use rule aggregate_baysor_no_prior_tile_assignments as aggregate_baysor_prior_tile_assignments with:
    input:
        lambda w : parsed.get_tiles_aggregation_input_files(w.dataset, w.rep_id, "baysor", w.id_code, areas=True, seg=w.seg),
        img = lambda w: parsed.get_replicate_file(w.dataset, 'images', int(w.rep_id)-1),
        mol = lambda w: parsed.get_replicate_file(w.dataset, 'molecules', int(w.rep_id)-1)
    output:
        '{results}/{dataset}/replicate{rep_id}/assignments_{seg}_baysor-{id_code}.csv',
        '{results}/{dataset}/replicate{rep_id}/areas_{seg}_baysor-{id_code}.csv'

use rule aggregate_baysor_no_prior_tile_assignments as aggregate_clustermap_tile_assignments with:
    input:
        lambda w : parsed.get_tiles_aggregation_input_files(w.dataset, w.rep_id, "clustermap", w.id_code, areas=False, seg=None),
        img = lambda w: parsed.get_replicate_file(w.dataset, 'images', int(w.rep_id)-1),
        mol = lambda w: parsed.get_replicate_file(w.dataset, 'molecules', int(w.rep_id)-1)
    output:
        '{results}/{dataset}/replicate{rep_id}/assignments_clustermap-{id_code}.csv',
        # no areas are aggregated (for clustermap they are calculated later)


# TODO: We can probably delete them, the full img case is captured with the n=1 tiling.
#rule baysor_prior:
#    threads: 22
#    resources:
#        mem_mb = lambda wildcards, attempt: 1000 # 200000 + 64000 * attempt
#    #conda:
#    #    "envs/base-env.yaml"
#    container:
#        "singularity_container/baysor_v0.6.2bin.sif"
#        #"docker://louisk92/txsim_baysor:v0.6.2bin"
#    input: 
#        '{results}/{dataset}/replicate{rep_id}/segments_{seg}.ome.tif',
#        mol = lambda w: parsed.get_replicate_file(w.dataset, 'molecules', int(w.rep_id)-1)
#    params:
#        hyper_params = lambda w: get_params('baysor', int(w.id_code), 'hyper_params'),
#        tmp = f"{config['TEMP']}"
#    output:
#        '{results}/{dataset}/replicate{rep_id}/assignments_{seg}_baysor-{id_code}.csv',
#        '{results}/{dataset}/replicate{rep_id}/areas_{seg}_baysor-{id_code}.csv'
#    shell:
#        "export JULIA_NUM_THREADS=20; python3 scripts/run_baysor.py "
#        "-m {input.mol} "
#        "-p \"{params.hyper_params}\" "
#        "-d {wildcards.results}/{wildcards.dataset}/replicate{wildcards.rep_id} "
#        "-id {wildcards.id_code} "
#        "-s {wildcards.seg} "
#        "--temp {params.tmp}/{wildcards.dataset}/rep{wildcards.rep_id}/{wildcards.seg}_baysor-{wildcards.id_code}"

#rule baysor_no_prior:
#    threads: 22
#    resources:
#        mem_mb = lambda wildcards, attempt: 1000 #200000 + 64000 * attempt
#    #conda:
#    #    "envs/base-env.yaml"
#    container:
#        "singularity_container/baysor_v0.6.2bin.sif"
#        #"docker://louisk92/txsim_baysor:v0.6.2bin"
#    input:
#        mol = lambda w: parsed.get_replicate_file(w.dataset, 'molecules', int(w.rep_id)-1)
#    params:
#        hyper_params = lambda w: get_params('baysor', int(w.id_code), 'hyper_params'),
#        tmp = f"{config['TEMP']}"
#    output:
#        '{results}/{dataset}/replicate{rep_id}/assignments_baysor-{id_code}.csv',
#        '{results}/{dataset}/replicate{rep_id}/areas_baysor-{id_code}.csv'
#    shell:
#        "export JULIA_NUM_THREADS=20; python3 scripts/run_baysor.py "
#        "-m {input.mol} "
#        "-p \"{params.hyper_params}\" "
#        "-d {wildcards.results}/{wildcards.dataset}/replicate{wildcards.rep_id} "
#        "-id {wildcards.id_code} "
#        "--temp {params.tmp}/{wildcards.dataset}/rep{wildcards.rep_id}/baysor-{wildcards.id_code}"


#################
# Normalization #
#################

rule normalize_total:
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: 64000 * attempt
    conda:
        "envs/txsim-env.yaml"
    input:
        '{results}/{dataset}/replicate{rep_id}/assignments_{assign}.csv',
        scd = '{results}/{dataset}/sc_normalized.h5ad'
    params:
        hyper_params = lambda w: get_params('total', int(w.id_code), 'hyper_params'),
        group_params = lambda w: get_params('total', int(w.id_code), 'group_params')
        # thr = lambda w: get_params('total', int(w.id_code), 'prior_threshold'),
        # ct = lambda w: get_params('total', int(w.id_code), 'ct_method'),
        # ctthresh = lambda w: get_params('total', int(w.id_code), 'ct_threshold'),
        # pergene = lambda w: get_params('total', int(w.id_code), 'per_gene_correction'),
	    # pergene_layer = lambda w: get_params('total', int(w.id_code), 'per_gene_layer')

    output:
        '{results}/{dataset}/replicate{rep_id}/normcounts_{assign}_total-{id_code}.h5ad'
    shell:
        "python3 scripts/gen_counts.py "
        "-as {wildcards.assign} "
        "--singlecell {input.scd} "
        "-d {wildcards.results}/{wildcards.dataset}/replicate{wildcards.rep_id} "
        "-n total "
        "-id {wildcards.id_code} "
        "-p \"{params.hyper_params}\" "
        "-g \"{params.group_params}\" "
        # "-t {params.thr} "
        # "-c {params.ct} "
        # "--ctcertthresh {params.ctthresh} "
        # "-g {params.pergene} "
        # "-l {params.pergene_layer}"

rule normalize_area:
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: 64000 * attempt
    conda:
        "envs/txsim-env.yaml"
    input:
        assign = '{results}/{dataset}/replicate{rep_id}/assignments_{method}.csv',
        scd = '{results}/{dataset}/sc_normalized.h5ad'
    params:
        hyper_params = lambda w: get_params('area', int(w.id_code), 'hyper_params'),
        group_params = lambda w: get_params('area', int(w.id_code), 'group_params')
        # thr = lambda w: get_params('total', int(w.id_code), 'prior_threshold'),
        # ct = lambda w: get_params('area', int(w.id_code), 'ct_method'),
        # ctthresh = lambda w: get_params('area', int(w.id_code), 'ct_threshold'),
        # pergene = lambda w: get_params('total', int(w.id_code), 'per_gene_correction'),
	    # pergene_layer = lambda w: get_params('total', int(w.id_code), 'per_gene_layer')

    output:
        '{results}/{dataset}/replicate{rep_id}/normcounts_{method}_area-{id_code}.h5ad'
    shell:
        "python3 scripts/gen_counts.py "
        "-as {wildcards.method} "
        "--singlecell {input.scd} "
        "-d {wildcards.results}/{wildcards.dataset}/replicate{wildcards.rep_id} "
        "-n area "
        "-id {wildcards.id_code} "
        "-p \"{params.hyper_params}\" "
        "-g \"{params.group_params}\" "
        # "-t {params.thr} "
        # "-c {params.ct} "
        # "--ctcertthresh {params.ctthresh} "
        # "-g {params.pergene} "
        # "-l {params.pergene_layer}"


########################
# Cell type annotation #
########################

# Rule for adding cell type annotations (from csvs) to the spatial adata and applying corrections according hparams

rule annotate_counts:
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: 64000 * attempt
    conda:
        "envs/txsim-env.yaml"
    input:
        counts = '{results}/{dataset}/replicate{rep_id}/normcounts_{method}.h5ad',
        scd = '{results}/{dataset}/sc_normalized.h5ad',
        ct_csv = '{results}/{dataset}/replicate{rep_id}/celltypes_{method}_{ct_method}-{id_code}.csv'
    params:
        hyper_params = lambda w: get_params(w.ct_method, int(w.id_code), 'hyper_params'),
        group_params = lambda w: get_params(w.ct_method, int(w.id_code), 'group_params')
    output:
        '{results}/{dataset}/replicate{rep_id}/counts_{method}_{ct_method}-{id_code}.h5ad'
    shell:
        "python3 scripts/annotate_counts.py "
        "-c {wildcards.method} "
        "--singlecell {input.scd} "
        "-n {input.ct_csv} "
        "-d {wildcards.results}/{wildcards.dataset}/replicate{wildcards.rep_id} "
        "-a {wildcards.ct_method} "
        "-id {wildcards.id_code} "
        "-p \"{params.hyper_params}\" "
        "-g \"{params.group_params}\" "


# csv generating methods

rule annotate_celltypes_txsim_methods:
    wildcard_constraints:
        ct_method="majority|ssam"
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: 64000 * attempt
    conda:
        "envs/txsim-env.yaml"
    input:
        counts = '{results}/{dataset}/replicate{rep_id}/normcounts_{method}.h5ad',
        scd = '{results}/{dataset}/sc_normalized.h5ad'
    params:
        hyper_params = lambda w: get_params(w.ct_method, int(w.id_code), 'hyper_params'),
        group_params = lambda w: get_params(w.ct_method, int(w.id_code), 'group_params')
    output:
        '{results}/{dataset}/replicate{rep_id}/celltypes_{method}_{ct_method}-{id_code}.csv'
    shell:
        "python3 scripts/annotate_celltypes.py "
        "-s {input.counts} "
        "-d {input.scd} "
        "-o {output} "
        "-m {wildcards.ct_method} "
        "-p \"{params.hyper_params}\" "
        "-g \"{params.group_params}\" "

rule annotate_celltypes_pciseqct: #NOTE: probably we'll never need hparams in this rule. Otherwise add this option.
    #NOTE: in output tried to use only {seg}_{assign} instead of {seg}_{assign}_{norm} with more flexible 
    # "assign" but it didn't work even though wildcard_constraints make sense to me...
    #wildcard_constraints:
    #    seg="^[^-_]*-[^-_]*$",
    #    assign=".*"
    #TODO: The current id_code of pciseqct and pciseq is the same. Probably we need a hyperparameter for pciseqct
    #      to refer to an id_code of pciseq.... Currently the assumption is we run these methods only once... code=0
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: 64000 * attempt
    conda:
        "envs/txsim-env.yaml"
    input:
        ann = '{results}/{dataset}/replicate{rep_id}/celltypes_{seg}_pciseq-{id_code}.csv',
        counts = '{results}/{dataset}/replicate{rep_id}/normcounts_{seg}_{assign}_{norm}.h5ad',
    output:
        '{results}/{dataset}/replicate{rep_id}/celltypes_{seg}_{assign}_{norm}_pciseqct-{id_code}.csv'
    shell:
        "python3 scripts/copy_celltypes_pciseqct.py -i {input.ann} -s {input.counts} -o {output}"

rule annotate_celltypes_mfishtools:
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: 64000 * attempt
    container:
        "docker://louisk92/txsim_mfishtools:2024-04-18" 
    input:
        counts = '{results}/{dataset}/replicate{rep_id}/normcounts_{method}.h5ad',
        scd = '{results}/{dataset}/sc_normalized.h5ad',
    output:
        '{results}/{dataset}/replicate{rep_id}/celltypes_{method}_mfishtools-{id_code}.csv'
    params:
        hyper_params = lambda w: get_params('mfishtools', int(w.id_code), 'hyper_params'),
        group_params = lambda w: get_params('mfishtools', int(w.id_code), 'group_params')
    shell:
        "Rscript scripts/annotate_celltypes_mfishtools.r " 
        "-s {input.counts} "
        "-d {input.scd} "
        "-o {output} "
        "-p \"{params.hyper_params}\" "
        "-g \"{params.group_params}\" "
        
rule annotate_celltypes_frmatch:
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: 64000 * attempt
    container:
        "docker://louisk92/txsim_frmatch:2024-04-18" 
    input:
        counts = '{results}/{dataset}/replicate{rep_id}/normcounts_{method}.h5ad',
        scd = '{results}/{dataset}/sc_normalized.h5ad',
    output:
        '{results}/{dataset}/replicate{rep_id}/celltypes_{method}_frmatch-{id_code}.csv'
    params:
        hyper_params = lambda w: get_params('frmatch', int(w.id_code), 'hyper_params'),
        group_params = lambda w: get_params('frmatch', int(w.id_code), 'group_params')
    shell:
        "Rscript scripts/annotate_celltypes_FRmatch.r " 
        "-s {input.counts} "
        "-d {input.scd} "
        "-o {output} "
        "-p \"{params.hyper_params}\" "
        "-g \"{params.group_params}\" "

rule annotate_celltypes_scrattchmapping:
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: 64000 * attempt
    container:
        "docker://louisk92/txsim_scrattchmapping:2024-04-19"
    input:
        counts = '{results}/{dataset}/replicate{rep_id}/normcounts_{method}.h5ad',
        scd = '{results}/{dataset}/sc_normalized.h5ad',
    output:
        '{results}/{dataset}/replicate{rep_id}/celltypes_{method}_scrattchmapping-{id_code}.csv'
    params:
        hyper_params = lambda w: get_params('scrattchmapping', int(w.id_code), 'hyper_params'),
        group_params = lambda w: get_params('scrattchmapping', int(w.id_code), 'group_params'),
        tmp = f"{config['TEMP']}"
    shell:
        "mkdir -p {params.tmp}/{wildcards.dataset}/rep{wildcards.rep_id}/{wildcards.method}_scrattchmapping-{wildcards.id_code} && "
        "Rscript scripts/annotate_celltypes_scrattchmapping.r " 
        "-s {input.counts} "
        "-d {input.scd} "
        "-o {output} "
        "-p \"{params.hyper_params}\" "
        "-g \"{params.group_params}\" "
        "-t {params.tmp}/{wildcards.dataset}/rep{wildcards.rep_id}/{wildcards.method}_scrattchmapping-{wildcards.id_code} "

rule annotate_celltypes_tangram:
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: 64000 * attempt
    conda:
        "envs/tangram-env.yaml"
    input:
        counts = '{results}/{dataset}/replicate{rep_id}/normcounts_{method}.h5ad',
        scd = '{results}/{dataset}/sc_normalized.h5ad',
    output:
        '{results}/{dataset}/replicate{rep_id}/celltypes_{method}_tangram-{id_code}.csv'
    params:
        hyper_params = lambda w: get_params('tangram', int(w.id_code), 'hyper_params'),
        group_params = lambda w: get_params('tangram', int(w.id_code), 'group_params')
    shell:
        "python3 scripts/annotate_celltypes_tangram.py " 
        "-s {input.counts} "
        "-d {input.scd} "
        "-o {output} "
        "-p \"{params.hyper_params}\" "
        "-g \"{params.group_params}\" "

rule annotate_celltypes_mapmycells:
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: 64000 * attempt
    conda:
        "envs/mapmycells-env.yaml"
    input:
        counts = '{results}/{dataset}/replicate{rep_id}/normcounts_{method}.h5ad',
        scd = '{results}/{dataset}/sc_normalized.h5ad',
    output:
        '{results}/{dataset}/replicate{rep_id}/celltypes_{method}_mapmycells-{id_code}.csv'
    params:
        hyper_params = lambda w: get_params('mapmycells', int(w.id_code), 'hyper_params'),
        group_params = lambda w: get_params('mapmycells', int(w.id_code), 'group_params')
    shell:
        "python3 scripts/annotate_celltypes_mapmycells.py " 
        "-s {input.counts} "
        "-d {input.scd} "
        "-o {output} "
        "-p \"{params.hyper_params}\" "
        "-g \"{params.group_params}\" "
        "-t {params.tmp}/{wildcards.dataset}/rep{wildcards.rep_id}/{wildcards.method}_mapmycells-{wildcards.id_code} "

# Consensus methods
#TODO: Probably you want to combine different consensus methods into one rule and have a hyperparameter.

def input_files_for_consensus_annotation(id_code, results, dataset, rep_id, method, consensus):
    """ Get cell type annotation input csvs for rule for consensus annotation
    
    Arguments
    ---------
    id_code: int
        ID that defines the parameters for the consensus annotation
    results: Results folder
    dataset: data set name
    rep_id: ID of replicate
    method: String that describes previously ran steps
    consensus: Type of consensus aggregation (supported: "nwconsensus", "gmconsensus")
    """
    hparams = get_params(consensus, int(id_code), 'hyper_params')
    ct_method_ids = hparams.get("ids")
    ct_method_ids = ct_method_ids.split('-')
    ct_methods = hparams.get("methods")
    ct_methods = ct_methods.split('-')
    file_paths = [
        f"{results}/{dataset}/replicate{rep_id}/celltypes_{method}_{ct_method}-{ct_method_id}.csv"
        for (ct_method_id, ct_method) in zip(ct_method_ids, ct_methods)
    ]    
    ## This is very much unnecessary, but I got a "python does not work anymore"-level error here and can't pinpoint it
    ## down. Think in future fresh envs this can be deleted. It worked before I had some deep conda cache issues.
    #file_paths = [f.replace(" ", "") for f in file_paths]
    return file_paths

rule annotate_celltypes_nwconsensus:
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: 64000 * attempt
    conda:
        "envs/txsim-env.yaml"
    input:
        lambda w: input_files_for_consensus_annotation(w.id_code, w.results, w.dataset, w.rep_id, w.method, "nwconsensus")
    output:
        '{results}/{dataset}/replicate{rep_id}/celltypes_{method}_nwconsensus-{id_code}.csv'
    params:
        hyper_params = lambda w: get_params('nwconsensus', int(w.id_code), 'hyper_params'),
        group_params = lambda w: get_params('nwconsensus', int(w.id_code), 'group_params')
    shell:
        "python3 scripts/annotate_celltypes_consensus_NWCS.py "
        "-i {input} "
        "-o {output} "
        
rule annotate_celltypes_gmconsensus:
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: 64000 * attempt
    conda:
        "envs/txsim-env.yaml"
    input:
        lambda w: input_files_for_consensus_annotation(w.id_code, w.results, w.dataset, w.rep_id, w.method, "gmconsensus")
    output:
        '{results}/{dataset}/replicate{rep_id}/celltypes_{method}_gmconsensus-{id_code}.csv'
    params:
        hyper_params = lambda w: get_params('gmconsensus', int(w.id_code), 'hyper_params'),
        group_params = lambda w: get_params('gmconsensus', int(w.id_code), 'group_params')
    shell:
        "python3 scripts/annotate_celltypes_consensus_GMCS.py "
        "-i {input} "
        "-o {output} "


#################################
# Prepare single cell reference #
#################################

rule normalize_sc:
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: 64000 * attempt
    conda:
        "envs/txsim-env.yaml"
    input:
        ref = lambda w: parsed.get_data_file(w.dataset, 'sc_data'),
        mol = lambda w: parsed.get_data_file_list(w.dataset, 'molecules')
    output:
        '{results}/{dataset}/sc_normalized.h5ad'
    shell:
        "python3 scripts/normalize_sc.py "
        "-sc {input.ref} "
        "-m {input.mol} "
        "-o {output} "


###########
# Metrics #
###########

rule metric:
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: 64000 * attempt
    conda:
        "envs/txsim-env.yaml"
    input:
        '{results}/{dataset}/{replicate}/counts_{methods}.h5ad',
        scd = '{results}/{dataset}/sc_normalized.h5ad'
    output:
        '{results}/{dataset}/{replicate}/metrics_{methods}.csv',
        '{results}/{dataset}/{replicate}/unfiltered/metrics_{methods}.csv'
    shell:
        "python3 scripts/calc_metrics.py "
        "-m {wildcards.methods} "
        "-d {wildcards.results}/{wildcards.dataset}/{wildcards.replicate} "
        "-sc {input.scd} "

rule quality_metric:
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: 64000 * attempt
    conda:
        "envs/txsim-env.yaml"
    input:
        '{results}/{dataset}/{replicate}/counts_{methods}.h5ad',
    output:
        '{results}/{dataset}/{replicate}/quality_metrics_{methods}.csv',
        '{results}/{dataset}/{replicate}/unfiltered/quality_metrics_{methods}.csv'
    shell:
        "python3 scripts/calc_quality_metrics.py "
        "-m {wildcards.methods} "
        "-d {wildcards.results}/{wildcards.dataset}/{wildcards.replicate} "

rule group_metrics:
    conda:
        "envs/txsim-env.yaml"
    input:
        parsed.get_metric_inputs, #function that can accept wildcards as input argument
        '{results}/{dataset}/group_metric_chunks.csv'
    output:
        '{results}/{dataset}/{replicate}/group_metrics-{id_code}.csv'
    shell:
        "python3 scripts/calc_group_metrics.py "
        "-d {wildcards.results}/{wildcards.dataset}/{wildcards.replicate} "
        "-id {wildcards.id_code}"


##############################
# Aggregate multiple samples #
##############################

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


###############################
# Metrics for aggregated data #
###############################

rule aggregate_metrics:
    conda:
        "envs/txsim-env.yaml"
    input:
        lambda w : parsed.get_aggregate_inputs(w, 'metrics') 
    output:
        '{results}/{dataset}/aggregated/aggregated_metrics_{method}.csv',
        '{results}/{dataset}/aggregated/unfiltered/aggregated_metrics_{method}.csv'

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
        '{results}/{dataset}/aggregated/aggregated_quality_metrics_{method}.csv',
        '{results}/{dataset}/aggregated/unfiltered/aggregated_quality_metrics_{method}.csv'
    shell:
        "python3 scripts/aggregate_metrics.py "
        "-m {wildcards.method} "
        "-d {wildcards.results}/{wildcards.dataset} "
        "-t quality_metrics "

rule aggregate_group_metrics:
    conda:
        "envs/txsim-env.yaml"
    input:
        lambda w : parsed.get_aggregate_inputs(w, 'group_metrics'),
        '{results}/{dataset}/aggregated/group_metrics-{id_code}.csv'
    output:
        '{results}/{dataset}/aggregated/aggregated_group_metrics-{id_code}.csv'
    shell:
        "python3 scripts/aggregate_group_metrics.py "
        "-d {wildcards.results}/{wildcards.dataset} "
        "-id {wildcards.id_code}"


