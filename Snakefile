import yaml
import itertools
import json
import csv
import pandas as pd
from TxsimConfig import *

parsed = ParsedConfig('configs/config.yaml')
final_files = parsed.gen_file_names()

#Main rule (will always be run)
rule all:
    input:
        final_files

rule label:
    output:
        '{results}/{dataset}/label_{seg}-{shp}.tif'
    params:
        img = lambda w: parsed.get_data_file(w.dataset, 'image'),
        exp = lambda w: parsed.get_method_params(w.seg, int(w.shp))['expand'],
        bry = lambda w: "-b " if parsed.get_method_params(w.seg, int(w.shp))['binary'] else ""
    shell:
        "python load_images.py "
        "-i {params.img} "
        "-o {wildcards.results}/{wildcards.dataset} "
        "-e {params.exp} "
        "-s {wildcards.seg} "
        "{params.bry}"

rule assign:
    input:
        '{results}/{dataset}/label_{seg}-{shp}.tif'
    params:
        mol = lambda w: parsed.get_data_file(w.dataset, 'molecules'),
        scd = lambda w: parsed.get_data_file(w.dataset, 'sc_data'),
        hyp = lambda w: parsed.get_method_params(w.assign, int(w.ahp))['p']
    output:
        '{results}/{dataset}/assignments_{seg}-{shp}_{assign}-{ahp}.csv'
    shell:
        "python run_{wildcards.assign}.py "
        "-m {params.mol} "
        "-p \"{params.hyp}\" "
        "-sc {params.scd} "
        "-d {wildcards.results}/{wildcards.dataset} "
        "-s {wildcards.seg} "

rule counts:
    input:
        '{results}/{dataset}/assignments_{seg}-{shp}_{assign}-{ahp}.csv'
    output:
        '{results}/{dataset}/counts_{seg}-{shp}_{assign}-{ahp}_{norm}-{nhp}.h5ad'
    shell:
        "python gen_counts.py "
        "-a {wildcards.assign} "
        "-d {wildcards.results}/{wildcards.dataset} "
        "-s {wildcards.seg} "
        "-n {wildcards.norm} "

rule metric:
    input:
        '{results}/{dataset}/counts_{seg}_{assign}_{norm}-{nhp}.h5ad'
    params:
        scd = lambda w: parsed.get_data_file(w.dataset, 'sc_data'),
    output:
        '{results}/{dataset}/metrics_{seg}_{assign}_{norm}-{nhp}.txt'
    shell:
        "python calc_metrics.py "
        "-a {wildcards.assign} "
        "-d {wildcards.results}/{wildcards.dataset} "
        "-s {wildcards.seg} "
        "-n {wildcards.norm} "
        "-sc {params.scd} "
