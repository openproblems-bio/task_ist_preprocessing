import yaml
import itertools
import json
import csv
import pandas as pd
from snakemake.io import expand
import os

class ParsedConfig:
    PROCESSES =  ['segmentation', 'assignment', 'normalize']

    def __init__(self, config):
        cfg = yaml.load(open(config, 'r'), Loader=yaml.FullLoader)
        self.cfg = cfg
        self.final_files = []
        #Maps all methods to list of parameters
        self.method_dict ={}

        #Creating all parameter combinations per batch
        for batch in self.cfg['PREPROCESSING']:
            #Run through each batch
            batch_combos = {}
            for group in ParsedConfig.PROCESSES:
                #Within each group of processes per batch, 
                #Generate batch parameter combinations for each method
                #Check if any of those combinations exist already
                #If not, add new entry to dictionary
                batch_combos[group] = []
                for method in self.cfg['PREPROCESSING'][batch][group]:
                    #List of all method parameter combinations
                    method_combos = []

                    group_params = self.cfg['PREPROCESSING'][batch].get(f'{group}_params')
                    params = self.cfg['PREPROCESSING'][batch][group][method]

                    #If there are params specific to a method, include them
                    if(params is not None):
                        #Create cartesian product of method-specific params
                        p_dict = list( dict(zip(params.keys(), x)) for x in itertools.product(*params.values()) )
                        
                        #Merge with cartesian product of overall parameters if necessary
                        if(group_params is not None):
                            temp_dict = { **{'p':p_dict} , **group_params}
                            combos = list(dict(zip(temp_dict.keys(), x)) for x in itertools.product(*temp_dict.values()))
                        else:
                            temp_dict = {'p':p_dict}
                            combos = list(dict(zip(temp_dict.keys(), x)) for x in itertools.product(*temp_dict.values()))
                        method_combos.extend(combos)
                        
                    #If no method params, only use overall params
                    else:
                        if(group_params is not None):
                            combos = list(dict(zip(group_params.keys(), x)) for x in itertools.product(*group_params.values()))
                            method_combos.extend(combos)
                        #If no method params and no overall params
                        else:
                            method_combos.extend({None})
                
                    #Check if any params match method_dict
                    #If not, new name->method_dict, add to batch_combos[group]
                    #Otherwise add method_dict name to batch_combos
                    if(self.method_dict.get(method) is None):
                        self.method_dict[method] = []
                    for new_params in method_combos:
                        unique = True
                        for i in range(len(self.method_dict[method])):
                            if self.method_dict[method][i] == new_params:
                                unique = False
                                idx = i
                        if unique:
                            idx = len(self.method_dict[method])
                            self.method_dict[method].append(new_params)
                        batch_combos[group].append(f'{method}-{idx}')

            #Generate the final files for the batch
            #Parameters such as 'seg' in format: '{method}-{id}'
            self.final_files.extend(
                expand(
                    cfg['ROOT']+'/{results}/{dataset}/metrics_{seg}_{assign}_{norm}.txt',
                    results=cfg['RESULTS'],
                    batch=batch,
                    dataset=cfg['PREPROCESSING'][batch]['dataset'],
                    seg=list(batch_combos['segmentation']),
                    assign=list(batch_combos['assignment']),
                    norm=list(batch_combos['normalize'])
                )
            )
        self.final_files = list(set(self.final_files))

        #Save the dictionary of method-id -> parameter combination
        output = {"names":[], "id":[], "params":[]}
        for method in self.method_dict:
            output["names"] += [method]  * len(self.method_dict[method])
            output["id"].extend(range(len(self.method_dict[method])))
            output["params"].extend(self.method_dict[method])
        
        df = pd.DataFrame.from_dict(output)
        output_folder = self.cfg['ROOT']+'/' + self.cfg['RESULTS']
        if not os.path.exists(output_folder):
            os.mkdir(output_folder)
        df.to_csv(output_folder + '/params_dict.csv')
        
    def gen_file_names(self):
        return self.final_files
    
    #TODO Type hints
    def get_data_file(self, dataset, file_name):
        return self.cfg['ROOT'] + '/' + self.cfg['DATA_SCENARIOS'][dataset][file_name]
    
    def get_method_params(self, method, id):
        return self.method_dict[method][id]
