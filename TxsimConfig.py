
import itertools
import json
import csv
import pandas as pd
from snakemake.io import expand
import os

class ParsedConfig:

    def __init__(
        self, 
        config: dict
    ):
        self.cfg = config
        self.final_files = []
        #Maps all methods to list of parameters
        self.method_dict ={}

        #Creating all parameter combinations per batch
        for batch in self.cfg['PREPROCESSING']:
            #Run through each batch
            batch_combos = {}
            #TODO decide if workflow should be in the config or generated from the config
            # would have to use regex if latter? 
            for group in self.cfg['PREPROCESSING'][batch]['workflow']:
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
            #TODO make the format more flexible 
            #Possible solution: 
            # create a string as the template file name based on contents of workflow
            # pass in a dictionary of group -> list(batch_combos) to expand
            # basically just **kwargs
            # ie replace seg= with segmentation
            

            file_template = os.path.join(self.cfg['RESULTS'], "{dataset}/metrics")
            wildcards = {'dataset': self.cfg['PREPROCESSING'][batch]['dataset']}

            for group in self.cfg['PREPROCESSING'][batch]['workflow']:
                file_template = file_template + "_{" + group + "}"
                wildcards[group] = list(batch_combos[group])

            self.final_files.extend( expand((file_template + ".txt"),**wildcards))
            
        self.final_files = list(set(self.final_files))

        #Save the dictionary of method-id -> parameter combination
        output = {"names":[], "id":[], "params":[]}
        for method in self.method_dict:
            output["names"] += [method]  * len(self.method_dict[method])
            output["id"].extend(range(len(self.method_dict[method])))
            output["params"].extend(self.method_dict[method])
        
        df = pd.DataFrame.from_dict(output)
        output_folder = self.cfg['RESULTS']
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        df.to_csv(output_folder + '/params_dict.csv')
        
    def gen_file_names(self):
        return self.final_files
    
    #TODO Type hints
    def get_data_file(self, dataset, file_name):
        return os.path.join(self.cfg['DATA_SCENARIOS'][dataset]['root_folder'] , self.cfg['DATA_SCENARIOS'][dataset][file_name])
    
    def get_method_params(self, method, id):
        if self.method_dict.get(method) is None:
            return None
        if id >= len(self.method_dict[method]):
            return None
        return self.method_dict[method][id]
