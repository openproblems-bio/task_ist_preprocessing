
import itertools
import json
import tifffile
import csv
import pandas as pd
import anndata as ad
import numpy as np
from snakemake.io import expand
import os
from collections import OrderedDict

class ParsedConfig:

    def __init__(
        self, 
        config: dict
    ):
        self.cfg = config
        self.final_files = []
        #Maps all methods to list of parameters
        self.method_dict ={}

        # Create params output folder
        output_folder = os.path.join(self.cfg['RESULTS'], 'params')
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        # Read from previous file if it exsists already
        elif os.path.exists( os.path.join(output_folder, 'params_dict.csv')  ):
            methods = pd.read_csv(os.path.join(output_folder, 'params_dict.csv'), index_col=0)   
            for method_name in pd.unique(methods['names']):
                for param_str in methods[methods['names'] == method_name]['params'].values:
                    try:
                        params = eval(param_str)
                        if(self.method_dict.get(method_name) is None):
                            self.method_dict[method_name] = [params]
                        else:
                            self.method_dict[method_name].append(params)
                    except:
                        if(self.method_dict.get(method_name) is None):
                            self.method_dict[method_name] = [None]
                        else:
                            self.method_dict[method_name].append(None)

        #Creating all parameter combinations per batch
        for batch in self.cfg['PREPROCESSING']:
            #Run through each batch
            batch_combos = {}
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

            file_template = os.path.join(self.cfg['RESULTS'], "{dataset}/metrics")
            wildcards = {'dataset': self.cfg['PREPROCESSING'][batch]['dataset']}

            for group in self.cfg['PREPROCESSING'][batch]['workflow']:
                file_template = file_template + "_{" + group + "}"
                wildcards[group] = list(batch_combos[group])

            self.final_files.extend( expand((file_template + ".csv"),**wildcards))
            
            # also for quality metrics
            
            file_template = os.path.join(self.cfg['RESULTS'], "{dataset}/quality_metrics")
            wildcards = {'dataset': self.cfg['PREPROCESSING'][batch]['dataset']}

            for group in self.cfg['PREPROCESSING'][batch]['workflow']:
                file_template = file_template + "_{" + group + "}"
                wildcards[group] = list(batch_combos[group])

            self.final_files.extend( expand((file_template + ".csv"),**wildcards))            
            
        self.final_files = list(set(self.final_files))

        #Save parameters

        #Save the dictionary of method-id -> parameter combination
        output = {"names":[], "id":[], "params":[]}
        for method in self.method_dict:
            output["names"] += [method]  * len(self.method_dict[method])
            output["id"].extend(range(len(self.method_dict[method])))
            output["params"].extend(self.method_dict[method])

            #Save each method to csv file
            m_df = pd.DataFrame()
            i = 0
            for d in self.method_dict[method]:
                if d is None:
                    i+=1
                    continue
                for key, value in d.items():
                    if key == 'p':
                        for k, v in d['p'].items():
                            if not k in m_df:
                                m_df.insert(len(m_df.columns), k, pd.Series( [v], index =[i]))
                            else:
                                m_df.loc[i,k] = v
                    else:
                        if not key in m_df:
                            m_df.insert(len(m_df.columns), key, pd.Series( [value], index =[i]))
                        else:
                            m_df.loc[i,key] = value
                i+=1
            m_df.to_csv(output_folder + f'/{method}_params.csv')
        
        df = pd.DataFrame.from_dict(output)
        df.to_csv(output_folder + '/params_dict.csv')

        readable_df = pd.DataFrame()
        for f in self.final_files:
            if "quality" in f: continue
            name = f.split('metrics_')[1].replace('.csv','')
            method_list = name.split('_')
            for m in method_list:
                method,idx = m.split('-'); idx = int(idx)
                d = self.method_dict[method][idx]
                if d is None:
                    continue
                for key, value in d.items():
                    if key == 'p':
                        for k, v in d['p'].items():
                            if not k in readable_df:
                                readable_df.insert(len(readable_df.columns), k, pd.Series( [v], index =[[name] ]))
                            else:
                                
                                readable_df.loc[ name ,k] = v
                    else:
                        if not key in readable_df:
                            readable_df.insert(len(readable_df.columns), key, pd.Series( [value], index =[ [name] ]))
                        else:
                            readable_df.loc[  name ,key] = value
            #pd.concat(readable_df, pd.Series([None]*len(readable_df.columns), index = [[name]]))
            if name not in readable_df.index:
                blank = pd.Series([np.nan]*len(readable_df.columns), index = readable_df.columns ).to_frame().T
                blank.index = [name]
                    
        
        readable_df.to_csv(output_folder + '/params_dict_readable.csv')

        #Checking data files
        for dataset in self.cfg['DATA_SCENARIOS']:
            for key in self.cfg['DATA_SCENARIOS'][dataset]:
                if key != 'root_folder':
                    if not os.path.exists(self.get_data_file(dataset,key)):
                        raise Exception(f"The following file was not found: {self.get_data_file(dataset,key)}")
            self.check_input_files(dataset)
        
    def gen_file_names(self):
        return self.final_files
    
    #`dataset` should be name of dataset, `file_name`` should be desired file, both as str
    def get_data_file(self, dataset, file_name):
        return os.path.join(self.cfg['DATA_SCENARIOS'][dataset]['root_folder'] , self.cfg['DATA_SCENARIOS'][dataset][file_name])
    
    #`method` should be name of method, `id_code` should be an int
    def get_method_params(self, method, id_code):
        if self.method_dict.get(method) is None:
            return None
        if id_code >= len(self.method_dict[method]):
            return None
        return self.method_dict[method][id_code]
    
    
    def check_input_files(self, dataset:str):
        """Test multiple requirements for the input files of a given dataset
        
        Arguments
        ---------
        dataset: str
            Name of dataset as given in the config.yaml
        
        """
        
        #############
        # Load data #
        #############
        
        sc_file = self.get_data_file(dataset,'sc_data')
        spots_file = self.get_data_file(dataset,'molecules')
        dapi_file = self.get_data_file(dataset,'image')
        
        adata = ad.read(sc_file)
        spots = pd.read_csv(spots_file)
        dapi = tifffile.imread(dapi_file)
        
        #########
        # Tests #
        #########
        
        # Test that there are no genes in the spatial data that we don't find in the sc data
        not_included = set(spots[spots.columns[0]]) - set(adata.var_names)
        if len(not_included) > 0:
            raise Exception(f"For dataset '{dataset}' the following genes are in the spatial data but not RNAseq: {list(not_included)}")
            
        # Test that the spots.csv header is correct
        for key in ["Gene", "x", "y"]:
            if key not in spots.columns:
                raise Exception(f"For dataset '{dataset}': '{key}' is not found in the header ({spots.columns}) in file {spots_file}")
        
        # Test that coordinates in spots.csv align with tif image pixels
        if spots[["x","y"]].min().min() < 0:
            raise Exception(f"For dataset '{dataset}': found negative coordinate values in {spots_file}")
        if (spots["y"].max() > dapi.shape[0]) or (spots["x"].max() > dapi.shape[1]):
            raise Exception(f"For dataset '{dataset}': Coordinates in \n\t{spots_file}\nexceed pixel range in \n\t{dapi_file}")
        border_th = 0.2
        large_border = spots["y"].min() > (border_th * dapi.shape[0])
        large_border |= spots["y"].max() < ( (1 - border_th) * dapi.shape[0])
        large_border |= spots["x"].min() > (border_th * dapi.shape[1])
        large_border |= spots["x"].max() < ( (1 - border_th) * dapi.shape[1])
        if large_border:
            raise Exception(
                f"For dataset '{dataset}': Either coordinates in \n\t{spots_file} \nare not aligned with pixel units in \n\t{dapi_file}\n" + 
                f"or there is a very large border in \n\t{dapi_file}\nwithout gene signals in \n\t{spots_file}"
            )
