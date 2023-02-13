
import itertools
import json
import tifffile
import csv
import pandas as pd
import anndata as ad
import numpy as np
from snakemake.io import expand
import os
import yaml
from collections import OrderedDict

class ParsedConfig:

    def __init__(
        self, 
        config: dict,
        defaults: str = None
    ):
        self.cfg = config
        self.final_files = []
        #Maps all methods to list of parameters
        self.method_dict ={}

        #Create (and read previous) params dict
        self.gen_params_dict()

        #Creating all parameter combinations per batch
        self.gen_combinations()

        #Save parameters
        self.save_parameters(defaults)

        #Calculate metrics across combinations
        self.gen_metrics()

        #Checking data files
        self.check_files()
    
    def gen_params_dict(self):
        # Create params output folder
        self.output_folder = os.path.join(self.cfg['RESULTS'], 'params')
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)

        # Read from previous file if it exists already
        elif os.path.exists( os.path.join(self.output_folder, 'params_dict.csv')  ):
            methods = pd.read_csv(os.path.join(self.output_folder, 'params_dict.csv'), index_col=0)   
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

    def gen_combinations(self):
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

            file_template = os.path.join(self.cfg['RESULTS'], "{dataset}/replicate{rep_id}/metrics") #dont use 'f' here
            wildcards = {"dataset": ['']}
            new_final_files = []

            for group in self.cfg['PREPROCESSING'][batch]['workflow']:
                file_template = file_template + "_{" + group + "}"
                wildcards[group] = list(batch_combos[group])


            for d_set in self.cfg['PREPROCESSING'][batch]['dataset']:
                #Number of replicates
                wildcards["rep_id"] = range(1, len( self.cfg['DATA_SCENARIOS'][d_set]['images'])+1)
                wildcards["dataset"]= [d_set]
                new_final_files.extend( expand((file_template + ".csv"),**wildcards))
                new_final_files.extend( expand((file_template.replace("replicate{rep_id}/", "aggregated/aggregated_") + ".csv"),**wildcards))
                new_final_files.extend( expand((file_template.replace("replicate{rep_id}/", "aggregated/") + ".csv"),**wildcards))
                

            # also for quality metrics and counts
            
            quality_metrics = new_final_files.copy()
            for final_file in quality_metrics:
                if '/metrics' in final_file:
                    new_final_files.append(final_file.replace('/metrics', '/quality_metrics'))
                if '_metrics' in final_file:
                    new_final_files.append(final_file.replace('_metrics', '_quality_metrics'))
                if 'aggregated_metrics' in final_file:
                    new_final_files.append(final_file.replace('aggregated_metrics', 'counts').replace(".csv",".h5ad"))
                    new_final_files.append(final_file.replace('aggregated_quality_metrics', 'quality_metrics'))
                    

            self.final_files.extend(new_final_files)

        self.final_files = list(set(self.final_files))
        #for f in self.final_files: print(f)

    def save_parameters(self, defaults: str):
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
            m_df.to_csv(self.output_folder + f'/{method}_params.csv')
        
        df = pd.DataFrame.from_dict(output)
        df.to_csv(self.output_folder + '/params_dict.csv')

        # Create a readable dictionary for every combination in final_files
        readable_df = pd.DataFrame()
        method_names = []
        for f in self.final_files:
            if "quality" in f: continue
            if "counts" in f: continue
            if "aggregated" in f: continue
            name = f.split('metrics_')[1].replace('.csv','')
            method_list = name.split('_')
            for m in method_list:
                method,idx = m.split('-'); idx = int(idx)
                d = self.method_dict[method][idx]
                if defaults is not None:
                    with open(defaults,'r') as def_file:
                        def_dict = yaml.safe_load(def_file).get(method)
                        if d is None:
                            d = def_dict
                        elif def_dict is not None:
                            def_dict.update(d)
                            d = def_dict
                if d is None: #This might not be necessary 
                    continue
                for key, value in d.items():
                    if key == 'p':
                        for k, v in d['p'].items():
                            if not k in readable_df:
                                readable_df.insert(len(readable_df.columns), k, pd.Series( [v], index =[ [name] ]))
                                method_names.append(method)
                            else:
                                readable_df.loc[name, k] = v
                    else:
                        if not key in readable_df:
                            readable_df.insert(len(readable_df.columns), key, pd.Series( [value], index =[ [name] ]))
                            method_names.append(method)
                        else:
                            readable_df.loc[name, key] = value
            #pd.concat(readable_df, pd.Series([None]*len(readable_df.columns), index = [[name]]))
            if name not in readable_df.index:
                blank = pd.Series([np.nan]*len(readable_df.columns), index = readable_df.columns ).to_frame().T
                blank.index = [name]
        
        #Add in the name of the method for each parameter. Note sure if there is a better way to do this
        readable_df.columns = pd.MultiIndex.from_tuples(list(zip(*np.array([readable_df.columns.to_numpy(), method_names]))) )            
        
        readable_df.to_csv(self.output_folder + '/params_dict_readable.csv')

    def gen_metrics(self):
        #TODO FIX HOW THE THINGY WORKS
        if not bool(self.cfg['METRICS']): return #Trick to see if its empty
        #Create dictionary of all runs for each dataset
        names = {}
        for f in self.final_files:
            if "/metrics" not in f : continue
            dataset = f.split('/metrics_')[0].split('/')[-2] #-2 since second to last = dataset, last = replicate
            name = f.split('metrics_')[1].replace('.csv','')
            if names.get(dataset) is None:
                names[dataset] = [name]
            else:
                names[dataset].append(name)
                
        self.metric_input_files = {}

        for metric_batch in self.cfg['METRICS']:
            dataset = self.cfg['METRICS'][metric_batch]['dataset']
            for rep in range(1, len( self.cfg['DATA_SCENARIOS'][dataset]['images'])+1):
                self.final_files.append( 
                        os.path.join(self.cfg['RESULTS'], f"{dataset}/replicate{rep}/group_metrics.csv"))
            self.final_files.append( os.path.join(self.cfg['RESULTS'], f"{dataset}/aggregated/group_metrics.csv"))
            self.final_files.append( os.path.join(self.cfg['RESULTS'], f"{dataset}/aggregated/aggregated_group_metrics.csv"))
            required_inputs = []
            for run_name in names[dataset]:
                required_inputs.append(
                    os.path.join(self.cfg['RESULTS'], f"{dataset}/metrics_{run_name}.csv"))
            self.metric_input_files[dataset] = required_inputs

    def check_files(self):
        for dataset in self.cfg['DATA_SCENARIOS']:
            num_replicates = 0
            for key in self.cfg['DATA_SCENARIOS'][dataset]:
                if 'images' in key or 'molecules' in key:
                    num_reps = len(self.cfg['DATA_SCENARIOS'][dataset][key])
                    if num_replicates > 0 and num_replicates != num_reps:
                        print(f"Warning: number of replicates in dataset '{dataset}' does not match ({num_reps} and {num_replicates})")
                    num_replicates = num_reps
                    for rep in range(0, ):
                        if not os.path.exists(self.get_replicate_file(dataset,key,rep)):
                            raise Exception(f"The following replicate was not found: {self.get_replicate_file(dataset,key,rep)}")
                elif key != 'root_folder':
                    if not os.path.exists(self.get_data_file(dataset,key)):
                        raise Exception(f"The following file was not found: {self.get_data_file(dataset,key)}")
            self.check_dataset(dataset)

    def gen_file_names(self):
        return self.final_files

    def get_metric_inputs(self, wildcards):
        output = []
        for input_file in self.metric_input_files[wildcards.dataset]:
            output.append(input_file.replace("/metrics", f"/{wildcards.replicate}/metrics"))
        return list(set(output))

    #`dataset` should be name of dataset, `file_name`` should be desired file, both as str
    def get_data_file(self, dataset, file_name):
        return os.path.join(self.cfg['DATA_SCENARIOS'][dataset]['root_folder'] , self.cfg['DATA_SCENARIOS'][dataset][file_name])

    def get_replicate_file(self, dataset, file_name, replicate_number):
        return os.path.join(self.cfg['DATA_SCENARIOS'][dataset]['root_folder'] , self.cfg['DATA_SCENARIOS'][dataset][file_name][replicate_number])
    
    def get_aggregate_inputs(self, wildcards, file_type):
        required_inputs = []
        
        #Prevent confusion between quality_metrics and metrics
        if file_type == 'metrics': file_type = "/metrics"

        for final_file in self.final_files:
            #Look for files with the same dataset and method signature as aggregated file
            if "aggregated/aggregated" in final_file:
                continue
            if f"{wildcards.results}/{wildcards.dataset}" not in final_file:
                continue
            if file_type == 'group_metrics':
                if 'group_metrics' in final_file and 'aggregated' not in final_file:
                    required_inputs.append(final_file)
                continue
            if f"{wildcards.method}" not in final_file:
                continue
            #Look for files with right type (counts, metrics, quality metrics)
            if file_type == 'counts' and '/metrics' in final_file:
                if 'aggregated' in final_file: continue # prevent cyclic graph
                required_inputs.append(final_file.replace('/metrics', '/counts').replace('.csv','.h5ad'))
            elif file_type != 'counts' and file_type in final_file:
                required_inputs.append(final_file)
        
        return list(set(required_inputs))

    #`method` should be name of method, `id_code` should be an int
    def get_method_params(self, method, id_code):
        if self.method_dict.get(method) is None:
            return None
        if id_code >= len(self.method_dict[method]):
            return None
        return self.method_dict[method][id_code]
    
    def check_dataset(self, dataset:str):
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
        spots_file = self.get_replicate_file(dataset,'molecules',0) #Just check first one in replicates
        dapi_file = self.get_replicate_file(dataset,'images',0)
        
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
