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
from typing import Optional

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
        self.gen_group_metrics(defaults)

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

                    # If there are params specific to a method, include them     
                    if params is not None:
                        #Create cartesian product of method-specific params
                        p_dict = list( dict(zip(params.keys(), x)) for x in itertools.product(*params.values()) )
                    else: p_dict = [{}]
                    
                    # If there are params for the group, include them     
                    if group_params is not None:
                        #Create cartesian product of group params
                        group_dict = list( dict(zip(group_params.keys(), x)) for x in itertools.product(*group_params.values()) )
                    else: group_dict = [{}]
                    
                    # Merge with cartesian product of overall parameters
                    temp_dict = { **{'hyper_params':p_dict} , **{'group_params':group_dict}}
                    combos = list(dict(zip(temp_dict.keys(), x)) for x in itertools.product(*temp_dict.values()))
                    method_combos.extend(combos)
        

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
                    if key == 'hyper_params':
                        for k, v in d['hyper_params'].items():
                            if not k in m_df:
                                m_df.insert(len(m_df.columns), k, pd.Series( [v], index =[i]))
                            else:
                                m_df.loc[i,k] = v
                    
                    elif key == 'group_params':
                        for k, v in d['group_params'].items():
                            if not k in m_df:
                                m_df.insert(len(m_df.columns), k, pd.Series( [v], index =[i]))
                            else:
                                m_df.loc[i,k] = v
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
                #Check each method based on the name of the final file
                method,idx = m.split('-'); idx = int(idx)
                if self.method_dict[method][idx] is not None:
                    d = self.method_dict[method][idx].copy()
                else:
                    d = None
                #Compare the config parameters with the default parameters
                with open(defaults,'r') as def_file:
                    def_dict = yaml.safe_load(def_file).get(method)
                    if def_dict is None:
                        raise Exception(f"No defaults found for method {method}")
                    if d is None:
                        d = def_dict
                    else:
                        if method == 'pciSeq' and d['hyper_params'].get('opts') is not None:
                            d['hyper_params'] = d['hyper_params']['opts']
                        #Check to see if unknown parameters
                        if d.get('hyper_params') is not None:
                            for key in d['hyper_params']:
                                if key not in def_dict:
                                    raise Exception(f"Unknown parameter '{key}' for method {method}")
                            #Update dictionary with defaults and save back into original
                            def_dict.update(d['hyper_params'])
                            d.pop('hyper_params')
                        #TODO check for the general parameters like segmentation_parameters
                        def_dict.update(d)
                        d = def_dict
                #Run through each key and insert it into a dataframe
                for key, value in d.items():
                    if key == 'group_params':
                        for k, v in d['group_params'].items():
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

            if name not in readable_df.index:
                blank = pd.Series([np.nan]*len(readable_df.columns), index = readable_df.columns ).to_frame().T
                blank.index = [name]
        
        #Add in the name of the method for each parameter. Not sure if there is a better way to do this
        readable_df.columns = pd.MultiIndex.from_tuples(list(zip(*np.array([readable_df.columns.to_numpy(), method_names]))) )            
        
        readable_df.to_csv(self.output_folder + '/params_dict_readable.csv') 

    def gen_group_metrics(self, defaults: str):
        #TODO fix when metrics are called on dataset not being run (is this even a bug?)
        if not bool(self.cfg['METRICS']): return #Trick to see if its empty
                
        self.metric_input_files = {}

        # pairs of runs per dataset that will be used for group metrics
        # 1 run = combination of methods (eg watershed-1_pciSeq-3_area-0)
        # for the pairs, type = frozenset
        self.list_of_group_pairs = {}
        chunk_size = {}

        for metric_batch in self.cfg['METRICS']:
            dataset = self.cfg['METRICS'][metric_batch]['dataset']
            
            names = [] #temp list of names in a given batch
            # 2 conditions: all or specific files, in either case, generate single list of all run names
            if (self.cfg['METRICS'][metric_batch]['files'] == 'all'):
                # go through every file for this dataset, and get all the runs 
                # from the regular metric file, grab the counts file
                for f in self.final_files:
                    #check if file is metric
                    if "/metrics" not in f : continue
                    #-2 since second to last = dataset, last = replicate
                    if (dataset != f.split('/metrics_')[0].split('/')[-2]): continue 
                    #convert to counts
                    name = f.split('metrics_')[1].replace('.csv','')
                    names.append(name)
            else:
                #TODO modify this check if input format is correct
                names.extend(self.cfg['METRICS'][metric_batch]['files'])

            # generate all combinations of names list, append the combinations to growing frozenset list
            # use frozenset so order doesn't matter AND it is hashable 
            # so we can have set of frozensets -> remove duplicate frozensets

            #create blank list if dataset doesn't have list yet
            if (self.list_of_group_pairs.get(dataset) is None):
                self.list_of_group_pairs[dataset] = []
            
            #gen all combinations
            names = list(set(names))
            for pair in itertools.combinations(names,2):
                self.list_of_group_pairs[dataset].append(frozenset(pair))
                
            # remove duplicate pairs
            self.list_of_group_pairs[dataset] = list(set(self.list_of_group_pairs[dataset]))
            
            #TODO make sure all chunk_size is same for given dataset
            # read a chunk size from metrics, use default otherwise 
            chunk_size[dataset] = self.cfg['METRICS'][metric_batch].get('chunk_size')
            if (chunk_size[dataset] is None): 
                with open(defaults,'r') as def_file:
                    chunk_size[dataset] = yaml.safe_load(def_file).get('metrics').get('chunk_size')

        for dataset in self.list_of_group_pairs.keys():
            #split the list into the right number of chunks -> create new list of numbers [0, 0, 0, 0, 0, 1, 1, 1...etc]
            len_pairs = len(self.list_of_group_pairs[dataset])
            chunk_id_list = np.int_(np.floor(np.array(range(0,len_pairs)) / chunk_size[dataset]))

            # create dataframe with 3 columns, run1, run2, and chunk_id
            tuple_list = [tuple(x) for x in self.list_of_group_pairs[dataset]]
            table_run_pairs  = pd.DataFrame(tuple_list, columns = ['run1','run2']).sort_values(by=['run1','run2'])
            table_run_pairs = table_run_pairs.reset_index(drop=True)
            table_run_pairs['chunk_id'] = chunk_id_list

            #TODO add reading matrix as input
            all_runs = list(np.unique(table_run_pairs[['run1','run2']].values))
            run_matrix = self.gen_group_metric_matrix(all_runs, dataset)

            #TODO check if existing is same, otherwise, write new table (unecessary?)
            update_table = True
            if os.path.exists( os.path.join(self.cfg['RESULTS'], f"{dataset}/group_metric_matrix.csv") ):
                old_matrix = pd.read_csv(os.path.join(self.cfg['RESULTS'], f"{dataset}/group_metric_matrix.csv"),index_col=0)
                old_matrix.sort_index(axis=0, inplace=True)
                old_matrix.sort_index(axis=1, inplace=True)
                run_matrix.sort_index(axis=0, inplace=True)
                run_matrix.sort_index(axis=1, inplace=True)
                if(old_matrix.equals(run_matrix)):
                    print("Not updating group metrics")
                    update_table = False

            if (update_table):
                print("Updating group metrics matrix and table")
                run_matrix.to_csv(os.path.join(self.cfg['RESULTS'], f"{dataset}/group_metric_matrix.csv"))
                table_run_pairs.to_csv(os.path.join(self.cfg['RESULTS'], f"{dataset}/group_metric_chunks.csv"), index=False)

            #Add final_file for each replicate, aggregate, and aggregated/aggregated
            for rep in range(1, len( self.cfg['DATA_SCENARIOS'][dataset]['images'])+1):
                for i in range(np.max(chunk_id_list)+1): 
                    self.final_files.append( 
                        os.path.join(self.cfg['RESULTS'], f"{dataset}/replicate{rep}/group_metrics-{i}.csv"))

            for i in range(np.max(chunk_id_list)+1): 
                self.final_files.append( os.path.join(self.cfg['RESULTS'], f"{dataset}/aggregated/group_metrics-{i}.csv"))

            for i in range(np.max(chunk_id_list)+1): 
                self.final_files.append( os.path.join(self.cfg['RESULTS'], f"{dataset}/aggregated/aggregated_group_metrics-{i}.csv"))
            
            # For each run in this dataset, make the counts an required input
            # Generic for dataset, add replicate in get_metric_inputs
            self.metric_input_files[dataset] = {}
            for i in range(np.max(chunk_id_list)+1):
                input_runs = list(table_run_pairs['run1'][table_run_pairs['chunk_id'] == i].values)
                input_runs.extend(table_run_pairs['run2'][table_run_pairs['chunk_id'] == i].values)
                input_runs = list(set(input_runs))
                required_inputs = []
                for run_name in input_runs:
                    required_inputs.append(
                        os.path.join(self.cfg['RESULTS'], f"{dataset}/counts_{run_name}.h5ad"))
                self.metric_input_files[dataset][str(i)] = required_inputs
        
    def gen_group_metric_matrix(self, all_runs: list, dataset: str):
        #Blank dataframe
        run_matrix = pd.DataFrame(np.zeros((len(all_runs), len(all_runs)),dtype=np.bool_), columns = all_runs, index = all_runs)

        #populate matrix based on pairings
        for p in self.list_of_group_pairs[dataset]:
            pair = tuple(p)
            run_matrix.loc[pair[0], pair[1]] = True
            run_matrix.loc[pair[1], pair[0]] = True

        # Set lower triangular to False (prevent duplicates) and save
        run_matrix.iloc[tuple(zip(np.tril_indices_from(run_matrix)))] = False
        return run_matrix

    def check_files(self):
        for dataset in self.cfg['DATA_SCENARIOS']:
            num_replicates = 0
            for key in self.cfg['DATA_SCENARIOS'][dataset]:
                
                #print("#x#x#x#x#x#x#x#x")
                #print(self.cfg['DATA_SCENARIOS'][dataset][key])
                #print(type(self.cfg['DATA_SCENARIOS'][dataset][key]))
                #print((type(self.cfg['DATA_SCENARIOS'][dataset][key]) == str))
                #print(isinstance(self.cfg['DATA_SCENARIOS'][dataset][key],str))
                
                #For replicate files
                if 'images' in key or 'molecules' in key:
                    num_reps = len(self.cfg['DATA_SCENARIOS'][dataset][key])
                    #Confirm number of replicates is equal for each key
                    if num_replicates > 0 and num_replicates != num_reps:
                        print(f"Warning: number of replicates in dataset '{dataset}' does not match ({num_reps} and {num_replicates})")
                    num_replicates = num_reps
                    for rep in range(0, num_replicates):
                        if not os.path.exists(self.get_replicate_file(dataset,key,rep)):
                            raise Exception(f"The following replicate was not found: {self.get_replicate_file(dataset,key,rep)}")
                
                #Check non-replicate files
                elif (key != 'root_folder') and (key != 'check_dataset') and (type(self.cfg['DATA_SCENARIOS'][dataset][key]) == str):
                    #print(f"key : {key}")
                    #print(f"key : {key}")
                    #print(f"type(self.cfg['DATA_SCENARIOS'][dataset][key]) == str : {type(self.cfg['DATA_SCENARIOS'][dataset][key]) == str}")
                    if not os.path.exists(self.get_data_file(dataset,key)):
                        raise Exception(f"The following file was not found: {self.get_data_file(dataset,key)}")
            
            #Confirm dataset needs to be checked (on by default)
            if 'check_dataset' in self.cfg['DATA_SCENARIOS'][dataset] and self.cfg['DATA_SCENARIOS'][dataset]['check_dataset'] is False:
                print(f"Not checking dataset '{dataset}' due to flag")
            else:
                #print(f"Checking dataset '{dataset}'...")
                self.check_dataset(dataset)

    def _file_names_for_tile_info(self):
        directories = [f.rsplit("/",1)[0] for f in self.final_files]
        directories = np.unique(directories)
        directories = [d for d in directories if ("replicate" in d.split("/")[-1])]
        return [d + "/tile_info.csv" for d in directories]

    def gen_file_names(self):
        if not self.cfg["PRE_COMPUTE_TILE_INFO"]:
            return self.final_files
        else:
            return self._file_names_for_tile_info()

    def get_metric_inputs(self, wildcards):
        #TODO modify so it adds file type to names
        output = []
        for input_file in self.metric_input_files[wildcards.dataset][wildcards.id_code]:
            output.append(input_file.replace("/counts", f"/{wildcards.replicate}/counts"))
        return list(set(output))

    #`dataset` should be name of dataset, `file_name`` should be desired file, both as str
    def get_data_file(self, dataset, file_name):
        #print("###############")
        #print(dataset, file_name)
        #print(self.cfg['DATA_SCENARIOS'][dataset]['root_folder'] , self.cfg['DATA_SCENARIOS'][dataset][file_name])
        return os.path.join(self.cfg['DATA_SCENARIOS'][dataset]['root_folder'] , self.cfg['DATA_SCENARIOS'][dataset][file_name])

    def get_data_file_list(self, dataset, file_list_name):
        assert type(self.cfg['DATA_SCENARIOS'][dataset][file_list_name]) == list;
        paths = [
            os.path.join(self.cfg['DATA_SCENARIOS'][dataset]['root_folder'],name) for name in self.cfg['DATA_SCENARIOS'][dataset][file_list_name]
        ]
        return paths

    def get_replicate_file(self, dataset, file_name, replicate_number):
        return os.path.join(self.cfg['DATA_SCENARIOS'][dataset]['root_folder'] , self.cfg['DATA_SCENARIOS'][dataset][file_name][replicate_number])
    
    def get_tiles_aggregation_input_files(
            self, dataset: str, rep_id: int, method: str, id_code: int, seg: Optional[str] = None, 
            areas: bool = True
        ):
        """Get input files of assignments (and areas) per tile for aggregation
        
        Arguments
        ---------
        dataset: str
        rep_id: int
        method: str
            Assignment method (baysor or clustermap)
        id_code: int
            Parameter id of assignment method
        seg: str
            Segmentation method with id code if prior segmentation is included
        areas: bool
            Whether to include csvs for cell areas of each tile
            
        Returns
        -------
        list:
            Assignments (and areas) csv paths
        
        """
        rep_dir = os.path.join(self.cfg['RESULTS'], f"{dataset}/replicate{rep_id}")
        tile_info_path = os.path.join(rep_dir,"tile_info.csv")
        if not os.path.exists(tile_info_path):
            raise ValueError("tile_info.csv doesn't exist. First run the pipeline once with 'PRE_COMPUTE_TILE_INFO: True' in the config.yaml")
        df = pd.read_csv(tile_info_path, index_col=0)
        nx = df.loc[method,"nx"]
        ny = df.loc[method,"ny"]
        extend_n_px = df.loc[method,"extend_n_px"]
        files = []
        seg_str = "" if seg is None else (seg+"_")
        for x in range(nx):
            for y in range(ny):
                files.append(
                    os.path.join(rep_dir, f"assignments_{seg_str}{method}-{id_code}_ny{ny}_nx{nx}_{y}_{x}_px{extend_n_px}.csv")
                )
                if areas:
                    files.append(
                        os.path.join(rep_dir, f"areas_{seg_str}{method}-{id_code}_ny{ny}_nx{nx}_{y}_{x}_px{extend_n_px}.csv")
                    )
        return files
    
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
                #TODO include chunks?
                if 'group_metrics' in final_file and 'aggregated' not in final_file and ('-'+wildcards.id_code) in final_file:
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
            #raise Exception(f"Cannot find method: {method}") 
        if id_code >= len(self.method_dict[method]):
            return None
            #raise Exception(f"Cannot find index {idx} in method '{method}'")
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
                
        # Test that all entries have type str in cell type expert annotations
        if "celltype" in spots.columns:
            for val in spots["celltype"].unique():
                assert isinstance(val, str), f"Cell type entry '{val}' in column 'celltype' in file {spots_file} is not a str."
        
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
