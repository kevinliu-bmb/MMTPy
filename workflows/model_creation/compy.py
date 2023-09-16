######BIG BOY MAIN PIPELINE

 

#required packages

from clean_community_file import clean_community
from com_biomass_file import com_biomass
from species_to_community_file import species_to_community
import cobra
from cobra import io
from cobra.io import load_matlab_model, load_json_model, save_json_model, save_matlab_model
import pandas as pd
from cobra import Model, Reaction, Metabolite
import csv

def compy(abunpath, modpath, respath, dietpath=None):
    
    #inspired by MICOM community building and mgpipe.m code
        
    """
    master pipeline which inputs the GEMs data and accesses the different functions: 


    INPUTS:
    abun_path: path to the species abundance .csv file
    formatting for the species abundance:
        the columns should have the names of the .mat files of the species you want to load
        (ex: normCoverageReduced)
        see file normCoverageReduced.csv in shared dropbox for example
        
    modpath: path to folder with all AGORA models 
        EX: '/Users/isabellagm/Desktop/research/MGH/AGORA-1.03/AGORA1mat/'
    respath: path where the community models will be output (defaults to same folder as )
        EX: "/Users/isabellagm/Desktop/working_compy_models/"
    dietpath: path to the AGORA compatible diet (for community model) .csv file
  
    OUTPUTS:
    all sample community models to a local folder
  
    """
    
    
    print("starting the ComPy pipeline!")
    #setting up importing modules:
    import cobra
    import pandas as pd
    from cobra import Model, Reaction, Metabolite
    import os

    print("loading in the abundance info!")
    #loading the abundance input file and formatting it correctly     
    sampleinfo = pd.read_csv(abunpath)
    sampleinfo.rename(columns = {list(sampleinfo)[0]:'species'}, inplace=True)
    sampleinfo.set_index('species', inplace=True)
    sample_list = sampleinfo.columns.tolist()
    
   
    #setting solver and model configurations 
    solver = [
                s
                for s in ["cplex", "gurobi", "osqp", "glpk"]
                if s in cobra.util.solver.solvers
            ][0]
    
    '''
    joining models 
    '''
    
####iterate through each column in sample info to look at each sample:
    
    for sample in sampleinfo.columns:

        print("working on sample: " + str(sample))
    #for 1 specific community sample
    
        #dic with model path as key and species abun as value
        model_path_abun = {}
        
        #finding the species that have an abundance above 0.001 and 
        #adding them to the model_path_abun_dic so its species name : abundace
        
        for num in range(len(sampleinfo[sample])):
            if sampleinfo[sample][num] > 0.001:
                species_name = sampleinfo.index[num]
                model_path_abun[modpath + species_name + '.mat'] = sampleinfo[sample][num]
            else:
                continue
        print("starting initial model")

        first_species = list(model_path_abun.keys())[0]
        print(first_species)
        first_model = cobra.io.load_matlab_model(first_species)
        print("nowwww")
        final_model = species_to_community(model=first_model, species_model_name=first_species)
        
        print("finished intial model!")
        
        #now we loop throop the remaining species model in this sample and merge them 
        #into the (combined) final_model:
        
        for species_model in list(model_path_abun.keys())[1:]:
            model = cobra.io.load_matlab_model(species_model)
            adding_model = species_to_community(model=model, species_model_name=species_model)
            rids = set(r.id for r in final_model.reactions)
            new = [r.id for r in adding_model.reactions if r.id not in rids]
            final_model.add_reactions(adding_model.reactions.get_by_any(new))
            print("finished adding " + str(species_model) +  " model!")
            
######now we need to add the diet + fecal compartments to the big boy model! 
        
        clean_model = clean_community(model = final_model)
    
    
#######now we need to add a community biomass:
        clean_model = com_biomass(model=clean_model, abunpath=abunpath, sample_com=sample)
    
            
        #ensuring the reversablity fits all compartments
        
        # for reac in clean_model.reactions:
            # clean_model.reactions.get_by_id(reac.id).lower_bound = -1000

        for reac in [i for i in clean_model.reactions if "DUt" in i.id]:
            clean_model.reactions.get_by_id(reac.id).lower_bound = 0
    
        for reac in [i for i in clean_model.reactions if "UFEt" in i.id]:
            clean_model.reactions.get_by_id(reac.id).lower_bound = 0

        #setting the given diet!!
        if dietpath is not None:
            diet_dic = {}
            with open(dietpath) as f:
                next(f)
                for key, *values in csv.reader(f):
                    diet_dic[key.strip()] = values[0]
                    
            for reac in [i.id for i in clean_model.reactions if "EX_" and "[d]" in i.id]:
                if reac in diet_dic.keys():
                    clean_model.reactions.get_by_id(reac).lower_bound = float(diet_dic[reac])
                    print("metabolite " + str(reac) + "has been added to diet!")
                else:
                    clean_model.reactions.get_by_id(reac).lower_bound = 0

        #setting the new community biomass!
        
        clean_model.objective = "communityBiomass"
        if not os.path.exists(respath):
            os.makedirs(respath)
        
        microbiome_model =  respath + sample + "_communitymodel_final.json"

        cobra.io.save_json_model(clean_model, microbiome_model)


    return clean_model 
