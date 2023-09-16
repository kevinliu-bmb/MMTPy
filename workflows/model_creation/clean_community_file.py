
import cobra
from cobra import io
from cobra.io import load_matlab_model, load_json_model, save_json_model, save_matlab_model
import pandas as pd
from cobra import Model, Reaction, Metabolite

def clean_community(model):
    

    """
    takes a combined community model and creates fecal and diet transport and exchange reactions 

  
        adds 4 types of reactions for every general metabolite in the lumen:
        (diet)
        EX_2omxyl[d]: 2omxyl[d] <=>
        DUt_2omxyl: 2omxyl[d] <=> 2omxyl[u]
        
        (fecal)
        UFEt_2omxyl: 2omxyl[u] <=> 2omxyl[fe]
        EX_2omxyl[fe]: 2omxyl[fe] <=>
        

    INPUTS:
    model: a .mat file of an AGORA single cell model
  
    OUTPUTS:
    model: updated model with fecal and diet compartments 
  
    """


####lets delete all of those EX_ reaction artifacts from the single cell models 
    #(reactions such as: EX_dad_2(e): dad_2[e] <=>, EX_thymd(e): thymd[e] <=> )
    
    while len([reac for reac in model.reactions if "_EX_" in reac.id or "(e)" in reac.id]) > 0:
        for reac in model.reactions:
            if "_EX_" in reac.id or "(e)" in reac.id:
                model.reactions.get_by_id(reac.id).remove_from_model()
                
                
 ###### now we can create the diet and fecal compartments (reactions and metabolites!)
    ###lets get all of our general extracellular metabolites:
    
    gen_mets = []
    
    for reac in model.reactions:
        if "IEX" in reac.id:
            gen_mets.append((model.reactions.get_by_id(reac.id).reaction).split(" <=> ")[0])

    gen_mets = set(gen_mets)
    
#####creating DIET compartments
    for met_name in gen_mets:
        
        #just the metabolite
        clean_met_name = met_name.split("[")[0]
        
        #metabolite[d]
        d_met_name =  str(clean_met_name) + "[d]"
        
        #reaction name DUt_metabolite
        DUt_name = "DUt_" + clean_met_name

        
###creating the d exchange reactions (EX_2omxyl[d]: 2omxyl[d] <=>)
    
        #making sure we don't double up!! (if user runs function more than once)
        if d_met_name not in [met.id for met in model.metabolites]:

            reac_name = "EX_" + str(d_met_name)
            reaction = Reaction(reac_name)
            reaction.name = str(d_met_name) + "diet exchange"
            reaction.subsystem = ' '
            reaction.lower_bound = -1000.  # This is the default
            reaction.upper_bound = 1000.  # This is the default

            model.add_reactions([reaction])


            model.add_metabolites([
            Metabolite(
            str(d_met_name),
            formula=' ',
            name='',
            compartment='d'
            ),
        ])


            new_dietreact = model.reactions.get_by_id(reac_name)
            new_dietreact.add_metabolites({model.metabolites.get_by_id(d_met_name): -1})
        
###creating the d transport reactions to lumen (DUt_4hbz: 4hbz[d] --> 4hbz[u])
        
    #making sure we don't double up!! (if user runs function more than once)
        if DUt_name not in [reac.id for reac in model.reactions]:
            
            DUt_formula = (f'{d_met_name} --> {met_name}')
            
            reac_name = DUt_name
            reaction = Reaction(reac_name)
            reaction.name = str(DUt_name) + "diet to lumen"
            reaction.subsystem = ' '
            reaction.lower_bound = 0.  # this is what mg pipe does, check on logic????
            reaction.upper_bound = 1000.  # This is the default

            model.add_reactions([reaction])

            #adding the correct d --> u formula to the reaction:
            reaction.reaction = DUt_formula
            
        
#####creating FECAL compartments

    for met_name in gen_mets:
        
        #just the metabolite
        clean_met_name = met_name.split("[")[0]
        
        #metabolite[d]
        fe_met_name =  str(clean_met_name) + "[fe]"
        
         #reaction name UFEt_metabolite
        UFEt_name = "UFEt_" + clean_met_name
        
###creating the fe exchange reactions EX_4abut[fe]: 4abut[fe] <=>)

        #making sure we don't double up!! (if user runs function more than once)
        if fe_met_name not in [met.id for met in model.metabolites]:
        
            reac_name = "EX_" + str(fe_met_name)
            reaction = Reaction(reac_name)
            reaction.name = str(fe_met_name) + "fecal exchange"
            reaction.subsystem = ' '
            reaction.lower_bound = -1000.  # This is the default
            reaction.upper_bound = 1000.  # This is the default

            model.add_reactions([reaction])


            model.add_metabolites([
            Metabolite(
            str(fe_met_name),
            formula=' ',
            name='',
            compartment='fe'
            ),
        ])


            new_fereact = model.reactions.get_by_id(reac_name)
            new_fereact.add_metabolites({model.metabolites.get_by_id(fe_met_name): -1})
    
###creating the fe transport reactions to lumen (UFEt_arabinoxyl: arabinoxyl[u] --> arabinoxyl[fe])

            
    #making sure we don't double up!! (if user runs function more than once)
        if UFEt_name not in [reac.id for reac in model.reactions]:
            
            UFEt_formula = (f'{met_name} --> {fe_met_name}')
            
            reac_name = UFEt_name
            reaction = Reaction(reac_name)
            reaction.name = str(UFEt_name) + "diet to lumen"
            reaction.subsystem = ' '
            reaction.lower_bound = 0.  # this is what mg pipe does, check on logic????
            reaction.upper_bound = 1000.  # This is the default

            model.add_reactions([reaction])

            #adding the correct d --> u formula to the reaction:
            reaction.reaction = UFEt_formula
                
            
    return model
