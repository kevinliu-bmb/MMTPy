#reactions get incorporated into the larger model

#should end up with all reactions tagged with this specific species, all should be added to com model:
import cobra
from cobra import io
from cobra.io import load_matlab_model, load_json_model, save_json_model, save_matlab_model
import pandas as pd
from cobra import Model, Reaction, Metabolite


def species_to_community(model, species_model_name):
    
    """
    takes a single cell AGORA GEM and changes its reaction and metabolite formatting so it 
    can be added to a community model in the Com_py pipeline 
  
        tags everything intracellular and intracellular to extracellular with the species name:
        (intracellular)
        tag[c] -> tag[c]

        (transport)
        tag[c] -> tag[u]

        (IEX reactions)
        tagged[e] -> general[e]
  
    INPUTS:
    model: a .mat file of an AGORA single cell model
    species_model_name: the species name to be tagged (extracted in the Com_py pipeline)
  
    OUTPUTS:
    model: updated model
  
    """
    
    #####TAGGING all reactions and metabolites in the cell with species name 
    #for each species model, we need to iterate through its reactions and add the species tag
    
    short_species_name = species_model_name.split("/")[-1]
    short_species_name = short_species_name.split(".")[0]
    
    
    
    #FIRST removing all exchange reactions, we will make our own!!!
    for r in model.reactions:
        if "EX_" in r.id and "biomass" not in r.id:
            model.remove_reactions(model.reactions.get_by_id(r.id))
    
    ##THE INTRA CELL reactions [c] -> [c]
    for r in model.reactions:
        if "[e]" not in r.reaction:
            if "[c]" in r.reaction:
                reaction_bounds = r.bounds
                reaction_name = r.id
                updated_r_name =  short_species_name + "_" + reaction_name
                model.reactions.get_by_id(reaction_name).id = updated_r_name


                #for each metabolite in that reaction!
                for m in r.metabolites:

                    #if that met already has the name tag skip it
                    if short_species_name not in str(m.id):
                        if "[c]" in str(m.id):
                            metabolite_name = m.id
                            model.metabolites.get_by_id(metabolite_name).compartment = "c"
                            no_c_name = metabolite_name.replace("[c]", "")
                            updated_m_name = short_species_name + "_" + no_c_name + "[c]"
                            model.metabolites.get_by_id(metabolite_name).id = updated_m_name
                            model.metabolites.get_by_id(updated_m_name).compartment = "c"
                        if "[p]" in str(m.id):
                            metabolite_name = m.id
                            model.metabolites.get_by_id(metabolite_name).compartment = ""
                            no_c_name = metabolite_name.replace("[p]", "")
                            updated_m_name = short_species_name + "_" + no_c_name + "[p]"
                            model.metabolites.get_by_id(metabolite_name).id = updated_m_name
                            model.metabolites.get_by_id(updated_m_name).compartment = "p"
                        
              
                #trying to get the reversiblity to carry over in my new models:
                model.reactions.get_by_id(updated_r_name).bounds = reaction_bounds
                print(model.reactions.get_by_id(updated_r_name).bounds)



    ##THE EXTRA CELL reactions [c] -> [e]
    
    for r in model.reactions:
        if "[e]" in r.reaction:
            reaction_name = r.id
            updated_r_name =  short_species_name + "_" + reaction_name
            model.reactions.get_by_id(reaction_name).id = updated_r_name
            
            #for each metabolite in that reaction!
            for m in r.metabolites:
                
                #if it is the cellular metabolite:
                
                if "[c]" in m.id:
                    #if that met already has the name tag skip it
                    if short_species_name not in str(m.id):
                        metabolite_name = m.id
                        model.metabolites.get_by_id(metabolite_name).compartment = "c"
                        no_c_name = metabolite_name.replace("[c]", "")
                        updated_m_name = short_species_name + "_" + no_c_name + "[c]"
                        model.metabolites.get_by_id(metabolite_name).id = updated_m_name
                        model.metabolites.get_by_id(updated_m_name).compartment = "c"
                
                if "[e]" in m.id:
                    if short_species_name not in str(m.id):
                        metabolite_name = m.id
                        #doing the compartments now changing to u
                        model.metabolites.get_by_id(metabolite_name).compartment = "u"
                        no_c_name = metabolite_name.replace("[e]", "")
                        updated_m_name = short_species_name + "_" + no_c_name + "[u]"
                        model.metabolites.get_by_id(metabolite_name).id = updated_m_name
                        model.metabolites.get_by_id(updated_m_name).compartment = "u"


                if "[p]" in m.id:
                    if short_species_name not in str(m.id):
                        metabolite_name = m.id
                        #doing the compartments now changing to u
                        model.metabolites.get_by_id(metabolite_name).compartment = "p"
                        no_c_name = metabolite_name.replace("[p]", "")
                        updated_m_name = short_species_name + "_" + no_c_name + "[p]"
                        model.metabolites.get_by_id(metabolite_name).id = updated_m_name
                        model.metabolites.get_by_id(updated_m_name).compartment = "p"
                        
                
                        
     ##THE IEX reactions tagged[e] -> general[e]
    
    #for all extracellular species specific metabolites
    for met in model.metabolites:
        if "[u]" in met.id:
            if short_species_name in met.id:
                replacing_specname = short_species_name + "_"
                general_name = met.id.replace(replacing_specname, '')
                IEX_formula = (f'{general_name} <=> {met.id}')

                IEX_reaction_name = short_species_name + "_IEX_" + general_name + "tr"

                reaction = Reaction(IEX_reaction_name)
                reaction.name = short_species_name + "_IEX"
                reaction.lower_bound = 0
                reaction.upper_bound = 1000
                reaction.subsystem = ""
                model.add_reactions([reaction])

                reaction.reaction = IEX_formula

    #making sure every reaction has the species tag on it:
    for r in model.reactions:
        if short_species_name not in r.id:
            reaction_name = r.id
            updated_r_name =  short_species_name + "_" + reaction_name
            model.reactions.get_by_id(reaction_name).id = updated_r_name
            
    #making sure every (c) compartment metabolite has the species tag:
    for m in model.metabolites:
        if "[c]" in m.id:
                if short_species_name not in m.id:
                    metabolite_name = m.id
                    model.metabolites.get_by_id(metabolite_name).compartment = "c"
                    no_c_name = metabolite_name.replace("[c]", "")
                    updated_m_name = short_species_name + "_" + no_c_name + "[c]"
                    model.metabolites.get_by_id(metabolite_name).id = updated_m_name
                    model.metabolites.get_by_id(updated_m_name).compartment = "c"
    
        #making sure every (p) compartment metabolite has the species tag:
    for m in model.metabolites:
        if "[p]" in m.id:
                if short_species_name not in m.id:
                    metabolite_name = m.id
                    model.metabolites.get_by_id(metabolite_name).compartment = "p"
                    no_p_name = metabolite_name.replace("[p]", "")
                    updated_m_name = short_species_name + "_" + no_p_name + "[p]"
                    model.metabolites.get_by_id(metabolite_name).id = updated_m_name
                    model.metabolites.get_by_id(updated_m_name).compartment = "p"

    
    #doing the individual biomass reaction!
    #each species will have a c --> c reaction that will then 
    for r in model.reactions:
        if short_species_name not in r.id:
            reaction_name = r.id
         # species_to_community.py   updated_r_name =  short_species_name + "_" + reaction_name
            model.reactions.get_by_id(reaction_name).id = updated_r_name
    
    return model
