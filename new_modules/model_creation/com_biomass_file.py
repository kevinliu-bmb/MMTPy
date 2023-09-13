
import cobra
from cobra import io
from cobra.io import load_matlab_model, load_json_model, save_json_model, save_matlab_model
import pandas as pd
from cobra import Model, Reaction, Metabolite

def com_biomass(model, abunpath, sample_com):
    
    """
    takes a combined community model and adds a community biomass formula 

        

    INPUTS:
    model: a .mat file of an AGORA single cell model
    abun_path: path to the species abundance .csv file
    sample_com: the sample name string (internal to the com_py pipeline)
  
    OUTPUTS:
    model: updated model with community biomass 
  
    """
    
###deleting all previous community biomass equations:
    while len([reac for reac in model.reactions if "Biomass" in reac.id]) > 0:
        for reac in model.reactions:
            if "Biomass" in reac.id:
                model.reactions.get_by_id(reac.id).remove_from_model()

    #extracting biomass metabolites from the different single cell models
    biomass_mets_list = []
    for mets in model.metabolites:
        if "biomass" in mets.id:
            biomass_mets_list.append(mets.id)

    #sort species alphabetically!
    biomass_mets_list = sorted(biomass_mets_list)

    #reading in the abundance file, sorting alphabetically:
    normCoverage = pd.read_csv(abunpath) 
    normCoverage = normCoverage.sort_values('X', ascending=True)
    normCoverage = normCoverage.reset_index(drop=True)
    normCoverage = normCoverage.loc[~((normCoverage[sample_com] < 0.0009))]
    normCoverage_list = normCoverage[sample_com].tolist()



    #creating the com bio reaction
    reaction = Reaction('communityBiomass')
    reaction.name = 'community biomass '
    reaction.subsystem = ' '
    reaction.lower_bound = 0  # This is the default
    reaction.upper_bound = 1000.  # This is the default
    #reaction.add_metabolites(com_biomass_dic)

    model.add_reactions([reaction])

    communityBiomass = model.reactions.communityBiomass

    counter = 0
    #adding the mets with their weighted abundances
    com_biomass_dic = {}

    for biomass in biomass_mets_list:
        com_biomass_dic[biomass] = -float(normCoverage_list[counter])
        counter += 1


    communityBiomass.add_metabolites(metabolites_to_add=com_biomass_dic, combine=True)

    #adding a microbeBiomass metabolite
    model.add_metabolites([
    Metabolite(
    'microbeBiomass[u]',
    formula=' ',
    name='product of community biomass',
    compartment='u'
    ),
])
    communityBiomass.add_metabolites({model.metabolites.get_by_id('microbeBiomass[u]'): 1})


####adding the EXCHANGE REACTION compartment:



    reac_name = "EX_" + "microbeBiomass[fe]"
    reaction = Reaction(reac_name)
    reaction.name = str(reac_name) + "fecal exchange"
    reaction.subsystem = ' '
    reaction.lower_bound = 0  # This is the default
    reaction.upper_bound = 1000.  # This is the default

    model.add_reactions([reaction])


 #adding a microbeBiomass fe metabolite
    model.add_metabolites([
    Metabolite(
    'microbeBiomass[fe]',
    formula=' ',
    name='product of community biomass',
    compartment='fe'
    ),
])


    new_fereact = model.reactions.get_by_id("EX_microbeBiomass[fe]")
    new_fereact.add_metabolites({model.metabolites.get_by_id("microbeBiomass[fe]"): -1})

####adding the UFEt REACTION :

    UFEt_formula = (f'microbeBiomass[u] --> microbeBiomass[fe]')

    reac_name = "UFEt_microbeBiomass"
    reaction = Reaction(reac_name)
    reaction.name = str("UFEt_microbeBiomass") + "diet to lumen"
    reaction.subsystem = ' '
    reaction.lower_bound = 0.  # this is what mg pipe does, check on logic????
    reaction.upper_bound = 1000.  # This is the default

    model.add_reactions([reaction])

    #adding the correct d --> u formula to the reaction:
    reaction.reaction = UFEt_formula  
        
        
    return model
