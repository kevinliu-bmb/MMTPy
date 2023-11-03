import os
import pandas as pd
import cobra
from cobra import Reaction, Metabolite


def get_communityBiomass(model_input, mgx_path):
    if isinstance(model_input, cobra.Model):
        model = model_input
    elif isinstance(model_input, str):
        if os.path.isfile(model_input):
            model = cobra.io.load_matlab_model(model_input)
        else:
            raise TypeError(
                "model_input is not a valid file path or cobra.Model object"
            )

    while [rxn for rxn in model.reactions if "Biomass" in rxn.id]:
        for rxn in model.reactions:
            if "Biomass" in rxn.id:
                model.reactions.get_by_id(rxn.id).remove_from_model()

    biomass_mets_list = [mets.id for mets in model.metabolites if "biomass" in mets.id]
    # sort species alphabetically!
    biomass_mets_list = sorted(biomass_mets_list)

    # reading in the abundance file, sorting alphabetically:
    normCoverage = pd.read_csv(mgx_path)
    normCoverage = normCoverage.sort_values("X", ascending=True)
    normCoverage = normCoverage.reset_index(drop=True)
    normCoverage = normCoverage.loc[~((normCoverage[model.name] < 0.0009))]
    normCoverage_list = normCoverage[model.name].tolist()

    # creating the com bio reaction
    reaction = Reaction("communityBiomass")
    _extracted_from_get_communityBiomass_(reaction, 0)
    reaction.name = "community biomass "
    # reaction.add_metabolites(com_biomass_dic)

    model.add_reactions([reaction])

    communityBiomass = model.reactions.communityBiomass

    com_biomass_dic = {
        biomass: -float(normCoverage_list[counter])
        for counter, biomass in enumerate(biomass_mets_list)
    }
    communityBiomass.add_metabolites(metabolites_to_add=com_biomass_dic, combine=True)

    # adding a microbeBiomass metabolite
    model.add_metabolites(
        [
            Metabolite(
                "microbeBiomass[u]",
                formula=" ",
                name="product of community biomass",
                compartment="u",
            ),
        ]
    )
    communityBiomass.add_metabolites(
        {model.metabolites.get_by_id("microbeBiomass[u]"): 1}
    )

    ####adding the EXCHANGE REACTION compartment:

    reac_name = "EX_" + "microbeBiomass[fe]"
    reaction = Reaction(reac_name)
    reaction.name = f"{str(reac_name)}fecal exchange"
    _extracted_from_get_communityBiomass_(reaction, 0)
    model.add_reactions([reaction])

    # adding a microbeBiomass fe metabolite
    model.add_metabolites(
        [
            Metabolite(
                "microbeBiomass[fe]",
                formula=" ",
                name="product of community biomass",
                compartment="fe",
            ),
        ]
    )

    new_fereact = model.reactions.get_by_id("EX_microbeBiomass[fe]")
    new_fereact.add_metabolites({model.metabolites.get_by_id("microbeBiomass[fe]"): -1})

    ####adding the UFEt REACTION :

    UFEt_formula = "microbeBiomass[u] --> microbeBiomass[fe]"

    reac_name = "UFEt_microbeBiomass"
    reaction = Reaction(reac_name)
    reaction.name = "UFEt_microbeBiomass" + "diet to lumen"
    _extracted_from_get_communityBiomass_(reaction, 0.0)
    model.add_reactions([reaction])

    # adding the correct d --> u formula to the reaction:
    reaction.reaction = UFEt_formula

    return model


# TODO Rename this here and in `get_communityBiomass`
def _extracted_from_get_communityBiomass_(reaction, arg1):
    reaction.subsystem = " "
    reaction.lower_bound = arg1
    reaction.upper_bound = 1000.0  # This is the default
