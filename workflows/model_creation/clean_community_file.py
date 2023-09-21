import cobra
from cobra import Metabolite, Model, Reaction


def clean_community(model: cobra.Model) -> cobra.Model:
    """
    Clean a combined community model and create fecal and diet transport and exchange reactions.

    For each general metabolite in the lumen, the function adds four types of reactions:
    1. Diet Exchange: EX_2omxyl[d]: 2omxyl[d] <=>
    2. Diet Transport: DUt_2omxyl: 2omxyl[d] <=> 2omxyl[u]
    3. Fecal Transport: UFEt_2omxyl: 2omxyl[u] <=> 2omxyl[fe]
    4. Fecal Exchange: EX_2omxyl[fe]: 2omxyl[fe] <=>

    Parameters:
    - model (Model): An AGORA single cell model.

    Returns:
    - Model: Updated model with fecal and diet compartments.
    """

    # Delete all of those EX_ reaction artifacts from the single cell models
    # (reactions such as: EX_dad_2(e): dad_2[e] <=>, EX_thymd(e): thymd[e] <=> )
    while (
        len([reac for reac in model.reactions if "_EX_" in reac.id or "(e)" in reac.id])
        > 0
    ):
        for reac in model.reactions:
            if "_EX_" in reac.id or "(e)" in reac.id:
                model.reactions.get_by_id(reac.id).remove_from_model()

    ## now we can create the diet and fecal compartments (reactions and metabolites!)
    ## lets get all of our general extracellular metabolites:
    gen_mets = []

    for reac in model.reactions:
        if "IEX" in reac.id:
            gen_mets.append(
                (model.reactions.get_by_id(reac.id).reaction).split(" <=> ")[0]
            )

    gen_mets = set(gen_mets)

    ## creating DIET compartments
    for met_name in gen_mets:
        # just the metabolite
        clean_met_name = met_name.split("[")[0]

        # metabolite[d]
        d_met_name = f"{clean_met_name}[d]"

        # reaction name DUt_metabolite
        DUt_name = f"DUt_{clean_met_name}"

        ## creating the d exchange reactions (EX_2omxyl[d]: 2omxyl[d] <=>)
        # making sure we don't double up!! (if user runs function more than once)
        if d_met_name not in [met.id for met in model.metabolites]:
            reac_name = f"EX_{d_met_name}"
            reaction = Reaction(reac_name)
            reaction.name = f"{d_met_name}diet exchange"
            reaction.subsystem = " "
            reaction.lower_bound = -1000.0  # This is the default
            reaction.upper_bound = 1000.0  # This is the default

            model.add_reactions([reaction])

            model.add_metabolites(
                [
                    Metabolite(str(d_met_name), formula=" ", name="", compartment="d"),
                ]
            )

            new_dietreact = model.reactions.get_by_id(reac_name)
            new_dietreact.add_metabolites({model.metabolites.get_by_id(d_met_name): -1})

        ## creating the d transport reactions to lumen (DUt_4hbz: 4hbz[d] --> 4hbz[u])
        # making sure we don't double up!! (if user runs function more than once)
        if DUt_name not in [reac.id for reac in model.reactions]:
            DUt_formula = f"{d_met_name} --> {met_name}"

            reac_name = DUt_name
            reaction = Reaction(reac_name)
            reaction.name = f"{DUt_name}diet to lumen"
            reaction.subsystem = " "
            reaction.lower_bound = 0.0
            reaction.upper_bound = 1000.0

            model.add_reactions([reaction])

            # adding the correct d --> u formula to the reaction:
            reaction.reaction = DUt_formula

    ## creating FECAL compartments
    for met_name in gen_mets:
        # just the metabolite
        clean_met_name = met_name.split("[")[0]

        # metabolite[d]
        fe_met_name = f"{clean_met_name}[fe]"

        # reaction name UFEt_metabolite
        UFEt_name = f"UFEt_{clean_met_name}"

        ## creating the fe exchange reactions EX_4abut[fe]: 4abut[fe] <=>)
        # making sure we don't double up!! (if user runs function more than once)
        if fe_met_name not in [met.id for met in model.metabolites]:
            reac_name = f"EX_{fe_met_name}"
            reaction = Reaction(reac_name)
            reaction.name = f"{fe_met_name}fecal exchange"
            reaction.subsystem = " "
            reaction.lower_bound = -1000.0  # This is the default
            reaction.upper_bound = 1000.0  # This is the default

            model.add_reactions([reaction])

            model.add_metabolites(
                [
                    Metabolite(
                        str(fe_met_name), formula=" ", name="", compartment="fe"
                    ),
                ]
            )

            new_fereact = model.reactions.get_by_id(reac_name)
            new_fereact.add_metabolites({model.metabolites.get_by_id(fe_met_name): -1})

        ###creating the fe transport reactions to lumen (UFEt_arabinoxyl: arabinoxyl[u] --> arabinoxyl[fe])

        # making sure we don't double up!! (if user runs function more than once)
        if UFEt_name not in [reac.id for reac in model.reactions]:
            UFEt_formula = f"{met_name} --> {fe_met_name}"

            reac_name = UFEt_name
            reaction = Reaction(reac_name)
            reaction.name = f"{UFEt_name}diet to lumen"
            reaction.subsystem = " "
            reaction.lower_bound = 0.0  # this is what mg pipe does, check on logic????
            reaction.upper_bound = 1000.0  # This is the default

            model.add_reactions([reaction])

            # adding the correct d --> u formula to the reaction:
            reaction.reaction = UFEt_formula

    return model
