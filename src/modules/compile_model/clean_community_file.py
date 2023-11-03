from cobra import Reaction, Metabolite


def change_lb(reaction, lb):
    reaction.subsystem = " "
    reaction.lower_bound = lb
    reaction.upper_bound = 1000.0


def create_reaction(model, rxn_name, lb, name):
    rxn = Reaction(rxn_name)
    change_lb(rxn, lb)
    rxn.name = name
    model.add_reactions([rxn])
    return rxn


def create_metabolite(model, met_id, compartment):
    met = Metabolite(id=met_id, formula=" ", name="", compartment=compartment)
    model.add_metabolites([met])
    return met


def rename_rxn(model):
    model.reactions = [
        rxn for rxn in model.reactions if "_EX_" not in rxn.id and "(e)" not in rxn.id
    ]

    met_comp_list = {
        (rxn.reaction).split(" <=> ")[0] for rxn in model.reactions if "_IEX_" in rxn.id
    }

    for met_id in met_comp_list:
        met_id = met_id.split("[")[0]
        met_id_comp = f"{met_id}[d]"
        DUt_name = f"DUt_{met_id}"

        if met_id_comp not in model.metabolites:
            EX_formula = f"EX_{met_id_comp}"
            rxn = create_reaction(
                model, EX_formula, -1000.0, f"{met_id_comp}diet exchange"
            )
            met = create_metabolite(model, met_id_comp, "d")
            rxn.add_metabolites({met: -1})

        if DUt_name not in model.reactions:
            DUt_formula = f"{met_id_comp} --> {met_id}"
            rxn = create_reaction(model, DUt_name, 0.0, f"{DUt_name}diet to lumen")
            rxn.reaction = DUt_formula

    for met_id in met_comp_list:
        met_id = met_id.split("[")[0]
        met_name_fe = f"{met_id}[fe]"
        UFEt_name = f"UFEt_{met_id}"

        if met_name_fe not in model.metabolites:
            EX_formula = f"EX_{met_name_fe}"
            rxn = create_reaction(
                model, EX_formula, -1000.0, f"{met_name_fe}fecal exchange"
            )
            met = create_metabolite(model, met_name_fe, "fe")
            rxn.add_metabolites({met: -1})

        if UFEt_name not in model.reactions:
            UFEt_formula = f"{met_id} --> {met_name_fe}"
            rxn = create_reaction(model, UFEt_name, 0.0, f"{UFEt_name}diet to lumen")
            rxn.reaction = UFEt_formula

    return model
