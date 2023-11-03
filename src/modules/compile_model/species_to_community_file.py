def update_name(model, old_name, new_name, compartment):
    model.get_by_id(old_name).id = new_name
    model.get_by_id(new_name).compartment = compartment


def remove_ex_rxn(model):
    for rxn in model.reactions:
        if "EX_" in rxn.id and "biomass" not in rxn.id:
            model.remove_reactions(model.reactions.get_by_id(rxn.id))


def create_IEX_reactions(model, short_species_name):
    for met in model.metabolites:
        if "[u]" in met.id and short_species_name in met.id:
            replacing_species_name = f"{short_species_name}_"
            general_name = met.id.replace(replacing_species_name, "")
            IEX_formula = f"{general_name} <=> {met.id}"

            IEX_reaction_name = f"{short_species_name}_IEX_{general_name}tr"

            reaction = Reaction(IEX_reaction_name)
            reaction.name = f"{short_species_name}_IEX"
            reaction.lower_bound = -1000.0
            reaction.upper_bound = 1000.0
            reaction.subsystem = ""
            model.add_reactions([reaction])

            reaction.reaction = IEX_formula


def rename_met_comp(model, species_name_short, comp):
    for met in model.metabolites:
        if f"[{comp}]" in met.id and species_name_short not in met.id:
            met_name = met.id
            met_name_comp_rm = met_name.replace(f"[{comp}]", "")
            met_name_update = f"{species_name_short}_{met_name_comp_rm}[{comp}]"
            update_name(model, met_name, met_name_update, comp)


def species_to_community(model, agora_model):
    species_name_short = agora_model.split("/")[-1]
    species_name_short = species_name_short.split(".")[0]

    remove_ex_rxn(model)

    for rxn in model.reactions:
        if "[e]" not in rxn.reaction and "[c]" in rxn.reaction:
            reaction_bounds = rxn.bounds
            rxn_name = rxn.id
            rxn_name_update = f"{species_name_short}_{rxn_name}"
            update_name(model, rxn_name, rxn_name_update, "c")

            for met in rxn.metabolites:
                if species_name_short not in str(met.id):
                    rename_met_comp(model, species_name_short, comp="c")
                    rename_met_comp(model, species_name_short, comp="p")

            model.reactions.get_by_id(rxn_name_update).bounds = reaction_bounds

    create_IEX_reactions(model, species_name_short)

    for rxn in model.reactions:
        if species_name_short not in rxn.id:
            rxn_name = rxn.id
            rxn_name_update = f"{species_name_short}_{rxn_name}"
            update_name(model, rxn_name, rxn_name_update, "c")

    rename_met_comp(model, species_name_short, comp="c")
    rename_met_comp(model, species_name_short, comp="p")

    return model
