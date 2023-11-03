import os
import csv

import cobra
import pandas as pd

from clean_community_file import rename_rxns
from modules.compile_model.get_communityBiomass import get_communityBiomass
from species_to_community_file import species_to_community


def set_lower_bound(reactions, lower_bound):
    for rxn in reactions:
        rxn.lower_bound = lower_bound


def make_community_model(mgx_path, model_path, output_dir, diet_path):
    mgx_df = pd.read_csv(mgx_path)
    mgx_df.rename(columns={list(mgx_df)[0]: "taxonomy"}, inplace=True)
    mgx_df.set_index("taxonomy", inplace=True)

    for sample in mgx_df.columns:
        model_path_abund_dict = {
            os.path.join(model_path, f"{mgx_df.index[i]}.mat"): mgx_df[sample][i]
            for i in range(len(mgx_df[sample]))
            if mgx_df[sample][i] > 0.001
        }

        agora_model = list(model_path_abund_dict.keys())[0]
        model = cobra.io.load_matlab_model(agora_model)
        model = species_to_community(model=model, agora_model=agora_model)

        for species_model in list(model_path_abund_dict.keys())[1:]:
            model = cobra.io.load_matlab_model(species_model)
            agora_model_add = species_to_community(
                model=model, agora_model=species_model
            )
            rxn_id_list = {r.id for r in model.reactions}
            new = [r.id for r in agora_model_add.reactions if r.id not in rxn_id_list]
            model.add_reactions(agora_model_add.reactions.get_by_any(new))

        model.name = sample
        model = rename_rxns(model_input=model)
        model = get_communityBiomass(model_input=model, mgx_path=mgx_path)

        set_lower_bound([i for i in model.reactions if "DUt" in i.id], 0)
        set_lower_bound([i for i in model.reactions if "UFEt" in i.id], 0)

        if diet_path is not None:
            with open(diet_path) as f:
                diet_dict = {key.strip(): values[0] for key, *values in csv.reader(f)}

            for rxn in [i.id for i in model.reactions if "EX_" and "[d]" in i.id]:
                if rxn in diet_dict:
                    model.reactions.get_by_id(rxn).lower_bound = float(diet_dict[rxn])
                else:
                    model.reactions.get_by_id(rxn).lower_bound = 0

        model.objective = "communityBiomass"

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        multispecies_model_output_path = os.path.join(
            output_dir, f"{sample}_multispecies.xml"
        )
        cobra.io.write_sbml_model(model, multispecies_model_output_path)

    return model
