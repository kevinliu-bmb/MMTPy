import os
import sys
from math import isnan

import cobra
import pandas as pd

from match_names_to_vmh import match_names_to_vmh


def print_logo(tool: str, tool_description: str, version: str):
    """
    Print the logo, tool name, and version for the tool.

    Parameters
    ----------
    tool : str
        The name of the tool.
    tool_description : str
        The description of the tool.
    version : str
        The version of the tool.

    Returns
    -------
    None
    """
    logo = r"""
     ____      _          _       _____           _     
    |__  /    | |    __ _| |__   |_   _|__   ___ | |___ 
      / / ____| |   / _` | '_ \    | |/ _ \ / _ \| / __|
     / / |____| |__| (_| | |_) |   | | (_) | (_) | \__ \
    /____|    |_____\__,_|_.__/    |_|\___/ \___/|_|___/
                                                        
    """

    tool_name = f"{tool} ({version})\n{tool_description}"

    output = f"{'#'*80}\n{logo}\n{tool_name}\n\n{'#'*80}\n"

    print(output)


def load_model(model_path: str) -> cobra.Model:
    """
    Load a multi-species model into memory given a path to a model file in a
    COBRApy supported format.

    Parameters
    ----------
    model_path : str
        Path to a multi-species model in any COBRApy supported format.

    Returns
    -------
    cobra.Model
        A COBRApy model loaded into memory.

    Raises
    ------
    ValueError
        If the model_path does not exist.
    ValueError
        If the model format is not supported.
    """

    if not os.path.exists(model_path):
        raise ValueError(f"Model path does not exist: {model_path}")

    print(f"\n[START] Loading model from {model_path}...")

    if model_path.endswith(".xml"):
        model = cobra.io.read_sbml_model(model_path)
    elif model_path.endswith(".json"):
        model = cobra.io.load_json_model(model_path)
    elif model_path.endswith(".yml"):
        model = cobra.io.load_yaml_model(model_path)
    elif model_path.endswith(".mat"):
        model = cobra.io.load_matlab_model(model_path)
    elif model_path.endswith(".sbml"):
        model = cobra.io.read_sbml_model(model_path)
    else:
        raise ValueError(f"Model format not supported for {model_path}")

    model.name = os.path.basename(model_path).split(".")[0]

    print(f"\n[DONE] {model.name} loaded.")

    return model


def set_default_bounds(model: cobra.Model) -> bool:
    """
    Set the bounds of the model's reactions according to default conventions;
    prints the changes and returns True if the bounds were different from the
    default state. Conventions are based on Heinken et al. (2022), mgPipe models.

    Parameters
    ----------
    model : cobra.Model
        The model whose reactions' bounds are to be set.

    Returns
    -------
    bool
        True if the bounds were different from the default state.

    Notes
    -----
    The conventions are as follows:
    1. Set the bounds of the fecal exchange (EX_met[fe]) reactions for metabolites to be (-1000., 1000000.)
    2. Set the bounds of the fecal exchange (EX_met[fe]) reaction for "microbeBiomass" to be (-10000., 1000000.)
    3. Set the bounds of the fecal transport (UFEt_met) reactions to be (0., 1000000.)
    4. Set the bounds of the microbe secretion/uptake (microbe_IEX_met[u]tr) reactions to be (-1000., 1000.)
    5. Set the bounds of the community biomass reaction to be (0.4, 1.)
    """

    print("\n[START] Setting default bounds...")

    saved_bounds = dict()
    new_bounds = dict()
    for rxn in model.reactions:
        # Set the bounds of the fecal exchange (EX_met[fe]) reactions to be (-1000., 1000000.)
        if (
            rxn.id.startswith("EX_")
            and rxn.id.endswith("[fe]")
            and "microbeBiomass" not in rxn.id
        ):
            saved_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
            model.reactions.get_by_id(rxn.id).bounds = (-1000.0, 1000000.0)
            new_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
        # Set the bounds of the fecal exchange (EX_met[fe]) reactions for the microbeBiomass to be (-10000., 1000000.)
        elif (
            rxn.id.startswith("EX_")
            and rxn.id.endswith("[fe]")
            and "microbeBiomass" in rxn.id
        ):
            saved_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
            model.reactions.get_by_id(rxn.id).bounds = (-10000.0, 1000000.0)
            new_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
        # Set the bounds of the fecal transport (UFEt_met) reactions to be (0., 1000000.)
        elif rxn.id.startswith("UFEt_"):
            saved_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
            model.reactions.get_by_id(rxn.id).bounds = (0.0, 1000000.0)
            new_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
        # Set the bounds of the microbe secretion/uptake (microbe_IEX_met[u]tr) reactions to be (-1000., 1000.)
        elif "IEX" in rxn.id and rxn.id.endswith("[u]tr"):
            saved_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
            model.reactions.get_by_id(rxn.id).bounds = (-1000.0, 1000.0)
            new_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
        # Set the bounds of the diet transport (DUt_met) reactions to be (0., 1000000.)
        elif rxn.id.startswith("DUt_"):
            saved_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
            model.reactions.get_by_id(rxn.id).bounds = (0.0, 1000000.0)
            new_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
        # Set the bounds of the community biomass reaction to be (0.4, 1.)
        elif rxn.id == "communityBiomass":
            saved_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
            model.reactions.get_by_id(rxn.id).bounds = (0.4, 1.0)
            new_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds

    n_changed_bounds = 0
    # Print out the changes
    for rxn, bounds in saved_bounds.items():
        if bounds != new_bounds[rxn]:
            print(f"Changed bounds for {rxn} from {bounds} to {new_bounds[rxn]}")
            n_changed_bounds += 1

    if n_changed_bounds > 0:
        bounds_changed = True
        print(f"\n[DONE] Changed bounds for {n_changed_bounds} reactions.")
    else:
        bounds_changed = False
        print("\n[DONE] No bounds were changed.")

    return bounds_changed


def convert_model_format(model_path, output_path):
    """
    Convert a mgPipe.m (Heinken et al., 2022) matlab model to a json model.

    Parameters
    ----------
    model_path : str
        Path to the model file.
    output_path : str
        Path to the output file.

    Returns
    -------
    None

    Notes
    -----
    If the metabolite charge is NaN, it is converted to a string.
    """
    model = load_model(model_path)

    print(f"\n[START] Converting model {model.name} to json format...")

    # Convert the metabolite charge to a string if it is NaN
    for metab in model.metabolites:
        if isnan(metab.charge):
            metab.charge = "NaN"

    if not os.path.exists(output_path):
        os.makedirs(output_path)

    if output_path.endswith("/"):
        output_path = output_path[:-1]

    converted_output_filepath = f"{output_path}/{model.name}.json"

    cobra.io.save_json_model(model, converted_output_filepath)

    print(f"\n[DONE] {model.name} converted to json format.")


def fetch_norm_sample_metabolomics_data(
    model_filepath: str,
    gcms_filepath: str,
    match_key_output_filepath: str,
    manual_matching_filepath: str = "data_dependencies/manually_matched_keys.txt",
) -> dict:
    """
    Generate a dictionary of VMH IDs and their corresponding normalized sample-specific metabolite values.

    Parameters
    ----------
    model_filepath : str
        Filepath to the model.
    gcms_filepath : str
        Filepath to the GC-MS data.
    match_key_output_filepath : str
        Filepath to the directory where the matched key file will be saved.
    manual_matching_filepath : str
        Filepath to the manually matched key file.

    Returns
    -------
    vmh_id_values : dict
        Dictionary of VMH IDs and their corresponding normalized sample-specific metabolite values.
    """

    model = load_model(model_filepath)

    print(f"\n[START] Fetching metabolomics data for {model.name}...")

    # Read metabolomics data
    metabolomics_data = pd.read_csv(gcms_filepath, sep=",", index_col=0)
    if match_key_output_filepath not in os.listdir():
        os.mkdir(match_key_output_filepath)

    if match_key_output_filepath[-1] != "/":
        match_key_output_filepath += "/"

    match_names_to_vmh(
        gcms_filepath=gcms_filepath,
        output_filepath=match_key_output_filepath,
        manual_matching_filepath=manual_matching_filepath,
    )

    matched_metabolite_names = {}
    with open(
        f"{match_key_output_filepath}{gcms_filepath.split('/')[-1].split('.')[-2]}_matched_key.txt",
        "r",
    ) as f:
        matches = f.readlines()
        if matches != "":
            matches = [match.strip().split("\t") for match in matches]
            for match in matches:
                matched_metabolite_names[match[0]] = match[1]

    idx_list = []
    # Given the matched metabolite names, find the index where the metabolite is found in the metabolomics data
    for col_name in metabolomics_data.columns:
        if col_name in matched_metabolite_names.keys():
            idx_list.append(metabolomics_data.columns.get_loc(col_name))

    sample_id = []
    # Given the model name, find which row the sample ID is in
    for index in metabolomics_data.index:
        if index in model.name:
            sample_id.append(index)
        elif index in model.name:
            sample_id.append(index)

    if len(sample_id) > 1:
        print("Multiple sample IDs found in metabolomics data. Please check the data.")
        sys.exit(1)
    elif len(sample_id) == 0:
        print("No sample ID found in metabolomics data. Please check the data.")
        sys.exit(1)
    else:
        sample_id = sample_id[0]

    # Given a sample ID, get the metabolomics data for that sample
    sample_metabolomic_data = metabolomics_data.loc[sample_id][min(idx_list) :]

    # Create a dictionary of metabolite names and their concentrations
    metabolite_raw_vals_dict = {
        name: float(conc) for name, conc in sample_metabolomic_data.items()
    }

    vmh_id_values = {}
    for vmh_name, vmh_id in matched_metabolite_names.items():
        for gcms_name, value in metabolite_raw_vals_dict.items():
            if vmh_name == gcms_name:
                vmh_id_values[vmh_id] = value

    # Normalize the values
    total = sum(vmh_id_values.values())
    normalized_vmh_id_values = {k: v / total for k, v in vmh_id_values.items()}

    assert sum(normalized_vmh_id_values.values()) == 1.0

    print(
        f"\nNumber of VMH ID-matched metabolites: {len(normalized_vmh_id_values)} of {len(metabolite_raw_vals_dict)}"
    )

    print(
        f"\n[DONE] Returning normalized sample-specific metabolomics values for {sample_id}."
    )

    return normalized_vmh_id_values
