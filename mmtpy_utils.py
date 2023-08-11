import os
import re
import sys
import time
from math import isnan

import cobra
import numpy as np
import pandas as pd
import pubchempy as pcp


def print_logo(tool: str, tool_description: str, version: str) -> None:
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


def load_model(model_path: str, simple_model_name: bool = True) -> cobra.Model:
    """
    Load a multi-species model into memory given a path to a model file in a COBRApy supported format.

    Parameters
    ----------
    model_path : str
        Path to a multi-species model in any COBRApy supported format.
    simple_model_name : bool
        If True, the model name will be set to the file name with only the sample identifier (e.g., Case_1) and not the full file name (e.g., Case_1.xml).

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

    print(f"\n[Loading model from {model_path}]")

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

    # Set the model name
    if simple_model_name:
        model_file_name = os.path.basename(model_path).split(".")[0]
        model_name_search_pattern = re.compile(r"(Case|Control)_\d+")
        model.name = model_name_search_pattern.search(model_file_name).group(0)
    else:
        model.name = os.path.basename(model_path).split(".")[0]

    print(f"\n[{model.name} loaded]")

    return model


def set_default_bounds(
    model: cobra.Model,
    source: str = "MMTpy",
    rxn_type: str = "all",
    silent: bool = False,
) -> bool:
    """
    Set the bounds of the model's reactions according to conventions;
    prints the changes and returns True if the bounds were different from the default state.
    Conventional bounds can either be set based on Heinken et al. (2022), mgPipe models or be set based on MMTpy conventions.

    Parameters
    ----------
    model : cobra.Model
        The model whose reactions' bounds are to be set.
    source : str, optional
        The definition of conventional bounds, by default "MMTpy"; options are of either "MMTpy" or "MATLAB".
    rxn_type : str, optional
        The type of reactions whose bounds are to be set, by default "all";
        options are of either "all", "FEX", "UFEt", "IEX", "DUt", or "commBiomass".
    silent : bool, optional
        Whether to print the changes, by default False.

    Returns
    -------
    bool
        True if the bounds were different from the default state.

    Notes
    -----
    The conventions are as follows based on Heinken et al. (2022), mgPipe models:
    1. Set the bounds of the fecal exchange (EX_met[fe]) reactions for metabolites to be (-1000., 1000000.)
    2. Set the bounds of the fecal exchange (EX_met[fe]) reaction for "microbeBiomass" to be (-10000., 1000000.)
    3. Set the bounds of the fecal transport (UFEt_met) reactions to be (0., 1000000.)
    4. Set the bounds of the microbe secretion/uptake (microbe_IEX_met[u]tr) reactions to be (-1000., 1000.)
    5. Set the bounds of the community biomass reaction to be (0.4, 1.)

    The default bounds for MMTpy are as follows:
    1. Set the bounds of the fecal exchange (EX_met[fe]) reactions for metabolites to be (0., 1000000.)
    2. Set the bounds of the fecal exchange (EX_met[fe]) reaction for "microbeBiomass" to be (0., 1000000.)
    3. Set the bounds of the fecal transport (UFEt_met) reactions to be (0., 1000000.)
    4. Set the bounds of the microbe secretion/uptake (microbe_IEX_met[u]tr) reactions to be (-1000., 1000.)
    5. Set the bounds of the community biomass reaction to be (0.4, 1.)
    """

    print(f"\n[Setting bounds for {rxn_type} reactions using {source} defaults]")

    saved_bounds = dict()
    new_bounds = dict()
    for rxn in model.reactions:
        # Set the bounds of the fecal exchange (EX_met[fe]) reactions to be (-1000., 1000000.)
        if (
            rxn.id.startswith("EX_")
            and rxn.id.endswith("[fe]")
            and "microbeBiomass" not in rxn.id
            and rxn_type in ["all", "FEX"]
        ):
            saved_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
            if source == "MMTpy":
                model.reactions.get_by_id(rxn.id).bounds = (0.0, 1000000.0)
            elif source == "MATLAB":
                model.reactions.get_by_id(rxn.id).bounds = (-1000.0, 1000000.0)
            else:
                raise ValueError(f"Source {source} not supported")
            new_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
        # Set the bounds of the fecal exchange (EX_met[fe]) reactions for the microbeBiomass to be (-10000., 1000000.)
        elif (
            rxn.id.startswith("EX_")
            and rxn.id.endswith("[fe]")
            and "microbeBiomass" in rxn.id
            and rxn_type in ["all", "FEX"]
        ):
            saved_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
            if source == "MMTpy":
                model.reactions.get_by_id(rxn.id).bounds = (0.0, 1000000.0)
            elif source == "MATLAB":
                model.reactions.get_by_id(rxn.id).bounds = (-10000.0, 1000000.0)
            else:
                raise ValueError(f"Source {source} not supported")
            new_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
        # Set the bounds of the fecal transport (UFEt_met) reactions to be (0., 1000000.)
        elif rxn.id.startswith("UFEt_") and rxn_type in ["all", "UFEt"]:
            saved_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
            model.reactions.get_by_id(rxn.id).bounds = (0.0, 1000000.0)
            new_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
        # Set the bounds of the microbe secretion/uptake (microbe_IEX_met[u]tr) reactions to be (-1000., 1000.)
        elif (
            "IEX" in rxn.id and rxn.id.endswith("[u]tr") and rxn_type in ["all", "IEX"]
        ):
            saved_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
            model.reactions.get_by_id(rxn.id).bounds = (-1000.0, 1000.0)
            new_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
        # Set the bounds of the diet transport (DUt_met) reactions to be (0., 1000000.)
        elif rxn.id.startswith("DUt_") and rxn_type in ["all", "DUt"]:
            saved_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
            model.reactions.get_by_id(rxn.id).bounds = (0.0, 1000000.0)
            new_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
        # Set the bounds of the community biomass reaction to be (0.4, 1.)
        elif rxn.id == "communityBiomass" and rxn_type in ["all", "commBiomass"]:
            saved_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
            model.reactions.get_by_id(rxn.id).bounds = (0.4, 1.0)
            new_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds

    n_changed_bounds = 0
    # Print out the changes
    for rxn, bounds in saved_bounds.items():
        if bounds != new_bounds[rxn]:
            if not silent:
                print(f"\tChanged bounds for {rxn} from {bounds} to {new_bounds[rxn]}")
            n_changed_bounds += 1

    if n_changed_bounds == 0:
        bounds_changed = False
        print("\tNo bounds were changed")
    else:
        bounds_changed = True
        print(
            f"\tChanged {n_changed_bounds}/{len(saved_bounds)} {rxn_type} reaction bounds for {model.name}"
        )

    return bounds_changed


def set_ba_diet_bounds(
    model: cobra.Model,
    diet_filepath: str = "data_dependencies/Heinken_2019_BA_diet.csv",
    silent: bool = False,
) -> None:
    """
    Changes select diet bounds according to Heinken et al. (2019) BA diet and tests model feasibility.
    Bounds are changed either if the BA diet lower bound is less than or equal to the model lower bound.

    Parameters
    ----------
    model : cobra.Model
        Model to be modified.
    diet_filepath : str, optional
        Path to csv file containing BA diet bounds. The default is "data_dependencies/Heinken_2019_BA_diet.csv".
    silent : bool, optional
        If True, no changes are printed. The default is False.

    Returns
    -------
    None
    """
    ba_diet_df = pd.read_csv(diet_filepath, skiprows=1)

    ba_diet_df["Reaction"] = ba_diet_df["Reaction"].str.replace("(e)", "", regex=False)

    ba_diet_df["Upper bound"] = ba_diet_df["Upper bound"].replace(1000.0, 10000.0)

    ba_diet_dict = {
        f"Diet_{ba_diet_df['Reaction'][i]}[d]": (
            ba_diet_df["Lower bound"][i],
            ba_diet_df["Upper bound"][i],
        )
        for i in range(len(ba_diet_df))
        if "EX_biomass[c]" not in ba_diet_df["Reaction"][i]
    }

    print(
        "\n[1/2] Changing select diet bounds according to Heinken et al. (2019) BA diet"
    )
    for diet_rxn, diet_bounds in ba_diet_dict.items():
        rxn_bounds = [
            (-0.1, 10000.0),
            (0.0, 10000.0),
            (-1000.0, 10000.0),
            (-50.0, 10000.0),
        ]
        if diet_rxn in model.reactions:
            if (
                diet_bounds in rxn_bounds
                and diet_bounds != model.reactions.get_by_id(diet_rxn).bounds
                and diet_bounds[0] <= model.reactions.get_by_id(diet_rxn).lower_bound
            ):
                if not silent:
                    print(
                        f"\t{diet_rxn}:\t bounds changed from {model.reactions.get_by_id(diet_rxn).bounds} to {ba_diet_dict[diet_rxn]}"
                    )
                model.reactions.get_by_id(diet_rxn).bounds = ba_diet_dict[diet_rxn]
        elif diet_rxn not in model.reactions:
            print(f"\tWarning: {diet_rxn}:\treaction not in model")
        else:
            print(f"\tWarning: {diet_rxn}:\tbounds not changed")

    print("\n[2/2] Testing model feasibility with BA diet")
    model.objective = 0
    ba_diet_sol = model.optimize()
    if ba_diet_sol.status == "optimal":
        print(
            "\tModel is feasible with relaxed BA diet, proceeding with remaining steps"
        )
    else:
        print("\tWarning: Model is infeasible with relaxed BA diet, terminating workflow")
        sys.exit()


def convert_model_format(
    model_path: str or cobra.Model, output_path: str = None
) -> None:
    """
    Convert a mgPipe.m (Heinken et al., 2022) MATLAB model to a json model.

    Parameters
    ----------
    model_path : str or cobra.Model
        Path to the model file or a COBRApy model loaded into memory.
    output_path : str
        Path to the output file.

    Returns
    -------
    None

    Raises
    ------
    ValueError
        If model_path is not a string or a COBRApy model.

    Notes
    -----
    If the metabolite charge is NaN, it is converted to a string.
    """
    if isinstance(model_path, str):
        model = load_model(model_path=model_path, simple_model_name=False)
    elif isinstance(model_path, cobra.Model):
        model = model_path
    else:
        raise ValueError(
            f"model_path must be a string or a COBRApy model, not {type(model_path)}"
        )

    print(f"\n[Converting model {model.name} to json format]")

    # Convert the metabolite charge to a string if it is NaN
    for metab in model.metabolites:
        if isnan(metab.charge):
            metab.charge = "nan"

    if not os.path.exists(output_path):
        os.makedirs(output_path)

    if output_path.endswith("/"):
        converted_output_filepath = f"{output_path}{model.name}.json"
    else:
        converted_output_filepath = f"{output_path}/{model.name}.json"

    cobra.io.save_json_model(model, converted_output_filepath)

    print(f"\n[{model.name} converted to json format]")


def convert_string(s: str) -> str:
    """
    Convert a string to a standard format for matching.

    Parameters
    ----------
    s : str
        String to be converted.

    Returns
    -------
    str
        Converted string.

    Notes
    -----
    The following operations are performed:
        1. Remove everything in parenthesis.
        2. Convert '.' to ',' between numbers.
        3. Add a hyphen between number and word if they are separated by space.
    """
    # Remove everything in parenthesis
    s = re.sub(r"\s*\(.*?\)", "", s)

    # Convert '.' to ',' between numbers
    s = re.sub(r"(\d)\.(\d)", r"\1,\2", s)

    # Add a hyphen between number and word if they are separated by space
    s = re.sub(r"(\d) (\w)", r"\1-\2", s)

    return s


def get_init_mbx_idx(df: pd.DataFrame) -> int:
    """
    Get the index of the first numerical column in a dataframe.

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe to be searched.

    Raises
    ------
    ValueError
        If no numerical columns are found.

    Returns
    -------
    int
        Index of the first numerical column.
    """
    for index, (col, dtype) in enumerate(df.dtypes.items()):
        if np.issubdtype(dtype, np.number):
            return index
    raise ValueError("No numerical columns found.")


def match_names_to_vmh(
    mbx_filepath: str,
    output_filepath: str,
    reuturn_matched_keys: bool,
    vmh_db_filepath: str = "data_dependencies/all_vmh_metabolites.tsv",
    manual_matching_filepath: str = "data_dependencies/manually_matched_keys.txt",
    silent: bool = False,
) -> dict:
    """
    Map the metabolite names detected by MBX to VMH identifiers for a given
    MBX dataset. The matching is performed in the following order:
        1. Direct matching of MBX names to VMH identifiers using the VMH database.
        2. Matching of MBX names to VMH identifiers via pubchempy.
        3. Manual matching of MBX names to VMH identifiers.

    Parameters
    ----------
    mbx_filepath : str
        Filepath to the MBX data.
    output_filepath : str
        Filepath (including .txt file name) for saving the matching keys.
    reuturn_matched_keys : bool
        Whether to return the matched keys as a dictionary.
    vmh_db_filepath : str
        Filepath to the VMH database of metabolites and their identifiers.
    manual_matching_filepath : str
        Filepath to the manually matched keys.
    silent : bool
        If True, no matchings are printed.

    Raises
    ------
    ValueError
        If model_path is not a string or a COBRApy model.

    Returns
    -------
    dict
        Dictionary of matched keys.

    Notes
    -----
    The metabolomics data must be a .csv file with samples as the column headers
    and metabolite names as the row indicies. The metabolite names must be canonical
    metabolite names (e.g. "L-Aspartate", "L-Asparagine", etc.). For matching via
    pubchempy, internet access is required; otherwise, the matching will fallback
    on direct and manual matching.
    """

    # Load the data
    vmh_db_df = pd.read_csv(vmh_db_filepath, index_col=0, header=0, delimiter="\t")
    mbx_data_df = pd.read_csv(mbx_filepath, index_col=0, header=0)

    print("\n[Matching MBX names to VMH identifiers]")

    print(
        "\n\t[1/3] Direct matching of MBX names to VMH identifiers using the VMH database"
    )

    # Get the index of the first numerical column in mbx_data_df
    first_mbx_idx = get_init_mbx_idx(mbx_data_df)

    # Create dictionaries for direct matching
    mbx_names_dict = {
        name: convert_string(name).lower()
        for name in mbx_data_df.columns[first_mbx_idx:].to_list()
    }
    vmh_names_dict = dict(zip(vmh_db_df["fullName"].index, vmh_db_df["fullName"]))

    # Perform direct matching
    direct_matching_dict = dict()
    for vmh_id, vmh_name in vmh_names_dict.items():
        for mbx_name, mbx_alt_name in mbx_names_dict.items():
            if mbx_alt_name.lower() == vmh_name.lower():
                direct_matching_dict[vmh_id] = mbx_name

    # Create a dict with key value pairs that remain unmatched for mbx_names_dict
    unmatched_dict = {
        vmh_id: name
        for vmh_id, name in mbx_names_dict.items()
        if name not in direct_matching_dict.values()
    }

    print("\n\t[2/3] Matching of MBX names to VMH identifiers via PubChemPy")
    # Match by pubchempy and vmh database
    # NOTE {vmh_id, matched_identifier}
    vmh_cid_dict = dict(zip(vmh_db_df["pubChemId"].index, vmh_db_df["pubChemId"]))
    vmh_inchikey_dict = dict(zip(vmh_db_df["inchiKey"].index, vmh_db_df["inchiKey"]))
    vmh_inchistring_dict = dict(
        zip(vmh_db_df["inchiString"].index, vmh_db_df["inchiString"])
    )
    vmh_smiles_dict = dict(zip(vmh_db_df["smile"].index, vmh_db_df["smile"]))

    # Empty dictionary to store standardized names
    # NOTE {MBX name, matched_identifier}
    iupac_names_dict = dict()
    cid_names_dict = dict()
    inchi_names_dict = dict()
    inchikey_names_dict = dict()
    smiles_names_dict = dict()

    # Iterate over each item in the compounds dictionary
    for mbx_name, compound in unmatched_dict.items():
        try:
            # Get the compound information from PubChem
            c = pcp.get_compounds(compound, "name")
            # If the compound was found, store its properties
            if c:
                if c[0].iupac_name in vmh_names_dict.values():
                    iupac_names_dict[mbx_name] = c[0].iupac_name
                if c[0].cid in vmh_cid_dict.values():
                    cid_names_dict[mbx_name] = int(c[0].cid)
                if c[0].inchi in vmh_inchistring_dict.values():
                    inchi_names_dict[mbx_name] = c[0].inchi
                if c[0].inchikey in vmh_inchikey_dict.values():
                    inchikey_names_dict[mbx_name] = c[0].inchikey
                if c[0].isomeric_smiles in vmh_smiles_dict.values():
                    smiles_names_dict[mbx_name] = c[0].isomeric_smiles
        except Exception as e:
            print(f"\t\tError getting info for '{compound}': {e}")
        time.sleep(0.5)

    pubchempy_matched_dict = dict()

    for vmh_id, vmh_inchikey in vmh_inchikey_dict.items():
        for mbx_name, pubchempy_inchikey in inchikey_names_dict.items():
            if vmh_inchikey != "nan" and vmh_inchikey == pubchempy_inchikey:
                if not silent:
                    print(f"\t\tMatched '{mbx_name}' to '{vmh_id}' using InChIKey")
                pubchempy_matched_dict[mbx_name] = vmh_id
    for vmh_id, vmh_cid in vmh_cid_dict.items():
        for mbx_name, pubchempy_cid in cid_names_dict.items():
            if vmh_cid != "nan" and vmh_cid == pubchempy_cid:
                if not silent:
                    print(f"\t\tMatched '{mbx_name}' to '{vmh_id}' using CID")
                pubchempy_matched_dict[mbx_name] = vmh_id
    for vmh_id, vmh_inchi in vmh_inchistring_dict.items():
        for mbx_name, pubchempy_inchi in inchi_names_dict.items():
            if vmh_inchi != "nan" and vmh_inchi == pubchempy_inchi:
                if not silent:
                    print(f"\t\tMatched '{mbx_name}' to '{vmh_id}' using InChI")
                pubchempy_matched_dict[mbx_name] = vmh_id

    # Combine the direct matching dictionary with the pubchempy matched dictionary
    pubchempy_matched_dict.update(
        {value: key for key, value in direct_matching_dict.items()}
    )

    if len(pubchempy_matched_dict) != len(mbx_data_df.index):
        for vmh_id, vmh_smiles in vmh_smiles_dict.items():
            for mbx_name, pubchempy_smiles in smiles_names_dict.items():
                if vmh_smiles != "nan" and vmh_smiles == pubchempy_smiles:
                    if not silent:
                        print(f"\t\tMatched '{mbx_name}' to '{vmh_id}' using SMILES")
                    pubchempy_matched_dict[mbx_name] = vmh_id

    print(
        "\n\t[3/3] Matching of MBX names to VMH identifiers via manual matching database"
    )
    manual_matching = dict()
    with open(manual_matching_filepath, "r") as f:
        for line in f:
            name, vmh = line.split("\t")
            if vmh.endswith("\n"):
                vmh = vmh[:-1]
            manual_matching[name] = vmh

    max_matched_dict = dict()
    for name, id in manual_matching.items():
        if id not in pubchempy_matched_dict.values():
            max_matched_dict[name] = id

    max_matched_dict.update(pubchempy_matched_dict)

    # If the output filepath does not exist, create it
    if not os.path.exists(output_filepath):
        os.makedirs(output_filepath)

    if output_filepath[-1] != "/":
        output_filepath += "/"

    key_output_filepath = (
        f"{output_filepath}{mbx_filepath.split('/')[-1].split('.')[-2]}_matched_key.txt"
    )

    # Write out the matched identifiers to a .txt file
    with open(key_output_filepath, "w") as f:
        for key, value in max_matched_dict.items():
            f.write(f"{key}:\t{value}\n")
    print(
        f"\t\t{len(max_matched_dict)} of {len(mbx_data_df.columns)-2} VMH identifiers matched to the MBX metabolite names"
    )

    print(
        f"\n\t[Matched MBX names to VMH identifiers and written to '{key_output_filepath}']"
    )

    if reuturn_matched_keys:
        return max_matched_dict


def fetch_norm_sample_mbx_data(
    model_input: cobra.Model or str,
    mbx_filepath: str,
    matched_mbx_names: dict,
) -> dict:
    """
    Generate a dictionary of VMH IDs and their corresponding normalized sample-specific metabolite values.

    Parameters
    ----------
    model_input : cobra.Model or str
        The COBRApy model loaded into memory or a file path to the model
    mbx_filepath : str
        Filepath to the MBX data.
    match_key_input : dict
        Dictionary of matched vmh_id and metabolite name keys.

    Returns
    -------
    dict
        Dictionary of VMH IDs and their corresponding normalized sample-specific metabolite values.

    Raises
    ------
    TypeError
        If model_input is not a cobra.Model or a filepath to a COBRApy model.
    ValueError
        If the model has multiple or no metabolomics data attributes.
    ValueError
        If the model mbx data values are not numeric after setting to float.
    ValueError
        If the matched_mbx_names dictionary is empty.
    """
    if isinstance(model_input, str):
        model = load_model(model_input)
    elif isinstance(model_input, cobra.Model):
        model = model_input
    else:
        raise TypeError(
            "model_input must be a cobra.Model or a filepath to a COBRApy model"
        )

    print(f"\n[Fetching MBX data for {model.name}]")

    # Read metabolomics data
    mbx_data = pd.read_csv(mbx_filepath, sep=",", index_col=0)

    idx_list = []
    # Given the matched metabolite names, find the index where the metabolite is found in the MBX data
    for col_name in mbx_data.columns:
        if col_name in matched_mbx_names.keys():
            idx_list.append(mbx_data.columns.get_loc(col_name))

    sample_id = []
    # Given the model name, find which row the sample ID is in
    for index in mbx_data.index:
        if index == model.name:
            sample_id.append(index)

    if len(sample_id) > 1:
        ValueError("Multiple sample IDs found in MBX data. Please check the data.")
    elif len(sample_id) == 0:
        ValueError("No sample ID found in MBX data. Please check the data.")
    else:
        sample_id = sample_id[0]

    # Given a sample ID, get the MBX data for that sample
    sample_mbx_data = mbx_data.loc[sample_id][min(idx_list) :]

    # Create a dictionary of metabolite names and their concentrations
    metab_raw_vals_dict = {name: conc for name, conc in sample_mbx_data.items()}

    # If the values of metab_raw_vals_dict are strings, remove any commas and convert them to floats
    for k, v in metab_raw_vals_dict.items():
        if isinstance(v, str):
            metab_raw_vals_dict[k] = float(v.replace(",", ""))

    # Check if all values in metab_raw_vals_dict are floats
    if not all(isinstance(v, float) for v in metab_raw_vals_dict.values()):
        raise ValueError("Not all values in metab_raw_vals_dict are floats.")

    # Raise ValueError if matched_mbx_names is empty
    if len(matched_mbx_names) == 0:
        raise ValueError("The VMH Identifier name-matching dictionary is empty.")

    # Create a dictionary of VMH IDs and their corresponding metabolite values if the metabolite is in the model
    vmh_id_values = dict()
    for vmh_name, vmh_id in matched_mbx_names.items():
        for mbx_name, value in metab_raw_vals_dict.items():
            if vmh_name == mbx_name and f"EX_{vmh_id}[fe]" in [
                rxn.id for rxn in model.reactions
            ]:
                vmh_id_values[vmh_id] = value

    # Normalize the values
    total_val = sum(vmh_id_values.values())
    norm_vmh_id_vals = {k: v / total_val for k, v in vmh_id_values.items()}

    print(
        f"\n\t[Number of name-matched and normalized metabolites in the model: {len(norm_vmh_id_vals)}/{len(metab_raw_vals_dict)} ({round(len(norm_vmh_id_vals)/len(metab_raw_vals_dict)*100, 2)}%)]"
    )
    print(f"\n[Returning normalized sample-specific MBX values for {sample_id}]")

    return norm_vmh_id_vals


def fetch_mbx_constr_list(model: cobra.Model, mbx_metab_norm_dict: dict) -> list:
    """
    Compute the MBX associated FEX reaction constraints for a sample and tests
    if the addition of the constraints gives a feasible solution.

    Parameters
    ----------
    model : cobra.Model
        The model to be constrained.
    mbx_metab_norm_dict : dict
        A dictionary of metabolite names and their normalized values.

    Returns
    -------
    list
        A list of the constraints added to the model.

    Raises
    ------
    ValueError
        If no constraints were added to the model.
    """
    print(f"\n[Fetching MBX constraints for {model.name}]")
    # Calculate the constraint for each metabolite's fecal exchange reactions
    metab_constrained_flux_expr = dict()
    for vmh_id, mbx_value in mbx_metab_norm_dict.items():
        if mbx_value != 0.0 and f"EX_{vmh_id}[fe]" in model.reactions:
            flux_expr = model.reactions.get_by_id(f"EX_{vmh_id}[fe]").flux_expression
            metab_constrained_flux_expr[
                f"EX_{vmh_id}[fe]"
            ] = flux_expr - mbx_value * sum(
                model.reactions.get_by_id(f"EX_{metab}[fe]").flux_expression
                for metab in mbx_metab_norm_dict.keys()
                if f"EX_{metab}[fe]" in model.reactions
            )

    # Add the constraints to a list for each metabolite with MBX data
    constraint_list = []
    for ex_rxn_id, constr_expr in metab_constrained_flux_expr.items():
        constraint = model.problem.Constraint(
            constr_expr, lb=0, ub=0, name=ex_rxn_id + "_constraint"
        )
        constraint_list.append(constraint)

    if len(constraint_list) == 0:
        raise ValueError("No constraints were added to the model.")

    print(
        f"\n\tAdding {len(constraint_list)} computed MBX constraints to {model.name} and testing if the solution is feasible"
    )
    # Add the constraints to the model
    model.add_cons_vars(constraint_list)
    model.solver.update()

    # Test if the solution is feasible
    model.objective = 0
    solution = model.optimize()

    if solution.status == "infeasible":
        model.remove_cons_vars(constraint_list)
        model.solver.update()
        print(
            "\n\tWarning: the solution is infeasible by introducing the constraints without slack variables"
        )
    else:
        print("\n\tThe solution is feasible without adding slack variables")
        model.remove_cons_vars(constraint_list)
        model.solver.update()

    return constraint_list


def solve_mbx_constraints(
    model: cobra.Model, constraints: list, parallel: bool = False
) -> list:
    """
    Add slack variables to infeasible constraints and test if the solution is feasible.

    Parameters
    ----------
    model : cobra.Model
        The model to be constrained.
    constraints : list
        A list of constraints to be added to the model.
    parallel : bool, optional
        If True, no detailed outputs will be printed, by default False

    Returns
    -------
    list
        A list of the constraints added to the model.
    list
        A list of the slack variables added to the model.
    """
    slack_constraints = []
    slack_variables = []
    print(f"\n[Introducing slack variables to all constraints]")
    for constraint in constraints:
        # Introduce a slack variable to the constraint
        slack_variable_pos = model.problem.Variable(
            constraint.name + "_slack_pos", lb=0
        )
        slack_variable_neg = model.problem.Variable(
            constraint.name + "_slack_neg", lb=0
        )

        slack_variables.extend([slack_variable_pos, slack_variable_neg])

        # Create a new constraint that includes the slack variable
        new_constraint = model.problem.Constraint(
            constraint.expression + slack_variable_pos - slack_variable_neg,
            lb=constraint.lb,
            ub=constraint.ub,
            name=constraint.name + "_slack",
        )

        slack_constraints.append(new_constraint)

    # Add the new constraints to the model
    model.add_cons_vars(slack_constraints)
    model.solver.update()

    # Set the objective to be the sum of the absolute values of the slack variables
    model.objective = model.problem.Objective(sum(slack_variables))

    # Run optimization
    solution = model.optimize(objective_sense="minimize")

    # Check if the solution is feasible
    if solution.status == "optimal":
        print(
            "\n\tThe solution is feasible with slack variables added to all constraints"
        )
    elif solution.status == "infeasible":
        print(
            "\n\tWarning: The solution is infeasible with slack variables added to all constraints"
        )

    print("\n[Simplifying the constraints based on the slack variable primal values]")

    # After optimizing the model, check the optimal values of the slack variables and set the constraints accordingly
    slack_pos_pos_val_names = []
    slack_pos_zero_val_names = []
    slack_neg_pos_val_names = []
    slack_neg_zero_val_names = []
    if not parallel:
        print("\n\tSlack variables and primals:")
    for var in slack_variables:
        if not parallel:
            print(f"\t\t{var.name}:\t{var.primal}")
        if var.primal > 0 and var.name.endswith("_slack_pos"):
            slack_pos_pos_val_names.append(var.name.replace("_pos", ""))
        elif var.primal > 0 and var.name.endswith("_slack_neg"):
            slack_neg_pos_val_names.append(var.name.replace("_neg", ""))
        elif var.primal == 0 and var.name.endswith("_slack_pos"):
            slack_pos_zero_val_names.append(var.name.replace("_pos", ""))
        elif var.primal == 0 and var.name.endswith("_slack_neg"):
            slack_neg_zero_val_names.append(var.name.replace("_neg", ""))

    # Find the intersection of the slack variable names
    intersection_gt = list(
        set(slack_pos_zero_val_names).intersection(slack_neg_pos_val_names)
    )
    intersection_lt = list(
        set(slack_pos_pos_val_names).intersection(slack_neg_zero_val_names)
    )
    intersection_eq = list(
        set(slack_pos_zero_val_names).intersection(slack_neg_zero_val_names)
    )

    # Adjust constraints based on slack variable values
    refined_constraints = []
    constraint_details_log = []
    if not parallel:
        print("\n\tSolved constraint details:")
    for constraint in slack_constraints:
        if constraint.name in intersection_gt:
            # Adjust constraint to v_EX_fecal_metab_i > x * summation_j_in_set_MBX_data(v_EX_fecal_metab_j)
            constraint.lb = 0
            constraint.ub = None
            refined_constraints.append(constraint)
            constraint_details = f"\t\t{constraint.name}:\t({constraint.lb}, {constraint.ub}),\t[Slacked, flux greater than constrained value]"
            constraint_details_log.append(constraint_details)
            if not parallel:
                print(constraint_details)
        elif constraint.name in intersection_lt:
            # Adjust constraint to v_EX_fecal_metab_i < x * summation_j_in_set_MBX_data(v_EX_fecal_metab_j)
            constraint.lb = None
            constraint.ub = 0
            refined_constraints.append(constraint)
            constraint_details = f"\t\t{constraint.name}:\t({constraint.lb}, {constraint.ub}),\t[Slacked, flux less than constrained value]"
            constraint_details_log.append(constraint_details)
            if not parallel:
                print(constraint_details)
        elif constraint.name in intersection_eq:
            # Adjust constraint to v_EX_fecal_metab_i = x * summation_j_in_set_MBX_data(v_EX_fecal_metab_j)
            constraint.lb = 0
            constraint.ub = 0
            refined_constraints.append(constraint)
            constraint_details = f"\t\t{constraint.name}:\t({constraint.lb}, {constraint.ub}),\t[Original, flux equal to constrained value]"
            constraint_details_log.append(constraint_details)
            if not parallel:
                print(constraint_details)

    # Remove the slack variables and constraints from the model
    model.remove_cons_vars(slack_variables)
    model.remove_cons_vars(slack_constraints)
    model.solver.update()

    print(f"\n[Returning solved constraints]")

    return refined_constraints, constraint_details_log
