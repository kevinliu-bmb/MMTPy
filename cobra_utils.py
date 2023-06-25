import os
import re
import sys
from math import isnan

import cobra
import pandas as pd
import pubchempy as pcp


version = "0.1.0"


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


def convert_model_format(model_path: str or cobra.Model, output_path: str = None):
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

    Notes
    -----
    If the metabolite charge is NaN, it is converted to a string.
    """
    if type(model_path) == str:
        model = load_model(model_path)
    else:
        model = model_path

    print(f"\n[START] Converting model {model.name} to json format...")

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

    print(f"\n[DONE] {model.name} converted to json format.")


def convert_string(s):
    """
    Convert a string to a standard format for matching.

    Parameters
    ----------
    s : str
        String to be converted.

    Returns
    -------
    s : str
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


def match_names_to_vmh(
    gcms_filepath: str,
    output_filepath: str,
    vmh_db_filepath: str = "data_dependencies/all_vmh_metabolites.tsv",
    manual_matching_filepath: str = "data_dependencies/manually_matched_keys.txt",
    show_logo: bool = False,
) -> None:
    """
    Map the metabolite names detected by GC-MS to VMH identifiers for a given
    GC-MS dataset. The matching is performed in the following order:
        1. Direct matching of GC-MS names to VMH identifiers.
        2. Matching of GC-MS names to VMH identifiers via pubchempy.
        3. Manual matching of GC-MS names to VMH identifiers.

    Parameters
    ----------
    gcms_filepath : str
        Filepath to the GC-MS data.
    output_filepath : str
        Filepath for saving the matching keys.
    vmh_db_filepath : str
        Filepath to the VMH database of metabolites and their identifiers.
    show_logo : bool
        Specification for printing the logo and function details.

    Returns
    -------
    None

    Notes
    -----
    The metabolomics data must be a .csv file with samples as the column headers
    and metabolite names as the row indicies. The metabolite names must be canonical
    metabolite names (e.g. "L-Aspartate", "L-Asparagine", etc.). For matching via
    pubchempy, internet access is required; otherwise, the matching will fallback
    on a combination of direct and manual matching.
    """
    # Define tool metadata
    tool = "match-names-to-vmh"
    tool_description = "Matching of GC-MS metabolite names to VMH identifiers."

    if show_logo:
        print_logo(tool, tool_description, version)

    # Load the data
    vmh_db = pd.read_csv(vmh_db_filepath, index_col=0, header=0, delimiter="\t")
    gcms_data = pd.read_csv(gcms_filepath, index_col=0, header=0)

    print("\n[START] Matching GC-MS names to VMH identifiers...")

    # Create dictionaries for direct matching
    gcms_names_dict = {
        name: convert_string(name).lower() for name in gcms_data.columns[2:].to_list()
    }
    vmh_names_dict = dict(zip(vmh_db["fullName"].index, vmh_db["fullName"]))

    # Perform direct matching
    direct_matching_dict = {}
    for vmh_id, vmh_name in vmh_names_dict.items():
        for gcms_name, gcms_alt_name in gcms_names_dict.items():
            if gcms_alt_name.lower() == vmh_name.lower():
                direct_matching_dict[vmh_id] = gcms_name

    # Create a dict with key value pairs that remain unmatched for gcms_names_dict
    unmatched_dict = {
        vmh_id: name
        for vmh_id, name in gcms_names_dict.items()
        if name not in direct_matching_dict.values()
    }

    # Match by pubchempy and vmh database
    # NOTE {vmh_id, matched_identifier}
    vmh_cid_dict = dict(zip(vmh_db["pubChemId"].index, vmh_db["pubChemId"]))
    vmh_inchikey_dict = dict(zip(vmh_db["inchiKey"].index, vmh_db["inchiKey"]))
    vmh_inchistring_dict = dict(zip(vmh_db["inchiString"].index, vmh_db["inchiString"]))
    vmh_smiles_dict = dict(zip(vmh_db["smile"].index, vmh_db["smile"]))

    # Empty dictionary to store standardized names
    # NOTE {GC-MS name, matched_identifier}
    iupac_names = {}
    cid_names = {}
    inchi_names = {}
    inchikey_names = {}
    smiles_names = {}

    # Iterate over each item in the compounds dictionary
    for gcms_name, compound in unmatched_dict.items():
        try:
            # Get the compound information from PubChem
            c = pcp.get_compounds(compound, "name")
            # If the compound was found, store its properties
            if c:
                if c[0].iupac_name in vmh_names_dict.values():
                    iupac_names[gcms_name] = c[0].iupac_name
                if c[0].cid in vmh_cid_dict.values():
                    cid_names[gcms_name] = int(c[0].cid)
                if c[0].inchi in vmh_inchistring_dict.values():
                    inchi_names[gcms_name] = c[0].inchi
                if c[0].inchikey in vmh_inchikey_dict.values():
                    inchikey_names[gcms_name] = c[0].inchikey
                if c[0].isomeric_smiles in vmh_smiles_dict.values():
                    smiles_names[gcms_name] = c[0].isomeric_smiles
        except Exception as e:
            print(f"\nError getting info for {compound}: {e}")

    pubchempy_matched_dict = {}

    for vmh_id, vmh_inchikey in vmh_inchikey_dict.items():
        for gcms_name, pubchempy_inchikey in inchikey_names.items():
            if vmh_inchikey != "nan" and vmh_inchikey == pubchempy_inchikey:
                print(f"\nMatched {gcms_name} to {vmh_id} using InChIKey")
                pubchempy_matched_dict[gcms_name] = vmh_id
    for vmh_id, vmh_cid in vmh_cid_dict.items():
        for gcms_name, pubchempy_cid in cid_names.items():
            if vmh_cid != "nan" and vmh_cid == pubchempy_cid:
                print(f"\nMatched {gcms_name} to {vmh_id} using CID")
                pubchempy_matched_dict[gcms_name] = vmh_id
    for vmh_id, vmh_inchi in vmh_inchistring_dict.items():
        for gcms_name, pubchempy_inchi in inchi_names.items():
            if vmh_inchi != "nan" and vmh_inchi == pubchempy_inchi:
                print(f"\nMatched {gcms_name} to {vmh_id} using InChI")
                pubchempy_matched_dict[gcms_name] = vmh_id

    # Combine the direct matching dictionary with the pubchempy matched dictionary
    pubchempy_matched_dict.update(
        {value: key for key, value in direct_matching_dict.items()}
    )

    if len(pubchempy_matched_dict) != len(gcms_data.index):
        for vmh_id, vmh_smiles in vmh_smiles_dict.items():
            for gcms_name, pubchempy_smiles in smiles_names.items():
                if vmh_smiles != "nan" and vmh_smiles == pubchempy_smiles:
                    print(f"\nMatched {gcms_name} to {vmh_id} using SMILES")
                    pubchempy_matched_dict[gcms_name] = vmh_id

    manual_matching = {}
    with open(manual_matching_filepath, "r") as f:
        for line in f:
            name, vmh = line.split("\t")
            if vmh.endswith("\n"):
                vmh = vmh[:-1]
            manual_matching[name] = vmh

    max_matched_dict = {}
    for name, id in manual_matching.items():
        if id not in pubchempy_matched_dict.values():
            max_matched_dict[name] = id

    max_matched_dict.update(pubchempy_matched_dict)

    print(
        f"\n{len(max_matched_dict)} of {len(gcms_data.columns)-2} VMH identifiers matched."
    )

    # If the output filepath does not exist, create it
    if not os.path.exists(output_filepath):
        os.makedirs(output_filepath)

    if output_filepath[-1] != "/":
        output_filepath += "/"

    key_output_filepath = f"{output_filepath}{gcms_filepath.split('/')[-1].split('.')[-2]}_matched_key.txt"

    # Write out the matched identifiers to a .txt file
    with open(key_output_filepath, "w") as f:
        for key, value in max_matched_dict.items():
            f.write(f"{key}\t{value}\n")

    print(
        f"\n[DONE] Matched GC-MS names to VMH identifiers and written to {key_output_filepath}"
    )


def fetch_norm_sample_metabolomics_data(
    model_input: cobra.Model or str,
    gcms_filepath: str,
    use_existing_matched_keys: bool = False,
    existing_keys_path: str = None,
    match_key_output_filepath: str = None,
    manual_matching_filepath: str = "data_dependencies/manually_matched_keys.txt",
    show_logo: bool = False,
) -> dict:
    """
    Generate a dictionary of VMH IDs and their corresponding normalized sample-specific metabolite values.

    Parameters
    ----------
    model_input : cobra.Model or str
        The COBRApy model loaded into memory or a file path to the model
    gcms_filepath : str
        Filepath to the GC-MS data.
    use_existing_matched_keys : bool
        Whether to use existing matched keys from match_names_to_vmh().
    existing_keys_path : str (optional)
        If use_existing_matched_keys is true, load the keys.
        Defaults to None.
    match_key_output_filepath : str (optional)
        Filepath to the directory where the matched key file will be saved.
        Defaults to None.
    manual_matching_filepath : str (optional)
        Filepath to the manually matched key file.
    show_logo : bool (optional)
        Specification for printing the logo and function details.

    Returns
    -------
    dict
        Dictionary of VMH IDs and their corresponding normalized sample-specific metabolite values.
    """
    tool = "fetch-norm-sample-metabomics"
    tool_description = "Gets the normalized metabolomics data for a sample"

    if show_logo:
        print_logo(tool, tool_description, version)

    if type(model_input) == str:
        model = load_model(model_input)

    print(f"\n[START] Fetching metabolomics data for {model.name}...")

    # Read metabolomics data
    metabolomics_data = pd.read_csv(gcms_filepath, sep=",", index_col=0)

    if use_existing_matched_keys:
        match_key_output_filepath = existing_keys_path
    else:
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
