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


def set_default_bounds(model: cobra.Model, source: str = "MMTpy") -> bool:
    """
    Set the bounds of the model's reactions according to conventions;
    prints the changes and returns True if the bounds were different from the
    default state.
    Conventional bounds can either be set based on Heinken et al. (2022), mgPipe
    models or be set based on MMTpy conventions.

    Parameters
    ----------
    model : cobra.Model
        The model whose reactions' bounds are to be set.
    source : str, optional
        The definition of conventional bounds, by default "MMTpy"; options are of
        either "MMTpy" or "MATLAB".

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

    print(f"\n[START] Setting default bounds using definitions based on {source}...")

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
            print(f"\tChanged bounds for {rxn} from {bounds} to {new_bounds[rxn]}")
            n_changed_bounds += 1

    if n_changed_bounds == 0:
        bounds_changed = False
        print("\n[DONE] No bounds were changed.")
    else:
        bounds_changed = True
        print(f"\n[DONE] Changed bounds for {n_changed_bounds} reactions.")

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
    mbx_filepath: str,
    output_filepath: str,
    vmh_db_filepath: str = "data_dependencies/all_vmh_metabolites.tsv",
    manual_matching_filepath: str = "data_dependencies/manually_matched_keys.txt",
    show_logo: bool = False,
) -> None:
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
    on direct and manual matching.
    """
    # Define tool metadata
    tool = "match-names-to-vmh"
    tool_description = "Matching of MBX metabolite names to VMH identifiers."

    if show_logo:
        print_logo(tool, tool_description, version)

    # Load the data
    vmh_db_df = pd.read_csv(vmh_db_filepath, index_col=0, header=0, delimiter="\t")
    mbx_data_df = pd.read_csv(mbx_filepath, index_col=0, header=0)

    print("\n[START] Matching MBX names to VMH identifiers...")

    print(
        "\n\t[1/3] Direct matching of MBX names to VMH identifiers using the VMH database."
    )
    # Create dictionaries for direct matching
    mbx_names_dict = {
        name: convert_string(name).lower() for name in mbx_data_df.columns[2:].to_list()
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

    print("\n\t[2/3] Matching of MBX names to VMH identifiers via PubChemPy.")
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

    pubchempy_matched_dict = dict()

    for vmh_id, vmh_inchikey in vmh_inchikey_dict.items():
        for mbx_name, pubchempy_inchikey in inchikey_names_dict.items():
            if vmh_inchikey != "nan" and vmh_inchikey == pubchempy_inchikey:
                print(f"\t\tMatched '{mbx_name}' to '{vmh_id}' using InChIKey")
                pubchempy_matched_dict[mbx_name] = vmh_id
    for vmh_id, vmh_cid in vmh_cid_dict.items():
        for mbx_name, pubchempy_cid in cid_names_dict.items():
            if vmh_cid != "nan" and vmh_cid == pubchempy_cid:
                print(f"\t\tMatched '{mbx_name}' to '{vmh_id}' using CID")
                pubchempy_matched_dict[mbx_name] = vmh_id
    for vmh_id, vmh_inchi in vmh_inchistring_dict.items():
        for mbx_name, pubchempy_inchi in inchi_names_dict.items():
            if vmh_inchi != "nan" and vmh_inchi == pubchempy_inchi:
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
                    print(f"\t\tMatched '{mbx_name}' to '{vmh_id}' using SMILES")
                    pubchempy_matched_dict[mbx_name] = vmh_id

    print(
        "\n\t[3/3] Matching of MBX names to VMH identifiers via manual matching database."
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
            f.write(f"{key}\t{value}\n")

    print(
        f"\n\t[DONE] Matched MBX names to VMH identifiers and written to '{key_output_filepath}'"
    )
    print(
        f"\t\t{len(max_matched_dict)} of {len(mbx_data_df.columns)-2} VMH identifiers matched to the MBX metabolite names."
    )


def fetch_norm_sample_mbx_data(
    model_input: cobra.Model or str,
    mbx_filepath: str,
    match_key_output_filepath: str,
    use_existing_matched_keys: bool = False,
    existing_keys_path: str = None,
    manual_matching_filepath: str = "data_dependencies/manually_matched_keys.txt",
    show_logo: bool = False,
) -> dict:
    """
    Generate a dictionary of VMH IDs and their corresponding normalized sample-specific metabolite values.

    Parameters
    ----------
    model_input : cobra.Model or str
        The COBRApy model loaded into memory or a file path to the model
    mbx_filepath : str
        Filepath to the MBX data.
    match_key_output_filepath : str
        Filepath to the directory where the matched key file will be saved. If the path is not supplied, the output will be saved in the working directory by default.
    use_existing_matched_keys : bool
        Whether to use existing matched keys from match_names_to_vmh().
    existing_keys_path : str (optional)
        If use_existing_matched_keys is true, load the keys; defaults to None.
    manual_matching_filepath : str (optional)
        Filepath to the manually matched key file; defaults to the work directory if a path is not supplied by the user.
    show_logo : bool (optional)
        Specification for printing the logo and function details.

    Returns
    -------
    dict
        Dictionary of VMH IDs and their corresponding normalized sample-specific metabolite values.

    Raises
    ------
    TypeError
        If model_input is not a cobra.Model or a filepath to a COBRApy model.
    ValueError
        If the model does not have a metabolomics data attribute.
    """
    tool = "fetch-norm-sample-mbx"
    tool_description = "Gets the normalized MBX data for a sample"

    if show_logo:
        print_logo(tool, tool_description, version)

    if type(model_input) == str:
        model = load_model(model_input)
    elif type(model_input) == cobra.Model:
        model = model_input
    else:
        raise TypeError(
            "model_input must be a cobra.Model or a filepath to a COBRApy model"
        )

    print(f"\n[START] Fetching MBX data for {model.name}...")

    # Read metabolomics data
    mbx_data = pd.read_csv(mbx_filepath, sep=",", index_col=0)

    if use_existing_matched_keys:
        if match_key_output_filepath == None:
            match_key_output_filepath = os.getcwd()
            print(
                f"\nNot using existing keys and path is not supplied; matched key file will be stored under the work directory: {match_key_output_filepath}"
            )
        else:
            match_key_output_filepath = existing_keys_path
    else:
        if match_key_output_filepath == None:
            match_key_output_filepath = os.getcwd()
            print(
                f"\nNot using existing keys and path is not supplied; matched key file will be stored under the work directory: {match_key_output_filepath}"
            )

        if match_key_output_filepath[-1] != "/":
            match_key_output_filepath += "/"

        match_names_to_vmh(
            mbx_filepath=mbx_filepath,
            output_filepath=match_key_output_filepath,
            manual_matching_filepath=manual_matching_filepath,
        )

    matched_mbx_names = dict()
    with open(
        f"{match_key_output_filepath}{mbx_filepath.split('/')[-1].split('.')[-2]}_matched_key.txt",
        "r",
    ) as f:
        matches = f.readlines()
        if matches != "":
            matches = [match.strip().split("\t") for match in matches]
            for match in matches:
                matched_mbx_names[match[0]] = match[1]

    idx_list = []
    # Given the matched metabolite names, find the index where the metabolite is found in the MBX data
    for col_name in mbx_data.columns:
        if col_name in matched_mbx_names.keys():
            idx_list.append(mbx_data.columns.get_loc(col_name))

    sample_id = []
    # Given the model name, find which row the sample ID is in
    for index in mbx_data.index:
        if index in model.name:
            sample_id.append(index)
        elif index in model.name:
            sample_id.append(index)

    if len(sample_id) > 1:
        print("Multiple sample IDs found in MBX data. Please check the data.")
        sys.exit(1)
    elif len(sample_id) == 0:
        print("No sample ID found in MBX data. Please check the data.")
        sys.exit(1)
    else:
        sample_id = sample_id[0]

    # Given a sample ID, get the MBX data for that sample
    sample_mbx_data = mbx_data.loc[sample_id][min(idx_list) :]

    # Create a dictionary of metabolite names and their concentrations
    metab_raw_vals_dict = {name: float(conc) for name, conc in sample_mbx_data.items()}

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

    print(f"\n[DONE] Returning normalized sample-specific MBX values for {sample_id}.")

    print(
        f"\tNumber of name-matched and normalized metabolites in the model: {len(norm_vmh_id_vals)}/{len(metab_raw_vals_dict)} ({round(len(norm_vmh_id_vals)/len(metab_raw_vals_dict)*100, 2)}%)"
    )

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
    """
    print(f"\n[START] Fetching MBX constraints for {model.name}...")
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

    print(
        f"\tTesting {len(constraint_list)} computed MBX constraints for {model.name}."
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
            "\n[WARNING] The solution is infeasible by introducing the constraints without slack variables."
        )
    else:
        print("\n[DONE] The solution is feasible without adding slack variables.")
        model.remove_cons_vars(constraint_list)
        model.solver.update()

    return constraint_list


def slack_constraints(model: cobra.Model, constraints: list) -> list:
    """
    Given a list of constraints, add slack variables to the constraints and test
    if the solution is feasible.

    Parameters
    ----------
    model : cobra.Model
        The model to be constrained.
    constraints : list
        A list of constraints to be added to the model.

    Returns
    -------
    list
        A list of the constraints added to the model.
    """
    # Add the constraints to the model and test for feasibility
    print("\n[START] Testing the feasibility of the model without slack variables...")
    model.add_cons_vars(constraints)
    model.solver.update()
    model.objective = 0
    solution = model.optimize()

    feasible_constraints = []

    if solution.status == "optimal":
        for constraint in constraints:
            print("\nFeasible model with constraints:")
            print(f"\t{constraint.name}:\t{constraint.primal}")
            feasible_constraints.append(constraint)
            print("\n[DONE] Returning feasible constraints.")
    else:
        # If the model is infeasible, loop through the constraints, remove them
        # from the model, and add a slack variable to a list.
        print("\n\tInfeasible model with constraints:")
        for constraint in constraints:
            print(f"\t{constraint.name}:\t{constraint.primal}")
        model.remove_cons_vars(constraints)
        model.solver.update()
        new_constraints = []

        print("\n\t[1/2] Adding slack variables to the constraints...")
        for constraint in constraints:
            # Introduce a slack variable to the constraint
            slack_variable_pos = model.problem.Variable(
                constraint.name + "_slack_pos", lb=0
            )
            slack_variable_neg = model.problem.Variable(
                constraint.name + "_slack_neg", lb=0
            )

            # Create a new constraint that includes the slack variable
            new_constraint = model.problem.Constraint(
                constraint.expression + slack_variable_pos - slack_variable_neg,
                lb=constraint.lb,
                ub=constraint.ub,
                name=constraint.name + "_slack",
            )

            # Add the new constraint to the list of new constraints
            new_constraints.append(new_constraint)

        # Add the new constraints to the model
        model.add_cons_vars(new_constraints)
        model.solver.update()

        print("\n\t[2/2] Testing the feasibility of the model with slack variables...")
        # Test for feasibility again
        model.objective = 0
        solution = model.optimize()
        if solution.status == "optimal":
            # If the solution is now feasible, add the new constraint to the list of feasible constraints
            print("\n\tFeasible model with slack variables:")
            for constraint in new_constraints:
                print(f"\t{constraint.name}:\t{constraint.primal}")
                feasible_constraints.append(constraint)
            # Remove the slack variables from the model
            model.remove_cons_vars(new_constraints)
            model.solver.update()
            print("\n[DONE] Returning feasible constraints.")
        else:
            # If the solution is still infeasible, print the values of the slack variables
            print("\n\tInfeasible model with slack variables:")
            for constraint in new_constraints:
                print(f"\t{constraint.name}:\t{constraint.primal}")
            # Remove the slack variables from the model
            model.remove_cons_vars(new_constraints)
            model.solver.update()
            print("\n[WARNING] Returning infeasible constraints.")

    return feasible_constraints
