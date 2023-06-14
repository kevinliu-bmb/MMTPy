import os
from math import isnan

import cobra
import pandas as pd




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
        if rxn.id.startswith("EX_") and rxn.id.endswith("[fe]") and "microbeBiomass" not in rxn.id:
            saved_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
            model.reactions.get_by_id(rxn.id).bounds = (-1000., 1000000.)
            new_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
        # Set the bounds of the fecal exchange (EX_met[fe]) reactions for the microbeBiomass to be (-10000., 1000000.)
        elif rxn.id.startswith("EX_") and rxn.id.endswith("[fe]") and "microbeBiomass" in rxn.id:
            saved_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
            model.reactions.get_by_id(rxn.id).bounds = (-10000., 1000000.)
            new_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
        # Set the bounds of the fecal transport (UFEt_met) reactions to be (0., 1000000.)
        elif rxn.id.startswith("UFEt_"):
            saved_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
            model.reactions.get_by_id(rxn.id).bounds = (0., 1000000.)
            new_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
        # Set the bounds of the microbe secretion/uptake (microbe_IEX_met[u]tr) reactions to be (-1000., 1000.)
        elif "IEX" in rxn.id and rxn.id.endswith("[u]tr"):
            saved_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
            model.reactions.get_by_id(rxn.id).bounds = (-1000., 1000.)
            new_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
        # Set the bounds of the diet transport (DUt_met) reactions to be (0., 1000000.)
        elif rxn.id.startswith("DUt_"):
            saved_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
            model.reactions.get_by_id(rxn.id).bounds = (0., 1000000.)
            new_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
        # Set the bounds of the community biomass reaction to be (0.4, 1.)
        elif rxn.id=="communityBiomass":
            saved_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
            model.reactions.get_by_id(rxn.id).bounds = (0.4, 1.)
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


def fetch_metabolomics(sample_id: str, gcms_filepath: str) -> dict:
    """
    Read in a .csv file of VMH name-matched GC-MS metabolomics data and returns
    a dictionary of the VMH metabolite identifiers and their corresponding
    normalized values.

    Parameters
    ----------
    sample_id : str
        The sample identifier of the model (e.g. Case_1, Control_1, etc.).
    gcms_filepath : str
        Path to the GC-MS metabolomics data file that is 

    Returns
    -------
    dict
        Dictionary of the VMH metabolite identifiers and their corresponding
        normalized values for the given model.

    Raises
    ------
    TypeError
        If the filepath is not a .csv file.

    Notes
    -----
    The GC-MS metabolomics data file must be a .csv file with samples as the
    column headers and metabolite names as the row headers. The sample identifiers
    must be in the format of "Case_1", "Control_1", etc. The metabolite names
    must be VMH identifiers. If the metabolite names are not VMH identifiers,
    see the match_names_to_vmh() function.

    TODO: match the metabolite names to VMH identifiers using the match_names_to_vmh() function
    and add them as the row headers. If the metabolite names cannot be fully matched,
    remove them and re-normalize the data.
    """

    print(f"\n[START] Fetching metabolomics data for {sample_id}...")

    # Read in the GC-MS sample file
    if isinstance(gcms_filepath, str) and gcms_filepath.endswith(".csv"):
        print("\nLoading the metabolomics data...")
        metab_data = pd.read_csv(gcms_filepath, sep = ",", index_col = 0)
    else:
        raise TypeError("The filepath must be a path to a .csv file.")

    # Normalizes the metabolomics data
    norm_metab_data = metab_data.div(metab_data.sum(axis=0), axis=1)

    # Check that the norm_metab_data sums to 1
    assert norm_metab_data.sum(axis=0).all() == 1.

    # Gets the list of metabolites that are non-zero for the sample
    non_zero_sample_metabs_list = []
    for i in range(len(norm_metab_data[sample_id])):
        if norm_metab_data[sample_id][i] != 0.:
            non_zero_sample_metabs_list.append(norm_metab_data[sample_id].index[i])

    # Create a dict of unique metabolite names and GC-MS measurements
    sample_metabs_dict = {}
    for met in non_zero_sample_metabs_list:
        sample_metabs_dict[met] = norm_metab_data[sample_id][met]
    for met, value in sample_metabs_dict.items():
        if isinstance(value, pd.Series):
            sample_metabs_dict[met] = value.sum()

    print(f"\n[DONE] Returning normalized sample-specific metabolomics values for {sample_id}.")
    
    return sample_metabs_dict
