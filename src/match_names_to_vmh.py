from utils import convert_string, get_init_mbx_idx


import pandas as pd
import pubchempy as pcp


import os
import time


def match_names_to_vmh(
    mbx_filepath: str,
    output_filepath: str,
    reuturn_matched_keys: bool,
    vmh_db_filepath: str = "workflows/optimization/data_dependencies/all_vmh_metabolites.tsv",
    manual_matching_filepath: str = "workflows/optimization/data_dependencies/manually_matched_keys.txt",
    silent: bool = False,
) -> dict:
    """Map the metabolite names detected by MBX to VMH identifiers for a given MBX dataset. The matching is performed in the following order:
        1. Direct matching of MBX names to VMH identifiers using the VMH database.
        2. Matching of MBX names to VMH identifiers via pubchempy.
        3. Manual matching of MBX names to VMH identifiers.

    Parameters:
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

    Raises:
    ValueError
        If model_path is not a string or a COBRApy model.

    Returns:
    dict
        Dictionary of matched keys.

    Notes:
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
    direct_matching_dict = {}
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
    iupac_names_dict = {}
    cid_names_dict = {}
    inchi_names_dict = {}
    inchikey_names_dict = {}
    smiles_names_dict = {}

    # Iterate over each item in the compounds dictionary
    for mbx_name, compound in unmatched_dict.items():
        try:
            if c := pcp.get_compounds(compound, "name"):
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

    pubchempy_matched_dict = {}

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
    pubchempy_matched_dict |= {
        value: key for key, value in direct_matching_dict.items()
    }

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

    max_matched_dict |= pubchempy_matched_dict

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
