import os
import re

import pandas as pd
import pubchempy as pcp

from cobra_utils import print_logo


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
) -> None:
    """
    Map the metabolite names detected by GC-MS to VMH identifiers for a given
    GC-MS dataset.

    Parameters
    ----------
    gcms_filepath : str
        Filepath to the GC-MS data.
    output_filepath : str
        Filepath for saving the matching keys.
    vmh_db_filepath : str
        Filepath to the VMH database of metabolites and their identifiers.

    Returns
    -------
    None

    Notes
    -----
    The metabolomics data must be a .csv file with samples as the column headers
    and metabolite names as the row indicies. The metabolite names must be canonical
    metabolite names (e.g. "L-Aspartate", "L-Asparagine", etc.). The matched keys
    are saved as a tab-delimited .txt file with the following format:
        {GC-MS name} {VMH identifier}
    """

    print_logo(
        tool="match-names-to-vmh",
        tool_description="Matches metabolite names to VMH identifiers via pubchempy.",
        version="0.1.0",
    )

    # Load the data
    vmh_db = pd.read_csv(vmh_db_filepath, index_col=0, header=0, delimiter="\t")
    gcms_data = pd.read_csv(gcms_filepath, index_col=0, header=0)

    print("\n[START] Matching GC-MS names to VMH identifiers...")

    # Create dictionaries for direct matching
    gcms_names_dict = {
        name: convert_string(name).lower() for name in gcms_data.columns[2:].to_list()
    }
    vmh_names_dict = dict(
        zip(vmh_db["fullName"].index, vmh_db["fullName"], strict=False)
    )

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
    vmh_cid_dict = dict(
        zip(vmh_db["pubChemId"].index, vmh_db["pubChemId"], strict=False)
    )
    vmh_inchikey_dict = dict(
        zip(vmh_db["inchiKey"].index, vmh_db["inchiKey"], strict=False)
    )
    vmh_inchistring_dict = dict(
        zip(vmh_db["inchiString"].index, vmh_db["inchiString"], strict=False)
    )
    vmh_smiles_dict = dict(zip(vmh_db["smile"].index, vmh_db["smile"], strict=False))

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
    pubchempy_matched_dict.update(direct_matching_dict)

    if len(pubchempy_matched_dict) != len(gcms_data.index):
        for vmh_id, vmh_smiles in vmh_smiles_dict.items():
            for gcms_name, pubchempy_smiles in smiles_names.items():
                if vmh_smiles != "nan" and vmh_smiles == pubchempy_smiles:
                    print(f"\nMatched {gcms_name} to {vmh_id} using SMILES")
                    pubchempy_matched_dict[gcms_name] = vmh_id

    print(
        f"\n{len(pubchempy_matched_dict)} of {len(gcms_data.columns)-2} VMH identifiers matched."
    )

    # If the output filepath does not exist, create it
    if not os.path.exists(output_filepath):
        os.makedirs(output_filepath)

    if output_filepath[-1] != "/":
        output_filepath += "/"

    key_output_filepath = f"{gcms_filepath.split('.')[-2]}_matched_key.txt"

    # Write out the matched identifiers to a .txt file
    with open(key_output_filepath, "w") as f:
        for key, value in pubchempy_matched_dict.items():
            f.write(f"{key}\t{value}\n")

    print(
        f"\n[DONE] Matched GC-MS names to VMH identifiers and written to {key_output_filepath}"
    )
