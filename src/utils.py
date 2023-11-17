import os
import re
import time
from math import isnan

import cobra

import numpy as np
import pandas as pd
import pubchempy as pcp


def load_model(model_path: str, simple_model_name: bool = True) -> cobra.Model:
    if not os.path.exists(model_path):
        raise ValueError(f"Model path does not exist: {model_path}")

    print(f"[Loading model from {model_path}]")

    extension_to_function = {
        ".xml": cobra.io.read_sbml_model,
        ".json": cobra.io.load_json_model,
        ".yml": cobra.io.load_yaml_model,
        ".mat": cobra.io.load_matlab_model,
        ".sbml": cobra.io.read_sbml_model,
    }
    file_extension = os.path.splitext(model_path)[1]
    if file_extension in extension_to_function:
        model = extension_to_function[file_extension](model_path)
    else:
        raise ValueError(f"Model format not supported for {model_path}")

    model_file_name = os.path.basename(model_path).split(".")[0]
    model.name = (
        re.compile(r"(Case|Control)_\d+").search(model_file_name)[0]
        if simple_model_name
        else model_file_name
    )

    print(f"[{model.name} loaded]")
    return model


def set_default_bounds(
    model: cobra.Model,
    source: str = "cobraGEMM",
    rxn_type: str = "all",
    silent: bool = False,
) -> bool:
    print(f"[Setting bounds for {rxn_type} reactions using {source} defaults]")

    saved_bounds = {}
    new_bounds = {}
    for rxn in model.reactions:
        # Set the bounds of the fecal exchange (EX_met[fe]) reactions to be (-1000., 1000000.)
        if (
            rxn.id.startswith("EX_")
            and rxn.id.endswith("[fe]")
            and "microbeBiomass" not in rxn.id
            and rxn_type in {"all", "FEX"}
        ):
            saved_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
            if source == "MATLAB":
                model.reactions.get_by_id(rxn.id).bounds = (-1000.0, 1000000.0)
            elif source == "cobraGEMM":
                model.reactions.get_by_id(rxn.id).bounds = (0.0, 1000000.0)
            else:
                raise ValueError(f"Source {source} not supported")
            new_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
        elif (
            rxn.id.startswith("EX_")
            and rxn.id.endswith("[fe]")
            and "microbeBiomass" in rxn.id
            and rxn_type in {"all", "FEX"}
        ):
            saved_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
            if source == "MATLAB":
                model.reactions.get_by_id(rxn.id).bounds = (-10000.0, 1000000.0)
            elif source == "cobraGEMM":
                model.reactions.get_by_id(rxn.id).bounds = (0.0, 1000000.0)
            else:
                raise ValueError(f"Source {source} not supported")
            new_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
        elif rxn.id.startswith("UFEt_") and rxn_type in {"all", "UFEt"}:
            saved_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
            model.reactions.get_by_id(rxn.id).bounds = (0.0, 1000000.0)
            new_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
        elif (
            "IEX" in rxn.id and rxn.id.endswith("[u]tr") and rxn_type in {"all", "IEX"}
        ):
            saved_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
            model.reactions.get_by_id(rxn.id).bounds = (-1000.0, 1000.0)
            new_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
        elif rxn.id.startswith("DUt_") and rxn_type in {"all", "DUt"}:
            saved_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
            model.reactions.get_by_id(rxn.id).bounds = (0.0, 1000000.0)
            new_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
        elif rxn.id == "communityBiomass" and rxn_type in {"all", "commBiomass"}:
            saved_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
            model.reactions.get_by_id(rxn.id).bounds = (0.4, 1.0)
            new_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds

    n_changed_bounds = 0
    # Print out the changes
    for rxn, bounds in saved_bounds.items():
        if bounds != new_bounds[rxn]:
            if not silent:
                print(f"Changed bounds for {rxn} from {bounds} to {new_bounds[rxn]}")
            n_changed_bounds += 1

    if n_changed_bounds == 0:
        bounds_changed = False
        print("No bounds were changed")
    else:
        bounds_changed = True
        print(
            f"Changed {n_changed_bounds}/{len(saved_bounds)} {rxn_type} reaction bounds for {model.name}"
        )

    return bounds_changed


def convert_model_format(model_path: str or cobra.Model, output_path: str = None):
    model = load_model(model_path) if isinstance(model_path, str) else model_path
    print(f"[Converting model {model.name} to json format]")

    for metab in model.metabolites:
        metab.charge = "nan" if isnan(metab.charge) else metab.charge

    os.makedirs(output_path, exist_ok=True)
    converted_output_filepath = os.path.join(output_path, f"{model.name}.json")
    cobra.io.save_json_model(model, converted_output_filepath)
    print(f"[{model.name} converted to json format]")


def convert_string(s: str) -> str:
    def replacement(match):
        if match.lastgroup == "dot":
            return match.group(1) + "," + match.group(2)
        elif match.lastgroup == "space":
            return match.group(1) + "-" + match.group(2)
        return ""

    pattern = re.compile(
        r"""
        \s*\(.*?\)|
        (?P<dot>(\d)\.(\d))|
        (?P<space>(\d) (\w))
    """,
        re.VERBOSE,
    )

    return pattern.sub(replacement, s)


def get_init_mbx_idx(df: pd.DataFrame) -> int:
    numerical_column = df.dtypes.apply(lambda x: np.issubdtype(x, np.number))
    if numerical_column.any():
        return numerical_column.idxmax()
    raise ValueError("No numerical columns found.")


def match_names_to_vmh(
    mbx_filepath: str,
    output_filepath: str,
    reuturn_matched_keys: bool,
    vmh_db_filepath: str,
    manual_matching_filepath: str,
    silent: bool = False,
) -> dict:
    # Load the data
    vmh_db_df = pd.read_csv(vmh_db_filepath, index_col=0, header=0, delimiter="")
    mbx_data_df = pd.read_csv(mbx_filepath, index_col=0, header=0)

    print("[Matching MBX names to VMH identifiers]")

    print(
        "[1/3] Direct matching of MBX names to VMH identifiers using the VMH database"
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

    print("[2/3] Matching of MBX names to VMH identifiers via PubChemPy")
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
            print(f"Error getting info for '{compound}': {e}")
        time.sleep(0.5)

    pubchempy_matched_dict = {}

    for vmh_id, vmh_inchikey in vmh_inchikey_dict.items():
        for mbx_name, pubchempy_inchikey in inchikey_names_dict.items():
            if vmh_inchikey != "nan" and vmh_inchikey == pubchempy_inchikey:
                if not silent:
                    print(f"Matched '{mbx_name}' to '{vmh_id}' using InChIKey")
                pubchempy_matched_dict[mbx_name] = vmh_id
    for vmh_id, vmh_cid in vmh_cid_dict.items():
        for mbx_name, pubchempy_cid in cid_names_dict.items():
            if vmh_cid != "nan" and vmh_cid == pubchempy_cid:
                if not silent:
                    print(f"Matched '{mbx_name}' to '{vmh_id}' using CID")
                pubchempy_matched_dict[mbx_name] = vmh_id
    for vmh_id, vmh_inchi in vmh_inchistring_dict.items():
        for mbx_name, pubchempy_inchi in inchi_names_dict.items():
            if vmh_inchi != "nan" and vmh_inchi == pubchempy_inchi:
                if not silent:
                    print(f"Matched '{mbx_name}' to '{vmh_id}' using InChI")
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
                        print(f"Matched '{mbx_name}' to '{vmh_id}' using SMILES")
                    pubchempy_matched_dict[mbx_name] = vmh_id

    print("[3/3] Matching of MBX names to VMH identifiers via manual matching database")
    manual_matching = {}
    with open(manual_matching_filepath, "r") as f:
        for line in f:
            name, vmh = line.split("")
            if vmh.endswith(""):
                vmh = vmh[:-1]
            manual_matching[name] = vmh

    max_matched_dict = {
        name: id
        for name, id in manual_matching.items()
        if id not in pubchempy_matched_dict.values()
    } | pubchempy_matched_dict
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
            f.write(f"{key}:{value}")
    print(
        f"{len(max_matched_dict)} of {len(mbx_data_df.columns)-2} VMH identifiers matched to the MBX metabolite names"
    )

    print(
        f"[Matched MBX names to VMH identifiers and written to '{key_output_filepath}']"
    )

    if reuturn_matched_keys:
        return max_matched_dict


def fetch_norm_sample_mbx_data(
    model_input: cobra.Model or str,
    mbx_filepath: str,
    matched_mbx_names: dict,
) -> dict:
    model = load_model(model_input) if isinstance(model_input, str) else model_input
    print(f"[Fetching MBX data for {model.name}]")

    mbx_data = pd.read_csv(mbx_filepath, sep=",", index_col=0)
    sample_id = next(
        iter([index for index in mbx_data.index if index == model.name]), None
    )
    if sample_id is None:
        raise ValueError("No sample ID found in MBX data. Please check the data.")

    sample_mbx_data = mbx_data.loc[sample_id]
    metab_raw_vals_dict = {
        k: float(v.replace(",", "")) if isinstance(v, str) else v
        for k, v in sample_mbx_data.items()
    }

    if not matched_mbx_names:
        raise ValueError("The VMH Identifier name-matching dictionary is empty.")

    vmh_id_values = {
        vmh_id: metab_raw_vals_dict[vmh_name]
        for vmh_name, vmh_id in matched_mbx_names.items()
        if f"EX_{vmh_id}[fe]" in model.reactions
    }

    total_val = sum(vmh_id_values.values())
    norm_vmh_id_vals = {k: v / total_val for k, v in vmh_id_values.items()}

    print(f"[Returning normalized sample-specific MBX values for {sample_id}]")
    return norm_vmh_id_vals


def fetch_mbx_constr_list(model: cobra.Model, mbx_metab_norm_dict: dict) -> list:
    print(f"[Fetching MBX constraints for {model.name}]")
    constraint_list = [
        model.problem.Constraint(
            model.reactions.get_by_id(f"EX_{vmh_id}[fe]").flux_expression
            - (
                mbx_value
                * sum(
                    model.reactions.get_by_id(f"EX_{metab}[fe]").flux_expression
                    for metab in mbx_metab_norm_dict
                    if f"EX_{metab}[fe]" in model.reactions
                )
            ),
            lb=0,
            ub=0,
            name=f"EX_{vmh_id}[fe]_constraint",
        )
        for vmh_id, mbx_value in mbx_metab_norm_dict.items()
        if mbx_value != 0.0 and f"EX_{vmh_id}[fe]" in model.reactions
    ]

    if not constraint_list:
        raise ValueError("No constraints were added to the model.")

    model.add_cons_vars(constraint_list)
    model.solver.update()
    model.objective = 0
    solution = model.optimize()
    if solution.status == "infeasible":
        model.remove_cons_vars(constraint_list)
        model.solver.update()
        print(
            "Warning: the solution is infeasible by introducing the constraints without slack variables"
        )
    else:
        print("The solution is feasible without adding slack variables")
        model.remove_cons_vars(constraint_list)
        model.solver.update()

    return constraint_list


def solve_mbx_constraints(
    model: cobra.Model, constraints: list, parallel: bool = False
) -> list:
    slack_variables = [
        model.problem.Variable(f"{constraint.name}_slack", lb=0, ub=None)
        for constraint in constraints
    ]
    slack_constraints = [
        model.problem.Constraint(
            constraint.expression + slack_variable,
            lb=constraint.lb,
            ub=constraint.ub,
            name=f"{constraint.name}_slack",
        )
        for constraint, slack_variable in zip(constraints, slack_variables)
    ]

    model.add_cons_vars(slack_constraints)
    model.solver.update()
    model.objective = model.problem.Objective(sum(slack_variables))
    solution = model.optimize(objective_sense="minimize")

    if solution.status == "optimal":
        print("The solution is feasible with slack variables added to all constraints")
    elif solution.status == "infeasible":
        print(
            "Warning: The solution is infeasible with slack variables added to all constraints"
        )

    clear_cons_vars(model, slack_constraints, slack_variables)
    return slack_constraints


def clear_cons_vars(model, slack_constraints, slack_variables):
    model.remove_cons_vars(slack_variables)
    model.remove_cons_vars(slack_constraints)
    model.solver.update()
