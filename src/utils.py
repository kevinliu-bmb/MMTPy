import os
import re
from math import isnan

import cobra
import numpy as np
import pandas as pd


def print_logo(tool: str, tool_description: str, version: str) -> None:
    """
    The `print_logo` function prints a logo, tool name, and version for a given
    tool.

    :param tool: The `tool` parameter is a string that represents the name of the
    tool. It is used to display the name of the tool in the logo and tool
    information
    :type tool: str
    :param tool_description: The `tool_description` parameter is a string that
    describes the tool. It provides information about what the tool does, its
    purpose, and any other relevant details
    :type tool_description: str
    :param version: The version parameter is a string that represents the version
    of the tool. It could be something like "1.0.0" or "v2.3.5"
    :type version: str
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


def reset_bounds(
    model: cobra.Model,
    source: str = "cg",
    rxn_type: str = "all",
) -> None:
    """
    The function `reset_bounds` sets bounds for specific reactions in a given model
    based on specified criteria.

    :param model: The `model` parameter is a `cobra.Model` object, which represents
    a metabolic model. It contains information about the reactions, metabolites,
    and constraints of the model
    :type model: cobra.Model
    :param source: The `source` parameter in the `reset_bounds` function is used to
    specify the source of default bounds to be used for setting the reaction
    bounds. It can take the following values:, defaults to cg
    :type source: str (optional)
    :param rxn_type: The `rxn_type` parameter in the `reset_bounds` function is
    used to specify the type of reactions for which the bounds should be reset. It
    can take the following values:, defaults to all
    :type rxn_type: str (optional)
    """
    print(f"\n[Setting bounds for {rxn_type} reactions using {source} defaults]")

    saved_bounds = {}
    new_bounds = {}
    for rxn in model.reactions:
        if (
            rxn.id.startswith("EX_")
            and rxn.id.endswith("[fe]")
            and "microbeBiomass" not in rxn.id
            and rxn_type in {"all", "fex"}
        ):
            saved_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
            if source == "cg-secrete":
                model.reactions.get_by_id(rxn.id).bounds = (0.0, 1000000.0)
            elif source in {"mmt", "cg"}:
                model.reactions.get_by_id(rxn.id).bounds = (-1000.0, 1000000.0)
            else:
                raise ValueError(f"Source {source} not supported")
            new_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
        elif (
            rxn.id.startswith("EX_")
            and rxn.id.endswith("[fe]")
            and "microbeBiomass" in rxn.id
            and rxn_type in {"all", "fex"}
        ):
            saved_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
            if source == "cg-secrete":
                model.reactions.get_by_id(rxn.id).bounds = (0.0, 1000000.0)
            elif source in {"mmt", "cg"}:
                model.reactions.get_by_id(rxn.id).bounds = (-10000.0, 1000000.0)
            else:
                raise ValueError(f"Source {source} not supported")
            new_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
        elif rxn.id.startswith("UFEt_") and rxn_type in {"all", "ufet"}:
            saved_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
            model.reactions.get_by_id(rxn.id).bounds = (0.0, 1000000.0)
            new_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
        elif (
            "iex" in rxn.id and rxn.id.endswith("[u]tr") and rxn_type in {"all", "iex"}
        ):
            saved_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
            model.reactions.get_by_id(rxn.id).bounds = (-1000.0, 1000.0)
            new_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
        elif rxn.id.startswith("DUt_") and rxn_type in {"all", "dut"}:
            saved_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
            model.reactions.get_by_id(rxn.id).bounds = (0.0, 1000000.0)
            new_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
        elif rxn.id == "communityBiomass" and rxn_type in {"all", "commBiomass"}:
            saved_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
            model.reactions.get_by_id(rxn.id).bounds = (0.4, 1.0)
            new_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds

    n_changed_bounds = 0
    # Print out the changes
    changed_bounds = [
        (rxn, bounds, new_bounds[rxn])
        for rxn, bounds in saved_bounds.items()
        if bounds != new_bounds[rxn]
    ]
    for rxn, bounds, new_bound in changed_bounds:
        print(f"\tChanged bounds for {rxn} from {bounds} to {new_bound}")
    n_changed_bounds = len(changed_bounds)

    if n_changed_bounds == 0:
        print("\tNo bounds were changed")
    else:
        print(
            f"\tChanged {n_changed_bounds}/{len(saved_bounds)} {rxn_type} reaction bounds for {model.name}"
        )


def convert_model_format(model_path: str or cobra.Model, output_path: str = None):
    """Convert a mgPipe.m (Heinken et al., 2022) MATLAB model to a json model.

    Parameters:
    model_path : str or cobra.Model
        Path to the model file or a COBRApy model loaded into memory.
    output_path : str
        Path to the output file.

    Raises:
    ValueError
        If model_path is not a string or a COBRApy model.

    Notes:
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
    """Convert a string to a standard format for matching.

    Parameters:
    s : str
        String to be converted.

    Returns:
    str
        Converted string.

    Notes:
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
    """Get the index of the first numerical column in a dataframe.

    Parameters:
    df : pd.DataFrame
        Dataframe to be searched.

    Raises:
    ValueError
        If no numerical columns are found.

    Returns:
    int
        Index of the first numerical column.
    """
    for index, (col, dtype) in enumerate(df.dtypes.items()):
        if np.issubdtype(dtype, np.number):
            return index
    raise ValueError("No numerical columns found.")


def fetch_norm_sample_mbx_data(
    model_input: cobra.Model or str,
    mbx_filepath: str,
    matched_mbx_names: dict,
) -> dict:
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

    idx_list = [
        mbx_data.columns.get_loc(col_name)
        for col_name in mbx_data.columns
        if col_name in matched_mbx_names
    ]
    sample_id = [index for index in mbx_data.index if index == model.name]
    if len(sample_id) > 1:
        ValueError("Multiple sample IDs found in MBX data. Please check the data.")
    elif not sample_id:
        ValueError("No sample ID found in MBX data. Please check the data.")
    else:
        sample_id = sample_id[0]

    # Given a sample ID, get the MBX data for that sample
    sample_mbx_data = mbx_data.loc[sample_id][min(idx_list) :]

    # Create a dictionary of metabolite names and their concentrations
    metab_raw_vals_dict = dict(sample_mbx_data.items())

    # If the values of metab_raw_vals_dict are strings, remove any commas and convert them to floats
    for k, v in metab_raw_vals_dict.items():
        if isinstance(v, str):
            metab_raw_vals_dict[k] = float(v.replace(",", ""))

    # Check if all values in metab_raw_vals_dict are floats
    if not all(isinstance(v, float) for v in metab_raw_vals_dict.values()):
        raise ValueError("Not all values in metab_raw_vals_dict are floats.")

    # Raise ValueError if matched_mbx_names is empty
    if not matched_mbx_names:
        raise ValueError("The VMH Identifier name-matching dictionary is empty.")

    # Create a dictionary of VMH IDs and their corresponding metabolite values if the metabolite is in the model
    vmh_id_values = {}
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
    """Compute the MBX associated FEX reaction constraints for a sample and tests
    if the addition of the constraints gives a feasible solution.

    Parameters:
    model : cobra.Model
        The model to be constrained.
    mbx_metab_norm_dict : dict
        A dictionary of metabolite names and their normalized values.

    Returns:
    list
        A list of the constraints added to the model.

    Raises:
    ValueError
        If no constraints were added to the model.
    """
    print(f"\n[Fetching MBX constraints for {model.name}]")
    # Calculate the constraint for each metabolite's fecal exchange reactions
    metab_constrained_flux_expr = {}
    for vmh_id, mbx_value in mbx_metab_norm_dict.items():
        if mbx_value != 0.0 and f"EX_{vmh_id}[fe]" in model.reactions:
            flux_expr = model.reactions.get_by_id(f"EX_{vmh_id}[fe]").flux_expression
            metab_constrained_flux_expr[f"EX_{vmh_id}[fe]"] = flux_expr - (
                mbx_value
                * sum(
                    model.reactions.get_by_id(f"EX_{metab}[fe]").flux_expression
                    for metab in mbx_metab_norm_dict
                    if f"EX_{metab}[fe]" in model.reactions
                )
            )

    # Add the constraints to a list for each metabolite with MBX data
    constraint_list = []
    for ex_rxn_id, constr_expr in metab_constrained_flux_expr.items():
        constraint = model.problem.Constraint(
            constr_expr, lb=0, ub=0, name=ex_rxn_id + "_constraint"
        )
        constraint_list.append(constraint)

    if not constraint_list:
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
    """Add slack variables to infeasible constraints and test if the solution is feasible.

    Parameters:
    model : cobra.Model
        The model to be constrained.
    constraints : list
        A list of constraints to be added to the model.
    parallel : bool, optional
        If True, no detailed outputs will be printed, by default False

    Returns:
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
