import os

import cobra

from utils import (
    fetch_mbx_constr_list,
    fetch_norm_sample_mbx_data,
    load_model,
    match_names_to_vmh,
    set_default_bounds,
    solve_mbx_constraints,
)
from solve_infeasible_model import adapt_diet_and_minimize_infeasibility


def optimize_model(
    model_input: cobra.Model or str,
    diet_path: str,
    output_path: str,
    silent: bool = False,
    return_outputs: bool = False,
    parallel: bool = False,
) -> dict:
    print(f"\n[STARTED] 'optimize_model' workflow for {model_input}")

    # Load the model
    if isinstance(model_input, cobra.Model):
        model = model_input
    elif isinstance(model_input, str):
        model = load_model(model_input, simple_model_name=True)
    else:
        raise ValueError(
            "The model_input must be a path to a model or a COBRApy model object."
        )

    # Set the default bounds
    set_default_bounds(model, source="cobraGEMM", rxn_type="FEX", silent=silent)

    # set the diet exchange reactions to the feasible bounds, if they exist
    model = adapt_diet_and_minimize_infeasibility(
        model=model,
        diet_file_path=diet_path,
        output_file=f"{output_path}/{model.name}_set_diet_flux.txt",
    )

    #########################################################
    # Part 1: maximize the flux through all FEX reactions
    #########################################################
    # Check the output path for existing FEX flux files, and if they exist, skip part 1
    FEX_flux_exists = False
    if os.path.exists(f"{output_path}/{model.name}_opt_FEX_flux.txt"):
        FEX_flux_exists = True
        maximized_FEX_flux_dict = {}

        # Read the existing maximized FEX fluxes
        with open(f"{output_path}/{model.name}_opt_FEX_flux.txt", "r") as f:
            for line in f:
                rxn_id = line.split(":\t")[0]
                flux = float(line.split(":\t")[1].strip("\n"))
                maximized_FEX_flux_dict[rxn_id] = flux

        FEX_rxn_list = list(maximized_FEX_flux_dict.keys())
        maximized_FEX_flux_list = list(maximized_FEX_flux_dict.values())

        print(
            f"\n[SKIPPED] Part 1: found existing maximized FEX fluxes for {model.name}"
        )

    if not FEX_flux_exists:
        if not parallel:
            print(f"\n[STARTED] Part 1: maximizing FEX fluxes for {model.name}")

        # Fetch all FEX reactions and store them in a list
        FEX_rxn_list = [
            rxn.id
            for rxn in model.reactions
            if rxn.id.startswith("EX_") and rxn.id.endswith("[fe]")
        ]

        counter_max = len(FEX_rxn_list)
        maximized_FEX_flux_list = []
        for counter, rxn in enumerate(FEX_rxn_list, start=1):
            if not parallel:
                print(
                    f"\n\tMaximizing FEX reaction {counter} of {counter_max} for {model.name}"
                )
            model.objective = rxn
            solution = model.optimize()
            maximized_FEX_flux_list.append(solution.objective_value)
            if not parallel:
                print(f"\t\t{rxn}:\t{solution.objective_value}")

        # Create a dictionary of the maximized FEX fluxes
        maximized_FEX_flux_dict = dict(zip(FEX_rxn_list, maximized_FEX_flux_list))

        if not parallel:
            print(f"\n[DONE] Part 1: maximization complete for {model.name}")

    #########################################################
    # Part 2: minimize the flux through all IEX reactions
    #########################################################
    if not parallel:
        print(f"\n[STARTED] Part 2: minimizing IEX fluxes for {model.name}")

    counter_max = len(FEX_rxn_list)
    minimized_IEX_flux_dict = {}
    for counter, i in enumerate(range(len(FEX_rxn_list)), start=1):
        if not parallel:
            print(
                f"\n\tMinimizing IEX reaction {counter} of {counter_max} for {model.name}"
            )
        if maximized_FEX_flux_list[i] != 0.0:
            # Store the old bounds for the FEX reaction
            saved_bounds = model.reactions.get_by_id(FEX_rxn_list[i]).bounds

            # Set the bounds for the FEX reaction to the calculated maximum
            model.reactions.get_by_id(FEX_rxn_list[i]).bounds = (
                maximized_FEX_flux_list[i],
                maximized_FEX_flux_list[i],
            )

            # Rename the FEX reaction to match the metabolite name
            metabolite = FEX_rxn_list[i].replace("EX_", "").replace("[fe]", "") + "[u]"

            # Iterate over all IEX reactions for each metabolite
            for rxn in model.metabolites.get_by_id(metabolite).reactions:
                # If it is an IEX reaction, minimize the reaction flux
                if "IEX" in rxn.id:
                    model.objective = model.reactions.get_by_id(rxn.id)
                    solution = model.optimize(objective_sense="minimize")
                    minimized_IEX_flux_dict[rxn.id] = solution.objective_value
                    if not parallel:
                        print(f"\t\t{rxn.id}:\t{solution.objective_value}")

            # Restore the bounds for the minimized IEX reaction
            model.reactions.get_by_id(FEX_rxn_list[i]).bounds = saved_bounds
        elif not parallel:
            print(
                f"\t\tSkipping the reaction '{FEX_rxn_list[i]}' because the maximized flux is 0.0"
            )

    model_rxn_bounds_dict = {rxn.id: rxn.bounds for rxn in model.reactions}
    # Write the maximized FEX fluxes and minimized IEX fluxes to files
    with open(f"{output_path}/{model.name}_opt_flux.txt", "w") as f:
        for rxn_id, value in minimized_IEX_flux_dict.items():
            f.write(f"{rxn_id}:\t{value}\n")

    if not FEX_flux_exists:
        with open(f"{output_path}/{model.name}_opt_FEX_flux.txt", "w") as f:
            for rxn_id in maximized_FEX_flux_dict.keys():
                f.write(f"{rxn_id}:\t{maximized_FEX_flux_dict[rxn_id]}\n")

    if not parallel:
        print(f"\n[DONE] Part 2: minimization complete for {model.name}")
    else:
        print(f"\n[DONE] 'optimize_model' workflow for {model_input}")

    if return_outputs:
        return maximized_FEX_flux_dict, minimized_IEX_flux_dict, model_rxn_bounds_dict


def optimize_model_mbx(
    model_input: str,
    diet_path: str,
    mbx_path: str,
    mbx_matched_keys_input: str or dict,
    output_path: str,
    silent: bool = False,
    verbose: bool = True,
    return_outputs: bool = False,
    parallel: bool = False,
) -> dict:
    print(f"\n[STARTED] 'optimize_model_mbx' workflow for {model_input}")

    # Load the model
    if isinstance(model_input, cobra.Model):
        model = model_input
    elif isinstance(model_input, str):
        model = load_model(model_input)
    else:
        raise ValueError(
            "The model_input must be a path to a model or a COBRApy model object."
        )

    # Set the solver interface
    model.solver = "gurobi"

    if isinstance(mbx_matched_keys_input, str):
        # If the input path does not exist, perform the matching; otherwise, load the dictionary
        if not os.path.exists(mbx_matched_keys_input):
            mbx_output_path = f"{output_path}/vmh_mbx_matched_keys.txt"
            mbx_matched_keys_dict = match_names_to_vmh(
                mbx_filepath=mbx_path,
                output_filepath=mbx_output_path,
                reuturn_matched_keys=True,
                silent=silent,
            )
        elif os.path.exists(mbx_matched_keys_input):
            mbx_matched_keys_dict = {}
            with open(mbx_matched_keys_input, "r") as f:
                for line in f:
                    mbx_matched_keys_dict[line.split(":\t")[0]] = line.split(":\t")[
                        1
                    ].strip("\n")
    elif isinstance(mbx_matched_keys_input, dict):
        if len(mbx_matched_keys_input) == 0:
            raise ValueError(
                "The VMH Identifier name-matching dictionary must not be empty."
            )
        mbx_matched_keys_dict = mbx_matched_keys_input
    else:
        raise ValueError(
            "The mbx_matched_keys_input must be a path to a matched file or a dictionary of the matched metabolomics names to the model metabolite."
        )

    # Get the model-specific metabolomics data
    mbx_metab_norm_dict = fetch_norm_sample_mbx_data(
        model_input=model,
        mbx_filepath=mbx_path,
        matched_mbx_names=mbx_matched_keys_dict,
    )

    # Set all fecal exchange reaction lower bounds to zero
    set_default_bounds(model=model, rxn_type="FEX", silent=silent)

    # set the diet exchange reactions to the feasible bounds, if they exist
    model = adapt_diet_and_minimize_infeasibility(
        model=model,
        diet_file_path=diet_path,
        output_file=f"{output_path}/{model.name}_set_diet_flux_mbx.txt",
    )

    # Fetch and test the constraint list
    mbx_constraints = fetch_mbx_constr_list(
        model=model, mbx_metab_norm_dict=mbx_metab_norm_dict
    )

    # Fetch the slack constraints if needed
    mbx_constr, constr_details = solve_mbx_constraints(
        model=model, constraints=mbx_constraints, parallel=parallel
    )

    # Add the constraints to the model
    model.add_cons_vars(mbx_constr)
    model.solver.update()

    # Optimize the model by minimizing the fluxes of reactions for each metabolite
    min_opt_solutions = {}
    if not parallel:
        print(f"\n[Minimizing the model {model.name} for each metabolite]")
    for rxn_id in [
        rxn.id
        for rxn in model.reactions
        if "_IEX_" in rxn.id and rxn.id.endswith("[u]tr")
    ]:
        model.objective = model.reactions.get_by_id(rxn_id)
        solution = model.optimize(objective_sense="minimize")
        if (solution.objective_value != 0.0 or verbose) and not parallel:
            print(
                f"\tMinimized: {model.reactions.get_by_id(rxn_id).id}:\t{solution.objective_value}"
            )
        min_opt_solutions[rxn_id] = solution.objective_value

    # # Optimize the model by maximizing the fluxes of reactions for each metabolite
    # max_opt_solutions = dict()
    # if not parallel:
    #     print(f"\n[Maximizing the model {model.name} for each metabolite]")
    # for rxn_id in [
    #     rxn.id
    #     for rxn in model.reactions
    #     if "_IEX_" in rxn.id and rxn.id.endswith("[u]tr")
    # ]:
    #     model.objective = model.reactions.get_by_id(rxn_id)
    #     solution = model.optimize(objective_sense="maximize")
    #     if solution.objective_value != 0.0 or verbose:
    #         if not parallel:
    #             print(
    #                 f"\tMaximized: {model.reactions.get_by_id(rxn_id).id}:\t{solution.objective_value}"
    #             )
    #     max_opt_solutions[rxn_id] = solution.objective_value

    # Save the results
    with open(f"{output_path}/{model.name}_slack_var_log.txt", "w") as f:
        for line in constr_details:
            f.write(f"{line}\n")
    with open(f"{output_path}/{model.name}_min_mbx_opt_flux.txt", "w") as f:
        for rxn_id in min_opt_solutions:
            f.write(f"{rxn_id}:\t{min_opt_solutions[rxn_id]}\n")
    # with open(f"{output_path}/{model.name}_max_mbx_opt_flux.txt", "w") as f:
    #     for rxn_id in max_opt_solutions:
    #         f.write(f"{rxn_id}:\t{max_opt_solutions[rxn_id]}\n")

    if parallel:
        print(f"\n[DONE] 'optimize_model_mbx' workflow for {model_input}")

    if return_outputs:
        return min_opt_solutions

