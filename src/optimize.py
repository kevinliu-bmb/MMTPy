import os

import cobra

def optimize_model_mbx(
    model_input: str,
    diet_path: str,
    met_path: str,
    met_matched_keys_input: str or dict,
    output_path: str,
    verbose: bool = True
) -> dict:
    """Optimize the model for each metabolite in the metabolomics data.

    Parameters:
    model_input : str or cobra.Model
        Path to the model file or a COBRApy model object.
    diet_path : str
        Path to the diet file used to generate the model.
    met_path : str
        Path to the metabolomics data file.
    met_matched_keys_input : str or dict
        Dictionary of the matched metabolomics names to the model metabolite or path to the matched file.
    output_path : str
        Path to the output directory.
    silent : bool
        If True, the function will not print the boundary changes and VMH name matching outputs.
    verbose : bool
        If True, the function will print the zero-valued optimization solution outputs.
    return_outputs : bool
        If True, the function will return the maximized and minimized IEX fluxes.
    parallel : bool
        If True, the function will not print the progress outputs.

    Returns:
    dict
        Dictionary of the optimized solutions for each metabolite.

    Raises:
    ValueError
        If the model_input is not a path to a model or a COBRApy model object.
    ValueError
        If the mbx_matched_keys_input is not a path to a matched file or a dictionary of the matched metabolomics names to the model metabolite.
    ValueError
        If the VMH Identifier name-matching dictionary is empty.
    """
    # Print the logo
    if not parallel:
        print_logo(
            tool="optimize_model_mbx",
            tool_description="Optimize a model for each metabolite given the MBX data.",
            version="0.1.0",
        )
    else:
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

    if isinstance(met_matched_keys_input, str):
        # If the input path does not exist, perform the matching; otherwise, load the dictionary
        if not os.path.exists(met_matched_keys_input):
            mbx_output_path = f"{output_path}/vmh_mbx_matched_keys.txt"
            mbx_matched_keys_dict = match_names_to_vmh(
                mbx_filepath=met_path,
                output_filepath=mbx_output_path,
                reuturn_matched_keys=True,
            )
        elif os.path.exists(met_matched_keys_input):
            mbx_matched_keys_dict = dict()
            with open(mbx_matched_keys_input, "r") as f:
                for line in f:
                    mbx_matched_keys_dict[line.split(":\t")[0]] = line.split(":\t")[
                        1
                    ].strip("\n")
    elif isinstance(met_matched_keys_input, dict):
        if len(met_matched_keys_input) == 0:
            raise ValueError(
                "The VMH Identifier name-matching dictionary must not be empty."
            )
        mbx_matched_keys_dict = met_matched_keys_input
    else:
        raise ValueError(
            "The mbx_matched_keys_input must be a path to a matched file or a dictionary of the matched metabolomics names to the model metabolite."
        )

    # Get the model-specific metabolomics data
    mbx_metab_norm_dict = fetch_norm_sample_mbx_data(
        model_input=model,
        mbx_filepath=met_path,
        matched_mbx_names=mbx_matched_keys_dict,
    )

    # Set all fecal exchange reaction lower bounds to zero
    set_default_bounds(model=model, rxn_type="FEX")

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
        model=model, constraints=mbx_constraints
    )

    # Add the constraints to the model
    model.add_cons_vars(mbx_constr)
    model.solver.update()

    # Optimize the model by minimizing the fluxes of reactions for each metabolite
    min_opt_solutions = {}
    print(f"\n[Minimizing the model {model.name} for each metabolite]")
    for rxn_id in [
        rxn.id
        for rxn in model.reactions
        if "_IEX_" in rxn.id and rxn.id.endswith("[u]tr")
    ]:
        model.objective = model.reactions.get_by_id(rxn_id)
        solution = model.optimize(objective_sense="minimize")
        if solution.objective_value != 0.0 or verbose:
            print(f"\tMinimized: {model.reactions.get_by_id(rxn_id).id}:\t{solution.objective_value}"
                )
        min_opt_solutions[rxn_id] = solution.objective_value

    max_opt_solutions = {}
    print(f"\n[Maximizing the model {model.name} for each metabolite]")
    for rxn_id in [
        rxn.id
        for rxn in model.reactions
        if "_IEX_" in rxn.id and rxn.id.endswith("[u]tr")
    ]:
        model.objective = model.reactions.get_by_id(rxn_id)
        solution = model.optimize(objective_sense="maximize")
        if solution.objective_value != 0.0 or verbose:
            print(f"\tMaximized: {model.reactions.get_by_id(rxn_id).id}:\t{solution.objective_value}"
                )
        max_opt_solutions[rxn_id] = solution.objective_value

    # Save the results
    with open(f"{output_path}/{model.name}_slack_var_log.txt", "w") as f:
        for line in constr_details:
            f.write(f"{line}\n")
    with open(f"{output_path}/{model.name}_min_mbx_opt_flux.txt", "w") as f:
        for rxn_id in min_opt_solutions:
            f.write(f"{rxn_id}:\t{min_opt_solutions[rxn_id]}\n")
    with open(f"{output_path}/{model.name}_max_mbx_opt_flux.txt", "w") as f:
        for rxn_id in max_opt_solutions:
            f.write(f"{rxn_id}:\t{max_opt_solutions[rxn_id]}\n")

    print(f"\n[DONE] 'optimize_model_mbx' workflow for {model_input}")

    return min_opt_solutions
