import cobra

from cobra_utils import (
    fetch_mbx_constr_list,
    fetch_norm_sample_mbx_data,
    load_model,
    print_logo,
    set_default_bounds,
    slack_constraints,
)


def optimize_model_mbx(
    model_input: str = "example_data/models/microbiota_model_diet_Case_1_18_month.json",
    mbx_path: str = "example_data/metabolomics_data.csv",
    output_path: str = "example_outputs",
):
    """
    Optimize the model for each metabolite in the metabolomics data.

    Parameters
    ----------
    model_input : str
        Path to the model file.
    mbx_path : str
        Path to the metabolomics data file.
    output_path : str
        Path to the output directory.

    Returns
    -------
    dict
        Dictionary of the optimized solutions for each metabolite.
    """
    # Print the logo
    print_logo(
        tool="optimize_model_mbx",
        tool_description="Optimize a model for each metabolite given the MBX data.",
        version="0.1.0-beta",
    )

    # Load the model
    if isinstance(model_input, cobra.Model):
        model = model_input
    elif isinstance(model_input, str):
        print("\nLoading the model...")
        model = load_model(model_input)
    else:
        raise ValueError(
            "The model_input must be a path to a model or a COBRApy model object."
        )

    # Set the solver interface
    model.solver = "gurobi"

    # Get the model-specific metabolomics data
    mbx_metab_norm_dict = fetch_norm_sample_mbx_data(
        model_input=model,
        mbx_filepath=mbx_path,
        match_key_output_filepath=output_path,
    )

    # Set all fecal exchange reaction lower bounds to zero
    set_default_bounds(model, rxn_type="FEX")

    # Fetch and test the constraint list
    mbx_constraints = fetch_mbx_constr_list(model, mbx_metab_norm_dict)

    # Fetch the slack constraints if needed
    feasible_constr, slacked_constr = slack_constraints(model, mbx_constraints)

    # Add the constraints to the model
    model.add_cons_vars(feasible_constr + slacked_constr)
    model.solver.update()

    opt_solutions = dict()
    # Optimize the model for each metabolite
    print(f"\n[Optimizing the model {model.name} for each metabolite]")
    for rxn_id in [
        rxn.id
        for rxn in model.reactions
        if rxn.id.startswith("EX_") and rxn.id.endswith("[fe]")
    ]:
        model.objective = model.reactions.get_by_id(rxn_id)
        solution = model.optimize()
        print(f"\t{model.reactions.get_by_id(rxn_id).id}:\t{solution.objective_value}")
        opt_solutions[rxn_id] = solution.objective_value

    # Save the results
    with open(f"{output_path}/optimized_mbx_model_results.txt", "w") as f:
        for rxn_id in opt_solutions:
            f.write(f"{rxn_id}\t{opt_solutions[rxn_id]}\n")

    return opt_solutions
