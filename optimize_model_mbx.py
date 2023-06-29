import cobra

from cobra_utils import load_model, fetch_norm_sample_mbx_data


def optimize_model_metabolomics(
    model_input: cobra.Model or str,
    gcms_filepath: str,
    reaction_type_prefix: str = "UFEt",
    reaction_type_suffix: str = "",
    opt_direction: str = "maximize",
):
    """
    Optimize the model reaction fluxes using metabolomics measurement data (e.g.,
    from GC-MS) to constrain the reaction types of interest.

    Parameters
    ----------
    model_input : cobra.Model or str
        COBRApy-supported model to be optimized.
    gcms_filepath : str
        File path to metabolomics data.
    reaction_type_prefix : str (defaults to "UFEt")
        Prefix of the desired optimization reaction type.
    reaction_type_suffix : str (defaults to "")
        Suffix of the desired optimization reaction type.
    opt_direction : str (defaults to "maximize")
        Direction for optimization; either "maximize" (default) or "minimize".

    Returns
    -------
    dict
        A dictionary containing the optimization objectives and solution values
        for each iteration.

    Raises
    ------
    None

    Notes
    -----
    TODO simulation testing using real dataset.
    """
    # Load the model
    if isinstance(model_input, str):
        model = load_model(model_input)
    elif isinstance(model_input, cobra.Model):
        model = model_input

    # Get the model-specific metabolomics data
    gcms_metab_norm = fetch_norm_sample_mbx_data(
        model_input=model, mbx_filepath=gcms_filepath
    )

    rxn_id_gcms_val = {
        f"{reaction_type_prefix}_{metab}{reaction_type_suffix}": gcms_value
        for metab, gcms_value in gcms_metab_norm.items()
        if f"{reaction_type_prefix}_{metab}{reaction_type_suffix}" in model.reactions
        and gcms_value != 0.0
    }

    # Ensure values are normalized
    total = sum(rxn_id_gcms_val.values())
    normalized_gcms_metab_vals = {k: v / total for k, v in rxn_id_gcms_val.items()}

    assert sum(normalized_gcms_metab_vals.values()) == 1.0

    # Calculate the constraint for each metabolite's specified reactions
    metab_constrained_flux_expr = dict()
    for vmh_id, gcms_value in gcms_metab_norm.items():
        if (
            gcms_value != 0.0
            and f"{reaction_type_prefix}_{vmh_id}{reaction_type_suffix}"
            in model.reactions
        ):
            metab_constrained_flux_expr[
                f"{reaction_type_prefix}_{vmh_id}{reaction_type_suffix}"
            ] = gcms_value * sum(
                model.reactions.get_by_id(
                    f"{reaction_type_prefix}_{metab}{reaction_type_suffix}"
                ).flux_expression
                for metab in gcms_metab_norm.keys()
                if f"{reaction_type_prefix}_{metab}{reaction_type_suffix}"
                in model.reactions
            )

    # Add the constraints to the model for each metabolite
    if reaction_type_prefix == "UFEt" and reaction_type_suffix == "":
        lowerb = 0.0
        upperb = 1000000.0
    elif reaction_type_prefix == "EX" and reaction_type_suffix == "[fe]":
        lowerb = -1000.0
        upperb = 1000000.0
    # elif IEX:
    #     lowerb = -1000.
    #     upperb = 1000.

    constraint_list = []
    for metab_rxn, flux_expr in metab_constrained_flux_expr.items():
        constraint = model.problem.Constraint(
            flux_expr, lb=lowerb, ub=upperb, name=metab_rxn
        )
        constraint_list.append(constraint)
    model.add_cons_vars(constraint_list)
    model.solver.update()

    opt_results_dict = dict()

    for met_rxn in metab_constrained_flux_expr.keys():
        model.objective = model.reactions.get_by_id(met_rxn)
        solution = model.optimize(objective_sense=opt_direction)
        print(f"{model.objective}:\t{solution.objective_value}")
        opt_results_dict[model.objective] = solution.objective_value

    return opt_results_dict
