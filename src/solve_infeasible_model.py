import pandas as pd
import cobra

from essential_metabolites import essential_metabolites


def adapt_diet_and_minimize_infeasibility(
    model: cobra.Model, diet_file_path: str, output_file: str
) -> cobra.Model:
    # read initial diet from the file
    diet_df = pd.read_csv(diet_file_path, sep="\t")
    original_diet = {
        row["Reaction"]: -abs(row["Flux Value"]) for _, row in diet_df.iterrows()
    }

    # add essential metabolites that are missing in the original diet
    for essential_met in essential_metabolites:
        diet_rxn_id = f"EX_{essential_met}[d]"
        if diet_rxn_id not in original_diet:
            original_diet[diet_rxn_id] = -0.1  # default value

    # apply the original diet to the model
    for rxn_id, flux_value in original_diet.items():
        if rxn_id in model.reactions:
            model.reactions.get_by_id(rxn_id).bounds = (flux_value, 0)

    initial_solution = set_commBiomass_obj(model)
    feasible_diet = {}
    if initial_solution.status == "optimal":
        print(
            "\n[DIET SOLVER] Model is feasible with given diet and essential metabolites."
        )
        for rxn_id, flux_value in original_diet.items():
            feasible_diet[rxn_id] = flux_value
    else:
        check_commBiomass_feasibility(original_diet, model, feasible_diet)
    # write adjusted diet to output file
    output_data = []
    for rxn_id, original_flux in original_diet.items():
        feasible_flux = feasible_diet.get(rxn_id, "N/A")
        is_essential = (
            "Yes"
            if rxn_id in [f"EX_{met}[d]" for met in essential_metabolites]
            else "No"
        )

        # Indicate if the metabolite is from the original diet or added
        from_original_diet = "Yes" if rxn_id in diet_df["Reaction"].values else "No"

        # Calculate the numerical difference
        if feasible_flux != "N/A":
            flux_difference = original_flux - feasible_flux
        else:
            flux_difference = "N/A"

        output_data.append(
            {
                "Reaction": rxn_id,
                "Flux Value": original_flux,
                "Flux Value Feasible": feasible_flux,
                "Is Essential": is_essential,
                "From Original Diet": from_original_diet,
                "Flux Difference": flux_difference,
            }
        )

    print("\n[DIET SOLVER] Writing adjusted diet to output file and returning model.")
    output_df = pd.DataFrame(output_data)
    output_df.to_csv(output_file, sep="\t", index=False)

    return model


def check_commBiomass_feasibility(original_diet, model, feasible_diet):
    # logic to adjust diet and re-check for feasibility
    print(
        "\n[DIET SOLVER] Model is infeasible with given diet and essential metabolites, relaxing bounds by 0.1."
    )
    for rxn_id in original_diet:
        if rxn_id in model.reactions:
            current_lower_bound = model.reactions.get_by_id(rxn_id).lower_bound
            model.reactions.get_by_id(rxn_id).bounds = (
                current_lower_bound - 0.1,
                0,
            )  # decrease lower bound by 0.1

    adjusted_solution = set_commBiomass_obj(model)
    if adjusted_solution.status != "optimal":
        raise ValueError(
            "\n[DIET SOLVER] Model is infeasible with relaxed diet and essential metabolites."
        )

    print(
        "\n[DIET SOLVER] Model is feasible with relaxed diet and essential metabolites."
    )
    for rxn_id in original_diet:
        feasible_diet[rxn_id] = model.reactions.get_by_id(rxn_id).lower_bound


def set_commBiomass_obj(model):
    # re-check for feasibility
    # set the communityBiomass reaction bounds to (0.4, 1)
    model.reactions.get_by_id("communityBiomass").bounds = (0.4, 1)
    # set the model objective to minimize the communityBiomass reaction
    model.objective = model.reactions.get_by_id("communityBiomass")
    return model.optimize()
