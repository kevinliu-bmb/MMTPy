import pandas as pd
import cobra


essential_metabolites = [
    '12dgr180', '26dap_M', '2dmmq8', '2obut', '3mop', '4abz', '4hbz', 'ac',
    'acgam', 'acmana', 'acnam', 'ade', 'adn', 'adocbl', 'ala_D', 'ala_L',
    'amet', 'amp', 'arab_D', 'arab_L', 'arg_L', 'asn_L', 'btn', 'ca2', 
    'cbl1', 'cgly', 'chor', 'chsterol', 'cit', 'cl', 'cobalt2', 'csn',
    'cu2', 'cys_L', 'cytd', 'dad_2', 'dcyt', 'ddca', 'dgsn', 'fald', 'fe2',
    'fe3', 'fol', 'for', 'gal', 'glc_D', 'gln_L', 'glu_L', 'gly', 'glyc',
    'glyc3p', 'gsn', 'gthox', 'gthrd', 'gua', 'h', 'h2o', 'h2s', 'his_L',
    'hxan', 'ile_L', 'k', 'lanost', 'leu_L', 'lys_L', 'malt', 'met_L', 'mg2',
    'mn2', 'mqn7', 'mqn8', 'nac', 'ncam', 'nmn', 'no2', 'ocdca', 'ocdcea',
    'orn', 'phe_L', 'pheme', 'pi', 'pnto_R', 'pro_L', 'ptrc', 'pydx', 'pydxn',
    'q8', 'rib_D', 'ribflv', 'ser_L', 'sheme', 'so4', 'spmd', 'thm', 'thr_L',
    'thymd', 'trp_L', 'ttdca', 'tyr_L', 'ura', 'val_L', 'xan', 'xyl_D', 'zn2'
]


def adapt_diet_and_minimize_infeasibility(model: cobra.Model, diet_file_path: str, output_file: str) -> cobra.Model:
    """Incorporate AGORA model-specific dietary components and check for model feasibility, reelaxing the lower bounds if necessary.

    Parameters:
    model: cobra.Model
        The model to be checked for feasibility.
    diet_file_path: str
        Path to the file containing the original diet.
    output_file: str
        Path to the file where the adjusted diet will be written.
        
    Returns:
    model: cobra.Model
        The model with the adjusted diet.
        
    Raises:
    ValueError: If the model is infeasible with the adjusted diet.
    """
    
    # read initial diet from the file
    diet_df = pd.read_csv(diet_file_path, sep='\t')
    original_diet = {row['Reaction']: -abs(row['Flux Value']) for _, row in diet_df.iterrows()}
    
    # add essential metabolites that are missing in the original diet
    for essential_met in essential_metabolites:
        diet_rxn_id = f"EX_{essential_met}[d]"
        if diet_rxn_id not in original_diet:
            original_diet[diet_rxn_id] = -0.1  # default value
            
    # apply the original diet to the model
    for rxn_id, flux_value in original_diet.items():
        if rxn_id in model.reactions:
            model.reactions.get_by_id(rxn_id).bounds = (flux_value, 0)
    
    # check for initial feasibility
    # set the communityBiomass reaction bounds to (0.4, 1)
    model.reactions.get_by_id('communityBiomass').bounds = (0.4, 1)
    # set the model objective to minimize the communityBiomass reaction
    model.objective = model.reactions.get_by_id('communityBiomass')
    initial_solution = model.optimize()
    feasible_diet = {}
    if initial_solution.status == 'optimal':
        print("\n[DIET SOLVER] Model is feasible with given diet and essential metabolites.")
        for rxn_id, flux_value in original_diet.items():
            feasible_diet[rxn_id] = flux_value
    else:
        # logic to adjust diet and re-check for feasibility
        print("\n[DIET SOLVER] Model is infeasible with given diet and essential metabolites, relaxing bounds by 0.1.")
        for rxn_id in original_diet.keys():
            if rxn_id in model.reactions:
                current_lower_bound = model.reactions.get_by_id(rxn_id).lower_bound
                model.reactions.get_by_id(rxn_id).bounds = (current_lower_bound - 0.1, 0) # decrease lower bound by 0.1
        
        # re-check for feasibility
        # set the communityBiomass reaction bounds to (0.4, 1)
        model.reactions.get_by_id('communityBiomass').bounds = (0.4, 1)
        # set the model objective to minimize the communityBiomass reaction
        model.objective = model.reactions.get_by_id('communityBiomass')
        adjusted_solution = model.optimize()
        if adjusted_solution.status == 'optimal':
            print("\n[DIET SOLVER] Model is feasible with relaxed diet and essential metabolites.")
            for rxn_id, flux_value in original_diet.items():
                feasible_diet[rxn_id] = model.reactions.get_by_id(rxn_id).lower_bound
        else:
            raise ValueError("\n[DIET SOLVER] Model is infeasible with relaxed diet and essential metabolites.")
            
    # write adjusted diet to output file
    output_data = []
    for rxn_id, original_flux in original_diet.items():
        feasible_flux = feasible_diet.get(rxn_id, "N/A")
        is_essential = "Yes" if rxn_id in [f"EX_{met}[d]" for met in essential_metabolites] else "No"
        
        # Indicate if the metabolite is from the original diet or added
        from_original_diet = "Yes" if rxn_id in diet_df['Reaction'].values else "No"

        # Calculate the numerical difference
        if feasible_flux != "N/A":
            flux_difference = original_flux - feasible_flux
        else:
            flux_difference = "N/A"
        
        output_data.append({
            "Reaction": rxn_id,
            "Flux Value": original_flux,
            "Flux Value Feasible": feasible_flux,
            "Is Essential": is_essential,
            "From Original Diet": from_original_diet,
            "Flux Difference": flux_difference
        })
        
    print("\n[DIET SOLVER] Writing adjusted diet to output file and returning model.")
    output_df = pd.DataFrame(output_data)
    output_df.to_csv(output_file, sep='\t', index=False)
    
    return model
