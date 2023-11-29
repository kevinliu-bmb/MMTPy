import time
import pandas as pd
import pubchempy as pcp
import cobra
from typing import Union
from utils.string_processing import StringProcessing
from utils.model_utils import ModelUtils
from utils.dataframe_analysis import DataframeAnalysis


class MetaboliteMatching:
    def __init__(self):
        self.string_processor = StringProcessing()
        self.dataframe_analyzer = DataframeAnalysis()

    def match_names_to_vmh(
        self,
        mbx_filepath: str,
        vmh_db_filepath: str = "data/all_vmh_metabolites.tsv",
        manual_matching_filepath: str = "data/manually_matched_keys.txt",
        silent: bool = False,
    ) -> dict:
        """Map the metabolite names detected by MBX to VMH identifiers."""
        # Load the data
        vmh_db_df = pd.read_csv(vmh_db_filepath, index_col=0, header=0, delimiter="\t")
        mbx_data_df = pd.read_csv(mbx_filepath, index_col=0, header=0)

        # Get the index of the first numerical column in mbx_data_df
        first_mbx_idx = self.dataframe_analyzer.get_init_mbx_idx(mbx_data_df)

        # Create dictionaries for direct matching
        mbx_names_dict = {
            name: self.string_processor.convert_string(name).lower()
            for name in mbx_data_df.columns[first_mbx_idx:].to_list()
        }
        vmh_names_dict = dict(zip(vmh_db_df["fullName"].index, vmh_db_df["fullName"]))

        # Perform direct matching
        direct_matching_dict = dict()
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

        # Matching by pubchempy and vmh database
        vmh_cid_dict = dict(zip(vmh_db_df["pubChemId"].index, vmh_db_df["pubChemId"]))
        vmh_inchikey_dict = dict(
            zip(vmh_db_df["inchiKey"].index, vmh_db_df["inchiKey"])
        )
        vmh_inchistring_dict = dict(
            zip(vmh_db_df["inchiString"].index, vmh_db_df["inchiString"])
        )
        vmh_smiles_dict = dict(zip(vmh_db_df["smile"].index, vmh_db_df["smile"]))

        iupac_names_dict = dict()
        cid_names_dict = dict()
        inchi_names_dict = dict()
        inchikey_names_dict = dict()
        smiles_names_dict = dict()

        # Iterate over each item in the compounds dictionary
        for mbx_name, compound in unmatched_dict.items():
            try:
                # Get the compound information from PubChem
                c = pcp.get_compounds(compound, "name")
                # If the compound was found, store its properties
                if c:
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

        pubchempy_matched_dict = dict()

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
        pubchempy_matched_dict.update(
            {value: key for key, value in direct_matching_dict.items()}
        )

        if len(pubchempy_matched_dict) != len(mbx_data_df.index):
            for vmh_id, vmh_smiles in vmh_smiles_dict.items():
                for mbx_name, pubchempy_smiles in smiles_names_dict.items():
                    if vmh_smiles != "nan" and vmh_smiles == pubchempy_smiles:
                        if not silent:
                            print(
                                f"\t\tMatched '{mbx_name}' to '{vmh_id}' using SMILES"
                            )
                        pubchempy_matched_dict[mbx_name] = vmh_id

        manual_matching = dict()
        with open(manual_matching_filepath, "r") as f:
            for line in f:
                name, vmh = line.split("\t")
                if vmh.endswith("\n"):
                    vmh = vmh[:-1]
                manual_matching[name] = vmh

        # Create and return the final matching dictionary
        max_matched_dict = dict()
        max_matched_dict.update(pubchempy_matched_dict)
        max_matched_dict.update(direct_matching_dict)

        return max_matched_dict

    def fetch_norm_sample_mbx_data(
        self,
        model_input: Union[cobra.Model, str],
        mbx_filepath: str,
        matched_mbx_names: dict,
    ) -> dict:
        """Fetch normalized sample-specific MBX data for a model."""
        if isinstance(model_input, str):
            model = ModelUtils().load_model(model_input)
        elif isinstance(model_input, cobra.Model):
            model = model_input
        else:
            raise TypeError(
                "model_input must be a cobra.Model or a filepath to a COBRApy model"
            )

        print(f"\n[Fetching MBX data for {model.name}]")

        # Read metolomics data
        mbx_data = pd.read_csv(mbx_filepath, sep=",", index_col=0)

        idx_list = []
        # Given the matched metabolite names, find the index where the metabolite is found in the MBX data
        for col_name in mbx_data.columns:
            if col_name in matched_mbx_names.keys():
                idx_list.append(mbx_data.columns.get_loc(col_name))

        sample_id = []
        # Given the model name, find which row the sample ID is in
        for index in mbx_data.index:
            if index == model.name:
                sample_id.append(index)

        if len(sample_id) > 1:
            raise ValueError(
                "Multiple sample IDs found in MBX data. Please check the data."
            )
        elif len(sample_id) == 0:
            raise ValueError("No sample ID found in MBX data. Please check the data.")
        else:
            sample_id = sample_id[0]

        # Given a sample ID, get the MBX data for that sample
        sample_mbx_data = mbx_data.loc[sample_id][min(idx_list) :]

        # Create a dictionary of metabolite names and their concentrations
        met_raw_vals_dict = {name: conc for name, conc in sample_mbx_data.items()}

        # If the values of met_raw_vals_dict are strings, remove any commas and convert them to floats
        for k, v in met_raw_vals_dict.items():
            if isinstance(v, str):
                met_raw_vals_dict[k] = float(v.replace(",", ""))

        # Check if all values in met_raw_vals_dict are floats
        if not all(isinstance(v, float) for v in met_raw_vals_dict.values()):
            raise ValueError("Not all values in met_raw_vals_dict are floats.")

        # Raise ValueError if matched_mbx_names is empty
        if len(matched_mbx_names) == 0:
            raise ValueError("The VMH Identifier name-matching dictionary is empty.")

        # Create a dictionary of VMH IDs and their corresponding metabolite values if the metabolite is in the model
        vmh_id_values = dict()
        for vmh_name, vmh_id in matched_mbx_names.items():
            for mbx_name, value in met_raw_vals_dict.items():
                if vmh_name == mbx_name and f"EX_{vmh_id}[fe]" in [
                    rxn.id for rxn in model.reactions
                ]:
                    vmh_id_values[vmh_id] = value

        # Normalize the values
        total_val = sum(vmh_id_values.values())
        norm_vmh_id_vals = {k: v / total_val for k, v in vmh_id_values.items()}

        print(
            f"\n\t[Number of name-matched and normalized metabolites in the model: {len(norm_vmh_id_vals)}/{len(met_raw_vals_dict)} ({round(len(norm_vmh_id_vals)/len(met_raw_vals_dict)*100, 2)}%)]"
        )
        print(f"\n[Returning normalized sample-specific MBX values for {sample_id}]")

        return norm_vmh_id_vals

    def fetch_mbx_constr_list(
        self, model: cobra.Model, mbx_met_norm_dict: dict
    ) -> list:
        """Compute the MBX associated FEX reaction constraints for a sample."""
        print(f"\n[Fetching MBX constraints for {model.name}]")
        # Calculate the constraint for each metabolite's fecal exchange reactions
        met_constrained_flux_expr = dict()
        for vmh_id, mbx_value in mbx_met_norm_dict.items():
            if mbx_value != 0.0 and f"EX_{vmh_id}[fe]" in model.reactions:
                flux_expr = model.reactions.get_by_id(
                    f"EX_{vmh_id}[fe]"
                ).flux_expression
                met_constrained_flux_expr[
                    f"EX_{vmh_id}[fe]"
                ] = flux_expr - mbx_value * sum(
                    model.reactions.get_by_id(f"EX_{met}[fe]").flux_expression
                    for met in mbx_met_norm_dict.keys()
                    if f"EX_{met}[fe]" in model.reactions
                )

        # Add the constraints to a list for each metabolite with MBX data
        constraint_list = []
        for ex_rxn_id, constr_expr in met_constrained_flux_expr.items():
            constraint = model.problem.Constraint(
                constr_expr, lb=0, ub=0, name=ex_rxn_id + "_constraint"
            )
            constraint_list.append(constraint)

        if len(constraint_list) == 0:
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
        self, model: cobra.Model, constraints: list, parallel: bool = False
    ) -> (list, list):
        """Add slack variables to infeasible constraints and test if the solution is feasible."""
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

        print(
            "\n[Simplifying the constraints based on the slack variable primal values]"
        )

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
                # Adjust constraint to v_EX_fecal_met_i > x * summation_j_in_set_MBX_data(v_EX_fecal_met_j)
                constraint.lb = 0
                constraint.ub = None
                refined_constraints.append(constraint)
                constraint_details = f"\t\t{constraint.name}:\t({constraint.lb}, {constraint.ub}),\t[Slacked, flux greater than constrained value]"
                constraint_details_log.append(constraint_details)
                if not parallel:
                    print(constraint_details)
            elif constraint.name in intersection_lt:
                # Adjust constraint to v_EX_fecal_met_i < x * summation_j_in_set_MBX_data(v_EX_fecal_met_j)
                constraint.lb = None
                constraint.ub = 0
                refined_constraints.append(constraint)
                constraint_details = f"\t\t{constraint.name}:\t({constraint.lb}, {constraint.ub}),\t[Slacked, flux less than constrained value]"
                constraint_details_log.append(constraint_details)
                if not parallel:
                    print(constraint_details)
            elif constraint.name in intersection_eq:
                # Adjust constraint to v_EX_fecal_met_i = x * summation_j_in_set_MBX_data(v_EX_fecal_met_j)
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
