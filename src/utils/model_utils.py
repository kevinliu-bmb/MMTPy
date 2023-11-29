import cobra
import os
import re
from math import isnan


class ModelUtils:
    model_name_search_pattern = re.compile(r"(Case|Control)_\d+")

    def load_model(
        self, model_path: str, simple_model_name: bool = True
    ) -> cobra.Model:
        """Load a multi-species model into memory given a path to a model file in a COBRApy supported format.

        Parameters:
        model_path : str
            Path to a multi-species model in any COBRApy supported format.
        simple_model_name : bool
            If True, the model name will be set to the file name with only the sample identifier (e.g., Case_1) and not the full file name (e.g., Case_1.xml).

        Returns:
        cobra.Model
            A COBRApy model loaded into memory.

        Raises:
        ValueError
            If the model_path does not exist.
        ValueError
            If the model format is not supported.
        """

        if not os.path.exists(model_path):
            raise ValueError(f"Model path does not exist: {model_path}")

        print(f"\n[Loading model from {model_path}]")

        if model_path.endswith(".xml"):
            model = cobra.io.read_sbml_model(model_path)
        elif model_path.endswith(".json"):
            model = cobra.io.load_json_model(model_path)
        elif model_path.endswith(".yml"):
            model = cobra.io.load_yaml_model(model_path)
        elif model_path.endswith(".mat"):
            model = cobra.io.load_matlab_model(model_path)
        elif model_path.endswith(".sbml"):
            model = cobra.io.read_sbml_model(model_path)
        else:
            raise ValueError(f"Model format not supported for {model_path}")

        if simple_model_name:
            model_file_name = os.path.basename(model_path).split(".")[0]
            model.name = self.model_name_search_pattern.search(model_file_name).group(0)
        else:
            model.name = os.path.basename(model_path).split(".")[0]

        print(f"\n[{model.name} loaded]")

        return model

    def set_default_bounds(
        self,
        model: cobra.Model,
        source: str = "cobraGEMM",
        rxn_type: str = "all",
        silent: bool = False,
    ) -> bool:
        """Set the bounds of the model's reactions according to conventions; prints the changes and returns True if the bounds were different from the default state. Conventional bounds can either be set based on Heinken et al. (2022), mgPipe models or be set based on cobraGEMM conventions.

        Parameters:
        model : cobra.Model
            The model whose reactions' bounds are to be set.
        source : str, optional
            The definition of conventional bounds, by default "cobraGEMM"; options are of either "cobraGEMM" or "MATLAB".
        rxn_type : str, optional
            The type of reactions whose bounds are to be set, by default "all";
            options are of either "all", "FEX", "UFEt", "IEX", "DUt", or "commBiomass".
        silent : bool, optional
            Whether to print the changes, by default False.

        Returns:
        bool
            True if the bounds were different from the default state.

        Notes
        -----
        The conventions are as follows based on Heinken et al. (2022), mgPipe models:
        1. Set the bounds of the fecal exchange (EX_met[fe]) reactions for metabolites to be (-1000., 1000000.)
        2. Set the bounds of the fecal exchange (EX_met[fe]) reaction for "microbeBiomass" to be (-10000., 1000000.)
        3. Set the bounds of the fecal transport (UFEt_met) reactions to be (0., 1000000.)
        4. Set the bounds of the microbe secretion/uptake (microbe_IEX_met[u]tr) reactions to be (-1000., 1000.)
        5. Set the bounds of the community biomass reaction to be (0.4, 1.)

        The default bounds for cobraGEMM are as follows:
        1. Set the bounds of the fecal exchange (EX_met[fe]) reactions for metabolites to be (0., 1000000.)
        2. Set the bounds of the fecal exchange (EX_met[fe]) reaction for "microbeBiomass" to be (0., 1000000.)
        3. Set the bounds of the fecal transport (UFEt_met) reactions to be (0., 1000000.)
        4. Set the bounds of the microbe secretion/uptake (microbe_IEX_met[u]tr) reactions to be (-1000., 1000.)
        5. Set the bounds of the community biomass reaction to be (0.4, 1.)
        """

        print(f"\n[Setting bounds for {rxn_type} reactions using {source} defaults]")

        saved_bounds = dict()
        new_bounds = dict()
        for rxn in model.reactions:
            # Set the bounds of the fecal exchange (EX_met[fe]) reactions to be (-1000., 1000000.)
            if (
                rxn.id.startswith("EX_")
                and rxn.id.endswith("[fe]")
                and "microbeBiomass" not in rxn.id
                and rxn_type in ["all", "FEX"]
            ):
                saved_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
                if source == "cobraGEMM":
                    model.reactions.get_by_id(rxn.id).bounds = (0.0, 1000000.0)
                elif source == "MATLAB":
                    model.reactions.get_by_id(rxn.id).bounds = (-1000.0, 1000000.0)
                else:
                    raise ValueError(f"Source {source} not supported")
                new_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
            # Set the bounds of the fecal exchange (EX_met[fe]) reactions for the microbeBiomass to be (-10000., 1000000.)
            elif (
                rxn.id.startswith("EX_")
                and rxn.id.endswith("[fe]")
                and "microbeBiomass" in rxn.id
                and rxn_type in ["all", "FEX"]
            ):
                saved_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
                if source == "cobraGEMM":
                    model.reactions.get_by_id(rxn.id).bounds = (0.0, 1000000.0)
                elif source == "MATLAB":
                    model.reactions.get_by_id(rxn.id).bounds = (-10000.0, 1000000.0)
                else:
                    raise ValueError(f"Source {source} not supported")
                new_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
            # Set the bounds of the fecal transport (UFEt_met) reactions to be (0., 1000000.)
            elif rxn.id.startswith("UFEt_") and rxn_type in ["all", "UFEt"]:
                saved_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
                model.reactions.get_by_id(rxn.id).bounds = (0.0, 1000000.0)
                new_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
            # Set the bounds of the microbe secretion/uptake (microbe_IEX_met[u]tr) reactions to be (-1000., 1000.)
            elif (
                "IEX" in rxn.id
                and rxn.id.endswith("[u]tr")
                and rxn_type in ["all", "IEX"]
            ):
                saved_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
                model.reactions.get_by_id(rxn.id).bounds = (-1000.0, 1000.0)
                new_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
            # Set the bounds of the diet transport (DUt_met) reactions to be (0., 1000000.)
            elif rxn.id.startswith("DUt_") and rxn_type in ["all", "DUt"]:
                saved_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
                model.reactions.get_by_id(rxn.id).bounds = (0.0, 1000000.0)
                new_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
            # Set the bounds of the community biomass reaction to be (0.4, 1.)
            elif rxn.id == "communityBiomass" and rxn_type in ["all", "commBiomass"]:
                saved_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds
                model.reactions.get_by_id(rxn.id).bounds = (0.4, 1.0)
                new_bounds[rxn.id] = model.reactions.get_by_id(rxn.id).bounds

        n_changed_bounds = 0
        # Print out the changes
        for rxn, bounds in saved_bounds.items():
            if bounds != new_bounds[rxn]:
                if not silent:
                    print(
                        f"\tChanged bounds for {rxn} from {bounds} to {new_bounds[rxn]}"
                    )
                n_changed_bounds += 1

        if n_changed_bounds == 0:
            bounds_changed = False
            print("\tNo bounds were changed")
        else:
            bounds_changed = True
            print(
                f"\tChanged {n_changed_bounds}/{len(saved_bounds)} {rxn_type} reaction bounds for {model.name}"
            )

        return bounds_changed

    def convert_model_format(
        self, model_path: str, output_path: str, output_format: str = "json"
    ) -> None:
        if isinstance(model_path, str):
            model = self.load_model(model_path=model_path, simple_model_name=False)
        elif isinstance(model_path, cobra.Model):
            model = model_path
        else:
            raise ValueError(
                f"model_path must be either a file path (str) or a COBRApy model, not {type(model_path)}"
            )

        print(f"\n[Converting model {model.name} to json format]")

        # Convert the metabolite charge to a string if it is NaN
        for met in model.metabolites:
            if isnan(met.charge):
                met.charge = "nan"

        if not os.path.exists(output_path):
            os.makedirs(output_path)

        if output_path.endswith("/"):
            converted_output_filepath = f"{output_path}{model.name}.json"
        else:
            converted_output_filepath = f"{output_path}/{model.name}.json"

        cobra.io.save_json_model(model, converted_output_filepath)

        print(f"\n[{model.name} converted to json format]")
