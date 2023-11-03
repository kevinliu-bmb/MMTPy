import os
import re
import tempfile

import cobra


def load_model(model_path: str, simple_model_name: bool = True) -> cobra.Model:
    """
    The function `load_model` loads a metabolic model from a specified file path
    and assigns a name to the model based on the file name.

    :param model_path: The `model_path` parameter is a string that specifies the
    path to the file containing the model. This can be an absolute path or a
    relative path
    :type model_path: str
    :param simple_model_name: The `simple_model_name` parameter is a boolean flag
    that determines whether the model name should be simplified or not. If
    `simple_model_name` is set to `True`, the model name will be extracted from the
    file name using a regular expression pattern. The regular expression pattern
    assumes that the file name, defaults to True
    :type simple_model_name: bool (optional)
    :return: a cobra.Model object.
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
        raise TypeError(f"Model format not supported for {model_path}")

    if simple_model_name:
        model_file_name = os.path.basename(model_path).split(".")[0]
        model_name_search_pattern = re.compile(r"(Case|Control)_\d+")
        model.name = model_name_search_pattern.search(model_file_name)[0]
    else:
        model.name = os.path.basename(model_path).split(".")[0]

    print(f"\n[{model.name} loaded]")

    return model


def validate_model(model: cobra.Model):
    """
    The function `validate_model` takes a cobra.Model object, writes it to a
    temporary file in SBML format, and then validates the SBML model.

    :param model: The `model` parameter is of type `cobra.Model`. This is a model
    object from the COBRApy package, which is used for constraint-based modeling of
    biological systems. It represents a mathematical representation of a biological
    system, typically a metabolic network
    :type model: cobra.Model
    """
    with tempfile.NamedTemporaryFile() as temp_file:
        cobra.io.write_sbml_model(model, temp_file.name)
        cobra.io.validate_sbml_model(temp_file.name)
