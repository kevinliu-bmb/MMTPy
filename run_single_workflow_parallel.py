import concurrent.futures
import os

from mmtpy_opt_workflows import optimize_model, optimize_model_mbx
from mmtpy_utils import convert_model_format

# Define paths
model_path = "example_data/models"
mbx_path = "example_data/metabolomics_data.csv"
output_path = "example_outputs"

# Add 1BA to the model
add_1ba = False

# Select workflow (either "optimize_model" or "optimize_model_mbx")
workflow = "optimize_model_mbx"


def main(
    model_path: str,
    mbx_path: str,
    output_path: str,
    workflow: str,
    add_1ba: bool,
):
    # Search through the model directory and find all the JSON and MATLAB files
    model_files = [f for f in os.listdir(model_path) if f.endswith(".json")]
    model_files_mat = [f for f in os.listdir(model_path) if f.endswith(".mat")]

    # If the number of JSON files is less than the number of MATLAB files, convert
    # all models to JSON format
    if len(model_files) < len(model_files_mat):
        with concurrent.futures.ProcessPoolExecutor() as executor:
            futures = []
            for model in model_files_mat:
                future = executor.submit(
                    convert_model_format, f"{model_path}/{model}", model_path
                )
                futures.append(future)
            for future in futures:
                future.result()

    # Search through the model directory and find all the JSON files
    model_files = [f for f in os.listdir(model_path) if f.endswith(".json")]

    # Initialize ProcessPoolExecutor
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = []
        for model_file in model_files:
            if workflow == "optimize_model":
                args = [
                    f"{model_path}/{model_file}",
                    output_path,
                    add_1ba,
                    True,
                    False,
                    True,
                ]
                future = executor.submit(optimize_model, *args)
            elif workflow == "optimize_model_mbx":
                args = [
                    f"{model_path}/{model_file}",
                    mbx_path,
                    output_path,
                    add_1ba,
                    True,
                    False,
                    False,
                    True,
                ]
                future = executor.submit(optimize_model_mbx, *args)
            else:
                raise ValueError(
                    "Invalid workflow. Choose either 'optimize_model' or 'optimize_model_mbx'"
                )

            futures.append(future)

        # Retrieve the results
        for future in futures:
            future.result()


if __name__ == "__main__":
    main(model_path, mbx_path, output_path, workflow, add_1ba)
