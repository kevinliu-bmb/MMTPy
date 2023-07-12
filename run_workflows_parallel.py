import concurrent.futures
import os

from mmtpy_opt_workflows import optimize_model, optimize_model_mbx
from mmtpy_utils import convert_model_format

# Define paths
model_path = "example_data/models"
mbx_path = "example_data/metabolomics_data_assert_non_large_value.csv"
output_path = "example_outputs"


def main(model_path, mbx_path, output_path):
    # Search through the model directory and find all the JSON and MATLAB files
    model_files = [f for f in os.listdir(model_path) if f.endswith(".json")]
    model_files_mat = [f for f in os.listdir(model_path) if f.endswith(".mat")]

    # If the number of JSON files is less than the number of MATLAB files, convert
    # all models to JSON format
    if len(model_files) < len(model_files_mat):
        for model in model_files_mat:
            convert_model_format(f"{model_path}/{model}", model_path)

    # Search through the model directory and find all the JSON files
    model_files = [f for f in os.listdir(model_path) if f.endswith(".json")]

    for model_file in model_files:
        # Arguments for the functions
        args_optimize_model = [
            f"{model_path}/{model_file}",
            output_path,
            False,
            True,
            False,
            True,
        ]
        args_optimize_model_mbx = [
            f"{model_path}/{model_file}",
            mbx_path,
            output_path,
            True,
            False,
            False,
            True,
        ]

        # Initialize ProcessPoolExecutor
        with concurrent.futures.ProcessPoolExecutor() as executor:
            # Submit the functions to the pool of workers
            future_optimize_model = executor.submit(
                optimize_model, *args_optimize_model
            )
            future_optimize_model_mbx = executor.submit(
                optimize_model_mbx, *args_optimize_model_mbx
            )

            # Retrieve the results
            future_optimize_model.result()
            future_optimize_model_mbx.result()


if __name__ == "__main__":
    main(model_path, mbx_path, output_path)
