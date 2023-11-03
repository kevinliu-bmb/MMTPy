import concurrent.futures
import os
from match_names_to_vmh import match_names_to_vmh

from optimization_workflows import optimize_model, optimize_model_mbx
from utils import convert_model_format


def main(
    model_path: str,
    mbx_path: str,
    mbx_matched_keys_input: str,
    output_path: str,
):
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
        # If the matched keys file does not exist, use the match_names_to_vmh function
        if isinstance(mbx_matched_keys_input, str) and not os.path.exists(
            mbx_matched_keys_input
        ):
            mbx_matched_keys_input = match_names_to_vmh(
                mbx_filepath=mbx_path,
                output_filepath=output_path,
                reuturn_matched_keys=True,
            )

        # Arguments for the functions
        args_optimize_model = [
            f"{model_path}/{model_file}",
            output_path,
            True,
            False,
            True,
        ]
        args_optimize_model_mbx = [
            f"{model_path}/{model_file}",
            mbx_path,
            mbx_matched_keys_input,
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
    # Define paths
    model_path = "workflows/optimization/example_data/models"
    mbx_path = "workflows/optimization/example_data/metabolomics_data_assert_non_large_value.csv"
    output_path = "workflows/optimization/example_outputs"
    mbx_matched_keys_input = (
        "workflows/optimization/example_outputs/metabolomics_data_matched_key.txt"
    )

    # Run the main function
    main(model_path, mbx_path, mbx_matched_keys_input, output_path)
