import concurrent.futures
import os

from optimization_workflows_min import (
    optimize_model,
    optimize_model_mbx,
)
from utils import convert_model_format, match_names_to_vmh


def main(
    model_path: str,
    diet_path: str,
    output_path: str,
    workflow: str,
    mbx_path: str = None,
    mbx_matched_keys_input: str = None,
):
    # Search through the model directory and find all the JSON and MATLAB files
    model_files = [f for f in os.listdir(model_path) if f.endswith(".json")]
    model_files_mat = [f for f in os.listdir(model_path) if f.endswith(".mat")]

    # If the number of JSON files is less than the number of MATLAB files, convert all models to JSON format
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
                    diet_path,
                    output_path,
                    True,
                    False,
                    True,
                ]
                future = executor.submit(optimize_model, *args)
            elif workflow == "optimize_model_mbx":
                # If the matched keys file does not exist, use the match_names_to_vmh function
                if isinstance(mbx_matched_keys_input, str) and not os.path.exists(
                    mbx_matched_keys_input
                ):
                    mbx_matched_keys_input = match_names_to_vmh(
                        mbx_filepath=mbx_path,
                        output_filepath=output_path,
                        reuturn_matched_keys=True,
                    )
                args = [
                    f"{model_path}/{model_file}",
                    diet_path,
                    mbx_path,
                    mbx_matched_keys_input,
                    output_path,
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
    # Define paths
    model_path = "/Users/kevinliu/Dropbox (Partners HealthCare)/Zomorrodi_Lab/Capstone_Project/workflow_data_io/cobraGEMM_data/italy_models_json"
    diet_path = "/Users/kevinliu/GitHub/cobraGEMM/workflows/optimization/example_data/mediterranean_diet.txt"
    # mbx_path = "/Users/kevinliu/Dropbox (Partners HealthCare)/Zomorrodi_Lab/Capstone_Project/workflow_data_io/cobraGEMM_data/metabolomics_data.csv"
    output_path = "/Users/kevinliu/Dropbox (Partners HealthCare)/Zomorrodi_Lab/Capstone_Project/workflow_data_io/cobraGEMM_data/output_min_italy_basic"
    # mbx_matched_keys_input = "/Users/kevinliu/Dropbox (Partners HealthCare)/Zomorrodi_Lab/Capstone_Project/workflow_data_io/cobraGEMM_data/output_maximize_us/metabolomics_data_matched_key.txt"

    # Select workflow (either "optimize_model" or "optimize_model_mbx")
    workflow = "optimize_model"

    # Run the main function
    main(model_path, diet_path, output_path, workflow)
