import concurrent.futures
import os
import argparse

from optimization_workflows import optimize_model, optimize_model_mbx
from utils import convert_model_format, match_names_to_vmh


def main(
    model_path: str,
    diet_path: str,
    mbx_path: str,
    mbx_matched_keys_input: str,
    output_path: str,
    workflow: str,
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
                        return_matched_keys=True,
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
    parser = argparse.ArgumentParser(
        description="Run a single workflow in parallel for all models in a directory."
    )
    parser.add_argument(
        "--model_path",
        type=str,
        help="Path to the directory containing the model files.",
    )
    parser.add_argument(
        "--diet_path",
        type=str,
        help="Path to the diet file.",
    )
    parser.add_argument(
        "--mbx_path",
        type=str,
        help="Path to the MBX file.",
    )
    parser.add_argument(
        "--mbx_matched_keys_input",
        type=str,
        help="Path to the matched keys file.",
    )
    parser.add_argument(
        "--output_path",
        type=str,
        help="Path to the output directory.",
    )
    parser.add_argument(
        "--workflow",
        type=str,
        help="Name of the workflow to run.",
    )
    parser.add_argument(
        "--n_workers",
        type=int,
        default=1,
        help="Number of workers to use.",
    )
    parser.add_argument(
        "--help",
        action="help",
        help="Show this help message and exit.",
    )

    args = parser.parse_args()

    main(
        model_path=args.model_path,
        diet_path=args.diet_path,
        mbx_path=args.mbx_path,
        mbx_matched_keys_input=args.mbx_matched_keys_input,
        output_path=args.output_path,
        workflow=args.workflow,
    )