from xopen import xopen
import time
import os
import math
import sys
import json
import signal
import socket
import traceback
import random
from pathlib import Path
from typing import Tuple
import errno


import warnings

warnings.filterwarnings("ignore", category=FutureWarning, module="pandas")

# COMPADRE imports
import ersa
from pop_classifier import run_new as run_pop_classifier

########################################


def signal_handler(signum, frame):
    try:
        if "server_socket" in globals() and server_socket:
            server_socket.close()
    except:
        pass
    os._exit(0)


def safe_print(*args, **kwargs):
    try:
        print(*args, flush=True, **kwargs)
    except BrokenPipeError:
        # Python flushes standard streams on exit, so redirect remaining output to devnull to avoid another BrokenPipeError at shutdown
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, sys.stdout.fileno())
        sys.exit(0)  # Python exits with error code 1 on EPIPE


def convert_value(value):
    """Convert string values of ERSA options to appropriate Python types"""
    # Try int first
    try:
        return int(value)
    except ValueError:
        pass

    # Try float
    try:
        return float(value)
    except ValueError:
        pass

    # Handle empty string as None for file paths
    if value == "None":
        return None

    # Return as string
    return value


def detect_min_cm(segment_dict):
    """
    Detect the minimum cM value from the loaded segments
    """
    if not segment_dict:
        return None

    min_cm_found = float("inf")
    for pair_segments in segment_dict.values():
        for segment in pair_segments:
            cm_value = segment[3]  # cM is at index 3
            if cm_value < min_cm_found:
                min_cm_found = cm_value

    return min_cm_found if min_cm_found != float("inf") else None


def update_min_cm(additional_options, segment_dict):
    """
    Update min_cm if it's still the default (2.5) and all segments are above that threshold

    This is to reflect the "true" cutoff that was used in segment detection but might not have been passed into COMPADRE via the ERSA min_cm flag
    """
    current_min_cm = additional_options["min_cm"]

    # Only auto-detect if still using the default value
    if current_min_cm != 2.5:
        safe_print(
            f"INFO: Current min_cm is {current_min_cm} and will be passed as such",
            file=sys.stderr,
        )
        return additional_options

    detected_min_cm = detect_min_cm(segment_dict)

    if detected_min_cm is not None and detected_min_cm > 2.5:
        safe_print(
            f"INFO: [COMPADRE] Auto-detected minimum cM: {detected_min_cm}",
            file=sys.stderr,
        )
        safe_print(
            "INFO: [COMPADRE] All segments are above default threshold (2.5)",
            file=sys.stderr,
        )
        safe_print(
            f"INFO: [COMPADRE] Updating min_cm for ERSA to: {detected_min_cm}",
            file=sys.stderr,
        )
        additional_options["min_cm"] = detected_min_cm

    return additional_options


def calculate_ersa_props(model_df):

    if model_df.empty:
        return 0, 0, 0, 1

    model_df["degree_of_relatedness"] = model_df["degree_of_relatedness"].astype(int)
    model_df["maxlnl2"] = model_df["maxlnl"].apply(lambda x: math.exp(x))

    model_df_2d = model_df[model_df["degree_of_relatedness"] == 2]
    model_df_2d_sum = model_df_2d["maxlnl2"].sum()
    model_df_3d = model_df[model_df["degree_of_relatedness"] == 3]
    model_df_3d_sum = model_df_3d["maxlnl2"].sum()

    # Newest version
    model_df_4d = model_df[
        (model_df["degree_of_relatedness"] >= 4)
        & (model_df["degree_of_relatedness"] < 8)
    ]
    model_df_4d_sum = model_df_4d["maxlnl2"].sum()
    model_df_un = model_df[
        (model_df["degree_of_relatedness"] >= 8)
        & (model_df["degree_of_relatedness"] <= 40)
    ]
    model_df_un_sum = model_df_un["maxlnl2"].sum()

    # re-proportion them
    total_prop = model_df_2d_sum + model_df_3d_sum + model_df_4d_sum + model_df_un_sum
    model_df_2d_prop = model_df_2d_sum / total_prop
    model_df_3d_prop = model_df_3d_sum / total_prop
    model_df_4d_prop = model_df_4d_sum / total_prop
    model_df_un_prop = model_df_un_sum / total_prop

    return model_df_2d_prop, model_df_3d_prop, model_df_4d_prop, model_df_un_prop


def clean_ersa_options(ersa_option_dict):

    # Check data types and revert to default values if they're not correct

    int_options = {
        "ascertained_position": -1,
        "mask_region_sim_count": 0,
        "mask_region_cross_length": 1000000,
        "mask_region_threshold": 4,
        "number_of_chromosomes": 22,
        "max_meioses": 40,
    }
    float_options = {
        "parent_offspring_zscore": 2.33,
        "pois_mean": 13.73,
        "exp_mean": 3.197036753,
        "rec_per_meioses": 35.2548101,
        "min_cm": 2.5,
        "confidence_level": 0.9,
        "max_cm": 10,
    }

    return {
        k: (
            v
            if isinstance(v, (int, float))
            else (
                int_options.get(k, v) if k in int_options else float_options.get(k, v)
            )
        )
        for k, v in ersa_option_dict.items()
    }


def create_socket() -> socket.socket:
    """function to wrap the creation logic for the socket"""
    created_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    created_socket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    created_socket.setsockopt(socket.SOL_SOCKET, socket.SO_KEEPALIVE, 1)
    created_socket.settimeout(None)
    return created_socket


PortResponse = Tuple[socket.socket | None, int | None, Exception | None]


def start_server(
    preferred_port: int,
    host="localhost",
    current_attempt: int = 0,
    max_attempts: int = 100,
) -> PortResponse:
    """finds an open port. First test the port provided by the user and then attempt to open:wrap
    a socket for another port"""

    if current_attempt == max_attempts:
        return (
            None,
            None,
            Exception(
                f"FATAL: Unable to find an open port after {max_attempts} attempts"
            ),
        )

    new_socket = create_socket()

    try:
        new_socket.bind((host, preferred_port))
        new_socket.listen(1)
        return new_socket, preferred_port, None
    except OSError as e:
        if e.errno == errno.EADDRINUSE:
            new_socket.close()
            safe_print(
                f"INFO: Port {preferred_port} is in use, searching for available port...",
                file=sys.stderr,
            )
            port = random.randint(4001, 8000)  # Try ports in this range
            return start_server(port, host, current_attempt + 1, max_attempts)
        else:
            return None, None, e


def parse_ersa_options(ersa_flag_str: str) -> dict:

    additional_options = {}

    # First step -- read in ERSA options from COMPADRE runtime flags
    ersa_flags = ersa_flag_str.strip('"').split("|")
    flag_names = ersa_flags[::2]
    flag_values = ersa_flags[1::2]

    for flag, value in zip(flag_names, flag_values):
        additional_options[flag] = convert_value(value)

    additional_options = clean_ersa_options(additional_options)

    return additional_options


########################################
def load_segment_information(
    segment_data_file: str, min_cm_options: int
) -> Tuple[dict, str]:
    """read in the shared IBD segmnet information into a dictionary that hte program can use

    Parameters
    ----------
    segment_data_file : str
        path to the shared IBD segments file

    min_cm_options : int
        threshold provided by the perl program to use as the
        mininimum segment threshold for ERSA

    Returns
    -------
    Tuple[dict, str]
        returns a tuple where the first element is the dictionary of
        shared pairwise IBD segments and the second element is whether
        or not there is IBD2 segment information in the file
    """
    segment_dict = {}
    ibd2_status = "false"

    if segment_data_file != "NA":
        safe_print(
            f"INFO: Reading in the pairwise IBD segment information from: {segment_data_file}",
            file=sys.stdout,
        )

        with xopen(
            segment_data_file, "r"
        ) as f:  # Open the large file here and populate dictionary that stays in system memory

            # Check number of columns to determine if it has IBD1/2 data or not
            header_line = f.readline().strip()
            num_columns = len(header_line.split())

            if (
                num_columns == 7
            ):  # File with segment-by-segment IBD status in 7th column
                ibd2_status = "true"

                if "chrom" in header_line:  # skip header if it exists
                    next(f)  # skip header line in this file version

                for line in f:
                    ls = line.strip().split("\t")
                    iid1, iid2, start, end, cmlen, chrom, ibd = (
                        ls[0],
                        ls[1],
                        int(ls[2]),
                        int(ls[3]),
                        round(float(ls[4]), 2),
                        int(ls[5]),
                        int(ls[6].strip()),
                    )
                    # If we move this block up here then we only have to make the tuple when it is needed
                    if cmlen >= min_cm_options:
                        # lets create a tuple of sorted ids that will make sure they are in the same order
                        key = ":".join(tuple(sorted([iid1, iid2])))

                        value = (chrom, start, end, cmlen, ibd)

                        segment_dict.setdefault(key, []).append(value)

            elif num_columns == 6:  # File without segment-by-segment IBD status

                if "chrom" in header_line:  # skip header if it exists
                    next(f)

                for line in f:
                    ls = line.split("\t")
                    if len(ls) == 11:  # germline1
                        iid1, iid2, start, end, cmlen, chrom = (
                            ls[0].split(" ")[0],
                            ls[1].split(" ")[0],
                            int(ls[3].split(" ")[0]),
                            int(ls[3].split(" ")[1]),
                            round(float(ls[6]), 2),
                            int(ls[2]),
                        )
                    else:  # germline2
                        iid1, iid2, start, end, cmlen, chrom = (
                            ls[0],
                            ls[1],
                            int(ls[2]),
                            int(ls[3]),
                            round(float(ls[4]), 2),
                            int(ls[5].strip()),
                        )
                    if cmlen >= min_cm_options:
                        key = ":".join(tuple(sorted([iid1, iid2])))
                        value = (
                            chrom,
                            start,
                            end,
                            cmlen,
                            "NA",
                        )  # no ibd1/2 data in this input
                        segment_dict.setdefault(key, []).append(value)

            else:
                safe_print(
                    "FATAL: [COMPADRE] Unrecognized segment file format. Please refer to the README (https://github.com/belowlab/compadre) for formatting guidelines."
                )
    return segment_dict, ibd2_status
    # Figure out the true min_cm to be passed into ERSA based on the contents of segment_dict


def main(
    segment_dict,
    additional_options: dict,
    ibd2_status: str,
    portnumber: int,
    output_directory: str,
):

    signal.signal(signal.SIGINT, signal_handler)  # CTRL+C
    signal.signal(signal.SIGTERM, signal_handler)  # Termination
    signal.signal(signal.SIGPIPE, signal_handler)  # Broken pipe
    global server_socket

    #########################

    socket_host = (
        os.environ["COMPADRE_HOST"] if "COMPADRE_HOST" in os.environ else "localhost"
    )

    # Find an available port
    server_socket, actual_port, error = start_server(portnumber, socket_host)

    if error:
        safe_print(
            f"FATAL: Encountered the following error while trying to start the compadre.py server: {str(error)}",
            file=sys.stderr,
        )
        raise error

    # If port changed, notify via stderr first (before stdout message)
    if actual_port != portnumber:
        safe_print(f"Port changed: {actual_port}", file=sys.stderr)

    safe_print(f"COMPADRE helper is ready on port {actual_port}")
    sys.stdout.flush()

    # The perl script forks a child process at this point so that it can
    # continue to read from the buffer between python and perl without losing
    # messages. This process is fast but might not be fast enough to get
    # every message. We are going to pause quickly to make sure the perl
    # program syncs with the python code
    time.sleep(0.01)

    if server_socket:
        try:
            while True:
                conn, _ = server_socket.accept()
                conn.settimeout(None)

                # NEW: Handle the client connection with better error handling
                should_shutdown = handle_client_connection(
                    conn,
                    segment_dict,
                    additional_options,
                    ibd2_status,
                    segment_data_file,
                    output_directory,
                )

                if should_shutdown:
                    safe_print(
                        "INFO: received shutdown signal from the handle_client_connection function",
                        file=sys.stdout,
                    )
                    break

        finally:
            server_socket.close()
            safe_print("PROGINFO: [COMPADRE] COMPADRE helper shutdown complete.")


def handle_client_connection(
    conn,
    segment_dict,
    additional_options,
    ibd2_status,
    segment_data_file,
    output_directory,
):
    try:
        with conn.makefile("r") as conn_input:
            for line in conn_input:

                # we need to check to see if the value sent to the server indicates that we should shut it down
                if line.strip() == "close":
                    response = json.dumps(
                        {"status": "success", "result": "", "message": "Closing server"}
                    )
                    conn.send((response + "\n").encode())
                    return True

                if line.strip() == "test":
                    response = json.dumps(
                        {"status": "success", "result": "", "message": "OK"}
                    )
                    conn.send((response + "\n").encode())
                    return False

                ms = line.strip().split("|")

                ########################################################
                # Population classifier
                if len(ms) >= 3 and ms[-1] == "pop_classifier":
                    try:
                        # safe_print("Processing pop_classifier request", file=sys.stdout)

                        eigenvec_file = ms[0]
                        pop_file = ms[1]

                        predictions = run_pop_classifier(eigenvec_file, pop_file)
                        predictions = "|".join(str(x) for x in predictions)

                        status = "success"
                        error_msg = ""

                    except Exception as e:
                        error_msg = f"POP_CLASSIFIER_ERROR: {str(e)}\n{traceback.format_exc()}\n"
                        status = "error"
                        predictions = ""
                        safe_print(error_msg, file=sys.stderr)
                    response = json.dumps(
                        {"status": status, "result": predictions, "message": error_msg}
                    )
                    conn.send((response + "\n").encode())

                ########################################################
                # PADRE
                elif len(ms) >= 1 and ms[-1] == "padre":
                    # make output directory
                    ersa_dir = f"{output_directory}/ersa"
                    if not os.path.exists(ersa_dir):
                        os.makedirs(ersa_dir, exist_ok=True)
                    ersa_outfile = f"{ersa_dir}/output_all_ersa"

                    try:
                        safe_print("PROGINFO: running ERSA for PADRE", file=sys.stderr)

                        # object to pass to ersa is the WHOLE dictionary
                        segment_obj = json.dumps(segment_dict)

                        # make output directory

                        ersa_options = {
                            "segment_dict": segment_obj,
                            "segment_files": segment_data_file,
                            "model_output_file": f"{ersa_outfile}.model",
                            "output_file": f"{ersa_outfile}.out",
                            "return_output": False,
                            "write_output": True,
                            "use_ibd2_siblings": ibd2_status,
                        }

                        # merge with additional options
                        ersa_options.update(additional_options)

                        ersa.runner(
                            ersa_options
                        )  # output written to file, not returned
                        error_msg = ""
                        status = "success"

                    except Exception as e:
                        error_msg = f"PADRE_ERROR: {str(e)}\n{traceback.format_exc()}\n"
                        ersa_outfile = ""
                        status = "error"
                        safe_print(error_msg, file=sys.stderr)

                    response = json.dumps(
                        {"status": status, "result": ersa_outfile, "message": error_msg}
                    )
                    conn.send((response + "\n").encode())

                ########################################################
                # Pairwise ERSA processing
                elif len(ms) >= 4:
                    id1, id2, vector_str, analysis_type = ms
                    try:
                        # Process the pairwise request
                        ersa_result = process_pairwise_ersa(
                            id1,
                            id2,
                            vector_str,
                            analysis_type,
                            segment_dict,
                            additional_options,
                            ibd2_status,
                            segment_data_file,
                            output_directory,
                        )

                        error_msg = ""
                        status = "success"
                    except ValueError as e:
                        error_msg = f"ERSA_ERROR: Malformed likelihood vector str. Expected 6 values. Provided str: {vector_str}\n{str(e)}\n{traceback.format_exc()}\n"
                        ersa_result = ""
                        status = "error"
                    except Exception as e:
                        error_msg = f"ERSA_ERROR: Failed to process {id1} <-> {id2}: {str(e)}\n{traceback.format_exc()}\n"
                        safe_print(error_msg, file=sys.stderr)
                        ersa_result = ""
                        status = "error"

                    result = json.dumps(
                        {"status": status, "result": ersa_result, "message": error_msg}
                    )
                    conn.send((result + "\n").encode())
                else:
                    error_msg = f"ERROR: During ERSA pairwise relatedness estimation received message with {len(ms)} parts: {ms}\n"
                    safe_print(error_msg, file=sys.stderr)
                    response = json.dumps(
                        {"status": "error", "result": "", "message": error_msg}
                    )
                    conn.send((response + "\n").encode())
    except KeyboardInterrupt:
        return True

    except Exception as e:
        error_msg = f"SOCKET_ERROR: {str(e)}\n{traceback.format_exc()}\n"
        safe_print(error_msg, file=sys.stderr)
        try:
            response = json.dumps(
                {"status": "error", "result": "", "message": error_msg}
            )
            conn.send((response + "\n").encode())
        except:
            pass  # Socket might be closed
    finally:
        conn.close()


def process_pairwise_ersa(
    id1,
    id2,
    vector_str,
    analysis_type,
    segment_dict,
    additional_options,
    ibd2_status,
    segment_data_file,
    output_directory,
) -> str:

    try:
        # Step 1: Parse inputs
        id1_temp = id1.split("_")[-1]  # just in case
        id2_temp = id2.split("_")[-1]

        idcombo = ":".join(tuple(sorted([id1_temp, id2_temp])))

        vector_arr = [float(x) for x in vector_str.split(",")]
        if len(vector_arr) != 6:
            raise ValueError(
                f"WARN: Expected 6 values in vector, got {len(vector_arr)}"
            )

        # Step 2: Check for segments
        if idcombo not in segment_dict.keys():
            # safe_print(f"No segments found for {idcombo} or {idcombo2}", file=sys.stderr)
            return "0,0,0,0,0,1"

        # Step 3: Prepare segment data
        try:
            segment_obj = {idcombo: segment_dict[idcombo]}

            segment_obj = json.dumps(segment_obj)

        except Exception as e:
            raise Exception(
                f"WARN: Failed to prepare segment data for {idcombo}: {str(e)}"
            )

        # Step 4: Set up ERSA options
        try:
            ersa_dir = Path(f"{output_directory}/ersa")
            if not ersa_dir.exists():
                ersa_dir.mkdir(exist_ok=True)
            ersa_outfile = f"{ersa_dir}/output_all_ersa"

            ersa_options = {
                "single_pair": idcombo,
                "segment_dict": segment_obj,
                "segment_files": segment_data_file,
                "model_output_file": f"{ersa_outfile}.model",
                "output_file": f"{ersa_outfile}.out",
                "return_output": True,
                "write_output": False,
                "use_ibd2_siblings": ibd2_status,
            }

            # merge with additional options
            ersa_options.update(additional_options)

        except Exception as e:
            raise Exception(f"FATAL: Failed to set up ERSA options: {str(e)}")

        # Step 5: Run ERSA
        try:
            # safe_print(f"Running ERSA for {key}", file=sys.stderr)
            output_model_df = ersa.runner(ersa_options)

        # except Exception as e:
        #     safe_print(f"ERSA failed for {key} with segments: {segment_dict[key]}", file=sys.stderr)
        #     raise Exception(f"ERSA runner failed: {str(e)}")

        except Exception as ersa_error:
            # Log the ERSA error but don't crash - return original vector
            safe_print(
                f"WARN: ERSA failed for {idcombo} with segments: {segment_dict[idcombo]}, reverting to original vector: {str(ersa_error)}",
                file=sys.stderr,
            )
            return vector_str

        # Step 6: Process results
        try:
            if output_model_df.empty:
                safe_print(
                    f"WARN: ERSA returned empty dataframe for {idcombo}",
                    file=sys.stderr,
                )
                return "0,0,0,0,0,1"

            if output_model_df["maxlnl"].isna().all():
                safe_print(
                    f"WARN: ERSA returned all NaN maxlnl values for {idcombo}",
                    file=sys.stderr,
                )
                return "0,0,0,0,0,1"

            ersa_props = calculate_ersa_props(output_model_df)
            prop02 = 1 - (vector_arr[0] + vector_arr[1])
            ersa_props_updated = tuple(x * prop02 for x in ersa_props)
            updated_vector = f"{vector_arr[0]},{vector_arr[1]},{ersa_props_updated[0]},{ersa_props_updated[1]},{ersa_props_updated[2]},{ersa_props_updated[3]}"

            # safe_print(f"ERSA completed successfully for {key}", file=sys.stderr)
            return updated_vector

        except Exception as e:
            raise Exception(f"FATAL: Failed to process ERSA results: {str(e)}")

    except Exception as e:
        # Re-raise with more context
        raise Exception(f"WARN: ERSA processing failed for {id1} <-> {id2}: {str(e)}")


if __name__ == "__main__":

    # This is how the script is called from bin/compadre_kickoff.pl

    segment_data_file = sys.argv[1]
    portnumber = int(sys.argv[2])
    ersa_flag_str = sys.argv[3]
    output_directory = sys.argv[4]

    additional_options = parse_ersa_options(ersa_flag_str)

    segment_dict, ibd2_status = load_segment_information(
        segment_data_file, additional_options["min_cm"]
    )

    safe_print(
        f"INFO: Number of pairs read in that have pairwise IBD data: {len(segment_dict)}",
        file=sys.stdout,
    )

    if segment_dict:  # Only if we have segments loaded
        additional_options = update_min_cm(
            additional_options, segment_dict
        )  # This will update min_cm passed into ersa if applicable

    main(segment_dict, additional_options, ibd2_status, portnumber, output_directory)
