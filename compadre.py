
import os, math, sys, json, signal, socket
import pandas as pd
import warnings
warnings.filterwarnings('ignore', category=FutureWarning, module='pandas')

# COMPADRE imports -- ersa function and the SVM population classifier

import ersa
from pop_classifier import run_new as run_pop_classifier

########################################

# def signal_handler(signum, frame):
#     sys.exit(0)

def signal_handler(signum, frame):
    try:
        if 'server_socket' in globals():
            server_socket.close()
    except:
        pass
    os._exit(0)

def safe_print(*args, **kwargs):
    try:
        print(*args, **kwargs)
        sys.stdout.flush()
    except BrokenPipeError:
        # Python flushes standard streams on exit, so redirect remaining output to devnull to avoid another BrokenPipeError at shutdown
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, sys.stdout.fileno())
        sys.exit(0)  # Python exits with error code 1 on EPIPE

def calculate_ersa_props(model_df):

    if model_df.empty:
        return 0, 0, 0, 1

    model_df['degree_of_relatedness'] = model_df['degree_of_relatedness'].astype(int)
    model_df['maxlnl2'] = model_df['maxlnl'].apply(lambda x: math.exp(x))
    
    model_df_2d = model_df[model_df['degree_of_relatedness'] == 2] 
    model_df_2d_sum = model_df_2d['maxlnl2'].sum()
    model_df_3d = model_df[model_df['degree_of_relatedness'] == 3] 
    model_df_3d_sum = model_df_3d['maxlnl2'].sum()

    ## Old version
    # model_df_4d = model_df[model_df['degree_of_relatedness'] == 4] 
    # model_df_4d_sum = model_df_4d['maxlnl2'].sum()
    # model_df_un = model_df[(model_df['degree_of_relatedness'] >= 5) & (model_df['degree_of_relatedness'] <= 40)] 
    # model_df_un_sum = model_df_un['maxlnl2'].sum()

    ## New version -- updated 7.2.24
    # model_df_4d = model_df[(model_df['degree_of_relatedness'] >= 4) & (model_df['degree_of_relatedness'] < 10)]
    # model_df_4d_sum = model_df_4d['maxlnl2'].sum()
    # model_df_un = model_df[(model_df['degree_of_relatedness'] >= 10) & (model_df['degree_of_relatedness'] <= 40)] 
    # model_df_un_sum = model_df_un['maxlnl2'].sum()

    # Newest version
    model_df_4d = model_df[(model_df['degree_of_relatedness'] >= 4) & (model_df['degree_of_relatedness'] < 8)]
    model_df_4d_sum = model_df_4d['maxlnl2'].sum()
    model_df_un = model_df[(model_df['degree_of_relatedness'] >= 8) & (model_df['degree_of_relatedness'] <= 40)] 
    model_df_un_sum = model_df_un['maxlnl2'].sum()

    # re-proportion them 
    total_prop = model_df_2d_sum + model_df_3d_sum + model_df_4d_sum + model_df_un_sum
    model_df_2d_prop = model_df_2d_sum/total_prop
    model_df_3d_prop = model_df_3d_sum/total_prop
    model_df_4d_prop = model_df_4d_sum/total_prop
    model_df_un_prop = model_df_un_sum/total_prop

    return model_df_2d_prop, model_df_3d_prop, model_df_4d_prop, model_df_un_prop

########################################

def main(segment_data_file, portnumber):

    #signal.signal(signal.SIGPIPE, signal_handler)

    signal.signal(signal.SIGINT, signal_handler)   # CTRL+C
    signal.signal(signal.SIGTERM, signal_handler)  # Termination
    signal.signal(signal.SIGPIPE, signal_handler)  # Broken pipe
    global server_socket

    segment_dict = {}
    ibd2_status = 'false'
    ibd2_counter = 0

    if segment_data_file != 'NA':

        with open (segment_data_file, 'r') as f: # Open the large file here and populate dictionary that stays in system memory

            # Check number of columns to determine if it has IBD1/2 data or not 
            header_line = f.readline().strip()
            num_columns = len(header_line.split())

            if num_columns == 7: # File with segment-by-segment IBD status in 7th column 
                ibd2_status = 'true'

                if 'chrom' in header_line: #skip header if it exists
                    next(f) # skip header line in this file version

                for line in f:
                    ls = line.strip().split('\t')
                    iid1, iid2, start, end, cmlen, chrom, ibd = ls[0], ls[1], int(ls[2]), int(ls[3]), round(float(ls[4]), 2), int(ls[5]), int(ls[6].strip())
                    key = f"{iid1}:{iid2}"
                    value = (chrom, start, end, cmlen, ibd)

                    # comment this out later
                    if ibd == 2:
                        ibd2_counter += 1

                    if cmlen >= 5.0: 
                        if key not in segment_dict:
                            segment_dict[key] = [value,]
                        else:
                            segment_dict[key] += [value,]

            elif num_columns == 6: # File without segment-by-segment IBD status 
                
                if 'chrom' in header_line: #skip header if it exists
                    next(f)

                for line in f:
                    ls = line.split('\t')
                    if len(ls) == 11: # germline1
                        iid1, iid2, start, end, cmlen, chrom = ls[0].split(' ')[0], ls[1].split(' ')[0], int(ls[3].split(' ')[0]), int(ls[3].split(' ')[1]), round(float(ls[6]), 2), int(ls[2])
                    else: # germline2
                        iid1, iid2, start, end, cmlen, chrom = ls[0], ls[1], int(ls[2]), int(ls[3]), round(float(ls[4]), 2), int(ls[5].strip())
                    key = f"{iid1}:{iid2}"
                    value = (chrom, start, end, cmlen, 'NA') # no ibd1/2 data in this input

                    if cmlen >= 5.0:
                        if key not in segment_dict:
                            segment_dict[key] = [value,]
                        else:
                            segment_dict[key] += [value,]

            else:
                safe_print('[COMPADRE] Unrecognized segment file format. Please refer to the README (https://github.com/belowlab/compadre) for formatting guidelines.')
                sys.exit(1)

    # else, no data provided, and the socket is only being set up for running the pop classifier

    ####################################################################################################
    # Everything above this is done ONCE -- when COMPADRE starts -- and kept in memory for easy access when new requests are made over the socket

    # Check for $COMPADRE_HOST ENV variable 
    socket_host = os.environ['COMPADRE_HOST'] if 'COMPADRE_HOST' in os.environ else 'localhost'

    server_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    server_socket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    server_socket.setsockopt(socket.SOL_SOCKET, socket.SO_KEEPALIVE, 1)
    server_socket.settimeout(None)

    try:
        server_socket.bind((socket_host, portnumber))
        server_socket.listen(1)
        safe_print(f"COMPADRE helper is ready.")
        sys.stdout.flush()

        # At this point, the script executes a while loop as a way to wait for incoming messages

        while True:

            try: # New try except for error handling socket issues 

                conn, address = server_socket.accept()
                #safe_print(f"Connection from: {address}")

                # Set timeout to None for client connections too
                conn.settimeout(None)

                msg = conn.recv(1024).decode()
                #safe_print(f"Received data: {msg}")
                
                if msg == 'close':
                    conn.send("Closing server".encode())
                    conn.close()
                    break

                ms = msg.strip().split('|')

        # except socket.error as e:
        #     safe_print(f"Socket error: {e}")
        #     continue
        # except Exception as e:
        #     safe_print(f"Error handling client: {e}")
        #     if 'conn' in locals():
        #         conn.close()
        #     continue

                ########################################################

                if ms[-1] == 'pop_classifier': # Run population classifier script and return success message

                    eigenvec_file = ms[0]
                    pop_file = ms[1]

                    predictions = run_pop_classifier(eigenvec_file, pop_file)
                    predictions = "|".join(str(x) for x in predictions)

                    conn.send(predictions.encode())
                    conn.close() 


                elif ms[-1] == 'padre': # Run PADRE
                    
                    # object to pass to ersa is the WHOLE dictionary
                    segment_obj = json.dumps(segment_dict)

                    # make output directory
                    ersa_dir = segment_data_file.split('/')[:-1]
                    ersa_dir = '/'.join(ersa_dir) + '/ersa'
                    if not os.path.exists(ersa_dir):
                        os.makedirs(ersa_dir, exist_ok=True)
                    ersa_outfile = f'{ersa_dir}/output_all_ersa'
                    
                    ersa_options = {
                        "segment_dict": segment_obj,
                        "segment_files": segment_data_file,
                        "model_output_file": f"{ersa_outfile}.model",
                        "output_file": f"{ersa_outfile}.out",
                        "return_output": False,
                        "write_output": True
                    }
                    ersa.runner(ersa_options)

                    conn.send(ersa_outfile.encode())
                    conn.close()


                else: # "Pairwise" normal behavior

                    id1, id2, vector_str, analysis_type = ms

                    id1_temp = id1.split('_')[-1] # just in case
                    id2_temp = id2.split('_')[-1]
                    
                    idcombo = f"{id1_temp}:{id2_temp}"
                    idcombo2 = f"{id2_temp}:{id1_temp}"

                    if idcombo not in segment_dict and idcombo2 not in segment_dict: # No segments with either id combo permutation
                        result = '0,0,0,0,0,1' 
                    
                    else:

                        if idcombo in segment_dict:
                            key = idcombo
                            segment_obj = {key : segment_dict[key]} 
                            segment_obj = json.dumps(segment_obj)

                        else:
                            key = idcombo2
                            segment_obj = {key : segment_dict[key]} 
                            segment_obj = json.dumps(segment_obj)

                        ersa_dir = segment_data_file.split('/')[:-1]
                        ersa_dir = '/'.join(ersa_dir) + '/ersa'
                        if not os.path.exists(ersa_dir):
                            os.makedirs(ersa_dir, exist_ok=True)
                        ersa_outfile = f'{ersa_dir}/output_new_{id1_temp}_{id2_temp}'

                        ersa_options = {
                            "single_pair": key,
                            "segment_dict": segment_obj,
                            "segment_files": segment_data_file,
                            "model_output_file": f"{ersa_outfile}.model",
                            "output_file": f"{ersa_outfile}.out",
                            "return_output": True,
                            "write_output": False
                        }
                        output_model_df = ersa.runner(ersa_options) # run ersa function mode 

                        if output_model_df.empty: 
                            result = '0,0,0,0,0,1' 

                        else:

                            if output_model_df['maxlnl'].isna().all():
                                result = '0,0,0,0,0,1'
                            
                            else:
                                ersa_props = calculate_ersa_props(output_model_df)
                                vector_arr = [float(x) for x in vector_str.split(',')]
                                prop02 = 1 - (vector_arr[0] + vector_arr[1])
                                ersa_props_updated = tuple(x * prop02 for x in ersa_props) 
                                updated_vector = f'{vector_arr[0]},{vector_arr[1]},{ersa_props_updated[0]},{ersa_props_updated[1]},{ersa_props_updated[2]},{ersa_props_updated[3]}'
                                result = updated_vector


                    # send whatever 'result' is at the end of this logic back to Perl
                    conn.send(result.encode())
                    conn.close()

            except KeyboardInterrupt:
                break
            except Exception as e:
                safe_print(f"Error: {e}")
                continue

    #large_file.close()

    except KeyboardInterrupt:
        pass
    finally:
        try:
            server_socket.close()
        except:
            pass
        safe_print("COMPADRE helper shutdown complete.")


if __name__ == '__main__':

    # This file isn't really meant to be ran in isolation, but you could instantiate the socket connection on your own by running it this way then send it messages from another script using the appropriate port

    # Pass segment data file and port number via positional args 
    segment_data_file = sys.argv[1]
    portnumber = int(sys.argv[2])

    main(segment_data_file, portnumber)