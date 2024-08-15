
import os, math, sys, json, signal, socket
import ersa

'''
This script is a new version of compadre_helper.py

- It's launched from script.pl (the equivalent of the bin/ PRIMUS .pl files )
- After launch, it opens a connection on port 6000 which Module.pm sends data through to the main() statement

'''

def signal_handler(signum, frame):
    sys.exit(0)

def safe_print(*args, **kwargs):
    try:
        print(*args, **kwargs)
        sys.stdout.flush()
    except BrokenPipeError:
        # Python flushes standard streams on exit; redirect remaining output to devnull to avoid another BrokenPipeError at shutdown
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, sys.stdout.fileno())
        sys.exit(0)  # Python exits with error code 1 on EPIPE

def calculate_ersa_props(model_df):

    model_df['degree_of_relatedness'] = model_df['degree_of_relatedness'].astype(int)
    model_df['maxlnl2'] = model_df['maxlnl'].apply(lambda x: math.exp(x))
    
    model_df_2d = model_df[model_df['degree_of_relatedness'] == 2] 
    model_df_2d_sum = model_df_2d['maxlnl2'].sum()
    model_df_3d = model_df[model_df['degree_of_relatedness'] == 3] 
    model_df_3d_sum = model_df_3d['maxlnl2'].sum()

    # Old version
    model_df_4d = model_df[model_df['degree_of_relatedness'] == 4] 
    model_df_4d_sum = model_df_4d['maxlnl2'].sum()
    model_df_un = model_df[(model_df['degree_of_relatedness'] >= 5) & (model_df['degree_of_relatedness'] <= 40)] 
    model_df_un_sum = model_df_un['maxlnl2'].sum()

    ## New version
    # model_df_4d = model_df[(model_df['degree_of_relatedness'] >= 4) & (model_df['degree_of_relatedness'] < 10)] # updated 7.2.24
    # model_df_4d_sum = model_df_4d['maxlnl2'].sum()
    # model_df_un = model_df[(model_df['degree_of_relatedness'] >= 10) & (model_df['degree_of_relatedness'] <= 40)] 
    # model_df_un_sum = model_df_un['maxlnl2'].sum()

    # re-proportion them 
    total_prop = model_df_2d_sum + model_df_3d_sum + model_df_4d_sum + model_df_un_sum
    model_df_2d_prop = model_df_2d_sum/total_prop
    model_df_3d_prop = model_df_3d_sum/total_prop
    model_df_4d_prop = model_df_4d_sum/total_prop
    model_df_un_prop = model_df_un_sum/total_prop

    return model_df_2d_prop, model_df_3d_prop, model_df_4d_prop, model_df_un_prop


def main():

    signal.signal(signal.SIGPIPE, signal_handler)

    matchfile = sys.argv[1]
    portnumber = int(sys.argv[2])
   
    # convert to dict and keep in sys memory
    segment_dict = {}
    with open (matchfile, 'r') as f: # Open the large file here -- ONCE
        for line in f:
            ls = line.split('\t')
            if len(ls) == 11: # germline1
                iid1, iid2, start, end, cmlen, chrom = ls[0].split(' ')[0], ls[1].split(' ')[0], int(ls[3].split(' ')[0]), int(ls[3].split(' ')[1]), round(float(ls[6]), 2), int(ls[2])
            else: # germline2
                iid1, iid2, start, end, cmlen, chrom = ls[0], ls[1], int(ls[2]), int(ls[3]), round(float(ls[4]), 2), int(ls[5].strip())
            key = f"{iid1}:{iid2}"
            value = (chrom, start, end, cmlen)
            if cmlen >= 5.0:
                if key not in segment_dict:
                    segment_dict[key] = [value,]
                else:
                    segment_dict[key] += [value,]

    # if len(segment_dict) != 0:
    #     safe_print(f'Successfully populated segment dict with {len(segment_dict)} keys')

    
    ####################################################################################################
    # Everything above this is done ONCE -- when COMPADRE starts -- and kept in memory for easy access when new requests are made over the socket

    server_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    server_socket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    server_socket.bind(('localhost', portnumber))
    server_socket.listen(1)
    safe_print("COMPADRE helper socket is ready")
    sys.stdout.flush()

    # At this point, the script executes a while loop as a way to wait for incoming messages.

    while True:
        conn, address = server_socket.accept()
        #safe_print(f"Connection from: {address}")
        msg = conn.recv(1024).decode()
        #safe_print(f"Received data: {msg}")
        
        if msg == 'close':
            conn.send("Closing server".encode())
            conn.close()
            break

        # Everything below here is the 'default' handling for COMPADRE requests

        msgsplit = msg.split('|')
        id1 = msgsplit[0]
        id2 = msgsplit[1]
        vector_str = msgsplit[2].strip()
        #ersa_input_file = msgsplit[3].strip() ## No longer needed since main() launches with the matchfile provided from Perl via sysarg

        id1_temp = id1.split('_')[-1] # just in case
        id2_temp = id2.split('_')[-1]

        # limit the dictionary to just the parts we want 

        idcombo = f"{id1_temp}:{id2_temp}"

        if idcombo in segment_dict.keys():

            segment_obj = {idcombo : segment_dict[idcombo]} # we're using this as ersa input now, NOT the match file (too big)
            segment_obj = json.dumps(segment_obj)

            # Old code that handles ersa function mode options 

            ersa_dir = matchfile.split('/')[:-1]
            ersa_dir = '/'.join(ersa_dir) + '/ersa'
            if not os.path.exists(ersa_dir):
                os.makedirs(ersa_dir, exist_ok=True)
            ersa_outfile = f'{ersa_dir}/output_new_{id1_temp}_{id2_temp}'

            ersa_options = {
                "single_pair": f"{id1_temp}:{id2_temp}",
                "segment_dict": segment_obj,
                "segment_files": matchfile,
                "model_output_file": f"{ersa_outfile}.model",
                "output_file": f"{ersa_outfile}.out",
                "return_output": True,
                "write_output": False
            }
            output_model_df = ersa.runner(ersa_options) # run ersa function mode 

            if len(output_model_df) == 0: # no ersa data outputted
                result = vector_str 

            else:
                ersa_props = calculate_ersa_props(output_model_df)
                vector_arr = [float(x) for x in vector_str.split(',')]
                prop02 = 1 - (vector_arr[0] + vector_arr[1])
                ersa_props_updated = tuple(x * prop02 for x in ersa_props) 
                updated_vector = f'{vector_arr[0]},{vector_arr[1]},{ersa_props_updated[0]},{ersa_props_updated[1]},{ersa_props_updated[2]},{ersa_props_updated[3]}'
                #print (updated_vector) 
                result = updated_vector

        else: # Combo isn't in dictionary because they share zero segments >= 5cM 
            result = vector_str

        # send whatever 'result' is at the end of this logic back to Perl
        conn.send(result.encode())
        conn.close()

    #large_file.close()

if __name__ == '__main__':

    main()