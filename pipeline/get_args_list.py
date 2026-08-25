import os


PROCESSED_PATH = "/home/leya/Code/Uni/desi-red-quasars/data/processing/processed_targetid_files"
TARGETID_PATH = "/home/leya/Code/Uni/desi-red-quasars/data/processing/targetid_list.txt"


def get_args_list(batch_size,n_worker):
    '''Create list of arguments for one round of processing, i.e. a list of tuples (process_id,list_of_targetids) for each worker.'''

    with open(TARGETID_PATH) as f:
        all_ids = [line.strip() for line in f if line.strip()]

    processed_files = os.listdir(PROCESSED_PATH)

    processed_ids = []

    for elem in processed_files:
        with open(f"{PROCESSED_PATH}/{elem}","r") as f:
            processed_ids.extend([line.strip() for line in f if line.strip()])

    unprocessed_ids = list(set(all_ids) - set(processed_ids)) #[elem for elem in all_ids if str(elem) not in processed_ids]

    args_list = []

    for w in range(n_worker):
        idxs = unprocessed_ids[w*batch_size:(w+1)*batch_size]
        args_list.append((w,idxs))
    
    return args_list