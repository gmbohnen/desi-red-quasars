from pipeline.pipeline import pipeline
from pipeline.parse_cli_args import parse_cli_args
from pipeline.get_args_list import get_args_list

import multiprocessing
import json
import sys
import time


def process_batch(args):
    process_idx = args[0]
    target_idxs = args[1]

    results_path = f"data/processing/batches/results_{process_idx}.jsonl"
    processed_path = f"data/processing/processed_targetid_files/processed_{process_idx}.txt"

    for idx in target_idxs:

        # run pipeline
        measurements = pipeline(idx)

        # dump results into results json
        with open(results_path, "a") as f:
            f.write(json.dumps(measurements) + "\n")
        
        # additionally store targetid in the processed list
        with open(process_batch, "a") as f:
            f.write(idx + "\n")
    
    return process_idx


if __name__ == "__main__":

    # parse arguments
    cli_args = sys.argv[1:]

    n_rounds, batch_size, n_worker = parse_cli_args(cli_args)
    n_batches = n_rounds * n_worker

    print(f"Start processing of {n_batches} batches of size {batch_size} with {n_worker} workers...")


    # for each round args_list is created again so the lists stored are not too large
    for _ in range(n_rounds):
        tic = time.time()

        args_list = get_args_list(batch_size=batch_size, n_worker=n_worker)

        # create pool, run parallel processing
        pool = multiprocessing.Pool(processes=n_worker)
        results = pool.map(process_batch, args_list)
        pool.close()
        pool.join()

        toc = time.time() - tic

        print(f"Round 1 finished in {toc/60:.2f} minutes.")

    print("--------------------\nProcessing finished.")