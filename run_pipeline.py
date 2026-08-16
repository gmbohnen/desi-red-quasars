from pipeline.pipeline import pipeline
from pipeline.parse_cli_args import parse_cli_args

import multiprocessing
import json
import sys

def process_batch(args):
    process_idx = args[0]
    target_idxs = args[1]

    # for idx in target_idxs:

        # run pipeline
        # TODO: pipeline() should only take idx as argument and load the arrays itself
        # TODO: should return dictionary of return values

        # dump results into results json

    pass


if __name__ == "__main__":

    ## PARSE COMMAND LINE ARGUMENTS

    cli_args = sys.argv[1:]

    n_rounds, batch_size, n_worker = parse_cli_args(cli_args)
    n_batches = n_rounds * n_worker

    print(f"Start processing of {n_batches} batches of size {batch_size} with {n_worker} workers...")



    # load idxs that go into each batch
    # args_list = []
    # for i in range(n_workers):
    #     args_list.append((i,test_idxs[30*i:30*(i+1)]))

    # set up workers

    # pool = multiprocessing.Pool(processes=n_workers)
    # results = pool.map(process_batch, args_list)
    # pool.close()
    # pool.join()
    # TODO do something with the results, or actually decide what the results should be lol

    pass