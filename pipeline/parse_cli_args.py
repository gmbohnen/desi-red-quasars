import getopt, sys
import multiprocessing
    
def parse_cli_args(args):
    options = "hn:s:w:"
    long_options = ["Help", "N_Rounds=", "Batch_Size=", "N_Worker="]

    # establish default values
    n_rounds = 1
    batch_size = 500
    n_worker = multiprocessing.cpu_count() - 4

    try:
        arguments, values = getopt.getopt(args, options, long_options)

        for currentArg, currentVal in arguments:

            if currentArg in ("-h", "--Help"):
                sys.exit("-n --N_Rounds\tSet how many batches each worker should process\n-s --Batch_Size\tSet the batch size\n-w --N_Worker\tSet how many workers should be used\n-h --Help\tShow help")

            elif currentArg in ("-n", "--N_Rounds"):
                n_rounds = int(currentVal)

            elif currentArg in ("-s", "--Batch_Size"):
                batch_size = int(currentVal)

            elif currentArg in ("-w", "--N_Worker"):
                n_worker = int(currentVal)

    except getopt.error as err:
        sys.exit(str(err))
    
    return n_rounds, batch_size, n_worker