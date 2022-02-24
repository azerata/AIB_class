from math import log2
import numpy as np
import pandas as pd
from timeit import default_timer as timer
import seaborn as sns
import sys

import global_linear 


def construct_test_data(sequences:list):

    out = [ [] for _ in range(10)]

    for i in range(10):
        for j in range(len(sequences)):
            tmp = len(sequences[j]) // 10
            out[i].append(sequences[j][:(i+1)*tmp])

    return out



def time_alg(test_data):

    ns, times = [], []
    

    for n in range(10):
        start = timer()
        optimal_alignment(test_data[n][0], test_data[n][1],gap_cost )
        end = timer()

        ns.append(len(test_data[n][0]))
        times.append(end - start)

    return pd.DataFrame({'n':ns, 'time':times})

cost = outer_cost()

sequences = get_sequences(sys.argv[2])

gap_cost = sys.argv[3]


test_data = construct_test_data(sequences)

print(time_alg(test_data))
