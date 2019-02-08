

import numpy as np

def weight(file):
    with open(file) as f:
        next(f) # skip header
        domain_sum = 0
        for line in f:
            str_list = line.replace('\n', '').split(',')
            protein = np.array(str_list[1:], dtype = int)
            domain_sum = domain_sum + protein

    # calculate weights
    whole_array = np.sum(domain_sum)
    weight = whole_array * (1/domain_sum)

    return(weight)

def probability(weight):
    '''
    weight [3.2, 2.6, 3.2] -> probability [0.35, 0.29, 0.35] (numpy.array)
    '''
    return(weight/weight.sum())
