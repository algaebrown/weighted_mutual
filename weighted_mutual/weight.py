

import numpy as np

def weight(ma):
    ''' Given whole binary matrix, return weight for each column '''
    # how many times domain occurs
    domain_sum = ma.sum(axis = 0)

    whole_array = np.sum(domain_sum)
    # the rarer has larger weight
    weight = whole_array * (1/domain_sum)
    
    # return normalized weight
    return weight/weight.sum()

