import numpy as np


def entropy(series_of_p):

    no_zero = series_of_p[series_of_p != 0] # remove zero so that calculate entropy will not have problems with log(0). It will not affect the result since we have 0*log(0) should = 0

    entropy = -(no_zero * np.log2(no_zero)).sum()

    return(entropy)


def weighted_entropy(row, p):
    '''
    input: row, p
    - row: contain 010101 of domain absence-presence of one protein; numpy.array
    - p: probability calculated from weight.probability; numpy.array
    output: weighted entropy of the row
    '''

    presence = (row * p).sum()
    absence = ((1-row)*p).sum()

    # if any of presence, absence = 0, then np.log2 will have problem
    # so the entropy function will remove 0 for us

    return(entropy(np.array([presence, absence])))

def weighted_joint_entropy(row1, row2, p):
    '''
    return joint weighted entropy of two proteins
    '''
    one_one = (row1 * row2 * p).sum() # 11
    one_zero = (row1 * (1-row2) * p).sum() #10
    zero_one = ((1-row1) * row2 * p).sum() #01
    zero_zero = ((1-row1) * (1-row2) * p).sum() #00

    return(entropy(np.array([one_one, one_zero, zero_zero, zero_one])))

def mutual_info(H_x, H_y, H_xy):
    return(H_x + H_y - H_xy)
