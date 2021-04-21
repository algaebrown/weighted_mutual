import numpy as np




def weighted_entropy(ma, weight):
    '''
    input: row, p
    - row: contain 010101 of domain absence-presence of one protein; numpy.array
    - p: probability calculated from weight.probability; numpy.array
    output: weighted entropy of the row
    '''
    presence = np.matmul(ma.toarray(), weight.T) # n_gene by 1
    absensce = np.matmul(1-ma.toarray(), weight.T)
    
    probs = np.stack([presence.flatten()[0], absensce.flatten()[0]]) # 2*n_genes

    log_probs = np.log2(probs)
    log_probs[log_probs == -np.inf] = 0 # replace those with 0

    entropy = -np.sum(np.array(probs)*np.array(log_probs), axis = 0)


    return entropy

def weighted_joint_entropy(row1, row2, weight):
    '''
    return joint weighted entropy of two proteins
    '''
    
    one_one = (row1 * row2 * weight).sum() # 11
    one_zero = (row1 * (1-row2) * weight).sum() #10
    zero_one = ((1-row1) * row2 * weight).sum() #01
    zero_zero = ((1-row1) * (1-row2) * weight).sum() #00

    probs = np.array([one_one, one_zero, zero_zero, zero_one])
    log_probs = np.log2(probs)
    log_probs[log_probs == -np.inf] = 0 # replace those with 0

    entropy = -np.sum(np.array(probs)*np.array(log_probs))

    return entropy

def mutual_info(H_x, H_y, H_xy):
    return(H_x + H_y - H_xy)
