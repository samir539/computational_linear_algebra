import numpy as np
import math


def operator_2_norm(A):
    """
    Given a real mxn matrix A, return the operator 2-norm.

    :param A: an mxn-dimensional numpy array

    :return o2norm: the norm
    """
    #biggest eigenvalue of ATA
    eigenvals =  np.linalg.eig(A.T@A)[0]
    operator_2norm = np.max(eigenvals)
    operator_2norm_sqr = np.sqrt(operator_2norm)
    # raise NotImplementedError

    return operator_2norm_sqr 


def cond(A):
    """
    Given a real mxn matrix A, return the condition number in the 2-norm.

    :return A: an mxn-dimensional numpy array

    :param ncond: the condition number
    """
    #find biggest and smallest eigenvalues of ATA and then divide them
    eigenvals =  np.linalg.eig(A.T@A)[0]
    biggest = np.max(eigenvals)
    smallest = np.min(eigenvals)
    ncond = biggest/smallest
    ncond = np.sqrt(ncond)
    # raise NotImplementedError

    return ncond
