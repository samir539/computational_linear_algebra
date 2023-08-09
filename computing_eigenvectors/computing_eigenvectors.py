import numpy as np

def Q1AQ1s(A):
    """
    For a matrix A, find the unitary matrix Q1 such that the first
    column of Q1*A has zeros below the diagonal. Then return A1 = Q1*A*Q1^*.

    :param A: an mxm numpy array

    :return A1: an mxm numpy array
    """
    m = A.shape[0]
    # do Q1*A 
    sign = lambda a: 1.0 if a >= 0 else -1.0    #make use of numpy sign
    #single loop
    x = A[0:m,0]
    e1 = np.zeros(m)  #could be better
    e1[0] = 1.0
    v = sign(x[0])*np.linalg.norm(x)*e1 + x
    v = v/np.linalg.norm(v)
    A[0:m,0:m] = A[0:m,0:m] - (2*np.outer(v,np.conj(v)))@A[0:m,0:m]     #LOOK AT MATRIX MATRIX MULTIPLICATION
    A[0:m,0:m] = A[0:m,0:m] - (2*A[0:m,0:m]@np.outer(v,np.conj(v)))
    A1 = A
    return A1
    # raise NotImplementedError


def hessenberg(A):
    """
    For a matrix A, transform to Hessenberg form H by Householder
    similarity transformations, in place.

    :param A: an mxm numpy array
    """
    m = A.shape[0]

    for i in range(m-1):
        sign = lambda a: 1.0 if a >= 0 else -1.0    #make use of numpy sign
        x = A[i+1:m,i]
        e1 = np.zeros(m-(i+1))  #could be better
        e1[0] = 1.0
        v = sign(x[0])*np.linalg.norm(x)*e1 + x
        v = v/np.linalg.norm(v)
        A[i+1:m,i:m] = A[i+1:m,i:m] - (2*np.outer(v,np.conj(v)))@A[i+1:m,i:m]     #LOOK AT MATRIX MATRIX MULTIPLICATION
        A[i:m,i+1:m] = A[i:m,i+1:m] - (2*A[i:m,i+1:m]@np.outer(v,np.conj(v)))
        

    


def hessenbergQ(A):
    """
    For a matrix A, transform to Hessenberg form H by Householder
    similarity transformations, in place, and return the matrix Q
    for which QHQ^* = A.

    :param A: an mxm numpy array
    
    :return Q: an mxm numpy array
    """
    #Append identity and do householder on both A and the extention 
    m = A.shape[0]
    append_I = np.eye(m)
    A_stack = np.hstack([A,append_I])
    for i in range(m-1):
        sign = lambda a: 1.0 if a >= 0 else -1.0    #make use of numpy sign
        x = A_stack[i+1:m,i]
        e1 = np.zeros(m-(i+1))  #could be better
        e1[0] = 1.0
        v = sign(x[0])*np.linalg.norm(x)*e1 + x
        v = v/np.linalg.norm(v)
        A_stack[i+1:m,i:m] = A_stack[i+1:m,i:m] - (2*np.outer(v,np.conj(v)))@A_stack[i+1:m,i:m]     #LOOK AT MATRIX MATRIX MULTIPLICATION
        A_stack[i:m,i+1:m] = A_stack[i:m,i+1:m] - (2*A_stack[i:m,i+1:m]@np.outer(v,np.conj(v)))
    #extract A
    A = A_stack[:,0:m]
    print("THIS IS A FROM HESSENBERG")
    #extract Q
    Q = A_stack[:,m:m+m]
    return Q


    raise NotImplementedError

def hessenberg_ev(H):
    """
    Given a Hessenberg matrix, return the eigenvectors.

    :param H: an mxm numpy array

    :return V: an mxm numpy array whose columns are the eigenvectors of H

    Do not change this function.
    """
    m, n = H.shape
    assert(m==n)
    assert(np.linalg.norm(H[np.tril_indices(m, -2)]) < 1.0e-6)
    _, V = np.linalg.eig(H)
    return V


def ev(A):
    """
    Given a matrix A, return the eigenvectors of A. This should
    be done by using your functions to reduce to upper Hessenberg
    form, before calling hessenberg_ev (which you should not edit!).

    :param A: an mxm numpy array

    :return V: an mxm numpy array whose columns are the eigenvectors of A
    """
    hessenberg(A)
    V = hessenberg_ev(A)
    return V

    # raise NotImplementedError
