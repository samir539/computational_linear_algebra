import numpy as np
import numpy.random as random
from cla_utils import *

def get_A100():
    """
    Return A100 matrix for investigating QR factoration.

    :return A: The 100x100 numpy array
    """
    m = 100
    random.seed(1111*m)
    A = random.randn(m, m) + 1j*random.randn(m, m)
    return A


def get_B100():
    """
    Return B100 matrix for investigating QR factoration.

    :return A: The 100x100 numpy array
    """
    m = 100
    random.seed(1111*m)
    A = random.randn(m, m) + 1j*random.randn(m, m)
    A[np.tril_indices(m, -2)] = 0
    return A


def get_C100():
    """
    Return C100 matrix for investigating QR factoration.

    :return A: The 100x100 numpy array
    """
    m = 100
    random.seed(1111*m)
    A = random.randn(m, m) + 1j*random.randn(m, m)
    A = 0.5*(A + np.conj(A).T)
    return A


def get_D100():
    """
    Return D100 matrix for investigating QR factoration.

    :return A: The 100x100 numpy array
    """
    m = 100
    random.seed(1111*m)
    A = random.randn(m, m) + 1j*random.randn(m, m)
    A = 0.5*(A + np.conj(A).T)
    A[np.tril_indices(m, -2)] = 0
    A[np.triu_indices(m, 2)] = 0
    return A


def get_A3():
    """
    Return A3 matrix for investigating power iteration.
    
    :return A3: a 3x3 numpy array.
    """

    return np.array([[ 0.76505141, -0.03865876,  0.42107996],
                     [-0.03865876,  0.20264378, -0.02824925],
                     [ 0.42107996, -0.02824925,  0.23330481]])


def get_B3():
    """
    Return B3 matrix for investigating power iteration.

    :return B3: a 3x3 numpy array.
    """
    return np.array([[ 0.76861909,  0.01464606,  0.42118629],
                     [ 0.01464606,  0.99907192, -0.02666057],
                     [ 0.42118629, -0.02666057,  0.23330798]])


def pow_it(A, x0, tol, maxit, store_iterations = False):
    """
    For a matrix A, apply the power iteration algorithm with initial
    guess x0, until either 

    ||r|| < tol where

    r = Ax - lambda*x,

    or the number of iterations exceeds maxit.

    :param A: an mxm numpy array
    :param x0: the starting vector for the power iteration
    :param tol: a positive float, the tolerance
    :param maxit: integer, max number of iterations
    :param store_iterations: if True, then return the entire sequence \
    of power iterates, instead of just the final iteration. Default is \
    False.

    :return x: an m dimensional numpy array containing the final iterate, or \
    if store_iterations, an mxmaxit dimensional numpy array containing all \
    the iterates.
    :return lambda0: the final eigenvalue.
    """
    m = A.shape[0]
    x = np.zeros((m,1))
    r_val = 0
    v = x0.reshape((m,1))
    count = 0
    lambda0 = 0
    while(r_val > tol and count < maxit):
        v = A@v
        v = v/np.linalg.norm(v)
        lambda0 = v.T@A@v
        r = A@v - lambda0*v
        r_val = np.linalg.norm(r)
        x = np.hstack((x,v))  
        count += 1

    if store_iterations == False:
        return x[:,-1], lambda0
    if store_iterations == True:
        return x[:,1:], lambda0
    # raise(NotImplementedError)
    # return x, lambda0


def inverse_it(A, x0, mu, tol, maxit, store_iterations = False):
    """
    For a Hermitian matrix A, apply the inverse iteration algorithm
    with initial guess x0, using the same termination criteria as
    for pow_it.
    :param A: an mxm numpy array
    :param mu: a floating point number, the shift parameter
    :param x0: the starting vector for the power iteration
    :param tol: a positive float, the tolerance
    :param maxit: integer, max number of iterations
    :param store_iterations: if True, then return the entire sequence \
    of inverse iterates, instead of just the final iteration. Default is \
    False.

    :return x: an m dimensional numpy array containing the final iterate, or \
    if store_iterations, an mxmaxit dimensional numpy array containing \
    all the iterates.
    :return l: a floating point number containing the final eigenvalue \
    estimate, or if store_iterations, a maxit dimensional numpy array containing \
    all the iterates.
    """
    m = A.shape[0]
    x = np.zeros((m,1))
    r_val = 10
    v = x0.reshape((m,1))
    count = 0
    lambda0 = 0
    while(r_val > tol and count < maxit):
        v = solve_LUP(A-mu*np.identity(m),v)
        v = v.reshape((m,1))
        v = v/np.linalg.norm(v)
        lambda0 = v.T@A@v
        r = A@v - lambda0*v
        r_val = np.linalg.norm(r)
        x = np.hstack((x,v))  
        count += 1
        # print(count)

    if store_iterations == False:
        return x[:,-1], lambda0
    if store_iterations == True:
        return x[:,1:], lambda0


    # raise NotImplementedError
    


def rq_it(A, x0, tol, maxit, store_iterations = False):
    """
    For a Hermitian matrix A, apply the Rayleigh quotient algorithm
    with initial guess x0, using the same termination criteria as
    for pow_it.

    :param A: an mxm numpy array
    :param x0: the starting vector for the power iteration
    :param tol: a positive float, the tolerance
    :param maxit: integer, max number of iterations
    :param store_iterations: if True, then return the entire sequence \
    of inverse iterates, instead of just the final iteration. Default is \
    False.

    :return x: an m dimensional numpy array containing the final iterate, or \
    if store_iterations, an mxmaxit dimensional numpy array containing \
    all the iterates.
    :return l: a floating point number containing the final eigenvalue \
    estimate, or if store_iterations, an m dimensional numpy array containing \
    all the iterates.
    """
    
    m = A.shape[0]
    x = np.zeros((m,1))
    x0 = x0/np.linalg.norm(x0) #normalise
    l = x0.T@A@x0
    v = x0.reshape((m,1))
    r_val = 0
    count = 0
    while(r_val > tol and count < maxit):
        v = solve_LUP(A-l*np.identity(m),v)
        v = v.reshape((m,1))
        v = v/np.linalg.norm(v)
        l = v.T@A@v
        r = A@v - l*v
        r_val = np.linalg.norm(r)
        x = np.hstack((x,v))  
        count += 1
    if store_iterations == False:
        return x[:,-1], l
    if store_iterations == True:
        return x[:,1:], l


    raise NotImplementedError


def pure_QR(A, maxit, tol):
    """
    For matrix A, apply the QR algorithm and return the result.

    :param A: an mxm numpy array
    :param maxit: the maximum number of iterations
    :param tol: termination tolerance
    :return Ak: the result
    """
    m = A.shape[0]
    # count = 0
    # for i in range(maxit):
    #     count = count + 1
    #     Q,R = householder_qr(A)
    #     A = R@Q
    #     if (np.linalg.norm(A[np.tril_indices(m, -1)])/m**2  < tol):
    #         return A
    # return A
    count = 0
    while((np.linalg.norm(A[np.tril_indices(m, -1)])/m**2  > tol) and count < maxit): #should be and
        count = count + 1
        Q,R = householder_qr(A)
        A = R@Q
    return A

# def pure_QR(A, maxit, tol):
#     """
#     For matrix A, apply the QR algorithm and return the result.

#     :param A: an mxm numpy array
#     :param maxit: the maximum number of iterations
#     :param tol: termination tolerance
#     :return Ak: the result
#     """
#     m = A.shape[0]
#     count = 0
#     for i in range(maxit):
#         count = count + 1
#         Q,R = householder_qr(A)
#         A = R@Q
#         if (np.linalg.norm(A[np.tril_indices(m, -1)])/m**2  < tol):
#             return A
#     return A



        


