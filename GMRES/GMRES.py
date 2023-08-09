import numpy as np
import numpy.random as random
from cla_utils import *


def arnoldi(A, b, k):
    """
    For a matrix A, apply k iterations of the Arnoldi algorithm,
    using b as the first basis vector.

    :param A: an mxm numpy array
    :param b: m dimensional numpy array, the starting vector
    :param k: integer, the number of iterations

    :return Q: an mx(k+1) dimensional numpy array containing the orthonormal basis
    :return H: a (k+1)xk dimensional numpy array containing the upper \
    Hessenberg matrix
    """
    m = A.shape[0]
    Q = np.zeros((m,k+1), dtype=A.dtype)
    H = np.zeros((k+1,k), dtype=A.dtype)
    q = b / np.sqrt(np.dot(b.conj(),b))
    Q[:,0] = q
    for i in range(k):
        v = A@q
        for j in range(i+1):  
            H[j,i] = np.dot(Q[:,j].conj(),v)
            v = v - H[j,i]*Q[:,j]
        H[i+1,i] = np.sqrt(np.dot(v.conj(),v))
        q = v / H[i+1,i]
        Q[:,i+1] = q
    return Q,H


# def GMRES(A, b, maxit, tol, x0=None, return_residual_norms=False,
#           return_residuals=False):
#     """
#     ss
#     For a matrix A, solve Ax=b using the basic GMRES algorithm.

#     :param A: an mxm numpy array
#     :param b: m dimensional numpy array
#     :param maxit: integer, the maximum number of iterations
#     :param tol: floating point number, the tolerance for termination
#     :param x0: the initial guess (if not present, use b)
#     :param return_residual_norms: logical
#     :param return_residuals: logical

#     :return x: an m dimensional numpy array, the solution
#     :return nits: if converged, the number of iterations required, otherwise \
#     equal to -1
#     :return rnorms: nits dimensional numpy array containing the norms of \
#     the residuals at each iteration
#     :return r: mxnits dimensional numpy array, column k contains residual \
#     at iteration k
#     """
#     m = A.shape[0]
#     q = b/np.linalg.norm(b)
#     if x0 is None:
#         x0 = b
#     r = np.linalg.norm(b-A@x0)
#     count = 1
#     for i in range(maxit):
#         Q,H = arnoldi(A,b,count)
#         e1 = np.zeros(count+1)
#         e1[0] = np.linalg.norm(b)
#         y = householder_ls(H,e1)
#         x_val = Q[:,:count]@y
#         r = np.linalg.norm(b-A@x_val)
#         count = count + 1
#         if r < tol:
#             return x_val,-1
#     return x_val,-1

# def GMRES(A, b, maxit, tol, x0=None, return_residual_norms=False,
#           return_residuals=False):
#     """
#     For a matrix A, solve Ax=b using the basic GMRES algorithm.

#     :param A: an mxm numpy array
#     :param b: m dimensional numpy array
#     :param maxit: integer, the maximum number of iterations
#     :param tol: floating point number, the tolerance for termination
#     :param x0: the initial guess (if not present, use b)
#     :param return_residual_norms: logical
#     :param return_residuals: logical

#     :return x: an m dimensional numpy array, the solution
#     :return nits: if converged, the number of iterations required, otherwise \
#     equal to -1
#     :return rnorms: nits dimensional numpy array containing the norms of \
#     the residuals at each iteration
#     :return r: mxnits dimensional numpy array, column k contains residual \
#     at iteration k
#     """
#     m = A.shape[0]
#     q = b/np.linalg.norm(b)
#     if (x0 is None):
#         x0 = b
#     r = b - A@x0    #get residual
#     count = 0
#     Q = np.zeros((m,maxit+1), dtype=A.dtype)
#     H = np.zeros((maxit+1,maxit), dtype=A.dtype)
#     q = b / np.sqrt(np.dot(b.conj(),b))
#     Q[:,0] = q
#     # res = tol*10 #set to number bigger than tolerance
#     for i in range (maxit):
#         #do arnoldi
#         v = A@q
#         for j in range(i+1):  
#             H[j,i] = np.dot(Q[:,j].conj(),v)
#             v = v - H[j,i]*Q[:,j]
#         H[i+1,i] = np.sqrt(np.dot(v.conj(),v))
#         q = v / H[i+1,i]
#         Q[:,i+1] = q
#         # b = np.zeros(maxit+1)
#         # b[0] = np.linalg.norm(r)
#         e1 = np.zeros(i+2)
#         e1[0] = 1
#         e1 = e1 * np.linalg.norm(b)
#         print("this is b",b[0])
#         y = householder_ls(H[:i+2,:i+1],e1)

#         print(y.shape[0])
#         reshpe = y.shape[0]
#         print("this is reshpe",reshpe)
#         print(H.shape)
#         y = y.reshape(reshpe,1)
#         print("y shape",y.shape)
#         print("h shape", H.shape)
#         r = np.linalg.norm(H[:i+2,:i+1]@y - b)#update residual
#         print("this is the shape of Q",Q.shape)
#         if r < tol:
#             return Q[:,:i+1]@y, 1
#         return Q[:,:i+1]@y,1

# def GMRES(A, b, maxit, tol, x0=None, return_residual_norms=True,
#           return_residuals=True):
#     """
#     For a matrix A, solve Ax=b using the basic GMRES algorithm.

#     :param A: an mxm numpy array
#     :param b: m dimensional numpy array
#     :param maxit: integer, the maximum number of iterations
#     :param tol: floating point number, the tolerance for termination
#     :param x0: the initial guess (if not present, use b)
#     :param return_residual_norms: logical
#     :param return_residuals: logical
#     :return x: an m dimensional numpy array, the solution
#     :return nits: if converged, the number of iterations required, otherwise \
#     equal to -1
#     :return rnorms: nits dimensional numpy array containing the norms of \
#     the residuals at each iteration
#     :return r: mxnits dimensional numpy array, column k contains residual \
#     at iteration k
#     """
#     m = A.shape[0]
#     if x0 is None: x0 = b
#     residual = b - A@x0
#     Q = np.zeros((m,maxit+1), dtype=A.dtype)
#     H = np.zeros((maxit+1,maxit), dtype=A.dtype)
#     q = b / np.sqrt(np.dot(b.conj(),b))
#     Q[:,0] = q
#     nits = 0
#     count = 1
#     if return_residual_norms is True:
#         residual_norms = np.zeros((maxit,1))
#         residual_norms[0] = np.linalg.norm(residual)
#     if return_residuals is True:
#         residuals = np.zeros((m,maxit))
#         residuals[:,0] = residual
#     for i in range(maxit):
#         v = A@Q[:,i] 
#         for j in range(i+1):
#             H[j,i] = np.dot(Q[:,j].conj(),v)
#             v = v - H[j,i]*Q[:,j]
#         H[i+1,i] = np.sqrt(np.dot(v.conj(),v))
#         Q[:,i+1] = v / H[i+1,i]
#         e1 = np.zeros(i+2)
#         e1[0] = 1
#         e1 = e1 * np.linalg.norm(b)
#         y = householder_ls(H[:i+2,:i+1],e1)
#         x = Q[:,:i+1]@y
#         residual = b-A@x
#         if return_residual_norms is True:
#             residual_norms[i] = np.linalg.norm(residual)
#         if return_residuals is True:
#             residuals[:,i] = residual
#         if np.linalg.norm(residual) < tol:
#             nits = i + 1
#             if return_residual_norms is True:
#                 residual_norms = residual_norms[~np.all(residual_norms==0, axis=1)]
#             if return_residuals is True:
#                 residuals = residuals[:,~np.all(residuals==0,axis=0)]
#             break
#         else:
#             nits = -1
#             continue
#     if return_residual_norms and return_residuals is True:
#         return x,nits, residual_norms,residuals
#     if return_residual_norms is True:
#         return x,nits,residual_norms
#     else:
#         return x,-1



def GMRES(A, b, maxit, tol, x0=None, return_residual_norms=False,
          return_residuals=False):
    """
    For a matrix A, solve Ax=b using the basic GMRES algorithm.

    :param A: an mxm numpy array
    :param b: m dimensional numpy array
    :param maxit: integer, the maximum number of iterations
    :param tol: floating point number, the tolerance for termination
    :param x0: the initial guess (if not present, use b)
    :param return_residual_norms: logical
    :param return_residuals: logical
    :return x: an m dimensional numpy array, the solution
    :return nits: if converged, the number of iterations required, otherwise \
    equal to -1
    :return rnorms: nits dimensional numpy array containing the norms of \
    the residuals at each iteration
    :return r: mxnits dimensional numpy array, column k contains residual \
    at iteration k
    """
    m = A.shape[0]
    if x0 is None: x0 = b
    residual = b - A@x0
    Q = np.zeros((m,maxit+1), dtype=A.dtype)
    H = np.zeros((maxit+1,maxit), dtype=A.dtype)
    q = b / np.sqrt(np.dot(b.conj(),b))
    Q[:,0] = q
    nits = 0
    count = 1
    if return_residual_norms is True:
        residual_norms = np.zeros((maxit,1))
        residual_norms[0] = np.linalg.norm(residual)
    if return_residuals is True:
        residuals = np.zeros((m,maxit))
        residuals[:,0] = residual
    for i in range(maxit):
        v = A@Q[:,i] 
        for j in range(i+1):
            H[j,i] = np.dot(Q[:,j].conj(),v)
            v = v - H[j,i]*Q[:,j]
        H[i+1,i] = np.sqrt(np.dot(v.conj(),v))
        Q[:,i+1] = v / H[i+1,i]
        e1 = np.zeros(i+2)
        e1[0] = 1
        e1 = e1 * np.linalg.norm(b)
        y = householder_ls(H[:i+2,:i+1],e1)
        x = Q[:,:i+1]@y
        residual = b-A@x
        if return_residual_norms is True:
            residual_norms[i] = np.linalg.norm(residual)
        if return_residuals is True:
            residuals[:,i] = residual
        if np.linalg.norm(residual) < tol:
            nits = i + 1
            if return_residual_norms is True:
                residual_norms = residual_norms[~np.all(residual_norms==0, axis=1)]
            if return_residuals is True:
                residuals = residuals[:,~np.all(residuals==0,axis=0)]
            break
        else:
            nits = -1
            continue


    if return_residual_norms and return_residuals is True:
        return x,nits, residual_norms,residuals
    if return_residual_norms is True:
        return x,nits,residual_norms
    else:
        return x,-1


def GMRES(A, b, maxit, tol, x0=None, return_residual_norms=False,
          return_residuals=False, callback=None):
    """
    For a matrix A, solve Ax=b using the basic GMRES algorithm.

    :param A: an mxm numpy array
    :param b: m dimensional numpy array
    :param maxit: integer, the maximum number of iterations
    :param tol: floating point number, the tolerance for termination
    :param x0: the initial guess (if not present, use b)
    :param return_residual_norms: logical
    :param return_residuals: logical
    :param callback
    :return x: an m dimensional numpy array, the solution
    :return nits: if converged, the number of iterations required, otherwise \
    equal to -1
    :return rnorms: nits dimensional numpy array containing the norms of \
    the residuals at each iteration
    :return r: mxnits dimensional numpy array, column k contains residual \
    at iteration k
    """
    m = A.shape[0]
    if x0 is None: x0 = b
    residual = b - A@x0
    Q = np.zeros((m,maxit+1), dtype=A.dtype)
    H = np.zeros((maxit+1,maxit), dtype=A.dtype)
    q = b / np.sqrt(np.dot(b.conj(),b))
    Q[:,0] = q
    nits = 0
    count = 1
    residual_norms = np.zeros((maxit,1))
    residual_norms[0] = np.linalg.norm(residual)
    residuals = np.zeros((m,maxit))
    residuals[:,0] = residual
    for i in range(maxit):
        v = A@Q[:,i] 
        for j in range(i+1):
            H[j,i] = np.dot(Q[:,j].conj(),v)
            v = v - H[j,i]*Q[:,j]
        H[i+1,i] = np.sqrt(np.dot(v.conj(),v))
        Q[:,i+1] = v / H[i+1,i]
        e1 = np.zeros(i+2)
        e1[0] = 1
        e1 = e1 * np.linalg.norm(b)
        y = householder_ls(H[:i+2,:i+1],e1)
        x = Q[:,:i+1]@y
        residual = b-A@x
        res_norm = np.linalg.norm(b-A@x)
        nits = nits + 1
        if return_residuals == True:
            residuals[:,i] = residual
        if return_residual_norms == True:
            residual_norms[i] = res_norm
        if np.linalg.norm(residual) < tol:
            break
    if return_residual_norms and return_residuals is True:
        return x,nits, residual_norms,residuals
    if return_residual_norms is True:
        return x,nits,residual_norms
    else:
        return x,-1



    

def get_AA100():
    """
    Get the AA100 matrix.

    :return A: a 100x100 numpy array used in exercises 10.
    """
    AA100 = np.fromfile('AA100.dat', sep=' ')
    AA100 = AA100.reshape((100, 100))
    return AA100


def get_BB100():
    """
    Get the BB100 matrix.

    :return B: a 100x100 numpy array used in exercises 10.
    """
    BB100 = np.fromfile('BB100.dat', sep=' ')
    BB100 = BB100.reshape((100, 100))
    return BB100


def get_CC100():
    """
    Get the CC100 matrix.

    :return C: a 100x100 numpy array used in exercises 10.
    """
    CC100 = np.fromfile('CC100.dat', sep=' ')
    CC100 = CC100.reshape((100, 100))
    return CC100
