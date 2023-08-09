from operator import truediv
import numpy as np


def householder(A, kmax=None, swap=None, reduced_tol=None):
    """
    Given a real mxn matrix A, find the reduction to upper triangular matrix R
    using Householder transformations. The reduction should be done "in-place",
    so that A is transformed to R.

    :param A: an mxn-dimensional numpy array
    :param kmax: an integer, the number of columns of A to reduce \
    to upper triangular. If not present, will default to n.
    :param swap:bool apply swapping procedure such that the rank of the matrix A is given by R_11 of the produced R matrix
    :param reduced_tol: float, if given the algorithm treats x values with norm below the threshold as if they were zero
                      

    :return R: an mxn-dimensional numpy array containing the upper \
    triangular matrix
    """
    
    A = A.astype('float64', copy=False)
    m, n = A.shape
    if kmax is None:
        kmax = n
    if swap is None:    #implment choice of swap and not swap
        swap = False
    if reduced_tol is None:
        reduced_tol = False

    sign = lambda a: 1.0 if a >= 0 else -1.0   
    # single loop standard householder
    for i in range(kmax):
        if swap == True:
            candidate_cols = np.linalg.norm(A[i:m,i:n], axis=0)
            if reduced_tol != False:
                #swap normally if bigger than thres
                candidate_cols[candidate_cols < reduced_tol] = 0
                test_bool = not np.any(candidate_cols)
                if test_bool == True:
                    #produce R'
                    for i in range(m):
                        for j in range(n):
                            if np.abs(A[i,j]) < 0.0000001:
                                A[i,j] = 0
                    R_prime = A[~np.all(A == 0, axis=1)]
                    print("this is R PRIME",R_prime)
                    return R_prime
                    break
            index_of_max = candidate_cols.argmax() + i
            incoming_col = np.copy(A[i:m,index_of_max])
            A[i:m,index_of_max] = A[i:m,i]
            A[i:m,i]  = incoming_col
        x = A[i:m,i]
        e1 = np.zeros(m-i)  
        e1[0] = 1.0
        v = sign(x[0])*np.linalg.norm(x)*e1 + x
        v = v/np.linalg.norm(v)
        A[i:m,i:n] = A[i:m,i:n] - (2*np.outer(v,np.conj(v)))@A[i:m,i:n]
    R=A
    return R



def solve_U(U, b):
    """
    Solve systems Ux_i=b_i for x_i with U upper triangular, i=1,2,...,k

    :param U: an mxm-dimensional numpy array, assumed upper triangular
    :param b: an mxk-dimensional numpy array, with ith column containing 
       b_i
    :return x: an mxk-dimensional numpy array, with ith column containing 
       the solution x_i

    """
    m,m = U.shape
    m,k = b.shape
    x = np.zeros((m,k))
    #single column case
    # x[m-1] = b[m-1]/U[m-1,m-1]
    # for i in range(m-1,-1,-1):
    #     x[i] = (b[i] - np.dot(U[i,i+1:m],x[i+1:m]))/U[i,i]
    # return x

    #k column case
    x[m-1,:] = b[m-1,:]/U[m-1,m-1]
    for i in range(m-1,-1,-1):
        x[i,:] = (b[i,:] - np.dot(U[i,i+1:m],x[i+1:m,:]))/U[i,i]
    return x
    #CHECK TYPES
 #   raise NotImplementedError


def householder_solve(A, b):
    """
    Given a real mxm matrix A, use the Householder transformation to solve
    Ax_i=b_i, i=1,2,...,k.

    :param A: an mxm-dimensional numpy array
    :param b: an mxk-dimensional numpy array whose columns are the \
    right-hand side vectors b_1,b_2,...,b_k.
    

    :return x: an mxk-dimensional numpy array whose columns are the \
    right-hand side vectors x_1,x_2,...,x_k.
    """
    m = b.shape[0]
    k = b.shape[1]

    x = np.zeros((m,k))
    A_ext = np.hstack([A,b])
    A_house = householder(A_ext,m) # look at minus 1
    x = solve_U(A_house[:,0:m], A_house[:,m:k+m]) 
   #raise NotImplementedError

    #own idea
    # m = b.shape[0]
    # k = b.shape[1]
    # x = np.zeros((m,k))
    # A_house = householder(A,)
    # Q_star_b = householder(b)
    # # print(A_house.shape, Q_star_b.shape)
    # x = solve_U(A_house,Q_star_b)

    return x


def householder_qr(A):
    """
    Given a real mxn matrix A, use the Householder transformation to find
    the full QR factorisation o f A.

    :param A: an mxn-dimensional numpy array

    :return Q: an mxm-dimensional numpy array
    :return R: an mxn-dimensional numpy array
    """
    m = A.shape[0]
    n = A.shape[1]
    append_I = np.eye(m)
    A_hat = np.hstack([A,append_I])
    A_house = householder(A_hat, n)
    R = A_house[:,0:n]
    Q_star = A_house[:,n:m+n]
    Q = Q_star.T

  # raise NotImplementedError
    return Q, R


def householder_ls(A, b):
    """
    Given a real mxn matrix A and an m dimensional vector b, find the
    least squares solution to Ax = b.

    :param A: an mxn-dimensional numpy array
    :param b: an m-dimensional numpy array

    :return x: an n-dimensional numpy array
    """
    # Make augmented arrayls
    m = A.shape[0]
    n = A.shape[1]
    b = b.reshape(m,1)
    A_ext = np.hstack([A,b])
    A_house = householder(A_ext, n)
    R = A_house[0:n,0:n]
    Q_star_b = A_house[0:n,-1]
    Q_star_b_fix = Q_star_b.reshape(n,1)
    x = solve_U(R, Q_star_b_fix).reshape(n,)

    return x
