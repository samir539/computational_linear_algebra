import numpy as np


def orthog_cpts(v, Q):
    """
    Given a vector v and an orthonormal set of vectors q_1,...q_n,
    compute v = r + u_1q_1 + u_2q_2 + ... + u_nq_n
    for scalar coefficients u_1, u_2, ..., u_n and
    residual vector r

    :param v: an m-dimensional numpy array
    :param Q: an mxn-dimensional numpy array whose columns are the \
    orthonormal vectors

    :return r: an m-dimensional numpy array containing the residual
    :return u: an n-dimensional numpy array containing the coefficients
    """
    #Find r
    r = np.array([]) #m dimensional 
    u = np.array([]) #n dimensional 
    dir_sum = 0
    for i in range(Q.shape[1]):
        dir_sum += np.dot(np.conjugate(Q[:,i]), v)*Q[:,i]
        u = np.append(u, np.dot(np. conjugate(Q[:,i]),v))
    r = v - dir_sum

    #Find coef of component


  #  raise NotImplementedError

    return r, u


def solveQ(Q, b):
    """
    Given a unitary mxm matrix Q and a vector b, solve Qx=b for x.

    :param Q: an mxm dimensional numpy array containing the unitary matrix
    :param b: the m dimensional array for the RHS

    :return x: m dimensional array containing the solution.
    """
    #Note
    #Qx=b
    #Q*Qx = Q*b
    #Ix = Q*b
    #x = Q*b
    x = np.array([])
    Q_star = np.conjugate(Q).T
    x = Q_star@b




   # raise NotImplementedError

    return x


def orthog_proj(Q):
    """
    Given a vector v and an orthonormal set of vectors q_1,...q_n,
    compute the orthogonal projector P that projects vectors onto
    the subspace spanned by those vectors.

    :param Q: an mxn-dimensional numpy array whose columns are the \
    orthonormal vectors

    :return P: an mxm-dimensional numpy array containing the projector
    """
    P = Q@np.conjugate(Q).T
   # raise NotImplementedError

    return P


def orthog_space(V):
    """
    Given set of vectors u_1,u_2,..., u_n, compute the
    orthogonal complement to the subspace U spanned by the vectors.

    :param V: an mxn-dimensional numpy array whose columns are the \
    vectors u_1,u_2,...,u_n.

    :return Q: an mxl-dimensional numpy array whose columns are an \
    orthonormal basis for the subspace orthogonal to U, for appropriate l.
    """
    q,r = np.linalg.qr(V)
    Q = q[:,V.shape[1]:V.shape[0]]

  # raise NotImplementedError

    return Q


def GS_classical(A):
    """
    Given an mxn matrix A, compute the QR factorisation by classical
    Gram-Schmidt algorithm, transforming A to Q in place and returning R.

    :param A: mxn numpy array

    :return R: nxn numpy array
    """
    #First try with double loop [works but fails the tests]
    # rows,cols = np.shape(A)
    # A = A.astype("clongdouble")
    # Q = np.empty([rows,cols], dtype="clongdouble")
    # for j in range(cols):
    #     v = A[:,j]
    #     for i in range(j):
    #         r_ij = np.dot(Q[:,i],A[:,j])
    #         v = v - r_ij*Q[:,i]
    #     r_jj = np.linalg.norm(v)
    #     Q[:,j] = v/r_jj

    # R = np.conjugate(Q).T@A
  #  raise NotImplementedError
  #single loop
    # m, n = np.shape(A)
    # R = np.zeros((n,n), dtype = "clongdouble")
    # Q = A.copy().astype('clongdouble')
    # for i in range(n):
    #     if i == 0:
    #         R[i,i] = np.linalg.norm(A[:,i])
    #         Q[:,i] = Q[:,i]/R[i,i]
    #     else:
    #         R[0:i,i] = (np.conjugate(Q).T@A)[0:i,i]
            
    #         v = Q[:,i] - Q@R[:,i] 
    #         R[i,i] = np.linalg.norm(v)
    #         Q[:,i] = v/R[i,i]

    # return R,Q

    # m, n = np.shape(A)
    # R = np.zeros((n,n), dtype = "clongdouble")
    # Q = A
    # for i in range(n):  
    #     if i == 0:
    #         R[i,i] = np.linalg.norm(A[:,i])
    #         Q[:,i] = Q[:,i]/R[i,i]
    #     else:
    #         R[0:i,i] = (np.conjugate(Q).T@A)[0:i,i]
            
    #         Q[:,i] = Q[:,i] - Q@R[:,i] 
    #         R[i,i] = np.linalg.norm(Q[:,i])
    #         Q[:,i] = Q[:,i]/R[i,i]

#     # return R,Q
# #NO IF STATEMENT IN FOR LOOP
#     A = A.astype(A.dtype, copy=False)

#     m,n = A.shape
#     R = np.zeros([n,n], dtype= A.dtype)
#     for i in range(n):
#         v_i  = A[:,i]
#         print(v_i, "this is v_i for iteration", i)
#         for j in range(i):
#             R[i,j] = np.dot(np.conj(np.transpose(A[:,i])),A[:,j]) 
#             print(R[i,j],"this is rij for ", i,j)
#             v_i = v_i - R[i,j]*A[:,j]
#             print("this is v_i", v_i)
#         print("THIS IS NORM",np.linalg.norm(v_i))
#         R[i,i] = np.linalg.norm(v_i)

#         A[:,i] = v_i/R[i,i]
#         print(A[:,i],"this is col",i, "for iteration", i)
#     return R.T

    #Single loop
    m,n = A.shape
    R = np.zeros([n,n], dtype = A.dtype)
    R[0,0] = np.linalg.norm(A[:,0])
    A[:,0] = A[:,0]/R[0,0]
    print("This is R[0,0]", R[0,0])
    print("this is a1", A[:,0])
    for i in range(1,n):
        v_i = A[:,i]
        print("this is v_i", v_i, "for step", i)
        R[0:i,i] = (A.T@A[:,i])[0:i]
        print("this is the R coloumn for the", i, "step", R[0:i,i])
        v_i = v_i - A@R[:,i]
        # print("THIS IS THE MAT MUL",R[:,i]@A[:,0:i])
        R[i,i] = np.linalg.norm(v_i)
        A[:,i] = v_i/R[i,i]
        print("THIS IS R", R)
    return R


def GS_modified(A):
    """
    Given an mxn matrix A, compute the QR factorisation by modified
    Gram-Schmidt algorithm, transforming A to Q in place and returning
    R.

    :param A: mxn numpy array

    :return R: nxn numpy array
    """
    m,n = np.shape(A)
    R = np.zeros([n,n], dtype = A.dtype)
    # for i in range(n):
    #     v_i = A[:,i]
    #     R[i,i] = np.linalg.norm(v_i)
    #     A[:,i] = v_i/R[i,i]
    #     # print("this is the col",A[:,i])
    #     for j in range(i+1,n):
    #         R[i,j]  = np.dot(np.conjugate(A[:,i].T),A[:,j])
    #         A[:,j] = A[:,j] - R[i,j]*A[:,i]

    for i in range(n):
        v_i = A[:,i]
        R[i,i] = np.linalg.norm(v_i)
        A[:,i] = v_i/R[i,i]
        R[i,i+1:] = (np.conjugate(A[:,i]))@A[:,i+1:n]
        A[:,i+1:n] = A[:,i+1:n] - np.outer(A[:,i],R[i,i+1:])
    return R,A  
   # raise NotImplementedError

   #single loop
    
    #return R


def GS_modified_get_R(A, k):
    """
    Given an mxn matrix A, with columns of A[:, 0:k] assumed orthonormal,
    return upper triangular nxn matrix R such that
    Ahat = A*R has the properties that
    1) Ahat[:, 0:k] = A[:, 0:k],
    2) A[:, k] is normalised and orthogonal to the columns of A[:, 0:k].

    :param A: mxn numpy array
    :param k: integer indicating the column that R should orthogonalise

    :return R: nxn numpy array
    """

    raise NotImplementedError

    return R

def GS_modified_R(A):
    """
    Implement the modified Gram Schmidt algorithm using the lower triangular
    formulation with Rs provided from GS_modified_get_R.

    :param A: mxn numpy array

    :return Q: mxn numpy array
    :return R: nxn numpy array
    """

    m, n = A.shape
    A = 1.0*A
    R = np.eye(n, dtype=A.dtype)
    for i in range(n):
        Rk = GS_modified_get_R(A, i)
        A[:,:] = np.dot(A, Rk)
        R[:,:] = np.dot(R, Rk)
    R = np.linalg.inv(R)
    return A, R
