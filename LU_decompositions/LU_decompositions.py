import numpy as np

def get_Lk(m, lvec):
    """Compute the lower triangular row operation mxm matrix L_k 
    which has ones on the diagonal, and below diagonal entries
    in column k given by lvec (k is inferred from the size of lvec).

    :param m: integer giving the dimensions of L.
    :param lvec: a m-k dimensional numpy array.
    :return Lk: an mxm dimensional numpy array.

    """
    lvec_len = lvec.shape[0]
    L_k = np.eye(m)
    L_k[-lvec_len+m:,-lvec_len+m-1] = -lvec
    return L_k



def LU_inplace(A):
    """Compute the LU factorisation of A, using the in-place scheme so
    that the strictly lower triangular components of the array contain
    the strictly lower triangular components of L, and the upper
    triangular components of the array contain the upper triangular
    components of U.

    :param A: an mxm-dimensional numpy array
    :return LU: A in form of LU factorisation with elements of L and U in one matrix

    """
    # no outer product
    m = A.shape[0]
    for i in range(0,m-1):
        for j in range(i+1,m):
            l_term = A[j,i]/A[i,i]
            A[j,i:m] = A[j,i:m] - l_term*A[i,i:m]
            A[j,i] = l_term
    
    return A

    # outer product implementation (one loop)
    m = A.shape[0]
    for i in range(0,m-1):
        l_vec = A[i+1:m,i]/A[i,i]
        A[i+1:m,i:m] = A[i+1:m,i:m] - np.outer(l_vec,A[i,i:m])
        A[i+1:m,i] = l_vec
    return A 




def solve_L(L, b):
    """
    Solve systems Lx_i=b_i for x_i with L lower triangular, i=1,2,...,k

    :param L: an mxm-dimensional numpy array, assumed lower triangular
    :param b: an mxk-dimensional numpy array, with ith column containing 
       b_i
    :return x: an mxk-dimensional numpy array, with ith column containing 
       the solution x_i

    """
    m,m = L.shape
    m,k = b.shape
    x = np.zeros((m,k))

    x[0,:] = b[0,:]/L[0,0]
    for i in range(1,m):
        x[i,:] = (b[i,:] - np.dot(L[i,:],x[:,:]))/L[i,i]
    return x


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




def inverse_LU(A):
    """
    Form the inverse of A via LU factorisation.

    :param A: an mxm-dimensional numpy array.

    :return Ainv: an mxm-dimensional numpy array.

    """
    m = A.shape[0]
    Identity = np.eye(m)

    #Decompose A in L and U
    LU_fused = LU_inplace(A)
    U = np.triu(LU_fused)
    L = np.tril(LU_fused)
    np.fill_diagonal(L,1)

    L_output = solve_L(L,Identity)
    A_inv = solve_U(U,L_output)

    return A_inv
                     
