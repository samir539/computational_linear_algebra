import numpy as np

def perm(p, i, j):
    """
    For p representing a permutation P, i.e. Px[i] = x[p[i]],
    replace with p representing the permutation P_{i,j}P, where
    P_{i,j} exchanges rows i and j.

    :param p: an m-dimensional numpy array of integers.
    """
    i = i
    j = j
    permutation_mat = p
    permutation_mat[i], permutation_mat[j] = permutation_mat[j], permutation_mat[i]
    return permutation_mat
    # raise NotImplementedError


def LUP_inplace(A,count_swaps = False):
    """
    Compute the LUP factorisation of A with partial pivoting, using the
    in-place scheme so that the strictly lower triangular components
    of the array contain the strictly lower triangular components of
    L, and the upper triangular components of the array contain the
    upper triangular components of U.

    :param A: an mxm-dimensional numpy array

    :return p: an m-dimensional integer array describing the permutation \
    i.e. (Px)[i] = x[p[i]]
    """

    #outer product implementation
    m = A.shape[0]
    p = np.arange(m) 
    swap_count = 0
    for i in range(m-1):
        pivot_index = np.abs(A[i:m,i]).argmax()
        pivot_index = int(pivot_index)
        if np.abs(A[i:m,i]).max() != A[i:m,i:m][0,0]:
            swap_count += 1
        A[i:m,i:m][[0,pivot_index]] = A[i:m,i:m][[pivot_index,0]]
        perm(p,i,pivot_index+i)    
        l_vec = A[i+1:m,i]/A[i,i]
        A[i+1:m,i:m] = A[i+1:m,i:m] - np.outer(l_vec,A[i,i:m])
        A[i+1:m,i] = l_vec
        for k in range(i):
            A[i:m,k][[0,pivot_index]] = A[i:m,k][[pivot_index,0]]
    if count_swaps == False:
        return p
    else:
        return p,swap_count

    

#######extra functions#######
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

    x[m-1,:] = b[m-1,:]/U[m-1,m-1]
    for i in range(m-1,-1,-1):
        x[i,:] = (b[i,:] - np.dot(U[i,i+1:m],x[i+1:m,:]))/U[i,i]
    return x

#############################


def solve_LUP(A, b):
    """
    Solve Ax=b using LUP factorisation.

    :param A: an mxm-dimensional numpy array
    :param b: an m-dimensional numpy array

    :return x: an m-dimensional numpy array
    """
    #reshape b
    b_dim = b.shape[0]
    b_reshaped = b.reshape(b_dim,1)

    #get permutation list
    permutation_list = LUP_inplace(A)

    #generate permutation matrix
    dim = permutation_list.shape[0]
    perm_mat = np.zeros((dim,dim))
    index = 0
    for i in permutation_list:
        perm_mat[index,i] = 1
        index = index + 1
    
    #get Pb
    Pb = perm_mat@b_reshaped
 
    #get L and U
    LU_fused = A
    U = np.triu(LU_fused)
    L = np.tril(LU_fused)
    np.fill_diagonal(L,1)

    #forward sub to find y: Ly = Pb = c
    y = solve_L(L,Pb)

    #backward sub to find x: Ux = y
    x = solve_U(U,y)
    x_dim = x.shape[0]
    x = x.reshape(x_dim,)
    
    return x
                     

def det_LUP(A):
    """
    Find the determinant of A using LUP factorisation.

    :param A: an mxm-dimensional numpy array

    :return detA: floating point number, the determinant.
    """
    #A = LU 
    #det(A) = det(L)det(U)
    permutation_list, no_of_swaps = LUP_inplace(A,count_swaps=True)

    #get L and U
    LU_fused = A
    U = np.triu(LU_fused)
    L = np.tril(LU_fused)
    np.fill_diagonal(L,1)

    #det(L)
    det_L = L.diagonal().prod()
    print("this is det_L",det_L)

    #det(U)
    det_U = U.diagonal().prod()
    det_A = det_U * det_L

    #Correct sign
    print("this is the number of swaps", no_of_swaps)
    det_A = det_A * (-1**(no_of_swaps))
    print(det_A)
    #find the number of swaps which occured
    # m = A.shape[0]
    # ref_vec = np.arange(m)
    # print("this is the ref_vec",ref_vec)
    # print("this is the permutation list",permutation_list)
    # # swap_count = ref_vec - permutation_list
    # # print("this is the difference",swap_count)
    # # print("the number of nonzero is", np.count_nonzero(swap_count))
    # # print("this is m",m)
    # # no_of_swaps = m - (m - np.count_nonzero(swap_count)) - 1
    # # print("this is no of swaps", no_of_swaps)
    # # print(det_U)
    # # if no_of_swaps == 0:
    # #     det_A = det_U * 1
    # #     print("this branch happened")
    # # else:
    # #     det_A = det_U * ((-1)**(no_of_swaps))
    # #     print("this ELSE branch happened")
    # no_of_swaps = 0
    # for i in range(m):
    #     if (permutation_list[i] != i):
    #         no_of_swaps = no_of_swaps + 1
    #         perm(permutation_list,permutation_list[i],i)

    return det_A


                     

