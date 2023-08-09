import numpy as np
import timeit
import numpy.random as random

# pre-construct a matrix in the namespace to use in tests
random.seed(1651)
A0 = random.randn(500, 500)
x0 = random.randn(500)


def basic_matvec(A, x):
    """
    Elementary matrix-vector multiplication.

    :param A: an mxn-dimensional numpy array
    :param x: an n-dimensional numpy array

    returns an m-dimensional numpy array which is the product of A with x

    This should be implemented using a double loop over the entries of A

    :return b: m-dimensional numpy array
    """
    b = np.array([])
    med = 0
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            med = med + A[i,j] * x[j]
        b = np.append(b,med)
        med = 0
    return b
 #   raise NotImplementedError


def column_matvec(A, x):
    """
    Matrix-vector multiplication using the representation of the product
    Ax as linear combinations of the columns of A, using the entries in 
    x as coefficients.


    :param A: an mxn-dimensional numpy array
    :param x: an n-dimensional numpy array

    :return b: an m-dimensional numpy array which is the product of A with x

    This should be implemented using a single loop over the entries of x
    """
    b = 0
    for i in range(len(x)):
        b = b + x[i]*A[:,i]
    return b
   # raise NotImplementedError


def timeable_basic_matvec():
    """
    Doing a matvec example with the basic_matvec that we can
    pass to timeit.
    """

    b = basic_matvec(A0, x0) # noqa


def timeable_column_matvec():
    """
    Doing a matvec example with the column_matvec that we can
    pass to timeit.
    """

    b = column_matvec(A0, x0) # noqa


def timeable_numpy_matvec():
    """
    Doing a matvec example with the builtin numpy matvec so that
    we can pass to timeit.
    """

    b = A0.dot(x0) # noqa


def time_matvecs():
    """
    Get some timings for matvecs.
    """

    print("Timing for basic_matvec")
    print(timeit.Timer(timeable_basic_matvec).timeit(number=1))
    print("Timing for column_matvec")
    print(timeit.Timer(timeable_column_matvec).timeit(number=1))
    print("Timing for numpy matvec")
    print(timeit.Timer(timeable_numpy_matvec).timeit(number=1))

def outerprod(u1,v1):
    """
    Returns the outter product of two vectors u1*v1^*

    : param u1: m-dimensional numpy array 
    : param v1: n-dimensional numpy array
    """
    v1_conj = np.conjugate(v1)
    output_mat = np.empty((0,v1.shape[0]))
    holding_row = np.array([])
    for i in range(u1.shape[0]):
        for j in range(v1.shape[0]):
            holding_row = np.append(holding_row,u1[i]*v1_conj[j])  
        output_mat = np.append(output_mat, np.array([holding_row]), axis=0)
        holding_row = np.array([])
    return output_mat

def rank2(u1, u2, v1, v2):
    """
    Return the rank2 matrix A = u1*v1^* + u2*v2^*.

    :param u1: m-dimensional numpy array
    :param u2: m-dimensional numpy array
    :param v1: n-dimensional numpy array
    :param v2: n-dimensional numpy array
    """
    # B = u1*v1^*, C = u2*v2^*

    v1_conj = np.conjugate(v1)
    v2_conj = np.conjugate(v2)
    B = outerprod(u1,v1_conj)
    C = outerprod(u2,v2_conj)


#   raise NotImplementedError

    A = B + C 
    #fails test, ask 
    return A


def rank1pert_inv(u, v):
    """
    Return the inverse of the matrix A = I + uv^*, where I
    is the mxm dimensional identity matrix, with

    :param u: m-dimensional numpy array
    :param v: m-dimensional numpy array
    """
    # Inverse is of form A^-1 = I + alpha*uv^*
    # alpha of form:  alpha = -1/(1+v*u)

    #v conj
    v_conj = np.conjugate(v)

    #Identity of correct dimension
    I = np.identity(u.shape[0])

    alpha = -1/(1+ np.inner(v_conj,u))
    Ainv = I + alpha*np.outer(u,v_conj)

#    raise NotImplementedError

    return Ainv


def ABiC(Ahat, xr, xi):
    """Return the real and imaginary parts of z = A*x, where A = B + iC
    with

    :param Ahat: an mxm-dimensional numpy array with Ahat[i,j] = B[i,j] \
    for i>=j and Ahat[i,j] = C[i,j] for i<j.

    :return zr: m-dimensional numpy arrays containing the real part of z.
    :return zi: m-dimensional numpy arrays containing the imaginary part of z.
    """

    raise NotImplementedError

    return zr, zi
