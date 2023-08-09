# Iterative Methods in Numerical Linear Algebra

Iterative methods are essential tools for approximating eigenvalues and eigenvectors of matrices. This README provides insights into key iterative methods used in numerical linear algebra, including power iteration, inverse iteration, Rayleigh quotient iteration, and pure QR iteration.

## Power Iteration

Power iteration is a fundamental iterative method for computing the dominant eigenvalue and its corresponding eigenvector of a matrix. The method involves iteratively multiplying the matrix by a vector and normalizing the result. The eigenvector converges to the dominant eigenvector, and the eigenvalue converges to the dominant eigenvalue. Power iteration is efficient for large matrices with a well-defined dominant eigenvalue.

## Inverse Iteration

Inverse iteration is an extension of power iteration that allows for the approximation of eigenvalues other than the dominant eigenvalue. Inverse iteration involves solving a linear system in each iteration using the matrix factorization. This method converges to the eigenvalue closest to the specified shift. Inverse iteration is particularly useful for finding eigenvalues in the vicinity of a given estimate.

## Rayleigh Quotient Iteration

Rayleigh quotient iteration is an improvement upon power iteration and inverse iteration. It utilizes the Rayleigh quotient—a measure of the eigenvalue relative to a given vector—as a shift in the iterations. This approach accelerates the convergence to the desired eigenvalue and enhances the accuracy of the approximation.

## Pure QR Iteration

Pure QR iteration is an iterative method for computing all eigenvalues and eigenvectors of a matrix. It involves iteratively applying QR factorization to the matrix and updating the factors. This process diagonalizes the matrix, revealing its eigenvalues on the diagonal. Pure QR iteration provides an efficient way to compute the full spectrum of eigenvalues and eigenvectors.

## Conclusion

Iterative methods are powerful techniques for approximating eigenvalues and eigenvectors, with applications in various fields of numerical linear algebra. Power iteration, inverse iteration, Rayleigh quotient iteration, and pure QR iteration offer flexible and efficient ways to extract valuable information from matrices, enabling the analysis of complex systems and phenomena.


