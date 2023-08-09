This repository hosts a collection of code files that implement various fundamental techniques and algorithms in computational linear algebra

# Topics Overview

## Basic Linear Algebra
Basic linear algebra involves the study of vector spaces, linear transformations, and matrices. It forms the foundation of various mathematical and computational techniques used in solving systems of linear equations, eigenvalue problems, and more complex mathematical operations.

## Gram-Schmidt
Gram-Schmidt is an orthogonalization process that transforms a set of linearly independent vectors into an orthogonal orthonormal basis. It's commonly used to simplify calculations and improve stability when working with vectors in a high-dimensional space.

## Householder Methods
Householder methods are numerical algorithms used for matrix transformations and solving linear systems. They utilize Householder reflections, which are orthogonal transformations that enable the reduction of matrices to more structured forms, such as upper triangular matrices.

## Backwards Stability
Backwards stability is a concept in numerical analysis that pertains to the behavior of an algorithm when it is subjected to small input perturbations. An algorithm is considered backward stable if its computed result is the exact result of a slightly perturbed input.

## Condition Number
The condition number is a numerical measure of the sensitivity of a mathematical problem to changes in its input. For linear systems, it indicates how much the solution can change due to small changes in the coefficients. A high condition number implies a problem that is ill-conditioned and susceptible to numerical instability.

## Computing Eigenvectors
Computing eigenvectors involves finding the vectors that remain unchanged in direction when a linear transformation (represented by a matrix) is applied. It's a crucial task in various scientific and engineering applications, such as understanding dynamic systems and analyzing data.

## Pure QR Methods
Pure QR methods are iterative techniques used to compute the eigenvalues and eigenvectors of a matrix. These methods rely on the QR factorization of the matrix and converge towards the eigenvalues and eigenvectors through iterative transformations.

## LU Methods
LU (Lower-Upper) decomposition is a matrix factorization technique where a given matrix is decomposed into the product of a lower triangular matrix and an upper triangular matrix. LU methods are used for solving systems of linear equations and can provide insight into the matrix's properties.

## LUP Factorisations
LUP factorization is an extension of LU decomposition that includes permutation matrices. It allows for factorizing a matrix into three components: a permutation matrix, a lower triangular matrix, and an upper triangular matrix. LUP factorization is useful for numerical stability and solving linear systems with partial pivoting.

## GMRES
GMRES (Generalized Minimal Residual) is an iterative method used to solve large, sparse systems of linear equations. It seeks to find a solution that minimizes the residual error in a least-squares sense. GMRES is particularly effective for systems that can't be efficiently solved using direct methods due to their size and sparsity.
