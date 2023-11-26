import numpy as np

def lu_decomposition(A):
    n = A.shape[0]
    L = np.eye(n)
    U = A.copy()

    for k in range(n-1):
        for i in range(k+1, n):
            factor = U[i, k] / U[k, k]
            L[i, k] = factor
            U[i, k:] -= factor * U[k, k:]

    return L, U

# Example usage:
A = np.array([[2, -1, 1], [3, 3, 9], [3, 3, 5]], dtype=float)
L, U = lu_decomposition(A)

print("Original matrix:")
print(A)

print("\nL matrix:")
print(L)

print("\nU matrix:")
print(U)

# Verify the decomposition by reconstructing A from L and U
A_reconstructed = np.dot(L, U)
print("\nReconstructed matrix from LU decomposition:")
print(A_reconstructed)

