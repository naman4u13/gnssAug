package com.gnssAug.helper.lambdaNew;

import com.gnssAug.helper.lambdaNew.models.*;
import org.apache.commons.math3.linear.*;

public class DecomposeLtDL {

	 /**
     * Perform an LtDL decomposition on the ambiguity variance-covariance matrix.
     *
     * @param Q_mat Variance-covariance matrix of the ambiguities (assumed symmetric positive-definite).
     * @return A result object containing L (lower unitriangular matrix) and d (diagonal elements vector).
     * @throws IllegalArgumentException if the matrix is not positive-definite.
     */
    public static LtDLDecompositionResult decomposeLtDL(RealMatrix Q_mat) {
        // Ensure Q_mat is square
        int n = Q_mat.getRowDimension();
        if (n != Q_mat.getColumnDimension()) {
            throw new IllegalArgumentException("Q_mat must be a square matrix.");
        }

        // Copy the input matrix to avoid modifying the original
        RealMatrix Q = Q_mat.copy();

        // Perform in-place iterations for the LtDL decomposition
        for (int k = n - 1; k >= 0; k--) {
            // Normalize the row Q[k, 0:k-1] by Q[k, k]
            double diagK = Q.getEntry(k, k);
            if (diagK <= 1e-12) {
                throw new IllegalArgumentException("Input matrix is not positive-definite or numerical errors occurred.");
            }

            for (int j = 0; j < k; j++) {
                Q.setEntry(k, j, Q.getEntry(k, j) / diagK);
            }

            // Update the upper-left block of the matrix
            for (int i = 0; i < k; i++) {
                for (int j = 0; j < k; j++) {
                    double value = Q.getEntry(i, j) - Q.getEntry(k, i) * diagK * Q.getEntry(k, j);
                    Q.setEntry(i, j, value);
                }
            }
        }

        // Extract the L matrix (lower unitriangular) and d vector (diagonal elements)
        RealMatrix L = MatrixUtils.createRealIdentityMatrix(n); // Start with identity matrix
        double[] d = new double[n];

        for (int i = 0; i < n; i++) {
            d[i] = Q.getEntry(i, i);
            for (int j = 0; j < i; j++) {
                L.setEntry(i, j, Q.getEntry(i, j));
            }
        }

        // Return results encapsulated in a result object
        return new LtDLDecompositionResult(L, new ArrayRealVector(d));
    }
	
}
