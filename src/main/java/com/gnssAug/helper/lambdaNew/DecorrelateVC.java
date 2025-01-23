package com.gnssAug.helper.lambdaNew;
import org.apache.commons.math3.linear.*; // For matrix operations

import com.gnssAug.helper.lambdaNew.models.*;

public class DecorrelateVC {

	/**
     * Decorrelate the ambiguity variance-covariance matrix by a Z-transformation.
     *
     * @param Qa_hat Variance-covariance matrix of the original ambiguities.
     * @param a_hat  Ambiguity float vector (optional). Pass null if not needed.
     * @return A result object containing the decorrelated matrices and vectors.
     */
	
	public static DecorrelateResult decorrelateVC(double[][] Qa_hat, double[] a_hat) {
		return decorrelateVC(MatrixUtils.createRealMatrix(Qa_hat), MatrixUtils.createRealVector(a_hat));
		
	}
    public static DecorrelateResult decorrelateVC(RealMatrix Qa_hat, RealVector a_hat) {
        // Step 1: Perform LtDL decomposition on Qa_hat
        LtDLDecompositionResult decomposition = DecomposeLtDL.decomposeLtDL(Qa_hat);
        RealMatrix La_mat = decomposition.getL();
        RealVector da_vec = decomposition.getD();

        // Step 2: Perform Z-transformation to reduce and order conditional variances
        TransformZResult zTransform = TransformZ.transformZ(La_mat, da_vec,null);
        RealMatrix Lz_mat = zTransform.getL();
        RealVector dz_vec = zTransform.getD();
        RealMatrix iZt_mat = zTransform.getIZt();

        // Step 3: Compute Qz_hat = Lz_mat' * diag(dz_vec) * Lz_mat
        RealMatrix Qz_hat = Lz_mat.transpose()
                .multiply(diagonalMatrixFromVector(dz_vec))
                .multiply(Lz_mat);

        // Optional outputs
        RealMatrix Z_mat = null;
        RealVector z_hat = null;

        if (a_hat != null) {
            // Retrieve the Z-transformation matrix: Z_mat = round(inv(iZt_mat.transpose()))
            Z_mat = new LUDecomposition(iZt_mat.transpose()).getSolver().getInverse();
            Z_mat = roundMatrix(Z_mat);

            // Transform the ambiguity float vector: z_hat = Z_mat' * a_hat
            z_hat = Z_mat.transpose().operate(a_hat);
        }

        // Return all outputs encapsulated in a result object
        return new DecorrelateResult(Qz_hat, Lz_mat, dz_vec, iZt_mat, Z_mat, z_hat);
    }

    /**
     * Helper method to create a diagonal matrix from a vector.
     *
     * @param vector The input vector.
     * @return A diagonal matrix with the vector elements as diagonal entries.
     */
    private static RealMatrix diagonalMatrixFromVector(RealVector vector) {
        int n = vector.getDimension();
        RealMatrix diagMatrix = MatrixUtils.createRealDiagonalMatrix(vector.toArray());
        return diagMatrix;
    }

    /**
     * Helper method to round all elements of a matrix.
     *
     * @param matrix The input matrix.
     * @return A matrix with all elements rounded to the nearest integer.
     */
    private static RealMatrix roundMatrix(RealMatrix matrix) {
        int rows = matrix.getRowDimension();
        int cols = matrix.getColumnDimension();
        RealMatrix roundedMatrix = MatrixUtils.createRealMatrix(rows, cols);

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                roundedMatrix.setEntry(i, j, Math.round(matrix.getEntry(i, j)));
            }
        }

        return roundedMatrix;
    }
}
