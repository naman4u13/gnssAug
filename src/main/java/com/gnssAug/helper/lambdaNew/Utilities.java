package com.gnssAug.helper.lambdaNew;

public class Utilities {

    /**
     * Checks the main inputs for the LAMBDA routine.
     * Validates the variance-covariance matrix (Qa_hat) and the float ambiguity vector (a_hat).
     *
     * @param Qa_hat Variance-covariance matrix of the original ambiguities.
     * @param a_hat  Ambiguity float vector (column vector).
     * @throws IllegalArgumentException If any of the checks fail.
     */
    public static void checkMainInputs(double[][] Qa_hat, double[] a_hat) {
        // Test #1: Validate the variance-covariance matrix
        validateVarianceCovarianceMatrix(Qa_hat);

        // Test #2: Validate the ambiguity vector (if provided)
        if (a_hat != null) {
            validateAmbiguityVector(Qa_hat, a_hat);
        }

        // If no errors are thrown, all tests are passed!
    }

    /**
     * Test #1: Validates the variance-covariance matrix.
     * The matrix must be symmetric and positive-definite.
     *
     * @param Qa_hat Variance-covariance matrix.
     * @throws IllegalArgumentException If the matrix is not symmetric or not positive-definite.
     */
    private static void validateVarianceCovarianceMatrix(double[][] Qa_hat) {
        // Test 1a: Check if the matrix is symmetric
        if (!isSymmetric(Qa_hat)) {
            throw new IllegalArgumentException("ATTENTION: Variance-covariance matrix needs to be symmetric!");
        }

        // Test 1b: Check if the matrix is positive-definite
        if (!isPositiveDefinite(Qa_hat)) {
            throw new IllegalArgumentException("ATTENTION: Variance-covariance matrix needs to be positive-definite!");
        }
    }

    /**
     * Test #2: Validates the ambiguity vector.
     * Ensures the vector is a column and is compatible with the variance-covariance matrix dimensions.
     *
     * @param Qa_hat Variance-covariance matrix.
     * @param a_hat  Ambiguity float vector.
     * @throws IllegalArgumentException If the vector is not a column or dimensions mismatch.
     */
    private static void validateAmbiguityVector(double[][] Qa_hat, double[] a_hat) {
        // Test 2a: Check if the vector is a column
        if (a_hat.length == 0) {
            throw new IllegalArgumentException("ATTENTION: Float ambiguity vector needs to be a column vector!");
        }

        // Test 2b: Check for dimensional compatibility
        if (Qa_hat.length != a_hat.length) {
            throw new IllegalArgumentException("ATTENTION: Dimension mismatch between float vector and its variance-covariance matrix!");
        }
    }

    /**
     * Helper method to check if a matrix is symmetric.
     *
     * @param matrix Input matrix.
     * @return True if the matrix is symmetric, false otherwise.
     */
    private static boolean isSymmetric(double[][] matrix) {
        int n = matrix.length;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (Math.abs(matrix[i][j] - matrix[j][i]) > 1e-12) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * Helper method to check if a matrix is positive-definite.
     * Uses Cholesky decomposition to verify.
     *
     * @param matrix Input matrix.
     * @return True if the matrix is positive-definite, false otherwise.
     */
    private static boolean isPositiveDefinite(double[][] matrix) {
        try {
            choleskyDecomposition(matrix);
            return true;
        } catch (IllegalArgumentException e) {
            return false;
        }
    }

    /**
     * Performs Cholesky decomposition on a matrix.
     * Throws an exception if the matrix is not positive-definite.
     *
     * @param matrix Input matrix.
     * @return Lower triangular matrix resulting from the decomposition.
     * @throws IllegalArgumentException If the matrix is not positive-definite.
     */
    private static double[][] choleskyDecomposition(double[][] matrix) {
        int n = matrix.length;
        double[][] L = new double[n][n];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j <= i; j++) {
                double sum = 0.0;

                for (int k = 0; k < j; k++) {
                    sum += L[i][k] * L[j][k];
                }

                if (i == j) {
                    double value = matrix[i][i] - sum;
                    if (value <= 0) {
                        throw new IllegalArgumentException("Matrix is not positive-definite!");
                    }
                    L[i][j] = Math.sqrt(value);
                } else {
                    L[i][j] = (matrix[i][j] - sum) / L[j][j];
                }
            }
        }

        return L;
    }

    /**
     * Rounds each element of a double array towards zero (equivalent to MATLAB's fix()).
     *
     * @param inputArray The array of doubles to be rounded towards zero.
     * @return A new array with each element rounded towards zero as double values.
     */
    public static double[] fix(double[] inputArray) {
        double[] roundedArray = new double[inputArray.length];
        for (int i = 0; i < inputArray.length; i++) {
            // Truncate towards zero
            roundedArray[i] = inputArray[i] < 0 ? Math.ceil(inputArray[i]) : Math.floor(inputArray[i]);
        }
        return roundedArray;
    }
    
}
