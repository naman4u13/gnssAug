package com.gnssAug.helper.lambdaNew;
/**
 * LAMBDA 4.0 | Check the main inputs for LAMBDA routine
 * This function checks the main inputs, i.e. variance-covariance matrix and
 * ambiguity vector, for the LAMBDA routine. The vc-matrix should be symmetric
 * positive-definite, while the vector dimensionality shall be compatible.
 *
 * -------------------------------------------------------------------------
 * _INPUTS:
 *   qaHat      Variance-covariance matrix of the original ambiguities
 *   aHat       Ambiguity float vector (column)
 *
 * _OUTPUTS:
 *   Error message is returned if one of the tests fails.
 *
 * _DEPENDENCIES:
 *   none
 *
 * _REFERENCES:
 *   none
 *
 * -------------------------------------------------------------------------
 * Copyright: Geoscience & Remote Sensing department @ TUDelft | 01/06/2024
 * Contact email:    LAMBDAtoolbox-CITG-GRS@tudelft.nl
 * -------------------------------------------------------------------------
 * Created by
 *   01/06/2024  - Lotfi Massarweh
 *       Implementation for LAMBDA 4.0 toolbox, based on LAMBDA 3.0
 *
 * Modified by
 *   dd/mm/yyyy  - Name Surname (author)
 *       >> Changes made in this new version
 * -------------------------------------------------------------------------
 */

import org.ejml.data.DMatrixRMaj;
import org.ejml.simple.SimpleMatrix;
import org.ejml.dense.row.CommonOps_DDRM;

import org.ejml.dense.row.factory.DecompositionFactory_DDRM;
import org.ejml.interfaces.decomposition.CholeskyDecomposition_F64;



public class Utilities {

	/**
     * Checks the inputs for the LAMBDA routine.
     *
     * @param QaHat Variance-covariance matrix of the original ambiguities (SimpleMatrix)
     * @param aHat  Ambiguity float vector (SimpleMatrix column vector)
     * @throws IllegalArgumentException if any of the checks fail
     */
    public static void checkMainInputs(SimpleMatrix QaHat, SimpleMatrix aHat) {
        // -------------------------------------------------------------------------
        // Test #1 on the variance-covariance matrix.
        // TEST 1a: Is the variance-covariance matrix "QaHat" symmetric?
        if (!isSymmetric(QaHat)) {
            throw new IllegalArgumentException("ATTENTION: Variance-covariance matrix needs to be symmetric!");
        }

        // TEST 1b: Is the variance-covariance matrix "QaHat" positive-definite?
        if (!isPositiveDefinite(QaHat)) {
            throw new IllegalArgumentException("ATTENTION: Variance-covariance matrix needs to be positive-definite!");
        }

        // No errors found? All tests #1 are passed!
        // -------------------------------------------------------------------------
        // Test #2 on the float ambiguity vector (if any).

        if (aHat != null) {
            // TEST 2a: Is the (float) ambiguity vector a column?
            if (aHat.numCols() != 1) {
                throw new IllegalArgumentException("ATTENTION: Float ambiguity vector needs to be a column vector!");
            }

            // TEST 2b: Do the ambiguity vector and variance-covariance matrix have compatible dimensions?
            if (aHat.numRows() != QaHat.numRows()) {
                throw new IllegalArgumentException("ATTENTION: Dimension mismatch between float vector and its variance-covariance matrix!");
            }
        }

        // No errors found? All tests #2 are passed!
        // -------------------------------------------------------------------------
    }

    /**
     * Checks if a matrix is symmetric.
     *
     * @param matrix The input matrix (SimpleMatrix)
     * @return true if the matrix is symmetric, false otherwise
     */
    private static boolean isSymmetric(SimpleMatrix matrix) {
        // Symmetry check: |Qa - Qa^T| < 1e-12
        SimpleMatrix difference = matrix.minus(matrix.transpose());
        return difference.elementMaxAbs() <= 1e-12;
    }

    /**
     * Checks if a matrix is positive-definite using Cholesky decomposition.
     *
     * @param matrix The input matrix (SimpleMatrix)
     * @return true if the matrix is positive-definite, false otherwise
     */
    private static boolean isPositiveDefinite(SimpleMatrix matrix) {
        // Create a copy of the matrix to prevent in-place modifications
        DMatrixRMaj denseMatrixCopy = matrix.copy().getDDRM();

        // Use EJML's factory to create a Cholesky decomposition object
        CholeskyDecomposition_F64<DMatrixRMaj> cholesky = DecompositionFactory_DDRM.chol(true);

        // Attempt to decompose the matrix; return true if successful
        return cholesky.decompose(denseMatrixCopy);
    }
}
