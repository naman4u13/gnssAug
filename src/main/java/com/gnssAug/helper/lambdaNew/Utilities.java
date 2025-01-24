package com.gnssAug.helper.lambdaNew;

import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.MatrixFeatures_DDRM;
import org.ejml.simple.SimpleMatrix;

/**
 * LAMBDA 4.0 | Check the main inputs for LAMBDA routine This class checks the
 * main inputs, i.e., variance-covariance matrix and ambiguity vector, for the
 * LAMBDA routine. The vc-matrix should be symmetric positive-definite, while
 * the vector dimensionality shall be compatible.
 *
 * -------------------------------------------------------------------------
 * INPUTS: qaHat Variance-covariance matrix of the original ambiguities aHat
 * Ambiguity float vector (column)
 *
 * OUTPUTS: Throws IllegalArgumentException if one of the tests fails.
 *
 * DEPENDENCIES: EJML library for matrix computations
 *
 * REFERENCES: none
 *
 * -------------------------------------------------------------------------
 * Copyright: Geoscience & Remote Sensing department @ TUDelft | 01/06/2024
 * Contact email: LAMBDAtoolbox-CITG-GRS@tudelft.nl
 * -------------------------------------------------------------------------
 * Created by 01/06/2024 - Lotfi Massarweh Implementation for LAMBDA 4.0
 * toolbox, based on LAMBDA 3.0
 *
 * Modified by dd/mm/yyyy - Name Surname (author) >> Changes made in this new
 * version
 * -------------------------------------------------------------------------
 */
public class Utilities {

	/**
	 * Checks the main inputs for the LAMBDA routine.
	 *
	 * @param qaHat Variance-covariance matrix of the original ambiguities
	 * @param aHat  Ambiguity float vector (column). Can be null if not provided.
	 * @throws IllegalArgumentException if any of the input checks fail
	 */
	public static void checkMainInputs(SimpleMatrix qaHat, SimpleMatrix aHat) {
		// --------------------------------------------------------------------------
		// Test #1 on the variance-covariance matrix.

		// TEST 1a: Is the variance-covariance matrix "qaHat" symmetric?
		if (!MatrixFeatures_DDRM.isSymmetric(qaHat.getDDRM(), 1e-12)) {
			throw new IllegalArgumentException("ATTENTION: variance-covariance matrix needs to be symmetric!");
		}

		// TEST 1b: Is the variance-covariance matrix "qaHat" positive-definite?
		org.ejml.dense.row.decomposition.chol.CholeskyDecompositionInner_DDRM cholesky = new org.ejml.dense.row.decomposition.chol.CholeskyDecompositionInner_DDRM();
		if (!cholesky.decompose(qaHat.getDDRM())) {
			throw new IllegalArgumentException("ATTENTION: variance-covariance matrix needs to be positive-definite!");
		}

		// >> No errors found? All tests #1 are passed!

		// --------------------------------------------------------------------------
		// Test #2 on the float ambiguity vector (if any).

		if (aHat != null) {
			// TEST 2a: Is the (float) ambiguity vector a column?
			if (aHat.numCols() != 1) {
				throw new IllegalArgumentException("ATTENTION: float ambiguity vector needs to be a column vector!");
			}

			// TEST 2b: Do the ambiguity vector & vc-matrix have compatible dimensions?
			if (aHat.numRows() != qaHat.numCols()) {
				throw new IllegalArgumentException(
						"ATTENTION: dimension mismatch between float vector and its vc-matrix");
			}
		}

		// >> No errors found? All tests #2 are passed!
	}

}
