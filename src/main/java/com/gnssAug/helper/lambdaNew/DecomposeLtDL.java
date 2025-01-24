package com.gnssAug.helper.lambdaNew;


import org.ejml.UtilEjml;
import org.ejml.simple.SimpleMatrix;

/**
 * LAMBDA 4.0 | Perform a LtDL-decomposition on the ambiguity vc-matrix
 * This class computes a LtDL decomposition given the ambiguity variance-
 * covariance matrix, which is assumed to be symmetric positive-definite.
 * 
 * -------------------------------------------------------------------------
 * COPYRIGHT: Geoscience & Remote Sensing department @ TUDelft | 01/06/2024
 * Contact email:    LAMBDAtoolbox-CITG-GRS@tudelft.nl
 * -------------------------------------------------------------------------
 * Created by
 *   01/06/2024  - Lotfi Massarweh
 *       Implementation for LAMBDA 4.0 toolbox
 * 
 * Modified by
 *   dd/mm/yyyy  - Name Surname author
 *       >> Changes made in this new version
 * -------------------------------------------------------------------------
 */
public class DecomposeLtDL {

    /**
     * Performs a LtDL-decomposition on the given covariance matrix.
     *
     * @param QMat Variance-covariance matrix of the ambiguities
     * @return An object containing the L matrix and d vector of the decomposition
     * @throws IllegalArgumentException if the input matrix is not positive-definite
     */
	public static DecompositionResult decomposeLtDL(SimpleMatrix QMat) {
	    // Problem dimensionality
	    int nn = QMat.numRows();

	    // In-place iterations (i.e., last-to-first) for the LtDL-decomposition
	    for (int kk = nn - 1; kk >= 0; kk--) {
	        double Q_kk_kk = QMat.get(kk, kk);
	        for (int j = 0; j < kk; j++) {
	            double value = QMat.get(kk, j) / Q_kk_kk;
	            QMat.set(kk, j, value);
	        }
	        for (int i = 0; i < kk; i++) {
	            for (int j = 0; j < kk; j++) {
	                double updatedValue = QMat.get(i, j) - QMat.get(kk, j) * Q_kk_kk * QMat.get(kk, i);
	                QMat.set(i, j, updatedValue);
	            }
	        }
	    }

	    // Extract main outputs
	    SimpleMatrix LMat = new SimpleMatrix(nn, nn);
	    for (int i = 0; i < nn; i++) {
	        for (int j = 0; j < nn; j++) {
	            if (i > j) {
	                LMat.set(i, j, QMat.get(i, j));
	            } else if (i == j) {
	                LMat.set(i, j, 1.0);
	            } else {
	                LMat.set(i, j, 0.0);
	            }
	        }
	    }

	    double[] dVec = new double[nn];
	    for (int i = 0; i < nn; i++) {
	        dVec[i] = QMat.get(i,i);
	    }

	    // Check positive-definiteness following the decomposition
	    for (int i = 0; i < nn; i++) {
	        if (dVec[i] < 1e-12) {
	            throw new IllegalArgumentException("ATTENTION: the input vc-matrix is not positive-definite or numerical errors affect the LtDL-decomposition!");
	        }
	    }

	    return new DecompositionResult(LMat, dVec);
	}

    /**
     * Class to hold the result of the LtDL-decomposition.
     */
    public static class DecompositionResult {
        private final SimpleMatrix LMat;
        private final double[] dVec;

        /**
         * Constructs a DecompositionResult with the given L matrix and d vector.
         *
         * @param LMat The lower unitriangular matrix L
         * @param dVec The diagonal elements vector D
         */
        public DecompositionResult(SimpleMatrix LMat, double[] dVec) {
            this.LMat = LMat;
            this.dVec = dVec;
        }

        /**
         * Gets the lower unitriangular matrix L.
         *
         * @return The L matrix
         */
        public SimpleMatrix getLMat() {
            return LMat;
        }

        /**
         * Gets the diagonal elements vector D.
         *
         * @return The d vector
         */
        public double[] getDVec() {
            return dVec;
        }
    }
}

