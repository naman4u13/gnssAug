package com.gnssAug.helper.lambdaNew.Estimators;

/**
 * LAMBDA 4.0 | Integer Bootstrapping (IB) estimator
 * This class computes a 'fixed' solution based on the Integer Bootstrapping 
 * (IB-)estimator, starting with a certain float ambiguity vector.
 *
 * -------------------------------------------------------------------------
 * _INPUTS:
 *   aHat       Ambiguity float vector (column)
 *   lMat       LtDL-decomposition matrix L (lower unitriangular)
 *
 * _OUTPUTS:
 *   aFix       Ambiguity fixed vector (column)
 *   aCond      Ambiguity "conditioned" float vector (column)    [OPTIONAL]
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
 *   dd/mm/yyyy  - Name Surname author - email address
 *       >> Changes made in this new version
 * -------------------------------------------------------------------------
 */

import org.ejml.data.DMatrixRMaj;
import org.ejml.simple.SimpleMatrix;

import java.util.Arrays;

/**
 * Main class for the Integer Bootstrapping (IB) estimator.
 */
public class EstimatorIB {

    /**
     * Computes a 'fixed' solution based on the Integer Bootstrapping (IB-)estimator,
     * starting with a certain float ambiguity vector.
     *
     * @param aHat Ambiguity float vector (column)
     * @param lMat LtDL-decomposition matrix L (lower unitriangular)
     * @return An EstimatorResult containing the fixed ambiguity vector and the conditioned float vector
     */
    public static SimpleMatrix estimatorIB(SimpleMatrix aHat, SimpleMatrix lMat) {
        // Problem dimensionality
        int nn = aHat.numRows();

        // Initialize main vectors
        double[] aFix = new double[nn];                  // Integer-fixed ambiguity vector
        Arrays.fill(aFix, Double.NaN);                   // Initialize with NaN

        double[] aCond = new double[nn];                 // Float (conditioned) ambiguity vector
        for (int i = 0; i < nn; i++) {
            aCond[i] = aHat.get(i, 0);
        }

        aFix[nn - 1] = Math.round(aCond[nn - 1]);        // Fixed (conditioned) ambiguity last component

        // Auxiliary vector used to compute float conditioned ambiguities "aCond"
        double[] sumCond = new double[nn];               // Initialized to 0

        // Iterative cycle, conditioning (last-to-first)
        for (int ii = nn - 2; ii >= 0; ii--) {
            // Compute the i-th ambiguity conditioned on previous ones (ii+1 to nn)
            for (int j = 0; j <= ii; j++) {
                sumCond[j] -= lMat.get(ii + 1, j) * (aCond[ii + 1] - aFix[ii + 1]);
            }

            // Conditioned ambiguity i-th component
            aCond[ii] = aHat.get(ii, 0) + sumCond[ii];
            //          = aHat[ii] + sum_{j=i+1...n} L_ji * ( aCond_j - aFix_j );

            // Sequentially rounded (conditioned) i-th component
            aFix[ii] = Math.round(aCond[ii]);
        }

        return new SimpleMatrix(aFix.length, 1, true, aFix);
    }

}