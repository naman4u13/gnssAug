package com.gnssAug.helper.lambdaNew.Estimators;

import org.ejml.simple.SimpleMatrix;

/*%% LAMBDA 4.0 | Integer Rounding (IR) estimator
 * This function computes a 'fixed' solution based on an Integer Rounding 
 * (IR-)estimator, starting with a certain float ambiguity vector.
 *
 * -------------------------------------------------------------------------
 *_INPUTS:
 *   a_hat       Ambiguity float vector (column)
 *
 *_OUTPUTS:
 *   a_fix       Ambiguity fixed vector (column)
 *
 *_DEPENDENCIES:
 *   none
 *
 *_REFERENCES:
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



public class EstimatorIR {
    /**
     * MAIN FUNCTION
     *--------------------------------------------------------------------------
     * Round all (unconditioned) float ambiguity components
     *--------------------------------------------------------------------------
     *
     * @param aHat Ambiguity float vector (column)
     * @return aFix Ambiguity fixed vector (column)
     */
    public static SimpleMatrix estimatorIR(SimpleMatrix aHat) {
        //--------------------------------------------------------------------------
        // Round all (unconditioned) float ambiguity components
        SimpleMatrix aFix = new SimpleMatrix(aHat.numRows(), aHat.numCols());
        for (int i = 0; i < aHat.getNumElements(); i++) {
            aFix.set(i, Math.round(aHat.get(i)));
        }
        //--------------------------------------------------------------------------
        return aFix;
    }
}
