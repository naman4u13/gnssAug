package com.gnssAug.helper.lambdaNew.Estimators;
/**
 * LAMBDA 4.0 | Integer Aperture (IA) using Fixed Failure-rate Ratio Test
 * This function computes the 'fixed' or 'float' solution based on the Integer
 * Aperture (IA) estimator, adopting a Fixed Failure-rate Ratio Test
 * (FFRT, [RD01]) based on lookup tables newly defined in [RD02], where
 * the maximum failure-rate threshold should be in the range [0.05-1%].
 *
 * -------------------------------------------------------------------------
 * INPUTS:
 *   aHat        Ambiguity float vector (column)
 *   lMat        Decomposition L matrix (lower unitriangular)
 *   dVec        Decomposition D matrix (vector of diagonal components)
 *   maxFR       Maximum failure-rate threshold, in the range [0.05-1%]
 *   muRatio     Ratio value used for an arbitrary ratio test    [OPTIONAL]
 *
 * OUTPUTS:
 *   aFix       IA solution based on FFRT (or an arbitrary ratio test)
 *   sqNorm     Squared-norm for the fixed solution
 *   nFixed     Number of fixed ambiguity components
 *
 * DEPENDENCIES:
 *   EstimatorILS.java
 *   ComputeSR_IBexact.java
 *   ComputeFFRTcoeff.java
 *
 * REFERENCES:
 *   [RD01] Verhagen, S., Teunissen, P.J.G. The ratio test for future GNSS 
 *       ambiguity resolution. GPS Solut 17, 535–548 (2013). 
 *       DOI: 10.1007/s10291-012-0299-z
 *
 *   [RD02] Hou, Yanqing, Sandra Verhagen, and Jie Wu. 2016. "An Efficient 
 *       Implementation of Fixed Failure-Rate Ratio Test for GNSS Ambiguity 
 *       Resolution" Sensors 16, no. 7: 945
 *       DOI: 10.3390/s16070945
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

import com.gnssAug.helper.lambdaNew.ComputeFFRTCoefficient;
import com.gnssAug.helper.lambdaNew.ComputeFFRTCoefficientOld;
import com.gnssAug.helper.lambdaNew.SuccessRate;
import com.gnssAug.helper.lambdaNew.SuccessRate.SRResult;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorILS.ILSResult;
import com.gnssAug.utility.Matrix;

public class EstimatorIA_FFRT {

    /**
     * Holds the results of the IA FFRT estimation.
     */
    public static class IAFFRTResult {
        private SimpleMatrix aFix;
        private double sqNorm;
        private int nFixed;

        public IAFFRTResult(SimpleMatrix aFix, double sqNorm, int nFixed) {
            this.aFix = aFix;
            this.sqNorm = sqNorm;
            this.nFixed = nFixed;
        }

        public SimpleMatrix getaFix() {
            return aFix;
        }

        public double getsqNorm() {
            return sqNorm;
        }

        public int getnFixed() {
            return nFixed;
        }
    }

    /**
     * Estimates the IA solution using Fixed Failure-rate Ratio Test.
     *
     * @param aHat    Ambiguity float vector (column)
     * @param lMat    Decomposition L matrix (lower unitriangular)
     * @param dVec    Decomposition D matrix (vector of diagonal components)
     * @param maxFR   Maximum failure-rate threshold, in the range [0.05-1%]
     * @param muRatio Ratio value used for an arbitrary ratio test [OPTIONAL]
     * @return EstimatorResult containing aFix, sqNorm, and nFixed
     * @throws IllegalArgumentException if the number of inputs is insufficient
     */
    public static IAFFRTResult estimatorIA_FFRT(SimpleMatrix aHat, SimpleMatrix lMat, double[] dVec, Double maxFR, Double muRatio) {
        // Problem dimensionality
        int nn = aHat.numRows();

        // Check the number of input arguments
        if (aHat == null || lMat == null || dVec == null) {
            throw new IllegalArgumentException("ATTENTION: number of inputs is insufficient!");
        }

        // Set default maximum Fixed Failure Rate if maxFR is not provided
        double effectiveMaxFR = (maxFR != null) ? maxFR : 0.1 / 100;

        // Instead of a maximum failure-rate criterion, consider an arbitrary ratio
        // test based on the input "muRatio", while setting "maxFR" to zero.
        if (muRatio != null) {
            effectiveMaxFR = 0.0; // Assure FR >= maxFR for arbitrary ratio tests
        }

        // Compute two best solutions used for the Fixed Failure-rate Ratio Test
        ILSResult ilsResult = EstimatorILS.estimatorILS(aHat, lMat, dVec, 2);
        SimpleMatrix aFixTemp = ilsResult.getAFix();
        double[] sqNormTemp = ilsResult.getSqNorm();

        // Compute Success Rate (SR) and Failure Rate (FR)
        SRResult srResult = SuccessRate.computeSR_IBexact(dVec);
        double sr = srResult.getSR();
        double fr = 1.0 - sr;

        // Check the Failure Rate (FR) with respect to current maximum FR threshold
        if (fr < effectiveMaxFR) {
            // FR is below the threshold, so return the ILS solution
            SimpleMatrix aFix = aFixTemp.extractVector(false, 0);
            double sqNorm = sqNormTemp[0];
            int nFixed = nn;

            return new IAFFRTResult(aFix, sqNorm, nFixed);
        } else {
            // FR is over the threshold, compute the mu-value for IA
            double muValue;
            if (muRatio != null) {
                muValue = muRatio; // Arbitrary ratio test
            } else {
                muValue = ComputeFFRTCoefficient.computeFFRTcoeff(effectiveMaxFR, fr, nn); // Fixed-FR ratio test
            }

            // Perform the Ratio Test based on the mu-value
            if (sqNormTemp[0] / sqNormTemp[1] > muValue) {
                SimpleMatrix aFix = aHat;
                double sqNorm = 0.0;
                int nFixed = 0;

                return new IAFFRTResult(aFix, sqNorm, nFixed);
            } else {
                SimpleMatrix aFix = aFixTemp.extractVector(false, 0);
                double sqNorm = sqNormTemp[0];
                int nFixed = nn;

                return new IAFFRTResult(aFix, sqNorm, nFixed);
            }
        }
    }
}

