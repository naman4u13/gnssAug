package com.gnssAug.helper.lambdaNew.Estimators;
import org.ejml.simple.SimpleMatrix;

import com.gnssAug.Android.constants.GnssDataConfig;
import com.gnssAug.helper.lambdaNew.ComputeSR_IBexact;
import com.gnssAug.helper.lambdaNew.ComputeVariance;
import com.gnssAug.helper.lambdaNew.ComputeVariance.VarianceResult;
import com.gnssAug.helper.lambdaNew.GammaIncompleteInverse;
import com.gnssAug.helper.lambdaNew.ComputeSR_IBexact.SR_IB;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorBIE.EstimatorBIEResult;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorILS.ILSResult;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorPAR_FFRT.PARResult_FFRT;

import java.util.Arrays;



/**
 * LAMBDA 4.0 | Partial Ambiguity Resolution (PAR) estimation based on ILS
 * This class provides a method to compute a partial 'fixed' solution based on the best integer 
 * least-squares solutions for the most precise subset, given a minimum success rate threshold.
 * Multiple best candidates can be selected for the integer-fixed subset, conditioning the remaining components accordingly.
 *
 * -------------------------------------------------------------------------
 * INPUTS:
 *   aHat        Ambiguity float vector (column)
 *   LMat        LtDL-decomposition matrix L (lower unitriangular)
 *   dVec        LtDL-decomposition matrix D (diagonal elements)
 *   nCands      Number of best integer solutions          [DEFAULT = 1]
 *   minSR       Minimum success rate threshold             [DEFAULT = 99.5%]
 *   alphaBIE    Use BIE estimator instead if alpha > 0    [DEFAULT = 0]
 *
 * OUTPUTS:
 *   aPAR        Partially 'fixed' solution given a minimum success rate
 *   nFixed      Number of fixed ambiguity (most precise) components
 *   SR_PAR      Success rate of ambiguity (most precise) subset 
 *
 * DEPENDENCIES:
 *   computeSR_IBexact
 *   estimatorILS
 *   estimatorBIE
 *
 * REFERENCES:
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
public class EstimatorPAR {

    /**
     * Computes a partially 'fixed' solution based on integer least-squares solutions.
     *
     * @param aHat    Ambiguity float vector (column)
     * @param LMat    LtDL-decomposition matrix L (lower unitriangular)
     * @param dVec    LtDL-decomposition matrix D (diagonal elements)
     * @param nCands  Number of best integer solutions          [DEFAULT = 1]
     * @param minSR   Minimum success rate threshold             [DEFAULT = 0.995]
     * @param alphaBIE Use BIE estimator instead if alpha > 0    [DEFAULT = 0]
     * @return        PARResult object containing aPAR, nFixed, SR_PAR
     * @throws IllegalArgumentException if number of inputs is insufficient
     */
    public static PARResult estimatorPAR(SimpleMatrix aHat, SimpleMatrix LMat, double[] dVec, 
                                         Integer nCands, Double minSR, Double alphaBIE) {
        // Problem dimensionality
        int nn = aHat.numRows();
        VarianceResult varRes = null;
        // Check number of input arguments and set default values if necessary
        if (aHat == null || LMat == null || dVec == null) {
            throw new IllegalArgumentException("ATTENTION: number of inputs is insufficient!");
        }
        
        if (nCands == null) {
            nCands = 1;          // Number of integer candidates for the partial fix
        }
        
        if (minSR == null) {
            minSR = 0.995;       // Default minimum success rate threshold
        }
        
        if (alphaBIE == null) {
            alphaBIE = 0.0;      // By default, use ILS estimator
        }
        
        // Compute success rate for IB (exact formulation)
        SR_IB srResult = ComputeSR_IBexact.computeSR_IBexact(dVec);
        double SR_IB = srResult.getSR();
        double[] SR_IB_cumul = srResult.getSR_cumul();
        
        int kk_PAR;
        double SR_PAR;
        int nFixed;
        SimpleMatrix aPAR;
        
        // Check largest subset above SR threshold
        if (SR_IB >= minSR) {
            // Full AR       (fixed solution)
            kk_PAR = 1;
            SR_PAR = SR_IB;
            nFixed = nn;
            SimpleMatrix qMat_subset = LMat.transpose().mult(SimpleMatrix.diag(dVec)).mult(LMat);
			varRes = ComputeVariance.computeVariance(qMat_subset, 1, 0, null,(int) GnssDataConfig.nSamplesMC);
        } else {
            // Find the first index where cumulative SR meets or exceeds minSR
            kk_PAR = -1;
            for (int i = 0; i < SR_IB_cumul.length; i++) {
                if (SR_IB_cumul[i] >= minSR) {
                    kk_PAR = i + 1; // MATLAB indices start at 1
                    break;
                }
            }
            
            if (kk_PAR != -1) {
                // Partial AR    (fixed solution)
                SR_PAR = SR_IB_cumul[kk_PAR - 1];
                nFixed = nn - (kk_PAR - 1);
            } else {
                // No AR         (float solution)
                aPAR  = aHat;
                SR_PAR = SR_IB;
                nFixed = 0;
                return new PARResult(aPAR, nFixed, SR_PAR,null);
            }
        }
        
        // Find fixed solution of subset {II} with sufficiently high success rate
        SimpleMatrix a_fix_PAR;
        if (alphaBIE > 0.0 && alphaBIE < 1.0) {
            double Chi2_BIE = 2.0 * GammaIncompleteInverse.gammaincinv(1 - alphaBIE, nFixed / 2.0);
            // NOTE: This uses the inverse of the regularized gamma function from Apache Commons Math
            // Ensure that Apache Commons Math library is included in the project dependencies
            
            // Call BIE-estimator (recursive implementation)
            SimpleMatrix aHat_subset = aHat.extractMatrix(kk_PAR - 1, nn, 0, 1);
            SimpleMatrix LMat_subset = LMat.extractMatrix(kk_PAR - 1, nn, kk_PAR - 1, nn);
            double[] dVec_subset = Arrays.copyOfRange(dVec, kk_PAR - 1, nn);
            EstimatorBIEResult bieResult = EstimatorBIE.estimatorBIE(aHat_subset, LMat_subset, dVec_subset, Chi2_BIE,null);
            a_fix_PAR = bieResult.getaBIE();
            // NOTE: this is an experimental PAR (BIE) approach still based on the 
            // SR criterion. We suggest to use "minSR = 0.50" & "alphaBIE = 1e-6", 
            // or to check the alternative implementation in 'estimatorPAR_BIE.m'
        } else {
            // Call ILS-estimator (search-and-shrink)
            SimpleMatrix aHat_subset = aHat.extractMatrix(kk_PAR - 1, nn, 0, 1);
            SimpleMatrix LMat_subset = LMat.extractMatrix(kk_PAR - 1, nn, kk_PAR - 1, nn);
            double[] dVec_subset = Arrays.copyOfRange(dVec, kk_PAR - 1, nn);
            ILSResult ilsResult = EstimatorILS.estimatorILS(aHat_subset, LMat_subset, dVec_subset, nCands);
            SimpleMatrix qMat_subset = LMat_subset.transpose().mult(SimpleMatrix.diag(dVec_subset)).mult(LMat_subset);
			varRes = ComputeVariance.computeVariance(qMat_subset, 1, 0, null,(int) GnssDataConfig.nSamplesMC);
            a_fix_PAR = ilsResult.getAFix();
        }
        
     // Extract float ambiguities before and after kkPAR
        SimpleMatrix aHat1 = aHat.extractMatrix(0, kk_PAR, 0, 1);          // a_hat(1:kk_PAR-1) in MATLAB
        SimpleMatrix aHat2 = aHat.extractMatrix(kk_PAR, aHat.numRows(), 0, 1); // a_hat(kk_PAR:end) in MATLAB

        SimpleMatrix QMat = LMat.transpose().mult(SimpleMatrix.diag(dVec)).mult(LMat);
        SimpleMatrix QMat_11 = QMat.extractMatrix(0, kk_PAR,0, kk_PAR);
        SimpleMatrix QMat_22 = QMat.extractMatrix(kk_PAR, nn,kk_PAR, nn);
        SimpleMatrix QMat_12 = QMat.extractMatrix(0, kk_PAR,kk_PAR, nn);
        SimpleMatrix QMat_21 = QMat_12.transpose();
        SimpleMatrix Q_fix_PAR = varRes.getVariance();
        
        SimpleMatrix a_cond_PAR = aHat1.minus(QMat_12.mult(QMat_22.invert()).mult(aHat2.minus(a_fix_PAR)));
       
        
        SimpleMatrix term1 = QMat_12.mult(QMat_22.invert()).mult(QMat_21);
        SimpleMatrix term2 = QMat_12.mult(QMat_22.invert()).mult(Q_fix_PAR).mult(QMat_22.invert()).mult(QMat_21);
        SimpleMatrix Q_cond_PAR = QMat_11.minus(term1).plus(term2);
        // Concatenate a_cond_PAR and a_fix_PAR vertically
        aPAR = new SimpleMatrix(nn, 1);

        // Set the first part of aPAR to aCondPAR
        aPAR.insertIntoThis(0, 0, a_cond_PAR);

        // Set the second part of aPAR to aFixPAR
        aPAR.insertIntoThis(a_cond_PAR.numRows(), 0, a_fix_PAR);

        SimpleMatrix QPAR = new SimpleMatrix(nn, nn);
        QPAR.insertIntoThis(0, 0, Q_cond_PAR);
        QPAR.insertIntoThis(kk_PAR, kk_PAR, Q_fix_PAR);
        
        // aPAR now contains the vertically concatenated result
    
		return new PARResult(aPAR, nFixed, SR_PAR,QPAR);    }
    

    
    /**
     * Class to hold the results of the PAR estimation.
     */
    public static class PARResult {
        private SimpleMatrix aPAR;
        private int nFixed;
        private double SR_PAR;
        private SimpleMatrix QPAR;
        
        public PARResult(SimpleMatrix aPAR, int nFixed, double SR_PAR,SimpleMatrix QPAR) {
            this.aPAR = aPAR;
            this.nFixed = nFixed;
            this.SR_PAR = SR_PAR;
        }

        public SimpleMatrix getaPAR() {
            return aPAR;
        }

        public int getnFixed() {
            return nFixed;
        }

        public double getSR_PAR() {
            return SR_PAR;
        }
        public SimpleMatrix getQPAR() {
			return QPAR;
		}
    }
    
    
}
