package com.gnssAug.helper.lambdaNew.Estimators;

import org.ejml.simple.SimpleMatrix;

import com.gnssAug.helper.lambdaNew.ComputeFFRTCoefficient;
import com.gnssAug.helper.lambdaNew.GammaIncompleteInverse;
import com.gnssAug.helper.lambdaNew.SuccessRate;
import com.gnssAug.helper.lambdaNew.SuccessRate.SRResult;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorBIE.EstimatorBIEResult;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorILS.ILSResult;

import java.util.Arrays;

/**
 * LAMBDA 4.0 | Partial Ambiguity Resolution (PAR) estimation based on ILS This
 * class provides a method to compute a partial 'fixed' solution based on the
 * best integer least-squares solutions for the most precise subset, given a
 * minimum success rate threshold. Multiple best candidates can be selected for
 * the integer-fixed subset, conditioning the remaining components accordingly.
 *
 * -------------------------------------------------------------------------
 * INPUTS: aHat Ambiguity float vector (column) LMat LtDL-decomposition matrix L
 * (lower unitriangular) dVec LtDL-decomposition matrix D (diagonal elements)
 * nCands Number of best integer solutions [DEFAULT = 1] minSR Minimum success
 * rate threshold [DEFAULT = 99.5%] alphaBIE Use BIE estimator instead if alpha
 * > 0 [DEFAULT = 0]
 *
 * OUTPUTS: aPAR Partially 'fixed' solution given a minimum success rate nFixed
 * Number of fixed ambiguity (most precise) components SR_PAR Success rate of
 * ambiguity (most precise) subset
 *
 * DEPENDENCIES: computeSR_IBexact estimatorILS estimatorBIE
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
 * Modified by dd/mm/yyyy - Name Surname author - email address >> Changes made
 * in this new version
 * -------------------------------------------------------------------------
 */
public class EstimatorPAR_FFRT {

	/**
	 * Computes a partially 'fixed' solution based on integer least-squares
	 * solutions.
	 *
	 * @param aHat     Ambiguity float vector (column)
	 * @param LMat     LtDL-decomposition matrix L (lower unitriangular)
	 * @param dVec     LtDL-decomposition matrix D (diagonal elements)
	 * @param nCands   Number of best integer solutions [DEFAULT = 1]
	 * @param minSR    Minimum success rate threshold [DEFAULT = 0.995]
	 * @param alphaBIE Use BIE estimator instead if alpha > 0 [DEFAULT = 0]
	 * @return PARResult object containing aPAR, nFixed, SR_PAR
	 * @throws IllegalArgumentException if number of inputs is insufficient
	 */
	public static PARResult_FFRT estimatorPAR_FFRT(SimpleMatrix aHat, SimpleMatrix LMat, double[] dVec, Integer nCands,
			Double minSR, Double alphaBIE) {
		// Problem dimensionality
		int nn = aHat.numRows();

		// Check number of input arguments and set default values if necessary
		if (aHat == null || LMat == null || dVec == null) {
			throw new IllegalArgumentException("ATTENTION: number of inputs is insufficient!");
		}

		if (nCands == null) {
			nCands = 1; // Number of integer candidates for the partial fix
		}

		if (minSR == null) {
			minSR = 0.995; // Default minimum success rate threshold
		}

		if (alphaBIE == null) {
			alphaBIE = 0.0; // By default, use ILS estimator
		}

		// Compute success rate for IB (exact formulation)
		SRResult srResult = SuccessRate.computeSR_IBexact(dVec);
		double SR_IB = srResult.getSR();
		double[] SR_IB_cumul = srResult.getSR_cumul();

		int kk_PAR;
		double SR_PAR;
		int nFixed;
		SimpleMatrix aPAR;

		
		kk_PAR = findFirstAboveThreshold_RatioTest(aHat, LMat, dVec, SR_IB_cumul, minSR);
		if (kk_PAR == -1) {
			// No AR
			return new PARResult_FFRT(aHat, 0, SR_IB);
		}
		SR_PAR = SR_IB_cumul[kk_PAR];
		nFixed = nn - kk_PAR;
		
		
		// Find fixed solution of subset {II} with sufficiently high success rate
		SimpleMatrix a_fix_PAR;
		if (alphaBIE > 0.0 && alphaBIE < 1.0) {
			double Chi2_BIE = 2.0 * GammaIncompleteInverse.gammaincinv(1 - alphaBIE, nFixed / 2.0);
			
			// Call BIE-estimator (recursive implementation)
			SimpleMatrix aHat_subset = aHat.extractMatrix(kk_PAR, nn, 0, 1);
			SimpleMatrix LMat_subset = LMat.extractMatrix(kk_PAR, nn, kk_PAR, nn);
			double[] dVec_subset = Arrays.copyOfRange(dVec, kk_PAR, nn);
			EstimatorBIEResult bieResult = EstimatorBIE.estimatorBIE(aHat_subset, LMat_subset, dVec_subset, Chi2_BIE,
					null);
			a_fix_PAR = bieResult.getaBIE();
			// NOTE: this is an experimental PAR (BIE) approach still based on the
			// SR criterion. We suggest to use "minSR = 0.50" & "alphaBIE = 1e-6",
			// or to check the alternative implementation in 'estimatorPAR_BIE.m'
		} else {
			// Call ILS-estimator (search-and-shrink)
			SimpleMatrix aHat_subset = aHat.extractMatrix(kk_PAR, nn, 0, 1);
			SimpleMatrix LMat_subset = LMat.extractMatrix(kk_PAR, nn, kk_PAR, nn);
			double[] dVec_subset = Arrays.copyOfRange(dVec, kk_PAR, nn);
			ILSResult ilsResult = EstimatorILS.estimatorILS(aHat_subset, LMat_subset, dVec_subset, nCands);
			a_fix_PAR = ilsResult.getAFix();
		}

		
		SimpleMatrix LMat_subset_transpose = LMat.extractMatrix(kk_PAR, nn, kk_PAR, nn).transpose();
		SimpleMatrix aHat_subset = aHat.extractMatrix(kk_PAR, nn, 0, 1);
		SimpleMatrix term = aHat_subset.minus(a_fix_PAR);

		SimpleMatrix inv_LMat = LMat.extractMatrix(kk_PAR, nn, kk_PAR, nn).transpose().invert();
		SimpleMatrix multiplication = LMat_subset_transpose.mult(inv_LMat).mult(term);
		SimpleMatrix a_cond_PAR = aHat.extractMatrix(0, kk_PAR, 0, 1).minus(multiplication);

		// Return PAR solution(s)
		aPAR = a_cond_PAR.combine(0, a_cond_PAR.numCols(), a_fix_PAR);

		return new PARResult_FFRT(aPAR, nFixed, SR_PAR);
	}

	private static int findFirstAboveThreshold_RatioTest(SimpleMatrix aHat, SimpleMatrix LMat, double[] dVec,
			double[] srCumul, double minSR) {

		double maxFR = 0.1 / 100.0;
		int n = aHat.numRows();
		for (int i = 0; i < n; i++) {
			if (srCumul[i] >= minSR) {
				return i;
			}
			ILSResult ilsResult = EstimatorILS.estimatorILS(aHat.extractMatrix(i, n, 0, 1),
					LMat.extractMatrix(i, n, i, n), Arrays.copyOfRange(dVec, i, n), 2);
			double mu_value = ComputeFFRTCoefficient.computeFFRTcoeff(maxFR, 1 - srCumul[i], n - 1); // Fixed-FR ratio

			double[] sqNorms = ilsResult.getSqNorm();
			// Step 4: Perform the Ratio Test based on the computed mu-value
			if (sqNorms[0] / sqNorms[1] < mu_value) {
				// If the ratio test passes, return the index
				return i;
			}

		}
		return -1; // Should not occur if threshold logic is correct
	}

	/**
	 * Class to hold the results of the PAR estimation.
	 */
	public static class PARResult_FFRT {
		private SimpleMatrix aPAR;
		private int nFixed;
		private double SR_PAR;

		public PARResult_FFRT(SimpleMatrix aPAR, int nFixed, double SR_PAR) {
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
	}

}
