package com.gnssAug.helper.lambdaNew.Estimators;

import org.ejml.simple.SimpleMatrix;

import com.gnssAug.helper.lambdaNew.ComputeFFRTCoefficient;
import com.gnssAug.helper.lambdaNew.ComputeSR_IBexact;
import com.gnssAug.helper.lambdaNew.ComputeSR_IBexact.SR_IB;
import com.gnssAug.helper.lambdaNew.ComputeVariance;
import com.gnssAug.helper.lambdaNew.EstimatorType;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorILS.ILSResult;
import com.gnssAug.Android.constants.GnssDataConfig;

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
	 * @throws Exception 
	 * @throws IllegalArgumentException if number of inputs is insufficient
	 */
	public static PARResult_FFRT estimatorPAR_FFRT(SimpleMatrix aHat, SimpleMatrix LMat, double[] dVec, Integer nCands,
			Double minSR, boolean estimateVar) throws Exception {
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

		// Compute success rate for IB (exact formulation)
		SR_IB srResult = ComputeSR_IBexact.computeSR_IBexact(dVec);
		double SR_IB = srResult.getSR();
		double[] SR_IB_cumul = srResult.getSR_cumul();

		int kk_PAR;
		double SR_PAR;
		int nFixed;
		SimpleMatrix aPAR;

		Object[] findFirstRes = findFirstAboveThreshold_RatioTest(aHat, LMat, dVec, SR_IB_cumul, minSR, estimateVar);

		if (findFirstRes == null) {
			// No AR
			SimpleMatrix qMat = LMat.transpose().mult(SimpleMatrix.diag(dVec))
					.mult(LMat);
			return new PARResult_FFRT(aHat, 0, SR_IB, new Object[] {qMat,0.0,1.0});
		}
		ILSResult ilsResult = (ILSResult) findFirstRes[0];
		Object[] stats = (Object[]) findFirstRes[1];
		kk_PAR = (int) findFirstRes[2];
		SR_PAR = SR_IB_cumul[kk_PAR];
		nFixed = nn - kk_PAR;

		// Find fixed solution of subset {II} with sufficiently high success rate
		SimpleMatrix a_fix_PAR;

		a_fix_PAR = ilsResult.getAFix().extractVector(false, 0);

		// Extract float ambiguities before and after kkPAR
		SimpleMatrix aHat1 = aHat.extractMatrix(0, kk_PAR, 0, 1); // a_hat(1:kk_PAR-1) in MATLAB
		SimpleMatrix aHat2 = aHat.extractMatrix(kk_PAR, aHat.numRows(), 0, 1); // a_hat(kk_PAR:end) in MATLAB

		SimpleMatrix QMat = LMat.transpose().mult(SimpleMatrix.diag(dVec)).mult(LMat);
		SimpleMatrix QMat_11 = QMat.extractMatrix(0, kk_PAR, 0, kk_PAR);
		SimpleMatrix QMat_22 = QMat.extractMatrix(kk_PAR, nn, kk_PAR, nn);
		SimpleMatrix QMat_12 = QMat.extractMatrix(0, kk_PAR, kk_PAR, nn);
		SimpleMatrix QMat_21 = QMat_12.transpose();
		SimpleMatrix Q_fix_PAR = stats != null ? (SimpleMatrix)stats[0] : SimpleMatrix.identity(nn - kk_PAR).scale(1e-10);

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
		if(stats==null)
		{
			stats = new Object[] {QPAR,0.0,1.0};
		}
		else
		{
			stats = new Object[] {QPAR,stats[1],stats[2]};
		}
		return new PARResult_FFRT(aPAR, nFixed, SR_PAR, stats);
	}

	private static Object[] findFirstAboveThreshold_RatioTest(SimpleMatrix aHat, SimpleMatrix LMat, double[] dVec,
			double[] srCumul, double minSR, boolean estimateVar) throws Exception {

		double maxFR = GnssDataConfig.PAR_FFRT_MAX_FR;
		int n = aHat.numRows();
		boolean flag = false;
		double muRatio = 1;
		for (int i = 0; i < n; i++) {
			ILSResult ilsResult = new EstimatorILS().estimatorILS(aHat.extractMatrix(i, n, 0, 1),
					LMat.extractMatrix(i, n, i, n), Arrays.copyOfRange(dVec, i, n), 2);
			if (srCumul[i] >= minSR) {
				flag = true;

			}
			if (!flag) {
				// Lambda.ratioinv(maxFR, 1 - srCumul[i], n-i);
				muRatio = ComputeFFRTCoefficient.computeFFRTcoeff(maxFR, 1 - srCumul[i], n - i); // Fixed-FR ratio
				double[] sqNorms = ilsResult.getSqNorm();

//				System.out.println("Combination count: " + (n - i));
//				System.out.println("failure rate :" + (1 - srCumul[i]));
//				System.out.println("MU :" + muRatio);
//				System.out.println("ratio :" + (sqNorms[0] / sqNorms[1]));

				// Step 4: Perform the Ratio Test based on the computed mu-value
				if (sqNorms[0] / sqNorms[1] < muRatio) {
					flag = true;
				}
			}
			if (flag) {
				SimpleMatrix LMat_subset = LMat.extractMatrix(i, n, i, n);
				double[] dVec_subset = Arrays.copyOfRange(dVec, i, n);
				SimpleMatrix qMat_subset = LMat_subset.transpose().mult(SimpleMatrix.diag(dVec_subset))
						.mult(LMat_subset);

				Object[] stats = null;
				if (estimateVar) {
					stats = ComputeVariance.computeVariance(qMat_subset, 2, 0, 1 / 100.0,
							(int) GnssDataConfig.nSamplesMC, null);
				}
//				if ((double)stats[1] == 0.0 && (double)stats[2] == 0.0) {
//					throw new Exception("PAR-FFRT Variance issue");
//				}
				// If the ratio test passes, return the index
				return new Object[] { ilsResult, stats, i };
			}

		}
		return null; // Should not occur if threshold logic is correct
	}

	/**
	 * Class to hold the results of the PAR estimation.
	 */
	public static class PARResult_FFRT {
		private SimpleMatrix aPAR;
		private int nFixed;
		private double SR_PAR;
		private Object[] stats;

		public PARResult_FFRT(SimpleMatrix aPAR, int nFixed, double SR_PAR, Object[] stats) {
			this.aPAR = aPAR;
			this.nFixed = nFixed;
			this.SR_PAR = SR_PAR;
			this.stats = stats;
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

		public Object[] getStats() {
			return stats;
		}
	}

}
