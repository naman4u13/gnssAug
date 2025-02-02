package com.gnssAug.helper.lambdaNew.Estimators;

import org.ejml.simple.SimpleMatrix;

import java.util.ArrayList;
import java.util.Arrays;

import org.ejml.data.DMatrix;
import org.ejml.data.DMatrixRMaj;
import org.ejml.simple.SimpleMatrix;
import org.ejml.simple.SimpleEVD;
import org.ejml.simple.SimpleSVD;

import com.gnssAug.helper.lambdaNew.GammaIncompleteInverse;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorILS.ILSResult;
import com.gnssAug.utility.Matrix;

public class EstimatorBIE {

	private ArrayList<SimpleMatrix> allIntCandidates;
	private ArrayList<Double> expTlist;

	public EstimatorBIE() {
		allIntCandidates = new ArrayList<SimpleMatrix>();
		expTlist = new ArrayList<Double>();
	}

	/**
	 * Best Integer Equivariant (BIE) estimator.
	 * 
	 * Computes the real-valued solution given by a Best Integer Equivariant (BIE)
	 * estimator assuming an underlying multivariate normal distribution. An initial
	 * search ellipsoid radius is provided to define an approximation to the
	 * infinite summation.
	 * 
	 * @param aHat Ambiguity float vector (column)
	 * @param LMat Decomposition L matrix (lower unitriangular)
	 * @param dVec Decomposition D matrix (vector of diagonal components)
	 * @param chi2 Radius of the initial search ellipsoid
	 * @param nMax Set maximum number of integer candidates
	 * @return EstimatorBIEResult containing aBIE and nInt
	 */
	public EstimatorBIEResult estimatorBIE(SimpleMatrix aHat, SimpleMatrix LMat, double[] dVec, Double chi2,
			Integer nMax,SimpleMatrix qHat) {
		// Problem dimensionality
		int nn = dVec.length;
		
		// Check number of input arguments and set defaults if necessary
		double chi2Value;
		double nMaxValue;
		if (chi2 == null) {
			chi2Value = 2 * GammaIncompleteInverse.gammaincinv(1 - 1e-6, nn / 2.0); // Finite search domain for BIE
			nMaxValue = Double.POSITIVE_INFINITY; // No limit to the number of integers used
		} else if (nMax == null) {
			chi2Value = chi2;
			nMaxValue = Double.POSITIVE_INFINITY; // No limit to the number of integers used
		} else {
			chi2Value = chi2;
			nMaxValue = nMax;
		}

		// Initialize main variables
		SimpleMatrix z_NUM = new SimpleMatrix(nn, 1); // Numerator vector for aBIE
		double z_DEN = 0.0; // Denominator scalar for aBIE
		int nInt = 0; // Number of integer vectors used

		// Auxiliary variable to store current integer vector used
		SimpleMatrix zVect = new SimpleMatrix(nn, 1);

		// Compute min-max values of the n-th component, i.e., zVect(nn)
		double aHatLast = aHat.get(nn - 1, 0);
		double dVecLast = dVec[nn - 1];
		int zMin =  (int) Math.ceil(aHatLast - Math.sqrt(dVecLast * chi2Value));
		int zMax =  (int) Math.floor(aHatLast + Math.sqrt(dVecLast * chi2Value));

		// Iterate over possible (if any) components "z_now" at the nth level
		for (int zNow = zMin; zNow <= zMax; zNow++) {
			zVect.set(nn - 1, 0, zNow); // Current n-th integer component

			// Fractional part of float nth component
			double zRest = aHatLast - zNow;

			// Compute argument of the exponential function (Gaussian distribution)
			double tValue = Math.pow(zRest, 2) / dVecLast;
			double expT = Math.exp(-0.5 * tValue);

			// Recursive nested function for lower levels kk = nn-1, nn-2, ..., 1.
			if (nn > 1) {
				// Update quantities for successive lower level nn-1
				double chi2New = chi2Value - tValue;
				SimpleMatrix zCond = aHat.extractMatrix(0, nn - 1, 0, 1)
						.minus(LMat.extractMatrix(nn - 1, nn, 0, nn - 1).transpose().scale(zRest));

				// Call function at lower level
				EstimatorBIENestedResult nestedResult = estimatorBIE_nested(zCond,
						LMat.extractMatrix(0, nn - 1, 0, nn - 1), Arrays.copyOfRange(dVec, 0, nn - 1), nn - 1, chi2New,
						nInt, zVect, nMaxValue,aHat,qHat);

				// Update main variables at current last level kk, i.e., kk = nn
				z_NUM = z_NUM.plus(nestedResult.zBIETemp.scale(expT));
				z_DEN += expT * nestedResult.zPdfTemp;

				nInt = nestedResult.nInt;
			} else {
				// Update main variables of this scalar problem, i.e., nn = 1
				z_NUM = z_NUM.plus(zVect.scale(expT));
				z_DEN += expT;
				nInt += 1;
				allIntCandidates.add(new SimpleMatrix(zVect));
				SimpleMatrix res = zVect.minus(aHat);
				double EXPT = Math.exp(-0.5*res.transpose().mult(qHat.invert()).mult(res).get(0));
				expTlist.add(EXPT);
			}

			// Check if exceeding the maximum number of integer candidates
			if (nInt > nMaxValue || nInt == 0) {
				nInt = 0; // An ILS-based BIE solution is computed
				break;
			}
		}

		SimpleMatrix aBIE;
		if (nInt == 0) {
			// EMPTY SET: force computation of BIE by using nearby integer candidates

			// Use a minimum number of integers based on neighbour pull-in regions
			int nIntegers = (int) (1 + (2 * (Math.pow(2, nn) - 1))); // More efficient in high dimensions

			// Call ILS estimator to find closest "nIntegers" candidates
			ILSResult ilsResult = new EstimatorILS().estimatorILS(aHat, LMat, dVec, nIntegers);
			for(int i=0;i<nIntegers;i++)
			{
				SimpleMatrix zVecttemp = new SimpleMatrix(ilsResult.getAFix().extractVector(false, i));
				allIntCandidates.add(zVecttemp);
				
			}
			// Define BIE weights (for Gaussian distribution)
			SimpleMatrix wBIE = new SimpleMatrix(ilsResult.getSqNorm().length, 1);
			for (int i = 0; i < ilsResult.getSqNorm().length; i++) {
				double sumExp = 0.0;
				for (int j = 0; j < ilsResult.getSqNorm().length; j++) {
					sumExp += Math.exp(-0.5 * (ilsResult.getSqNorm()[i] - ilsResult.getSqNorm()[j]));
				}
				wBIE.set(i, 0, 1.0 / sumExp);
				expTlist.add(1/sumExp);
			}

			// Return BIE solution
			aBIE = ilsResult.getAFix().mult(wBIE);
		} else {
			// Return BIE solution
			aBIE = z_NUM.divide(z_DEN);
		}

		EstimatorBIEResult result = new EstimatorBIEResult(aBIE, nInt);

		return result;
	}

	public SimpleMatrix getAllIntCandidates() {
		SimpleMatrix _allIntCandidates = new SimpleMatrix(allIntCandidates.get(0).numRows(), allIntCandidates.size());
		for(int i=0;i<allIntCandidates.size();i++)
		{
			_allIntCandidates.setColumn(i, 0, Matrix.matrix2ArrayVec(allIntCandidates.get(i)));
		}
		return _allIntCandidates;
	}

	public SimpleMatrix getExpTlist() {
		double[] expArray = expTlist.stream().mapToDouble(i->i).toArray();
		SimpleMatrix expMatrix = new SimpleMatrix(expArray.length,1,true,expArray);
		return expMatrix;
	}

	/**
	 * Nested function for BIE estimator.
	 * 
	 * @param aHat  Ambiguity float vector at current level
	 * @param LMat  Decomposition L matrix at current level
	 * @param dVec  Decomposition D matrix at current level
	 * @param kk    Current level
	 * @param chi2  Updated Chi2 value
	 * @param nInt  Current number of integer candidates
	 * @param zVect Current integer vector
	 * @param nMax  Maximum number of integer candidates
	 * @return EstimatorBIENestedResult containing zBIETemp, zPdfTemp, and nInt
	 */
	private EstimatorBIENestedResult estimatorBIE_nested(SimpleMatrix aHat, SimpleMatrix LMat, double[] dVec,
			int kk, double chi2, int nInt, SimpleMatrix zVect, double nMax,SimpleMatrix AHAT,SimpleMatrix QHAT) {
		// Initialize main variables
		SimpleMatrix zBIE = new SimpleMatrix(zVect.numRows(), zVect.numCols());
		double zPdf = 0.0;

		// Define integers' range at current k-th level
		double aHatK = aHat.get(kk - 1, 0);
		double dVecK = dVec[kk - 1];
		int zMin =  (int) Math.ceil(aHatK - Math.sqrt(dVecK * chi2));
		int zMax =  (int) Math.floor(aHatK + Math.sqrt(dVecK * chi2));

		// Iterate over possible (if any) components "z_now" at this k-th level
		for (int zNow = zMin; zNow <= zMax; zNow++) {
			zVect.set(kk - 1, 0, zNow); // Current k-th integer component

			// Fractional part of float k-th component
			double zRest = aHatK - zNow;

			// Compute argument of the exponential function (Gaussian distribution)
			double tValue = Math.pow(zRest, 2) / dVecK;
			double expT = Math.exp(-0.5 * tValue);

			// Recursive nested function for lower levels kk-1, kk-2, ..., 1.
			if (kk > 1) {
				// Update quantities for next lower level kk-1
				double chi2New = chi2 - tValue;
				SimpleMatrix zCond = aHat.extractMatrix(0, kk - 1, 0, 1)
						.minus(LMat.extractMatrix(kk - 1, kk, 0, kk - 1).transpose().scale(zRest));

				// Call function at lower level kk-1
				EstimatorBIENestedResult nestedResult = estimatorBIE_nested(zCond,
						LMat.extractMatrix(0, kk - 1, 0, kk - 1), Arrays.copyOfRange(dVec, 0, kk - 1), kk - 1, chi2New,
						nInt, zVect, nMax,AHAT,QHAT);

				// Update variables at current k-th level
				zBIE = zBIE.plus(nestedResult.zBIETemp.scale(expT));
				zPdf += expT * nestedResult.zPdfTemp;

				nInt = nestedResult.nInt;
			} else {
				// Update the main variables once reaching level kk = 1
				zBIE = zBIE.plus(zVect.scale(expT));
				zPdf += expT;
				nInt += 1;
				allIntCandidates.add(new SimpleMatrix(zVect));
				SimpleMatrix res = zVect.minus(AHAT);
				double EXPT = Math.exp(-0.5*res.transpose().mult(QHAT.invert()).mult(res).get(0));
				expTlist.add(EXPT);
			}

			// Check if exceeding the maximum number of integer candidates
			if (nInt > nMax || nInt == 0) {
				nInt = 0; // An ILS-based BIE solution is computed
				break;
			}
		}

		EstimatorBIENestedResult result = new EstimatorBIENestedResult();
		result.zBIETemp = zBIE;
		result.zPdfTemp = zPdf;
		result.nInt = nInt;
		return result;
	}

	/**
	 * Class to hold the results of the estimatorBIE function.
	 */
	public class EstimatorBIEResult {
		private SimpleMatrix aBIE; // Ambiguity fixed (real-valued) solution using BIE
		private int nInt; // Number of integer candidates included in the computation

		public EstimatorBIEResult(SimpleMatrix aBIE, int nInt) {
			super();
			this.aBIE = aBIE;
			this.nInt = nInt;
		}

		public SimpleMatrix getaBIE() {
			return aBIE;
		}

		public int getnInt() {
			return nInt;
		}

	}

	/**
	 * Class to hold the results of the nested estimatorBIE_nested function.
	 */
	private class EstimatorBIENestedResult {
		public SimpleMatrix zBIETemp; // Temporary zBIE
		public double zPdfTemp; // Temporary zPDF
		public int nInt; // Number of integer candidates

	}

}
