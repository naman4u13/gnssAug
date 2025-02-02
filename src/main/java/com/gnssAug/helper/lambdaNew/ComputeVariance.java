package com.gnssAug.helper.lambdaNew;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.CholeskyDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.ejml.simple.SimpleMatrix;
import com.gnssAug.helper.lambdaNew.DecomposeLtDL.DecompositionResult;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorBIE;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorBIE.EstimatorBIEResult;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorIA_FFRT;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorIA_FFRT.IAFFRTResult;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorILS;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorILS.ILSResult;
import com.gnssAug.utility.Matrix;
import java.util.HashMap;
import java.util.HashSet;

import java.util.Random;

public class ComputeVariance {

	/**
	 * Computes the success rate based on numerical simulations.
	 *
	 * @param QMat Variance-Covariance matrix of the ambiguities
	 * @return Result object containing Success Rate (SR), Failure Rate (FR), and
	 *         Computational Time (timeCPU)
	 * @throws IllegalArgumentException if not enough inputs are provided or if an
	 *                                  invalid estimator is selected
	 */
	public static VarianceResult computeVariance(SimpleMatrix QMat, int estimator, int decorr, Double maxFR,
			Integer nSamples, Double muRatio) {
		// Problem dimensionality
		int nn = QMat.numCols();
		SimpleMatrix LMat;
		double[] dVec;

		// Z-transformation (optional) for decorrelating the ambiguities
		if (decorr == 1) {
			DecorrelateVCResult decorResult = DecorrelateVC.decorrelateVC(QMat.copy(), null);
			QMat = decorResult.getQzHat();
			LMat = decorResult.getLzMat();
			dVec = decorResult.getDzVec();

		} else {
			// Retrieve LtDL-decomposition of the original ambiguity vc-matrix
			DecompositionResult ltldlResult = DecomposeLtDL.decomposeLtDL(QMat.copy());
			LMat = ltldlResult.getLMat();
			dVec = ltldlResult.getDVec();
		}
		// Decorrelate the ambiguity vc-matrix and return a LtDL-decomposition

		// Generation of numerical samples
		if (nSamples == 0 || nSamples == null) {
			// Get an approximative value of success rate based on IB formulation
			double P0 = ComputeSR_IBexact.computeSR_IBexact(dVec).getSR();

			// Compute the number of samples to be used
			nSamples = Utilities.computeNumSamples(P0);
		}
		// ESTIMATORS: numerical simulations for computing the success rate

		// Initialize random number generator
		Random rand = new Random();
		// Check if the matrix is symmetric
		RealMatrix _QMat = new Array2DRowRealMatrix(Matrix.matrix2Array(QMat));
		if (!isSymmetric(_QMat, 1e-20)) {

			_QMat = makeSymmetric(_QMat);
		}
		RealMatrix _cholQMat = new CholeskyDecomposition(_QMat).getL();

		SimpleMatrix cholQMat = new SimpleMatrix(_cholQMat.getData());

		// Now cholQMat contains the lower triangular matrix from the Cholesky
		// decomposition
		SimpleMatrix aHatAll = cholQMat.transpose().mult(generateRandn(nn, nSamples, rand));

		// Initialize all the ambiguity fixed vectors
		SimpleMatrix aFixAll = new SimpleMatrix(nn, nSamples);
		for (int i = 0; i < aFixAll.getNumElements(); i++) {
			aFixAll.set(i, Double.NaN);
		}

		// ESTIMATORS: numerical simulations for computing the success rate

		switch (estimator) {
		case 1: // Use ILS
			for (int ii = 0; ii < nSamples; ii++) {
				SimpleMatrix aHat = aHatAll.extractVector(false, ii);
				ILSResult ilsResult = new EstimatorILS().estimatorILS(aHat, LMat, dVec, 1);
				SimpleMatrix aFix = ilsResult.getAFix();
				aFixAll.insertIntoThis(0, ii, aFix);

			}
			break;

		case 2: // Use IA-FFRT (ILS w/ Fixed Failure-rate Ratio Test)
			for (int ii = 0; ii < nSamples; ii++) {
				SimpleMatrix aHat = aHatAll.extractVector(false, ii);
				IAFFRTResult IAFFRTresult = new EstimatorIA_FFRT().estimatorIA_FFRT(aHat, LMat, dVec, maxFR, muRatio);
				SimpleMatrix aFix = IAFFRTresult.getaFix();
				if (IAFFRTresult.getnFixed() != 0) {
					aFixAll.insertIntoThis(0, ii, aFix);

				}

			}
			break;
			
		case 3: // Use IA-FFRT (ILS w/ Fixed Failure-rate Ratio Test)
			for (int ii = 0; ii < nSamples; ii++) {
				SimpleMatrix aHat = aHatAll.extractVector(false, ii);
				EstimatorBIEResult bieResult = new EstimatorBIE().estimatorBIE(aHat, LMat, dVec, null, null, QMat);
				SimpleMatrix aFix = bieResult.getaBIE();
				aFixAll.insertIntoThis(0, ii, aFix);

			}
			break;

		default:
			throw new IllegalArgumentException("ATTENTION: the estimator selected is not available! Use 1-2.");
		}
		Object[] varCalRes = OptimizedVarCalc2.calculateVariance(aFixAll, nSamples);
		SimpleMatrix variance = (SimpleMatrix) varCalRes[0];
		double approxSR = (double) varCalRes[1];
		double approxFR = (double) varCalRes[2];
		System.out.println("Approximate Success Rate : " + (approxSR*100));
		System.out.println("1 - Approximate Success Rate : " + (1 - approxSR)*100);
		System.out.println("Approximate Failure Rate : " + approxFR*100);
		return new VarianceResult(variance, 0);
	}

	// Helper method to generate a matrix of random Gaussian numbers
	private static SimpleMatrix generateRandn(int rows, int cols, Random rand) {
		double[] data = new double[rows * cols];
		for (int i = 0; i < data.length; i++) {
			data[i] = rand.nextGaussian();
		}
		return new SimpleMatrix(rows, cols, true, data);
	}

	// Result class to hold output values
	public static class VarianceResult {
		private final SimpleMatrix variance;
		private final double timeCPU;

		public VarianceResult(SimpleMatrix variance, double timeCPU) {
			this.variance = variance;
			this.timeCPU = timeCPU;
		}

		public SimpleMatrix getVariance() {
			return variance;
		}

		public double getTimeCPU() {
			return timeCPU;
		}

	}

	private static boolean isSymmetric(RealMatrix matrix, double tolerance) {
		for (int i = 0; i < matrix.getRowDimension(); i++) {
			for (int j = 0; j < matrix.getColumnDimension(); j++) {
				if (Math.abs(matrix.getEntry(i, j) - matrix.getEntry(j, i)) > tolerance) {
					return false;
				}
			}
		}
		return true;
	}

	private static RealMatrix makeSymmetric(RealMatrix matrix) {
		int n = matrix.getRowDimension();
		RealMatrix symmetricMatrix = matrix.copy();

		for (int i = 0; i < n; i++) {
			for (int j = i + 1; j < n; j++) {
				double avg = (matrix.getEntry(i, j) + matrix.getEntry(j, i)) / 2.0;
				symmetricMatrix.setEntry(i, j, avg);
				symmetricMatrix.setEntry(j, i, avg);
			}
		}
		return symmetricMatrix;
	}

}
