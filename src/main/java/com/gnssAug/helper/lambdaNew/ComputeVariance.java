package com.gnssAug.helper.lambdaNew;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.CholeskyDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.decomposition.chol.CholeskyDecompositionLDL_DDRM;
import org.ejml.dense.row.factory.DecompositionFactory_DDRM;
import org.ejml.interfaces.decomposition.CholeskyDecomposition_F64;
import org.ejml.simple.SimpleMatrix;

import com.gnssAug.helper.lambdaNew.DecomposeLtDL.DecompositionResult;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorIA_FFRT;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorIA_FFRT.IAFFRTResult;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorIB;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorILS;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorILS.ILSResult;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorIR;
import com.gnssAug.utility.Matrix;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.stream.IntStream;

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
			Integer nSamples) {
		// Problem dimensionality
		int nn = QMat.numCols();

		HashSet<String> allCandidatesKey = new HashSet<String>();
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
				ILSResult ilsResult = EstimatorILS.estimatorILS(aHat, LMat, dVec, 1);
				SimpleMatrix aFix = ilsResult.getAFix();
				for (SimpleMatrix res : ilsResult.getAllCandidates()) {
					String key = serializeMatrix(res, 1e-9);
					allCandidatesKey.add(key);
				}
				for (int row = 0; row < aFix.numRows(); row++) {
					aFixAll.set(row, ii, aFix.get(row));
				}
			}
			break;

		case 2: // Use IA-FFRT (ILS w/ Fixed Failure-rate Ratio Test)
			for (int ii = 0; ii < nSamples; ii++) {
				SimpleMatrix aHat = aHatAll.extractVector(false, ii);
				IAFFRTResult IAFFRTresult = EstimatorIA_FFRT.estimatorIA_FFRT(aHat, LMat, dVec, maxFR, null);
				SimpleMatrix aFix = IAFFRTresult.getaFix();
				for (SimpleMatrix res : IAFFRTresult.getAllCandidates()) {
					String key = serializeMatrix(res, 1e-9);
					allCandidatesKey.add(key);

				}
				for (int row = 0; row < aFix.numRows(); row++) {
					aFixAll.set(row, ii, aFix.get(row));
				}
			}
			break;
		case 3: // Use ILS
			try {
				Object[] parallelThreadILSsoln = parallelThreadILS(cholQMat, nn, nSamples, LMat, dVec, rand);
				aFixAll = (SimpleMatrix) parallelThreadILSsoln[0];
				allCandidatesKey = (HashSet<String>) parallelThreadILSsoln[1];
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			break;

		default:
			throw new IllegalArgumentException("ATTENTION: the estimator selected is not available! Use 1-7.");
		}
		Object[] varCalRes = OptimizedVarCalc.calculateVariance(allCandidatesKey, aFixAll, nSamples, 1e-9);
		SimpleMatrix variance = (SimpleMatrix) varCalRes[0];
		double approxSR = (double) varCalRes[1];
		double approxFR = (double) varCalRes[2];
		System.out.println("1 - Approximate Success Rate : " + (1 - approxSR));
		System.out.println("Approximate Failure Rate : " + approxFR);
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

	private static Object[] parallelThreadILS(SimpleMatrix cholQMat, int nn, int nSamples, SimpleMatrix LMat,
			double[] dVec, Random rand) throws InterruptedException {

		// Generate aHatAll matrix
		SimpleMatrix aHatAll = cholQMat.transpose().mult(generateRandn(nn, nSamples, rand));

		// Thread-safe map for storing aFix results (using column index as key)
		ConcurrentMap<Integer, SimpleMatrix> aFixMap = new ConcurrentHashMap<>();

		// Thread-safe set implementation using ConcurrentHashMap for storing all
		// candidate keys
		ConcurrentHashMap<String, Boolean> allCandidatesKeySet = new ConcurrentHashMap<>();

		// Parallel processing using IntStream
		IntStream.range(0, nSamples).parallel().forEach(ii -> {
			// Extract column vector aHat for the current sample
			SimpleMatrix aHat = aHatAll.extractVector(false, ii);

			// Perform the ILS estimation
			ILSResult ilsResult = EstimatorILS.estimatorILS(aHat, LMat, dVec, 1);

			// Get the fixed solution (aFix)
			SimpleMatrix aFix = ilsResult.getAFix();

			// Store all candidate keys in a thread-safe set
			for (SimpleMatrix res : ilsResult.getAllCandidates()) {
				String key = serializeMatrix(res, 1e-9);
				allCandidatesKeySet.put(key, Boolean.TRUE); // Simulating a set
			}

			// Store the fixed solution (aFix) in the thread-safe map
			aFixMap.put(ii, aFix);
		});

		// Merge results back into aFixAll
		SimpleMatrix aFixAll = new SimpleMatrix(nn, nSamples);
		aFixMap.forEach((colIndex, aFix) -> {
			for (int row = 0; row < aFix.numRows(); row++) {
				aFixAll.set(row, colIndex, aFix.get(row));
			}
		});

		HashSet<String> allCandidatesKey = new HashSet<String>(allCandidatesKeySet.keySet());
		return new Object[] { aFixAll, allCandidatesKey };
	}

	/**
	 * Serializes a SimpleMatrix to a key with rounding to avoid numerical precision
	 * issues.
	 *
	 * @param matrix    The matrix to serialize (assumes column vector).
	 * @param tolerance Numerical tolerance for rounding.
	 * @return A serialized string key representing the matrix.
	 */
	private static String serializeMatrix(SimpleMatrix matrix, double tolerance) {
		StringBuilder keyBuilder = new StringBuilder();
		for (int i = 0; i < matrix.numRows(); i++) {
			// Round each element to the specified tolerance
			double roundedValue = Math.round(matrix.get(i, 0) / tolerance) * tolerance;
			keyBuilder.append(roundedValue).append(",");
		}
		return keyBuilder.toString();
	}

}
