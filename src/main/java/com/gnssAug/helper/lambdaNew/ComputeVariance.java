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
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

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
        RealMatrix _cholQMat  = new CholeskyDecomposition(_QMat).getL();
		
		SimpleMatrix cholQMat = new SimpleMatrix(_cholQMat.getData());
		long startTime = System.nanoTime();
		// Now cholQMat contains the lower triangular matrix from the Cholesky
		// decomposition
		SimpleMatrix aHatAll = cholQMat.transpose().mult(generateRandn(nn, nSamples, rand));
		double timeCPU = (System.nanoTime() - startTime) / 1e9;
		System.out.println("Random Function time: "+timeCPU);
		// Initialize all the ambiguity fixed vectors
		SimpleMatrix aFixAll = new SimpleMatrix(nn, nSamples);
		for (int i = 0; i < aFixAll.getNumElements(); i++) {
			aFixAll.set(i, Double.NaN);
		}

		// ESTIMATORS: numerical simulations for computing the success rate
		startTime = System.nanoTime();

		switch (estimator) {
		case 1: // Use ILS
			for (int ii = 0; ii < nSamples; ii++) {
				SimpleMatrix aHat = aHatAll.extractVector(false, ii);
				ILSResult ilsResult = EstimatorILS.estimatorILS(aHat, LMat, dVec, 1);
				SimpleMatrix aFix = ilsResult.getAFix();
				for(SimpleMatrix res:ilsResult.getAllCandidates())
				{
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
				for(SimpleMatrix res:IAFFRTresult.getAllCandidates())
				{
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
				aFixAll = computeILSResults( cholQMat,  nn,  nSamples,  LMat,  dVec,  rand);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			break;

		default:
			throw new IllegalArgumentException("ATTENTION: the estimator selected is not available! Use 1-7.");
		}

		timeCPU = (System.nanoTime() - startTime) / 1e9;
		System.out.println("Estimation time : " + timeCPU);

		startTime = System.nanoTime();
		SimpleMatrix variance = new SimpleMatrix(nn, nn);
		 // Step 2: Process each candidate and compute contributions
        for (String candidateKey : allCandidatesKey) {
            
            int count = 0; 
            SimpleMatrix candidate = deserializeMatrix(candidateKey, nn);

			for (int jj = 0; jj < nSamples; jj++) {
				if (Matrix.areMatricesEqual(candidate, aFixAll.extractVector(false, jj), 1e-9)) {
					count++;
				}
			}
			variance = variance.plus(candidate.mult(candidate.transpose()).scale((count * 1.0) / nSamples));

		}
//		 SimpleMatrix variance =
//		 OptimizedVarianceCalculator.calculateVariance(allCandidates, aFixAll,
//		 nSamples);
//		 SimpleMatrix variance =
//				 VarianceCalculator.calculateVariance(allCandidatesKey, aFixAll,
//		 nSamples,1e-9);
		
		timeCPU = (System.nanoTime() - startTime) / 1e9;
		System.out.println("Computing Variance time : " + timeCPU);
		return new VarianceResult(variance, timeCPU);
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
    
    public static SimpleMatrix computeILSResults(SimpleMatrix cholQMat, int nn, int nSamples, SimpleMatrix LMat, double[] dVec, Random rand) throws InterruptedException {
        // Step 1: Generate all samples using Cholesky decomposition
        SimpleMatrix aHatAll = cholQMat.transpose().mult(generateRandn(nn, nSamples, rand));

        // Step 2: Preallocate result matrix for fixed solutions
        SimpleMatrix aFixAll = new SimpleMatrix(nn, nSamples);

        // Step 3: Use ExecutorService for parallel processing
        int numThreads = Runtime.getRuntime().availableProcessors(); // Number of CPU cores
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);

        // Step 4: Process each sample in parallel
        List<Future<Void>> futures = new ArrayList<>();
        for (int threadId = 0; threadId < numThreads; threadId++) {
            int start = (threadId * nSamples) / numThreads; // Start index for this thread
            int end = ((threadId + 1) * nSamples) / numThreads; // End index for this thread

            Future<Void> future = executor.submit(() -> {
                for (int ii = start; ii < end; ii++) {
                    // Extract aHat for the current sample
                    SimpleMatrix aHat = aHatAll.extractVector(false, ii);

                    // Perform ILS estimation
                    ILSResult ilsResult = EstimatorILS.estimatorILS(aHat, LMat, dVec, 1);
                    SimpleMatrix aFix = ilsResult.getAFix();

                    // Store the fixed results in the shared matrix
                    synchronized (aFixAll) { // Synchronize access to shared resource
                        for (int row = 0; row < aFix.numRows(); row++) {
                            aFixAll.set(row, ii, aFix.get(row));
                        }
                    }
                }
                return null;
            });

            futures.add(future);
        }

        // Step 5: Wait for all threads to complete
        for (Future<Void> future : futures) {
            try {
				future.get();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (ExecutionException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} // This ensures the computation completes
        }

        // Step 6: Shutdown the executor
        executor.shutdown();

        return aFixAll;
    }
    
    /**
     * Serializes a SimpleMatrix to a key with rounding to avoid numerical precision issues.
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
    private static SimpleMatrix deserializeMatrix(String serializedMatrix, int numRows) {
        // Split the serialized string into individual values using a comma as the delimiter
        String[] values = serializedMatrix.split(",");

        // Ensure the number of rows matches the number of values in the serialized matrix
        if (values.length != numRows) {
            throw new IllegalArgumentException("The number of rows does not match the serialized matrix length.");
        }

        // Create a new SimpleMatrix with the specified number of rows and 1 column
        SimpleMatrix matrix = new SimpleMatrix(numRows, 1);

        // Parse each value and set it into the matrix
        for (int i = 0; i < numRows; i++) {
            double value = Double.parseDouble(values[i]); // Convert the string value back to a double
            matrix.set(i, 0, value); // Set the value in the matrix
        }

        return matrix; // Return the reconstructed SimpleMatrix
    }
}
