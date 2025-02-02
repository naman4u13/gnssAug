package com.gnssAug.helper.lambdaNew;

import org.ejml.simple.SimpleMatrix;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

public class OptimizedVarCalc2 {

	/**
	 * Calculates the variance matrix based on the provided candidates and samples.
	 *
	 * @param allCandidatesKey List of candidate vectors (column vectors as
	 *                         SimpleMatrix objects).
	 * @param aFixAll          Matrix containing all fixed vectors (column-wise
	 *                         SimpleMatrix).
	 * @param nSamples         Number of samples (columns in aFixAll).
	 * @param tolerance        Numerical tolerance for comparing matrices.
	 * @return Variance matrix (SimpleMatrix).
	 */
	public static Object[] calculateVariance(SimpleMatrix aFixAll, int nSamples) {
		int nn = aFixAll.numRows();
		SimpleMatrix variance = new SimpleMatrix(nn, nn);

		// Input validation to avoid edge cases
		if (aFixAll == null || aFixAll.numCols() == 0 || nSamples <= 0) {
			throw new IllegalArgumentException("Invalid inputs: Samples cannot be empty, and nSamples must be > 0.");
		}
		double failRate = 0;
		double successRate = 0;

		for (int jj = 0; jj < nSamples; jj++) {
			SimpleMatrix sample = aFixAll.extractVector(false, jj); // Extract column vector
            double sum =sample.elementSum();
            if(Double.isNaN(sum)) {
            	continue;
            }
			variance = variance.plus(sample.mult(sample.transpose()).scale((1.0) / nSamples));
			if (areAllElementsIntegers(sample)) {
				if (sum == 0.0) {
					successRate++;
				} else {
					failRate++;
				}
			}

		}
		// Step 3: Enforce symmetry on the variance matrix (optional but useful for
		// numerical stability)
		variance = enforceSymmetry(variance);
		failRate = failRate / nSamples;
		successRate = successRate/nSamples;
		return new Object[] { variance, successRate, failRate };
	}

	/**
	 * Enforces symmetry on a matrix by averaging it with its transpose.
	 *
	 * @param matrix The input matrix (SimpleMatrix).
	 * @return A symmetric version of the matrix.
	 */
	private static SimpleMatrix enforceSymmetry(SimpleMatrix matrix) {
		return matrix.plus(matrix.transpose()).scale(0.5);
	}

	public static boolean areAllElementsIntegers(SimpleMatrix matrix) {
		int rows = matrix.numRows();
		int cols = matrix.numCols();

		// Iterate over all elements and check if they are integers
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				double element = matrix.get(i, j);
				// Check if the element is an integer (i.e., it's equal to its integer cast)
				if (element != (int) element) { // If element is not equal to its integer part
					return false;
				}
			}
		}
		return true; // All elements are integers
	}
}