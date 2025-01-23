package com.gnssAug.helper.lambdaNew;

import java.util.Arrays;

import org.apache.commons.math3.linear.*;

public class Estimators {

	/**
	 * Computes a fixed solution based on the Integer Rounding (IR) estimator.
	 *
	 * @param aHat The float ambiguity vector (RealVector).
	 * @return The fixed ambiguity vector (RealVector) obtained by rounding.
	 */
	public static RealVector estimatorIR(RealVector aHat) {
		// Create a new vector to store the fixed values
		RealVector aFix = aHat.copy();

		// Round each component of the input vector to the nearest integer
		for (int i = 0; i < aFix.getDimension(); i++) {
			aFix.setEntry(i, Math.round(aFix.getEntry(i)));
		}

		return aFix;
	}

	/**
	 * Computes a fixed solution based on the Integer Bootstrapping (IB) estimator.
	 *
	 * @param aHat The ambiguity float vector as a RealVector.
	 * @param LMat The LtDL-decomposition matrix L as a RealMatrix.
	 * @return A 2-element array containing: - [0]: The fixed ambiguity vector
	 *         (aFix) as a RealVector. - [1]: The conditioned ambiguity vector
	 *         (aCond) as a RealVector.
	 */
	public static RealVector estimatorIB(RealVector aHat, RealMatrix LMat) {
		int n = aHat.getDimension();

		// Initialize main vectors
		RealVector aFix = new ArrayRealVector(n); // Integer-fixed ambiguity vector
		RealVector aCond = aHat.copy(); // Float (conditioned) ambiguity vector
		RealVector sumCond = new ArrayRealVector(n); // Auxiliary vector for conditioned ambiguities

		// Fix the last component of the ambiguity
		aFix.setEntry(n - 1, Math.round(aCond.getEntry(n - 1)));

		// Iterative cycle, conditioning (last-to-first)
		for (int i = n - 2; i >= 0; i--) {
			// Compute the i-th ambiguity conditioned on the previous ones (i+1 to n)
			for (int j = i + 1; j < n; j++) {
				double correction = LMat.getEntry(j, i) * (aCond.getEntry(j) - aFix.getEntry(j));
				sumCond.addToEntry(i, -correction);
			}

			// Conditioned ambiguity for the i-th component
			aCond.setEntry(i, aHat.getEntry(i) + sumCond.getEntry(i));

			// Sequentially round the conditioned i-th component
			aFix.setEntry(i, Math.round(aCond.getEntry(i)));
		}

		// Return both aFix and aCond
		return aFix;
	}

	/**
	 * Integer Least-Squares (ILS) estimator by search-and-shrink.
	 *
	 * @param aHat   Ambiguity float vector (RealVector).
	 * @param LMat   LtDL-decomposition matrix L (lower unitriangular, RealMatrix).
	 * @param dVec   LtDL-decomposition matrix D (diagonal elements, double array).
	 * @param nCands Number of best integer solutions (default is 1).
	 * @return A 2-element array containing: - [0]: Fixed ambiguity vector
	 *         (RealVector). - [1]: Squared norm of the solution as a double.
	 */
	public static Object[] estimatorILS(RealVector aHat, RealMatrix LMat, double[] dVec, int nCands) {
		int n = aHat.getDimension();

		if (nCands <= 0) {
			throw new IllegalArgumentException("Number of candidates (nCands) must be greater than zero.");
		}

		// Initialization of output variables
		RealVector[] aFix = new RealVector[nCands];
		double[] sqNorm = new double[nCands];
		Arrays.fill(sqNorm, Double.POSITIVE_INFINITY); // Initialize to infinity

		// Start from an ellipsoid with infinite radius
		double maxChi2 = Double.POSITIVE_INFINITY;

		// Variables for the iterative process
		RealVector aCond = aHat.copy();
		RealVector zCond = new ArrayRealVector(n);
		RealVector left = new ArrayRealVector(n);
		RealVector step = new ArrayRealVector(n);

		// Initialization at the last level
		aCond.setEntry(n - 1, aHat.getEntry(n - 1));
		zCond.setEntry(n - 1, Math.round(aCond.getEntry(n - 1)));
		left.setEntry(n - 1, aCond.getEntry(n - 1) - zCond.getEntry(n - 1));
		step.setEntry(n - 1, Math.signum(left.getEntry(n - 1)));

		// Ensure step is positive to avoid stalls
		if (step.getEntry(n - 1) == 0) {
			step.setEntry(n - 1, 1);
		}

		RealMatrix S = LMat.copy();
		double[] dist = new double[n + 1]; // Distance metric
		Arrays.fill(dist, 0);

		int k = n - 1; // Start from the last ambiguity component
		boolean endSearch = false;

		while (!endSearch) {
			double newDist = dist[k] + Math.pow(left.getEntry(k), 2) / dVec[k];

			while (newDist < maxChi2) {
				if (k > 0) {
					// Move down to the next level
					k--;
					dist[k] = newDist;

					// Update S matrix
					for (int j = k + 1; j < n; j++) {
						S.setEntry(j - 1, k, S.getEntry(j, k) - left.getEntry(j) * LMat.getEntry(j, k));
					}

					// Update conditioned ambiguity
					aCond.setEntry(k, aHat.getEntry(k) + S.getColumnVector(k).getSubVector(k, n - k).dotProduct(left));
					zCond.setEntry(k, Math.round(aCond.getEntry(k)));
					left.setEntry(k, aCond.getEntry(k) - zCond.getEntry(k));
					step.setEntry(k, Math.signum(left.getEntry(k)));

					// Ensure step is positive
					if (step.getEntry(k) == 0) {
						step.setEntry(k, 1);
					}
				} else {
					// Store the candidate found
					RealVector candidate = zCond.copy();
					double candidateSqNorm = newDist;

					if (sqNorm[nCands - 1] > candidateSqNorm) {
						sqNorm[nCands - 1] = candidateSqNorm;
						aFix[nCands - 1] = candidate;
						// Sort solutions by squared norm
						for (int i = nCands - 1; i > 0; i--) {
							if (sqNorm[i] < sqNorm[i - 1]) {
								double tempNorm = sqNorm[i];
								sqNorm[i] = sqNorm[i - 1];
								sqNorm[i - 1] = tempNorm;

								RealVector tempVec = aFix[i];
								aFix[i] = aFix[i - 1];
								aFix[i - 1] = tempVec;
							}
						}
					}

					// Try next valid integer
					k = 0;
					zCond.setEntry(k, zCond.getEntry(k) + step.getEntry(k));
					left.setEntry(k, aCond.getEntry(k) - zCond.getEntry(k));
					step.setEntry(k, -step.getEntry(k) - Math.signum(step.getEntry(k)));
				}
				newDist = dist[k] + Math.pow(left.getEntry(k), 2) / dVec[k];
			}

			// Move up or end search
			while (newDist >= maxChi2) {
				if (k == n - 1) {
					endSearch = true;
					break;
				}
				k++;
				zCond.setEntry(k, zCond.getEntry(k) + step.getEntry(k));
				left.setEntry(k, aCond.getEntry(k) - zCond.getEntry(k));
				step.setEntry(k, -step.getEntry(k) - Math.signum(step.getEntry(k)));
				newDist = dist[k] + Math.pow(left.getEntry(k), 2) / dVec[k];
			}
		}

		return new Object[] { aFix, sqNorm };
	}

	public static Object[] estimatorBIE(RealVector a_hat, RealMatrix L_mat, RealVector d_vec, double Chi2) {
		return estimatorBIE(a_hat, L_mat, d_vec, Chi2, Integer.MAX_VALUE);
	}

	public static Object[] estimatorBIE(RealVector a_hat, RealMatrix L_mat, RealVector d_vec, double Chi2, int N_max) {
		int nn = (int) d_vec.getDimension();

		// Default input arguments
		if (Chi2 <= 0) {
			Chi2 = 2 * GammaIncompleteInverse.gammaincinv(1 - 1e-6, nn / 2.0);
			N_max = Integer.MAX_VALUE;
		}

		// Initialize main variables
		RealVector z_NUM = new ArrayRealVector(nn); // Numerator for a_BIE
		double z_DEN = 0.0; // Denominator for a_BIE
		int N_int = 0; // Count of integer vectors used

		// Auxiliary variable to store current integer vector
		RealVector z_vect = new ArrayRealVector(nn);

		// Compute min-max values for the n-th component
		int z_min = (int) Math.ceil(a_hat.getEntry(nn - 1) - Math.sqrt(d_vec.getEntry(nn - 1) * Chi2));
		int z_max = (int) Math.floor(a_hat.getEntry(nn - 1) + Math.sqrt(d_vec.getEntry(nn - 1) * Chi2));

		// Iterate over possible integer components
		for (int z_now = z_min; z_now <= z_max; z_now++) {
			z_vect.setEntry(nn - 1, z_now); // Set the n-th component

			// Fractional part
			double z_rest = a_hat.getEntry(nn - 1) - z_vect.getEntry(nn - 1);

			// Compute Gaussian weight
			double t_value = Math.pow(z_rest, 2) / d_vec.getEntry(nn - 1);
			double exp_t = Math.exp(-0.5 * t_value);

			if (nn > 1) {
				// Update for lower levels (nn-1)
				double chi2_new = Chi2 - t_value;
				RealVector z_cond = a_hat.getSubVector(0, nn - 1)
						.subtract(L_mat.getColumnVector(nn - 1).getSubVector(0, nn - 1).mapMultiply(z_rest));

				// Recursive call to the nested estimatorBIE
				Object[] nestedResult = estimatorBIE_nested(z_cond, L_mat.getSubMatrix(0, nn - 2, 0, nn - 2),
						d_vec.getSubVector(0, nn - 1), nn - 1, chi2_new, N_int, z_vect, N_max);

				// Update main variables at the last level
				z_NUM = z_NUM.add(((RealVector) nestedResult[0]).mapMultiply(exp_t));
				z_DEN += exp_t * ((double) nestedResult[1]);
				N_int = (int) nestedResult[2];

			} else {
				// Scalar case: directly update main variables
				z_NUM = z_NUM.add(z_vect.mapMultiply(exp_t));
				z_DEN += exp_t;
				N_int++;
			}

			// Check if we exceed the max number of candidates
			if (N_int > N_max || N_int == 0) {
				N_int = 0; // Force fallback to ILS-based BIE
				break;
			}
		}

		// Handle empty candidate set
		if (N_int == 0) {
			int nIntegers = 1 + 2 * (nn * nn - 1); // Minimum number of candidates
			Object[] ilsOutput = estimatorILS(a_hat, L_mat, d_vec.toArray(), nIntegers);

			RealVector[] z_fixes = (RealVector[]) ilsOutput[0]; // Change to RealVector[]
			double[] sqnorm_fixes = (double[]) ilsOutput[1]; // Change to double[]

			RealVector w_BIE = new ArrayRealVector(z_fixes.length);
			for (int i = 0; i < w_BIE.getDimension(); i++) {
				w_BIE.setEntry(i, Math.exp(-0.5 * sqnorm_fixes[i]));
			}

			double sumWeights = Arrays.stream(w_BIE.toArray()).sum();
			w_BIE.mapDivideToSelf(sumWeights);

			// Compute the BIE solution as a weighted sum of z_fixes
			RealVector a_BIE = new ArrayRealVector(z_fixes[0].getDimension());
			for (int i = 0; i < z_fixes.length; i++) {
				a_BIE = a_BIE.add(z_fixes[i].mapMultiply(w_BIE.getEntry(i)));
			}

			return new Object[] { a_BIE, 0 };
		} else {
			// Return BIE solution
			return new Object[] { z_NUM.mapDivide(z_DEN), N_int };
		}
	}

	private static Object[] estimatorBIE_nested(RealVector a_hat, RealMatrix L_mat, RealVector d_vec, int kk,
			double Chi2, int N_int, RealVector z_vect, int N_max) {
		RealVector z_BIE = new ArrayRealVector(z_vect.getDimension());
		double z_PDF = 0;

		int z_min = (int) Math.ceil(a_hat.getEntry(kk - 1) - Math.sqrt(d_vec.getEntry(kk - 1) * Chi2));
		int z_max = (int) Math.floor(a_hat.getEntry(kk - 1) + Math.sqrt(d_vec.getEntry(kk - 1) * Chi2));

		for (int z_now = z_min; z_now <= z_max; z_now++) {
			z_vect.setEntry(kk - 1, z_now);

			double z_rest = a_hat.getEntry(kk - 1) - z_vect.getEntry(kk - 1);
			double t_value = Math.pow(z_rest, 2) / d_vec.getEntry(kk - 1);
			double exp_t = Math.exp(-0.5 * t_value);

			if (kk > 1) {
				double chi2_new = Chi2 - t_value;
				RealVector z_cond = a_hat.getSubVector(0, kk - 1)
						.subtract(L_mat.getColumnVector(kk - 1).getSubVector(0, kk - 1).mapMultiply(z_rest));

				Object[] nestedResult = estimatorBIE_nested(z_cond, L_mat.getSubMatrix(0, kk - 2, 0, kk - 2),
						d_vec.getSubVector(0, kk - 1), kk - 1, chi2_new, N_int, z_vect, N_max);

				z_BIE = z_BIE.add(((RealVector) nestedResult[0]).mapMultiply(exp_t));
				z_PDF += exp_t * ((double) nestedResult[1]);
				N_int = (int) nestedResult[2];
			} else {
				z_BIE = z_BIE.add(z_vect.mapMultiply(exp_t));
				z_PDF += exp_t;
				N_int++;
			}

			if (N_int > N_max || N_int == 0) {
				return new Object[] { z_BIE, z_PDF, 0 };
			}
		}

		return new Object[] { z_BIE, z_PDF, N_int };
	}

	/**
	 * Partial Ambiguity Resolution (PAR) estimation based on ILS.
	 *
	 * @param aHat     Ambiguity float vector (RealVector).
	 * @param LMat     LtDL-decomposition matrix L (lower unitriangular,
	 *                 RealMatrix).
	 * @param dVec     LtDL-decomposition matrix D (diagonal elements, double
	 *                 array).
	 * @param nCands   Number of best integer solutions (default = 1).
	 * @param minSR    Minimum success rate threshold (default = 0.995).
	 * @param alphaBIE Use BIE estimator if alpha > 0 (default = 0).
	 * @return A result object containing: - aPAR: Partially resolved ambiguity
	 *         vector (RealVector). - nFixed: Number of fixed ambiguity components.
	 *         - srPAR: Success rate of the ambiguity subset.
	 */
	
	public static Object[] estimatorPAR(RealVector aHat, RealMatrix LMat, double[] dVec, int nCands, double minSR,
			double alphaBIE) {
		int n = aHat.getDimension();

		// Default parameter values
		if (nCands <= 0)
			nCands = 1;
		if (minSR <= 0)
			minSR = 0.995;
		if (alphaBIE < 0 || alphaBIE > 1)
			alphaBIE = 0;

		// Compute success rates
		Object[] srResult = SuccessRate.computeSR_IBexact(dVec);
		double srIB = (double) srResult[0];
		double[] srCumul = (double[]) srResult[1];

		// Determine largest subset above success rate threshold
		int kkPAR;
		double srPAR;
		int nFixed;
		if (srIB >= minSR) {
			// Full AR
			kkPAR = 0;
			srPAR = srIB;
			nFixed = n;

		} else if (srCumul[srCumul.length - 1] >= minSR) {
			// Partial AR
			kkPAR = findFirstAboveThreshold(srCumul, minSR);
			srPAR = srCumul[kkPAR];
			nFixed = n - kkPAR;

		} else {
			// No AR
			return new Object[] { aHat, 0, srIB };
		}

		// Find fixed solution for subset with sufficiently high success rate
		RealVector aFixPAR;
		if (alphaBIE > 0 && alphaBIE < 1) {
			double chi2BIE = 2 * GammaIncompleteInverse.gammaincinv(1 - alphaBIE, nFixed / 2.0);
			aFixPAR = (RealVector) estimatorBIE(aHat.getSubVector(kkPAR, n - kkPAR),
					LMat.getSubMatrix(kkPAR, n - 1, kkPAR, n - 1),
					new ArrayRealVector(Arrays.copyOfRange(dVec, kkPAR, n)), chi2BIE)[0];
		} else {
			aFixPAR = (RealVector) estimatorILS(aHat.getSubVector(kkPAR, n - kkPAR),
					LMat.getSubMatrix(kkPAR, n - 1, kkPAR, n - 1), Arrays.copyOfRange(dVec, kkPAR, n), nCands)[0];
		}

		// Compute conditioned float solution for the unresolved subset
		RealMatrix lSub = LMat.getSubMatrix(kkPAR, n - 1, 0, kkPAR - 1);
		RealMatrix lDiag = LMat.getSubMatrix(kkPAR, n - 1, kkPAR, n - 1);
		RealVector residual = aHat.getSubVector(kkPAR, n - kkPAR).subtract(aFixPAR);

		DecompositionSolver solver = new LUDecomposition(lDiag.transpose()).getSolver();
		RealVector adjustment = solver.solve(residual);
		RealVector aCondPAR = aHat.getSubVector(0, kkPAR).subtract(lSub.transpose().operate(adjustment));

		// Combine conditioned and fixed solutions
		RealVector aPAR = new ArrayRealVector(n);
		aPAR.setSubVector(0, aCondPAR);
		aPAR.setSubVector(kkPAR, aFixPAR);

		return new Object[] { aPAR, nFixed, srPAR };
	}

	/**
	 * Finds the first index where the value in the array exceeds the threshold.
	 */
	private static int findFirstAboveThreshold(double[] array, double threshold) {
		for (int i = 0; i < array.length; i++) {
			if (array[i] >= threshold) {
				return i;
			}
		}
		return -1; // Should not occur if threshold logic is correct
	}
	
	
	public static Object[] estimatorIAFFRT(RealVector a_hat, RealMatrix L_mat, RealVector d_vec, double maxFR, Double mu_RATIO) {
        int nn = a_hat.getDimension(); // Problem dimensionality

        // Check if maxFR is not set, use the default 0.1%
        if (maxFR <= 0) {
            maxFR = 0.1 / 100.0; // Default max failure rate is 0.1%
        }

        // If mu_RATIO is provided, consider an arbitrary ratio test
        if (mu_RATIO != null) {
            maxFR = 0; // Arbitrary ratio tests ignore the failure rate
        }

        // Step 1: Compute the two best solutions using ILS estimator
        Object[] ilsOutput = estimatorILS(a_hat, L_mat, d_vec.toArray(), 2);
        RealVector[] a_fix_temp = (RealVector[]) ilsOutput[0]; // Array of RealVectors
        double[] sqnorm_temp = (double[]) ilsOutput[1]; // Squared norms of the solutions

        // Step 2: Compute Success Rate (SR) and Failure Rate (FR)
        double SR = (double) SuccessRate.computeSR_IBexact(d_vec.toArray())[0];
        double FR = 1.0 - SR;

        // Step 3: Check the Failure Rate (FR) against the maximum threshold
        if (FR < maxFR) {
            // If FR is below the threshold, return the best ILS solution
            return new Object[] {a_fix_temp[0], sqnorm_temp[0], nn};
        } else {
            // If FR exceeds the threshold, compute the mu-value
            double mu_value;
            if (mu_RATIO != null) {
                mu_value = mu_RATIO; // Arbitrary ratio test
            } else {
                mu_value = ComputeFFRTCoefficient.computeFFRTcoeff(maxFR, FR, nn); // Fixed-FR ratio test
            }

            // Step 4: Perform the Ratio Test based on the computed mu-value
            if (sqnorm_temp[0] / sqnorm_temp[1] > mu_value) {
                // If the ratio test fails, return the float solution
                return new Object[] {a_hat, 0, 0};
            } else {
                // If the ratio test passes, return the best ILS solution
                return new Object[] {a_fix_temp[0], sqnorm_temp[0], nn};
            }
        }
    }

}
