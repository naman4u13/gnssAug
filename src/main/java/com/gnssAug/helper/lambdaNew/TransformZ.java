package com.gnssAug.helper.lambdaNew;
import org.apache.commons.math3.linear.*;

import com.gnssAug.helper.lambdaNew.models.TransformZResult;

public class TransformZ {
	/**
	 * Decorrelate ambiguities using an admissible Z-transformation.
	 *
	 * @param L_mat  Old LtDL-decomposition matrix L (lower unitriangular).
	 * @param d_vec  Old LtDL-decomposition diagonal elements.
	 * @param iZt_mat Old inverse transpose Z-transformation matrix (unimodular).
	 *                If null, an identity matrix is used.
	 * @return A result object containing the updated L, d, and iZt matrices.
	 */
	public static TransformZResult transformZ(RealMatrix L_mat, RealVector d_vec, RealMatrix iZt_mat) {
		int n = d_vec.getDimension();

		// Check if iZt_mat is null, initialize it as an identity matrix
		if (iZt_mat == null) {
			iZt_mat = MatrixUtils.createRealIdentityMatrix(n);
		}

		// Iterative loop for swapping and decorrelating adjacent components
		int k = n - 2;
		while (k >= 0) {
			int kp1 = k + 1;

			// Check current pairs {k, k+1} and a correlation-like term L_mat(k+1, k)
			double corr = L_mat.getEntry(kp1, k);
			double mu = Math.round(corr);
			if (mu != 0) {
				corr -= mu;
			}

			// Condition for swapping adjacent ambiguities
			double delta = d_vec.getEntry(k) + corr * corr * d_vec.getEntry(kp1);
			if (delta < d_vec.getEntry(kp1)) {
				// Perform decorrelation for L_mat(k+1, k) if needed
				if (mu != 0) {
					for (int i = kp1; i < n; i++) {
						L_mat.setEntry(i, k, L_mat.getEntry(i, k) - mu * L_mat.getEntry(i, kp1));
					}
					for (int i = 0; i < n; i++) {
						iZt_mat.setEntry(i, kp1, iZt_mat.getEntry(i, kp1) + mu * iZt_mat.getEntry(i, k));
					}

					// Reduce the entire column L_mat(kk+1:nn,kk) for better stability
					for (int i = kp1 + 1; i < n; i++) {
						mu = Math.round(L_mat.getEntry(i, k));
						if (mu != 0) {
							for (int j = i; j < n; j++) {
								L_mat.setEntry(j, k, L_mat.getEntry(j, k) - mu * L_mat.getEntry(j, i));
							}
							for (int j = 0; j < n; j++) {
								iZt_mat.setEntry(j, i, iZt_mat.getEntry(j, i) + mu * iZt_mat.getEntry(j, k));
							}
						}
					}
				}

				// Compute auxiliary variables for performing the adjacent swapping
				double lambda = L_mat.getEntry(kp1, k) * d_vec.getEntry(kp1) / delta;
				double eta = d_vec.getEntry(k) / delta;

				// STEP I: Adjacent swapping operation
				double lkk = -L_mat.getEntry(kp1, k);
				double lkp1 = 1.0;

				for (int j = 0; j < k; j++) {
					double temp1 = lkk * L_mat.getEntry(k, j) + lkp1 * L_mat.getEntry(kp1, j);
					double temp2 = eta * L_mat.getEntry(k, j) + lambda * L_mat.getEntry(kp1, j);
					L_mat.setEntry(k, j, temp1);
					L_mat.setEntry(kp1, j, temp2);
				}

				// STEP II: Update decomposition in the specific swapped block
				L_mat.setEntry(kp1, k, lambda);
				d_vec.setEntry(k, eta * d_vec.getEntry(kp1));
				d_vec.setEntry(kp1, delta);

				// STEP III: Update decomposition in the other conditioned block
				for (int i = kp1 + 1; i < n; i++) {
					double temp1 = L_mat.getEntry(i, k);
					double temp2 = L_mat.getEntry(i, kp1);
					L_mat.setEntry(i, k, temp2);
					L_mat.setEntry(i, kp1, temp1);
				}

				for (int i = 0; i < n; i++) {
					double temp1 = iZt_mat.getEntry(i, k);
					double temp2 = iZt_mat.getEntry(i, kp1);
					iZt_mat.setEntry(i, k, temp2);
					iZt_mat.setEntry(i, kp1, temp1);
				}

				// If a swap took place at lower levels, move up
				if (k < n - 2) {
					k++;
				}
			} else {
				// No swap took place, so move one level down
				k--;
			}
		}

		// Assure that all ambiguity components are decorrelated
		computeIGTRow(L_mat, iZt_mat,null,null);

		return new TransformZResult(L_mat, d_vec, iZt_mat);
	}

	/**
	 * Ensure all ambiguity components are decorrelated by row operations.
	 *
	 * @param L_mat   The L matrix.
	 * @param iZt_mat The inverse transpose Z-transformation matrix.
	 */
	public static void computeIGTRow(RealMatrix L_mat, RealMatrix iZt_mat, Integer ii_min, Integer ii_max) {
	    int n = L_mat.getColumnDimension(); // Dimensionality of the matrix

	    // Initialize iZt_mat as identity if null
	    if (iZt_mat == null) {
	        iZt_mat = MatrixUtils.createRealIdentityMatrix(n);
	    }

	    // Default values for ii_min and ii_max
	    if (ii_min == null || ii_min < 2) {
	        ii_min = 2;
	    }
	    if (ii_max == null || ii_max > n) {
	        ii_max = n;
	    }

	    // Ensure indices are valid
	    if (ii_min > ii_max || ii_min < 2 || ii_max > n) {
	        throw new IllegalArgumentException("Invalid values for ii_min and ii_max.");
	    }

	    // Iterate over each row from ii_min to ii_max
	    for (int ii = ii_min - 1; ii < ii_max; ii++) { // Adjust for 0-based indexing in Java
	        // Round elements of the current row up to ii-1
	        double[] mu_vect = new double[ii];
	        for (int j = 0; j < ii; j++) {
	            mu_vect[j] = Math.round(L_mat.getEntry(ii, j));
	        }

	        // Process non-zero elements in mu_vect
	        for (int jj = 0; jj < mu_vect.length; jj++) {
	            if (mu_vect[jj] != 0) {
	                double mu = mu_vect[jj];

	                // Update L_mat for the ii-th row and below
	                for (int k = ii; k < n; k++) {
	                    L_mat.setEntry(k, jj, L_mat.getEntry(k, jj) - mu * L_mat.getEntry(k, ii));
	                }

	                // Update iZt_mat
	                for (int k = 0; k < n; k++) {
	                    iZt_mat.setEntry(k, ii, iZt_mat.getEntry(k, ii) + mu * iZt_mat.getEntry(k, jj));
	                }
	            }
	        }
	    }
	}
}
