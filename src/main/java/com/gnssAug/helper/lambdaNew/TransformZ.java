package com.gnssAug.helper.lambdaNew;

import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.simple.SimpleMatrix;

/**
 * LAMBDA 4.0 | Decorrelate ambiguities by an admissible Z-transformation
 * 
 * This class provides a method to decorrelate ambiguities by reduction and ordering of 
 * conditional variances. The Z-transformation matrix (unimodular) is then
 * obtained conventionally as the inverse transpose of Z.
 * 
 * -------------------------------------------------------------------------
 * Inputs:
 *   L_mat       Old LtDL-decomposition matrix L (lower unitriangular)
 *   d_vec       Old LtDL-decomposition matrix D (diagonal elements)
 *   iZt_mat     Old inverse transpose Z-transformation matrix (unimodular)
 * 
 * Outputs:
 *   L_mat       New LtDL-decomposition matrix L (lower unitriangular)
 *   d_vec       New LtDL-decomposition matrix D (diagonal elements)
 *   iZt_mat     New inverse transpose Z-transformation matrix (unimodular)
 * 
 * Dependencies:
 *   computeIGT_row
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
public class TransformZ {

    /**
     * Decorrelates ambiguities by reduction and ordering of conditional variances.
     * 
     * @param LMat      Old LtDL-decomposition matrix L (lower unitriangular)
     * @param dVec      Old LtDL-decomposition matrix D (diagonal elements)
     * @param iZtMat    Old inverse transpose Z-transformation matrix (unimodular). If null, identity matrix is assumed.
     * @return          A TransformResult object containing the new L_mat, d_vec, and iZt_mat
     */
    public static TransformResult transformZ(SimpleMatrix LMat, double[] dVec, SimpleMatrix iZtMat) {
        // Problem dimensionality
        int nn = dVec.length;

        // Check number of inputs
        if (LMat == null || dVec == null) {
            throw new IllegalArgumentException("ATTENTION: number of inputs is insufficient!");
        }

        if (iZtMat == null) {
            // Case without any a priori Z-transformation matrix
            iZtMat = SimpleMatrix.identity(nn);
        }

        // NOTE: we use iZt = inv(Z'), assuming a transformation z_hat = Z' * a_hat

        // ALGORITHM: matrix reduction with conditional variances ordering

        // Iterative loop for swapping & decorrelate adjacent components
        int kk = nn - 2;
        while (kk >= 0) {
            int kp1 = kk + 1;
            // ----------------------------------------------------------------------
            // Check current pairs {k,k+1} and a correlation-like term L_mat(kp1,kk)
            double CORR = LMat.get(kp1, kk);
            long muLong = Math.round(CORR);
            double mu = (double) muLong;
            if (muLong != 0) {
                CORR = CORR - mu;
            }
            // ----------------------------------------------------------------------
            // Condition for swapping adjacent ambiguities
            double delta = dVec[kk] + Math.pow(CORR, 2) * dVec[kp1];
            if (delta < dVec[kp1]) {
                // ------------------------------------------------------------------
                // Check if decorrelation for L_mat(kk+1,kk) was needed 
                if (muLong != 0) {
                    // L_mat(kp1:nn, kk) = L_mat(kp1:nn, kk) - mu * L_mat(kp1:nn, kp1);
                    for (int ii = kp1; ii < nn; ii++) {
                        double updatedValue = LMat.get(ii, kk) - mu * LMat.get(ii, kp1);
                        LMat.set(ii, kk, updatedValue);
                    }

                    // iZt_mat(:,kp1) = iZt_mat(:,kp1) + mu * iZt_mat(:,kk);
                    for (int ii = 0; ii < nn; ii++) {
                        double updatedValue = iZtMat.get(ii, kp1) + mu * iZtMat.get(ii, kk);
                        iZtMat.set(ii, kp1, updatedValue);
                    }

                    // Reduce entire column L_mat(kk+1:nn, kk) -> better stability
                    for (int ii = kp1 + 1; ii < nn; ii++) {
                        long muInnerLong = Math.round(LMat.get(ii, kk));
                        double muInner = (double) muInnerLong;
                        if (muInnerLong != 0) {
                            // L_mat(ii:nn, kk) = L_mat(ii:nn, kk) - mu * L_mat(ii:nn, ii);
                            for (int jj = ii; jj < nn; jj++) {
                                double updatedValue = LMat.get(jj, kk) - muInner * LMat.get(jj, ii);
                                LMat.set(jj, kk, updatedValue);
                            }
                            // iZt_mat(:,ii) = iZt_mat(:,ii) + mu * iZt_mat(:,kk);
                            for (int jj = 0; jj < nn; jj++) {
                                double updatedValue = iZtMat.get(jj, ii) + muInner * iZtMat.get(jj, kk);
                                iZtMat.set(jj, ii, updatedValue);
                            }
                        }
                    }
                }
                // ------------------------------------------------------------------
                // Compute auxiliary variables for performing the adjacent swapping
                double lambda = LMat.get(kp1, kk) * dVec[kp1] / delta;    // Auxiliary #1
                double eta = dVec[kk] / delta;                           // Auxiliary #2

                // STEP I: adjacent swapping operation
                // Creating swapMatrix = [ -L_mat(kp1,kk)    1 ;
                //                           eta          lambda ];
                double[] swapRow1 = { -LMat.get(kp1, kk), 1.0 };
                double[] swapRow2 = { eta, lambda };

                // Perform swapMatrix * LMat([kk kp1],1:kk-1)
                for (int col = 0; col < kk; col++) {
                    double val1 = swapRow1[0] * LMat.get(kk, col) + swapRow1[1] * LMat.get(kp1, col);
                    double val2 = swapRow2[0] * LMat.get(kk, col) + swapRow2[1] * LMat.get(kp1, col);
                    LMat.set(kk, col, val1);
                    LMat.set(kp1, col, val2);
                }

                // STEP II: update decomposition in the specific swapped block
                LMat.set(kp1, kk, lambda);
                dVec[kk] = eta * dVec[kp1];
                dVec[kp1] = delta;

                // STEP III: update decomposition in the other conditioned block
                for (int i = kk + 2; i < nn; i++) {
                    double temp1 = LMat.get(i, kk);
                    double temp2 = LMat.get(i, kp1);
                    LMat.set(i, kk, temp2);
                    LMat.set(i, kp1, temp1);
                }
                for (int i = 0; i < nn; i++) {
                    double temp = iZtMat.get(i, kk);
                    double tempNext = iZtMat.get(i, kp1);
                    iZtMat.set(i, kk, tempNext);
                    iZtMat.set(i, kp1, temp);
                }
                // ------------------------------------------------------------------
                // If a swap took place at lower levels, we move up
                if (kk < nn - 2) {
                    kk += 1;
                }
            } else {
                // No swap took place, so we move one level down
                kk -= 1;
            }
            // ----------------------------------------------------------------------
        }

        // Assure that all the ambiguity components are ultimately decorrelated
        // [L_mat,iZt_mat] = computeIGT_row(L_mat,iZt_mat);
        // Assuming computeIGT_row is another method to be implemented.
        IGTResult result = computeIGTRow(LMat, iZtMat);

        return new TransformResult(result.getLMat(),dVec,result.iZtMat);
    }
    
    /**
     * A helper class to hold the result of the computeIGTRow method.
     */
    public static class IGTResult {
        private SimpleMatrix lMat;
        private SimpleMatrix iZtMat;

        public IGTResult(SimpleMatrix lMat, SimpleMatrix iZtMat) {
            this.lMat = lMat;
            this.iZtMat = iZtMat;
        }

        public SimpleMatrix getLMat() {
            return lMat;
        }

        public SimpleMatrix getIZtMat() {
            return iZtMat;
        }
    }

    /**
     * Computes the Integer Gauss Transformations over a range of matrix rows.
     * 
     * @param lMat  Old LtDL-decomposition matrix L (lower unitriangular)
     * @param iZtMat Old inverse transpose of Z-transformation matrix
     * @param iiMin Minimum index of rows to be processed (0-based)
     * @param iiMax Maximum index of rows to be processed (0-based)
     * @return IGTResult containing the new L matrix and inverse transpose Z matrix
     * @throws IllegalArgumentException if iiMin and iiMax are out of valid range
     */
    public static IGTResult computeIGTRow(SimpleMatrix lMat, SimpleMatrix iZtMat, int iiMin, int iiMax) {
        // Problem dimensionality
        int nn = lMat.numCols();

        // Check that input "iiMin" and "iiMax" are correct
        if (iiMin > iiMax || iiMin < 0 || iiMax >= nn) {
            throw new IllegalArgumentException("ATTENTION: something is wrong with \"iiMin\" and \"iiMax\"!");
        }

        // Iterate over each row from "iiMin" till "iiMax" (up -> down)
        for (int ii = iiMin; ii <= iiMax; ii++) {
            // Round elements of current row "ii"
            double[] muVect = new double[ii];
            for (int j = 0; j < ii; j++) {
                muVect[j] = Math.round(lMat.get(ii, j));
            }

            // Find indices where muVect is not zero
            int[] indexMu = findNonZeroIndices(muVect);

            // At the ii-th row, process columns defined in "indexMu"
            for (int jj : indexMu) {
                double mu = muVect[jj];
                for (int row = ii; row < nn; row++) {
                    double updatedValue = lMat.get(row, jj) - mu * lMat.get(row, ii);
                    lMat.set(row, jj, updatedValue);
                }
                for (int row = 0; row < iZtMat.numRows(); row++) {
                    double updatedValue = iZtMat.get(row, ii) + mu * iZtMat.get(row, jj);
                    iZtMat.set(row, ii, updatedValue);
                }
            }
        }

        return new IGTResult(lMat, iZtMat);
    }

    /**
     * Finds the indices of the non-zero elements in the given vector.
     * 
     * @param vector The input vector
     * @return An array of indices where the vector elements are non-zero
     */
    private static int[] findNonZeroIndices(double[] vector) {
        return java.util.stream.IntStream.range(0, vector.length)
                .filter(i -> vector[i] != 0)
                .toArray();
    }

    /**
     * Overloaded method to computeIGTRow when no previous Z-transformation was performed.
     * Initializes iZtMat as the identity matrix and sets default iiMin and iiMax.
     * 
     * @param lMat Old LtDL-decomposition matrix L (lower unitriangular)
     * @return IGTResult containing the new L matrix and inverse transpose Z matrix
     */
    public static IGTResult computeIGTRow(SimpleMatrix lMat) {
        int nn = lMat.numCols();
        SimpleMatrix iZtMat = SimpleMatrix.identity(nn);
        int iiMin = 0; // Changed to zero-based
        int iiMax = nn - 1; // Changed to zero-based
        return computeIGTRow(lMat, iZtMat, iiMin, iiMax);
    }

    /**
     * Overloaded method to computeIGTRow when only one of iiMin or iiMax is missing.
     * Sets the missing parameters to default values.
     * 
     * @param lMat  Old LtDL-decomposition matrix L (lower unitriangular)
     * @param iZtMat Old inverse transpose of Z-transformation matrix
     * @param iiMin Minimum index of rows to be processed (if provided, 0-based)
     * @return IGTResult containing the new L matrix and inverse transpose Z matrix
     */
    public static IGTResult computeIGTRow(SimpleMatrix lMat, SimpleMatrix iZtMat) {
        int nn = lMat.numCols();
        int iiMin = 0; // Changed to zero-based
        int iiMax = nn - 1; // Changed to zero-based
        return computeIGTRow(lMat, iZtMat, iiMin, iiMax);
    }
    
    /**
     * A class to hold the results of the transformZ method.
     */
    public static class TransformResult {
        private SimpleMatrix Lmat;
        private double[] dVec;
        private SimpleMatrix iZtMat;

        public TransformResult(SimpleMatrix Lmat, double[] dVec, SimpleMatrix iZtMat) {
            this.Lmat = Lmat;
            this.dVec = dVec;
            this.iZtMat = iZtMat;
        }

		public SimpleMatrix getLmat() {
			return Lmat;
		}

		public double[] getdVec() {
			return dVec;
		}

		public SimpleMatrix getiZtMat() {
			return iZtMat;
		}
    }
}
