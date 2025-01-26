package com.gnssAug.helper.lambdaNew.Estimators;

import org.ejml.data.DMatrixRMaj;
import org.ejml.simple.SimpleMatrix;
import java.util.Arrays;
import java.util.Comparator;

/**
 * LAMBDA 4.0 | Integer Least-Squares (ILS) estimator by search-and-shrink
 * This class computes a 'fixed' solution based on Integer Least-Squares
 * (ILS-)estimator using the search-and-shrink approach. The latter adopts 
 * the algorithm proposed by Ghasemmehdi and Agrell (2011, [RD01]).
 *
 * -------------------------------------------------------------------------
 * INPUTS
 *   aHat        Ambiguity float vector (column)
 *   LMat        LtDL-decomposition matrix L (lower unitriangular)
 *   dVec        LtDL-decomposition matrix D (diagonal elements)
 *   nCands      Number of best integer solutions [DEFAULT = 1]
 *
 * OUTPUTS
 *  aFix        Ambiguity fixed vector (column) - Best ILS solution(s)
 *  sqNorm      Squared norm for each ILS solution
 *
 * DEPENDENCIES:
 *   none
 *
 * REFERENCES
 *   [RD01] A. Ghasemmehdi and E. Agrell, "Faster Recursions in Sphere 
 *       Decoding" in IEEE Transactions on Information Theory, vol. 57, 
 *       no. 6, pp. 3530-3536, June 2011. 
 *       DOI: 10.1109/TIT.2011.2143830.
 *
 * -------------------------------------------------------------------------
 * Copyright: Geoscience & Remote Sensing department @ TUDelft | 01/06/2024
 * Contact email:    LAMBDAtoolbox-CITG-GRS@tudelft.nl
 * -------------------------------------------------------------------------
 * Created by
 *   01/06/2024  - Lotfi Massarweh
 *       Implementation for LAMBDA 4.0 toolbox
 *
 * Modified by
 *   dd/mm/yyyy  - Name Surname author - email address
 *       >> Changes made in this new version
 * -------------------------------------------------------------------------
 */
public class EstimatorILS_old {

    /**
     * Container class for the ILS estimation results.
     */
    public static class ILSResult {
        private final SimpleMatrix aFix;
        private final double[] sqNorm;

        public ILSResult(SimpleMatrix aFix, double[] sqNorm) {
            this.aFix = aFix;
            this.sqNorm = sqNorm;
        }

        public SimpleMatrix getAFix() {
            return aFix;
        }

        public double[] getSqNorm() {
            return sqNorm;
        }
    }

    /**
     * Computes a 'fixed' solution based on Integer Least-Squares (ILS) estimator using the search-and-shrink approach.
     *
     * @param aHat   Ambiguity float vector (column)
     * @param LMat   LtDL-decomposition matrix L (lower unitriangular)
     * @param dVec   LtDL-decomposition matrix D (diagonal elements)
     * @param nCands Number of best integer solutions [DEFAULT = 1]
     * @return ILSResult containing the best ILS solution(s) and their corresponding squared norms
     * @throws IllegalArgumentException if the number of input arguments is insufficient
     */
    public static ILSResult estimatorILS(SimpleMatrix aHat, SimpleMatrix LMat, double[] dVec, Integer nCands) {
        // Problem dimensionality
        int nn = aHat.numRows();

        // Check number of input arguments
        if (aHat == null || LMat == null || dVec == null) {
            throw new IllegalArgumentException("ATTENTION: number of inputs is insufficient!");
        }

        // Set default number of candidates if not provided
        int numCands = (nCands != null) ? nCands : 1;

        // Determine which level to move to after zCond[0] is chosen at level 1.
        int k0;
        if (numCands == 1 && nn > 1) {
            k0 = 1; // Cannot find a better candidate, so directly try level 2 (index 1)
        } else {
            k0 = 0; // Try to further improve candidates at level 1 (index 0)
        }

        // Initialization output variables 
        SimpleMatrix aFix = new SimpleMatrix(nn, numCands);
        double[] sqNorm = new double[numCands];

        // Initial number of candidate solutions
        int intCount = 0;
        int intMax = numCands;

        // Start from an ellipsoid with infinite radius
        double maxChi2 = Double.POSITIVE_INFINITY;

        // Initialization variables used
        SimpleMatrix aCond = new SimpleMatrix(nn, 1);
        double[] zCond = new double[nn];
        double[] left = new double[nn];
        int[] step = new int[nn];

        // Initialization at the n-th level (last index nn-1)
        aCond.set(nn - 1, 0, aHat.get(nn - 1, 0));
        zCond[nn - 1] = Math.round(aCond.get(nn - 1, 0));
        left[nn - 1] = aCond.get(nn - 1, 0) - zCond[nn - 1];
        step[nn - 1] = Integer.signum((int) left[nn - 1]);

        // NOTE: very rarely, we need a positive step to avoid stall in the case an
        // exact integer value is provided for "aHat(nn)"
        if (step[nn - 1] == 0) {
            step[nn - 1] = 1;
        }

        // Used to compute conditional ambiguities
        SimpleMatrix S = new SimpleMatrix(nn, nn);

        // Initializing the variable "dist[k] = sum_{j=k+1}^{n} (a_j - aCond_j)^2 / d_j"
        double[] dist = new double[nn + 1];
        Arrays.fill(dist, 0.0);

        // Additional variables needed for keeping track of the conditional update
        int[] path = new int[nn];
        Arrays.fill(path, nn - 1);

        // Algorithm GHAH (Ghasemmehdi and Agrell, 2011; [RD01])
        int kk = nn - 1; // Start main search-loop from the last ambiguity component

        // Iterative search
        boolean endSearch = false;
        while (!endSearch) {

            // Current (partial) distance of a candidate solution
            double newDist = dist[kk] + (left[kk] * left[kk]) / dVec[kk];

            // Keep moving down if current (partial) distance is smaller than radius
            while (newDist < maxChi2) {
                if (kk != 0) {
                    // Move down to level "k-1"
                    kk = kk - 1;
                    dist[kk] = newDist;

                    // Conditionally update recalling previous updates by "path"
                    for (int jj = path[kk]; jj >= kk + 1; jj--) {
                        S.set(jj - 1, kk, S.get(jj, kk) - left[jj] * LMat.get(jj, kk));
                    }

                    aCond.set(kk, 0, aHat.get(kk, 0) + S.get(kk, kk));
                    zCond[kk] = Math.round(aCond.get(kk, 0));
                    left[kk] = aCond.get(kk, 0) - zCond[kk];
                    step[kk] = Integer.signum((int) left[kk]);

                    // NOTE: very rarely, we need a positive step to avoid stall in 
                    // the case an exact integer value is found for "aCond(kk)"
                    if (step[kk] == 0) {
                        step[kk] = 1;
                    }
                } else {
                    // Store the candidate found and try next valid integer
                    if (intCount < numCands - 1) {
                        intCount = intCount + 1;
                        for (int i = 0; i < nn; i++) {
                            aFix.set(i, intCount-1, zCond[i]);
                        }
                        sqNorm[intCount-1] = newDist;
                    } else {
                        for (int i = 0; i < nn; i++) {
                            aFix.set(i, intMax - 1, zCond[i]);
                        }
                        sqNorm[intMax - 1] = newDist;
                        double currentMaxChi2 = sqNorm[intMax - 1];
                        for (int i = 0; i < sqNorm.length; i++) {
                            if (sqNorm[i] > currentMaxChi2) {
                                currentMaxChi2 = sqNorm[i];
                                intMax = i + 1;
                            }
                        }
                        maxChi2 = currentMaxChi2;
                    }

                    // Next valid integer (kk+1 level)
                    kk = k0;
                    zCond[kk] = zCond[kk] + step[kk];
                    left[kk] = aCond.get(kk, 0) - zCond[kk];
                    step[kk] = -step[kk] - Integer.signum(step[kk]);

                }
                newDist = dist[kk] + (left[kk] * left[kk]) / dVec[kk];
            }
            int iLevel = kk;

            // Exit or move up
            while (newDist >= maxChi2) {
                if (kk == nn - 1) {
                    endSearch = true;
                    break;
                }
                kk = kk + 1; // Move up to level "kk+1"
                zCond[kk] = zCond[kk] + step[kk]; // Next valid integer
                left[kk] = aCond.get(kk, 0) - zCond[kk];
                step[kk] = -step[kk] - Integer.signum(step[kk]);
                newDist = dist[kk] + (left[kk] * left[kk]) / dVec[kk];
            }

            if (!endSearch) {
                // Define "path" for the successive conditional update
                for (int j = iLevel; j < kk; j++) {
                    path[j] = kk;
                }
                for (int jj = iLevel - 1; jj >= 0; jj--) {
                    if (path[jj] < kk) {
                        path[jj] = kk;
                    } else {
                        break; // Exit from this for-loop
                    }
                }
            }
        }

        // Sort the solutions by their corresponding residuals' squared norm
        // Create an array of indices
        Integer[] indices = new Integer[sqNorm.length];
        for (int i = 0; i < indices.length; i++) {
            indices[i] = i;
        }

        // Sort indices based on sqNorm values
        Arrays.sort(indices, Comparator.comparingDouble(i -> sqNorm[i]));

        // Create sorted sqNorm and aFix
        double[] sortedSqNorm = new double[sqNorm.length];
        SimpleMatrix sortedAFix = new SimpleMatrix(nn, sqNorm.length);
        for (int i = 0; i < indices.length; i++) {
            sortedSqNorm[i] = sqNorm[indices[i]];
            for (int j = 0; j < nn; j++) {
                sortedAFix.set(j, i, aFix.get(j, indices[i]));
            }
        }

        // If only nCands solutions are needed, truncate
        if (numCands < sortedSqNorm.length) {
            sortedSqNorm = Arrays.copyOfRange(sortedSqNorm, 0, numCands);
            sortedAFix = sortedAFix.extractMatrix(0, nn, 0, numCands);
        }

        return new ILSResult(sortedAFix, sortedSqNorm);
    }
}

