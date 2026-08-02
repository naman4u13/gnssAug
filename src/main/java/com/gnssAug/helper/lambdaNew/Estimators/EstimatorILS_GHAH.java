package com.gnssAug.helper.lambdaNew.Estimators;

import org.ejml.simple.SimpleMatrix;

/**
 * LAMBDA 4.0 | ILS estimator using the GHAH (Ghasemmehdi & Agrell 2011) algorithm.
 *
 * Extends EstimatorILS and overrides estimatorILS() with two optimisations
 * over the plain search-and-shrink in the parent class:
 *
 *   1. k0 shortcut: when ncands==1, after a candidate is found at level 1,
 *      jump back to level 2 instead of level 1.  The nearest integer at level 1
 *      is always round(acond[1]), so there is no benefit in trying others.
 *
 *   2. path array + column-wise S update: instead of recomputing the full S-row
 *      on every descent, only the single diagonal element S[kk-1][kk-1] is updated
 *      by accumulating contributions column-by-column through the path pointer.
 *      This avoids O(n) work per descent step when levels are skipped during ascent.
 *
 * Reference: A. Ghasemmehdi and E. Agrell, "Faster Recursions in Sphere Decoding",
 *            IEEE Trans. Inf. Theory, vol. 57, no. 6, pp. 3530-3536, June 2011.
 *            DOI: 10.1109/TIT.2011.2143830
 *
 * All inputs, outputs, and the ILSResult type are identical to EstimatorILS.
 * To revert to the original algorithm, replace "EstimatorILS_GHAH" with
 * "EstimatorILS" at every call site (no other changes required).
 */
public class EstimatorILS_GHAH extends EstimatorILS {

    @Override
    public ILSResult estimatorILS(SimpleMatrix ahat, SimpleMatrix L, double[] D, Integer ncands) {
        int n = ahat.numRows();
        SimpleMatrix afixed = new SimpleMatrix(n, ncands);
        double[] sqnorm = new double[ncands];

        // k0: after storing a candidate at level 1, jump back here.
        // For ncands==1, level 1's best integer is always round(acond[1]) — no need to try others.
        int k0 = (ncands == 1 && n > 1) ? 2 : 1;

        int intCount = 0;
        int intMax = ncands;
        double maxChi2 = Double.MAX_VALUE;

        // Arrays indexed as [kk-1] for 1-indexed level kk
        double[] acond = new double[n];
        double[] zcond = new double[n];
        double[] left  = new double[n];
        double[] step  = new double[n];

        // dist[kk-1] = partial squared norm: sum_{j=kk+1}^{n} left[j-1]^2 / D[j-1]
        // dist[n-1] = 0 always (top level has no contribution from above).
        // Java initialises double[] to 0, so dist[n-1] is already correct.
        double[] dist = new double[n];

        // S matrix: only the diagonal element S[kk-1][kk-1] is read per level.
        // GHAH populates one column at a time via the path array.
        SimpleMatrix S = new SimpleMatrix(n, n);

        // path[kk-1] = the highest 1-indexed level from which we last descended through level kk.
        // Initialised to n so the first full descent propagates from the top.
        int[] path = new int[n];
        for (int i = 0; i < n; i++) path[i] = n;

        // Initialise at the top level
        int kk = n;
        acond[kk - 1] = ahat.get(kk - 1, 0);
        zcond[kk - 1] = Math.round(acond[kk - 1]);
        left[kk - 1]  = acond[kk - 1] - zcond[kk - 1];
        step[kk - 1]  = Math.signum(left[kk - 1]);
        if (step[kk - 1] == 0.0) step[kk - 1] = 1.0;

        boolean endSearch = false;

        while (!endSearch) {
            double newDist = dist[kk - 1] + left[kk - 1] * left[kk - 1] / D[kk - 1];

            // ---- Inner descent loop: keep moving down while inside the ellipsoid ----
            while (newDist < maxChi2) {
                if (kk != 1) {
                    // Move down to level kk-1
                    kk--;
                    dist[kk - 1] = newDist;

                    // GHAH S-column update.
                    // MATLAB (1-indexed): for jj = path(kk):-1:kk+1
                    //                        S(jj-1,kk) = S(jj,kk) - left(jj)*L_mat(jj,kk)
                    // Java  (0-indexed): S[jj-2][kk-1] = S[jj-1][kk-1] - left[jj-1]*L[jj-1][kk-1]
                    for (int jj = path[kk - 1]; jj >= kk + 1; jj--) {
                        S.set(jj - 2, kk - 1,
                              S.get(jj - 1, kk - 1) - left[jj - 1] * L.get(jj - 1, kk - 1));
                    }

                    acond[kk - 1] = ahat.get(kk - 1, 0) + S.get(kk - 1, kk - 1);
                    zcond[kk - 1] = Math.round(acond[kk - 1]);
                    left[kk - 1]  = acond[kk - 1] - zcond[kk - 1];
                    step[kk - 1]  = Math.signum(left[kk - 1]);
                    if (step[kk - 1] == 0.0) step[kk - 1] = 1.0;

                } else {
                    // kk == 1: a complete integer candidate has been found — store it
                    if (intCount < ncands - 1) {
                        intCount++;
                        SimpleMatrix zcondMatrix = new SimpleMatrix(n, 1, true, zcond);
                        afixed.insertIntoThis(0, intCount - 1, zcondMatrix);
                        sqnorm[intCount - 1] = newDist;
                    } else {
                        SimpleMatrix zcondMatrix = new SimpleMatrix(n, 1, true, zcond);
                        afixed.insertIntoThis(0, intMax - 1, zcondMatrix);
                        sqnorm[intMax - 1] = newDist;
                        maxChi2 = sqnorm[0];
                        intMax = 1;
                        for (int i = 1; i < sqnorm.length; i++) {
                            if (sqnorm[i] > maxChi2) {
                                maxChi2 = sqnorm[i];
                                intMax = i + 1;
                            }
                        }
                    }

                    // GHAH k0 shortcut: jump to k0 and enumerate the next integer there
                    kk = k0;
                    zcond[kk - 1] += step[kk - 1];
                    left[kk - 1]   = acond[kk - 1] - zcond[kk - 1];
                    step[kk - 1]   = -step[kk - 1] - Math.signum(step[kk - 1]);
                }

                newDist = dist[kk - 1] + left[kk - 1] * left[kk - 1] / D[kk - 1];
            }

            int iLevel = kk;

            // ---- Inner ascent loop: move up until back inside the ellipsoid ----
            while (newDist >= maxChi2) {
                if (kk == n) {
                    endSearch = true;
                    break;
                }
                kk++;
                zcond[kk - 1] += step[kk - 1];
                left[kk - 1]   = acond[kk - 1] - zcond[kk - 1];
                step[kk - 1]   = -step[kk - 1] - Math.signum(step[kk - 1]);
                newDist = dist[kk - 1] + left[kk - 1] * left[kk - 1] / D[kk - 1];
            }

            // Update path: levels iLevel..kk-1 now source their S-column update from kk.
            // Then propagate downward to any earlier levels whose path is stale (< kk).
            // MATLAB: path(iLevel:kk-1) = kk
            for (int i = iLevel; i <= kk - 1; i++) {
                path[i - 1] = kk;
            }
            // MATLAB: for jj = iLevel-1:-1:1; if path(jj)<kk; path(jj)=kk; else break; end
            for (int jj = iLevel - 1; jj >= 1; jj--) {
                if (path[jj - 1] < kk) {
                    path[jj - 1] = kk;
                } else {
                    break;
                }
            }
        }

        // Sort candidates by ascending squared norm (reuses parent's arraySort)
        int[] order = arraySort(sqnorm);
        SimpleMatrix reorderedAFixed = new SimpleMatrix(n, ncands);
        for (int i = 0; i < ncands; i++) {
            SimpleMatrix column = afixed.extractVector(false, order[i]);
            reorderedAFixed.insertIntoThis(0, i, column);
        }
        afixed = reorderedAFixed;

        return new ILSResult(afixed, sqnorm);
    }
}
