package com.gnssAug.helper.lambdaNew.Estimators;

import org.ejml.simple.SimpleMatrix;

import com.gnssAug.helper.lambdaNew.Estimators.EstimatorILS_old.ILSResult;

//Integer ambiguity vector search by employing the search-and-shrink technique.
public class EstimatorILS {
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

    public static ILSResult estimatorILS(SimpleMatrix ahat, SimpleMatrix L, double[] D, Integer ncands) {
        int n = ahat.numRows();
        SimpleMatrix afixed = new SimpleMatrix(n,ncands);
        double[] sqnorm = new double[ncands];

        double Chi2 = Double.MAX_VALUE;
        double[] dist = new double[n];
        boolean endsearch = false;
        int count = 0;

        double[] acond = new double[n];
        acond[n-1] = ahat.get(n-1,0);
        double[] zcond = new double[n];
        zcond[n-1] = Math.round(acond[n-1]);
        double left = acond[n-1] - zcond[n-1];
        double[] step = new double[n];
        step[n-1] = Math.signum(left);
        if(step[n-1] == 0.0) step[n-1] = 1.0;

        int imax = ncands;

        SimpleMatrix S = new SimpleMatrix(n,n);
        int k = n;

        while (!endsearch){
            double newdist = dist[k-1] + left*left/D[k-1];
            if(newdist < Chi2){
                if (k != 1){
                    k = k - 1;
                    dist[k-1] = newdist;
                 // Update the row (k-1) of matrix S from row k of S and L
                    SimpleMatrix rowK_S = S.extractMatrix(k, k + 1, 0, k); // Row k of S (1xk matrix)
                    SimpleMatrix rowK_L = L.extractMatrix(k, k + 1, 0, k); // Row k of L (1xk matrix)
                    double multiplier = zcond[k] - acond[k];               // Scalar multiplier

                    // Compute the updated row (k-1)
                    SimpleMatrix updatedRow = rowK_S.plus(rowK_L.scale(multiplier));

                    // Set the updated row (k-1) back into S
                    S.insertIntoThis(k - 1, 0, updatedRow);
                    acond[k-1] = ahat.get(k-1,0) + S.get(k-1,k-1);
                    zcond[k-1] = Math.round(acond[k-1]);
                    left = acond[k-1] - zcond[k-1];
                    step[k-1] = Math.signum(left);
                    if(step[k-1] == 0.0) step[k-1] = 1.0;
                }else{
                    if(count < ncands - 1){
                        count = count + 1;
                     // Convert zcond to a column SimpleMatrix
                        SimpleMatrix zcondMatrix = new SimpleMatrix(n, 1, true, zcond);

                        // Set column (count - 1) of aFixed with zcondMatrix
                        afixed.insertIntoThis(0, count - 1, zcondMatrix);
                        sqnorm[count - 1] = newdist;
                    }else{
                    	// Convert zcond to a column SimpleMatrix
                    	SimpleMatrix zcondMatrix = new SimpleMatrix(n, 1, true, zcond);

                    	// Set column (imax - 1) of aFixed with zcondMatrix
                    	afixed.insertIntoThis(0, imax - 1, zcondMatrix);
                        sqnorm[imax-1] = newdist;
                        Chi2 = sqnorm[0];
                        imax = 1;
                        for(int i = 1;i < sqnorm.length;i++){
                            if(sqnorm[i] > Chi2){
                                Chi2 = sqnorm[i];
                                imax = i + 1;
                            }
                        }
                    }
                    zcond[0] = zcond[0] + step[0];
                    left = acond[0] - zcond[0];
                    step[0] = -step[0] - Math.signum(step[0]);
                }
            }else{
                if(k == n){
                    endsearch = true;
                }else{
                    k += 1;
                    zcond[k-1] = zcond[k-1] + step[k-1];
                    left = acond[k-1] - zcond[k-1];
                    step[k-1] = -step[k-1] - Math.signum(step[k-1]);
                }
            }
        }
        int[] order = arraySort(sqnorm);
        // Reorder the columns of afixed based on the "order" array
        SimpleMatrix reorderedAFixed = new SimpleMatrix(n, ncands);
        for (int i = 0; i < ncands; i++) {
            // Extract the column corresponding to the "order[i]" index
            SimpleMatrix column = afixed.extractVector(false, order[i]);
            // Set the reordered column into the "reorderedAFixed" matrix
            reorderedAFixed.insertIntoThis(0, i, column);
        }

        // Update afixed with the reordered matrix
        afixed = reorderedAFixed;
        return new ILSResult(afixed, sqnorm);
    }

    private static int[] arraySort(double[] arr) {
        double temp;
        int index;
        int k=arr.length;
        int[]Index= new int[k];
        for(int i=0;i<k;i++) {
            Index[i]=i;
        }

        for(int i=0;i<arr.length;i++) {
            for(int j=0;j<arr.length-i-1;j++) {
                if(arr[j]>arr[j+1]) {
                    temp = arr[j];
                    arr[j] = arr[j+1];
                    arr[j+1] = temp;

                    index=Index[j];
                    Index[j] = Index[j+1];
                    Index[j+1] = index;
                }
            }
        }
        return Index;
    }

}
