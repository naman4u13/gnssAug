package com.gnssAug.helper.lambdaNew;
import org.apache.commons.math3.special.Erf;
public class SuccessRate {
	
	public static Object[] computeSR_IBexact(double[] dVec) {
        int n = dVec.length;

        double[] srVect = new double[n];
        double[] srCumul = new double[n];

        // Success rate of each conditioned component
        for (int i = 0; i < n; i++) {
            srVect[i] = Erf.erf(1.0 / Math.sqrt(8.0 * dVec[i]));
        }

        // Cumulative success rate (Partial AR) for incremental subsets (last-to-first order)
        srCumul[n - 1] = srVect[n - 1];
        for (int i = n - 2; i >= 0; i--) {
            srCumul[i] = srCumul[i + 1] * srVect[i];
        }

        // Success rate for Full AR (FAR)
        double sr = srCumul[0];

        return new Object[] {sr,srCumul,srVect};
    }

}
