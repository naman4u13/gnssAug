package com.gnssAug.helper.lambdaNew;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.special.Erf;

public class ComputeSR_IBexact {

	/**
	 * Container class for the success rate results.
	 */
	public static class SR_IB {
		private double SR;
		private double[] SR_cumul;
		private double[] SR_vect;

		/**
		 * Constructor for SRResult.
		 *
		 * @param SR       Success rate for Full AR (FAR)
		 * @param SR_cumul Success rate for Partial AR (PAR) with incremental subsets
		 * @param SR_vect  Success rate for each individual (conditioned) component
		 */
		public SR_IB(double SR, double[] SR_cumul, double[] SR_vect) {
			this.SR = SR;
			this.SR_cumul = SR_cumul;
			this.SR_vect = SR_vect;
		}

		/**
		 * Gets the success rate for Full AR (FAR).
		 *
		 * @return SR
		 */
		public double getSR() {
			return SR;
		}

		/**
		 * Gets the success rate for Partial AR (PAR) with incremental subsets.
		 *
		 * @return SR_cumul
		 */
		public double[] getSR_cumul() {
			return SR_cumul;
		}

		/**
		 * Gets the success rate for each individual (conditioned) component.
		 *
		 * @return SR_vect
		 */
		public double[] getSR_vect() {
			return SR_vect;
		}
	}

	/**
	 * Computes the success rate based on an Integer Bootstrapping analytical
	 * formula (exact).
	 *
	 * @param dVec Conditional variances vector
	 * @return SRResult containing SR, SR_cumul, and SR_vect
	 */
	public static SR_IB computeSR_IBexact(double[] dVec) {
		// Problem dimensionality
		int nn = dVec.length;

		// Success rate of each (conditioned) component using IB analytical formula
		double[] SR_vect = new double[nn];
		for (int i = 0; i < nn; i++) {
			SR_vect[i] = Erf.erf(1.0 / Math.sqrt(8.0 * dVec[i]));
			// SR_vect[i] = 2 * normCDF(0.5 / Math.sqrt(dVec[i])) - 1; // Slower & needs
			// toolbox!
		}

		// Success rate (Partial AR) for incremental subsets assuming last-to-first
		double[] SR_cumul = new double[nn];
		double cumulativeProduct = 1.0;
		for (int i = nn - 1; i >= 0; i--) {
			cumulativeProduct *= SR_vect[i];
			SR_cumul[i] = cumulativeProduct;
		}

		// Success-rate for the Full AR (FAR)
		double SR = SR_cumul[0];
		
		
		
		
		double Ps = 1;
        NormalDistribution normalDistribution = new NormalDistribution();
        for (int i = 0;i < dVec.length;i++){
            double cdf = normalDistribution.probability(-Double.MAX_VALUE, 0.5/Math.sqrt(dVec[i]));
            Ps *= (2 * cdf - 1);
        }

		return new SR_IB(SR, SR_cumul, SR_vect);
	}

}