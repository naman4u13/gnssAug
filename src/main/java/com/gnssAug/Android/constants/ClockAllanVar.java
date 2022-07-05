package com.gnssAug.Android.constants;

public class ClockAllanVar {

	private static final double SpeedofLight = 299792458;
	private static final double c2 = SpeedofLight * SpeedofLight;

	public class TCXO_low_qaulity {
		// Typical Allan Variance Coefficients for TCXO (low quality)
		public static final double h0 = 2E-19;
		public static final double h_2 = 2E-20;
		public static final double sf = c2 * h0 / 2;
		public static final double sg = 2 * Math.PI * Math.PI * c2 * h_2;
	}
}
