package com.gnssAug.Android.constants;

public class ImuDataSheets {

	/*
	 * In-Run bias is provided in IMU datasheets as zero-g bias(acc) or zero-g
	 * rate(gyro)
	 */

	// The nominal "average" value at Earth's surface, known as standard gravity

	private final static double g = 9.80665;

	public class Pixel4 {
		// in m/s
		public static double accInRunBias = 150 * g * 1e-3;
		// in rad/s
		public static double gyroInRunBias = Math.toRadians(3);
		// in m/s
		public static double accTurnOnBias = 150 * g * 1e-3;
		// in rad/s
		public static double gyroTurnOnBias = Math.toRadians(3);
		// in m/(s^1.5)
		public static double vrw = 180 * g * 1e-6;
		// in rad/(s^0.5)
		public static double arw = Math.toRadians(0.007);
		// in sec
		public static double accBiasCorrTime = 300;
		// in sec
		public static double gyroBiasCorrTime = 300;

	}

}
