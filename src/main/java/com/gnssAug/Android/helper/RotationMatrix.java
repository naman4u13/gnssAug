package com.gnssAug.Android.helper;

public class RotationMatrix {

	public static double[] dcm2euler(double[][] dcm) {

		double yaw = Math.atan2(dcm[1][0], dcm[0][0]);
		double pitch = -Math.asin(dcm[2][0]);
		double roll = Math.atan2(dcm[2][1], dcm[2][2]);
		return new double[] { yaw, pitch, roll };
	}

	public static double[][] euler2dcm(double[] euler) {
		double yaw = euler[0];
		double pitch = euler[1];
		double roll = euler[2];

		double[][] dcm = new double[3][3];
		double cp = Math.cos(pitch);
		double sp = Math.sin(pitch);
		double sr = Math.sin(roll);
		double cr = Math.cos(roll);
		double sy = Math.sin(yaw);
		double cy = Math.cos(yaw);

		dcm[0][0] = cp * cy;
		dcm[1][0] = (sr * sp * cy) - (cr * sy);
		dcm[2][0] = (cr * sp * cy) + (sr * sy);
		dcm[0][1] = cp * sy;
		dcm[1][1] = (sr * sp * sy) + (cr * cy);
		dcm[2][1] = (cr * sp * sy) - (sr * cy);
		dcm[0][2] = -sp;
		dcm[1][2] = sr * cp;
		dcm[2][2] = cr * cp;

		return dcm;
	}

	public static void check(double x, double y) {
		double angle = Math.atan2(y, x);
		double angle2 = Math.abs(Math.atan(y / x));
		angle2 += y > 0 ? (x < 0 ? Math.PI - (2 * angle2) : 0) : (x < 0 ? -Math.PI : -(2 * angle2));
		if (Math.abs(angle - angle2) > 0.0001f) {
			System.err.println("ERROR");
		}
	}
}
