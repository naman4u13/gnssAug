package com.gnssAug.Android.utility;

public class Matrix {

	// Only for 3*3 matrix
	public static double[][] multiply(double[][] a, double[][] b) {
		double[][] c = new double[3][3];

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				c[i][j] = 0;
				for (int k = 0; k < 3; k++) {
					c[i][j] += a[i][k] * b[k][j];
				}
			}
		}
		return c;
	}

	public static double[][] getSkewSymMat(double[] a) {
		return new double[][] { { 0, -a[2], a[1] }, { a[2], 0, -a[0] }, { -a[1], a[0], 0 } };
	}
}
