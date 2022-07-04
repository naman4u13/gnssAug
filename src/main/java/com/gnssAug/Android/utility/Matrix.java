package com.gnssAug.Android.utility;

import java.util.Arrays;

import org.ejml.simple.SimpleMatrix;

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
		return getSkewSymMat(a, false);
	}

	public static double[][] getSkewSymMat(double[] a, boolean isNeg) {
		if (isNeg) {
			Arrays.stream(a).forEach(i -> i = -i);
		}
		return new double[][] { { 0, -a[2], a[1] }, { a[2], 0, -a[0] }, { -a[1], a[0], 0 } };
	}

	public static double[][] matrix2Array(SimpleMatrix matrix) {
		double[][] array = new double[matrix.numRows()][matrix.numCols()];
		for (int r = 0; r < matrix.numRows(); r++) {
			for (int c = 0; c < matrix.numCols(); c++) {
				array[r][c] = matrix.get(r, c);
			}
		}
		return array;
	}
}
