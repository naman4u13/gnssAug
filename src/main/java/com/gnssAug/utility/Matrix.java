package com.gnssAug.utility;

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
		// Make a deep copy, so that source array does not get impacted
		double[] b = new double[] { a[0], a[1], a[2] };
		if (isNeg) {
			for (int i = 0; i < 3; i++) {
				b[i] *= -1;
			}

		}
		return new double[][] { { 0, -b[2], b[1] }, { b[2], 0, -b[0] }, { -b[1], b[0], 0 } };
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

	public static SimpleMatrix getProjection(SimpleMatrix A, SimpleMatrix W) {
		SimpleMatrix At = A.transpose();
		SimpleMatrix P_A = A.mult((At.mult(W).mult(A)).invert()).mult(At).mult(W);
		return P_A;
	}

	public static SimpleMatrix getPerpendicularProjection(SimpleMatrix A, SimpleMatrix W) {

		SimpleMatrix P_A = getProjection(A, W);
		int n = P_A.numRows();
		SimpleMatrix P_A_perpendicular = SimpleMatrix.identity(n).minus(P_A);
		return P_A_perpendicular;
	}

	public static double getNorm(SimpleMatrix A, SimpleMatrix B) {
		double c = A.transpose().mult(B.invert()).mult(A).get(0);
		return c;
	}

	public static SimpleMatrix getPseudoInv(SimpleMatrix A, SimpleMatrix W) {
		SimpleMatrix At = A.transpose();
		SimpleMatrix Aplus = ((At.mult(W).mult(A)).invert()).mult(At).mult(W);
		return Aplus;
	}
}
