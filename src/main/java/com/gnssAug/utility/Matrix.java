package com.gnssAug.utility;

import java.util.ArrayList;
import java.util.List;

import org.ejml.data.Complex_F64;
import org.ejml.simple.SimpleEVD;
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

	// Only for 3*3 matrix
	public static double[][] scale(double[][] a, double scaleFact) {
		int n = a.length;
		int m = a[1].length;
		double[][] b = new double[n][m];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				b[i][j] = a[i][j] * scaleFact;
			}
		}
		return b;
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
	
	public static double[] matrix2ArrayVec(SimpleMatrix matrix) {
		double[] array = new double[matrix.numRows()];
		for (int r = 0; r < matrix.numRows(); r++) {
			
				array[r] = matrix.get(r, 0);
			
		}
		return array;
	}
	
	public static SimpleMatrix ArrayList2Vector(ArrayList<Integer> list)
	{
		int n = list.size();
		double[] data = list.stream().mapToDouble(i -> i).toArray();
		SimpleMatrix vec = new SimpleMatrix(n,1,true,data);
		return vec;
	}

	public static SimpleMatrix getProjection(SimpleMatrix A, SimpleMatrix W) {
		SimpleMatrix At = A.transpose();
		SimpleMatrix P_A = null;
		P_A = A.mult((At.mult(W).mult(A)).invert()).mult(At).mult(W);
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
	
	public static boolean areMatricesEqual(SimpleMatrix matrixA, SimpleMatrix matrixB, double tolerance) {
        // Check dimensions first
        if (matrixA.numRows() != matrixB.numRows() || matrixA.numCols() != matrixB.numCols()) {
            return false; // Matrices are not of the same size
        }
        // Use isIdentical to check element-wise equality within the tolerance
        return matrixA.isIdentical(matrixB, tolerance);
    }
	
	/**
     * Method to compute the correlation coefficient matrix from a covariance matrix
     * @param covMatrix The covariance matrix (as a SimpleMatrix)
     * @return The correlation coefficient matrix
     */
    public static SimpleMatrix computeCorrelationMatrix(SimpleMatrix covMatrix) {
        // Get the number of variables (i.e., the number of rows/columns)
        int n = covMatrix.numRows();

        // Create a new matrix to store the correlation coefficient matrix
        SimpleMatrix corrMatrix = new SimpleMatrix(n, n);

        // Calculate the standard deviations (square roots of diagonal elements)
        SimpleMatrix stdDevs = new SimpleMatrix(n, 1);
        for (int i = 0; i < n; i++) {
            stdDevs.set(i, 0, Math.sqrt(covMatrix.get(i, i))); // Square root of diagonal element
        }

        // Compute the correlation matrix
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double cov = covMatrix.get(i, j); // Covariance value
                double stdDevProduct = stdDevs.get(i, 0) * stdDevs.get(j, 0); // Product of standard deviations
                double corr = cov / stdDevProduct; // Correlation coefficient
                corrMatrix.set(i, j, corr);
            }
        }

        return corrMatrix;
    }
    
    public static double computeDecorrelationNumber(SimpleMatrix varianceMatrix) {
        // Perform eigenvalue decomposition
        List<Complex_F64> eigenvalues = new SimpleEVD(varianceMatrix.getMatrix()).getEigenvalues();

        
        // Find the maximum and minimum eigenvalues
        double lambdaMax = Double.NEGATIVE_INFINITY;
        double lambdaMin = Double.POSITIVE_INFINITY;

        for ( int i=0;i<eigenvalues.size();i++) {
        	Double eigenvalue = eigenvalues.get(i).getMagnitude();
            if (eigenvalue > lambdaMax) {
                lambdaMax = eigenvalue;
            }
            if (eigenvalue < lambdaMin) {
                lambdaMin = eigenvalue;
            }
        }

        // Compute and return the decorrelation number (ratio of max to min eigenvalue)
        return lambdaMax / lambdaMin;
    }
}
