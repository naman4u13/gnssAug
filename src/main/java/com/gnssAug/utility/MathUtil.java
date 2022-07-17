package com.gnssAug.utility;

public class MathUtil {

	public static double getEuclidean(double[] x, double[] y) {
		int n = x.length;
		double dist = 0;
		for (int i = 0; i < 3; i++) {
			dist += Math.pow(x[i] - y[i], 2);
		}
		dist = Math.sqrt(dist);
		return dist;

	}
}
