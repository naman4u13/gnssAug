package com.gnssAug.utility;

import java.util.ArrayList;

public class MathUtil {

	public static double getEuclidean(double[] x, double[] y) {

		double dist = 0;
		for (int i = 0; i < 3; i++) {
			dist += Math.pow(x[i] - y[i], 2);
		}
		dist = Math.sqrt(dist);
		return dist;

	}

	public static long getFact(long x) {
		long fact = 1;
		while (x != 0) {
			fact *= x;
			x--;
		}
		return fact;
	}

	public static long getCombCount(int n, int m) {
		long count = getFact(n) / (getFact(m) * getFact(n - m));
		return count;
	}

	public static double RMS(ArrayList<Double> list) {
		return Math.sqrt(list.stream().mapToDouble(x -> x * x).average().orElse(Double.NaN));
	}

	public static double MAE(ArrayList<Double> list) {
		return list.stream().mapToDouble(x -> x).average().orElse(Double.NaN);
	}
}
