package com.gnssAug.utility;

public class Vector {

	public static double[] crossProd(double[] a, double[] b) {
		double[] c = new double[3];
		c[0] = (a[1] * b[2]) - (a[2] * b[1]);
		c[1] = -((a[0] * b[2]) - (a[2] * b[0]));
		c[2] = (a[0] * b[1]) - (a[1] * b[0]);
		return c;
	}

	public static double dotProd(double[] a, double[] b) {
		double dotProd = 0;
		for (int i = 0; i < 3; i++) {
			dotProd += a[i] * b[i];
		}
		return dotProd;
	}

	public static double[] add(double[] a, double[] b) {
		int n = a.length;
		double[] sum = new double[n];
		for (int i = 0; i < n; i++) {
			sum[i] = a[i] + b[i];
		}
		return sum;
	}

	public static double[] subtract(double[] a, double[] b) {
		int n = a.length;
		double[] diff = new double[n];
		for (int i = 0; i < n; i++) {
			diff[i] = a[i] - b[i];
		}
		return diff;
	}

	public static double[] scale(double[] a, double b) {
		int n = a.length;

		for (int i = 0; i < n; i++) {
			a[i] *= b;
		}
		return a;
	}

	public static double mod(double[] a) {
		int n = a.length;
		double mod = 0;
		for (int i = 0; i < n; i++) {
			mod += a[i] * a[i];
		}
		mod = Math.sqrt(mod);
		return mod;
	}
	

    public static double magnitude(double[] a) {
        return Math.sqrt((a[0] * a[0]) + (a[1] * a[1]) + (a[2] * a[2]));
    }

    public static double[] normalize(double[] a) {
        double mag = magnitude(a);
        if (mag == 0) 
        	{	
        		return  new double[] {0, 0, 0};
        	}
        	
        return new double[] {a[0] / mag, a[1] / mag, a[2] / mag};
    }
}
