package com.gnssAug.utility;

import java.util.ArrayList;
import java.util.stream.IntStream;

import org.ejml.simple.SimpleMatrix;

import com.gnssAug.Android.models.Satellite;

public class Weight {

	public static double[][] computeCovInvMat2(ArrayList<Satellite> satList) {

		int SVcount = satList.size();
		double[][] covInvMat = new double[SVcount][SVcount];
		IntStream.range(0, SVcount).forEach(i -> covInvMat[i][i] = 1
				/ Weight.computeCoVariance(satList.get(i).getCn0DbHz(), satList.get(i).getElevAzm()[0]));
//		IntStream.range(0, SVcount)
//				.forEach(i -> covInvMat[i][i] = 1 / Math.pow(Math.sin(satList.get(i).getElevAzm()[0]), 2));

		return covInvMat;
	}

	public static double[][] computeCovInvMat(ArrayList<com.gnssAug.Rinex.models.Satellite> satList) {

		int SVcount = satList.size();
		double[][] covInvMat = new double[SVcount][SVcount];
		IntStream.range(0, SVcount).forEach(i -> covInvMat[i][i] = 1
				/ Weight.computeCoVariance(satList.get(i).getCNo(), satList.get(i).getElevAzm()[0]));
//		IntStream.range(0, SVcount)
//				.forEach(i -> covInvMat[i][i] = 1 / Math.pow(Math.sin(satList.get(i).getElevAzm()[0]), 2));

		return covInvMat;
	}

	public static SimpleMatrix getNormCyy(ArrayList<Satellite> satList, double priorVarOfUnitW) {
		int n = satList.size();
		double[][] weight = Weight.computeCovInvMat2(satList);
		SimpleMatrix Cyy = null;

		double[][] cov = new double[n][n];
		double max = Double.MIN_VALUE;
		for (int i = 0; i < n; i++) {
			max = Math.max(weight[i][i], max);
		}
		for (int i = 0; i < n; i++) {
			weight[i][i] = weight[i][i] / max;
		}

		for (int i = 0; i < n; i++) {
			cov[i][i] = priorVarOfUnitW / weight[i][i];
		}
		Cyy = new SimpleMatrix(cov);
		return Cyy;
	}

	public static double[][] normalize(double[][] W) {

		int n = W.length;
		double rms = Math
				.sqrt(IntStream.range(0, n).mapToDouble(i -> W[i][i] * W[i][i]).reduce(0, (i, j) -> i + j) / n);

		IntStream.range(0, n).forEach(i -> W[i][i] = W[i][i] / rms);
		return W;
	}

	public static double computeCoVariance(double CNo, double ElevAng) {
		double var = Math.pow(10, -(CNo / 10)) / Math.pow(Math.sin(ElevAng), 2);
		return var;
	}

}
