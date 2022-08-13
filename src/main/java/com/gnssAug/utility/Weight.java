package com.gnssAug.utility;

import java.util.ArrayList;
import java.util.stream.IntStream;

import com.gnssAug.Android.models.Satellite;

public class Weight {

	public static double[][] computeCovMat(ArrayList<Satellite> satList) {

		int SVcount = satList.size();
		double[][] covMat = new double[SVcount][SVcount];
		IntStream.range(0, SVcount).forEach(i -> covMat[i][i] = Weight.computeCoVariance(satList.get(i).getCn0DbHz(),
				satList.get(i).getElevAzm()[0]));
		// double[][] normWeight = Weight.normalize(covMat);
		return covMat;
	}

	public static double[][] computeCovMat2(ArrayList<com.gnssAug.Rinex.models.Satellite> satList) {

		int SVcount = satList.size();
		double[][] covMat = new double[SVcount][SVcount];
		IntStream.range(0, SVcount).forEach(i -> covMat[i][i] = 1
				/ Weight.computeCoVariance(satList.get(i).getCNo(), satList.get(i).getElevAzm()[0]));
		// double[][] normWeight = Weight.normalize(covMat);
		return covMat;
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
