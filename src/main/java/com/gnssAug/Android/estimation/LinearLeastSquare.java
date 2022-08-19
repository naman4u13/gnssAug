package com.gnssAug.Android.estimation;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.stream.IntStream;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.ejml.dense.row.MatrixFeatures_DDRM;
import org.ejml.simple.SimpleMatrix;

import com.gnssAug.Android.models.Satellite;
import com.gnssAug.utility.LatLonUtil;
import com.gnssAug.utility.Weight;

public class LinearLeastSquare {
	private final static double SpeedofLight = 299792458;
	private static double[] residual = null;
	private static double postVarOfUnitW;
	private static SimpleMatrix Cxx_hat;
	private static double[] dop;
	private static ArrayList<Satellite> testedSatList;

	public static double[] process(ArrayList<Satellite> satList, boolean isWLS) throws Exception {
		return process(satList, isWLS, false, false);
	}

	public static double[] process(ArrayList<Satellite> satList, boolean isWLS, boolean doAnalyze, boolean doTest)
			throws Exception {
		// Satellite count
		int n = satList.size();
		// Weight matrix
		double[][] weight = new double[n][n];
		boolean useAndroidW = false;
		testedSatList = satList;
		/*
		 * If 'isWLS' flag is true, the estimation method is WLS and weight matrix will
		 * be based on elevation angle otherwise identity matrix will assigned for LS
		 */
		if (isWLS) {
			//
			if (useAndroidW) {
				for (int i = 0; i < n; i++) {

					weight[i][i] = 1
							/ Math.pow(satList.get(i).getReceivedSvTimeUncertaintyNanos() * SpeedofLight * 1e-9, 2);
				}
			} else {
				weight = Weight.computeCovInvMat2(satList);
			}

		} else {
			for (int i = 0; i < n; i++) {
				weight[i][i] = 1;
			}
		}
		double[] estEcefClk = regress(satList, weight);
		if (doAnalyze) {
			HashSet<Integer> indexSet = qualityControl(weight, estEcefClk, satList, useAndroidW, doTest);
			if (!indexSet.isEmpty()) {
				testedSatList = new ArrayList<Satellite>();
				int j = 0;
				int _n = n - indexSet.size();
				double[][] _weight = new double[_n][_n];

				for (int i = 0; i < n; i++) {
					Satellite sat = satList.get(i);
					if (!indexSet.contains(i)) {

						testedSatList.add(sat);
						_weight[j][j] = weight[i][i];
						j++;

					}
				}
				estEcefClk = regress(testedSatList, _weight);
				qualityControl(_weight, estEcefClk, testedSatList, useAndroidW, false);
			}
		}
		return estEcefClk;
	}

	public static double[] regress(ArrayList<Satellite> satList, double[][] weight) throws Exception {

		int n = satList.size();
		// variable to store estimated Rx position and clk offset
		double[] estEcefClk = new double[] { 0, 0, 0, 0 };
		/*
		 * Error variable based on norm value deltaX vector, intially assigned a big
		 * value
		 */
		double error = Double.MAX_VALUE;
		// Threshold to stop iteration or regression
		double threshold = 1e-3;
		// Minimum 4 satellite are required to proceed
		if (n >= 4) {

			while (error >= threshold) {

				// Misclosure vector
				double[][] deltaPR = new double[n][1];
				// Jacobian or Design Matrix
				double[][] h = new double[n][4];
				// Iterate through each satellite, to compute LOS vector and Approx pseudorange
				for (int i = 0; i < n; i++) {

					Satellite sat = satList.get(i);
					// Its not really a ECI, therefore don't get confused
					double[] satEcef = sat.getSatEci();
					double PR = sat.getPseudorange();
					// Approx Geometric Range
					double approxGR = Math.sqrt(IntStream.range(0, 3).mapToDouble(j -> satEcef[j] - estEcefClk[j])
							.map(j -> Math.pow(j, 2)).reduce((j, k) -> j + k).getAsDouble());
					// Approx Pseudorange Range
					double approxPR = approxGR + (SpeedofLight * estEcefClk[3]);
					deltaPR[i][0] = approxPR - PR;
					int index = i;
					IntStream.range(0, 3).forEach(j -> h[index][j] = (satEcef[j] - estEcefClk[j]) / approxGR);
					h[i][3] = 1;
				}
				// Least Squares implementation
				SimpleMatrix H = new SimpleMatrix(h);
				SimpleMatrix Ht = H.transpose();
				SimpleMatrix W = new SimpleMatrix(weight);
				SimpleMatrix HtWHinv = null;

				HtWHinv = (Ht.mult(W).mult(H)).invert();

				SimpleMatrix DeltaPR = new SimpleMatrix(deltaPR);
				SimpleMatrix DeltaX = HtWHinv.mult(Ht).mult(W).mult(DeltaPR);
				// updating Rx state vector, by adding deltaX vector
				IntStream.range(0, 3).forEach(i -> estEcefClk[i] = estEcefClk[i] + DeltaX.get(i, 0));
				estEcefClk[3] += (-DeltaX.get(3, 0)) / SpeedofLight;
				// Recomputing error - norm of deltaX vector
				error = Math.sqrt(IntStream.range(0, 3).mapToDouble(i -> Math.pow(DeltaX.get(i, 0), 2)).reduce(0,
						(i, j) -> i + j));

			}

			/*
			 * Regression is completed, error is below threshold, successfully estimated Rx
			 * Position and Clk Offset
			 */

			return estEcefClk;
		}

		throw new Exception("Satellite count is less than 4, can't compute user position");

	}

	public static HashSet<Integer> qualityControl(double[][] weight, double[] estEcefClk, ArrayList<Satellite> satList,
			boolean useAndroidW, boolean doTest) throws Exception {
		HashSet<Integer> indexSet = new HashSet<Integer>();
		int n = satList.size();
		residual = new double[n];
		double[][] h = new double[n][4];
		for (int i = 0; i < n; i++) {
			Satellite sat = satList.get(i);
			// Its not really a ECI, therefore don't get confused
			double[] satEcef = sat.getSatEci();
			double pr = sat.getPseudorange();

			// Approx Geometric Range
			double gr_hat = Math.sqrt(IntStream.range(0, 3).mapToDouble(j -> satEcef[j] - estEcefClk[j])
					.map(j -> Math.pow(j, 2)).reduce((j, k) -> j + k).getAsDouble());
			// Approx Pseudorange Range
			double pr_hat = gr_hat + (SpeedofLight * estEcefClk[3]);
			residual[i] = pr_hat - pr;
			int _i = i;
			IntStream.range(0, 3).forEach(j -> h[_i][j] = -(satEcef[j] - estEcefClk[j]) / gr_hat);
			h[i][3] = 1;
		}
		double priorVarOfUnitW = 1;
		SimpleMatrix Cyy = null;
		if (useAndroidW) {
			Cyy = new SimpleMatrix(weight).invert();
		} else {
			priorVarOfUnitW = 4;
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
		}

		SimpleMatrix e_hat = new SimpleMatrix(n, 1, true, residual);
		SimpleMatrix Cyy_inv = Cyy.invert();
		double detectT = e_hat.transpose().mult(Cyy_inv).mult(e_hat).get(0);
		postVarOfUnitW = detectT * priorVarOfUnitW / (n - 4);
		SimpleMatrix H = new SimpleMatrix(h);
		SimpleMatrix Ht = H.transpose();
		Cxx_hat = (Ht.mult(Cyy_inv).mult(H)).invert();

		boolean isSingleOut = true;
		if (doTest && n > 5) {
			ChiSquaredDistribution csd = new ChiSquaredDistribution(n - 4);
			double alpha = 1;
			if (detectT == 0) {
				throw new Exception("Error: T stat is zero");
			}
			if ((1 - csd.cumulativeProbability(detectT)) < alpha) {

				SimpleMatrix Cyy_hat = H.mult(Cxx_hat).mult(Ht);
				SimpleMatrix Cee_hat = Cyy.minus(Cyy_hat);
				if (!MatrixFeatures_DDRM.isDiagonalPositive(Cyy.getMatrix())) {

					throw new Exception("PositiveDiagonal Cyy test Failed");
				}
				if (!MatrixFeatures_DDRM.isDiagonalPositive(Cyy_hat.getMatrix())) {

					throw new Exception("PositiveDiagonal Cyy_hat test Failed");
				}
				if (!MatrixFeatures_DDRM.isDiagonalPositive(Cee_hat.getMatrix())) {

					throw new Exception("PositiveDiagonal Cee_hat test Failed");
				}

				ArrayList<Integer> indexList = new ArrayList<Integer>();
				ArrayList<Double> wList = new ArrayList<Double>();
				double alpha2 = 0.05;
				for (int i = 0; i < n; i++) {
					SimpleMatrix c = new SimpleMatrix(n, 1);
					c.set(i, 1);
					SimpleMatrix ct = c.transpose();
					double w = (ct.mult(Cyy_inv).mult(e_hat).get(0))
							/ (Math.sqrt(ct.mult(Cyy_inv).mult(Cee_hat).mult(Cyy_inv).mult(c).get(0)));
					NormalDistribution norm = new NormalDistribution();
					if (1 - norm.cumulativeProbability(Math.abs(w)) < (alpha2 / 2)) {
						indexList.add(i);
						if (isSingleOut) {
							wList.add(Math.abs(w));

						}
					}
				}
				if (isSingleOut) {
					double maxW = Double.MIN_VALUE;
					int index = -1;
					for (int i = 0; i < wList.size(); i++) {
						if (wList.get(i) > maxW) {
							index = indexList.get(i);
							maxW = wList.get(i);
						}
					}
					if (index != -1) {
						indexSet.add(index);
					}
				} else {

					double minT = Double.MAX_VALUE;
					int m = indexList.size();

					int len = Math.min(m, n - 4 - 1);
					// Normalized t stat with number of sat, this is a personal theorized stat,
					// unsure about its validity
					double normT = detectT / (n - 4);
					for (int i = 1; i <= len; i++) {
						Iterator<int[]> iterator = CombinatoricsUtils.combinationsIterator(m, i);
						while (iterator.hasNext()) {
							int[] combination = iterator.next();
							HashSet<Integer> _indexSet = new HashSet<Integer>();
							for (int j = 0; j < combination.length; j++) {
								_indexSet.add(indexList.get(combination[j]));
							}
							ArrayList<Satellite> _satList = new ArrayList<Satellite>();
							int k = 0;
							int _n = n - _indexSet.size();
							double[][] _Cyy_inv = new double[_n][_n];
							for (int j = 0; j < n; j++) {
								if (!_indexSet.contains(j)) {
									_satList.add(satList.get(j));
									_Cyy_inv[k][k] = Cyy_inv.get(j, j);
									k++;
								}
							}
							double t = 0;
							t = getNormT(_satList, _Cyy_inv);
							if (t < normT && t < minT) {
								if (t == 0) {
									throw new Exception("Error: T stat is zero");
								}
								minT = t;
								indexSet = new HashSet<Integer>(_indexSet);
							}
						}
					}

				}
				if (indexSet.isEmpty()) {
					System.err.println("Quality Control implementation failed: No Identification despite detection");
					// throw new Exception("Quality Control implementation failed");
				}
			}
		}
		// Convert to ENU frame
		SimpleMatrix R = new SimpleMatrix(4, 4);
		R.insertIntoThis(0, 0, new SimpleMatrix(LatLonUtil.getEcef2EnuRotMat(estEcefClk)));
		R.set(3, 3, 1);
		Cxx_hat = R.mult(Cxx_hat).mult(R.transpose());
		SimpleMatrix _dop = R.mult((Ht.mult(H)).invert()).mult(R.transpose());
		dop = new double[] { _dop.get(0, 0), _dop.get(1, 1), _dop.get(2, 2), _dop.get(3, 3) };
		return indexSet;
	}

	// Get sum of weighted squared residuals based test statistic
	public static double getNormT(ArrayList<Satellite> satList, double[][] _Cyy_inv) throws Exception {
		int n = satList.size();
		double[] res = new double[n];
		double[] estEcefClk = regress(satList, _Cyy_inv);
		for (int i = 0; i < n; i++) {
			Satellite sat = satList.get(i);
			// Its not really a ECI, therefore don't get confused
			double[] satEcef = sat.getSatEci();
			double pr = sat.getPseudorange();
			final double[] _estEcefClk = estEcefClk;
			// Approx Geometric Range
			double gr_hat = Math.sqrt(IntStream.range(0, 3).mapToDouble(j -> satEcef[j] - _estEcefClk[j])
					.map(j -> Math.pow(j, 2)).reduce((j, k) -> j + k).getAsDouble());
			// Approx Pseudorange Range
			double pr_hat = gr_hat + (SpeedofLight * estEcefClk[3]);
			res[i] = pr_hat - pr;

		}
		SimpleMatrix e_hat = new SimpleMatrix(n, 1, true, res);
		SimpleMatrix Cyy_inv = new SimpleMatrix(_Cyy_inv);
		double detectT = e_hat.transpose().mult(Cyy_inv).mult(e_hat).get(0);

		return (detectT / (n - 4));
	}

	public static double[] getEstVel(ArrayList<Satellite> satList, double[] estEcefClk) {
		// Satellite count
		int n = satList.size();
		// Weight matrix
		double[][] weight = new double[n][n];
		for (int i = 0; i < n; i++) {
			double elevAngle = satList.get(i).getElevAzm()[0];
			double CNo = satList.get(i).getCn0DbHz();
			double var = Math.pow(10, -(CNo / 10)) / Math.pow(Math.sin(elevAngle), 2);
			weight[i][i] = 1 / Math.pow(satList.get(i).getReceivedSvTimeUncertaintyNanos() * SpeedofLight * 1e-9, 2);
		}
		double[] estVel = new double[4];
		double[] z = new double[n];

		// Jacobian or Design Matrix
		double[][] h = new double[n][4];
		// Minimum 4 satellite are required to proceed

		if (n >= 4) {
			// Iterate through each satellite, to compute LOS vector and Approx pseudorange
			for (int i = 0; i < n; i++) {

				Satellite sat = satList.get(i);
				// Its not really a ECI, therefore don't get confused
				double[] satEcef = sat.getSatEci();

				// Approx Geometric Range
				double approxGR = Math.sqrt(IntStream.range(0, 3).mapToDouble(j -> satEcef[j] - estEcefClk[j])
						.map(j -> Math.pow(j, 2)).reduce((j, k) -> j + k).getAsDouble());
				int index = i;
				IntStream.range(0, 3).forEach(j -> h[index][j] = -(satEcef[j] - estEcefClk[j]) / approxGR);
				h[i][3] = 1;
				double rangeRate = sat.getRangeRate();
				SimpleMatrix satVel = new SimpleMatrix(3, 1, true, sat.getSatVel());
				SimpleMatrix A = new SimpleMatrix(1, 3, true, new double[] { -h[i][0], -h[i][1], -h[i][2] });

				/*
				 * Observable derived from doppler and satellite velocity, refer Kaplan and
				 * personal notes
				 */
				double dopplerDerivedObs = rangeRate - A.mult(satVel).get(0);
				z[i] = dopplerDerivedObs;

			}
			// Least Squares implementation
			SimpleMatrix H = new SimpleMatrix(h);
			SimpleMatrix Ht = H.transpose();
			SimpleMatrix W = new SimpleMatrix(weight);
			SimpleMatrix HtWHinv = (Ht.mult(W).mult(H)).invert();
			SimpleMatrix Z = new SimpleMatrix(n, 1, true, z);
			SimpleMatrix X = HtWHinv.mult(Ht).mult(W).mult(Z);

			// Velocity
			IntStream.range(0, 3).forEach(i -> estVel[i] = X.get(i, 0));
			estVel[3] = X.get(3, 0);
		}

		return estVel;
	}

	public static double getPostVarOfUnitW() {
		return postVarOfUnitW;
	}

	public static double[] getResidual() {
		return residual;
	}

	public static SimpleMatrix getCxx_hat() {
		return Cxx_hat;
	}

	public static double[] getDop() {
		return dop;
	}

	public static ArrayList<Satellite> getTestedSatList() {
		return testedSatList;
	}

}
