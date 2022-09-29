package com.gnssAug.Android.estimation;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.stream.IntStream;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.ejml.simple.SimpleMatrix;

import com.gnssAug.Android.models.Satellite;
import com.gnssAug.utility.LatLonUtil;
import com.gnssAug.utility.Matrix;
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
		boolean adaptAugState = false;
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
				double[][] _weight = weight;
				testedSatList = satList;
				if (adaptAugState) {
					estEcefClk = regress(satList, weight, indexSet);
				} else {
					testedSatList = new ArrayList<Satellite>();
					int j = 0;
					int _n = n - indexSet.size();
					_weight = new double[_n][_n];
					for (int i = 0; i < n; i++) {
						Satellite sat = satList.get(i);
						if (!indexSet.contains(i)) {
							testedSatList.add(sat);
							_weight[j][j] = weight[i][i];
							j++;
						}
					}
					estEcefClk = regress(testedSatList, _weight);

				}
				qualityControl(_weight, estEcefClk, testedSatList, useAndroidW, false);
			}
		}
		return estEcefClk;
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
		double priorVarOfUnitW = 7;
		SimpleMatrix Cyy = null;
		if (useAndroidW) {
			Cyy = new SimpleMatrix(weight).invert();

		} else {
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
		double globalTq = e_hat.transpose().mult(Cyy_inv).mult(e_hat).get(0);
		postVarOfUnitW = globalTq * priorVarOfUnitW / (n - 4);
		if (n == 4) {
			postVarOfUnitW = -1;
		}
		SimpleMatrix H = new SimpleMatrix(h);
		SimpleMatrix Ht = H.transpose();
		Cxx_hat = (Ht.mult(Cyy_inv).mult(H)).invert();
		if (doTest && n > 5) {
			ChiSquaredDistribution csd = new ChiSquaredDistribution(n - 4);
			double alpha = 0.05;
			if (globalTq == 0) {
				throw new Exception("Error: T stat is zero");
			}
			// Detection
			double globalPVal = 1 - csd.cumulativeProbability(globalTq);
			if (globalPVal < alpha) {
				double pVal_min = Double.MAX_VALUE;
				double fd_test_max = Double.MIN_VALUE;
				int len = n - 5;
				for (int i = 1; i <= len; i++) {
					Iterator<int[]> iterator = CombinatoricsUtils.combinationsIterator(n, i);
					csd = new ChiSquaredDistribution(i);
					while (iterator.hasNext()) {
						HashSet<Integer> _indexSet = new HashSet<Integer>();
						int[] combination = iterator.next();
						SimpleMatrix C = new SimpleMatrix(n, i);
						for (int j = 0; j < i; j++) {
							C.set(combination[j], j, 1);
							_indexSet.add(combination[j]);
						}
						SimpleMatrix P_H_perpendicular = Matrix.getPerpendicularProjection(H, Cyy_inv);
						SimpleMatrix C_ = P_H_perpendicular.mult(C);
						SimpleMatrix P_C_ = Matrix.getProjection(C_, Cyy_inv);
						double Tq = Matrix.getNorm(P_C_.mult(e_hat), Cyy);
						double pVal = 1 - csd.cumulativeProbability(Tq);
						if (pVal < pVal_min) {
//							if (_indexSet.size() > indexSet.size()) {
//								System.out.print("");
//							}
							pVal_min = pVal;
							fd_test_max = Tq / i;
							indexSet = new HashSet<Integer>(_indexSet);
						} else if (pVal == 0 && pVal_min == 0) {
							double fd_test = Tq / i;
							if (fd_test > fd_test_max) {
//								if (_indexSet.size() > indexSet.size()) {
//									System.out.print("");
//								}
								fd_test_max = fd_test;
								indexSet = new HashSet<Integer>(_indexSet);
							}
						}

					}
				}
				if (indexSet.isEmpty()) {
					throw new Exception("IndexSet cannot be empty: Impossible to have detection but no identification");
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

	public static double[] regress(ArrayList<Satellite> satList, double[][] weight) throws Exception {

		return regress(satList, weight, null);

	}

	public static double[] regress(ArrayList<Satellite> satList, double[][] weight, HashSet<Integer> indexSet)
			throws Exception {

		boolean augBiasState = indexSet != null && !indexSet.isEmpty();
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
			SimpleMatrix C = new SimpleMatrix(0, 0);
			if (augBiasState) {
				ArrayList<Integer> indexes = new ArrayList<Integer>(indexSet);
				Collections.sort(indexes);
				int q = indexes.size();
				C = new SimpleMatrix(n, q);
				for (int j = 0; j < q; j++) {
					int i = indexes.get(j);
					C.set(i, j, 1);
				}
			}
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
					deltaPR[i][0] = PR - approxPR;
					int index = i;
					IntStream.range(0, 3).forEach(j -> h[index][j] = -(satEcef[j] - estEcefClk[j]) / approxGR);
					h[i][3] = 1;
				}
				// Least Squares implementation
				SimpleMatrix H = new SimpleMatrix(h);
				H = H.concatColumns(C);
				SimpleMatrix Ht = H.transpose();
				SimpleMatrix W = new SimpleMatrix(weight);
				SimpleMatrix HtWHinv = null;

				HtWHinv = (Ht.mult(W).mult(H)).invert();

				SimpleMatrix DeltaPR = new SimpleMatrix(deltaPR);
				SimpleMatrix DeltaX = HtWHinv.mult(Ht).mult(W).mult(DeltaPR);
				// updating Rx state vector, by adding deltaX vector
				IntStream.range(0, 3).forEach(i -> estEcefClk[i] = estEcefClk[i] + DeltaX.get(i, 0));
				estEcefClk[3] += DeltaX.get(3, 0) / SpeedofLight;
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
