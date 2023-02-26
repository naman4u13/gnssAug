package com.gnssAug.Android.estimation;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.stream.IntStream;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.ejml.simple.SimpleMatrix;

import com.gnssAug.Android.constants.Measurement;
import com.gnssAug.Android.models.Satellite;
import com.gnssAug.utility.LatLonUtil;
import com.gnssAug.utility.MathUtil;
import com.gnssAug.utility.Matrix;
import com.gnssAug.utility.Weight;

public class LinearLeastSquare {
	private final static double SpeedofLight = 299792458;
	private static HashMap<Measurement, double[]> residualMap = new HashMap<Measurement, double[]>();
	private static HashMap<Measurement, Double> postVarOfUnitWMap = new HashMap<Measurement, Double>();
	private static HashMap<Measurement, SimpleMatrix> Cxx_hat_Map = new HashMap<Measurement, SimpleMatrix>();

	private static double[] dop;
	private static HashMap<Measurement, ArrayList<Satellite>> testedSatListMap = new HashMap<Measurement, ArrayList<Satellite>>();
	final private static double pseudorange_priorVarOfUnitW = 401;
	final private static double doppler_priorVarOfUnitW = 4;

	public static double[] getEstPos(ArrayList<Satellite> satList, boolean isWLS, boolean useIGS) throws Exception {
		return process(satList, isWLS, false, false,false, Measurement.Pseudorange, null, useIGS);
	}

	public static double[] getEstPos(ArrayList<Satellite> satList, boolean isWLS, boolean doAnalyze, boolean doTest,boolean outlierAnalyze,
			boolean useIGS) throws Exception {
		return process(satList, isWLS, doAnalyze, doTest,outlierAnalyze, Measurement.Pseudorange, null, useIGS);
	}

	public static double[] getEstVel(ArrayList<Satellite> satList, boolean isWLS, double[] refPos, boolean useIGS)
			throws Exception {
		return process(satList, isWLS, false, false,false, Measurement.Doppler, refPos, useIGS);
	}

	public static double[] getEstVel(ArrayList<Satellite> satList, boolean isWLS, boolean doAnalyze, boolean doTest,boolean outlierAnalyze,
			double[] refPos, boolean useIGS) throws Exception {
		return process(satList, isWLS, doAnalyze, doTest,outlierAnalyze, Measurement.Doppler, refPos, useIGS);
	}
	
	public static double[] process(ArrayList<Satellite> satList, boolean isWLS,
			boolean doAnalyze, boolean doTest, boolean outlierAnalyze, Measurement type, double[] refPos,boolean useIGS)
			throws Exception {
		// Satellite count
		int n = satList.size();
		int DIA_type = 2;
		// Weight matrix
		double[][] weight = new double[n][n];
		boolean useAndroidW = false;
		for (int i = 0; i < n; i++) {
			satList.get(i).setOutlier(false);
		}
		ArrayList<Satellite> testedSatList = new ArrayList<Satellite>(satList);

		double scale = 1;
		if (type == Measurement.Doppler) {
			scale = 100;
		}
		/*
		 * If 'isWLS' flag is true, the estimation method is WLS and weight matrix will
		 * be based on elevation angle otherwise identity matrix will assigned for LS
		 */
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
			weight = Matrix.scale(weight, scale);
		} else {
			for (int i = 0; i < n; i++) {
				weight[i][i] = 1;
			}
		}
		String[] obsvCodeList = useIGS ? (String[]) findObsvCodeSet(satList) : null;
		
		double[] estState = estimate(satList, weight, null, refPos, type, useIGS, obsvCodeList);
		if (doAnalyze) {
			switch (DIA_type) {
			case 1:
				HashSet<Integer> indexSet = qualityControl(weight, estState, satList, useAndroidW, doTest, type, refPos,
						useIGS, obsvCodeList);

				if (!indexSet.isEmpty()) {
					testedSatList = new ArrayList<Satellite>();
					int j = 0;
					int _n = n - indexSet.size();
					double[][] _weight = new double[_n][_n];

					for (int i = 0; i < n; i++) {
						Satellite sat = satList.get(i);
						if (!indexSet.contains(i)) {
							sat.setOutlier(false);
							testedSatList.add(sat);
							_weight[j][j] = weight[i][i];
							j++;

						} else {
							sat.setOutlier(true);
						}
					}
					obsvCodeList = useIGS ? (String[]) findObsvCodeSet(testedSatList) : null;
					estState = estimate(testedSatList, _weight, null, refPos, type, useIGS, obsvCodeList);
					if (!outlierAnalyze) {
						qualityControl(_weight, estState, testedSatList, useAndroidW, false, type, refPos, useIGS,
								obsvCodeList);
					}
				}
				break;
			case 2:
				if (outlierAnalyze) {
					qualityControl(weight, estState, satList, useAndroidW, false, type, refPos,
							useIGS, obsvCodeList);

				}
				double[][] _weight = weight;
				if (doTest) {

					int index = 0;
					int _n = n;
					while (index != -1) {

						index = qualityControl2(_weight, estState, testedSatList, useAndroidW, type, refPos, useIGS,
								obsvCodeList);
						if (index != -1) {
							Satellite sat = testedSatList.remove(index);
							satList.get(satList.indexOf(sat)).setOutlier(true);
							_n = _n - 1;
							int j = 0;
							double[][] _tempweight = new double[_n][_n];
							for (int i = 0; i < _n + 1; i++) {
								if (i != index) {
									_tempweight[j][j] = _weight[i][i];
									j++;
								}
							}
							_weight = _tempweight;
							obsvCodeList = findObsvCodeSet(testedSatList);

							estState = estimate(testedSatList, _weight, null, refPos, type, useIGS, obsvCodeList);
						}
					}
				}
				if (!outlierAnalyze) {
					qualityControl(_weight, estState, testedSatList, useAndroidW, false, type, refPos, useIGS,
							obsvCodeList);
				}
				break;
			}
		}
		if (outlierAnalyze) {
			testedSatListMap.put(type, satList);
		} else {
			testedSatListMap.put(type, testedSatList);
		}
		return estState;

	}

	private static HashSet<Integer> qualityControl(double[][] weight, double[] estState, ArrayList<Satellite> satList,
			boolean useAndroidW, boolean doTest, Measurement type, double[] refPos, boolean useIGS,
			String[] obsvCodeList) throws Exception {
		HashSet<Integer> indexSet = new HashSet<Integer>();
		int n = satList.size();
		int m = 1;
		if (useIGS) {
			m = obsvCodeList.length;
		}
		int l = 3 + m;
		double[] residual = new double[n];
		double[][] h = new double[n][l];
		double[] userEcef = estState;
		if (type == Measurement.Doppler) {
			userEcef = refPos;
		}
		for (int i = 0; i < n; i++) {
			Satellite sat = satList.get(i);
			String obsvCode = sat.getObsvCode();
			// Its not really a ECI, therefore don't get confused
			double[] satEcef = sat.getSatEci();
			double gr_hat = MathUtil.getEuclidean(satEcef, userEcef);
			for (int j = 0; j < 3; j++) {
				h[i][j] = -(satEcef[j] - userEcef[j]) / gr_hat;
			}
			double y = 0;
			double y_hat = 0;
			if (useIGS) {
				for (int j = 0; j < m; j++) {
					if (obsvCode.equals(obsvCodeList[j])) {
						h[i][3 + j] = 1;
						y_hat += estState[3 + j];
					}
				}
			} else {
				h[i][3] = 1;
				y_hat += estState[3];
			}

			if (type == Measurement.Pseudorange) {
				y = sat.getPseudorange();
				// Approx Pseudorange
				y_hat += gr_hat;

			} else if (type == Measurement.Doppler) {
				double rangeRate = sat.getRangeRate();
				SimpleMatrix satVel = new SimpleMatrix(3, 1, true, sat.getSatVel());
				SimpleMatrix A = new SimpleMatrix(1, 3, true, new double[] { -h[i][0], -h[i][1], -h[i][2] });
				/*
				 * Observable derived from doppler and satellite velocity, refer Kaplan and
				 * personal notes
				 */
				double dopplerDerivedObs = rangeRate - A.mult(satVel).get(0);
				y = dopplerDerivedObs;
				// Approx DopplerDerivedObs
				y_hat += -(A.get(0) * estState[0]) - (A.get(1) * estState[1]) - (A.get(2) * estState[2]);

			}
			residual[i] = y - y_hat;
		}
		double priorVarOfUnitW = 1;
		SimpleMatrix Cyy = null;
		if (useAndroidW) {
			Cyy = new SimpleMatrix(weight).invert();

		} else {
			if (type == Measurement.Pseudorange) {
				priorVarOfUnitW = pseudorange_priorVarOfUnitW;
			} else if (type == Measurement.Doppler) {
				priorVarOfUnitW = doppler_priorVarOfUnitW;
			}
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
		double postVarOfUnitW = globalTq * priorVarOfUnitW / (n - l);
		if (n == l) {
			postVarOfUnitW = -1;
		}
		SimpleMatrix H = new SimpleMatrix(h);
		SimpleMatrix Ht = H.transpose();
		SimpleMatrix Cxx_hat = (Ht.mult(Cyy_inv).mult(H)).invert();

		if (doTest && n > (3 + m + 1)) {
			ChiSquaredDistribution csd = new ChiSquaredDistribution(n - l);
			double alpha = 0.01;
			if (globalTq == 0) {
				throw new Exception("Error: T stat is zero");
			}
			// Detection
			double globalPVal = 1 - csd.cumulativeProbability(globalTq);
			if (globalPVal < alpha) {
				double pVal_min = Double.MAX_VALUE;
				double fd_test_max = Double.MIN_VALUE;
				int len = n - (l + 1);
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
						SimpleMatrix P_C_ = null;

						P_C_ = Matrix.getProjection(C_, Cyy_inv);

						double Tq = Matrix.getNorm(P_C_.mult(e_hat), Cyy);
						double pVal = 1 - csd.cumulativeProbability(Tq);
						if (pVal < pVal_min) {
							pVal_min = pVal;
							fd_test_max = Tq / i;
							indexSet = new HashSet<Integer>(_indexSet);
						} else if (pVal == 0 && pVal_min == 0) {
							double fd_test = Tq / i;
							if (fd_test > fd_test_max) {
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
		SimpleMatrix R = new SimpleMatrix(l, l);
		R.insertIntoThis(0, 0, new SimpleMatrix(LatLonUtil.getEcef2EnuRotMat(estState)));
		for (int i = 0; i < m; i++) {
			R.set(3 + i, 3 + i, 1);
		}
		Cxx_hat = R.mult(Cxx_hat).mult(R.transpose());
		SimpleMatrix _dop = R.mult((Ht.mult(H)).invert()).mult(R.transpose());
		dop = new double[] { _dop.get(0, 0), _dop.get(1, 1), _dop.get(2, 2), _dop.get(3, 3) };
		Cxx_hat_Map.put(type, Cxx_hat);
		postVarOfUnitWMap.put(type, postVarOfUnitW);
		residualMap.put(type, residual);
		return indexSet;
	}

	// Baarda's Iterative Data Snooping
	private static int qualityControl2(double[][] weight, double[] estState, ArrayList<Satellite> satList,
			boolean useAndroidW, Measurement type, double[] refPos, boolean useIGS, String[] obsvCodeList)
			throws Exception {

		int n = satList.size();
		int m = 1;
		if (useIGS) {
			m = obsvCodeList.length;
		}
		int l = 3 + m;
		double[] residual = new double[n];
		double[][] h = new double[n][l];
		double[] userEcef = estState;
		if (type == Measurement.Doppler) {
			userEcef = refPos;
		}
		for (int i = 0; i < n; i++) {
			Satellite sat = satList.get(i);
			String obsvCode = sat.getObsvCode();
			// Its not really a ECI, therefore don't get confused
			double[] satEcef = sat.getSatEci();
			double gr_hat = MathUtil.getEuclidean(satEcef, userEcef);
			for (int j = 0; j < 3; j++) {
				h[i][j] = -(satEcef[j] - userEcef[j]) / gr_hat;
			}
			double y = 0;
			double y_hat = 0;
			if (useIGS) {
				for (int j = 0; j < m; j++) {
					if (obsvCode.equals(obsvCodeList[j])) {
						h[i][3 + j] = 1;
						y_hat += estState[3 + j];
					}
				}
			} else {
				h[i][3] = 1;
				y_hat += estState[3];
			}

			if (type == Measurement.Pseudorange) {
				y = sat.getPseudorange();
				// Approx Pseudorange
				y_hat += gr_hat;

			} else if (type == Measurement.Doppler) {
				double rangeRate = sat.getRangeRate();
				SimpleMatrix satVel = new SimpleMatrix(3, 1, true, sat.getSatVel());
				SimpleMatrix A = new SimpleMatrix(1, 3, true, new double[] { -h[i][0], -h[i][1], -h[i][2] });
				/*
				 * Observable derived from doppler and satellite velocity, refer Kaplan and
				 * personal notes
				 */
				double dopplerDerivedObs = rangeRate - A.mult(satVel).get(0);
				y = dopplerDerivedObs;
				// Approx DopplerDerivedObs
				y_hat += -(A.get(0) * estState[0]) - (A.get(1) * estState[1]) - (A.get(2) * estState[2]);

			}
			residual[i] = y - y_hat;
		}
		double priorVarOfUnitW = 1;
		SimpleMatrix Cyy = null;
		if (useAndroidW) {
			Cyy = new SimpleMatrix(weight).invert();

		} else {
			if (type == Measurement.Pseudorange) {
				priorVarOfUnitW = pseudorange_priorVarOfUnitW;
			} else if (type == Measurement.Doppler) {
				priorVarOfUnitW = doppler_priorVarOfUnitW;
			}
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
		int index = -1;
		if (n > (l + 1)) {
			ChiSquaredDistribution csd = new ChiSquaredDistribution(n - l);
			double alpha = 0.01;
			if (globalTq == 0) {
				throw new Exception("Error: T stat is zero");
			}
			// Detection
			double globalPVal = 1 - csd.cumulativeProbability(globalTq);
			if (globalPVal < alpha) {
				// Identification
				double max_w = Double.MIN_VALUE;
				SimpleMatrix H = new SimpleMatrix(h);
				SimpleMatrix P_H_perpendicular = Matrix.getPerpendicularProjection(H, Cyy_inv);
				SimpleMatrix Cee_hat = P_H_perpendicular.mult(Cyy).mult(P_H_perpendicular.transpose());
				for (int j = 0; j < n; j++) {
					double w = Math.abs(e_hat.get(j) / Math.sqrt(Cee_hat.get(j, j)));
					if (w > max_w) {
						max_w = w;
						index = j;
					}
				}

			}
		}

		return index;
	}

	private static double[] estimate(ArrayList<Satellite> satList, double[][] weight, HashSet<Integer> indexSet,
			double[] refPos, Measurement type, boolean useIGS, String[] obsvCodeList) throws Exception {

		if (type == Measurement.Pseudorange) {
			return estimatePos(satList, weight, null, useIGS, obsvCodeList);
		} else if (type == Measurement.Doppler) {
			return estimateVel(satList, weight, refPos, useIGS, obsvCodeList);
		} else {
			throw new Exception("Fatal Error: Wrong Measurement Type chosen");
		}
	}

	private static double[] estimatePos(ArrayList<Satellite> satList, double[][] weight, HashSet<Integer> indexSet,
			boolean useIGS, String[] obsvCodeList) throws Exception {

		boolean augBiasState = indexSet != null && !indexSet.isEmpty();
		int n = satList.size();
		int m = 1;
		if (useIGS) {
			m = obsvCodeList.length;
		}
		// variable to store estimated Rx position and clk offset
		double[] estEcefClk = new double[3 + m];
		/*
		 * Error variable based on norm value deltaX vector, intially assigned a big
		 * value
		 */
		double error = Double.MAX_VALUE;
		// Threshold to stop iteration or regression
		double threshold = 1e-3;
		// Minimum 3+m satellites are required to proceed
		if (n >= 3 + m) {
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
				double[][] h = new double[n][3 + m];
				// Iterate through each satellite, to compute LOS vector and Approx pseudorange
				for (int i = 0; i < n; i++) {

					Satellite sat = satList.get(i);
					String obsvCode = sat.getObsvCode();
					// Its not really a ECI, therefore don't get confused
					double[] satEcef = sat.getSatEci();
					double PR = sat.getPseudorange();
					// Approx Geometric Range
					double approxGR = Math.sqrt(IntStream.range(0, 3).mapToDouble(j -> satEcef[j] - estEcefClk[j])
							.map(j -> Math.pow(j, 2)).reduce((j, k) -> j + k).getAsDouble());
					// Approx Pseudorange Range
					double approxPR = approxGR;
					int index = i;
					IntStream.range(0, 3).forEach(j -> h[index][j] = -(satEcef[j] - estEcefClk[j]) / approxGR);
					if (useIGS) {
						for (int j = 0; j < m; j++) {
							if (obsvCode.equals(obsvCodeList[j])) {
								h[i][3 + j] = 1;
								approxPR += estEcefClk[3 + j];
							}
						}
					} else {
						h[i][3] = 1;
						approxPR += estEcefClk[3];
					}
					deltaPR[i][0] = PR - approxPR;

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
				IntStream.range(0, 3+m).forEach(i -> estEcefClk[i] = estEcefClk[i] + DeltaX.get(i, 0));
				

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

		throw new Exception("Satellite count is less than " + (3 + m) + ", can't compute user position");

	}

	private static double[] estimateVel(ArrayList<Satellite> satList, double[][] weight, double[] estEcefClk,
			boolean useIGS, String[] obsvCodeList) throws Exception {
		// Satellite count
		int n = satList.size();
		int m = 1;
		if (useIGS) {
			m = obsvCodeList.length;
		}
		double[] estVel = new double[3 + m];
		double[] z = new double[n];

		// Jacobian or Design Matrix
		double[][] h = new double[n][3 + m];
		// Minimum 4 satellite are required to proceed

		if (n >= 3 + m) {
			// Iterate through each satellite, to compute LOS vector and Approx pseudorange
			for (int i = 0; i < n; i++) {

				Satellite sat = satList.get(i);
				String obsvCode = sat.getObsvCode();
				// Its not really a ECI, therefore don't get confused
				double[] satEcef = sat.getSatEci();

				// Approx Geometric Range
				double approxGR = Math.sqrt(IntStream.range(0, 3).mapToDouble(j -> satEcef[j] - estEcefClk[j])
						.map(j -> Math.pow(j, 2)).reduce((j, k) -> j + k).getAsDouble());
				int index = i;
				IntStream.range(0, 3).forEach(j -> h[index][j] = -(satEcef[j] - estEcefClk[j]) / approxGR);
				if (useIGS) {
					for (int j = 0; j < m; j++) {
						if (obsvCode.equals(obsvCodeList[j])) {
							h[i][3 + j] = 1;
						}
					}
				} else {
					h[i][3] = 1;
				}
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
			IntStream.range(0, 3 + m).forEach(i -> estVel[i] = X.get(i, 0));
			return estVel;
		}

		throw new Exception("Satellite count is less than " + (3 + m) + " , can't compute user position");

	}

	public static double getPostVarOfUnitW(Measurement type) {
		return postVarOfUnitWMap.get(type);
	}

	public static double[] getResidual(Measurement type) {
		return residualMap.get(type);
	}

	public static SimpleMatrix getCxx_hat(Measurement type) {
		return Cxx_hat_Map.get(type);
	}

	public static double[] getDop() {
		return dop;
	}

	public static ArrayList<Satellite> getTestedSatList(Measurement type) {
		return testedSatListMap.get(type);
	}

	private static String[] findObsvCodeSet(ArrayList<Satellite> satList) {
		LinkedHashSet<String> obsvCodeSet = new LinkedHashSet<String>();
		for (int i = 0; i < satList.size(); i++) {
			obsvCodeSet.add(satList.get(i).getObsvCode());
		}
		return obsvCodeSet.toArray(new String[0]);

	}

}
