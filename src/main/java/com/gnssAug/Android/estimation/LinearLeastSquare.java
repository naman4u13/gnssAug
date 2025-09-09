package com.gnssAug.Android.estimation;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.stream.IntStream;

import org.apache.commons.collections.set.ListOrderedSet;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.ejml.simple.SimpleMatrix;

import com.gnssAug.Android.constants.GnssDataConfig;
import com.gnssAug.Android.constants.Measurement;
import com.gnssAug.Android.models.Satellite;
import com.gnssAug.utility.LatLonUtil;
import com.gnssAug.utility.MathUtil;
import com.gnssAug.utility.Matrix;
import com.gnssAug.utility.SatUtil;
import com.gnssAug.utility.Weight;

public class LinearLeastSquare {
	private final static double SpeedofLight = 299792458;
	private static HashMap<Measurement, double[]> residualMap = new HashMap<Measurement, double[]>();
	private static HashMap<Measurement, Double> postVarOfUnitWMap = new HashMap<Measurement, Double>();
	private static HashMap<Measurement, HashMap<String, SimpleMatrix>> Cxx_hat_Map = new HashMap<Measurement, HashMap<String, SimpleMatrix>>();
	private static HashMap<Measurement, SimpleMatrix> Cyy_Map = new HashMap<Measurement, SimpleMatrix>();
	private static HashMap<Measurement, HashMap<String, SimpleMatrix>> Cxx_hat_updated_Map = new HashMap<Measurement, HashMap<String, SimpleMatrix>>();
	private static HashMap<Measurement, SimpleMatrix> Cyy_updated_Map = new HashMap<Measurement, SimpleMatrix>();

	private static double[] dop;
	private static HashMap<Measurement, ArrayList<Satellite>> testedSatListMap = new HashMap<Measurement, ArrayList<Satellite>>();
	

	public static double[] getEstPos(ArrayList<Satellite> satList, boolean isWLS, boolean useIGS) throws Exception {
		return process(satList, isWLS, false, false, false, Measurement.Pseudorange, null, useIGS);
	}

	public static double[] getEstPos(ArrayList<Satellite> satList, boolean isWLS, boolean doAnalyze, boolean doTest,
			boolean outlierAnalyze, boolean useIGS) throws Exception {
		return process(satList, isWLS, doAnalyze, doTest, outlierAnalyze, Measurement.Pseudorange, null, useIGS);
	}

	public static double[] getEstVel(ArrayList<Satellite> satList, boolean isWLS, double[] refPos, boolean useIGS)
			throws Exception {
		return process(satList, isWLS, false, false, false, Measurement.Doppler, refPos, useIGS);
	}

	public static double[] getEstVel(ArrayList<Satellite> satList, boolean isWLS, boolean doAnalyze, boolean doTest,
			boolean outlierAnalyze, double[] refPos, boolean useIGS) throws Exception {
		return process(satList, isWLS, doAnalyze, doTest, outlierAnalyze, Measurement.Doppler, refPos, useIGS);
	}

	public static double[] process(ArrayList<Satellite> satList, boolean isWLS, boolean doAnalyze, boolean doTest,
			boolean outlierAnalyze, Measurement type, double[] refPos, boolean useIGS) throws Exception {
		// Satellite count
		int n = satList.size();

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
		String[] obsvCodeList = useIGS ? (String[]) SatUtil.findObsvCodeArray(satList) : null;
		ListOrderedSet ssiSet = SatUtil.findSSIset(obsvCodeList);
		
		double[] estState = estimate(satList, weight, null, refPos, type, useIGS, ssiSet);
		if (doAnalyze) {

			if (outlierAnalyze) {
				analyze(weight, estState, satList, useAndroidW, type, refPos, useIGS, ssiSet);

			}
			double[][] _weight = weight;
			if (doTest) {

				int index = 0;
				int _n = n;
				while (index != -1) {

					index = qualityControl(_weight, estState, testedSatList, useAndroidW, type, refPos, useIGS,
							ssiSet);
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
						obsvCodeList = SatUtil.findObsvCodeArray(testedSatList);
						ssiSet = SatUtil.findSSIset(obsvCodeList);
						estState = estimate(testedSatList, _weight, null, refPos, type, useIGS, ssiSet);
					}
				}
			}
			if (!outlierAnalyze) {
				analyze(_weight, estState, testedSatList, useAndroidW, type, refPos, useIGS, ssiSet);

			}

		}
		if (outlierAnalyze) {
			testedSatListMap.put(type, satList);
		} else {
			testedSatListMap.put(type, testedSatList);
		}
		return estState;

	}

	private static void analyze(double[][] weight, double[] estState, ArrayList<Satellite> satList, boolean useAndroidW,
			Measurement type, double[] refPos, boolean useIGS, ListOrderedSet ssiSet) throws Exception {

		int n = satList.size();
		int m = 1;
		if (useIGS) {
			m = ssiSet.size();
		}
		int l = 3 + m;
		int[] constellCount = new int[m];
		double[] residual = new double[n];
		double[][] h = new double[n][l];
		double[] userEcef = estState;
		if (type == Measurement.Doppler) {
			userEcef = refPos;
		}
		for (int i = 0; i < n; i++) {
			Satellite sat = satList.get(i);
			char ssi = sat.getObsvCode().charAt(0);
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
					if (ssi == (char)ssiSet.get(j)) {
						h[i][3 + j] = 1;
						y_hat += estState[3 + j];
						constellCount[j] +=1; 
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
				priorVarOfUnitW = GnssDataConfig.pseudorange_priorVarOfUnitW;
			} else if (type == Measurement.Doppler) {
				priorVarOfUnitW = GnssDataConfig.doppler_priorVarOfUnitW;
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

		SimpleMatrix H = new SimpleMatrix(h);
		SimpleMatrix Ht = H.transpose();
		SimpleMatrix Cxx_hat= (Ht.mult(Cyy_inv).mult(H)).invert();
		
		SimpleMatrix P_H_perpendicular = Matrix.getPerpendicularProjection(H, Cyy_inv);
		SimpleMatrix Cee_hat = P_H_perpendicular.mult(Cyy).mult(P_H_perpendicular.transpose());
		SimpleMatrix redunMatrix = Cee_hat.mult(Cyy_inv);
		// Cyy updated is computed based on the theory of posteriori variance of unit
		// weight
		SimpleMatrix Cyy_updated = null;
		SimpleMatrix Cxx_hat_updated = null;
		if (n == l) {
			postVarOfUnitW = -1;
		}
		boolean flag = false;
		for(int i=0;i<m;i++)
		{
			if(constellCount[i]==1)
			{
				flag = true;
			}
		}
		if ((n <= l)||flag) {
			Cyy_updated = Cyy;
			Cxx_hat_updated = Cxx_hat;
			System.err.println("Redundancy for a particular measuement is zero, cannot calculate Cyy and Cxx hat");
		} else {
			Cyy_updated = new SimpleMatrix(n, n);
			double sum = 0;
			for (int i = 0; i < n; i++) {
				Cyy_updated.set(i, i, Math.pow(e_hat.get(i), 2) / redunMatrix.get(i, i));
				sum += redunMatrix.get(i, i);
			}
			Cxx_hat_updated = (Ht.mult(Cyy_updated.invert()).mult(H)).invert();
			
			if (Math.abs(sum - (n - (3 + m))) > 0.1) {
				throw new Exception("Error in LS redundancy computation");
			}
		}
//		if(type==Measurement.Doppler)
//		{
//			Cyy_updated = Cyy_updated.scale(25);
//		}
//		else
//		{
//			Cyy_updated = Cyy_updated.scale(4);
//		}

//		for (int i = 0; i < n; i++) {
//			residual[i] = e_hat.get(i) / Math.sqrt(Cee_hat.get(i, i));
//		}

		// Convert to ENU frame
		SimpleMatrix R = new SimpleMatrix(l, l);
		R.insertIntoThis(0, 0, new SimpleMatrix(LatLonUtil.getEcef2EnuRotMat(estState)));
		for (int i = 0; i < m; i++) {
			R.set(3 + i, 3 + i, 1);
		}
		Cxx_hat_Map.computeIfAbsent(type, k -> new HashMap<String, SimpleMatrix>()).put("ECEF", Cxx_hat);
		Cxx_hat_updated_Map.computeIfAbsent(type, k -> new HashMap<String, SimpleMatrix>()).put("ECEF",
				Cxx_hat_updated);
		Cxx_hat = R.mult(Cxx_hat).mult(R.transpose());
		
		Cxx_hat_updated = R.mult(Cxx_hat_updated).mult(R.transpose());
		
		Cxx_hat_Map.computeIfAbsent(type, k -> new HashMap<String, SimpleMatrix>()).put("ENU", Cxx_hat);
		Cxx_hat_updated_Map.computeIfAbsent(type, k -> new HashMap<String, SimpleMatrix>()).put("ENU", Cxx_hat_updated);
		SimpleMatrix _dop = R.mult((Ht.mult(H)).invert()).mult(R.transpose());
		dop = new double[] { _dop.get(0, 0), _dop.get(1, 1), _dop.get(2, 2), _dop.get(3, 3) };
		postVarOfUnitWMap.put(type, postVarOfUnitW);
		residualMap.put(type, residual);
		Cyy_Map.put(type, Cyy);
		Cyy_updated_Map.put(type, Cyy_updated);

	}

	// Baarda's Iterative Data Snooping
	private static int qualityControl(double[][] weight, double[] estState, ArrayList<Satellite> satList,
			boolean useAndroidW, Measurement type, double[] refPos, boolean useIGS,ListOrderedSet ssiSet)
			throws Exception {

		int n = satList.size();
		int m = 1;
		if (useIGS) {
			m = ssiSet.size();
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
			char ssi = sat.getObsvCode().charAt(0);
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
					if (ssi == (char)ssiSet.get(j)) {
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
				priorVarOfUnitW = GnssDataConfig.pseudorange_priorVarOfUnitW;
			} else if (type == Measurement.Doppler) {
				priorVarOfUnitW = GnssDataConfig.doppler_priorVarOfUnitW;
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
			double[] refPos, Measurement type, boolean useIGS, ListOrderedSet ssiSet) throws Exception {

		if (type == Measurement.Pseudorange) {
			return estimatePos(satList, weight, null, useIGS, ssiSet);
		} else if (type == Measurement.Doppler) {
			return estimateVel(satList, weight, refPos, useIGS, ssiSet);
		} else {
			throw new Exception("Fatal Error: Wrong Measurement Type chosen");
		}
	}

	private static double[] estimatePos(ArrayList<Satellite> satList, double[][] weight, HashSet<Integer> indexSet,
			boolean useIGS, ListOrderedSet ssiSet) throws Exception {

		boolean augBiasState = indexSet != null && !indexSet.isEmpty();
		int n = satList.size();
		int m = 1;
		if (useIGS) {
			m = ssiSet.size();
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
					char ssi = sat.getObsvCode().charAt(0);
					
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
							if (ssi == (char)ssiSet.get(j)) {
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
				SimpleMatrix HtWHinv= null;
				try {
				HtWHinv = (Ht.mult(W).mult(H)).invert();
				}
				catch (Exception e) {
					// TODO: handle exception
					System.out.println();
				}
				
				SimpleMatrix DeltaPR = new SimpleMatrix(deltaPR);
				SimpleMatrix DeltaX = HtWHinv.mult(Ht).mult(W).mult(DeltaPR);
				// updating Rx state vector, by adding deltaX vector
				IntStream.range(0, 3 + m).forEach(i -> estEcefClk[i] = estEcefClk[i] + DeltaX.get(i, 0));

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
			boolean useIGS, ListOrderedSet ssiSet) throws Exception {
		// Satellite count
		int n = satList.size();
		int m = 1;
		if (useIGS) {
			m = ssiSet.size();
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
				char ssi = sat.getObsvCode().charAt(0);
				// Its not really a ECI, therefore don't get confused
				double[] satEcef = sat.getSatEci();

				// Approx Geometric Range
				double approxGR = Math.sqrt(IntStream.range(0, 3).mapToDouble(j -> satEcef[j] - estEcefClk[j])
						.map(j -> Math.pow(j, 2)).reduce((j, k) -> j + k).getAsDouble());
				int index = i;
				IntStream.range(0, 3).forEach(j -> h[index][j] = -(satEcef[j] - estEcefClk[j]) / approxGR);
				if (useIGS) {
					for (int j = 0; j < m; j++) {
						if (ssi == (char)ssiSet.get(j)) {
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

	public static SimpleMatrix getCyy(Measurement type) {
		return Cyy_Map.get(type);
	}

	public static SimpleMatrix getCxx_hat(Measurement type, String frame) {
		return Cxx_hat_Map.get(type).get(frame);
	}

	public static SimpleMatrix getCyy_updated(Measurement type) {
		return Cyy_updated_Map.get(type);
	}

	public static SimpleMatrix getCxx_hat_updated(Measurement type, String frame) {
		return Cxx_hat_updated_Map.get(type).get(frame);
	}

	public static double[] getDop() {
		return dop;
	}

	public static ArrayList<Satellite> getTestedSatList(Measurement type) {
		return testedSatListMap.get(type);
	}

	

}
