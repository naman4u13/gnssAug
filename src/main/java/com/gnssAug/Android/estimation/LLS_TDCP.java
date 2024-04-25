package com.gnssAug.Android.estimation;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.stream.IntStream;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.ejml.simple.SimpleMatrix;

import com.gnssAug.Android.constants.GnssDataConfig;
import com.gnssAug.Android.constants.Measurement;
import com.gnssAug.Android.models.Satellite;
import com.gnssAug.Android.models.TDCP;
import com.gnssAug.utility.LatLonUtil;
import com.gnssAug.utility.MathUtil;
import com.gnssAug.utility.Matrix;
import com.gnssAug.utility.SatUtil;
import com.gnssAug.utility.Weight;

public class LLS_TDCP {

	private final static double SpeedofLight = 299792458;
	private static double[] residual = null;
	private static double postVarOfUnitW;
	private static HashMap<String, SimpleMatrix> Cxx_hat_Map = new HashMap<String, SimpleMatrix>();
	private static SimpleMatrix Cyy = null;
	private static HashMap<String, SimpleMatrix> Cxx_hat_updated = new HashMap<String, SimpleMatrix>();
	private static SimpleMatrix Cyy_updated = null;
	private static int commonSatCount;
	private static double[] dop;
	private static ArrayList<TDCP> _testedTdcpList = null;

	public static double[] getEstVel(ArrayList<Satellite> currentSatList, ArrayList<Satellite> prevSatList,
			boolean isWLS, double[] refPos, boolean useIGS) throws Exception {
		return process(currentSatList, prevSatList, isWLS, false, false, false, refPos, useIGS,false);
	}

	public static double[] getEstVel(ArrayList<Satellite> currentSatList, ArrayList<Satellite> prevSatList,
			boolean isWLS, boolean doAnalyze, boolean doTest, boolean outlierAnalyze, double[] refPos, boolean useIGS,boolean useGFCSD)
			throws Exception {
		return process(currentSatList, prevSatList, isWLS, doAnalyze, doTest, outlierAnalyze, refPos, useIGS,useGFCSD);
	}

	public static double[] process(ArrayList<Satellite> currentsatList, ArrayList<Satellite> prevSatList, boolean isWLS,
			boolean doAnalyze, boolean doTest, boolean outlierAnalyze, double[] refPos, boolean useIGS,boolean useGFCSD)
			throws Exception {

		ArrayList<TDCP> tdcpList = new ArrayList<TDCP>();

		// Satellite count
		int n_curr = currentsatList.size();
		int n_prev = prevSatList.size();
		for (int i = 0; i < n_curr; i++) {
			Satellite current_sat = currentsatList.get(i);
			String satID = current_sat.getObsvCode() + current_sat.getSvid();
			for (int j = 0; j < n_prev; j++) {
				Satellite prev_sat = prevSatList.get(j);
				String prev_satID = prev_sat.getObsvCode() + prev_sat.getSvid();
				if (satID.equals(prev_satID)) {
					SimpleMatrix satEci = new SimpleMatrix(3, 1, true, current_sat.getSatEci());
					SimpleMatrix prev_satEci = new SimpleMatrix(3, 1, true, prev_sat.getSatEci());
					SimpleMatrix unitLOS = new SimpleMatrix(1, 3, true, SatUtil.getUnitLOS(current_sat.getSatEci(), refPos));
					double ionoRate = current_sat.getIonoErr()-prev_sat.getIonoErr();
					double tropoRate = current_sat.getTropoErr()-prev_sat.getTropoErr();
					double dopplerDR = ((current_sat.getRangeRate()+prev_sat.getRangeRate())/2)-tropoRate+ionoRate;
					double phaseDR = current_sat.getPhase() - prev_sat.getPhase()+ionoRate;
					double satVelCorr = unitLOS.mult(satEci.minus(prev_satEci)).get(0);
					double wavelength = SpeedofLight/current_sat.getCarrierFrequencyHz();
					double approxCS = Math.abs(phaseDR-dopplerDR);
					if(useGFCSD)
					{
						if(approxCS<5*wavelength)
						{
							
							tdcpList.add(new TDCP(current_sat, phaseDR, satVelCorr, unitLOS,wavelength));
						}
						
					}
					else
					{
						tdcpList.add(new TDCP(current_sat, phaseDR, satVelCorr, unitLOS,wavelength));
					}
					
					
				}
			}
		}
		int n = tdcpList.size();
		commonSatCount = n;
		ArrayList<Satellite> satList = new ArrayList<Satellite>();
		for(int i=0;i<n;i++)
		{
			satList.add(tdcpList.get(i).getSat());
		}
		String[] obsvCodeList = useIGS ? (String[]) SatUtil.findObsvCodeArray(satList) : null;
		// Weight matrix
		double[][] weight = new double[n][n];
		boolean useAndroidW = false;
		for (int i = 0; i < n; i++) {
			tdcpList.get(i).setOutlier(false);
		}
		ArrayList<TDCP> testedTdcpList = new ArrayList<TDCP>(tdcpList);

		/*
		 * If 'isWLS' flag is true, the estimation method is WLS and weight matrix will
		 * be based on elevation angle otherwise identity matrix will assigned for LS
		 */
		/*
		 * If 'isWLS' flag is true, the estimation method is WLS and weight matrix will
		 * be based on elevation angle otherwise identity matrix will assigned for LS
		 */
		if (isWLS) {

			weight = Weight.computeCovInvMat2(satList);

		} else {
			for (int i = 0; i < n; i++) {
				weight[i][i] = 1;
			}
		}

		double[] estState = estimate(tdcpList, weight, null, refPos, useIGS, obsvCodeList);
		if (doAnalyze) {

			if (outlierAnalyze) {
				analyze(weight, estState, tdcpList, useAndroidW, refPos, useIGS, obsvCodeList);

			}
			double[][] _weight = weight;
			if (doTest) {

				int index = 0;
				int _n = n;
				while (index != -1) {

					index = qualityControl(_weight, estState, testedTdcpList, useAndroidW, refPos, useIGS,
							obsvCodeList);
					if (index != -1) {
						TDCP tdcp = testedTdcpList.remove(index);
						tdcpList.get(tdcpList.indexOf(tdcp)).setOutlier(true);
						
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
						satList.clear();
						for(int i=0;i<testedTdcpList.size();i++)
						{
							satList.add(testedTdcpList.get(i).getSat());
						}
						obsvCodeList = SatUtil.findObsvCodeArray(satList);

						estState = estimate(testedTdcpList, _weight, null, refPos, useIGS, obsvCodeList);
					}
				}
			}
			if (!outlierAnalyze) {
				analyze(_weight, estState, testedTdcpList, useAndroidW, refPos, useIGS, obsvCodeList);

			}

		}
		if (outlierAnalyze) {
			_testedTdcpList = new ArrayList<TDCP>(tdcpList);
		} else {
			_testedTdcpList = new ArrayList<TDCP>(testedTdcpList);
		}
		return estState;

	}

	private static void analyze(double[][] weight, double[] estState, ArrayList<TDCP> tdcpList, boolean useAndroidW,
			double[] refPos, boolean useIGS, String[] obsvCodeList) throws Exception {

		int n = tdcpList.size();
		int m = 1;
		if (useIGS) {
			m = obsvCodeList.length;
		}
		int l = 3 + m;
		int[] constellCount = new int[m];
		residual = new double[n];
		double[][] h = new double[n][l];
		double[] userEcef = refPos;

		for (int i = 0; i < n; i++) {
			TDCP tdcp = tdcpList.get(i);
			Satellite sat = tdcp.getSat();
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
						constellCount[j] += 1;
					}
				}
			} else {
				h[i][3] = 1;
				y_hat += estState[3];
			}

			double deltaRange = tdcp.getDeltaRange();
			SimpleMatrix A = new SimpleMatrix(1, 3, true, new double[] { -h[i][0], -h[i][1], -h[i][2] });
			y = deltaRange-tdcp.getSatVelCorr();
			// Approx DopplerDerivedObs
			y_hat += -(A.get(0) * estState[0]) - (A.get(1) * estState[1]) - (A.get(2) * estState[2]);
			residual[i] = y - y_hat;
		}
		double priorVarOfUnitW = 1;
		if (useAndroidW) {
			Cyy = new SimpleMatrix(weight).invert();

		} else {
			
			priorVarOfUnitW = GnssDataConfig.tdcp_priorVarOfUnitW;
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
		postVarOfUnitW = globalTq * priorVarOfUnitW / (n - l);

		SimpleMatrix H = new SimpleMatrix(h);
		SimpleMatrix Ht = H.transpose();
		SimpleMatrix Cxx_hat = (Ht.mult(Cyy_inv).mult(H)).invert();

		SimpleMatrix P_H_perpendicular = Matrix.getPerpendicularProjection(H, Cyy_inv);
		SimpleMatrix Cee_hat = P_H_perpendicular.mult(Cyy).mult(P_H_perpendicular.transpose());
		SimpleMatrix redunMatrix = Cee_hat.mult(Cyy_inv);
		
		if (n == l) {
			postVarOfUnitW = -1;
		}
		boolean flag = false;
		for (int i = 0; i < m; i++) {
			if (constellCount[i] == 1) {
				flag = true;
			}
		}

		// Convert to ENU frame
		SimpleMatrix R = new SimpleMatrix(l, l);
		R.insertIntoThis(0, 0, new SimpleMatrix(LatLonUtil.getEcef2EnuRotMat(estState)));
		for (int i = 0; i < m; i++) {
			R.set(3 + i, 3 + i, 1);
		}
		Cxx_hat_Map.put("ECEF", Cxx_hat);
		
		Cxx_hat = R.mult(Cxx_hat).mult(R.transpose());
		Cxx_hat_Map.put("ENU", Cxx_hat);
		
		SimpleMatrix _dop = R.mult((Ht.mult(H)).invert()).mult(R.transpose());
		dop = new double[] { _dop.get(0, 0), _dop.get(1, 1), _dop.get(2, 2), _dop.get(3, 3) };
		for (int i = 0; i < n; i++) {
			TDCP tdcp = tdcpList.get(i);
			residual[i] = residual[i]/tdcp.getWavelength();
		}
		
	}

	// Baarda's Iterative Data Snooping
	private static int qualityControl(double[][] weight, double[] estState, ArrayList<TDCP> tdcpList,
			boolean useAndroidW, double[] refPos, boolean useIGS, String[] obsvCodeList)
			throws Exception {

		int n = tdcpList.size();
		int m = 1;
		if (useIGS) {
			m = obsvCodeList.length;
		}
		int l = 3 + m;
		int[] constellCount = new int[m];
		residual = new double[n];
		double[][] h = new double[n][l];
		double[] userEcef = refPos;

		for (int i = 0; i < n; i++) {
			TDCP tdcp = tdcpList.get(i);
			Satellite sat = tdcp.getSat();
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
						constellCount[j] += 1;
					}
				}
			} else {
				h[i][3] = 1;
				y_hat += estState[3];
			}

			double deltaRange = tdcp.getDeltaRange();
			SimpleMatrix A = new SimpleMatrix(1, 3, true, new double[] { -h[i][0], -h[i][1], -h[i][2] });
			y = deltaRange-tdcp.getSatVelCorr();
			// Approx DopplerDerivedObs
			y_hat += -(A.get(0) * estState[0]) - (A.get(1) * estState[1]) - (A.get(2) * estState[2]);
			residual[i] = y - y_hat;
		}
		double priorVarOfUnitW = 1;
		if (useAndroidW) {
			Cyy = new SimpleMatrix(weight).invert();

		} else {
			
			priorVarOfUnitW = GnssDataConfig.tdcp_priorVarOfUnitW;
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

	private static double[] estimate(ArrayList<TDCP> tdcpList, double[][] weight, HashSet<Integer> indexSet,
			double[] refPos, boolean useIGS, String[] obsvCodeList) throws Exception {

		return estimateVel(tdcpList, weight, refPos, useIGS, obsvCodeList);

	}

	private static double[] estimateVel(ArrayList<TDCP> tdcpList, double[][] weight, double[] estEcefClk,
			boolean useIGS, String[] obsvCodeList) throws Exception {
		// Satellite count
		int n = tdcpList.size();
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

				TDCP tdcp = tdcpList.get(i);
				Satellite sat = tdcp.getSat();
				String obsvCode = sat.getObsvCode();
				SimpleMatrix unitLOS = tdcp.getUnitLOS();
				SimpleMatrix negUnitLOS = unitLOS.scale(-1);
				int index = i;
				IntStream.range(0, 3).forEach(j -> h[index][j] = negUnitLOS.get(j));
				if (useIGS) {
					for (int j = 0; j < m; j++) {
						if (obsvCode.equals(obsvCodeList[j])) {
							h[i][3 + j] = 1;
						}
					}
				} else {
					h[i][3] = 1;
				}
				double deltaRange = tdcp.getDeltaRange();
				z[i] = deltaRange - tdcp.getSatVelCorr();

			}
			// Least Squares implementation
			SimpleMatrix H = new SimpleMatrix(h);
			SimpleMatrix Ht = H.transpose();
			SimpleMatrix W = new SimpleMatrix(weight);
			SimpleMatrix HtWHinv= (Ht.mult(W).mult(H)).invert();
			
			SimpleMatrix Z = new SimpleMatrix(n, 1, true, z);
			SimpleMatrix X = HtWHinv.mult(Ht).mult(W).mult(Z);

			// Velocity
			IntStream.range(0, 3 + m).forEach(i -> estVel[i] = X.get(i, 0));
			return estVel;
		}

		throw new Exception("Satellite count is less than " + (3 + m) + " , can't compute user position");

	}

	public static double getPostVarOfUnitW() {
		return postVarOfUnitW;
	}

	public static double[] getResidual() {
		return residual;
	}

	public static SimpleMatrix getCyy() {
		return Cyy;
	}

	public static SimpleMatrix getCxx_hat(String frame) {
		return Cxx_hat_Map.get(frame);
	}

	public static int getCommSatCount()
	{
		return commonSatCount;
	}
	public static double[] getDop() {
		return dop;
	}

	public static ArrayList<TDCP> getTestedTdcpList() {
		return _testedTdcpList;
	}
}
