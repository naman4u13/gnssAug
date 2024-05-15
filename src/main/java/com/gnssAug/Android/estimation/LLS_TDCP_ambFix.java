package com.gnssAug.Android.estimation;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.stream.IntStream;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.ejml.simple.SimpleMatrix;
import com.gnssAug.Android.constants.GnssDataConfig;
import com.gnssAug.Android.models.CycleSlipDetect;
import com.gnssAug.Android.models.Satellite;
import com.gnssAug.helper.lambda.Lambda;
import com.gnssAug.utility.LatLonUtil;
import com.gnssAug.utility.MathUtil;
import com.gnssAug.utility.Matrix;
import com.gnssAug.utility.SatUtil;
import com.gnssAug.utility.Weight;

public class LLS_TDCP_ambFix {

	private final static double SpeedofLight = 299792458;
	private static double[] residual = null;
	private static double postVarOfUnitW;
	private static HashMap<String, SimpleMatrix> Cxx_hat_Map = new HashMap<String, SimpleMatrix>();
	private static SimpleMatrix Cyy = null;
	private static int commonSatCount;
	private static double[] dop;
	private static ArrayList<CycleSlipDetect> testedCsdList = null;
	private static long ambDetectedCount = 0;
	private static long ambRepairedCount = 0;
	protected static TreeMap<Long, Integer> ambDetectedCountMap = new TreeMap<Long, Integer>();
	protected static TreeMap<Long, Integer> ambRepairedCountMap = new TreeMap<Long, Integer>();

	public static double[] getEstVel(ArrayList<Satellite> currentSatList, ArrayList<Satellite> prevSatList,
			boolean isWLS, double[] refPos, boolean useIGS,long currentTime) throws Exception {
		return process(currentSatList, prevSatList, isWLS, false, refPos, useIGS,currentTime);
	}

	public static double[] getEstVel(ArrayList<Satellite> currentSatList, ArrayList<Satellite> prevSatList,
			boolean isWLS, boolean doAnalyze, double[] refPos, boolean useIGS,long currentTime) throws Exception {
		return process(currentSatList, prevSatList, isWLS, doAnalyze, refPos, useIGS,currentTime);
	}

	public static double[] process(ArrayList<Satellite> currentsatList, ArrayList<Satellite> prevSatList, boolean isWLS,
			boolean doAnalyze, double[] refPos, boolean useIGS,long currentTime) throws Exception {

		ArrayList<CycleSlipDetect> csdList = new ArrayList<CycleSlipDetect>();

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
					SimpleMatrix unitLOS = new SimpleMatrix(1, 3, true,
							SatUtil.getUnitLOS(current_sat.getSatEci(), refPos));
					double ionoRate = current_sat.getIonoErr() - prev_sat.getIonoErr();
					double tropoRate = current_sat.getTropoErr() - prev_sat.getTropoErr();
					double dopplerDR = ((current_sat.getRangeRate() + prev_sat.getRangeRate()) / 2) - tropoRate
							+ ionoRate;
					double phaseDR = current_sat.getPhase() - prev_sat.getPhase() + ionoRate;
					double satVelCorr = unitLOS.mult(satEci.minus(prev_satEci)).get(0);
					double wavelength = SpeedofLight / current_sat.getCarrierFrequencyHz();
					double approxCS = Math.abs(phaseDR - dopplerDR);

					if (approxCS < 5 * wavelength) {
						csdList.add(new CycleSlipDetect(current_sat, dopplerDR, phaseDR, ionoRate, false, wavelength,
								satVelCorr, unitLOS,currentTime));
					} else if (approxCS < 100 * wavelength) {
						csdList.add(new CycleSlipDetect(current_sat, dopplerDR, phaseDR, ionoRate, true, wavelength,
								satVelCorr, unitLOS,currentTime));
					}
				}
			}
		}
		int n = csdList.size();
		commonSatCount = n;
		ArrayList<Satellite> testedSatList = new ArrayList<Satellite>();
		ArrayList<CycleSlipDetect> testedCsdList = new ArrayList<CycleSlipDetect>();
		for (int i = 0; i < n; i++) {
			if (!csdList.get(i).isCS()) {
				testedSatList.add(csdList.get(i).getSat());
				testedCsdList.add(csdList.get(i));
			}
		}
		n = testedSatList.size();
		// Weight matrix
		double[][] weight = new double[n][n];
		boolean useAndroidW = false;

		/*
		 * If 'isWLS' flag is true, the estimation method is WLS and weight matrix will
		 * be based on elevation angle otherwise identity matrix will assigned for LS
		 */
		/*
		 * If 'isWLS' flag is true, the estimation method is WLS and weight matrix will
		 * be based on elevation angle otherwise identity matrix will assigned for LS
		 */
		if (isWLS) {

			weight = Weight.computeCovInvMat2(testedSatList);

		} else {
			for (int i = 0; i < n; i++) {
				weight[i][i] = 1;
			}
		}

		double[] estState = null;

		int index = 0;

		int ctr = 0;
		while (index != -1) {

			estState = estimateVel(testedCsdList, new SimpleMatrix(weight), refPos, useIGS);
			index = qualityControl(weight, estState, testedCsdList, useAndroidW, refPos, useIGS);
			if (index != -1) {
				CycleSlipDetect csd = testedCsdList.remove(index);
				csdList.get(csdList.indexOf(csd)).setCS(true);
				n = n - 1;
				int j = 0;
				double[][] _tempweight = new double[n][n];
				for (int i = 0; i < n + 1; i++) {
					if (i != index) {
						_tempweight[j][j] = weight[i][i];
						j++;
					}
				}
				weight = _tempweight;
				ctr++;
			}
		}
		estState = estimateVel(csdList, isWLS, refPos, useIGS, doAnalyze,currentTime);

		return estState;

	}

	// Baarda's Iterative Data Snooping
	private static int qualityControl(double[][] weight, double[] estState, ArrayList<CycleSlipDetect> csdList,
			boolean useAndroidW, double[] refPos, boolean useIGS) throws Exception {

		int n = csdList.size();
		ArrayList<Satellite> satList = new ArrayList<Satellite>();
		for (int i = 0; i < n; i++) {

			satList.add(csdList.get(i).getSat());

		}
		String[] obsvCodeList = useIGS ? (String[]) SatUtil.findObsvCodeArray(satList) : null;
		int m = obsvCodeList.length;
		int l = 3 + m;
		int[] constellCount = new int[m];
		residual = new double[n];
		double[][] h = new double[n][l];
		double[] userEcef = refPos;

		for (int i = 0; i < n; i++) {
			CycleSlipDetect csd = csdList.get(i);
			Satellite sat = csd.getSat();
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

			double phaseDR = csd.getCarrierPhaseDR();
			SimpleMatrix A = new SimpleMatrix(1, 3, true, new double[] { -h[i][0], -h[i][1], -h[i][2] });
			y = phaseDR - csd.getSatVelCorr();
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

	private static double[] estimateVel(ArrayList<CycleSlipDetect> csdList, boolean isWLS, double[] estEcefClk,
			boolean useIGS, boolean doAnalyze,long currentTime) throws Exception {
		// Satellite count
		int n = csdList.size();
		int ambCount = 0;
		ArrayList<Satellite> satList = new ArrayList<Satellite>();
		for (int i = 0; i < n; i++) {
			satList.add(csdList.get(i).getSat());
			if (csdList.get(i).isCS()) {
				ambCount++;
			}
		}
		ambDetectedCountMap.put(currentTime,ambCount);
		ambDetectedCount += ambCount;
		String[] obsvCodeList = useIGS ? (String[]) SatUtil.findObsvCodeArray(satList) : null;
		int m = 1;
		if (useIGS) {
			m = obsvCodeList.length;
		}
		double[] estVel = new double[3 + m];
		SimpleMatrix z = new SimpleMatrix(2 * n, 1);

		// Jacobian or Design Matrix
		SimpleMatrix H = new SimpleMatrix(2 * n, 3 + m + ambCount);
		SimpleMatrix unitLOS = new SimpleMatrix(SatUtil.getUnitLOS(satList, estEcefClk));
		H.insertIntoThis(0, 0, unitLOS.scale(-1));
		H.insertIntoThis(n, 0, unitLOS.scale(-1));

		SimpleMatrix W = new SimpleMatrix(2 * n, 2 * n);
		if (isWLS) {

			SimpleMatrix doppler_Cyy = Weight.getNormCyy(satList, GnssDataConfig.doppler_priorVarOfUnitW);
			SimpleMatrix tdcp_Cyy = Weight.getNormCyy(satList, GnssDataConfig.tdcp_priorVarOfUnitW);
			W.insertIntoThis(0, 0, tdcp_Cyy.invert());
			W.insertIntoThis(n, n, doppler_Cyy.invert());

		} else {

			SimpleMatrix doppler_Cyy = SimpleMatrix.identity(n).scale(GnssDataConfig.doppler_priorVarOfUnitW);
			SimpleMatrix tdcp_Cyy = SimpleMatrix.identity(n).scale(GnssDataConfig.tdcp_priorVarOfUnitW);
			W.insertIntoThis(0, 0, tdcp_Cyy.invert());
			W.insertIntoThis(n, n, doppler_Cyy.invert());
		}

		int ctr = 3 + m;
		// Minimum 4 satellite are required to proceed
		if ((2 * n) >= 3 + m + ambCount) {
			// Iterate through each satellite, to compute LOS vector and Approx pseudorange
			for (int i = 0; i < n; i++) {

				CycleSlipDetect csdObj = csdList.get(i);

				z.set(i, csdObj.getCarrierPhaseDR() - csdObj.getSatVelCorr());
				z.set(i + n, csdObj.getDopplerDR() - csdObj.getSatVelCorr());
				String obsvCode = satList.get(i).getObsvCode();
				double wavelength = csdObj.getWavelength();
				for (int j = 0; j < m; j++) {
					if (obsvCodeList[j].equals(obsvCode)) {
						H.set(i, 3 + j, 1);
						H.set(n + i, 3 + j, 1);
						if (csdObj.isCS()) {
							H.set(i, ctr, wavelength);
							// Cyy_phase.set(i, i, 1e8);
							ctr++;
						}

					}
				}

			}

			SimpleMatrix Ht = H.transpose();
			SimpleMatrix HtWHinv = (Ht.mult(W).mult(H)).invert();
			
			SimpleMatrix x = HtWHinv.mult(Ht).mult(W).mult(z);
			if (ambCount > 0) {
				SimpleMatrix floatAmb = x.extractMatrix(3 + m, 3 + m + ambCount, 0, 1);
				SimpleMatrix floatAmbCov = HtWHinv.extractMatrix(3 + m, 3 + m + ambCount, 3 + m, 3 + m + ambCount);
				System.out.println("Float Ambiguity");
				System.out.println(floatAmb.toString());
				System.out.println("Float Ambiguity Covariance");
				System.out.println(floatAmbCov.toString());
				
				Jama.Matrix ahat = new Jama.Matrix(Matrix.matrix2Array(floatAmb));
				Jama.Matrix Qahat = new Jama.Matrix(Matrix.matrix2Array(floatAmbCov));
				SimpleMatrix afixed = new SimpleMatrix(floatAmb);
				
				Lambda lmd = new Lambda(ahat, Qahat, 6,"MU",(1/3.0),"NCANDS",10);
				int nFixed = lmd.getNfixed();
				double Ps = lmd.getPs();
				if(nFixed==0&&ambCount>1)
				{
					lmd = new Lambda(ahat, Qahat, 5,"MU",(1/3.0),"NCANDS",10);
					Ps = lmd.getPs();
					afixed = new SimpleMatrix(lmd.getafixed().getArray());
					nFixed = lmd.getNfixed();
					
				}
				else
				{
					afixed = new SimpleMatrix(lmd.getafixed().getArray());
				}
				
				if(nFixed!=0)
				{
					SimpleMatrix P = HtWHinv;
					SimpleMatrix Cba = P.extractMatrix(0, 3+m, 3+m, 3+m+ambCount);
					SimpleMatrix Cbb_hat = P.extractMatrix(0, 3+m, 0, 3+m);
					SimpleMatrix b_hat  = x.extractMatrix(0,3+m,0,1);
					SimpleMatrix Caa_inv = floatAmbCov.invert();
					SimpleMatrix a_hat = new SimpleMatrix(floatAmb);
					SimpleMatrix a_inv_hat = afixed.extractMatrix(0, ambCount, 0, 1);
					
					SimpleMatrix b_inv_hat = b_hat.minus(Cba.mult(Caa_inv).mult(a_hat.minus(a_inv_hat))); 
					SimpleMatrix Cbb_inv_hat = Cbb_hat.minus(Cba.mult(Caa_inv).mult(Cba.transpose()));
					
					x = new SimpleMatrix(3+m+ambCount,1);
					x.insertIntoThis(0, 0, b_inv_hat);
					x.insertIntoThis(3+m,0,a_inv_hat);
					P = new SimpleMatrix(Cbb_inv_hat);
					System.out.println("Fixed Ambiguity Sequence");
					System.out.println(a_inv_hat.toString());
					System.out.println(" N Fixed : "+nFixed );
					System.out.println(" Failure Rate : "+(1-Ps));
					ambRepairedCountMap.put(currentTime, nFixed);
					ambRepairedCount += nFixed;
					
					
				}
				
			}
			if (doAnalyze) {
				SimpleMatrix e_hat = (z.minus(H.mult(x))).extractMatrix(0, n, 0, 1);
				residual = Matrix.matrix2ArrayVec(e_hat);
				Cyy = W.invert();
				SimpleMatrix P_H_perpendicular = Matrix.getPerpendicularProjection(H, W);
				SimpleMatrix Cee_hat = P_H_perpendicular.mult(Cyy).mult(P_H_perpendicular.transpose());
				SimpleMatrix redunMatrix = Cee_hat.mult(W);
				double phase_redun = 0;
				for (int i = 0; i < n; i++) {
					phase_redun += redunMatrix.get(i, i);
					residual[i] = residual[i]/csdList.get(i).getWavelength();
				}
				double globalTq = e_hat.transpose().mult(W.extractMatrix(0, n, 0, n)).mult(e_hat).get(0);
				postVarOfUnitW = globalTq * GnssDataConfig.tdcp_priorVarOfUnitW / phase_redun;
				SimpleMatrix R = new SimpleMatrix(3 + m, 3 + m);
				R.insertIntoThis(0, 0, new SimpleMatrix(LatLonUtil.getEcef2EnuRotMat(estEcefClk)));
				for (int i = 0; i < m; i++) {
					R.set(3 + i, 3 + i, 1);
				}
				SimpleMatrix Cxx_hat = HtWHinv.extractMatrix(0, 3 + m, 0, 3 + m);
				Cxx_hat_Map.put("ECEF", Cxx_hat);

				Cxx_hat = R.mult(Cxx_hat).mult(R.transpose());
				Cxx_hat_Map.put("ENU", Cxx_hat);
				testedCsdList = new ArrayList<CycleSlipDetect>(csdList);
			}

			// Velocity
			for(int i=0;i<3+m;i++)
			{
				estVel[i] = x.get(i, 0);
			}
			
			return estVel;
		}

		throw new Exception("Satellite count is less than " + (3 + m) + " , can't compute user position");

	}

	private static double[] estimateVel(ArrayList<CycleSlipDetect> csdList, SimpleMatrix W, double[] estEcefClk,
			boolean useIGS) throws Exception {

		// Satellite count
		int n = csdList.size();
		ArrayList<Satellite> satList = new ArrayList<Satellite>();
		for (int i = 0; i < n; i++) {

			satList.add(csdList.get(i).getSat());

		}
		String[] obsvCodeList = useIGS ? (String[]) SatUtil.findObsvCodeArray(satList) : null;
		int m = obsvCodeList.length;
		double[] estVel = new double[3 + m];
		SimpleMatrix z = new SimpleMatrix(n, 1);
		// Jacobian or Design Matrix
		SimpleMatrix H = new SimpleMatrix(n, 3 + m);

		// Minimum 4 satellite are required to proceed
		if (n >= 3 + m) {
			// Iterate through each satellite, to compute LOS vector and Approx pseudorange
			for (int i = 0; i < n; i++) {

				CycleSlipDetect csdObj = csdList.get(i);
				z.set(i, csdObj.getCarrierPhaseDR() - csdObj.getSatVelCorr());
				String obsvCode = satList.get(i).getObsvCode();
				SimpleMatrix unitLOS = csdObj.getUnitLOS();
				SimpleMatrix negUnitLOS = unitLOS.scale(-1);
				int index = i;
				IntStream.range(0, 3).forEach(j -> H.set(index, j, negUnitLOS.get(j)));
				for (int j = 0; j < m; j++) {
					if (obsvCodeList[j].equals(obsvCode)) {
						H.set(i, 3 + j, 1);

					}
				}

			}

			SimpleMatrix Ht = H.transpose();
			SimpleMatrix HtWHinv = (Ht.mult(W).mult(H)).invert();
			SimpleMatrix x = HtWHinv.mult(Ht).mult(W).mult(z);
			// Velocity
			IntStream.range(0, 3 + m).forEach(i -> estVel[i] = x.get(i, 0));
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

	public static int getCommSatCount() {
		return commonSatCount;
	}

	public static double[] getDop() {
		return dop;
	}

	public static ArrayList<CycleSlipDetect> getCsdList() {
		return testedCsdList;
	}

	public static long getAmbDetectedCount() {
		return ambDetectedCount;
	}

	public static long getAmbRepairedCount() {
		return ambRepairedCount;
	}

	public static TreeMap<Long, Integer> getAmbDetectedCountMap() {
		return ambDetectedCountMap;
	}

	public static TreeMap<Long, Integer> getAmbRepairedCountMap() {
		return ambRepairedCountMap;
	}
	
	
}
