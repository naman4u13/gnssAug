package com.gnssAug.Android.estimation.KalmanFilter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.stream.IntStream;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.ejml.dense.row.MatrixFeatures_DDRM;
import org.ejml.simple.SimpleMatrix;
import com.gnssAug.Android.constants.GnssDataConfig;
import com.gnssAug.Android.estimation.LinearLeastSquare;
import com.gnssAug.Android.estimation.KalmanFilter.Models.KFconfig;
import com.gnssAug.Android.models.CycleSlipDetect;
import com.gnssAug.Android.models.Satellite;
import com.gnssAug.helper.lambdaNew.EstimatorType;
import com.gnssAug.helper.lambdaNew.LAMBDA;
import com.gnssAug.helper.lambdaNew.LAMBDA.LambdaResult;
import com.gnssAug.helper.lambdaNew.LAMBDA_all;
import com.gnssAug.helper.lambdaNew.LAMBDA_all.LambdaAllResult;
import com.gnssAug.utility.LatLonUtil;
import com.gnssAug.utility.MathUtil;
import com.gnssAug.utility.Matrix;
import com.gnssAug.utility.SatUtil;
import com.gnssAug.utility.Weight;

public class EKF_TDCP_ambFix_allEst extends EKFParent {
	private long ambDetectedCount = 0;
	private HashMap<EstimatorType, Long> ambRepairedCount = null;
	private HashMap<EstimatorType, ArrayList<Object[]>> srfrMap = new HashMap<EstimatorType, ArrayList<Object[]>>();
	public EKF_TDCP_ambFix_allEst() {
		kfObj = new KFconfig();
	}

	private TreeMap<Long, Integer> ambDetectedCountMap;
	private TreeMap<Long, HashMap<EstimatorType, Integer>> ambRepairedCountMap;
	private TreeMap<Long, HashMap<EstimatorType, ArrayList<CycleSlipDetect>>> csdListMap;

	public TreeMap<Long, double[]> process(TreeMap<Long, ArrayList<Satellite>> SatMap, ArrayList<Long> timeList,
			boolean useIGS, String[] obsvCodeList, boolean doTest, boolean outlierAnalyze,
			ArrayList<double[]> truePosEcefList) throws Exception {
		boolean isWeighted = false;
		int m = obsvCodeList.length;
		int n = 3 + m;
		SimpleMatrix x = new SimpleMatrix(n, 1);
		SimpleMatrix P = new SimpleMatrix(n, n);
		/*
		 * state XYZ is intialized WLS generated Rx position estimated using first epoch
		 * data, a-priori estimate error covariance matrix for state XYZ is therefore
		 * assigned 25 m^2 value. Other state variables are assigned infinite(big)
		 * variance
		 */
		double[] refPos = LinearLeastSquare.getEstPos(SatMap.firstEntry().getValue(), true, useIGS);
		double[] intialVel = LinearLeastSquare.getEstVel(SatMap.firstEntry().getValue(), true, refPos, useIGS);
		IntStream.range(0, n).forEach(i -> x.set(i, intialVel[i]));

		// Total State
		IntStream.range(0, n).forEach(i -> P.set(i, i, 100));
		kfObj.setState_ProcessCov(x, P);

		ambRepairedCount = new HashMap<EstimatorType, Long>();
		csdListMap = new TreeMap<Long, HashMap<EstimatorType, ArrayList<CycleSlipDetect>>>();
		for (EstimatorType est : new EstimatorType[] { EstimatorType.ILS, EstimatorType.PAR, EstimatorType.IA_FFRT,
				EstimatorType.BIE, EstimatorType.PAR_FFRT }) {
			ambRepairedCount.put(est, (long) 0);

		}
		ambDetectedCountMap = new TreeMap<Long, Integer>();
		ambRepairedCountMap = new TreeMap<Long, HashMap<EstimatorType, Integer>>();
		return iterate(SatMap, timeList, useIGS, obsvCodeList, doTest, outlierAnalyze, isWeighted, truePosEcefList);
	}

	TreeMap<Long, double[]> iterate(TreeMap<Long, ArrayList<Satellite>> SatMap, ArrayList<Long> timeList,
			boolean useIGS, String[] obsvCodeList, boolean doTest, boolean outlierAnalyze, boolean isWeighted,
			ArrayList<double[]> truePosEcefList) throws Exception {
		TreeMap<Long, double[]> estStateMap = new TreeMap<Long, double[]>();
		int m = obsvCodeList.length;
		int x_size = 3 + m;
		long prevTime = timeList.get(0);
		double[] prevTruePos = truePosEcefList.get(0);
		// Start from 2nd epoch
		for (int i = 1; i < timeList.size(); i++) {

			System.out.println("\n\n Epoch : " + i);

			long currentTime = timeList.get(i);
			double deltaT = (currentTime - prevTime) / 1e3;
			ArrayList<Satellite> currSatList = SatMap.get(currentTime);
			ArrayList<Satellite> prevSatList = SatMap.get(prevTime);
			ArrayList<CycleSlipDetect> csdList = new ArrayList<CycleSlipDetect>();
			double[] currentTruePos = truePosEcefList.get(i);
			double[] refPos = LinearLeastSquare.getEstPos(currSatList, true, useIGS);
			int n_curr = currSatList.size();
			int n_prev = prevSatList.size();
			for (int j = 0; j < n_curr; j++) {
				Satellite current_sat = currSatList.get(j);
				String satID = current_sat.getObsvCode() + current_sat.getSvid();
				for (int k = 0; k < n_prev; k++) {
					Satellite prev_sat = prevSatList.get(k);
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
						double trueDR = MathUtil.getEuclidean(currentTruePos, current_sat.getSatEci())
								- MathUtil.getEuclidean(prevTruePos, prev_sat.getSatEci());

						if (approxCS < 5 * wavelength) {
							csdList.add(new CycleSlipDetect(current_sat, dopplerDR, phaseDR, ionoRate, false,
									wavelength, satVelCorr, unitLOS, currentTime - timeList.get(0), trueDR));
						} else if (approxCS < 100 * wavelength) {
							csdList.add(new CycleSlipDetect(current_sat, dopplerDR, phaseDR, ionoRate, false,
									wavelength, satVelCorr, unitLOS, currentTime - timeList.get(0), trueDR));
						}

					}
				}
			}

			HashMap<EstimatorType, ArrayList<CycleSlipDetect>> csdMap = runFilter(currentTime, deltaT, csdList,
					obsvCodeList, doTest, outlierAnalyze, useIGS, i, isWeighted, refPos);
			SimpleMatrix x = kfObj.getState();
			SimpleMatrix P = kfObj.getCovariance();

			double[] estState = new double[x_size];
			IntStream.range(0, x_size).forEach(j -> estState[j] = x.get(j));
			// Add position estimate to the list
			estStateMap.put(currentTime, estState);

			// Convert to ENU frame
			SimpleMatrix R = new SimpleMatrix(x_size, x_size);
			R.insertIntoThis(0, 0, new SimpleMatrix(LatLonUtil.getEcef2EnuRotMat(estState)));
			IntStream.range(3, x_size).forEach(j -> R.set(j, j, 1));
			csdListMap.put(currentTime, csdMap);

			if (!MatrixFeatures_DDRM.isPositiveDefinite(P.getMatrix())) {
				throw new Exception("PositiveDefinite test Failed");
			}
			prevTime = currentTime;
			prevTruePos = currentTruePos;
		}
		return estStateMap;
	}

	// Innovation Based Testing
	private HashMap<EstimatorType, ArrayList<CycleSlipDetect>> runFilter(long currentTime, double deltaT,
			ArrayList<CycleSlipDetect> csdList, String[] obsvCodeList, boolean doTest, boolean outlierAnalyze,
			boolean useIGS, int ct, boolean isWeighted, double[] refPos) throws Exception {

		// Satellite count
		int n = csdList.size();
		ArrayList<Satellite> satList = new ArrayList<Satellite>();
		for (int i = 0; i < n; i++) {
			satList.add(csdList.get(i).getSat());
		}
		int m = obsvCodeList.length;
		// Last update VC-matrix
		SimpleMatrix priorP = new SimpleMatrix(kfObj.getCovariance());

		// Assign Q and F matrix
		kfObj.configTDCP(deltaT, m, refPos);

		kfObj.predict();
		SimpleMatrix x = kfObj.getState();
		SimpleMatrix z = new SimpleMatrix(n, 1);

		SimpleMatrix H = new SimpleMatrix(n, 3 + m);
		SimpleMatrix unitLOS = new SimpleMatrix(SatUtil.getUnitLOS(satList, refPos));
		H.insertIntoThis(0, 0, unitLOS.scale(-1));
		for (int i = 0; i < n; i++) {
			CycleSlipDetect csdObj = csdList.get(i);
			z.set(i, csdObj.getDopplerDR() - csdObj.getSatVelCorr());
			String obsvCode = satList.get(i).getObsvCode();
			for (int j = 0; j < m; j++) {
				if (obsvCodeList[j].equals(obsvCode)) {
					H.set(i, 3 + j, 1);
				}
			}
		}

		SimpleMatrix ze = H.mult(x);
		innovation = Matrix.matrix2ArrayVec(z.minus(ze));

		SimpleMatrix doppler_Cyy = null;
		// Measurement Noise
		if (true) {
			doppler_Cyy = Weight.getNormCyy(satList, GnssDataConfig.doppler_priorVarOfUnitW);
		} else {
			doppler_Cyy = SimpleMatrix.identity(n).scale(GnssDataConfig.doppler_priorVarOfUnitW);
		}

		SimpleMatrix R = new SimpleMatrix(doppler_Cyy);

		ArrayList<Satellite> testedSatList = new ArrayList<Satellite>(satList);
		if (doTest && !outlierAnalyze) {
			Object[] params = performTesting(R, H, n, m, satList, testedSatList, z, ze, csdList, false);
			R = (SimpleMatrix) params[0];
			H = (SimpleMatrix) params[1];
			z = (SimpleMatrix) params[2];
			ze = (SimpleMatrix) params[3];
		}
		// Perform Update Step
		kfObj.update(z, R, ze, H);
		x = kfObj.getState();
		SimpleMatrix P = kfObj.getCovariance();

		// TDCP integration begins

		int ambCount = 0;

		testedSatList = new ArrayList<Satellite>();
		ArrayList<CycleSlipDetect> testedCsdList = new ArrayList<CycleSlipDetect>();

		for (int i = 0; i < n; i++) {
			if (!csdList.get(i).isCS()) {
				testedSatList.add(csdList.get(i).getSat());
				testedCsdList.add(csdList.get(i));
			}
		}
		int tested_n = testedSatList.size();

		z = new SimpleMatrix(tested_n, 1);
		H = new SimpleMatrix(tested_n, 3 + m);
		SimpleMatrix testedUnitLOS = new SimpleMatrix(SatUtil.getUnitLOS(testedSatList, refPos));
		H.insertIntoThis(0, 0, testedUnitLOS.scale(-1));
		for (int i = 0; i < tested_n; i++) {
			CycleSlipDetect csdObj = testedCsdList.get(i);
			z.set(i, csdObj.getCarrierPhaseDR() - csdObj.getSatVelCorr());
			String obsvCode = testedSatList.get(i).getObsvCode();
			for (int j = 0; j < m; j++) {
				if (obsvCodeList[j].equals(obsvCode)) {
					H.set(i, 3 + j, 1);
				}
			}
		}

		ze = H.mult(x);
		innovation = Matrix.matrix2ArrayVec(z.minus(ze));
		// Measurement Noise
		SimpleMatrix tested_tdcp_Cyy = null;
		if (isWeighted) {
			tested_tdcp_Cyy = Weight.getNormCyy(testedSatList, GnssDataConfig.tdcp_priorVarOfUnitW);
		} else {
			tested_tdcp_Cyy = SimpleMatrix.identity(testedSatList.size()).scale(GnssDataConfig.tdcp_priorVarOfUnitW);
		}

		R = new SimpleMatrix(tested_tdcp_Cyy);

		// Testing for CS
		performTesting(R, H, tested_n, m, satList, testedSatList, z, ze, csdList, true);

		// Resume full SatList or CSDList
		for (int i = 0; i < n; i++) {
			if (csdList.get(i).isCS()) {
				ambCount++;
			}
		}

		ambDetectedCountMap.put(currentTime, ambCount);
		ambDetectedCount += ambCount;
		SimpleMatrix x_new = new SimpleMatrix(3 + m + ambCount, 1);
		x_new.insertIntoThis(0, 0, x);
		P = kfObj.getCovariance();
		SimpleMatrix P_new = new SimpleMatrix(3 + m + ambCount, 3 + m + ambCount);
		P_new.insertIntoThis(0, 0, P);
		kfObj.setState_ProcessCov(x_new, P_new);
		for (int i = 0; i < ambCount; i++) {
			P_new.set(3 + m + i, 3 + m + i, 1e12);
		}
		H = new SimpleMatrix(n, 3 + m + ambCount);
		H.insertIntoThis(0, 0, unitLOS.scale(-1));
		z = new SimpleMatrix(n, 1);
		SimpleMatrix tdcp_Cyy = null;
		if (isWeighted) {
			tdcp_Cyy = Weight.getNormCyy(satList, GnssDataConfig.tdcp_priorVarOfUnitW);
		} else {
			tdcp_Cyy = SimpleMatrix.identity(n).scale(GnssDataConfig.tdcp_priorVarOfUnitW);
		}
		R = new SimpleMatrix(tdcp_Cyy);
		double[] sysout_var = new double[ambCount];
		String[] sysout_svid = new String[ambCount];
		int ctr = 3 + m;
		for (int i = 0; i < n; i++) {
			CycleSlipDetect csdObj = csdList.get(i);
			z.set(i, csdObj.getCarrierPhaseDR() - csdObj.getSatVelCorr());
			String obsvCode = satList.get(i).getObsvCode();
			double wavelength = csdObj.getWavelength();
			for (int j = 0; j < m; j++) {
				if (obsvCodeList[j].equals(obsvCode)) {
					H.set(i, 3 + j, 1);
				}
			}
			if (csdObj.isCS()) {
				sysout_var[ctr - 3 - m] = R.get(i, i) / Math.pow(wavelength, 2);
				sysout_svid[ctr - 3 - m] = obsvCode.charAt(0) + "" + satList.get(i).getSvid();
				H.set(i, ctr, wavelength);
				ctr++;
			}
		}
		ze = H.mult(x_new);
		innovation = Matrix.matrix2ArrayVec(z.minus(ze));

		// Perform Update Step
		kfObj.update(z, R, ze, H);
		x = kfObj.getState();
		P = kfObj.getCovariance();

		HashMap<EstimatorType, ArrayList<CycleSlipDetect>> csdMap = new HashMap<EstimatorType, ArrayList<CycleSlipDetect>>();
		for (EstimatorType est : new EstimatorType[] { EstimatorType.ILS, EstimatorType.PAR, EstimatorType.IA_FFRT,
				EstimatorType.BIE, EstimatorType.PAR_FFRT }) {
			ArrayList<CycleSlipDetect> newCsdList = new ArrayList<CycleSlipDetect>();
			for (int i = 0; i < n; i++) {
				newCsdList.add(new CycleSlipDetect(csdList.get(i)));
			}
			csdMap.put(est, newCsdList);
		}
		SimpleMatrix ilsBCaron = x.extractMatrix(0, 3 + m, 0, 1);
		SimpleMatrix ilsCbbCaron = P.extractMatrix(0, 3 + m, 0, 3 + m);
		if (ambCount > 0) {
			SimpleMatrix floatAmb = x.extractMatrix(3 + m, 3 + m + ambCount, 0, 1);
			SimpleMatrix floatAmbCov = P.extractMatrix(3 + m, 3 + m + ambCount, 3 + m, 3 + m + ambCount);
			System.out.println("Float Ambiguity");
			System.out.println(floatAmb.toString());
			System.out.println("Float Ambiguity Covariance");
			System.out.println(floatAmbCov.toString());
			System.out.println("SVIDs:");
			System.out.println(Arrays.toString(sysout_svid));
			System.out.println("TDCP variance:");
			System.out.println(Arrays.toString(sysout_var));
			SimpleMatrix a_hat = new SimpleMatrix(floatAmb);
			SimpleMatrix Q_ahat = new SimpleMatrix(floatAmbCov);
			SimpleMatrix afixed = new SimpleMatrix(floatAmb);
			boolean estimateVar = true;
			// System.out.println(" Failure Rate : " + (1 - Ps));
			HashMap<EstimatorType, LambdaAllResult> lmdMap = LAMBDA_all.computeLambda(a_hat, Q_ahat, estimateVar);

			for (EstimatorType est : lmdMap.keySet()) {
				LambdaAllResult lmd = lmdMap.get(est);
				int nFixed = lmd.getnFixed();
				double Ps = lmd.getSr();
				afixed = lmd.getaFix();
				SimpleMatrix qFixed = lmd.getqFix();

				if (nFixed != 0) {
					SimpleMatrix Cba = P.extractMatrix(0, 3 + m, 3 + m, 3 + m + ambCount);

					SimpleMatrix Cbb_hat = P.extractMatrix(0, 3 + m, 0, 3 + m);
					SimpleMatrix b_hat = x.extractMatrix(0, 3 + m, 0, 1);
					SimpleMatrix Caa_hat_inv = floatAmbCov.invert();
					SimpleMatrix Caa_caron = new SimpleMatrix(qFixed);
					SimpleMatrix a_caron = new SimpleMatrix(afixed);
					SimpleMatrix b_caron = b_hat.minus(Cba.mult(Caa_hat_inv).mult(a_hat.minus(a_caron)));
					SimpleMatrix fixedTermContri = Cba.mult(Caa_hat_inv).mult(Caa_caron).mult(Caa_hat_inv)
							.mult(Cba.transpose());
					SimpleMatrix Cbb_caron = Cbb_hat.minus(Cba.mult(Caa_hat_inv).mult(Cba.transpose()))
							.plus(fixedTermContri);
					if (est == EstimatorType.PAR_FFRT) {

						ilsBCaron = new SimpleMatrix(b_caron);
						x.insertIntoThis(3 + m, 0, a_caron);
						ilsCbbCaron = new SimpleMatrix(Cbb_caron);
					}
					System.out.println("Fixed Ambiguity Sequence : " + est.toString());
					System.out.println(a_caron.toString());
					System.out.println(" N Fixed : " + nFixed);
					ambRepairedCountMap.computeIfAbsent(currentTime, k -> new HashMap<EstimatorType, Integer>())
							.put(est, nFixed);
					ambRepairedCount.put(est, ambRepairedCount.get(est) + nFixed);
					int count = 0;

					for (int i = 0; i < n; i++) {
						CycleSlipDetect csdObj = csdMap.get(est).get(i);
						if (csdObj.isCS()) {
							csdObj.setRepaired(true);
							csdObj.setIntAmb(a_caron.get(count));
							csdObj.setIntAmbCov(Caa_caron.get(count, count));
							csdObj.setFloatAmb(a_hat.get(count));
							csdObj.setFloatAmbCov(Q_ahat.get(count,count));
							count++;
						}

					}
					if(estimateVar)
					{
						srfrMap.computeIfAbsent(est, k->new ArrayList<Object[]>()).add(new Object[] {ambCount,nFixed,lmd.getApproxSR(),lmd.getApproxFR()});
					}

				}

			}

		}

		x = new SimpleMatrix(ilsBCaron);
		P = new SimpleMatrix(ilsCbbCaron);
		kfObj.setState_ProcessCov(x, P);
		return csdMap;
	}

	private Object[] performTesting(SimpleMatrix R, SimpleMatrix H, int n, int m, ArrayList<Satellite> satList,
			ArrayList<Satellite> testedSatList, SimpleMatrix z, SimpleMatrix ze, ArrayList<CycleSlipDetect> csdList,
			boolean isTDCP) throws Exception {
		// Pre-fit residual/innovation
		SimpleMatrix v = new SimpleMatrix(n, 1, true, innovation);
		SimpleMatrix P = kfObj.getCovariance();
		SimpleMatrix Cvv = ((H.mult(P).mult(H.transpose())).plus(R));
		SimpleMatrix Cvv_inv = Cvv.invert();
		double globalTq = v.transpose().mult(Cvv_inv).mult(v).get(0);
		ChiSquaredDistribution csd = new ChiSquaredDistribution(n);
		double alpha = 0.01;
		if (isTDCP) {
			alpha = 0.01;
		}
		if (globalTq == 0) {
			throw new Exception("Error: T stat is zero");
		}
		// Detection
		double globalPVal = 1 - csd.cumulativeProbability(globalTq);
		int _n = testedSatList.size();

		while (globalPVal < alpha && _n > (n / 2)) {

			double max_w = Double.MIN_VALUE;
			int index = -1;
			for (int i = 0; i < _n; i++) {
				SimpleMatrix cv = new SimpleMatrix(_n, 1);
				cv.set(i, 1);
				double w = Math.abs(cv.transpose().mult(Cvv_inv).mult(v).get(0))
						/ Math.sqrt(cv.transpose().mult(Cvv_inv).mult(cv).get(0));
				if (w > max_w) {
					max_w = w;
					index = i;
				}
			}
			int rem_index = satList.indexOf(testedSatList.remove(index));
			if (isTDCP) {

				satList.get(rem_index).setOutlier(true);
				csdList.get(rem_index).setCS(true);
			}
			_n = testedSatList.size();
			SimpleMatrix R_ = new SimpleMatrix(_n, _n);
			SimpleMatrix z_ = new SimpleMatrix(_n, 1);
			SimpleMatrix ze_ = new SimpleMatrix(_n, 1);
			SimpleMatrix H_ = new SimpleMatrix(_n, 3 + m);
			SimpleMatrix v_ = new SimpleMatrix(_n, 1);
			int j = 0;
			for (int i = 0; i < _n + 1; i++) {
				if (i != index) {
					R_.set(j, j, R.get(i, i));
					v_.set(j, v.get(i));
					z_.set(j, z.get(i));
					ze_.set(j, ze.get(i));
					for (int k = 0; k < 3 + m; k++) {
						H_.set(j, k, H.get(i, k));
					}
					j++;
				}
			}
			R = new SimpleMatrix(R_);
			H = new SimpleMatrix(H_);
			z = new SimpleMatrix(z_);
			ze = new SimpleMatrix(ze_);
			v = new SimpleMatrix(v_);

			Cvv = ((H.mult(P).mult(H.transpose())).plus(R));
			Cvv_inv = Cvv.invert();
			globalTq = v.transpose().mult(Cvv_inv).mult(v).get(0);
			if (globalTq == 0) {
				throw new Exception("Error: T stat is zero");
			}
			csd = new ChiSquaredDistribution(_n);
			globalPVal = 1 - csd.cumulativeProbability(globalTq);
		}
		return new Object[] { R, H, z, ze };

	}

	public long getAmbDetectedCount() {
		return ambDetectedCount;
	}

	public HashMap<EstimatorType, Long> getAmbRepairedCount() {
		return ambRepairedCount;
	}

	public TreeMap<Long, Integer> getAmbDetectedCountMap() {
		return ambDetectedCountMap;
	}

	public TreeMap<Long, HashMap<EstimatorType, Integer>> getAmbRepairedCountMap() {
		return ambRepairedCountMap;
	}

	public TreeMap<Long, HashMap<EstimatorType, ArrayList<CycleSlipDetect>>> getCsdListMap() {
		return csdListMap;
	}

	private boolean notNaNOrInfinity(double d) {
		return !(Double.isNaN(d) || Double.isInfinite(d));
	}

	public HashMap<EstimatorType, ArrayList<Object[]>> getSrfrMap() {
		return srfrMap;
	}
	
}
