package com.gnssAug.Rinex.estimation;

import java.time.ZonedDateTime;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;
import java.util.stream.IntStream;
import org.apache.commons.collections.set.ListOrderedSet;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.ejml.dense.row.MatrixFeatures_DDRM;
import org.ejml.simple.SimpleMatrix;
import com.gnssAug.Rinex.constants.GnssDataConfig;
import com.gnssAug.Android.constants.Measurement;
import com.gnssAug.Android.estimation.KalmanFilter.EKFParent;
import com.gnssAug.Android.estimation.KalmanFilter.Models.KFconfig;
import com.gnssAug.Android.models.CycleSlipDetect;
import com.gnssAug.Rinex.models.Satellite;
import com.gnssAug.helper.ComputeSolidEarthTide;
import com.gnssAug.helper.lambdaNew.EstimatorType;
import com.gnssAug.helper.lambdaNew.LAMBDA;
import com.gnssAug.helper.lambdaNew.LAMBDA.LambdaResult;
import com.gnssAug.utility.LatLonUtil;
import com.gnssAug.utility.MathUtil;
import com.gnssAug.utility.Matrix;
import com.gnssAug.utility.SatUtil;
import com.gnssAug.utility.Time;
import com.gnssAug.utility.Vector;
import com.gnssAug.utility.Weight;

public class EKF_PPP extends EKFParent {

	private TreeMap<Long, Map<Measurement, double[]>> innovationMap;
	private TreeMap<Long, Map<Measurement, double[]>> residualMap;
	private TreeMap<Long, Map<Measurement, Double>> postVarOfUnitWMap;
	private TreeMap<Long, Map<Measurement, Double>> redundancyNoMap;
	private TreeMap<Long, HashMap<String, Double>> ambMap;
	private TreeMap<Long, HashMap<String, Double>> ionoMap;
	private TreeMap<Long, double[]> dopMap;
	private TreeMap<Long, Double> tropoMap;
	private TreeMap<Long, double[]> clkOffMap;
	private TreeMap<Long, double[]> clkDriftMap;
	private long csDetectedCount = 0;
	private long csRepairedCount = 0;
	private long ambFixedCount = 0;
	private KFconfig vel_kfObj;
	private KFconfig pos_kfObj;
	protected TreeMap<Long, ArrayList<Satellite>> satListMap;
	private HashMap<String, Integer> consecutiveCSmap;
	final private int csThresh = 100;
	final private int consecutiveSlips = 1;
	private int nonZeroArtificalCScount = 0;


	
	public EKF_PPP() {
		vel_kfObj = new KFconfig();
		pos_kfObj = new KFconfig();
	}

	private TreeMap<Long, Integer> csDetectedCountMap;
	private TreeMap<Long, Integer> csRepairedCountMap;
	private TreeMap<Long, Integer> ambFixedCountMap;
	private HashMap<String, int[]> cycleSlipCount;
	private TreeMap<Long, ArrayList<CycleSlipDetect>> csdListMap;

	

	public TreeMap<Long, double[]> process(TreeMap<Long, ArrayList<Satellite>> SatMap, HashMap<String, double[]> rxPCO,
			ArrayList<Long> timeList, boolean doAnalyze, boolean doTest, String[] obsvCodeList, double[] rxARP,
			boolean isStatic, boolean repairCS, boolean fixAmb) throws Exception {
		boolean predictPhaseClock = false;
		int clkOffNum = obsvCodeList.length;
		ListOrderedSet ssiSet = new ListOrderedSet();
		for (int i = 0; i < clkOffNum; i++) {
			ssiSet.add(obsvCodeList[i].charAt(0));
		}
		int clkDriftNum = ssiSet.size();

		SimpleMatrix x_vel = new SimpleMatrix(3 + (2 * clkDriftNum), 1);
		SimpleMatrix P_vel = new SimpleMatrix(3 + (2 * clkDriftNum), 3 + (2 * clkDriftNum));
		SimpleMatrix x_pos = new SimpleMatrix(6 + clkOffNum + clkDriftNum + 1, 1);
		SimpleMatrix P_pos = new SimpleMatrix(6 + clkOffNum + clkDriftNum + 1, 6 + clkOffNum + clkDriftNum + 1);

		// Initializing position state
		double[] intialPos = LinearLeastSquare.getEstPos(SatMap.firstEntry().getValue(), rxPCO, true);
		for (int i = 0; i < 3; i++) {
			x_pos.set(i, intialPos[i]);
			P_pos.set(i, i, 1e4);
		}
		// Initializing clock offset state
		IntStream.range(3, 3 + clkOffNum).forEach(i -> P_pos.set(i, i, 1e4));

		double[] intialVel = null;
		if (isStatic) {
			intialVel = new double[3 + clkDriftNum];
			IntStream.range(0, 3).forEach(i -> P_vel.set(i, i, 1e-16));
			IntStream.range(3 + clkOffNum, 6 + clkOffNum).forEach(i -> P_pos.set(i, i, 1e-16));
		} else {

			intialVel = LinearLeastSquare.getEstVel(SatMap.firstEntry().getValue(), rxPCO, true, false, false, false,
					intialPos);
			IntStream.range(0, 3).forEach(i -> P_vel.set(i, i, 100));
			IntStream.range(3 + clkOffNum, 3 + clkOffNum + 3).forEach(i -> P_pos.set(i, i, 100));
		}
		// Initializing velocity state
		for (int i = 0; i < 3; i++) {
			x_vel.set(i, intialVel[i]);
			x_pos.set(3 + clkOffNum + i, intialVel[i]);
		}

		// Initializing clock drift state and Tropo
		IntStream.range(6 + clkOffNum, 6 + clkOffNum + clkDriftNum + 1).forEach(i -> P_pos.set(i, i, 1e4));
		IntStream.range(3, 3 + (2 * clkDriftNum)).forEach(i -> P_vel.set(i, i, 1e4));

		vel_kfObj.setState_ProcessCov(x_vel, P_vel);
		pos_kfObj.setState_ProcessCov(x_pos, P_pos);
		if (doAnalyze) {
			innovationMap = new TreeMap<Long, Map<Measurement, double[]>>();
			errCovMap = new TreeMap<Long, SimpleMatrix>();
			residualMap = new TreeMap<Long, Map<Measurement, double[]>>();
			postVarOfUnitWMap = new TreeMap<Long, Map<Measurement, Double>>();
			redundancyNoMap = new TreeMap<Long, Map<Measurement, Double>>();
			redundancyList = new ArrayList<double[]>();
			satCountMap = new TreeMap<Long, Long>();
			satListMap = new TreeMap<Long, ArrayList<Satellite>>();
			measNoiseMap = new TreeMap<Long, double[]>();
			csdListMap = new TreeMap<Long, ArrayList<CycleSlipDetect>>();
			ambMap = new TreeMap<Long, HashMap<String, Double>>();
			ionoMap = new TreeMap<Long, HashMap<String, Double>>();
			tropoMap = new TreeMap<Long, Double>();
			clkOffMap = new TreeMap<Long, double[]>();
			clkDriftMap = new TreeMap<Long, double[]>();
			dopMap = new TreeMap<Long, double[]>();

		}
		cycleSlipCount = new HashMap<String, int[]>();
		csDetectedCountMap = new TreeMap<Long, Integer>();
		csRepairedCountMap = new TreeMap<Long, Integer>();
		ambFixedCountMap = new TreeMap<Long, Integer>();
		consecutiveCSmap = new HashMap<String, Integer>();
		return iterate(SatMap, rxPCO, timeList, ssiSet, doAnalyze, doTest, rxARP, obsvCodeList, repairCS, fixAmb,predictPhaseClock);

	}

	private TreeMap<Long, double[]> iterate(TreeMap<Long, ArrayList<Satellite>> SatMap, HashMap<String, double[]> rxPCO,
			ArrayList<Long> timeList, ListOrderedSet ssiSet, boolean doAnalyze, boolean doTest, double[] rxARP,
			String[] obsvCodeList, boolean repairCS, boolean fixAmb,boolean predictPhaseClock) throws Exception {

		TreeMap<Long, double[]> estStateMap = new TreeMap<Long, double[]>();
		HashMap<String, Integer> ambMap = new HashMap<String, Integer>();
		HashMap<String, Integer> ionoMap = new HashMap<String, Integer>();
		long prevTime = timeList.get(0);
		double[] slip = new double[] {1,2,3,4,5,6,7};
		// Start from 2nd epoch
		for (int i = 1; i < timeList.size(); i++) {
			System.out.println("\n\n Epoch : " + i);
			long currentTime = timeList.get(i);
			double deltaT = (currentTime - prevTime) / 1e3;
			ArrayList<Satellite> currSatList = SatMap.get(currentTime);
			ArrayList<Satellite> prevSatList = SatMap.get(prevTime);
			int csIncrement = 0;
			if(i%4==0)
			{
				csIncrement =((int)(i/4))%10;
				if(csIncrement!=0)
				{
					nonZeroArtificalCScount++;
				}
				slip[0] += csIncrement;
				System.out.println(currSatList.get(0).getSVID()+" : CS value added = "+slip[0]);
			}
			if(i%3==0)
			{
				csIncrement =((int)(i/3))%10;
				if(csIncrement!=0)
				{
					nonZeroArtificalCScount++;
				}
				slip[1] += csIncrement;
				System.out.println(currSatList.get(1).getSVID()+" : CS value added = "+slip[1]);
			}
			if(i%6==0)
			{
				csIncrement =((int)(i/6))%10;
				if(csIncrement!=0)
				{
					nonZeroArtificalCScount++;
				}
				slip[2] += csIncrement;
				System.out.println(currSatList.get(2).getSVID()+" : CS value added = "+slip[2]);
			}
			if(i%8==0)
			{
				csIncrement =((int)(i/8))%10;
				if(csIncrement!=0)
				{
					nonZeroArtificalCScount++;
				}
				slip[3] += csIncrement;
				System.out.println(currSatList.get(3).getSVID()+" : CS value added = "+slip[3]);
			}
			if(i%4==0)
			{
				csIncrement =((int)(i/4))%10;
				if(csIncrement!=0)
				{
					nonZeroArtificalCScount++;
				}
				slip[4] += csIncrement;
				System.out.println(currSatList.get(4).getSVID()+" : CS value added = "+slip[4]);
			}
			if(i%6==0)
			{
				csIncrement =((int)(i/6))%10;
				if(csIncrement!=0)
				{
					nonZeroArtificalCScount++;
				}
				slip[5] += csIncrement;
				System.out.println(currSatList.get(5).getSVID()+" : CS value added = "+slip[5]);
			}
////			if(i%3==0)
////			{
////				slip[6] += ((int)(i/3))%10;
////				System.out.println(currSatList.get(6).getSVID()+" : CS value added = "+slip[6]);
////			}
			
			Satellite sat_temp = currSatList.get(0);
			sat_temp.setPhase(sat_temp.getPhase()+(sat_temp.getCarrier_wavelength()*slip[0]));
			sat_temp = currSatList.get(1);
			sat_temp.setPhase(sat_temp.getPhase()+(sat_temp.getCarrier_wavelength()*slip[1]));
			sat_temp = currSatList.get(2);
			sat_temp.setPhase(sat_temp.getPhase()+(sat_temp.getCarrier_wavelength()*slip[2]));
			sat_temp = currSatList.get(3);
			sat_temp.setPhase(sat_temp.getPhase()+(sat_temp.getCarrier_wavelength()*slip[3]));
			sat_temp = currSatList.get(4);
			sat_temp.setPhase(sat_temp.getPhase()+(sat_temp.getCarrier_wavelength()*slip[4]));
			sat_temp = currSatList.get(5);
			sat_temp.setPhase(sat_temp.getPhase()+(sat_temp.getCarrier_wavelength()*slip[5]));
//			sat_temp = currSatList.get(6);
//			sat_temp.setPhase(sat_temp.getPhase()+(sat_temp.getCarrier_wavelength()*slip[6]));

			ArrayList<CycleSlipDetect> csdList = new ArrayList<CycleSlipDetect>();
			double[] refPos = LinearLeastSquare.getEstPos(currSatList, rxPCO, true);
			int n_curr = currSatList.size();
			int n_prev = prevSatList.size();
			ZonedDateTime zdTime = Time.convertUsingToInstant(currSatList.get(0).getTime());
			double[] timeVaryingTides = ComputeSolidEarthTide.calculateTimeVaryingTides(refPos, false, zdTime);
			double[] permanentTide = ComputeSolidEarthTide.getMeanTideCorrection(refPos);
			SimpleMatrix earthTide = new SimpleMatrix(3, 1, true, Vector.add(timeVaryingTides, permanentTide));
			for (int j = 0; j < n_curr; j++) {
				Satellite current_sat = currSatList.get(j);
				String satID = current_sat.getObsvCode() + current_sat.getSVID();
				for (int k = 0; k < n_prev; k++) {
					Satellite prev_sat = prevSatList.get(k);
					String prev_satID = prev_sat.getObsvCode() + prev_sat.getSVID();
					if (satID.equals(prev_satID)) {

						SimpleMatrix satEci = new SimpleMatrix(3, 1, true, current_sat.getSatEci());
						SimpleMatrix prev_satEci = new SimpleMatrix(3, 1, true, prev_sat.getSatEci());
						SimpleMatrix unitLOS = new SimpleMatrix(1, 3, true,
								SatUtil.getUnitLOS(current_sat.getSatEci(), refPos));
						double earthTide_range = unitLOS.mult(earthTide).get(0);
						current_sat.setPseudorange(current_sat.getPseudorange() + earthTide_range);
						current_sat.setPhase(current_sat.getPhase() + earthTide_range);
						double ionoRate = current_sat.getIonoErr() - prev_sat.getIonoErr();
						double tropoRate = current_sat.getTropoErr() - prev_sat.getTropoErr();
						double dopplerDR = ((current_sat.getPseudoRangeRate() + prev_sat.getPseudoRangeRate()) / 2)
								- tropoRate + ionoRate;
						double phaseDR = current_sat.getPhase() - prev_sat.getPhase() + ionoRate;
						double prDR = current_sat.getPseudorange() - prev_sat.getPseudorange() - ionoRate;
						double satVelCorr = unitLOS.mult(satEci.minus(prev_satEci)).get(0);
						double wavelength = SpeedofLight / current_sat.getCarrier_frequency();
						double approxCS = Math.abs(phaseDR - dopplerDR);
						double[] pco = rxPCO.get(current_sat.getObsvCode());
						double trueDR = MathUtil.getEuclidean(Vector.add(rxARP, pco), current_sat.getSatEci())
								- MathUtil.getEuclidean(Vector.add(rxARP, pco), prev_sat.getSatEci());

						if (approxCS < csThresh * wavelength) {
							csdList.add(new CycleSlipDetect(current_sat, dopplerDR, phaseDR, ionoRate, false,
									wavelength, satVelCorr, unitLOS, currentTime - timeList.get(0), trueDR, prDR));
						}

					}
				}
			}
			runFilter_vel(currentTime, deltaT, csdList, ssiSet, true, i, refPos, repairCS, doTest);

			runFilter_pos(currentTime, deltaT, csdList, ssiSet, ambMap, ionoMap, obsvCodeList, rxPCO, fixAmb, doAnalyze,
					doTest,predictPhaseClock);
			SimpleMatrix x = pos_kfObj.getState();
			SimpleMatrix P = pos_kfObj.getCovariance();

			double[] estState = new double[3];
			for (int j = 0; j < 3; j++) {
				estState[j] = x.get(j);
			}
			
			// Add position estimate to the list
			estStateMap.put(currentTime, estState);
			prevTime = currentTime;

		}
		System.out.println("Artificial CS count : "+nonZeroArtificalCScount);
		return estStateMap;

	}

	private void runFilter_vel(long currentTime, double deltaT, ArrayList<CycleSlipDetect> csdList,
			ListOrderedSet ssiSet, boolean useIGS, int ct, double[] refPos, boolean repairCS, boolean doTest)
			throws Exception {

		// Satellite count
		int n = csdList.size();
		ArrayList<Satellite> satList = new ArrayList<Satellite>();
		for (int i = 0; i < n; i++) {
			satList.add(csdList.get(i).getIgs_sat());
		}
		int m = ssiSet.size();

		// Assign Q and F matrix
		vel_kfObj.configTDCP(deltaT, m, refPos);
		vel_kfObj.predict();
		SimpleMatrix x = vel_kfObj.getState();
		SimpleMatrix z = new SimpleMatrix(n, 1);

		SimpleMatrix H = new SimpleMatrix(n, 3 + (2 * m));
		SimpleMatrix unitLOS = null;
		unitLOS = new SimpleMatrix(SatUtil.igs_getUnitLOS(satList, refPos));
		H.insertIntoThis(0, 0, unitLOS.scale(-1));
		for (int i = 0; i < n; i++) {
			CycleSlipDetect csdObj = csdList.get(i);
			z.set(i, csdObj.getDopplerDR() - csdObj.getSatVelCorr());
			char ssi = satList.get(i).getObsvCode().charAt(0);

			for (int j = 0; j < m; j++) {
				if ((char) ssiSet.get(j) == ssi) {
					H.set(i, 3 + j, 1);
				}
			}
		}

		SimpleMatrix ze = H.mult(x);
		innovation = Matrix.matrix2ArrayVec(z.minus(ze));

		SimpleMatrix doppler_Cyy = null;
		// Measurement Noise
		doppler_Cyy = SimpleMatrix.identity(satList.size()).scale(GnssDataConfig.doppler_priorVarOfUnitW);

		SimpleMatrix R = new SimpleMatrix(doppler_Cyy);

		ArrayList<Satellite> testedSatList = new ArrayList<Satellite>(satList);
		if (doTest) {
			Object[] params = performTesting(R, H, n, m, satList, testedSatList, z, ze, csdList, false, vel_kfObj);
			R = (SimpleMatrix) params[0];
			H = (SimpleMatrix) params[1];
			z = (SimpleMatrix) params[2];
			ze = (SimpleMatrix) params[3];
		}
		// Perform Update Step
		vel_kfObj.update(z, R, ze, H);
		x = vel_kfObj.getState();
		SimpleMatrix P = vel_kfObj.getCovariance();

		// TDCP integration begins

		int ambCount = 0;

		testedSatList = new ArrayList<Satellite>();
		ArrayList<CycleSlipDetect> testedCsdList = new ArrayList<CycleSlipDetect>();

		for (int i = 0; i < n; i++) {
			if (!csdList.get(i).isCS()) {
				testedSatList.add(csdList.get(i).getIgs_sat());
				testedCsdList.add(csdList.get(i));
			}
		}
		int tested_n = testedSatList.size();

		z = new SimpleMatrix(tested_n, 1);
		H = new SimpleMatrix(tested_n, 3 + (2 * m));
		SimpleMatrix testedUnitLOS = new SimpleMatrix(SatUtil.igs_getUnitLOS(testedSatList, refPos));
		H.insertIntoThis(0, 0, testedUnitLOS.scale(-1));
		for (int i = 0; i < tested_n; i++) {
			CycleSlipDetect csdObj = testedCsdList.get(i);
			z.set(i, csdObj.getCarrierPhaseDR() - csdObj.getSatVelCorr());

			char ssi = testedSatList.get(i).getObsvCode().charAt(0);
			for (int j = 0; j < m; j++) {
				if ((char) ssiSet.get(j) == ssi) {
					H.set(i, 3 + m + j, 1);
				}
			}
		}

		ze = H.mult(x);
		innovation = Matrix.matrix2ArrayVec(z.minus(ze));
		// Measurement Noise
		SimpleMatrix tested_tdcp_Cyy = SimpleMatrix.identity(testedSatList.size())
				.scale(GnssDataConfig.tdcp_priorVarOfUnitW);

		R = new SimpleMatrix(tested_tdcp_Cyy);
		int excludeCount = 0;
		// Testing for CS
		performTesting(R, H, tested_n, m, satList, testedSatList, z, ze, csdList, true, vel_kfObj);
		HashMap<String, Integer> currConsecutiveCSmap = new HashMap<String, Integer>();
		// Resume full SatList or CSDList
		for (int i = 0; i < n; i++) {
			Satellite sat = csdList.get(i).getIgs_sat();
			String satID = sat.getObsvCode() + "" + sat.getSVID();
			int[] record = cycleSlipCount.computeIfAbsent(satID, k -> new int[2]);
			record[1] += 1;
			if (csdList.get(i).isCS()) {
				record[0] += 1;
				cycleSlipCount.put(satID, record);
				if (consecutiveCSmap.containsKey(satID)) {
					int count = consecutiveCSmap.get(satID) + 1;
					currConsecutiveCSmap.put(satID, count);
					if (count > consecutiveSlips) {
						csdList.get(i).setReset(true);
						excludeCount++;
						continue;
					}
				} else {
					currConsecutiveCSmap.put(satID, 1);
				}
				ambCount++;
			}
			
		}
		consecutiveCSmap = new HashMap<String, Integer>(currConsecutiveCSmap);
		csDetectedCountMap.put(currentTime, ambCount);
		csDetectedCount += ambCount;

		SimpleMatrix x_new = new SimpleMatrix(3 + (2 * m) + ambCount, 1);
		x_new.insertIntoThis(0, 0, x);
		P = vel_kfObj.getCovariance();
		SimpleMatrix P_new = new SimpleMatrix(3 + (2 * m) + ambCount, 3 + (2 * m) + ambCount);
		P_new.insertIntoThis(0, 0, P);
		vel_kfObj.setState_ProcessCov(x_new, P_new);
		for (int i = 0; i < ambCount; i++) {
			P_new.set(3 + (2 * m) + i, 3 + (2 * m) + i, 1e12);
		}
		H = new SimpleMatrix(n-excludeCount, 3 + (2 * m) + ambCount);
		
		z = new SimpleMatrix(n-excludeCount, 1);
		SimpleMatrix tdcp_Cyy = SimpleMatrix.identity(n-excludeCount).scale(GnssDataConfig.tdcp_priorVarOfUnitW);

		R = new SimpleMatrix(tdcp_Cyy);
		String[] sysout_svid = new String[ambCount];
		int ctr = 3 + (2 * m);
		int _i = 0;
		for (int i = 0; i < n; i++) {
			CycleSlipDetect csdObj = csdList.get(i);
			if (!csdObj.isReset()) {
				z.set(_i, csdObj.getCarrierPhaseDR() - csdObj.getSatVelCorr());
				H.insertIntoThis(_i, 0, unitLOS.scale(-1).extractMatrix(i, i+1, 0, 3));
				double wavelength = csdObj.getWavelength();
				char ssi = satList.get(i).getObsvCode().charAt(0);
				for (int j = 0; j < m; j++) {
					if ((char) ssiSet.get(j) == ssi) {
						H.set(_i, 3 + m + j, 1);
					}
				}

				if (csdObj.isCS()) {
					sysout_svid[ctr - 3 - (2 * m)] = satList.get(i).getObsvCode()+""+satList.get(i).getSVID();
					H.set(_i, ctr, wavelength);
					ctr++;
				}
				_i++;

			}
		}

		ze = H.mult(x_new);
		innovation = Matrix.matrix2ArrayVec(z.minus(ze));

		// Perform Update Step
		vel_kfObj.update(z, R, ze, H);
		x = vel_kfObj.getState();
		P = vel_kfObj.getCovariance();

		if (ambCount > 0 && repairCS) {
			SimpleMatrix floatAmb = x.extractMatrix(3 + (2 * m), 3 + (2 * m) + ambCount, 0, 1);
			SimpleMatrix floatAmbCov = P.extractMatrix(3 + (2 * m), 3 + (2 * m) + ambCount, 3 + (2 * m),
					3 + (2 * m) + ambCount);
			System.out.println("Float Cyle Slip");
			System.out.println(floatAmb.toString());
			System.out.println("Float Cyle Slip Covariance");
//			System.out.println(floatAmbCov.toString());
//			System.out.println("Float Cyle Slip Correlation Matrix");
			System.out.println(Matrix.computeCorrelationMatrix(floatAmbCov).toString());
			System.out.println("SVIDs:");
			System.out.println(Arrays.toString(sysout_svid));
//			System.out.println("TDCP variance:");
//			System.out.println(Arrays.toString(sysout_var));

			SimpleMatrix a_hat = new SimpleMatrix(floatAmb);
			SimpleMatrix Q_ahat = new SimpleMatrix(floatAmbCov);
			Q_ahat = Q_ahat.plus(Q_ahat.transpose()).scale(0.5);
			SimpleMatrix afixed = new SimpleMatrix(floatAmb);
			boolean estimateVar = false;
			LambdaResult lmd = LAMBDA.computeLambda(a_hat, Q_ahat, EstimatorType.PAR_FFRT, estimateVar);
			int nFixed = lmd.getnFixed();
			double Ps = lmd.getSr();
			afixed = lmd.getaFix();
			SimpleMatrix QaFixed = lmd.getQaFix();
			System.out.println(" Failure Rate : " + (1 - Ps));
			if (nFixed != 0) {
				SimpleMatrix Cba = P.extractMatrix(0, 3 + (2 * m), 3 + (2 * m), 3 + (2 * m) + ambCount);
				SimpleMatrix Cbb_hat = P.extractMatrix(0, 3 + (2 * m), 0, 3 + (2 * m));
				SimpleMatrix b_hat = x.extractMatrix(0, 3 + (2 * m), 0, 1);
				SimpleMatrix Caa_hat_inv = floatAmbCov.invert();
				SimpleMatrix Caa_caron = new SimpleMatrix(QaFixed);
				SimpleMatrix a_caron = new SimpleMatrix(afixed);
				SimpleMatrix b_caron = b_hat.minus(Cba.mult(Caa_hat_inv).mult(a_hat.minus(a_caron)));
				SimpleMatrix fixedTermContri = Cba.mult(Caa_hat_inv).mult(Caa_caron).mult(Caa_hat_inv)
						.mult(Cba.transpose());
				SimpleMatrix Cbb_caron = Cbb_hat.minus(Cba.mult(Caa_hat_inv).mult(Cba.transpose()))
						.plus(fixedTermContri);
				x = new SimpleMatrix(3 + (2 * m) + ambCount, 1);
				x.insertIntoThis(0, 0, b_caron);
				x.insertIntoThis(3 + (2 * m), 0, a_caron);
				P = new SimpleMatrix(Cbb_caron);
				System.out.println("Fixed Cyle Slip Sequence");
				System.out.println(a_caron.toString());
				System.out.println("Fixed Cyle Slip Variance");
				System.out.println(Caa_caron.toString());
				System.out.println(" N Fixed : " + nFixed);

				csRepairedCountMap.put(currentTime, nFixed);
				csRepairedCount += nFixed;

				int count = 0;
				for (int i = 0; i < n; i++) {
					CycleSlipDetect csdObj = csdList.get(i);
					if (!csdObj.isReset() && csdObj.isCS()) {

						csdObj.setRepaired(true);
						csdObj.setIntAmb(a_caron.get(count));
						csdObj.setIntAmbCov(Caa_caron.get(count, count));
						csdObj.setFloatAmb(a_hat.get(count));
						csdObj.setFloatAmbCov(Q_ahat.get(count, count));
						count++;
					}
				}

			}

		}

		x = x.extractMatrix(0, 3 + (2 * m), 0, 1);
		P = P.extractMatrix(0, 3 + (2 * m), 0, 3 + (2 * m));
		vel_kfObj.setState_ProcessCov(x, P);

	}

	private void runFilter_pos(long currentTime, double deltaT, ArrayList<CycleSlipDetect> csdList,
			ListOrderedSet ssiSet, HashMap<String, Integer> ambMap, HashMap<String, Integer> ionoMap,
			String[] obsvCodeList, HashMap<String, double[]> rxPCO, boolean fixAmb, boolean doAnalyze,
			boolean doTest,boolean predictPhaseClock) throws Exception {
		HashMap<String, Integer> new_ambMap = new HashMap<String, Integer>();
		HashMap<String, Integer> new_ionoMap = new HashMap<String, Integer>();
		ListOrderedSet uniqSat = new ListOrderedSet();
		ArrayList<double[]> ionoParams = new ArrayList<double[]>();
		ArrayList<String> list_ = new ArrayList<String>();
		// Satellite count
		int n = csdList.size();
		ArrayList<Satellite> satList = new ArrayList<Satellite>();
		for (int i = 0; i < n; i++) {
			CycleSlipDetect csdObj = csdList.get(i);
			Satellite sat = csdObj.getIgs_sat();
			String satID = sat.getObsvCode().charAt(0) + "" + sat.getSVID();
			satList.add(sat);
			list_.add(satID);
			if (!uniqSat.contains(satID)) {
				uniqSat.add(satID);
				double freq2 = Math.pow(sat.getCarrier_frequency(), 2);
				double ionoCoeff = (40.3 * 1e16) / freq2;
				ionoParams.add(new double[] { sat.getElevAzm()[0], sat.getIonoErr() / ionoCoeff });
			}
		}
		int clkOffNum = obsvCodeList.length;
		int clkDriftNum = ssiSet.size();
		int ionoParamNum = uniqSat.size();
		int coreStateNum = 6 + clkOffNum + clkDriftNum + 1;
		int totalStateNum = coreStateNum + n + ionoParamNum;

		ArrayList<String> uniqSatList = new ArrayList<String>(uniqSat);

		SimpleMatrix x = pos_kfObj.getState();
		SimpleMatrix P = pos_kfObj.getCovariance();
		SimpleMatrix _x = new SimpleMatrix(totalStateNum, 1);
		SimpleMatrix _P = new SimpleMatrix(totalStateNum, totalStateNum);
		_x.insertIntoThis(0, 0, x.extractMatrix(0, coreStateNum, 0, 1));
		_P.insertIntoThis(0, 0, P.extractMatrix(0, coreStateNum, 0, coreStateNum));
		for (int j = coreStateNum; j < totalStateNum - ionoParamNum; j++) {
			Satellite sat = satList.get(j - coreStateNum);
			CycleSlipDetect csdObj = csdList.get(j - coreStateNum);
			String satObsvCode = sat.getObsvCode() + "" + sat.getSVID();
			boolean flag1 = csdObj.isCS() == false;
			boolean flag2 = csdObj.isCS() == true && csdObj.isRepaired() == true && csdObj.isReset() == false;
			boolean flag = flag1 || flag2;
			if (ambMap.containsKey(satObsvCode) && flag) {
				int k = ambMap.get(satObsvCode);
				double repairedCSval = 0;
				double repairedCSVar = 0;
				if (flag2) {
					repairedCSval = csdObj.getIntAmb();
					repairedCSVar = csdObj.getIntAmbCov();
				}
				_x.set(j, x.get(k) + repairedCSval);
				_P.set(j, j, P.get(k, k) + repairedCSVar);
				for (int l = 0; l < coreStateNum; l++) {
					_P.set(l, j, P.get(l, k));
					_P.set(j, l, P.get(k, l));
				}
				for (int l = j + 1; l < coreStateNum + n; l++) {
					Satellite _sat = satList.get(l - coreStateNum);
					String _satObsvCode = _sat.getObsvCode() + "" + _sat.getSVID();
					CycleSlipDetect _csdObj = csdList.get(l - coreStateNum);

					boolean _flag1 = _csdObj.isCS() == false;
					boolean _flag2 = _csdObj.isCS() == true && _csdObj.isRepaired() == true&&_csdObj.isReset() == false;;
					boolean _flag = _flag1 || _flag2;
					if (ambMap.containsKey(_satObsvCode) && _flag) {
						int _k = ambMap.get(_satObsvCode);
						_P.set(j, l, P.get(k, _k));
						_P.set(l, j, P.get(_k, k));

					}
				}
				for (int l = coreStateNum + n; l < totalStateNum; l++) {

					String satID = uniqSatList.get(l - (coreStateNum + n));
					if (ionoMap.containsKey(satID)) {
						int _k = ionoMap.get(satID);
						_P.set(j, l, P.get(k, _k));
						_P.set(l, j, P.get(_k, k));

					}
				}

			} else {
				
				double wl = SpeedofLight / sat.getCarrier_frequency();
				_x.set(j, (sat.getPhase() - (sat.getPseudorange())) / wl);
				_P.set(j, j, 1e16);

			}

			new_ambMap.put(satObsvCode, j);
			System.out.println(satObsvCode + " value: " + _x.get(j) + "  variance  " + _P.get(j, j));
		}
		System.out.println("\n");
		ambMap.clear();
		ambMap.putAll(new_ambMap);

		for (int j = coreStateNum + n; j < totalStateNum; j++) {

			String satID = uniqSatList.get(j - (coreStateNum + n));
			if (ionoMap.containsKey(satID)) {
				int k = ionoMap.get(satID);
				_x.set(j, x.get(k));
				_P.set(j, j, P.get(k, k));
				for (int l = 0; l < coreStateNum; l++) {
					_P.set(l, j, P.get(l, k));
					_P.set(j, l, P.get(k, l));
				}
				for (int l = j + 1; l < totalStateNum; l++) {
					String _satID = uniqSatList.get(l - (coreStateNum + n));
					if (ionoMap.containsKey(_satID)) {
						int _k = ionoMap.get(_satID);
						_P.set(j, l, P.get(k, _k));
						_P.set(l, j, P.get(_k, k));

					}
				}

			} else {

				_x.set(j, 0);
				_P.set(j, j, 1e6);

			}

			new_ionoMap.put(satID, j);
		}
		ionoMap.clear();
		ionoMap.putAll(new_ionoMap);
		pos_kfObj.setState_ProcessCov(_x, _P);
		// Assign Q and F matrix
		pos_kfObj.configPPP(deltaT, clkOffNum, clkDriftNum, totalStateNum, ionoParams, false,predictPhaseClock,false);
		pos_kfObj.predict();
		x = pos_kfObj.getState();
		SimpleMatrix R = new SimpleMatrix((3 * n) + ionoParamNum, (3 * n) + ionoParamNum);
		SimpleMatrix Cyy_pseudorange = SimpleMatrix.identity(n).scale(GnssDataConfig.pseudorange_priorVarOfUnitW);// Weight.igs_getNormCyy(satList,GnssDataConfig.pseudorange_priorVarOfUnitW);
		SimpleMatrix Cyy_phase = SimpleMatrix.identity(n).scale(GnssDataConfig.phase_priorVarOfUnitW);
		SimpleMatrix Cyy_doppler = SimpleMatrix.identity(n).scale(GnssDataConfig.doppler_priorVarOfUnitW);// Weight.igs_getNormCyy(satList,GnssDataConfig.doppler_priorVarOfUnitW);
		SimpleMatrix Cyy_GIM_iono = SimpleMatrix.identity(ionoParamNum).scale(GnssDataConfig.GIM_TECU_variance);
		R.insertIntoThis(0, 0, Cyy_pseudorange);
		R.insertIntoThis(n, n, Cyy_phase);
		R.insertIntoThis(2 * n, 2 * n, Cyy_doppler);
		R.insertIntoThis(3 * n, 3 * n, Cyy_GIM_iono);
		Object[] z_ze_H = get_z_ze_H(x, coreStateNum, clkOffNum, clkDriftNum, ionoParamNum, n, totalStateNum, satList,
				obsvCodeList, rxPCO, ionoParams, csdList, uniqSatList, ssiSet, false, currentTime);
		SimpleMatrix z = (SimpleMatrix) z_ze_H[0];
		SimpleMatrix ze = (SimpleMatrix) z_ze_H[1];
		SimpleMatrix H = (SimpleMatrix) z_ze_H[2];
		SimpleMatrix e_prior_hat = z.minus(ze);
		if (doAnalyze) {

			double[] pseudorange_e_prior_hat = Matrix.matrix2ArrayVec(e_prior_hat.extractMatrix(0, n, 0, 1));
			SimpleMatrix  _phase_e_prior_hat= e_prior_hat.extractMatrix(n, 2 * n, 0, 1);
			for (int i = 0; i < n; i++) {

				Satellite sat = satList.get(i);
				double wl = sat.getCarrier_wavelength();
				_phase_e_prior_hat.set(i,_phase_e_prior_hat.get(i)/wl);
			}
			double[] phase_e_prior_hat = Matrix.matrix2ArrayVec(_phase_e_prior_hat);
			double[] doppler_e_prior_hat = Matrix.matrix2ArrayVec(e_prior_hat.extractMatrix(2 * n, 3 * n, 0, 1));
			double[] iono_e_prior_hat = Matrix
					.matrix2ArrayVec(e_prior_hat.extractMatrix(3 * n, (3 * n) + ionoParamNum, 0, 1));
			Map<Measurement, double[]> measMap = Map.of(Measurement.Pseudorange, pseudorange_e_prior_hat,
					Measurement.CarrierPhase, phase_e_prior_hat, Measurement.Doppler, doppler_e_prior_hat,
					Measurement.GIM_Iono, iono_e_prior_hat);
			innovationMap.put(currentTime, measMap);
		}
		if (doTest) {
			P = pos_kfObj.getCovariance();
			SimpleMatrix R_pseudorange = performObservableTesting(e_prior_hat.extractMatrix(0, n, 0, 1), P,
					R.extractMatrix(0, n, 0, n), H.extractMatrix(0, n, 0, totalStateNum), n);
			R.insertIntoThis(0, 0, R_pseudorange);
			SimpleMatrix R_doppler = performObservableTesting(e_prior_hat.extractMatrix(2 * n, 3 * n, 0, 1), P,

					R.extractMatrix(2 * n, 3 * n, 2 * n, 3 * n), H.extractMatrix(2 * n, 3 * n, 0, totalStateNum), n);
			R.insertIntoThis(2 * n, 2 * n, R_doppler);
		}
		pos_kfObj.update(z, R, ze, H);
		if (doAnalyze) {
			x = pos_kfObj.getState();
//			P = pos_kfObj.getCovariance();
//			System.out.println("Design Matrix");
//			System.out.println(H.toString());
//			SimpleMatrix mainVar = P.diag();
//			System.out.println("Position Variance + " + mainVar.extractMatrix(0, 3, 0, 1));
//			System.out.println("Clock Offsets Variance + " + mainVar.extractMatrix(3, 3 + clkOffNum, 0, 1));
//			System.out.println("Velocity Variance + " + mainVar.extractMatrix(3 + clkOffNum, 3 + clkOffNum + 3, 0, 1));
//			System.out.println("Clock Drifts Variance + " + mainVar.extractMatrix(3 + clkOffNum + 3, 3 + clkOffNum + 3 + clkDriftNum, 0, 1));
//			System.out.println("Troposphere Variance + " + mainVar.extractMatrix(3 + clkOffNum + 3 + clkDriftNum, 3 + clkOffNum + 3 + clkDriftNum + 1, 0, 1));
//			System.out.println("Ambiguities Variance + " + mainVar.extractMatrix(coreStateNum, coreStateNum + n, 0, 1));
//			System.out.println("Ionosphere Variance + " + mainVar.extractMatrix(coreStateNum + n, totalStateNum, 0, 1));
			SimpleMatrix RotMat = new SimpleMatrix(3 + clkOffNum, 3 + clkOffNum);
			RotMat.insertIntoThis(0, 0, new SimpleMatrix(
					LatLonUtil.getEcef2EnuRotMat(Matrix.matrix2ArrayVec(x.extractMatrix(0, 3, 0, 1)))));
			for (int i = 0; i < clkOffNum; i++) {
				RotMat.set(3 + i, 3 + i, 1);
			}
			SimpleMatrix H_dop = H.extractMatrix(0, n, 0, 3 + clkOffNum);
			SimpleMatrix _dop = RotMat.mult((H_dop.transpose().mult(H_dop)).invert()).mult(RotMat.transpose());
			double[] dop = new double[] { _dop.get(0, 0), _dop.get(1, 1), _dop.get(2, 2), _dop.get(3, 3) };
			dopMap.put(currentTime, dop);
			z_ze_H = get_z_ze_H(x, coreStateNum, clkOffNum, clkDriftNum, ionoParamNum, n, totalStateNum, satList,
					obsvCodeList, rxPCO, ionoParams, csdList, uniqSatList, ssiSet, true, currentTime);
			z = (SimpleMatrix) z_ze_H[0];
			ze = (SimpleMatrix) z_ze_H[1];
			performAnalysis(z, ze, satList, R, H, currentTime, pos_kfObj);

		}

		// Ambiguity fixing for position filter
		if (n > 0 && fixAmb) {
			x = pos_kfObj.getState();
			P = pos_kfObj.getCovariance();
			SimpleMatrix floatAmb = x.extractMatrix(coreStateNum, coreStateNum + n, 0, 1);
			SimpleMatrix floatAmbCov = P.extractMatrix(coreStateNum, coreStateNum + n, coreStateNum, coreStateNum + n);
			System.out.println("Float Ambiguity");
			System.out.println(floatAmb.toString());
			System.out.println("Float Ambiguity Covariance");
			System.out.println(floatAmbCov.toString());

			// ADD BSD HERE

			// Group satellites by system
			HashMap<String, ArrayList<Integer>> sysGroups = new HashMap<>();
			for (int i = 0; i < n; i++) {
				String groupKey = satList.get(i).getObsvCode(); // e.g., "G1C" for system + freq
				sysGroups.computeIfAbsent(groupKey, k -> new ArrayList<>()).add(i);
			}
			boolean skipAR = false;
			// Select reference satellite per system based on elevation and CN0
			HashMap<String, Integer> refs = new HashMap<>();
			for (Map.Entry<String, ArrayList<Integer>> entry : sysGroups.entrySet()) {
				String sys = entry.getKey();
				ArrayList<Integer> idxs = entry.getValue();
				int refIdx = -1;
				double maxWeight = -1.0;
				for (int idx : idxs) {
					CycleSlipDetect csdObj = csdList.get(idx);
					if (csdObj.isCS()) {
						continue;
					}
					Satellite sat = satList.get(idx);
					double weight = 1.0 / Weight.computeCoVariance(sat.getCNo(), sat.getElevAzm()[0]);
					if (weight > maxWeight) {
						maxWeight = weight;
						refIdx = idx;
					}
				}
				if (refIdx != -1) {
					refs.put(sys, refIdx);
				} else {
					skipAR = true;
					break;
				}
			}

			int sd_n = 0;
			for (Map.Entry<String, ArrayList<Integer>> entry : sysGroups.entrySet()) {
				String sys = entry.getKey();
				if (refs.containsKey(sys)) {
					sd_n += entry.getValue().size() - 1;
				}
			}
			if (sd_n <= 0 || skipAR) {
				// No AR possible, skip or handle as float
				System.out.println("Insufficient satellites for BSD AR");
			} else {
				// Build transformation matrix Z (sd_n x n)
				SimpleMatrix Z = new SimpleMatrix(sd_n, n);
				int row = 0;
				for (Map.Entry<String, ArrayList<Integer>> entry : sysGroups.entrySet()) {
					String sys = entry.getKey();
					if (!refs.containsKey(sys)) {
						continue;
					}
					int ref = refs.get(sys);
					ArrayList<Integer> idxs = entry.getValue();
					for (int idx : idxs) {
						if (idx == ref) {
							continue;
						}
						Z.set(row, idx, 1.0);
						Z.set(row, ref, -1.0);
						row++;
					}
				}

				// Compute SD float ambiguities and covariance
				SimpleMatrix sd_a = Z.mult(floatAmb);
				SimpleMatrix sd_Q = Z.mult(floatAmbCov).mult(Z.transpose());
				sd_Q = sd_Q.plus(sd_Q.transpose()).scale(0.5);

				// Perform LAMBDA on SD
				boolean estimateVar = false;
				LambdaResult lmd = LAMBDA.computeLambda(sd_a, sd_Q, EstimatorType.PAR, estimateVar, 1, 0.99999);
				int nFixed = lmd.getnFixed();
				double Ps = lmd.getSr();
				SimpleMatrix sd_fix = lmd.getaFix();
				SimpleMatrix sd_Qfix = lmd.getQaFix();
				System.out.println(" Failure Rate : " + (1 - Ps));
				if (nFixed != 0) {
					System.out.println(" N Fixed (SD): " + nFixed);
					// Compute pseudoinverse of SD cov for stability (handles rank issues)
					SimpleMatrix sd_Q_inv = sd_Q.pseudoInverse(); // Use EJML's pseudoInverse()
					// Optional regularization if sd_Q is ill-conditioned
					double epsilon = 1e-8;
					sd_Q_inv = (sd_Q.plus(SimpleMatrix.identity(sd_n).scale(epsilon))).invert();
					// Recovery: UD fixed mean adjustment
					SimpleMatrix delta_sd = sd_fix.minus(Z.mult(floatAmb));
					SimpleMatrix adjustment = floatAmbCov.mult(Z.transpose()).mult(sd_Q_inv).mult(delta_sd);
					SimpleMatrix a_caron = floatAmb.plus(adjustment);
					// Recovery: UD fixed covariance (conditional variance)
					SimpleMatrix tempMat = floatAmbCov.mult(Z.transpose()).mult(sd_Q_inv).mult(Z).mult(floatAmbCov);
					SimpleMatrix Qa_caron = floatAmbCov.minus(tempMat); // Base conditional cov
					// Add back fixed SD cov contribution (for partial fixes)
					SimpleMatrix fixed_contrib = floatAmbCov.mult(Z.transpose()).mult(sd_Q_inv).mult(sd_Qfix)
							.mult(sd_Q_inv).mult(Z).mult(floatAmbCov);
					Qa_caron = Qa_caron.plus(fixed_contrib);
					Qa_caron = Qa_caron.plus(Qa_caron.transpose()).scale(0.5); // Symmetrize

					// Now use a_caron and Qa_caron in the conditional update
					SimpleMatrix a_hat = new SimpleMatrix(floatAmb);
					SimpleMatrix delta_a = a_hat.minus(a_caron);

					int bSize = coreStateNum;
					int aSize = n;
					int cSize = ionoParamNum;
					int ionoStart = coreStateNum + n;
					int ionoEnd = ionoStart + cSize;
					SimpleMatrix Cba = P.extractMatrix(0, bSize, coreStateNum, coreStateNum + aSize);
					SimpleMatrix Cca = P.extractMatrix(ionoStart, ionoEnd, coreStateNum, coreStateNum + aSize);
					SimpleMatrix Cbc = P.extractMatrix(0, bSize, ionoStart, ionoEnd);
					SimpleMatrix Cbb_hat = P.extractMatrix(0, bSize, 0, bSize);
					SimpleMatrix Ccc_hat = P.extractMatrix(ionoStart, ionoEnd, ionoStart, ionoEnd);
					SimpleMatrix b_hat = x.extractMatrix(0, bSize, 0, 1);
					SimpleMatrix c_hat = x.extractMatrix(ionoStart, ionoEnd, 0, 1);
					SimpleMatrix Caa_hat_inv = floatAmbCov.invert();
					SimpleMatrix b_caron = b_hat.minus(Cba.mult(Caa_hat_inv).mult(delta_a));
					SimpleMatrix c_caron = c_hat.minus(Cca.mult(Caa_hat_inv).mult(delta_a));
					SimpleMatrix temp = Caa_hat_inv.mult(Qa_caron).mult(Caa_hat_inv);
					SimpleMatrix Pbb_c = Cbb_hat.minus(Cba.mult(Caa_hat_inv).mult(Cba.transpose()))
							.plus(Cba.mult(temp).mult(Cba.transpose()));
					SimpleMatrix Pcc_c = Ccc_hat.minus(Cca.mult(Caa_hat_inv).mult(Cca.transpose()))
							.plus(Cca.mult(temp).mult(Cca.transpose()));
					SimpleMatrix Pbc_c = Cbc.minus(Cba.mult(Caa_hat_inv).mult(Cca.transpose()))
							.plus(Cba.mult(temp).mult(Cca.transpose()));
					SimpleMatrix Pba_c = Cba.mult(Caa_hat_inv).mult(Qa_caron);
					SimpleMatrix Pca_c = Cca.mult(Caa_hat_inv).mult(Qa_caron);
					SimpleMatrix x_new = new SimpleMatrix(totalStateNum, 1);
					x_new.insertIntoThis(0, 0, b_caron);
					x_new.insertIntoThis(coreStateNum, 0, a_caron);
					x_new.insertIntoThis(ionoStart, 0, c_caron);
					SimpleMatrix P_new = new SimpleMatrix(totalStateNum, totalStateNum);
					P_new.insertIntoThis(0, 0, Pbb_c);
					P_new.insertIntoThis(coreStateNum, coreStateNum, Qa_caron);
					P_new.insertIntoThis(ionoStart, ionoStart, Pcc_c);
					P_new.insertIntoThis(0, coreStateNum, Pba_c);
					P_new.insertIntoThis(coreStateNum, 0, Pba_c.transpose());
					P_new.insertIntoThis(0, ionoStart, Pbc_c);
					P_new.insertIntoThis(ionoStart, 0, Pbc_c.transpose());
					P_new.insertIntoThis(coreStateNum, ionoStart, Pca_c.transpose());
					P_new.insertIntoThis(ionoStart, coreStateNum, Pca_c);
					P_new = P_new.plus(P_new.transpose()).scale(0.5);
					pos_kfObj.setState_ProcessCov(x_new, P_new);
					System.out.println("Fixed Ambiguity Sequence");
					System.out.println(a_caron.toString());
					System.out.println("Fixed Ambiguity Variance");
					System.out.println(Qa_caron.toString());
					System.out.println(" N Fixed : " + nFixed);

					ambFixedCountMap.put(currentTime, nFixed);
					ambFixedCount += nFixed;
				}
			}
		}

	}

	private Object[] performTesting(SimpleMatrix R, SimpleMatrix H, int n, int m, ArrayList<Satellite> satList,
			ArrayList<Satellite> testedSatList, SimpleMatrix z, SimpleMatrix ze, ArrayList<CycleSlipDetect> csdList,
			boolean isTDCP, KFconfig kf_Obj) throws Exception {
		// Pre-fit residual/innovation
		SimpleMatrix v = new SimpleMatrix(n, 1, true, innovation);
		SimpleMatrix P = kf_Obj.getCovariance();
		SimpleMatrix Cvv = ((H.mult(P).mult(H.transpose())).plus(R));
		SimpleMatrix Cvv_inv = Cvv.invert();
		double globalTq = v.transpose().mult(Cvv_inv).mult(v).get(0);
		ChiSquaredDistribution csd = new ChiSquaredDistribution(n);
		double alpha = 0.01;
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
			SimpleMatrix H_ = new SimpleMatrix(_n, 3 + (2 * m));
			SimpleMatrix v_ = new SimpleMatrix(_n, 1);
			int j = 0;
			for (int i = 0; i < _n + 1; i++) {
				if (i != index) {
					R_.set(j, j, R.get(i, i));
					v_.set(j, v.get(i));
					z_.set(j, z.get(i));
					ze_.set(j, ze.get(i));
					for (int k = 0; k < 3 + (2 * m); k++) {
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

	private SimpleMatrix performObservableTesting(SimpleMatrix v, SimpleMatrix P, SimpleMatrix R, SimpleMatrix H,
			int n) throws Exception {
		// Pre-fit residual/innovation
		SimpleMatrix Ht = H.transpose();
		SimpleMatrix Cvv = ((H.mult(P).mult(Ht)).plus(R));
		SimpleMatrix Cvv_inv = Cvv.invert();
		SimpleMatrix K = P.mult(Ht).mult(Cvv_inv);
		SimpleMatrix HK = H.mult(K);
		int totalN = HK.numRows();
		double redunNo = SimpleMatrix.identity(totalN).minus(HK).trace();
		double globalTq = v.transpose().mult(Cvv_inv).mult(v).get(0);
		ChiSquaredDistribution csd = new ChiSquaredDistribution(redunNo);
		double alpha = 0.01;
		if (globalTq == 0) {
			throw new Exception("Error: T stat is zero");
		}
		// Detection
		double globalPVal = 1 - csd.cumulativeProbability(globalTq);
		while (globalPVal < alpha && redunNo > (n / 2.0)) {
			double max_w = Double.MIN_VALUE;
			int index = -1;
			for (int i = 0; i < n; i++) {
				SimpleMatrix cv = new SimpleMatrix(n, 1);
				cv.set(i, 1);
				double w = Math.abs(cv.transpose().mult(Cvv_inv).mult(v).get(0))
						/ Math.sqrt(cv.transpose().mult(Cvv_inv).mult(cv).get(0));
				if (w > max_w) {
					max_w = w;
					index = i;
				}
			}
			R.set(index, index, 1e16);
			Cvv = ((H.mult(P).mult(H.transpose())).plus(R));
			Cvv_inv = Cvv.invert();
			K = P.mult(Ht).mult(Cvv_inv);
			HK = H.mult(K);
			totalN = HK.numRows();
			redunNo = SimpleMatrix.identity(totalN).minus(HK).trace();
			globalTq = v.transpose().mult(Cvv_inv).mult(v).get(0);
			if (globalTq == 0) {
				throw new Exception("Error: T stat is zero");
			}
			csd = new ChiSquaredDistribution(redunNo);
			globalPVal = 1 - csd.cumulativeProbability(globalTq);
		}
		return R;

	}

	private Object[] get_z_ze_H(SimpleMatrix x, int coreStateNum, int clkOffNum, int clkDriftNum, int ionoParamNum,
			int n, int totalStateNum, ArrayList<Satellite> satList, String[] obsvCodeList,
			HashMap<String, double[]> rxPCO, ArrayList<double[]> ionoParams, ArrayList<CycleSlipDetect> csdList,
			ArrayList<String> uniqSatList, ListOrderedSet ssiSet, boolean doAnalyze, long currentTime) {
		double[] estPos = new double[] { x.get(0), x.get(1), x.get(2) };
		double[] rxClkOff = new double[clkOffNum];// in meters
		for (int i = 0; i < clkOffNum; i++) {
			rxClkOff[i] = x.get(i + 3);
		}
		double[] estVel = new double[] { x.get(3 + clkOffNum), x.get(4 + clkOffNum), x.get(5 + clkOffNum) };
		double[] rxClkDrift = new double[clkDriftNum];// in meters
		for (int i = 0; i < clkDriftNum; i++) {
			rxClkDrift[i] = x.get(6 + clkOffNum + i);
		}
		double estTropo = x.get(coreStateNum - 1);
		SimpleMatrix estAmb = x.extractMatrix(coreStateNum, coreStateNum + n, 0, 1);
		SimpleMatrix estIonoTec = x.extractMatrix(coreStateNum + n, totalStateNum, 0, 1);

		if (doAnalyze) {
			tropoMap.put(currentTime, estTropo);
			clkOffMap.put(currentTime, rxClkOff);
			clkDriftMap.put(currentTime, rxClkDrift);
		}

		SimpleMatrix H = new SimpleMatrix((3 * n) + ionoParamNum, totalStateNum);
		SimpleMatrix z = new SimpleMatrix((3 * n) + ionoParamNum, 1);
		SimpleMatrix ze = new SimpleMatrix((3 * n) + ionoParamNum, 1);
		SimpleMatrix unitLOS = new SimpleMatrix(SatUtil.igs_getUnitLOS(satList, estPos));
		H.insertIntoThis(0, 0, unitLOS.scale(-1));
		H.insertIntoThis(n, 0, unitLOS.scale(-1));
		H.insertIntoThis(2 * n, 3 + clkOffNum, unitLOS.scale(-1));
		HashMap<String, Double> _ionoMap = new HashMap<String, Double>();
		HashMap<String, Double> _ambMap = new HashMap<String, Double>();
		for (int i = 0; i < n; i++) {

			Satellite sat = satList.get(i);
			CycleSlipDetect csdObj = csdList.get(i);
			String obsvCode = sat.getObsvCode();
			char ssi = obsvCode.charAt(0);
			String satID = sat.getObsvCode().charAt(0) + "" + sat.getSVID();
			double wavelength = SpeedofLight / sat.getCarrier_frequency();
			double freq2 = Math.pow(sat.getCarrier_frequency(), 2);
			double ionoCoeff = (40.3 * 1e16) / freq2;
			int ionoIndex = uniqSatList.indexOf(satID);
			double[] pco = rxPCO.get(obsvCode);
			z.set(i, sat.getPseudorange());
			z.set(i + n, sat.getPhase());
			z.set(i + (2 * n), csdObj.getDopplerDR() - csdObj.getSatVelCorr());
			double geometricRange = Math
					.sqrt(IntStream.range(0, 3).mapToDouble(j -> estPos[j] + pco[j] - sat.getSatEci()[j])
							.map(j -> j * j).reduce(0, (j, k) -> j + k));

			double estPR = geometricRange + (sat.getWetMF() * estTropo) + (ionoCoeff * estIonoTec.get(ionoIndex));
			double estCP = geometricRange + (sat.getWetMF() * estTropo) - (ionoCoeff * estIonoTec.get(ionoIndex))
					+ (wavelength * estAmb.get(i));
			ze.set(i, estPR);
			ze.set(i + n, estCP);
			ze.set(i + (2 * n), H.extractMatrix(i, i + 1, 0, 3).mult(new SimpleMatrix(3, 1, true, estVel)).get(0));
			ze.set(i + n, ze.get(i + n) + rxClkOff[0]);
			H.set(i + n, 3, 1);
			for (int j = 0; j < clkOffNum; j++) {
				double clkOffVal = rxClkOff[0];
				if (obsvCode.equals(obsvCodeList[j])) {
					if (j == 0) {
						clkOffVal = 0;
					} else {
						H.set(i, 3, 1);
					}
					ze.set(i, ze.get(i) + rxClkOff[j] + clkOffVal);
					H.set(i, 3 + j, 1);
				}
			}
			for (int j = 0; j < clkDriftNum; j++) {
				if ((char) ssiSet.get(j) == ssi) {
					H.set(i + (2 * n), 6 + clkOffNum + j, 1);
					ze.set(i + (2 * n), ze.get(i + (2 * n)) + rxClkDrift[j]);
				}
			}
			H.set(i, coreStateNum - 1, sat.getWetMF());
			H.set(i + n, coreStateNum - 1, sat.getWetMF());

			H.set(i + n, coreStateNum + i, wavelength);

			H.set(i, coreStateNum + n + ionoIndex, ionoCoeff);
			H.set(i + n, coreStateNum + n + ionoIndex, -ionoCoeff);

			if (doAnalyze) {
				_ionoMap.put(satID, estIonoTec.get(ionoIndex));
				_ambMap.put(sat.getObsvCode() + "" + sat.getSVID(), estAmb.get(i));
			}

		}
		for (int i = 0; i < ionoParamNum; i++) {
			z.set(i + (3 * n), ionoParams.get(i)[1]);
			ze.set(i + (3 * n), x.get(coreStateNum + n + i));
			H.set(i + (3 * n), coreStateNum + n + i, 1);

		}
		if (doAnalyze) {
			ionoMap.put(currentTime, _ionoMap);
			ambMap.put(currentTime, _ambMap);
		}
		return new Object[] { z, ze, H };

	}

	private void performAnalysis(SimpleMatrix z, SimpleMatrix ze, ArrayList<Satellite> satList, SimpleMatrix R,
			SimpleMatrix H, long currentTime, KFconfig kf_Obj) {

		int n = satList.size();
		SimpleMatrix e_post_hat = z.minus(ze);
		// Post-fit residual
		SimpleMatrix Cyy_inv = R.invert();
		// Compute Redundancies
		SimpleMatrix K = kf_Obj.getKalmanGain();
		SimpleMatrix HK = H.mult(K);
		int totalN = HK.numRows();
		SimpleMatrix redunMat = SimpleMatrix.identity(totalN).minus(HK);
		redunMat = redunMat.plus(redunMat.transpose()).scale(0.5);
		double pseudorange_redundancyNo = redunMat.extractMatrix(0, n, 0, n).trace();
		double phase_redundancyNo = redunMat.extractMatrix(n, 2 * n, n, 2 * n).trace();
		double doppler_redundancyNo = redunMat.extractMatrix(2 * n, 3 * n, 2 * n, 3 * n).trace();
		double iono_redundancyNo = redunMat.extractMatrix(3 * n, totalN, 3 * n, totalN).trace();

		SimpleMatrix pseudorange_e_post_hat = e_post_hat.extractMatrix(0, n, 0, 1);
		SimpleMatrix phase_e_post_hat = e_post_hat.extractMatrix(n, 2 * n, 0, 1);
		SimpleMatrix doppler_e_post_hat = e_post_hat.extractMatrix(2 * n, 3 * n, 0, 1);
		SimpleMatrix iono_e_post_hat = e_post_hat.extractMatrix(3 * n, totalN, 0, 1);

		double pseudorange_postVarOfUnitW = pseudorange_e_post_hat.transpose().mult(Cyy_inv.extractMatrix(0, n, 0, n))
				.mult(pseudorange_e_post_hat).get(0) / pseudorange_redundancyNo;
		double phase_postVarOfUnitW = phase_e_post_hat.transpose().mult(Cyy_inv.extractMatrix(n, 2 * n, n, 2 * n))
				.mult(phase_e_post_hat).get(0) / phase_redundancyNo;
		double doppler_postVarOfUnitW = doppler_e_post_hat.transpose()
				.mult(Cyy_inv.extractMatrix(2 * n, 3 * n, 2 * n, 3 * n)).mult(doppler_e_post_hat).get(0)
				/ doppler_redundancyNo;
		double iono_postVarOfUnitW = iono_e_post_hat.transpose()
				.mult(Cyy_inv.extractMatrix(3 * n, totalN, 3 * n, totalN)).mult(iono_e_post_hat).get(0)
				/ iono_redundancyNo;

		Map<Measurement, double[]> measResMap = Map.of(Measurement.Pseudorange,
				Matrix.matrix2ArrayVec(pseudorange_e_post_hat), Measurement.CarrierPhase,
				Matrix.matrix2ArrayVec(phase_e_post_hat), Measurement.Doppler,
				Matrix.matrix2ArrayVec(doppler_e_post_hat), Measurement.GIM_Iono,
				Matrix.matrix2ArrayVec(iono_e_post_hat));
		residualMap.put(currentTime, measResMap);
		Map<Measurement, Double> measPostVarMap = Map.of(Measurement.Pseudorange, pseudorange_postVarOfUnitW,
				Measurement.CarrierPhase, phase_postVarOfUnitW, Measurement.Doppler, doppler_postVarOfUnitW,
				Measurement.GIM_Iono, iono_postVarOfUnitW);
		postVarOfUnitWMap.put(currentTime, measPostVarMap);
		Map<Measurement, Double> measRedunNoMap = Map.of(Measurement.Pseudorange, pseudorange_redundancyNo,
				Measurement.CarrierPhase, phase_redundancyNo, Measurement.Doppler, doppler_redundancyNo,
				Measurement.GIM_Iono, iono_redundancyNo);
		redundancyNoMap.put(currentTime, measRedunNoMap);
		satListMap.put(currentTime, satList);
		satCountMap.put(currentTime, (long) n);

	}

	public long getCsDetectedCount() {
		return csDetectedCount;
	}

	public long getCsRepairedCount() {
		return csRepairedCount;
	}

	public TreeMap<Long, Integer> getCsDetectedCountMap() {
		return csDetectedCountMap;
	}

	public TreeMap<Long, Integer> getCsRepairedCountMap() {
		return csRepairedCountMap;
	}

	public TreeMap<Long, ArrayList<CycleSlipDetect>> getCsdListMap() {
		return csdListMap;
	}

	public long getAmbFixedCount() {
		return ambFixedCount;
	}

	public TreeMap<Long, Integer> getAmbFixedCountMap() {
		return ambFixedCountMap;
	}

	public HashMap<String, int[]> getCycleSlipCount() {
		return cycleSlipCount;
	}

	public TreeMap getInnovationMap() {
		return innovationMap;
	}

	public TreeMap getResidualMap() {
		return residualMap;
	}

	public TreeMap getPostVarOfUnitWMap() {
		return postVarOfUnitWMap;
	}

	public TreeMap getSatListMap() {
		return satListMap;
	}

	public TreeMap<Long, HashMap<String, Double>> getAmbMap() {
		return ambMap;
	}

	public TreeMap<Long, HashMap<String, Double>> getIonoMap() {
		return ionoMap;
	}

	public TreeMap<Long, Double> getTropoMap() {
		return tropoMap;
	}

	public TreeMap<Long, double[]> getClkOffMap() {
		return clkOffMap;
	}

	public TreeMap<Long, double[]> getClkDriftMap() {
		return clkDriftMap;
	}

	public TreeMap<Long, Map<Measurement, Double>> getRedundancyNoMap() {
		return redundancyNoMap;
	}

	public TreeMap<Long, double[]> getDopMap() {
		return dopMap;
	}

}
