package com.gnssAug.Android.estimation.KalmanFilter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeMap;
import java.util.stream.IntStream;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.ejml.dense.row.MatrixFeatures_DDRM;
import org.ejml.simple.SimpleMatrix;
import org.hipparchus.linear.AbstractRealMatrix;
import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.RealMatrix;
import org.orekit.estimation.measurements.gnss.IntegerLeastSquareSolution;
import org.orekit.estimation.measurements.gnss.LambdaMethod;
import org.orekit.estimation.measurements.gnss.SimpleRatioAmbiguityAcceptance;

import com.gnssAug.Android.constants.GnssDataConfig;
import com.gnssAug.Android.constants.Measurement;
import com.gnssAug.Android.estimation.LinearLeastSquare;
import com.gnssAug.Android.estimation.KalmanFilter.Models.KFconfig;
import com.gnssAug.Android.models.CycleSlipDetect;
import com.gnssAug.Android.models.Satellite;
import com.gnssAug.helper.IntegerLeastSquares;
import com.gnssAug.utility.Combination;
import com.gnssAug.utility.LatLonUtil;
import com.gnssAug.utility.Matrix;
import com.gnssAug.utility.SatUtil;
import com.gnssAug.utility.Vector;
import com.gnssAug.utility.Weight;

public class EKF_PPP extends KFDopplerParent {

	private final double SpeedofLight = 299792458;
	private final double ionoErrCoeff = 40.3 * (1E16) * 20;
	public EKF_PPP() {
		kfObj = new KFconfig();
	}

	public TreeMap<Long, double[]> process(TreeMap<Long, ArrayList<Satellite>> SatMap, ArrayList<Long> timeList,
			boolean useIGS, String[] obsvCodeList, boolean doAnalyze, boolean doTest, boolean outlierAnalyze)
			throws Exception {
		boolean isWeighted = true;
		int m = obsvCodeList.length;
		int n = 3 + m;
		double[][] _X = new double[n][1];
		double[][] P = new double[n][n];
		/*
		 * state XYZ is intialized WLS generated Rx position estimated using first epoch
		 * data, a-priori estimate error covariance matrix for state XYZ is therefore
		 * assigned 25 m^2 value. Other state variables are assigned infinite(big)
		 * variance
		 */
		double[] intialECEF = LinearLeastSquare.getEstPos(SatMap.firstEntry().getValue(), true, useIGS);
		IntStream.range(0, n).forEach(i -> _X[i][0] = intialECEF[i]);

		// Total State
		SimpleMatrix X = new SimpleMatrix(_X);
		IntStream.range(0, n).forEach(i -> P[i][i] = 100);
		// Error State intialized as zero
		double[][] x = new double[n][1];
		kfObj.setState_ProcessCov(x, P);
		if (doAnalyze) {
			innovationMap = new TreeMap<Long, double[]>();
			errCovMap = new TreeMap<Long, SimpleMatrix>();
			residualMap = new TreeMap<Long, double[]>();
			postVarOfUnitWMap = new TreeMap<Long, Double>();
			redundancyList = new ArrayList<double[]>();
			satCountMap = new TreeMap<Long, Long>();
			satListMap = new TreeMap<Long, ArrayList<Satellite>>();
			measNoiseMap = new TreeMap<Long, double[]>();
		}
		return iterate(X, SatMap, timeList, useIGS, obsvCodeList, doAnalyze, doTest, outlierAnalyze,isWeighted);

	}

	TreeMap<Long, double[]> iterate(SimpleMatrix X, TreeMap<Long, ArrayList<Satellite>> SatMap,
			ArrayList<Long> timeList, boolean useIGS, String[] obsvCodeList, boolean doAnalyze, boolean doTest,
			boolean outlierAnalyze,boolean isWeighted) throws Exception {
		TreeMap<Long, double[]> estStateMap = new TreeMap<Long, double[]>();
		int m = obsvCodeList.length;
		int x_size = 3 + m;
		long time = timeList.get(0);
		HashMap<String, Integer> satBookKeep = new HashMap<String, Integer>();
		prevVel = LinearLeastSquare.getEstVel(SatMap.get(time), isWeighted, true, doTest, false,
				new double[] { X.get(0), X.get(1), X.get(2) }, useIGS);
		prev_Cxx_dot_hat = LinearLeastSquare.getCxx_hat(Measurement.Doppler, "ECEF");
		// Start from 2nd epoch
		for (int i = 1; i < timeList.size(); i++) {
			long currentTime = timeList.get(i);
			double deltaT = (currentTime - time) / 1e3;
			ArrayList<Satellite> satList = SatMap.get(currentTime);
			//CSDR2(SatMap.get(time), satList, useIGS);
			CSDR_Orekit(SatMap.get(time), satList, useIGS);
//			int n = satList.size(); 
//			SimpleMatrix Cxx_dot_hat = predictTotalState(X, satList, deltaT, useIGS, doTest,isWeighted,obsvCodeList);
//			SimpleMatrix x = kfObj.getState();
//			SimpleMatrix P = kfObj.getCovariance();
//			HashMap<Integer, Integer> newSatBookKeep = new HashMap<Integer, Integer>();
//			SimpleMatrix _x = new SimpleMatrix(3 + m + (2*n), 1);
//			SimpleMatrix _P = new SimpleMatrix(3 + m + (2*n), 3 + m + (2*n));
//			_x.insertIntoThis(0, 0, x.extractMatrix(0, 3+m, 0, 1));
//			_P.insertIntoThis(0, 0, P.extractMatrix(0, 3+m, 0, 3+m));
//			for (int j = 3 + m + 0; j < 3 + m + n; j++) {
//				Satellite sat = satList.get(j - 3);
//				String obsvCode = sat.getObsvCode()+sat.getSvid();
//				if (satBookKeep.containsKey(obsvCode) && sat.isPhaseLocked() == true) {
//					int k = map.get(prn);
//					_x.set(j, x.get(k));
//					_P.set(0, j, P.get(0, k));
//					_P.set(1, j, P.get(1, k));
//					_P.set(2, j, P.get(2, k));
//					_P.set(j, 0, P.get(k, 0));
//					_P.set(j, 1, P.get(k, 1));
//					_P.set(j, 2, P.get(k, 2));
//					_P.set(j, j, P.get(k, k));
//					for (int l = j + 1; l < 3 + n; l++) {
//						Observation _obs = baselineObsvs.get(l - 3);
//						int _prn = _obs.getPrn();
//						if (map.containsKey(_prn) && _obs.isPhaseLocked() == true) {
//							int m = map.get(_prn);
//							_P.set(j, l, P.get(k, m));
//							_P.set(l, j, P.get(m, k));
//
//						}
//					}
//
//				} else {
//					_x.set(j, (obs.getPhaseL1() - obs.getPseudorange()) / wavelengthL1);
//					_P.set(j, j, 1e8);
//
//				}
//
//				newMap.put(prn, j);
//			}
//			
//			
//			
//			runFilter(X, currentTime, deltaT, satList, obsvCodeList, Cxx_dot_hat, doAnalyze, doTest, outlierAnalyze,
//					useIGS, i,isWeighted);
//			SimpleMatrix P = kfObj.getCovariance();
//			double[] estState = new double[x_size];
//			IntStream.range(0, x_size).forEach(j -> estState[j] = X.get(j));
//			// Add position estimate to the list
//			estStateMap.put(currentTime, estState);
//			if (doAnalyze) {
//				// Convert to ENU frame
//				SimpleMatrix R = new SimpleMatrix(x_size, x_size);
//				R.insertIntoThis(0, 0, new SimpleMatrix(LatLonUtil.getEcef2EnuRotMat(estState)));
//				IntStream.range(3, x_size).forEach(j -> R.set(j, j, 1));
//
//				SimpleMatrix errCov = R.mult(P).mult(R.transpose());
//				errCovMap.put(currentTime, errCov);
//				innovationMap.put(currentTime, innovation);
//			}
//			if (!MatrixFeatures_DDRM.isPositiveSemidefinite(P.getMatrix())) {
//				throw new Exception("PositiveDefinite test Failed");
//			}
			time = currentTime;
		}
		return estStateMap;
	}

	private SimpleMatrix predictTotalState(SimpleMatrix X, ArrayList<Satellite> satList, double deltaT, boolean useIGS,
			boolean doTest,boolean isWeighted,String[] obsvCodeList) throws Exception {
		double[] vel = LinearLeastSquare.getEstVel(SatUtil.createCopy(satList), isWeighted, true, doTest, false,
				new double[] { X.get(0), X.get(1), X.get(2) }, useIGS);
		SimpleMatrix Cxx_dot_hat = LinearLeastSquare.getCxx_hat(Measurement.Doppler, "ECEF");
		Object[] resettedVar = SatUtil.resetVar(Measurement.Doppler, obsvCodeList, vel, Cxx_dot_hat);
		vel = (double[]) resettedVar[0];
		Cxx_dot_hat = (SimpleMatrix) resettedVar[1];
		double[] avg_vel = new double[vel.length];
		for (int i = 0; i < vel.length; i++) {
			avg_vel[i] = (vel[i] + prevVel[i]) * 0.5;
			X.set(i, X.get(i) + (prevVel[i] * deltaT));
		}
		
		SimpleMatrix avg_Cxx_dot_hat = Cxx_dot_hat.plus(prev_Cxx_dot_hat).scale(0.5);
		prevVel = Arrays.copyOf(vel, vel.length);
		SimpleMatrix temp = prev_Cxx_dot_hat;
		prev_Cxx_dot_hat = new SimpleMatrix(Cxx_dot_hat);
		return temp;
	}

	// Innovation Based Testing
	private void runFilter(SimpleMatrix X, long currentTime, double deltaT, ArrayList<Satellite> satList,
			String[] obsvCodeList, SimpleMatrix Cxx_dot_hat, boolean doAnalyze, boolean doTest, boolean outlierAnalyze,
			boolean useIGS, int ct,boolean isWeighted) throws Exception {

		boolean useAndroidW = false;
		// Satellite count
		int n = satList.size();
		int m = obsvCodeList.length;

		// Last update VC-matrix
		SimpleMatrix priorP = new SimpleMatrix(kfObj.getCovariance());
		SimpleMatrix priorX = new SimpleMatrix(X);
		// Assign Q and F matrix
		kfObj.configDoppler(deltaT, Cxx_dot_hat, m, X);
		kfObj.predict();

		double[] estPos = new double[] { X.get(0), X.get(1), X.get(2) };
		double[] rxClkOff = new double[m];// in meters
		for (int i = 0; i < m; i++) {
			rxClkOff[i] = X.get(i + 3);
		}

		/*
		 * H is the Jacobian matrix of partial derivatives Observation StateModel(h) of
		 * with respect to x
		 */
		double[][] _H = new double[n][3 + m];
		// Measurement vector
		double[][] z = new double[n][1];
		// Estimated Measurement vector
		double[][] ze = new double[n][1];
		// Measurement Noise
		double[][] _R = new double[n][n];
		innovation = new double[n];
		temp_innovation = new double[n];

		if (isWeighted) {
			if (useAndroidW) {
				for (int i = 0; i < n; i++) {
					_R[i][i] = Math.pow(satList.get(i).getReceivedSvTimeUncertaintyNanos() * SpeedofLight * 1e-9, 2);
				}
			} else {
//				LinearLeastSquare.getEstPos(satList, true, true, false, false, useIGS);
//				SimpleMatrix Cyy = LinearLeastSquare.getCyy_updated(Measurement.Pseudorange);
//				_R = Matrix.matrix2Array(Cyy);
				SimpleMatrix Cyy = Weight.getNormCyy(satList, GnssDataConfig.pseudorange_priorVarOfUnitW);
				_R = Matrix.matrix2Array(Cyy);
			}
		} else {
			for (int i = 0; i < n; i++) {
				_R[i][i] = GnssDataConfig.pseudorange_priorVarOfUnitW;
			}
		}
		/*
		 * Notice measurement 'z' is difference b/w PR and estimated PR, because it's a
		 * complimentary filter predicted error state is zero and therefore Hx = 0 or
		 * estimated measurement(ze) is zero
		 */
		for (int i = 0; i < n; i++) {
			Satellite sat = satList.get(i);
			// Its not really a ECI, therefore don't get confused
			String obsvCode = sat.getObsvCode();
			double[] satEcef = sat.getSatEci();
			double PR = sat.getPseudorange();
			double[] LOS = IntStream.range(0, 3).mapToDouble(j -> satEcef[j] - estPos[j]).toArray();
			// Approx Geometric Range
			double approxGR = Vector.mod(LOS);
			// Approx Pseudorange Range
			double approxPR = approxGR;
			for (int j = 0; j < m; j++) {
				if (obsvCode.equals(obsvCodeList[j])) {
					approxPR += rxClkOff[j];
					_H[i][3 + j] = 1;
				}
			}

			z[i][0] = PR - approxPR;
			final int _i = i;
			IntStream.range(0, 3).forEach(j -> _H[_i][j] = -LOS[j] / approxGR);
			innovation[i] = z[i][0] - ze[i][0];
		}
		SimpleMatrix R = new SimpleMatrix(_R);
		SimpleMatrix H = new SimpleMatrix(_H);
		ArrayList<Satellite> testedSatList = new ArrayList<Satellite>(satList);
		if (doTest && !outlierAnalyze) {
			Object[] params = performTesting(R, H, n, m, satList, testedSatList, z, ze, ct);
			R = (SimpleMatrix) params[0];
			H = (SimpleMatrix) params[1];
			z = (double[][]) params[2];
			ze = (double[][]) params[3];
		}

		kfObj.update(z, R, ze, H);
		SimpleMatrix x = kfObj.getState();
		for (int i = 0; i < 3 + m; i++) {
			X.set(i, X.get(i) + x.get(i));
			x.set(i, 0);
		}
		if (doAnalyze) {
			performAnalysis(X, testedSatList, satList, R, H, priorP, currentTime, n, obsvCodeList, doTest,
					outlierAnalyze, ct);
		}
		if (doTest && outlierAnalyze) {
			kfObj.setState_ProcessCov(new SimpleMatrix(3 + m, 1), priorP);
			for (int i = 0; i < 3 + m; i++) {
				X.set(i, priorX.get(i));
			}
			kfObj.predict();
			Object[] params = performTesting(R, H, n, m, satList, testedSatList, z, ze, ct);
			R = (SimpleMatrix) params[0];
			H = (SimpleMatrix) params[1];
			z = (double[][]) params[2];
			ze = (double[][]) params[3];
			kfObj.update(z, R, ze, H);
			x = kfObj.getState();
			for (int i = 0; i < 3 + m; i++) {
				X.set(i, X.get(i) + x.get(i));
				x.set(i, 0);
			}
		}
	}

	private void performAnalysis(SimpleMatrix X, ArrayList<Satellite> testedSatList, ArrayList<Satellite> satList,
			SimpleMatrix R, SimpleMatrix H, SimpleMatrix priorP, long currentTime, int n, String[] obsvCodeList,
			boolean doTest, boolean outlierAnalyze, int ct) {

		int _n = testedSatList.size();
		double[] measNoise = new double[_n];
		double[] residual = (double[]) get_z_ze_res(X, testedSatList, obsvCodeList)[2];
		// Post-fit residual
		SimpleMatrix e_post_hat = new SimpleMatrix(_n, 1, true, residual);
		SimpleMatrix Cyy_inv = R.invert();

		// Compute Redundancies
		SimpleMatrix K = kfObj.getKalmanGain();
		SimpleMatrix HK = H.mult(K);
		SimpleMatrix phi = kfObj.getPhi();
		SimpleMatrix Cvv = kfObj.getCvv();
//		if (ct > 1994&&ct<2000) {
//			System.out.println("Cvv epoch " + ct);
//			System.out.println(Cvv);
//		}
		SimpleMatrix HtCvvInvH = H.transpose().mult(Cvv.invert()).mult(H);
		SimpleMatrix Q = kfObj.getQ();
		double rX = phi.mult(priorP).mult(phi.transpose()).mult(HtCvvInvH).trace();
		double rW = Q.mult(HtCvvInvH).trace();
		double rZ = SimpleMatrix.identity(HK.numRows()).minus(HK).trace();
		double rSum = rX + rW + rZ;
		if (_n - rSum > 0.01) {
			System.err.println("FATAL ERROR: Redundancy sum is wrong");
		}
		double postVarOfUnitW = e_post_hat.transpose().mult(Cyy_inv).mult(e_post_hat).get(0) / rZ;
		redundancyList.add(new double[] { _n, rSum, rX, rW, rZ });
		postVarOfUnitWMap.put(currentTime, postVarOfUnitW);
		residualMap.put(currentTime, residual);
		if (doTest == true) {
			_n = n - _n;
		}
		satCountMap.put(currentTime, (long) _n);
		if (outlierAnalyze) {
			satListMap.put(currentTime, satList);
		} else {
			satListMap.put(currentTime, testedSatList);
		}
		measNoiseMap.put(currentTime, measNoise);

	}

	private Object[] performTesting(SimpleMatrix R, SimpleMatrix H, int n, int m, ArrayList<Satellite> satList,
			ArrayList<Satellite> testedSatList, double[][] z, double[][] ze, int ct) throws Exception {
		// Pre-fit residual/innovation
		SimpleMatrix v = new SimpleMatrix(n, 1, true, innovation);
		SimpleMatrix P = kfObj.getCovariance();
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
//		if(globalPVal < alpha)
//		{
//			System.out.println("\nOutlier Detection Epoch "+ct);
//			System.out.println("Size "+_n);
//		}
		
		
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
			//System.out.print(" " + index + " ");
			satList.get(satList.indexOf(testedSatList.remove(index))).setOutlier(true);
			_n = testedSatList.size();
			SimpleMatrix R_ = new SimpleMatrix(_n, _n);
			double[][] z_ = new double[_n][1];
			double[][] ze_ = new double[_n][1];
			SimpleMatrix H_ = new SimpleMatrix(_n, 3 + m);
			SimpleMatrix v_ = new SimpleMatrix(_n, 1);
			int j = 0;
			for (int i = 0; i < _n + 1; i++) {
				if (i != index) {
					R_.set(j, j, R.get(i, i));
					v_.set(j, v.get(i));
					z_[j][0] = z[i][0];
					ze_[j][0] = ze[i][0];
					for (int k = 0; k < 3 + m; k++) {
						H_.set(j, k, H.get(i, k));
					}
					j++;
				}
			}
			R = new SimpleMatrix(R_);
			H = new SimpleMatrix(H_);
			z = new double[_n][1];
			ze = new double[_n][1];
			v = new SimpleMatrix(v_);
			for (int i = 0; i < _n; i++) {
				z[i] = z_[i];
				ze[i] = ze_[i];
			}
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

	private Object[] get_z_ze_res(SimpleMatrix X, ArrayList<Satellite> satList, String[] obsvCodeList) {
		int n = satList.size();
		int m = obsvCodeList.length;
		double[][] z = new double[n][1];
		double[][] ze = new double[n][1];
		double[] residual = new double[n];
		double[] estPos = new double[] { X.get(0), X.get(1), X.get(2) };
		double[] rxClkOff = new double[m];// in meters
		for (int i = 0; i < m; i++) {
			rxClkOff[i] = X.get(i + 3);
		}
		for (int i = 0; i < n; i++) {
			Satellite sat = satList.get(i);
			String obsvCode = sat.getObsvCode();
			double approxPR = Math.sqrt(IntStream.range(0, 3).mapToDouble(j -> estPos[j] - sat.getSatEci()[j])
					.map(j -> j * j).reduce(0, (j, k) -> j + k));
			for (int j = 0; j < m; j++) {
				if (obsvCode.equals(obsvCodeList[j])) {
					approxPR += rxClkOff[j];
				}
			}
			z[i][0] = sat.getPseudorange() - approxPR;
			residual[i] = z[i][0] - ze[i][0];
		}

		return new Object[] { z, ze, residual };
	}

	private void CSDR_old(ArrayList<Satellite> prevSatList,ArrayList<Satellite> satList,boolean useIGS) throws Exception
	{
		// Geometry Free CSD
		int n = satList.size();
		int n_prev = prevSatList.size();
		ArrayList<CycleSlipDetect> csdList = new ArrayList<CycleSlipDetect>();
		ArrayList<Satellite> csdSatList = new ArrayList<Satellite>();
		ArrayList<CycleSlipDetect> testing = new ArrayList<CycleSlipDetect>();
		int initialAmbCount = 0;
		double[] userXYZ = LinearLeastSquare.getEstPos(satList, true, useIGS);
		// Iterating to find common satellites b/w two epochs
		for(int i=0;i<n;i++)
		{
			Satellite sat = satList.get(i);
			String satID = sat.getObsvCode()+sat.getSvid();
			for(int j=0;j<n_prev;j++)
			{
				Satellite prev_sat = prevSatList.get(j);
				String prev_satID = prev_sat.getObsvCode()+prev_sat.getSvid();
				if(satID.equals(prev_satID))
				{
					sat.setPhaseLocked(true);
					SimpleMatrix satEci = new SimpleMatrix(3, 1, true, sat.getSatEci());
					SimpleMatrix prev_satEci = new SimpleMatrix(3, 1, true, prev_sat.getSatEci());
					SimpleMatrix unitLOS = new SimpleMatrix(1,3,true,SatUtil.getUnitLOS(sat.getSatEci(), userXYZ));
					double satVelCorr = unitLOS.mult(satEci.minus(prev_satEci)).get(0);
					double tropoRate = sat.getTropoErr()-prev_sat.getTropoErr();
					double dopplerDR = ((sat.getRangeRate()+prev_sat.getRangeRate())/2)-tropoRate;
					double carrierPhaseDR = sat.getPhase()-prev_sat.getPhase();
					double ionoRate = sat.getIonoErr()-prev_sat.getIonoErr();
					double wavelength = SpeedofLight/sat.getCarrierFrequencyHz();
					
					boolean isCS = false;
					double approxCS = Math.abs(carrierPhaseDR-dopplerDR);
					if(approxCS>wavelength)
					{
						if(approxCS<100*wavelength)
						{
							isCS = true;
							initialAmbCount++;
							testing.add(new CycleSlipDetect(sat, dopplerDR, carrierPhaseDR, ionoRate, isCS,wavelength,satVelCorr,i,approxCS/wavelength,unitLOS));
						}
						else
						{
							sat.setPhaseLocked(false);
							break;
						}

					}
					CycleSlipDetect csdObj = new CycleSlipDetect(sat, dopplerDR, carrierPhaseDR, ionoRate, isCS,wavelength,satVelCorr,i,approxCS,unitLOS);
					csdList.add(csdObj);
					csdSatList.add(sat);
					break;
				}
			}
		}
		
		int l = csdList.size();
		String[] obsvCodeList= SatUtil.findObsvCodeArray(csdSatList);
		int m = obsvCodeList.length;
		SimpleMatrix z = new SimpleMatrix(3*l,1);
		//SimpleMatrix x = new SimpleMatrix(3+m+initialAmbCount+l,1);
		SimpleMatrix H = new SimpleMatrix(3*l,3+m+initialAmbCount+l);
		SimpleMatrix unitLOS = new SimpleMatrix(SatUtil.getUnitLOS(csdSatList, userXYZ));
		SimpleMatrix identity = SimpleMatrix.identity(l);
		
		
		
		H.insertIntoThis(0, 0, unitLOS.scale(-1));
		H.insertIntoThis(l, 0, unitLOS.scale(-1));
		H.insertIntoThis(0, 3+m+initialAmbCount, identity);
		H.insertIntoThis(l, 3+m+initialAmbCount, identity);
		H.insertIntoThis(2*l, 3+m+initialAmbCount, identity);
		SimpleMatrix Cyy_doppler = Weight.getNormCyy(csdSatList, GnssDataConfig.doppler_priorVarOfUnitW);
		SimpleMatrix Cyy_phase = new SimpleMatrix(Cyy_doppler);
		SimpleMatrix Cyy_ionoRate = new SimpleMatrix(l,l);
		
		
		int ctr = 3+m;
		for(int i =0;i<l;i++)
		{
			CycleSlipDetect csdObj = csdList.get(i);
			
			z.set(i, csdObj.getCarrierPhaseDR()-csdObj.getSatVelCorr());
			z.set(i+l, csdObj.getDopplerDR()-csdObj.getSatVelCorr());
			z.set(i+l+l, -csdObj.getIonoRate());
			
			String obsvCode = csdSatList.get(i).getObsvCode();
			double wavelength = csdObj.getWavelength();
			for(int j=0;j<m;j++)
			{
				if(obsvCodeList[j].equals(obsvCode))
				{
					H.set(i, 3+j, 1);
					H.set(l+i, 3+j, 1);
					if(csdObj.isCS())
					{
						H.set(i, ctr, wavelength);
						//Cyy_phase.set(i, i, 1e8);
						ctr++;
					}
					Cyy_ionoRate.set(i,i,ionoErrCoeff/ Math.pow(SpeedofLight/wavelength, 2));
				}
			}
		}
		SimpleMatrix W = new SimpleMatrix(3*l,3*l);
		W.insertIntoThis(0, 0,Cyy_phase.invert());
		W.insertIntoThis(l, l,Cyy_doppler.invert());
		W.insertIntoThis(2*l, 2*l,Cyy_ionoRate.invert());
		System.out.println(H);
		
		SimpleMatrix Ht = H.transpose();
		SimpleMatrix HtWHinv = (Ht.mult(W).mult(H)).invert();
		SimpleMatrix x = HtWHinv.mult(Ht).mult(W).mult(z);
		SimpleMatrix floatAmb = x.extractMatrix(3+m, 3+m+initialAmbCount,0, 1 );
		SimpleMatrix floatAmbCov = HtWHinv.extractMatrix(3+m, 3+m+initialAmbCount, 3+m, 3+m+initialAmbCount);
		
		IntegerLeastSquares ILS = new IntegerLeastSquares(floatAmb, floatAmbCov);
		boolean isFixed = ILS.process();
		System.out.println();
		if(isFixed)
		{
//			SimpleMatrix intAmb = ILS.getIntAmb();
//			for(int i=0;i<n;i++)
//			{
//				Satellite sat = csdSatList.get(i);
//				sat.setPhase(sat.getPhase()+intAmb.get(i));
//				
//			}
		}
		else
		{
//			for(int i=0;i<n;i++)
//			{
//				Satellite sat = csdSatList.get(i);
//				
//				
//			}
		}
		
	}
	
	private void CSDR_WLS(ArrayList<Satellite> prevSatList,ArrayList<Satellite> satList,boolean useIGS) throws Exception
	{
		// Geometry Free CSD
		int n = satList.size();
		int n_prev = prevSatList.size();
		ArrayList<CycleSlipDetect> csdList = new ArrayList<CycleSlipDetect>();
		ArrayList<CycleSlipDetect> testing = new ArrayList<CycleSlipDetect>();
		ArrayList<Satellite> csdSatList = new ArrayList<Satellite>();
		int initialAmbCount = 0;
		double[] userXYZ = LinearLeastSquare.getEstPos(satList, true, useIGS);
		// Iterating to find common satellites b/w two epochs
		for(int i=0;i<n;i++)
		{
			Satellite sat = satList.get(i);
			String satID = sat.getObsvCode()+sat.getSvid();
			for(int j=0;j<n_prev;j++)
			{
				Satellite prev_sat = prevSatList.get(j);
				String prev_satID = prev_sat.getObsvCode()+prev_sat.getSvid();
				if(satID.equals(prev_satID))
				{
					sat.setPhaseLocked(true);
					SimpleMatrix satEci = new SimpleMatrix(3, 1, true, sat.getSatEci());
					SimpleMatrix prev_satEci = new SimpleMatrix(3, 1, true, prev_sat.getSatEci());
					SimpleMatrix unitLOS = new SimpleMatrix(1,3,true,SatUtil.getUnitLOS(sat.getSatEci(), userXYZ));
					double ionoRate = sat.getIonoErr()-prev_sat.getIonoErr();
					double satVelCorr = unitLOS.mult(satEci.minus(prev_satEci)).get(0);
					double tropoRate = sat.getTropoErr()-prev_sat.getTropoErr();
					double dopplerDR = ((sat.getRangeRate()+prev_sat.getRangeRate())/2)-tropoRate+ionoRate;
					double carrierPhaseDR = sat.getPhase()-prev_sat.getPhase()+ionoRate;
					
					double wavelength = SpeedofLight/sat.getCarrierFrequencyHz();
					
					boolean isCS = false;
					double approxCS = Math.abs(carrierPhaseDR-dopplerDR);
					if(approxCS>wavelength)
					{
						if(approxCS<25*wavelength)
						{
							isCS = true;
							initialAmbCount++;
							testing.add(new CycleSlipDetect(sat, dopplerDR, carrierPhaseDR, ionoRate, isCS,wavelength,satVelCorr,i,approxCS/wavelength,unitLOS));
						}
						else
						{
							sat.setPhaseLocked(false);
							break;
						}

					}
					CycleSlipDetect csdObj = new CycleSlipDetect(sat, dopplerDR, carrierPhaseDR, ionoRate, isCS,wavelength,satVelCorr,i,approxCS,unitLOS);
					csdList.add(csdObj);
					csdSatList.add(sat);
					break;
				}
			}
		}
		
		int l = csdList.size();
		String[] obsvCodeList= SatUtil.findObsvCodeArray(csdSatList);
		int m = obsvCodeList.length;
		SimpleMatrix z = new SimpleMatrix(2*l,1);
		//SimpleMatrix x = new SimpleMatrix(3+m+initialAmbCount+l,1);
		SimpleMatrix H = new SimpleMatrix(2*l,3+m+initialAmbCount);
		SimpleMatrix unitLOS = new SimpleMatrix(SatUtil.getUnitLOS(csdSatList, userXYZ));
		
		
		
		
		H.insertIntoThis(0, 0, unitLOS.scale(-1));
		H.insertIntoThis(l, 0, unitLOS.scale(-1));
		
		SimpleMatrix Czz_doppler = Weight.getNormCyy(csdSatList, GnssDataConfig.doppler_priorVarOfUnitW);
		SimpleMatrix Czz_phase = new SimpleMatrix(Czz_doppler);
		
		
		
		int ctr = 3+m;
		for(int i =0;i<l;i++)
		{
			CycleSlipDetect csdObj = csdList.get(i);
			
			z.set(i, csdObj.getCarrierPhaseDR()-csdObj.getSatVelCorr());
			z.set(i+l, csdObj.getDopplerDR()-csdObj.getSatVelCorr());
			
			
			String obsvCode = csdSatList.get(i).getObsvCode();
			double wavelength = csdObj.getWavelength();
			for(int j=0;j<m;j++)
			{
				if(obsvCodeList[j].equals(obsvCode))
				{
					H.set(i, 3+j, 1);
					H.set(l+i, 3+j, 1);
					if(csdObj.isCS())
					{
						H.set(i, ctr, wavelength);
						//Cyy_phase.set(i, i, 1e8);
						ctr++;
					}
					
				}
			}
		}
		SimpleMatrix W = new SimpleMatrix(2*l,2*l);
		W.insertIntoThis(0, 0,Czz_phase.invert());
		W.insertIntoThis(l, l,Czz_doppler.invert());
		
		System.out.println(H);
		
		SimpleMatrix Ht = H.transpose();
		SimpleMatrix HtWHinv = (Ht.mult(W).mult(H)).invert();
		SimpleMatrix x = HtWHinv.mult(Ht).mult(W).mult(z);
		SimpleMatrix floatAmb = x.extractMatrix(3+m, 3+m+initialAmbCount,0, 1 );
		SimpleMatrix floatAmbCov = HtWHinv.extractMatrix(3+m, 3+m+initialAmbCount, 3+m, 3+m+initialAmbCount);
		
		SimpleMatrix adaptW = computeAdaptiveW(H, x, z, W, HtWHinv);
		SimpleMatrix adaptCxx = (Ht.mult(adaptW).mult(H)).invert();
		SimpleMatrix adaptFloatAmbCov = adaptCxx.extractMatrix(3+m, 3+m+initialAmbCount, 3+m, 3+m+initialAmbCount);
		
		IntegerLeastSquares ILS = new IntegerLeastSquares(floatAmb, adaptFloatAmbCov);
		boolean isFixed = ILS.process();
		System.out.println();
		if(isFixed)
		{
//			SimpleMatrix intAmb = ILS.getIntAmb();
//			for(int i=0;i<n;i++)
//			{
//				Satellite sat = csdSatList.get(i);
//				sat.setPhase(sat.getPhase()+intAmb.get(i));
//				
//			}
		}
		else
		{
//			for(int i=0;i<n;i++)
//			{
//				Satellite sat = csdSatList.get(i);
//				
//				
//			}
		}
		
	}
	
	private void CSDR_LS(ArrayList<Satellite> prevSatList,ArrayList<Satellite> satList,boolean useIGS) throws Exception
	{
		// Geometry Free CSD
		int n = satList.size();
		int n_prev = prevSatList.size();
		ArrayList<CycleSlipDetect> csdList = new ArrayList<CycleSlipDetect>();
		ArrayList<CycleSlipDetect> testing = new ArrayList<CycleSlipDetect>();
		ArrayList<Satellite> csdSatList = new ArrayList<Satellite>();
		int initialAmbCount = 0;
		double[] userXYZ = LinearLeastSquare.getEstPos(satList, true, useIGS);
		// Iterating to find common satellites b/w two epochs
		for(int i=0;i<n;i++)
		{
			Satellite sat = satList.get(i);
			String satID = sat.getObsvCode()+sat.getSvid();
			for(int j=0;j<n_prev;j++)
			{
				Satellite prev_sat = prevSatList.get(j);
				String prev_satID = prev_sat.getObsvCode()+prev_sat.getSvid();
				if(satID.equals(prev_satID))
				{
					sat.setPhaseLocked(true);
					SimpleMatrix satEci = new SimpleMatrix(3, 1, true, sat.getSatEci());
					SimpleMatrix prev_satEci = new SimpleMatrix(3, 1, true, prev_sat.getSatEci());
					SimpleMatrix unitLOS = new SimpleMatrix(1,3,true,SatUtil.getUnitLOS(sat.getSatEci(), userXYZ));
					double ionoRate = sat.getIonoErr()-prev_sat.getIonoErr();
					double satVelCorr = unitLOS.mult(satEci.minus(prev_satEci)).get(0);
					double tropoRate = sat.getTropoErr()-prev_sat.getTropoErr();
					double dopplerDR = ((sat.getRangeRate()+prev_sat.getRangeRate())/2)-tropoRate+ionoRate;
					double carrierPhaseDR = sat.getPhase()-prev_sat.getPhase()+ionoRate;
					
					double wavelength = SpeedofLight/sat.getCarrierFrequencyHz();
					
					boolean isCS = false;
					double approxCS = Math.abs(carrierPhaseDR-dopplerDR);
					if(approxCS>wavelength)
					{
						if(approxCS<100*wavelength)
						{
							isCS = true;
							initialAmbCount++;
							testing.add(new CycleSlipDetect(sat, dopplerDR, carrierPhaseDR, ionoRate, isCS,wavelength,satVelCorr,i,approxCS/wavelength,unitLOS));
						}
						else
						{
							sat.setPhaseLocked(false);
							break;
						}

					}
					CycleSlipDetect csdObj = new CycleSlipDetect(sat, dopplerDR, carrierPhaseDR, ionoRate, isCS,wavelength,satVelCorr,i,approxCS,unitLOS);
					csdList.add(csdObj);
					csdSatList.add(sat);
					break;
				}
			}
		}
		
		int l = csdList.size();
		String[] obsvCodeList= SatUtil.findObsvCodeArray(csdSatList);
		int m = obsvCodeList.length;
		SimpleMatrix z = new SimpleMatrix(2*l,1);
		//SimpleMatrix x = new SimpleMatrix(3+m+initialAmbCount+l,1);
		SimpleMatrix H = new SimpleMatrix(2*l,3+m+initialAmbCount);
		SimpleMatrix unitLOS = new SimpleMatrix(SatUtil.getUnitLOS(csdSatList, userXYZ));
		
		
		
		
		H.insertIntoThis(0, 0, unitLOS.scale(-1));
		H.insertIntoThis(l, 0, unitLOS.scale(-1));
		
		SimpleMatrix Czz_doppler = Weight.getNormCyy(csdSatList, GnssDataConfig.doppler_priorVarOfUnitW);
		SimpleMatrix Czz_phase = new SimpleMatrix(Czz_doppler);
		
		
		
		int ctr = 3+m;
		for(int i =0;i<l;i++)
		{
			CycleSlipDetect csdObj = csdList.get(i);
			
			z.set(i, csdObj.getCarrierPhaseDR()-csdObj.getSatVelCorr());
			z.set(i+l, csdObj.getDopplerDR()-csdObj.getSatVelCorr());
			
			
			String obsvCode = csdSatList.get(i).getObsvCode();
			double wavelength = csdObj.getWavelength();
			for(int j=0;j<m;j++)
			{
				if(obsvCodeList[j].equals(obsvCode))
				{
					H.set(i, 3+j, 1);
					H.set(l+i, 3+j, 1);
					if(csdObj.isCS())
					{
						H.set(i, ctr, wavelength);
						//Cyy_phase.set(i, i, 1e8);
						ctr++;
					}
					
				}
			}
		}
//		SimpleMatrix W = new SimpleMatrix(2*l,2*l);
//		W.insertIntoThis(0, 0,Cyy_phase.invert());
//		W.insertIntoThis(l, l,Cyy_doppler.invert());
		SimpleMatrix W = SimpleMatrix.identity(2*l);
		System.out.println(H);
		
		SimpleMatrix Ht = H.transpose();
		SimpleMatrix HtWHinv = (Ht.mult(W).mult(H)).invert();
		SimpleMatrix x = HtWHinv.mult(Ht).mult(W).mult(z);
		SimpleMatrix floatAmb = x.extractMatrix(3+m, 3+m+initialAmbCount,0, 1 );
		SimpleMatrix floatAmbCov = HtWHinv.extractMatrix(3+m, 3+m+initialAmbCount, 3+m, 3+m+initialAmbCount);
		
		SimpleMatrix adaptW = computeAdaptiveW(H, x, z, W, HtWHinv);
		SimpleMatrix adaptCxx = (Ht.mult(adaptW).mult(H)).invert();
		SimpleMatrix adaptFloatAmbCov = adaptCxx.extractMatrix(3+m, 3+m+initialAmbCount, 3+m, 3+m+initialAmbCount);
		
		
		IntegerLeastSquares ILS = new IntegerLeastSquares(floatAmb, adaptFloatAmbCov);
		boolean isFixed = ILS.process();
		System.out.println();
		if(isFixed)
		{
//			SimpleMatrix intAmb = ILS.getIntAmb();
//			for(int i=0;i<n;i++)
//			{
//				Satellite sat = csdSatList.get(i);
//				sat.setPhase(sat.getPhase()+intAmb.get(i));
//				
//			}
		}
		else
		{
//			for(int i=0;i<n;i++)
//			{
//				Satellite sat = csdSatList.get(i);
//				
//				
//			}
		}
		
	}
	
	private void CSDR_Orekit(ArrayList<Satellite> prevSatList,ArrayList<Satellite> satList,boolean useIGS) throws Exception
	{
		// Geometry Free CSD
		int n = satList.size();
		int n_prev = prevSatList.size();
		ArrayList<CycleSlipDetect> csdList = new ArrayList<CycleSlipDetect>();
		ArrayList<CycleSlipDetect> testing = new ArrayList<CycleSlipDetect>();
		ArrayList<Satellite> csdSatList = new ArrayList<Satellite>();
		int initialAmbCount = 0;
		double[] userXYZ = LinearLeastSquare.getEstPos(satList, true, useIGS);
		// Iterating to find common satellites b/w two epochs
		for(int i=0;i<n;i++)
		{
			Satellite sat = satList.get(i);
			String satID = sat.getObsvCode()+sat.getSvid();
			for(int j=0;j<n_prev;j++)
			{
				Satellite prev_sat = prevSatList.get(j);
				String prev_satID = prev_sat.getObsvCode()+prev_sat.getSvid();
				if(satID.equals(prev_satID))
				{
					sat.setPhaseLocked(true);
					SimpleMatrix satEci = new SimpleMatrix(3, 1, true, sat.getSatEci());
					SimpleMatrix prev_satEci = new SimpleMatrix(3, 1, true, prev_sat.getSatEci());
					SimpleMatrix unitLOS = new SimpleMatrix(1,3,true,SatUtil.getUnitLOS(sat.getSatEci(), userXYZ));
					double ionoRate = sat.getIonoErr()-prev_sat.getIonoErr();
					double satVelCorr = unitLOS.mult(satEci.minus(prev_satEci)).get(0);
					double tropoRate = sat.getTropoErr()-prev_sat.getTropoErr();
					double dopplerDR = ((sat.getRangeRate()+prev_sat.getRangeRate())/2)-tropoRate+ionoRate;
					double carrierPhaseDR = sat.getPhase()-prev_sat.getPhase()+ionoRate;
					
					double wavelength = SpeedofLight/sat.getCarrierFrequencyHz();
					
					boolean isCS = false;
					double approxCS = Math.abs(carrierPhaseDR-dopplerDR);
					if(approxCS>3*wavelength)
					{
						if(approxCS<100*wavelength)
						{
							isCS = true;
							initialAmbCount++;
							testing.add(new CycleSlipDetect(sat, dopplerDR, carrierPhaseDR, ionoRate, isCS,wavelength,satVelCorr,i,approxCS/wavelength,unitLOS));
						}
						else
						{
							sat.setPhaseLocked(false);
							break;
						}

					}
					CycleSlipDetect csdObj = new CycleSlipDetect(sat, dopplerDR, carrierPhaseDR, ionoRate, isCS,wavelength,satVelCorr,i,approxCS,unitLOS);
					csdList.add(csdObj);
					csdSatList.add(sat);
					break;
				}
			}
		}
		
		int l = csdList.size();
		String[] obsvCodeList= SatUtil.findObsvCodeArray(csdSatList);
		int m = obsvCodeList.length;
		SimpleMatrix z = new SimpleMatrix(2*l,1);
		//SimpleMatrix x = new SimpleMatrix(3+m+initialAmbCount+l,1);
		SimpleMatrix H = new SimpleMatrix(2*l,3+m+initialAmbCount);
		SimpleMatrix unitLOS = new SimpleMatrix(SatUtil.getUnitLOS(csdSatList, userXYZ));
		
		
		
		
		H.insertIntoThis(0, 0, unitLOS.scale(-1));
		H.insertIntoThis(l, 0, unitLOS.scale(-1));
		
		SimpleMatrix Czz_doppler = Weight.getNormCyy(csdSatList, GnssDataConfig.doppler_priorVarOfUnitW);
		SimpleMatrix Czz_phase = new SimpleMatrix(Czz_doppler).scale(0.1);
		
		
		
		int ctr = 3+m;
		for(int i =0;i<l;i++)
		{
			CycleSlipDetect csdObj = csdList.get(i);
			
			z.set(i, csdObj.getCarrierPhaseDR()-csdObj.getSatVelCorr());
			z.set(i+l, csdObj.getDopplerDR()-csdObj.getSatVelCorr());
			
			
			String obsvCode = csdSatList.get(i).getObsvCode();
			double wavelength = csdObj.getWavelength();
			for(int j=0;j<m;j++)
			{
				if(obsvCodeList[j].equals(obsvCode))
				{
					H.set(i, 3+j, 1);
					H.set(l+i, 3+j, 1);
					if(csdObj.isCS())
					{
						H.set(i, ctr, wavelength);
						//Cyy_phase.set(i, i, 1e8);
						ctr++;
					}
					
				}
			}
		}
		SimpleMatrix W = new SimpleMatrix(2*l,2*l);
		W.insertIntoThis(0, 0,Czz_phase.invert());
		W.insertIntoThis(l, l,Czz_doppler.invert());
		
		System.out.println(H);
		
		SimpleMatrix Ht = H.transpose();
		SimpleMatrix HtWHinv = (Ht.mult(W).mult(H)).invert();
		SimpleMatrix x = HtWHinv.mult(Ht).mult(W).mult(z);
		SimpleMatrix floatAmb = x.extractMatrix(3+m, 3+m+initialAmbCount,0, 1 );
		SimpleMatrix floatAmbCov = HtWHinv.extractMatrix(3+m, 3+m+initialAmbCount, 3+m, 3+m+initialAmbCount);
		
		
		
		
		
		RealMatrix Cxx_floatAmb_hat =  new Array2DRowRealMatrix(Matrix.matrix2Array(floatAmbCov));
		LambdaMethod lm = new LambdaMethod();
		
		// Full ambiguity Resolution
		
		int[] indirection = IntStream.range(0, initialAmbCount).toArray();
		IntegerLeastSquareSolution[] ILSsol = lm.solveILS(5, Matrix.matrix2ArrayVec(floatAmb), indirection, Cxx_floatAmb_hat);
		SimpleRatioAmbiguityAcceptance ratioTest = new SimpleRatioAmbiguityAcceptance(1.0/3.0);
		IntegerLeastSquareSolution acceptedILS = ratioTest.accept(ILSsol);
		
		
		if(acceptedILS==null)
		{
			
			int count = initialAmbCount-1;
			int[] arr = IntStream.range(0, initialAmbCount).toArray();
			while(acceptedILS==null&&count>=0)
			{
				ArrayList<int[]> combs = Combination.getCombination(arr,initialAmbCount,count);
				for(int i=0;i<combs.size();i++)
				{
					int[] comb = combs.get(i);
					double[] newFloatAmb = IntStream.range(0, count).mapToDouble(j->floatAmb.get(comb[j])).toArray();
					ILSsol = lm.solveILS(5, newFloatAmb,comb ,Cxx_floatAmb_hat);
					acceptedILS = ratioTest.accept(ILSsol);
					if(acceptedILS!=null)
					{
						break;
					}
				}
				count--;
			}
			

		}
			
		
		System.out.println();

		
	}
	
	private SimpleMatrix computeAdaptiveW(SimpleMatrix H, SimpleMatrix x, SimpleMatrix z, SimpleMatrix W,SimpleMatrix HtWHinv)
	{
		SimpleMatrix residual = z.minus(H.mult(x));
		int n = z.getNumElements();
		SimpleMatrix  identity = SimpleMatrix.identity(n);
		SimpleMatrix redunMat = identity.minus(H.mult(HtWHinv).mult(H.transpose()).mult(W)); 
		SimpleMatrix adaptW = new SimpleMatrix(n,n);
		int n_2 = n/2;
		for(int i=0;i<n_2;i++)
		{
			adaptW.set(i,i,1/(Math.pow(residual.get(i+n_2), 2)/redunMat.get(i+n_2,i+n_2)));
			adaptW.set(i+n_2,i+n_2,1/(Math.pow(residual.get(i+n_2), 2)/redunMat.get(i+n_2,i+n_2)));
		}
		return adaptW;
	}
	
	

	
}
