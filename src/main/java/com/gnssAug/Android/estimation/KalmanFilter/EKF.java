package com.gnssAug.Android.estimation.KalmanFilter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.stream.IntStream;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.ejml.dense.row.MatrixFeatures_DDRM;
import org.ejml.simple.SimpleMatrix;

import com.gnssAug.Android.estimation.LinearLeastSquare;
import com.gnssAug.Android.estimation.KalmanFilter.Models.Flag;
import com.gnssAug.Android.estimation.KalmanFilter.Models.KFconfig;
import com.gnssAug.Android.models.Satellite;
import com.gnssAug.utility.LatLonUtil;
import com.gnssAug.utility.Matrix;
import com.gnssAug.utility.SatUtil;
import com.gnssAug.utility.Weight;
import com.gnssAug.Android.constants.GnssDataConfig;
import com.gnssAug.Android.constants.Measurement;
import com.gnssAug.Android.constants.State;

public class EKF {
	private final double SpeedofLight = 299792458;
	private KFconfig kfObj;
	private double[] innovation;

	private TreeMap<Long, HashMap<Measurement, double[]>> innovationMap;
	private TreeMap<Long, HashMap<Measurement, double[]>> residualMap;
	private TreeMap<Long, HashMap<State, double[]>> measNoiseMap;
	private TreeMap<Long, HashMap<Measurement, Double>> postVarOfUnitWMap;
	// Posteriori Err Cov
	private TreeMap<Long, HashMap<State, SimpleMatrix>> errCovMap;
	private HashMap<Measurement, ArrayList<double[]>> redundancyMap;

	// Satellite Count
	private TreeMap<Long, HashMap<Measurement, Long>> satCountMap;
	private TreeMap<Long, HashMap<Measurement, ArrayList<Satellite>>> satListMap;

	public EKF() {
		kfObj = new KFconfig();
	}

	public TreeMap<Long, double[]> process(TreeMap<Long, ArrayList<Satellite>> SatMap, ArrayList<Long> timeList,
			Flag flag, boolean useDoppler, boolean useIGS, String[] obsvCodeList, boolean doAnalyze, boolean doTest,
			boolean outlierAnalyze) throws Exception {

		return process(SatMap, timeList, flag, useDoppler, useIGS, obsvCodeList, doAnalyze, doTest, outlierAnalyze,
				false, false);

	}

	public TreeMap<Long, double[]> process(TreeMap<Long, ArrayList<Satellite>> SatMap, ArrayList<Long> timeList,
			Flag flag, boolean useDoppler, boolean useIGS, String[] obsvCodeList, boolean doAnalyze, boolean doTest,
			boolean outlierAnalyze, boolean complementary) throws Exception {

		return process(SatMap, timeList, flag, useDoppler, useIGS, obsvCodeList, doAnalyze, doTest, outlierAnalyze,
				complementary, false);

	}

	public TreeMap<Long, double[]> process(TreeMap<Long, ArrayList<Satellite>> SatMap, ArrayList<Long> timeList,
			Flag flag, boolean useDoppler, boolean useIGS, String[] obsvCodeList, boolean doAnalyze, boolean doTest,
			boolean outlierAnalyze, boolean complementary, boolean useEstVel) throws Exception {

		boolean isWeighted = true;
		int n = 0;
		int m = obsvCodeList.length;
		/* constant position model - state vector(n=5) -> (x,y,z,cdt,cdt_dot) */
		/*
		 * constant velocity model - state vector(n=8) ->
		 * (x,y,z,cdt,x_dot,y_dot,z_dot,cdt_dot)
		 */
		if (flag == Flag.POSITION) {
			n = 3 + (2 * m);
		} else if (flag == Flag.VELOCITY) {
			n = 6 + (2 * m);
		}
		double[][] x = new double[n][1];
		double[][] P = new double[n][n];
		/*
		 * state XYZ is intialized WLS generated Rx position estimated using first epoch
		 * data, a-priori estimate error covariance matrix for state XYZ is therefore
		 * assigned 25 m^2 value. Other state variables are assigned infinite(big)
		 * variance
		 */
		double[] intialECEF = LinearLeastSquare.getEstPos(SatMap.firstEntry().getValue(), true, useIGS);
		
		IntStream.range(0, 3 + m).forEach(i -> x[i][0] = intialECEF[i]);
		IntStream.range(0, 3 + m).forEach(i -> P[i][i] = 100);
		if (flag == Flag.POSITION) {
			IntStream.range(3 + m, 3 + (2 * m)).forEach(i -> P[i][i] = 1e13);
		} else {
			double[] intialVel = LinearLeastSquare.getEstVel(SatUtil.createCopy(SatMap.firstEntry().getValue()), isWeighted,
					true, doTest, false, intialECEF, useIGS);
			SimpleMatrix intialVelCov =  LinearLeastSquare.getCxx_hat(Measurement.Doppler, "ECEF");
			IntStream.range(3 + m, 6 + (2 * m)).forEach(i -> x[i][0] = intialVel[i - (3 + m)]);
			for(int i=3+m;i<6+(2*m);i++)
			{
				for(int j=3+m;j<6+(2*m);j++)
				{
					P[i][j] = intialVelCov.get(i-(3+m), j-(3+m));
				}
			}
			//IntStream.range(3 + m, 6 + m).forEach(i -> P[i][i] = 4);
			IntStream.range(6 + m, 6 + (2 * m)).forEach(i -> P[i][i] += 1e5);
		}

		kfObj.setState_ProcessCov(x, P);
		if (doAnalyze) {
			innovationMap = new TreeMap<Long, HashMap<Measurement, double[]>>();
			errCovMap = new TreeMap<Long, HashMap<State, SimpleMatrix>>();
			residualMap = new TreeMap<Long, HashMap<Measurement, double[]>>();
			postVarOfUnitWMap = new TreeMap<Long, HashMap<Measurement, Double>>();
			redundancyMap = new HashMap<Measurement, ArrayList<double[]>>();
			satCountMap = new TreeMap<Long, HashMap<Measurement, Long>>();
			satListMap = new TreeMap<Long, HashMap<Measurement, ArrayList<Satellite>>>();
			measNoiseMap = new TreeMap<Long, HashMap<State, double[]>>();
		}
		// Begin iteration or recursion
		return iterate(SatMap, timeList, flag, useDoppler, obsvCodeList, doAnalyze, doTest, outlierAnalyze, useIGS,
				complementary, useEstVel,isWeighted);

	}

	private TreeMap<Long, double[]> iterate(TreeMap<Long, ArrayList<Satellite>> SatMap, ArrayList<Long> timeList,
			Flag flag, boolean useDoppler, String[] obsvCodeList, boolean doAnalyze, boolean doTest,
			boolean outlierAnalyze, boolean useIGS, boolean complementary, boolean useEstVel,boolean isWeighted) throws Exception {
		TreeMap<Long, double[]> estStateMap = new TreeMap<Long, double[]>();
		Measurement[] measArr = useDoppler ? new Measurement[] { Measurement.Pseudorange, Measurement.Doppler }
				: new Measurement[] { Measurement.Pseudorange };
		int m = obsvCodeList.length;
		long time = timeList.get(0);
		// Start from 2nd epoch
		for (int i = 1; i < timeList.size(); i++) {
			long currentTime = timeList.get(i);
			ArrayList<Satellite> satList = SatMap.get(currentTime);
			double deltaT = (currentTime - time) / 1e3;
//			if (i<150) {
//				System.out.println("P(+)  \n"+kfObj.getCovariance().extractMatrix(0, 3+m, 0, 3+m));
//			}
			// Perform Predict and Update
			runFilter(deltaT, satList, flag, useDoppler, obsvCodeList, useIGS, currentTime, doAnalyze, doTest,
					outlierAnalyze, complementary, measArr, useEstVel, i,isWeighted);
			// Fetch Posteriori state estimate and estimate error covariance matrix
			SimpleMatrix x = kfObj.getState();
			SimpleMatrix P = kfObj.getCovariance();
			int n = x.numRows();
			double[] estState = new double[n];
			for (int j = 0; j < n; j++) {
				estState[j] = x.get(j);
			}
			// Add position estimate to the list
			estStateMap.put(currentTime, estState);
			if (doAnalyze) {

				// Convert to ENU frame
				SimpleMatrix R = new SimpleMatrix(n, n);
				SimpleMatrix rotMat = new SimpleMatrix(
						LatLonUtil.getEcef2EnuRotMat(new double[] { estState[0], estState[1], estState[2] }));
				R.insertIntoThis(0, 0, rotMat);
				IntStream.range(3, 3 + m).forEach(j -> R.set(j, j, 1));
				if (flag == Flag.VELOCITY) {
					R.insertIntoThis(3 + m, 3 + m, rotMat);
					IntStream.range(6 + m, 6 + (2 * m)).forEach(j -> R.set(j, j, 1));
				} else {
					IntStream.range(3 + m, 3 + (2 * m)).forEach(j -> R.set(j, j, 1));
				}
				SimpleMatrix errCov = R.mult(P).mult(R.transpose());
				errCovMap.computeIfAbsent(currentTime, j -> new HashMap<State, SimpleMatrix>()).put(State.Position,
						errCov.extractMatrix(0, 3 + m, 0, 3 + m));
				int satCount = satList.size();
				innovationMap.computeIfAbsent(currentTime, j -> new HashMap<Measurement, double[]>())
						.put(Measurement.Pseudorange, Arrays.copyOfRange(innovation, 0, satCount));
				if (flag == Flag.VELOCITY) {
					errCovMap.get(currentTime).put(State.Velocity,
							errCov.extractMatrix(3 + m, 6 + (2 * m), 3 + m, 6 + (2 * m)));
				}
				if (useDoppler) {
					innovationMap.get(currentTime).put(Measurement.Doppler,
							Arrays.copyOfRange(innovation, satCount, 2 * satCount));
				}

			}
			/*
			 * Check whether estimate error covariance matrix is positive semidefinite
			 * before further proceeding
			 */
			if (!MatrixFeatures_DDRM.isPositiveSemidefinite(P.getMatrix())) {

				throw new Exception("PositiveDefinite test Failed");
			}
			time = currentTime;

		}

		return estStateMap;
	}

	private void runFilter(double deltaT, ArrayList<Satellite> satList, Flag flag, boolean useDoppler,
			String[] obsvCodeList, boolean useIGS, long currentTime, boolean doAnalyze, boolean doTest,
			boolean outlierAnalyze, boolean complementary, Measurement[] measArr, boolean useEstVel, int ct, boolean isWeighted)
			throws Exception {

		boolean useAndroidW = false;
		// Satellite count
		int n = satList.size();
		int m = obsvCodeList.length;

		SimpleMatrix priorP = new SimpleMatrix(kfObj.getCovariance());
		SimpleMatrix priorX = new SimpleMatrix(kfObj.getState());

		// Assign Q and F matrix
		kfObj.config(deltaT, flag, m, useDoppler, complementary, useEstVel);
		kfObj.predict();
		if (useEstVel && complementary) {
			SimpleMatrix P = kfObj.getCovariance();
			for (int i = 0; i < 3 + m; i++) {
				for (int j = 0; j < 3 + m; j++) {
					P.set(i, j + 3 + m, 0);
					P.set(i + 3 + m, j, 0);
					if (i != j) {
						P.set(i + 3 + m, j + 3 + m, 0);
					}
				}
			}
			kfObj.setProcessCov(P);
		}

		SimpleMatrix x = kfObj.getState();
		double[] estPos = new double[] { x.get(0), x.get(1), x.get(2) };
		double[] estVel_Z = null;
		SimpleMatrix Czz_estVel = null;
		int size = useDoppler ? 2 * n : n;
		if (useEstVel) {
			estVel_Z = LinearLeastSquare.getEstVel(SatUtil.createCopy(satList), isWeighted, true, doTest, false, estPos,
					useIGS);
			Czz_estVel = LinearLeastSquare.getCxx_hat(Measurement.Doppler, "ECEF");
			size = n + 3 + m;
		}
		/*
		 * H is the Jacobian matrix of partial derivatives Observation StateModel(h) of
		 * with respect to x
		 */
		SimpleMatrix H = getJacobian(satList, estPos, flag, useDoppler, obsvCodeList, useEstVel);
		Object[] z_ze_res = get_z_ze_res(x, satList, obsvCodeList, H, useDoppler, useEstVel, estVel_Z);

		// Measurement vector
		double[][] z = (double[][]) z_ze_res[0];
		// Estimated Measurement vector
		double[][] ze = (double[][]) z_ze_res[1];
		// Innovation vector
		innovation = (double[]) z_ze_res[2];

		// Measurement Noise
		SimpleMatrix R = new SimpleMatrix(size, size);

		if (isWeighted) {
			if (useAndroidW) {
				for (int i = 0; i < n; i++) {
					R.set(i, i, Math.pow(satList.get(i).getReceivedSvTimeUncertaintyNanos() * SpeedofLight * 1e-9, 2));
					if (useDoppler) {
						R.set(i + n, i + n, Math.pow(satList.get(i).getPseudorangeRateUncertaintyMetersPerSecond(), 2));
					}
				}
			} else {
//				LinearLeastSquare.getEstPos(satList, true, true, false, false, useIGS);
//				SimpleMatrix Cyy = LinearLeastSquare.getCyy_updated(Measurement.Pseudorange);
//				if (useDoppler) {
//					LinearLeastSquare.getEstVel(satList, true, true, false, false, estPos, useIGS);
//					Cyy.concatRows(LinearLeastSquare.getCyy_updated(Measurement.Doppler));
//				}
//				R = Cyy;
				SimpleMatrix Cyy = Weight.getNormCyy(satList, GnssDataConfig.pseudorange_priorVarOfUnitW);
				R.insertIntoThis(0, 0, Cyy);
				if(useDoppler)
				{
					Cyy = Weight.getNormCyy(satList, GnssDataConfig.doppler_priorVarOfUnitW);
					R.insertIntoThis(n, n, Cyy);
				}
				if(useEstVel)
				{
					R.insertIntoThis(n, n, Czz_estVel);
					for (int i = 3; i < 3 + m; i++) {
						R.set(i + n, i + n, R.get(i + n, i + n) + GnssDataConfig.clkDriftVar);
					}
				}
				
			}
		} else {
			for (int i = 0; i < n; i++) {
				R.set(i, i, GnssDataConfig.pseudorange_priorVarOfUnitW);
				if (useDoppler) {
					R.set(i + n, i + n, GnssDataConfig.doppler_priorVarOfUnitW);
				}
			}
			if (useEstVel) {
				for (int i = 0; i < 3 + m; i++) {
					for (int j = 0; j < 3 + m; j++) {
						R.set(i + n, j + n, Czz_estVel.get(i, j));
					}
				}
				for (int i = 3; i < 3 + m; i++) {
					R.set(i + n, i + n, R.get(i + n, i + n) + GnssDataConfig.clkDriftVar);
				}
			}
		}
		HashMap<Measurement, ArrayList<Satellite>> satMap = new HashMap<Measurement, ArrayList<Satellite>>();
		HashMap<Measurement, ArrayList<Satellite>> testedSatMap = new HashMap<Measurement, ArrayList<Satellite>>();
		for (Measurement meas : measArr) {
			satMap.put(meas, SatUtil.createCopy(satList));
			testedSatMap.put(meas, new ArrayList<Satellite>(satMap.get(meas)));
		}
		if (doTest && !outlierAnalyze) {
			Object[] params = performTesting(R, H, size, m, testedSatMap, satMap, z, ze, useDoppler, useEstVel,ct);
			R = (SimpleMatrix) params[0];
			H = (SimpleMatrix) params[1];
			z = (double[][]) params[2];
			ze = (double[][]) params[3];
		}
		// Perform Update Step
		kfObj.update(z, R, ze, H);
		if (doAnalyze) {
			performAnalysis(testedSatMap, satMap, R, H, priorP, currentTime, n, obsvCodeList, doTest, outlierAnalyze,
					useDoppler, useEstVel, ct);
		}
		if (doTest && outlierAnalyze) {
			kfObj.setState_ProcessCov(priorX, priorP);
			kfObj.predict();
			if (useEstVel && complementary) {
				SimpleMatrix P = kfObj.getCovariance();
				for (int i = 0; i < 3 + m; i++) {
					for (int j = 0; j < 3 + m; j++) {
						P.set(i, j + 3 + m, 0);
						P.set(i + 3 + m, j, 0);
						if (i != j) {
							P.set(i + 3 + m, j + 3 + m, 0);
						}
					}
				}
				kfObj.setProcessCov(P);

			}
			Object[] params = performTesting(R, H, size, m, testedSatMap, satMap, z, ze, useDoppler, useEstVel,ct);
			R = (SimpleMatrix) params[0];
			H = (SimpleMatrix) params[1];
			z = (double[][]) params[2];
			ze = (double[][]) params[3];
			kfObj.update(z, R, ze, H);
			
		}

	}

	private SimpleMatrix getJacobian(ArrayList<Satellite> satList, double[] estECEF, Flag flag, boolean useDoppler,
			String[] obsvCodeList, boolean useEstVel) {
		int n = satList.size();
		int m = obsvCodeList.length;
		int rows = n;
		int stateN = 3 + (2 * m);
		if (flag == Flag.VELOCITY) {
			stateN = 6 + (2 * m);
			if (useDoppler) {
				rows = 2 * n;
			} else if (useEstVel) {
				rows = n + 3 + m;
			}
		}
		double[][] H = new double[rows][stateN];

		for (int i = 0; i < n; i++) {
			Satellite sat = satList.get(i);
			String obsvCode = sat.getObsvCode();
			// Line of Sight vector
			// Its not really a ECI, therefore don't get confused
			double[] LOS = IntStream.range(0, 3).mapToDouble(j -> sat.getSatEci()[j] - estECEF[j]).toArray();
			// Geometric Range
			double GR = Math.sqrt(Arrays.stream(LOS).map(j -> j * j).reduce(0.0, (j, k) -> j + k));
			// Converting LOS to unit vector
			final int _i = i;
			IntStream.range(0, 3).forEach(j -> H[_i][j] = -LOS[j] / GR);
			for (int j = 0; j < m; j++) {
				if (obsvCode.equals(obsvCodeList[j])) {
					H[i][3 + j] = 1;
					if (useDoppler) {
						H[i + n][6 + m + j] = 1;
						// H[i + n][stateN-1] = 1;
					}
				}
			}
			if (useDoppler) {
				IntStream.range(0, 3).forEach(j -> H[_i + n][j + 3 + m] = -LOS[j] / GR);

			}
		}
		if (useEstVel) {
			for (int i = 0; i < 3 + m; i++) {
				H[i + n][i + 3 + m] = 1;
			}
		}
		return new SimpleMatrix(H);

	}

	private void performAnalysis(HashMap<Measurement, ArrayList<Satellite>> testedSatMap,
			HashMap<Measurement, ArrayList<Satellite>> satMap, SimpleMatrix R, SimpleMatrix H, SimpleMatrix priorP,
			long currentTime, int n, String[] obsvCodeList, boolean doTest, boolean outlierAnalyze, boolean useDoppler,
			boolean useEstVel, int ct) {

		SimpleMatrix x = kfObj.getState();
		int n_pr = testedSatMap.get(Measurement.Pseudorange).size();
		int n_doppler = 0;
		int size = n_pr;
		double[] residual = (double[]) get_z_ze_res(x, testedSatMap, obsvCodeList, H, useDoppler)[2];
		double[] pr_res = Arrays.copyOfRange(residual, 0, n_pr);
		// Post-fit residual
		SimpleMatrix e_post_hat_pr = new SimpleMatrix(n_pr, 1, true, pr_res);
		SimpleMatrix Cyy_inv_pr = R.extractMatrix(0, n_pr, 0, n_pr).invert();

		// Compute Redundancies

		SimpleMatrix K = kfObj.getKalmanGain();
		SimpleMatrix HK = H.mult(K);
		SimpleMatrix phi = kfObj.getPhi();
		SimpleMatrix Cvv = kfObj.getCvv();
//		if (ct<5) {
//			System.out.println("Cvv epoch " + ct);
//			System.out.println(Cvv);
//		}
		SimpleMatrix HtCvvInvH = H.transpose().mult(Cvv.invert()).mult(H);
		SimpleMatrix Q = kfObj.getQ();
		double rX = phi.mult(priorP).mult(phi.transpose()).mult(HtCvvInvH).trace();
		double rW = Q.mult(HtCvvInvH).trace();
		int numRows = HK.numRows();
		SimpleMatrix I_HK = SimpleMatrix.identity(numRows).minus(HK);
		double rZ_pr = I_HK.extractMatrix(0, n_pr, 0, n_pr).trace();
		double postVarOfUnitW_pr = e_post_hat_pr.transpose().mult(Cyy_inv_pr).mult(e_post_hat_pr).get(0) / rZ_pr;
		double rZ_doppler = 0;
		double rZ_estVel = 0;
		if (useDoppler) {
			n_doppler = testedSatMap.get(Measurement.Doppler).size();
			size = size + n_doppler;
			rZ_doppler = I_HK.extractMatrix(n_pr, size, n_pr, size).trace();
			double[] doppler_res = Arrays.copyOfRange(residual, n_pr, size);
			SimpleMatrix e_post_hat_doppler = new SimpleMatrix(n_doppler, 1, true, doppler_res);
			SimpleMatrix Cyy_inv_doppler = R.extractMatrix(n_pr, size, n_pr, size).invert();
			double postVarOfUnitW_doppler = e_post_hat_doppler.transpose().mult(Cyy_inv_doppler)
					.mult(e_post_hat_doppler).get(0) / rZ_doppler;
			redundancyMap.computeIfAbsent(Measurement.Doppler, k -> new ArrayList<double[]>())
					.add(new double[] { n_doppler, 0, rX, rW, rZ_doppler });
			postVarOfUnitWMap.computeIfAbsent(currentTime, k -> new HashMap<Measurement, Double>())
					.put(Measurement.Doppler, postVarOfUnitW_doppler);
			residualMap.computeIfAbsent(currentTime, k -> new HashMap<Measurement, double[]>()).put(Measurement.Doppler,
					doppler_res);

		} else if (useEstVel) {
			size += 3 + obsvCodeList.length;
			rZ_estVel = I_HK.extractMatrix(n_pr, size, n_pr, size).trace();
		}
		double rZ = rZ_pr + rZ_doppler + rZ_estVel;
		double rSum = rX + rW + rZ;
		if (Math.abs(numRows - rSum) > 0.01) {
			System.err.println("FATAL ERROR: Redundancy sum is wrong");
		}
		redundancyMap.computeIfAbsent(Measurement.Pseudorange, k -> new ArrayList<double[]>())
				.add(new double[] { n_pr, 0, rX, rW, rZ_pr });

		postVarOfUnitWMap.computeIfAbsent(currentTime, k -> new HashMap<Measurement, Double>())
				.put(Measurement.Pseudorange, postVarOfUnitW_pr);
		residualMap.computeIfAbsent(currentTime, k -> new HashMap<Measurement, double[]>()).put(Measurement.Pseudorange,
				pr_res);

		if (doTest == true) {
			n_pr = n - n_pr;
			n_doppler = n - n_doppler;
		}
		HashMap<Measurement, Long> count = new HashMap<Measurement, Long>();
		count.put(Measurement.Pseudorange, (long) n_pr);
		if (useDoppler) {
			count.put(Measurement.Doppler, (long) n_doppler);
		}
		satCountMap.put(currentTime, count);
		if (outlierAnalyze) {
			satListMap.put(currentTime, satMap);
		} else {
			satListMap.put(currentTime, testedSatMap);
		}

	}

	private Object[] performTesting(SimpleMatrix R, SimpleMatrix H, int size, int m,
			HashMap<Measurement, ArrayList<Satellite>> testedSatMap, HashMap<Measurement, ArrayList<Satellite>> satMap,
			double[][] z, double[][] ze, boolean useDoppler, boolean useEstVel,int ct) throws Exception {
		int n_pr = testedSatMap.get(Measurement.Pseudorange).size();
		int n_doppler = 0;
		int n_estVel = 0;
		double rZ_estVel = 0;
		int cols = H.numCols();
		if (useDoppler) {
			n_doppler = testedSatMap.get(Measurement.Doppler).size();
		}
		int _size = n_pr + n_doppler;
		SimpleMatrix P = kfObj.getCovariance();
		SimpleMatrix Cvv = ((H.mult(P).mult(H.transpose())).plus(R));
		if (useEstVel) {
			Cvv = ((H.extractMatrix(0, _size, 0, 3+m).mult(P.extractMatrix(0, 3+m, 0, 3+m)).mult(H.extractMatrix(0, _size, 0, 3+m).transpose())).plus(R.extractMatrix(0, _size, 0, _size)));
			//Cvv = Cvv.extractMatrix(0, _size, 0, _size);
			n_estVel = 3 + m;
//			SimpleMatrix K = P.mult(H.transpose()).mult(Cvv.invert());
//			SimpleMatrix HK = H.mult(K);
//			SimpleMatrix I_HK = SimpleMatrix.identity(HK.numRows()).minus(HK);
//			
//			rZ_estVel = I_HK.extractMatrix(n_pr, n_pr + n_estVel, n_pr, n_pr + n_estVel).trace();
			// Cvv = Cvv.extractMatrix(0, _size, 0, _size);
		}
		// Pre-fit residual/innovation
		SimpleMatrix v = new SimpleMatrix(_size, 1, true, innovation);
		SimpleMatrix Cvv_inv = Cvv.invert();
		double globalTq = v.transpose().mult(Cvv_inv).mult(v).get(0);
		ChiSquaredDistribution csd = new ChiSquaredDistribution(_size);
		double alpha = 0.01;
		if (globalTq == 0) {
			throw new Exception("Error: T stat is zero");
		}
		// Detection
		double globalPVal = 1 - csd.cumulativeProbability(globalTq);
//		if(globalPVal < alpha)
//		{
//			System.out.println("\nOutlier Detection Epoch "+ct);
//			System.out.println("Size "+_size);
//		}
		while (globalPVal < alpha && (_size > ((size - n_estVel) / 2))) {
			double max_w = Double.MIN_VALUE;
			int index = -1;
			for (int i = 0; i < _size; i++) {
				SimpleMatrix cv = new SimpleMatrix(_size, 1);
				cv.set(i, 1);
				double w = Math.abs(cv.transpose().mult(Cvv_inv).mult(v).get(0))
						/ Math.sqrt(cv.transpose().mult(Cvv_inv).mult(cv).get(0));
				if (w > max_w) {
					max_w = w;
					index = i;
				}
			}
//			System.out.print(" " + index + " ");
			Measurement meas = null;
			int _index = index;
			if (index >= n_pr) {
				meas = Measurement.Doppler;
				_index -= n_pr;
				n_doppler -= 1;
			} else {
				meas = Measurement.Pseudorange;
				n_pr -= 1;
			}
			ArrayList<Satellite> satList = satMap.get(meas);
			ArrayList<Satellite> testedSatList = testedSatMap.get(meas);
			satList.get(satList.indexOf(testedSatList.remove(_index))).setOutlier(true);
			satMap.put(meas, satList);
			testedSatMap.put(meas, testedSatList);
			_size = n_pr + n_doppler + n_estVel;
			SimpleMatrix R_ = new SimpleMatrix(_size, _size);
			double[][] z_ = new double[_size][1];
			double[][] ze_ = new double[_size][1];
			SimpleMatrix H_ = new SimpleMatrix(_size, cols);
			SimpleMatrix v_ = new SimpleMatrix(_size - n_estVel, 1);
			int j = 0;
			for (int i = 0; i < _size + 1; i++) {
				if (i != index) {
					if (useEstVel) {
						if (i <= n_pr) {
							v_.set(j, v.get(i));
						}
						int l = 0;
						for (int k = 0; k < _size + 1; k++) {
							if (k != index) {
								R_.set(j, l, R.get(i, k));
								l++;
							}	
						}
					} else {
						v_.set(j, v.get(i));
						R_.set(j, j, R.get(i, i));
					}
					z_[j][0] = z[i][0];
					ze_[j][0] = ze[i][0];
					for (int k = 0; k < cols; k++) {
						H_.set(j, k, H.get(i, k));
					}
					j++;
				}
			}
			R = new SimpleMatrix(R_);
			H = new SimpleMatrix(H_);
			z = new double[_size][1];
			ze = new double[_size][1];
			v = new SimpleMatrix(v_);
			for (int i = 0; i < _size; i++) {
				z[i] = z_[i];
				ze[i] = ze_[i];
			}
			_size -= n_estVel;
			Cvv = ((H.mult(P).mult(H.transpose())).plus(R));
			if(useEstVel)
			{
				Cvv = ((H.extractMatrix(0, _size, 0, 3+m).mult(P.extractMatrix(0, 3+m, 0, 3+m)).mult(H.extractMatrix(0, _size, 0, 3+m).transpose())).plus(R.extractMatrix(0, _size, 0, _size)));
			}
			//Cvv = Cvv.extractMatrix(0, _size, 0, _size);
			Cvv_inv = Cvv.invert();
			globalTq = v.transpose().mult(Cvv_inv).mult(v).get(0);
			if (globalTq == 0) {
				throw new Exception("Error: T stat is zero");
			}
			csd = new ChiSquaredDistribution(_size);
			globalPVal = 1 - csd.cumulativeProbability(globalTq);
		}
//		if(ct<150)
//		{
//			System.out.println("\n Epoch "+ct);
//			System.out.println("P  \n"+P.extractMatrix(0, 3+m, 0, 3+m));
//			System.out.println("H  \n"+H.extractMatrix(0, n_pr, 0, 3+m));
//			System.out.println("R  \n"+R.extractMatrix(0, n_pr, 0, n_pr));
//			System.out.println("Cvv  \n"+Cvv);
//			double[] residual = new double[n_pr];
//			for(int i=0;i<n_pr;i++)
//			{
//				residual[i] = z[i][0]-ze[i][0];
//			}
//			
//			System.out.println("Residual  \n"+Arrays.toString(residual));
//
//		}
		return new Object[] { R, H, z, ze };

	}

	private Object[] get_z_ze_res(SimpleMatrix x, HashMap<Measurement, ArrayList<Satellite>> satMap,
			String[] obsvCodeList, SimpleMatrix H, boolean useDoppler) {

		int n_pr = satMap.get(Measurement.Pseudorange).size();
		int n_doppler = 0;
		int m = obsvCodeList.length;
		int size = n_pr;
		double[] estVel = null;
		double[] rxClkDrift = null;

		double[] estPos = new double[] { x.get(0), x.get(1), x.get(2) };
		double[] rxClkOff = new double[m];// in meters
		for (int i = 0; i < m; i++) {
			rxClkOff[i] = x.get(i + 3);
		}
		if (useDoppler) {
			n_doppler = satMap.get(Measurement.Doppler).size();
			size = size + n_doppler;
			estVel = new double[] { x.get(3 + m), x.get(4 + m), x.get(5 + m) };
			rxClkDrift = new double[m];// in meters
			for (int i = 0; i < m; i++) {
				rxClkDrift[i] = x.get(6 + m + i);
			}
		}
		double[][] z = new double[size][1];
		double[][] ze = new double[size][1];
		double[] residual = new double[size];

		for (int i = 0; i < n_pr; i++) {
			Satellite sat = satMap.get(Measurement.Pseudorange).get(i);
			String obsvCode = sat.getObsvCode();
			z[i][0] = sat.getPseudorange();
			ze[i][0] = Math.sqrt(IntStream.range(0, 3).mapToDouble(j -> estPos[j] - sat.getSatEci()[j]).map(j -> j * j)
					.reduce(0, (j, k) -> j + k));
			for (int j = 0; j < m; j++) {
				if (obsvCode.equals(obsvCodeList[j])) {
					ze[i][0] += rxClkOff[j];
				}
			}
			residual[i] = z[i][0] - ze[i][0];
		}
		if (useDoppler) {
			for (int i = 0; i < n_doppler; i++) {
				Satellite sat = satMap.get(Measurement.Doppler).get(i);
				String obsvCode = sat.getObsvCode();
				double rangeRate = sat.getRangeRate();
				SimpleMatrix satVel = new SimpleMatrix(3, 1, true, sat.getSatVel());
				SimpleMatrix A = new SimpleMatrix(1, 3, true, new double[] { -H.get(i + n_pr, 0 + 3 + m),
						-H.get(i + n_pr, 1 + 3 + m), -H.get(i + n_pr, 2 + 3 + m) });

				/*
				 * Observable derived from doppler and satellite velocity, refer Kaplan and
				 * personal notes
				 */
				double dopplerDerivedObs = rangeRate - A.mult(satVel).get(0);
				z[i + n_pr][0] = dopplerDerivedObs;
				// Est doppler derived observable
				double estDopplerDerivedObs = -A.mult(new SimpleMatrix(3, 1, true, estVel)).get(0);
				for (int j = 0; j < m; j++) {
					if (obsvCode.equals(obsvCodeList[j])) {
						estDopplerDerivedObs += rxClkDrift[j];
					}
				}
				ze[i + n_pr][0] = estDopplerDerivedObs;
				residual[i + n_pr] = z[i + n_pr][0] - ze[i + n_pr][0];
			}
		}

		return new Object[] { z, ze, residual };

	}

	private Object[] get_z_ze_res(SimpleMatrix x, ArrayList<Satellite> satList, String[] obsvCodeList, SimpleMatrix H,
			boolean useDoppler, boolean useEstVel, double[] estVel_Z) {
		int n = satList.size();
		int m = obsvCodeList.length;
		int size = n;
		double[] estVel = null;
		double[] rxClkDrift = null;
		double[] estPos = new double[] { x.get(0), x.get(1), x.get(2) };
		double[] rxClkOff = new double[m];// in meters
		for (int i = 0; i < m; i++) {
			rxClkOff[i] = x.get(i + 3);
		}
		if (useDoppler) {
			size = 2 * n;
		} else if (useEstVel) {
			size = n + 3 + m;
		}
		if (useDoppler || useEstVel) {
			estVel = new double[] { x.get(3 + m), x.get(4 + m), x.get(5 + m) };
			rxClkDrift = new double[m];// in meters
			for (int i = 0; i < m; i++) {
				rxClkDrift[i] = x.get(6 + m + i);
			}
		}
		double[][] z = new double[size][1];
		double[][] ze = new double[size][1];
		double[] residual = new double[size];

		for (int i = 0; i < n; i++) {
			final int _i = i;
			Satellite sat = satList.get(i);
			String obsvCode = sat.getObsvCode();
			z[i][0] = sat.getPseudorange();
			ze[i][0] = Math.sqrt(IntStream.range(0, 3).mapToDouble(j -> estPos[j] - satList.get(_i).getSatEci()[j])
					.map(j -> j * j).reduce(0, (j, k) -> j + k));
			for (int j = 0; j < m; j++) {
				if (obsvCode.equals(obsvCodeList[j])) {
					ze[i][0] += rxClkOff[j];
				}
			}
			residual[i] = z[i][0] - ze[i][0];

			if (useDoppler) {
				double rangeRate = sat.getRangeRate();
				SimpleMatrix satVel = new SimpleMatrix(3, 1, true, sat.getSatVel());
				SimpleMatrix A = new SimpleMatrix(1, 3, true,
						new double[] { -H.get(i, 0), -H.get(i, 1), -H.get(i, 2) });

				/*
				 * Observable derived from doppler and satellite velocity, refer Kaplan and
				 * personal notes
				 */
				double dopplerDerivedObs = rangeRate - A.mult(satVel).get(0);
				z[i + n][0] = dopplerDerivedObs;
				// Est doppler derived observable
				double estDopplerDerivedObs = -A.mult(new SimpleMatrix(3, 1, true, estVel)).get(0);
				for (int j = 0; j < m; j++) {
					if (obsvCode.equals(obsvCodeList[j])) {
						estDopplerDerivedObs += rxClkDrift[j];
					}
				}
				ze[i + n][0] = estDopplerDerivedObs;
				residual[i + n] = z[i + n][0] - ze[i + n][0];

			}
		}
		if (useEstVel) {
			for (int i = 0; i < 3; i++) {
				z[i + n][0] = estVel_Z[i];
				ze[i + n][0] = estVel[i];
				residual[i + n] = z[i + n][0] - ze[i + n][0];
			}
			for (int i = 0; i < m; i++) {
				z[i + 3 + n][0] = estVel_Z[3 + i];
				ze[i + 3 + n][0] = rxClkDrift[i];
				residual[i + 3 + n] = z[i + 3 + n][0] - ze[i + 3 + n][0];
			}
		}

		return new Object[] { z, ze, residual };
	}

	public TreeMap<Long, HashMap<Measurement, double[]>> getInnovationMap() {
		return innovationMap;
	}

	public TreeMap<Long, HashMap<Measurement, double[]>> getResidualMap() {
		return residualMap;
	}

	public TreeMap<Long, HashMap<Measurement, Double>> getPostVarOfUnitWMap() {
		return postVarOfUnitWMap;
	}

	public TreeMap<Long, HashMap<State, SimpleMatrix>> getErrCovMap() {
		return errCovMap;
	}

	public ArrayList<double[]> getRedundancyList(Measurement meas) {
		return redundancyMap.get(meas);
	}

	public TreeMap<Long, HashMap<Measurement, Long>> getSatCountMap() {
		return satCountMap;
	}

	public TreeMap<Long, HashMap<Measurement, ArrayList<Satellite>>> getSatListMap() {
		return satListMap;
	}
}
