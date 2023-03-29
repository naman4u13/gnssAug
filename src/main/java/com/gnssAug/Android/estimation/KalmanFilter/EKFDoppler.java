package com.gnssAug.Android.estimation.KalmanFilter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeMap;
import java.util.stream.IntStream;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.ejml.dense.row.MatrixFeatures_DDRM;
import org.ejml.simple.SimpleMatrix;

import com.gnssAug.Android.constants.Measurement;
import com.gnssAug.Android.estimation.LinearLeastSquare;
import com.gnssAug.Android.estimation.KalmanFilter.Models.KFconfig;
import com.gnssAug.Android.models.Satellite;
import com.gnssAug.utility.LatLonUtil;
import com.gnssAug.utility.Matrix;
import com.gnssAug.utility.Vector;
import com.gnssAug.utility.Weight;

public class EKFDoppler {


	private final double SpeedofLight = 299792458;
	private KFconfig kfObj;

	private double[] innovation;
	private double[] temp_innovation;
	private TreeMap<Long, double[]> innovationMap;
	private TreeMap<Long, double[]> residualMap;
	private TreeMap<Long, double[]> measNoiseMap;
	private TreeMap<Long, Double> postVarOfUnitWMap;
	// Posteriori Err Cov
	private TreeMap<Long, SimpleMatrix> errCovMap;
	private ArrayList<double[]> redundancyList;
	// Satellite Count
	private TreeMap<Long, Long> satCountMap;
	private TreeMap<Long, ArrayList<Satellite>> satListMap;
	final static private double priorVarOfUnitW = Math.pow(13.9,2);
	
	static private double[] prevVel;
	static private SimpleMatrix prev_Cxx_dot_hat;
	public EKFDoppler() {
		kfObj = new KFconfig();
	}

	public TreeMap<Long, double[]> process(TreeMap<Long, ArrayList<Satellite>> SatMap, ArrayList<Long> timeList,
			boolean useIGS, String[] obsvCodeList, boolean doAnalyze, boolean doTest, boolean outlierAnalyze)
			throws Exception {

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
		return iterate(X, SatMap, timeList, useIGS, obsvCodeList, doAnalyze, doTest, outlierAnalyze);

	}

	TreeMap<Long, double[]> iterate(SimpleMatrix X, TreeMap<Long, ArrayList<Satellite>> SatMap,
			ArrayList<Long> timeList, boolean useIGS, String[] obsvCodeList, boolean doAnalyze, boolean doTest,
			boolean outlierAnalyze) throws Exception {
		TreeMap<Long, double[]> estStateMap = new TreeMap<Long, double[]>();
		int m = obsvCodeList.length;
		int x_size = 3 + m;
		long time = timeList.get(0);
		prevVel = LinearLeastSquare.getEstVel(SatMap.get(time), false, true, true, false,
				new double[] { X.get(0), X.get(1), X.get(2) }, useIGS);
		prev_Cxx_dot_hat = LinearLeastSquare.getCxx_hat(Measurement.Doppler,"ECEF");
		// Start from 2nd epoch
		for (int i = 1; i < timeList.size(); i++) {
			long currentTime = timeList.get(i);
			double deltaT = (currentTime - time) / 1e3;
			ArrayList<Satellite> satList = SatMap.get(currentTime);
			SimpleMatrix Cxx_dot_hat = predictTotalState(X, satList, deltaT, useIGS);
			runFilter(X, currentTime, deltaT, satList, obsvCodeList,Cxx_dot_hat, doAnalyze, doTest, outlierAnalyze,useIGS);
			SimpleMatrix P = kfObj.getCovariance();
			double[] estState = new double[x_size];
			IntStream.range(0, x_size).forEach(j -> estState[j] = X.get(j));
			// Add position estimate to the list
			estStateMap.put(currentTime, estState);
			if (doAnalyze) {
				// Convert to ENU frame
				SimpleMatrix R = new SimpleMatrix(x_size, x_size);
				R.insertIntoThis(0, 0, new SimpleMatrix(LatLonUtil.getEcef2EnuRotMat(estState)));
				IntStream.range(3, x_size).forEach(j -> R.set(j, j, 1));

				SimpleMatrix errCov = R.mult(P).mult(R.transpose());
				errCovMap.put(currentTime, errCov);
				innovationMap.put(currentTime, innovation);
			}
			if (!MatrixFeatures_DDRM.isPositiveSemidefinite(P.getMatrix())) {
				throw new Exception("PositiveDefinite test Failed");
			}
			time = currentTime;
		}
		return estStateMap;
	}

	private SimpleMatrix predictTotalState(SimpleMatrix X, ArrayList<Satellite> satList, double deltaT, boolean useIGS)
			throws Exception {

		
		double[] vel = LinearLeastSquare.getEstVel(satList, false, true, true, false,
				new double[] { X.get(0), X.get(1), X.get(2) }, useIGS);
		double[] avg_vel = new double[vel.length];
		for (int i = 0; i < vel.length; i++) {
			avg_vel[i] = (vel[i]+prevVel[i])*0.5;
			X.set(i, X.get(i) + (avg_vel[i]* deltaT));
		}
		SimpleMatrix Cxx_dot_hat = LinearLeastSquare.getCxx_hat(Measurement.Doppler,"ECEF");
		SimpleMatrix avg_Cxx_dot_hat =  Cxx_dot_hat.plus(prev_Cxx_dot_hat).scale(0.25);
		prevVel = Arrays.copyOf(vel, vel.length);
		prev_Cxx_dot_hat = new SimpleMatrix(Cxx_dot_hat);
		return avg_Cxx_dot_hat;
	}
	
	// Innovation Based Testing
	private void runFilter(SimpleMatrix X, long currentTime, double deltaT, ArrayList<Satellite> satList,
			String[] obsvCodeList,SimpleMatrix Cxx_dot_hat, boolean doAnalyze, boolean doTest, boolean outlierAnalyze,boolean useIGS) throws Exception {

		boolean isWeighted = false;
		boolean useAndroidW = false;
		// Satellite count
		int n = satList.size();
		int m = obsvCodeList.length;
		
		// Last update VC-matrix
		SimpleMatrix priorP = new SimpleMatrix(kfObj.getCovariance());
		SimpleMatrix priorX = new SimpleMatrix(X);
		// Assign Q and F matrix
		kfObj.configDoppler(deltaT, Cxx_dot_hat, m,X);
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
				LinearLeastSquare.getEstPos(satList, true, true, false, false, useIGS);
				SimpleMatrix Cyy = LinearLeastSquare.getCyy_updated(Measurement.Pseudorange);
				_R = Matrix.matrix2Array(Cyy);
//				SimpleMatrix Cyy = Weight.getNormCyy(satList, priorVarOfUnitW);
//				_R = Matrix.matrix2Array(Cyy);
			}
		} else {
			for (int i = 0; i < n; i++) {
				_R[i][i] = priorVarOfUnitW;
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
			Object[] params = performTesting(R, H, n, m, satList, testedSatList, z, ze);
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
			performAnalysis(X, testedSatList, satList, R, H, priorP, currentTime,  n, obsvCodeList,
					doTest, outlierAnalyze);
		}
		if (doTest && outlierAnalyze) {
			kfObj.setState_ProcessCov(new SimpleMatrix(3 + m, 1), priorP);
			X = new SimpleMatrix(priorX);
			kfObj.predict();
			Object[] params = performTesting(R, H, n, m, satList, testedSatList, z, ze);
			R = (SimpleMatrix) params[0];
			H = (SimpleMatrix) params[1];
			z = (double[][]) params[2];
			ze = (double[][]) params[3];
			kfObj.update(z, R, ze, H);
			for (int i = 0; i < 3 + m; i++) {
				X.set(i, X.get(i) + x.get(i));
				x.set(i, 0);
			}
		}

	}
	
	// Residual Based Testing
		private void runFilter2(SimpleMatrix X, long currentTime, double deltaT, ArrayList<Satellite> satList,
				String[] obsvCodeList, boolean doAnalyze, boolean doTest, boolean outlierAnalyze,boolean useIGS) throws Exception {

			boolean isWeighted = true;
			boolean useAndroidW = false;
			// Satellite count
			int n = satList.size();
			int m = obsvCodeList.length;
			SimpleMatrix Cxx_dot_hat = LinearLeastSquare.getCxx_hat_updated(Measurement.Doppler,"ECEF");
			// Last update VC-matrix
			SimpleMatrix priorP = new SimpleMatrix(kfObj.getCovariance());
			SimpleMatrix priorX = new SimpleMatrix(X);
			// Assign Q and F matrix
			kfObj.configDoppler(deltaT, Cxx_dot_hat, m,X);
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
					LinearLeastSquare.getEstPos(satList, true, true, false, false, useIGS);
					SimpleMatrix Cyy = LinearLeastSquare.getCyy_updated(Measurement.Pseudorange);
					_R = Matrix.matrix2Array(Cyy);
//					SimpleMatrix Cyy = Weight.getNormCyy(satList, priorVarOfUnitW);
//					_R = Matrix.matrix2Array(Cyy);
				}
			} else {
				for (int i = 0; i < n; i++) {
					//_R[i][i] = priorVarOfUnitW;
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
			
			kfObj.update(z, R, ze, H);
			if (doAnalyze && outlierAnalyze) {
				performAnalysis(X, testedSatList, satList, R, H, priorP, currentTime,  n, obsvCodeList,
						doTest, outlierAnalyze);
			}
			
			
			

		}

	private void performAnalysis(SimpleMatrix X, ArrayList<Satellite> testedSatList, ArrayList<Satellite> satList,
			SimpleMatrix R, SimpleMatrix H, SimpleMatrix priorP, long currentTime,  int n,
			String[] obsvCodeList, boolean doTest, boolean outlierAnalyze) {

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
			ArrayList<Satellite> testedSatList, double[][] z, double[][] ze) throws Exception {
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
					for (int k = 0; k < 3 +  m; k++) {
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
			csd = new ChiSquaredDistribution(_n - 1);
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

	public TreeMap<Long, double[]> getInnovationMap() {
		return innovationMap;
	}

	public TreeMap<Long, SimpleMatrix> getErrCovMap() {
		return errCovMap;
	}

	public TreeMap<Long, double[]> getResidualMap() {
		return residualMap;
	}

	public TreeMap<Long, Double> getPostVarOfUnitWMap() {
		return postVarOfUnitWMap;
	}

	public ArrayList<double[]> getRedundancyList() {
		return redundancyList;
	}

	public TreeMap<Long, Long> getSatCountMap() {
		return satCountMap;
	}

	public TreeMap<Long, ArrayList<Satellite>> getSatListMap() {
		return satListMap;
	}

	public TreeMap<Long, double[]> getMeasNoiseMap() {
		return measNoiseMap;
	}
}
