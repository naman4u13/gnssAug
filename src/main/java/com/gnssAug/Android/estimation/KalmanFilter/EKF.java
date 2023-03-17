package com.gnssAug.Android.estimation.KalmanFilter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.stream.IntStream;

import org.ejml.dense.row.MatrixFeatures_DDRM;
import org.ejml.simple.SimpleMatrix;

import com.gnssAug.Android.estimation.LinearLeastSquare;
import com.gnssAug.Android.estimation.KalmanFilter.Models.Flag;
import com.gnssAug.Android.estimation.KalmanFilter.Models.KFconfig;
import com.gnssAug.Android.models.Satellite;
import com.gnssAug.utility.LatLonUtil;
import com.gnssAug.Android.constants.Measurement;
import com.gnssAug.Android.constants.State;

public class EKF {
	private final double SpeedofLight = 299792458;
	private KFconfig kfObj;
	private double[] innovation;
	private double[] temp_innovation;
	private TreeMap<Long,HashMap<State,double[]>> innovationMap;
	private TreeMap<Long, HashMap<Measurement,double[]>> residualMap;
	private TreeMap<Long, HashMap<State,double[]>> measNoiseMap;
	private TreeMap<Long, HashMap<Measurement,Double>> postVarOfUnitWMap;
	// Posteriori Err Cov
	private TreeMap<Long, HashMap<State,SimpleMatrix>> errCovMap;
	private HashMap<Measurement,ArrayList<double[]>> redundancyMap;
	// Satellite Count
	private TreeMap<Long, Long> satCountMap;
	private TreeMap<Long, ArrayList<Satellite>> satListMap;

	public EKF() {
		kfObj = new KFconfig();
	}

	public TreeMap<Long, double[]> process(TreeMap<Long, ArrayList<Satellite>> SatMap, ArrayList<Long> timeList,
			Flag flag, boolean useDoppler, boolean useIGS,String[] obsvCodeList, boolean doAnalyze, boolean doTest, boolean outlierAnalyze) throws Exception {

		int n = 0;
		int m = obsvCodeList.length;
		/* constant position model - state vector(n=5) -> (x,y,z,cdt,cdt_dot) */
		/*
		 * constant velocity model - state vector(n=8) ->
		 * (x,y,z,cdt,x_dot,y_dot,z_dot,cdt_dot)
		 */
		if (flag == Flag.POSITION) {
			n = 3+(2*m);
		} else if (flag == Flag.VELOCITY) {
			n = 6+(2*m);
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
		double[] intialVel = LinearLeastSquare.getEstVel(SatMap.firstEntry().getValue(), true, intialECEF, useIGS);
		IntStream.range(0, 3+m).forEach(i -> x[i][0] = intialECEF[i]);
		IntStream.range(0, 3+m).forEach(i -> P[i][i] = 100);
		if (n ==  3+(2*m)) {
			IntStream.range( 3+m,3+(2*m)).forEach(i -> P[i][i] = 1e5);
		} else {
			IntStream.range(3+m,6+(2*m) ).forEach(i -> x[i][0] = intialVel[i-(3+m)]);
			IntStream.range(3+m, 6+m).forEach(i -> P[i][i] = 1);
			IntStream.range(6+m, 6+(2*m)).forEach(i -> P[i][i] = 1e5);
		}

		kfObj.setState_ProcessCov(x, P);
		if (doAnalyze) {
			innovationMap = new TreeMap<Long, HashMap<State,double[]>>();
			errCovMap = new TreeMap<Long, HashMap<State,SimpleMatrix>>();
			residualMap = new TreeMap<Long, HashMap<Measurement,double[]>>();
			postVarOfUnitWMap = new TreeMap<Long, HashMap<Measurement,Double>>();
			redundancyMap = new HashMap<Measurement,ArrayList<double[]>>();
			satCountMap = new TreeMap<Long, Long>();
			satListMap = new TreeMap<Long, ArrayList<Satellite>>();
			measNoiseMap = new TreeMap<Long, HashMap<State,double[]>>();
		}
		// Begin iteration or recursion
		return iterate(SatMap, timeList, flag, useDoppler,obsvCodeList,doAnalyze,doTest,outlierAnalyze);

	}

	private TreeMap<Long, double[]> iterate(TreeMap<Long, ArrayList<Satellite>> SatMap, ArrayList<Long> timeList,
			Flag flag, boolean useDoppler,String[] obsvCodeList, boolean doAnalyze, boolean doTest, boolean outlierAnalyze) throws Exception {
		TreeMap<Long, double[]> estStateMap = new TreeMap<Long, double[]>();

		int m = obsvCodeList.length;
		long time = timeList.get(0);
		// Start from 2nd epoch
		for (int i = 1; i < timeList.size(); i++) {
			long currentTime = timeList.get(i);
			ArrayList<Satellite> satList = SatMap.get(currentTime);
			double deltaT = (currentTime - time) / 1e3;
			// Perform Predict and Update
			runFilter(deltaT, satList, flag, useDoppler,obsvCodeList);
			// Fetch Posteriori state estimate and estimate error covariance matrix
			SimpleMatrix x = kfObj.getState();
			SimpleMatrix P = kfObj.getCovariance();
			int n = x.numRows();
			double[] estState = new double[n];
			for(int j=0;j<n;j++)
			{
				estState[j] = x.get(j);
			}
			// Add position estimate to the list
			estStateMap.put(currentTime, estState);
			if (doAnalyze) {
			
				// Convert to ENU frame
				SimpleMatrix R = new SimpleMatrix(n, n);
				SimpleMatrix rotMat = new SimpleMatrix(LatLonUtil.getEcef2EnuRotMat(new double[] {estState[0],estState[1],estState[2]}));
				R.insertIntoThis(0, 0, rotMat);
				IntStream.range(3, 3+m).forEach(j -> R.set(j, j, 1));
				if(flag==Flag.VELOCITY)
				{
					R.insertIntoThis(3+m, 3+m, rotMat);
					IntStream.range(3+(2*m), 6+(2*m)).forEach(j -> R.set(j, j, 1));
				}
				else
				{
					IntStream.range(3+m, 3+(2*m)).forEach(j -> R.set(j, j, 1));
				}
				SimpleMatrix errCov = R.mult(P).mult(R.transpose());
				errCovMap.computeIfAbsent(currentTime, j->new HashMap<State,SimpleMatrix>()).put(State.Position, errCov.extractMatrix(0, 3+m, 0, 3+m));
				int satCount = satList.size();
				innovationMap.computeIfAbsent(currentTime, j->new HashMap<State,double[]>()).put(State.Position, Arrays.copyOfRange(temp_innovation, 0, satCount));
				if(flag==Flag.VELOCITY)
				{
					errCovMap.get(currentTime).put(State.Velocity, errCov.extractMatrix(3+m, 6+(2*m), 3+m, 6+(2*m)));
					innovationMap.get(currentTime).put(State.Velocity, Arrays.copyOfRange(temp_innovation, satCount,2*satCount));
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

	private void runFilter(double deltaT, ArrayList<Satellite> satList, Flag flag, boolean useDoppler,String[] obsvCodeList)
			throws Exception {

		// Satellite count
		int n = satList.size();
		int m = obsvCodeList.length;
		// Assign Q and F matrix
		kfObj.config(deltaT, flag,m);
		kfObj.predict();

		SimpleMatrix x = kfObj.getState();
		double[] estPos = new double[] { x.get(0), x.get(1), x.get(2) };
		double[] rxClkOff = new double[m];// in meters
		for (int i = 0; i < m; i++) {
			rxClkOff[i] = x.get(i + 3);
		}
		double[] estVel = null;
		double[] rxClkDrift = null;
		/*
		 * H is the Jacobian matrix of partial derivatives Observation StateModel(h) of
		 * with respect to x
		 */
		double[][] H = getJacobian(satList, estPos, flag, useDoppler,obsvCodeList);
		// Measurement vector
		double[][] z = null;
		// Estimated Measurement vector
		double[][] ze = null;
		// Measurement Noise
		double[][] R = null;
		if (useDoppler) {
			z = new double[2 * n][1];
			ze = new double[2 * n][1];
			R = new double[2 * n][2 * n];
			estVel = new double[] { x.get(3+m), x.get(4+m), x.get(5+m) };
			rxClkDrift = new double[m];// in meters
			for (int i = 0; i < m; i++) {
				rxClkDrift[i] = x.get(6+m+i);
			}
		} else {
			z = new double[n][1];
			ze = new double[n][1];
			R = new double[n][n];
		}

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
			R[i][i] = Math.pow(sat.getReceivedSvTimeUncertaintyNanos() * SpeedofLight * 1e-9, 2);

			if (useDoppler) {
				double rangeRate = sat.getRangeRate();
				SimpleMatrix satVel = new SimpleMatrix(3, 1, true, sat.getSatVel());
				SimpleMatrix A = new SimpleMatrix(1, 3, true, new double[] { -H[i][0], -H[i][1], -H[i][2] });

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
				R[i + n][i + n] = Math.pow(sat.getReceivedSvTimeUncertaintyNanos() * SpeedofLight * 1e-9, 2) / 10;
			}

		}

		// Perform Update Step
		kfObj.update(z, R, ze, H);

	}

	private double[][] getJacobian(ArrayList<Satellite> satList, double[] estECEF,Flag flag, boolean useDoppler,String[] obsvCodeList) {
		int n = satList.size();
		int m = obsvCodeList.length;
		int rows = n;
		int stateN = 3+(2*m);
		if (flag==Flag.VELOCITY && useDoppler) {
			rows = 2 * n;
			stateN = 6+(2*m);
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
					if (useDoppler)
					{
						H[i+n][6 + m + j] = 1;
					}
				}
			}
			if (useDoppler) {
				IntStream.range(0, 3).forEach(j -> H[_i + n][j + 3 + m] = -LOS[j] / GR);
				
			}
		}

		return H;

	}
	
	private void performAnalysis(SimpleMatrix X, ArrayList<Satellite> testedSatList, ArrayList<Satellite> satList,
			SimpleMatrix R, SimpleMatrix H, SimpleMatrix priorP, long currentTime,  int n,
			String[] obsvCodeList, boolean doTest, boolean outlierAnalyze,Flag flag) {

		int _n = testedSatList.size();
		double[] residual = (double[]) get_z_ze_res(X, testedSatList, obsvCodeList)[2];
		double[] pr_res = Arrays.copyOfRange(residual, 0, _n);
		// Post-fit residual
		SimpleMatrix e_post_hat_pr = new SimpleMatrix(_n, 1, true, pr_res);
		SimpleMatrix Cyy_inv_pr = R.extractMatrix(0, _n, 0, _n).invert();
		
		// Compute Redundancies
		SimpleMatrix K = kfObj.getKalmanGain();
		SimpleMatrix HK = H.mult(K);
		SimpleMatrix phi = kfObj.getPhi();
		SimpleMatrix Cvv = kfObj.getCvv();
		SimpleMatrix HtCvvInvH = H.transpose().mult(Cvv.invert()).mult(H);
		SimpleMatrix Q = kfObj.getQ();
		double rX = phi.mult(priorP).mult(phi.transpose()).mult(HtCvvInvH).trace();
		double rW = Q.mult(HtCvvInvH).trace();
		int m = HK.numRows();
		SimpleMatrix tempMat = SimpleMatrix.identity(m).minus(HK);
		double rZ_pr = tempMat.extractMatrix(0, m/2, 0, m/2).trace();
		double rZ_doppler = tempMat.extractMatrix(m/2,m, m/2,m).trace();
		double rZ = SimpleMatrix.identity(HK.numRows()).minus(HK).trace();
		double rSum = rX + rW + rZ;
		if ((2*_n) - rSum > 0.01) {
			System.err.println("FATAL ERROR: Redundancy sum is wrong");
		}
		double postVarOfUnitW_pr = e_post_hat_pr.transpose().mult(Cyy_inv_pr).mult(e_post_hat_pr).get(0) / rZ_pr;
		
		if(flag == Flag.VELOCITY) {
			double[] doppler_res = Arrays.copyOfRange(residual, _n, 2*_n);
			SimpleMatrix e_post_hat_doppler = new SimpleMatrix(_n, 1, true, doppler_res);
			SimpleMatrix Cyy_inv_doppler = R.extractMatrix(_n, 2*_n, _n,2*_n).invert();
			double postVarOfUnitW_doppler = e_post_hat_doppler.transpose().mult(Cyy_inv_doppler).mult(e_post_hat_doppler).get(0) / rZ_doppler;
			postVarOfUnitWMap.computeIfAbsent(currentTime, k->new HashMap<Measurement,Double>()).put(Measurement.Doppler, postVarOfUnitW_doppler);
			residualMap.computeIfAbsent(currentTime,  k->new HashMap<Measurement,double[]>()).put(Measurement.Doppler, doppler_res);
			}
		
		redundancyMap.computeIfAbsent(Measurement.Pseudorange, k->new ArrayList<double[]>()).add(new double[] { _n, rSum-rZ+rZ_pr, rX, rW, rZ_pr });
		redundancyMap.computeIfAbsent(Measurement.Doppler, k->new ArrayList<double[]>()).add(new double[] { _n, rSum-rZ+rZ_doppler, rX, rW, rZ_doppler });
		postVarOfUnitWMap.computeIfAbsent(currentTime, k->new HashMap<Measurement,Double>()).put(Measurement.Pseudorange, postVarOfUnitW_pr);
		residualMap.computeIfAbsent(currentTime,  k->new HashMap<Measurement,double[]>()).put(Measurement.Pseudorange, pr_res);
		
		if (doTest == true) {
			_n = n - _n;
		}
		satCountMap.put(currentTime, (long) _n);
		if (outlierAnalyze) {
			satListMap.put(currentTime, satList);
		} else {
			satListMap.put(currentTime, testedSatList);
		}
		

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
}
