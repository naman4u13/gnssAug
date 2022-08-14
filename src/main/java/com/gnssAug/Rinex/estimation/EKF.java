package com.gnssAug.Rinex.estimation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.stream.IntStream;

import org.ejml.dense.row.MatrixFeatures_DDRM;
import org.ejml.simple.SimpleMatrix;

import com.gnssAug.Android.estimation.KalmanFilter.Models.KFconfig;
import com.gnssAug.Rinex.models.Satellite;

public class EKF {

	private final double SpeedofLight = 299792458;
	private KFconfig kfObj;
	private double prObsNoiseVar;

	public EKF() {
		kfObj = new KFconfig();
	}

	public TreeMap<Long, double[]> process(TreeMap<Long, ArrayList<Satellite>> SatMap, HashMap<String, double[]> rxPCO,
			ArrayList<Long> timeList) throws Exception {

		int n = 5;

		/*
		 * state XYZ is intialized WLS generated Rx position estimated using first epoch
		 * data, a-priori estimate error covariance matrix for state XYZ is therefore
		 * assigned 25 m^2 value. Other state variables are assigned infinite(big)
		 * variance
		 */
		double[][] x = new double[n][1];
		double[][] P = new double[n][n];
		double[] intialECEF = LinearLeastSquare.process(SatMap.firstEntry().getValue(), rxPCO, true);
		IntStream.range(0, 3).forEach(i -> x[i][0] = intialECEF[i]);
		x[3][0] = SpeedofLight * intialECEF[3];
		IntStream.range(0, 4).forEach(i -> P[i][i] = 10);
		P[4][4] = 1e5;

		kfObj.setState_ProcessCov(x, P);
		// Begin iteration or recursion
		return iterate(SatMap, rxPCO, timeList);

	}

	private TreeMap<Long, double[]> iterate(TreeMap<Long, ArrayList<Satellite>> SatMap, HashMap<String, double[]> rxPCO,
			ArrayList<Long> timeList) throws Exception {
		TreeMap<Long, double[]> estStateMap = new TreeMap<Long, double[]>();

		long time = timeList.get(0);
		// Start from 2nd epoch
		for (int i = 1; i < timeList.size(); i++) {
			long currentTime = timeList.get(i);
			ArrayList<Satellite> satList = SatMap.get(currentTime);
			double deltaT = (currentTime - time) / 1e3;
			// Perform Predict and Update
			runFilter(deltaT, satList, rxPCO);
			// Fetch Posteriori state estimate and estimate error covariance matrix
			SimpleMatrix x = kfObj.getState();
			SimpleMatrix P = kfObj.getCovariance();

			double[] estState = new double[] { x.get(0), x.get(1), x.get(2) };

			// Add position estimate to the list
			estStateMap.put(currentTime, estState);

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

	private void runFilter(double deltaT, ArrayList<Satellite> satList, HashMap<String, double[]> rxPCO)
			throws Exception {

		// Satellite count
		int n = satList.size();

		// Assign Q and F matrix
		kfObj.configIGS(deltaT);
		kfObj.predict();

		SimpleMatrix x = kfObj.getState();
		double[] estPos = new double[] { x.get(0), x.get(1), x.get(2) };
		double rxClkOff = x.get(3);// in meters

		/*
		 * H is the Jacobian matrix of partial derivatives Observation StateModel(h) of
		 * with respect to x
		 */
		double[][] H = getJacobian(satList, estPos);
		// Measurement vector
		double[][] z = null;
		// Estimated Measurement vector
		double[][] ze = null;
		// Measurement Noise
		double[][] R = null;

		z = new double[n][1];
		ze = new double[n][1];
		R = new double[n][n];

		for (int i = 0; i < n; i++) {
			final int _i = i;
			Satellite sat = satList.get(i);
			String obsvCode = sat.getSSI() + "" + sat.getFreqID() + "C";
			double[] pco = rxPCO.get(obsvCode);
			double[] rxAPC = IntStream.range(0, 3).mapToDouble(j -> estPos[j] + pco[j]).toArray();
			z[i][0] = sat.getPseudorange();
			ze[i][0] = Math.sqrt(IntStream.range(0, 3).mapToDouble(j -> rxAPC[j] - satList.get(_i).getSatEci()[j])
					.map(j -> j * j).reduce(0, (j, k) -> j + k)) + rxClkOff;
			R[i][i] = 4;
		}

		// Perform Update Step
		kfObj.update(z, R, ze, H);

	}

	private double[][] getJacobian(ArrayList<Satellite> satList, double[] estECEF) {
		int n = satList.size();
		int rows = n;
		double[][] H = new double[rows][5];

		for (int i = 0; i < n; i++) {
			Satellite sat = satList.get(i);
			// Line of Sight vector
			// Its not really a ECI, therefore don't get confused
			double[] LOS = IntStream.range(0, 3).mapToDouble(j -> sat.getSatEci()[j] - estECEF[j]).toArray();
			// Geometric Range
			double GR = Math.sqrt(Arrays.stream(LOS).map(j -> j * j).reduce(0.0, (j, k) -> j + k));
			// Converting LOS to unit vector
			final int _i = i;
			IntStream.range(0, 3).forEach(j -> H[_i][j] = -LOS[j] / GR);
			H[i][3] = 1;

		}

		return H;

	}

}
