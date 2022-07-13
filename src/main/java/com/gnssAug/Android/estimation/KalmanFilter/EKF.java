package com.gnssAug.Android.estimation.KalmanFilter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeMap;
import java.util.stream.IntStream;

import org.ejml.dense.row.MatrixFeatures_DDRM;
import org.ejml.simple.SimpleMatrix;

import com.gnssAug.Android.estimation.LinearLeastSquare;
import com.gnssAug.Android.estimation.KalmanFilter.Models.Flag;
import com.gnssAug.Android.estimation.KalmanFilter.Models.KFconfig;
import com.gnssAug.Android.models.Satellite;

public class EKF {
	private final double SpeedofLight = 299792458;
	private KFconfig kfObj;
	private double prObsNoiseVar;

	public EKF() {
		kfObj = new KFconfig();
	}

	public TreeMap<Long, double[]> process(TreeMap<Long, ArrayList<Satellite>> SatMap, ArrayList<Long> timeList,
			Flag flag, boolean useDoppler) throws Exception {

		int n = 0;
		/* constant position model - state vector(n=5) -> (x,y,z,cdt,cdt_dot) */
		/*
		 * constant velocity model - state vector(n=8) ->
		 * (x,y,z,cdt,x_dot,y_dot,z_dot,cdt_dot)
		 */
		if (flag == Flag.POSITION) {
			n = 5;
		} else if (flag == Flag.VELOCITY) {
			n = 8;
		}
		double[][] x = new double[n][1];
		double[][] P = new double[n][n];
		/*
		 * state XYZ is intialized WLS generated Rx position estimated using first epoch
		 * data, a-priori estimate error covariance matrix for state XYZ is therefore
		 * assigned 25 m^2 value. Other state variables are assigned infinite(big)
		 * variance
		 */
		double[] intialECEF = LinearLeastSquare.process(SatMap.firstEntry().getValue(), true);
		IntStream.range(0, 3).forEach(i -> x[i][0] = intialECEF[i]);
		x[3][0] = SpeedofLight * intialECEF[3];
		IntStream.range(0, 4).forEach(i -> P[i][i] = 100);
		if (n == 5) {
			P[4][4] = 1e13;
		} else {
			IntStream.range(4, 7).forEach(i -> P[i][i] = 1);
			P[7][7] = 1e3;
		}

		kfObj.setState_ProcessCov(x, P);
		// Begin iteration or recursion
		return iterate(SatMap, timeList, flag, useDoppler);

	}

	private TreeMap<Long, double[]> iterate(TreeMap<Long, ArrayList<Satellite>> SatMap, ArrayList<Long> timeList,
			Flag flag, boolean useDoppler) throws Exception {
		TreeMap<Long, double[]> ecefMap = new TreeMap<Long, double[]>();
		long time = timeList.get(0);
		// Start from 2nd epoch
		for (int i = 1; i < timeList.size(); i++) {
			long currentTime = timeList.get(i);
			ArrayList<Satellite> satList = SatMap.get(currentTime);
			double deltaT = (currentTime - time) / 1e3;
			// Perform Predict and Update
			runFilter(deltaT, satList, flag, useDoppler);
			// Fetch Posteriori state estimate and estimate error covariance matrix
			SimpleMatrix x = kfObj.getState();
			SimpleMatrix P = kfObj.getCovariance();
			double[] estECEF = new double[] { x.get(0), x.get(1), x.get(2) };
			// Add position estimate to the list
			ecefMap.put(currentTime, estECEF);
			/*
			 * Check whether estimate error covariance matrix is positive semidefinite
			 * before further proceeding
			 */
			if (!MatrixFeatures_DDRM.isPositiveSemidefinite(P.getMatrix())) {

				throw new Exception("PositiveDefinite test Failed");
			}
			time = currentTime;

		}

		return ecefMap;
	}

	private void runFilter(double deltaT, ArrayList<Satellite> satList, Flag flag, boolean useDoppler) {

		// Satellite count
		int n = satList.size();

		// Assign Q and F matrix
		kfObj.config(deltaT, flag);
		kfObj.predict();

		SimpleMatrix x = kfObj.getState();
		double[] estPos = new double[] { x.get(0), x.get(1), x.get(2) };
		double rxClkOff = x.get(3);// in meters
		double[] estVel = new double[] { x.get(4), x.get(5), x.get(6) };
		double rxClkDrift = x.get(7);// in meters

		/*
		 * H is the Jacobian matrix of partial derivatives Observation StateModel(h) of
		 * with respect to x
		 */
		double[][] H = getJacobian(satList, estPos, x.numRows(), useDoppler);
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
		} else {
			z = new double[n][1];
			ze = new double[n][1];
			R = new double[n][n];
		}

		for (int i = 0; i < n; i++) {
			final int _i = i;
			Satellite sat = satList.get(i);
			z[i][0] = sat.getPseudorange();
			ze[i][0] = Math.sqrt(IntStream.range(0, 3).mapToDouble(j -> estPos[j] - satList.get(_i).getSatEci()[j])
					.map(j -> j * j).reduce(0, (j, k) -> j + k)) + rxClkOff;
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
				double estDopplerDerivedObs = -A.mult(new SimpleMatrix(3, 1, true, estVel)).get(0) + rxClkDrift;
				ze[i + n][0] = estDopplerDerivedObs;
				R[i + n][i + n] = Math.pow(sat.getPseudorangeRateUncertaintyMetersPerSecond(), 2);
			}

		}

		// Perform Update Step
		kfObj.update(z, R, ze, H);

	}

	private double[][] getJacobian(ArrayList<Satellite> satList, double[] estECEF, int stateN, boolean useDoppler) {
		int n = satList.size();
		int rows = n;
		if (stateN == 8 && useDoppler) {
			rows = 2 * n;
		}
		double[][] H = new double[rows][stateN];

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
			if (useDoppler) {
				IntStream.range(0, 3).forEach(j -> H[_i + n][j + 4] = -LOS[j] / GR);
				H[i + n][7] = 1;
			}
		}

		return H;

	}
}
