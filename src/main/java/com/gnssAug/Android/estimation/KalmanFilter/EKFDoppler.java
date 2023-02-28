package com.gnssAug.Android.estimation.KalmanFilter;

import java.util.ArrayList;
import java.util.TreeMap;
import java.util.stream.IntStream;

import org.ejml.dense.row.MatrixFeatures_DDRM;
import org.ejml.simple.SimpleMatrix;

import com.gnssAug.Android.constants.Measurement;
import com.gnssAug.Android.estimation.LinearLeastSquare;
import com.gnssAug.Android.estimation.KalmanFilter.Models.KFconfig;
import com.gnssAug.Android.models.Satellite;
import com.gnssAug.utility.Matrix;
import com.gnssAug.utility.Vector;
import com.gnssAug.utility.Weight;

public class EKFDoppler {

	private final double pseudorange_priorVarOfUnitW =7.4662;
	private final double SpeedofLight = 299792458;
	private KFconfig kfObj;

	public EKFDoppler() {
		kfObj = new KFconfig();
	}

	public TreeMap<Long, double[]> process(TreeMap<Long, ArrayList<Satellite>> SatMap, ArrayList<Long> timeList,
			boolean useIGS, String[] obsvCodeList) throws Exception {

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
		IntStream.range(0, 3 + m).forEach(i -> _X[i][0] = intialECEF[i]);

		// Total State
		SimpleMatrix X = new SimpleMatrix(_X);
		IntStream.range(0, 3 + m).forEach(i -> P[i][i] = 100);
		// Error State intialized as zero
		double[][] x = new double[n][1];
		kfObj.setState_ProcessCov(x, P);
		return iterate(X, SatMap, timeList, useIGS, obsvCodeList);

	}

	TreeMap<Long, double[]> iterate(SimpleMatrix X, TreeMap<Long, ArrayList<Satellite>> SatMap,
			ArrayList<Long> timeList, boolean useIGS, String[] obsvCodeList) throws Exception {
		TreeMap<Long, double[]> estStateMap = new TreeMap<Long, double[]>();
		int m = obsvCodeList.length;
		int x_size = 3 +  m;
		long time = timeList.get(0);
		// Start from 2nd epoch
		for (int i = 1; i < timeList.size(); i++) {
			long currentTime = timeList.get(i);
			ArrayList<Satellite> satList = SatMap.get(currentTime);
			double deltaT = (currentTime - time) / 1e3;
			predictTotalState(X, satList, deltaT, useIGS);
			runFilter(X, deltaT, satList, obsvCodeList);
			double[] estState = new double[x_size];
			IntStream.range(0, x_size).forEach(j -> estState[j] = X.get(j));
			// Add position estimate to the list
			estStateMap.put(currentTime, estState);
			if (!MatrixFeatures_DDRM.isPositiveSemidefinite(kfObj.getCovariance().getMatrix())) {
				throw new Exception("PositiveDefinite test Failed");
			}
			time = currentTime;
		}
		return estStateMap;
	}

	private void predictTotalState(SimpleMatrix X, ArrayList<Satellite> satList, double deltaT, boolean useIGS)
			throws Exception {

		double[] vel = LinearLeastSquare.getEstVel(satList, true, true, false, false,
				new double[] { X.get(0), X.get(1), X.get(2) }, useIGS);
//		for(int i=3;i<vel.length;i++)
//		{
//			vel[i] -=122;
//		}
		for (int i = 0; i < vel.length; i++) {
			X.set(i, X.get(i) + (vel[i] * deltaT));
		}
	}

	private void runFilter(SimpleMatrix X, double deltaT, ArrayList<Satellite> satList, String[] obsvCodeList)
			throws Exception {

		boolean isWeighted = true;
		boolean useAndroidW = false;
		// Satellite count
		int n = satList.size();
		int m = obsvCodeList.length;
		SimpleMatrix Cxx_dot_hat = LinearLeastSquare.getCxx_hat(Measurement.Doppler);
		// Assign Q and F matrix
		kfObj.configDoppler(deltaT, Cxx_dot_hat, m);
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
		double[][] H = new double[n][3+m];
		// Measurement vector
		double[][] z = new double[n][1];
		// Estimated Measurement vector
		double[][] ze = new double[n][1];
		// Measurement Noise
		double[][] R = new double[n][n];

		double priorVarOfUnitW = pseudorange_priorVarOfUnitW;
		if (isWeighted) {
			if (useAndroidW) {
				for (int i = 0; i < n; i++) {
					R[i][i] = Math.pow(satList.get(i).getReceivedSvTimeUncertaintyNanos() * SpeedofLight * 1e-9, 2);
				}
			} else {
				SimpleMatrix Cyy = Weight.getNormCyy(satList, priorVarOfUnitW);
				R = Matrix.matrix2Array(Cyy);
			}
		} else {
			for (int i = 0; i < n; i++) {
				R[i][i] = priorVarOfUnitW;
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
					H[i][3 + j] = 1;
				}
			}

			z[i][0] = PR - approxPR;
			final int _i = i;
			IntStream.range(0, 3).forEach(j -> H[_i][j] = -LOS[j] / approxGR);
		}

		// Perform Update Step
		kfObj.update(z, R, ze, H);
		SimpleMatrix x = kfObj.getState();
		for (int i = 0; i < 3 + m; i++) {
			X.set(i, X.get(i) + x.get(i));
			x.set(i, 0);
		}

	}

}
