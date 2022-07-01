package com.gnssAug.Android.estimation.KalmanFilter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.stream.IntStream;

import org.ejml.simple.SimpleMatrix;

import com.gnssAug.Android.constants.AndroidSensor;
import com.gnssAug.Android.constants.ImuDataSheets;
import com.gnssAug.Android.estimation.LinearLeastSquare;
import com.gnssAug.Android.estimation.KalmanFilter.Models.Covariance;
import com.gnssAug.Android.estimation.KalmanFilter.Models.State;
import com.gnssAug.Android.helper.Rotation;
import com.gnssAug.Android.models.IMUsensor;
import com.gnssAug.Android.models.Satellite;
import com.gnssAug.Android.utility.LatLonUtil;
import com.gnssAug.Android.utility.Matrix;

public class INSfusion {

	public static void process(TreeMap<Long, HashMap<AndroidSensor, IMUsensor>> imuMap,
			ArrayList<ArrayList<Satellite>> SVlist, ArrayList<Long> timeList, double[][] dcm) throws Exception {

		// Convert DCM from Body Frame to ENU to Body frame to NED
		dcm = LatLonUtil.enu_ned_convert(dcm);
		double[] ecef0 = LinearLeastSquare.process(SVlist.get(0), true);
		double[] llh0 = LatLonUtil.ecef2lla(ecef0);
		IntStream.range(0, 2).forEach(i -> llh0[i] = Math.toRadians(llh0[i]));
		// Velocity in ENU frame, zero initially
		double[] vel0 = new double[3];
		State x = new State(llh0[0], llh0[1], llh0[2], vel0[0], vel0[1], vel0[2], dcm, 0, 0, 0, 0, 0, 0);
		// Attitude std deviation is 20 degree, values mentioned below in covariance
		// matrix is in radians
		double accBiasCov = Math.pow(ImuDataSheets.Pixel4.accTurnOnBias, 2);
		double gyroBiasCov = Math.pow(ImuDataSheets.Pixel4.gyroTurnOnBias, 2);
		Covariance p = new Covariance(100, 100, 100, 0.1, 0.1, 0.1, 0.1156, 0.1156, 0.1156, accBiasCov, accBiasCov,
				accBiasCov, gyroBiasCov, gyroBiasCov, gyroBiasCov);
		int n = SVlist.size();
		if (n != timeList.size()) {
			System.err.println("FATAL ERROR in INSfusion.java");
			throw new Exception("size of SVlist and timelist does not match");
		}
		int i = 0;
		do {
			Iterator<Entry<Long, HashMap<AndroidSensor, IMUsensor>>> iterator = imuMap.entrySet().iterator();
			while (iterator.hasNext()) {
				long time = iterator.next().getKey();
				if (time == timeList.get(i)) {

				} else {
					if (time > timeList.get(i)) {
						i++;
					}
					HashMap<AndroidSensor, IMUsensor> imuSensor = iterator.next().getValue();

				}

			}
		} while (i < n);

	}

	private static void updateTotalState(State x, HashMap<AndroidSensor, IMUsensor> imuSensor, double tau) {
		double[] obsAcc = LatLonUtil.enu_ned_convert(imuSensor.get(AndroidSensor.Accelerometer).getVal());
		double[] estAcc = IntStream.range(0, 3).mapToDouble(j -> obsAcc[j] - x.getBiasAcc()[j]).toArray();
		double[] obsGyro = imuSensor.get(AndroidSensor.Gyroscope).getVal();
		double[] estGyro = IntStream.range(0, 3).mapToDouble(j -> obsGyro[j] - x.getBiasGyro()[j]).toArray();

		double lat = x.getP()[0];
		double alt = x.getP()[2];
		double[] vel = x.getV();
		double Re = LatLonUtil.getNormalEarthRadius(lat);
		double earthAngularRate = LatLonUtil.omega_ie;
		// Attitude update
		SimpleMatrix oldDcm = new SimpleMatrix(x.getDcm());
		SimpleMatrix omega_b_ib = new SimpleMatrix(Matrix.getSkewSymMat(estGyro));
		double[] _omega_n_ie = new double[] { earthAngularRate * Math.cos(lat), 0, -earthAngularRate * Math.sin(lat) };
		SimpleMatrix omega_n_ie = new SimpleMatrix(Matrix.getSkewSymMat(_omega_n_ie));
		double[] _omega_n_en = new double[] { vel[0], -vel[1], -vel[0] * Math.tan(lat) };
		Arrays.stream(_omega_n_en).forEach(i -> i = i / (Re + alt));
		SimpleMatrix omega_n_en = new SimpleMatrix(Matrix.getSkewSymMat(_omega_n_en));
		SimpleMatrix I = SimpleMatrix.identity(3);
		SimpleMatrix newDcm = (oldDcm.mult(I.plus(omega_b_ib.scale(tau))))
				.minus((omega_n_ie.plus(omega_n_en)).mult(oldDcm.scale(tau)));
		newDcm = Rotation.reorthonormDcm(newDcm);

		// Specific-Force Frame Transformation
		SimpleMatrix f_b_ib = new SimpleMatrix(3, 1, true, estAcc);
		SimpleMatrix f_n_ib = f_b_ib.scale(0.5).mult((oldDcm.plus(newDcm)));

		// Velocity Update
		SimpleMatrix oldVel = new SimpleMatrix(3, 1, true, vel);
		SimpleMatrix g_n_b = new SimpleMatrix(3, 1, true, new double[] { 0, 0, -LatLonUtil.getGravity(lat, alt) });
		SimpleMatrix newVel = oldVel
				.plus((f_n_ib.plus(g_n_b).minus(oldVel.mult(omega_n_en.plus((omega_n_ie.scale(2)))))).scale(tau));

	}

}
