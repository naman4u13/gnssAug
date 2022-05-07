package com.gnssAug.Android.helper.INS;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.TreeMap;

import com.gnssAug.Android.constants.AndroidSensor;
import com.gnssAug.Android.estimation.LinearLeastSquare;
import com.gnssAug.Android.helper.GeomagneticField;
import com.gnssAug.Android.helper.RotationMatrix;
import com.gnssAug.Android.models.IMUsensor;
import com.gnssAug.Android.models.Satellite;
import com.gnssAug.Android.utility.LatLonUtil;

public class StateInitialization {

	public static void initialize(TreeMap<Long, HashMap<AndroidSensor, IMUsensor>> imuMap,
			ArrayList<ArrayList<Satellite>> SVlist) throws Exception {

		double[] acc0 = new double[3];
		double[] mag0 = new double[3];
		double[] ecef0 = new double[3];
		Iterator<Entry<Long, HashMap<AndroidSensor, IMUsensor>>> iterator = imuMap.entrySet().iterator();
		// time period at the beginning of the data, for which acc and mag data will be
		// averaged for. It is assumed that IMU is stationary during this period
		int t = 15;
		for (int i = 1; i <= t; i++) {
			HashMap<AndroidSensor, IMUsensor> imu = iterator.next().getValue();
			double[] acc = imu.get(AndroidSensor.Accelerometer).getVal();
			double[] accBias = imu.get(AndroidSensor.Accelerometer).getBias();
			double[] mag = imu.get(AndroidSensor.Magnetometer).getVal();
			double[] magBias = imu.get(AndroidSensor.Magnetometer).getBias();
			ArrayList<Satellite> satList = SVlist.get(i - 1);
			double[] ecef = LinearLeastSquare.process(satList, true);
			for (int j = 0; j < 3; j++) {
				acc0[j] = (acc0[j] * (i - 1) / i) + ((acc[j] - accBias[j]) / i);
				mag0[j] = (mag0[j] * (i - 1) / i) + ((mag[j] - magBias[j]) / i);
				ecef0[j] = (ecef0[j] * (i - 1) / i) + (ecef[j] / i);
			}

		}
		double[] llh0 = LatLonUtil.ecef2lla(ecef0);
		long utcTimeMilli = SVlist.get(0).get(0).getUtcTimeMillis();
		double g = norm2(acc0);
		double m = norm2(mag0);
		// Rotation from body frame to ENU (where N is magnetic north)
		double[][] dcm = getRotationMatrix(acc0, mag0);
		// Require Rotation matrix from body frame to ENU (where N is geographic North)
		// Calculate declination angle at initial position and time
		GeomagneticField gmf = new GeomagneticField((float) llh0[0], (float) llh0[1], (float) llh0[2], utcTimeMilli);
		double declinationAngle = Math.toRadians(gmf.getDeclination());
		double[] euler = RotationMatrix.dcm2euler(dcm);
		// Compensate for declination angle to obtain euler angles for geographic north
		// frame
		euler[0] = euler[0] + declinationAngle;
		dcm = RotationMatrix.euler2dcm(euler);

		System.out.print("");
	}

	private static double norm2(double[] vector) {
		return Math.sqrt(Arrays.stream(vector).map(i -> i * i).sum());
	}

	public static double[][] getRotationMatrix(double[] gravity, double[] magnetic) throws Exception {

		double[][] R = new double[3][3];

		// TODO: move this to native code for efficiency
		double Ax = gravity[0];
		double Ay = gravity[1];
		double Az = gravity[2];
		final double normsqA = (Ax * Ax + Ay * Ay + Az * Az);
		final double g = 9.81f;
		final double freeFallGravitySquared = 0.01f * g * g;
		if (normsqA < freeFallGravitySquared) {
			throw new Exception("Error in IMU stateIntialization: gravity less than 10% of normal value");
		}
		final double Ex = magnetic[0];
		final double Ey = magnetic[1];
		final double Ez = magnetic[2];
		double Hx = Ey * Az - Ez * Ay;
		double Hy = Ez * Ax - Ex * Az;
		double Hz = Ex * Ay - Ey * Ax;
		final double normH = (float) Math.sqrt(Hx * Hx + Hy * Hy + Hz * Hz);
		if (normH < 0.1f) {
			throw new Exception(
					"Error in IMU stateIntialization: device is close to free fall (or in space?), or close to magnetic north pole. Typical values are > 100.");
		}
		final double invH = 1.0f / normH;
		Hx *= invH;
		Hy *= invH;
		Hz *= invH;
		final double invA = 1.0f / (float) Math.sqrt(Ax * Ax + Ay * Ay + Az * Az);
		Ax *= invA;
		Ay *= invA;
		Az *= invA;
		final double Mx = Ay * Hz - Az * Hy;
		final double My = Az * Hx - Ax * Hz;
		final double Mz = Ax * Hy - Ay * Hx;

		R[0][0] = Hx;
		R[0][1] = Hy;
		R[0][2] = Hz;
		R[1][0] = Mx;
		R[1][1] = My;
		R[1][2] = Mz;
		R[2][0] = Ax;
		R[2][1] = Ay;
		R[2][2] = Az;

		return R;
	}
}
