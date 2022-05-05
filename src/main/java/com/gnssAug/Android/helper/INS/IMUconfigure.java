package com.gnssAug.Android.helper.INS;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeMap;

import com.gnssAug.Android.constants.AndroidSensor;
import com.gnssAug.Android.models.IMUsensor;
import com.gnssAug.Android.utility.Interpolator;

public class IMUconfigure {

	private static TreeMap<Long, HashMap<AndroidSensor, IMUsensor>> imuMap;

	public static TreeMap<Long, HashMap<AndroidSensor, IMUsensor>> configure(long t0, int sampleRate,
			ArrayList<IMUsensor> imuList) throws Exception {
		ArrayList<IMUsensor> acc = new ArrayList<IMUsensor>();
		ArrayList<IMUsensor> gyro = new ArrayList<IMUsensor>();
		ArrayList<IMUsensor> mag = new ArrayList<IMUsensor>();
		imuMap = new TreeMap<Long, HashMap<AndroidSensor, IMUsensor>>();
		int n = imuList.size();
		for (int i = 0; i < n; i++) {
			IMUsensor imu = imuList.get(i);
			if (imu.getType() == AndroidSensor.Accelerometer) {
				acc.add(imu);
			} else if (imu.getType() == AndroidSensor.Gyroscope) {
				gyro.add(imu);
			} else if (imu.getType() == AndroidSensor.Magnetometer) {
				mag.add(imu);
			}
		}
		if (sampleRate % 1000 == 0) {
			System.err.println("Erroneous imu sampling rate");
			throw new Exception("Erroneous imu sampling rate");
		}
		int diff = 1000 / sampleRate;
		downsample(t0, diff, acc);
		downsample(t0, diff, gyro);
		downsample(t0, diff, mag);
		imuMap.remove(imuMap.lastKey());
		imuMap.remove(imuMap.firstKey());
		return imuMap;
	}

	private static void downsample(long t, int diff, ArrayList<IMUsensor> imuList) throws Exception {

		long tRx0 = imuList.get(0).gettRx();
		if (t < tRx0) {
			t += Math.ceil((tRx0 - t) / (double) diff) * diff;
		} else if (tRx0 < t) {
			t -= Math.floor((t - tRx0) / (double) diff) * diff;
		}

		if (tRx0 > t) {
			System.err.println("ERROR in downsampling IMU");
			throw new Exception("ERROR in downsampling IMU");
		}
		int n = imuList.size();
		for (int i = 0; i < n; i++) {
			IMUsensor imu = imuList.get(i);
			long tRx = imu.gettRx();
			if (t <= tRx) {
				if (t < tRx) {
					IMUsensor imu_prev = imuList.get(i - 1);
					IMUsensor _imu = interpolateIMU(imu_prev.gettRx(), tRx, t, imu_prev, imu);
					imuMap.computeIfAbsent(t, k -> new HashMap<AndroidSensor, IMUsensor>()).put(_imu.getType(), _imu);

				} else {
					imuMap.computeIfAbsent(t, k -> new HashMap<AndroidSensor, IMUsensor>()).put(imu.getType(), imu);

				}
				t += diff;
				i--;
			} else if (tRx < t) {
				continue;
			}
		}

	}

	private static IMUsensor interpolateIMU(long x0, long x1, long x, IMUsensor y0, IMUsensor y1) {
		long[] X = new long[] { x0, x1 };
		double xVal = Interpolator.linear(X, new double[] { y0.getX(), y1.getX() }, x);
		double yVal = Interpolator.linear(X, new double[] { y0.getY(), y1.getY() }, x);
		double zVal = Interpolator.linear(X, new double[] { y0.getZ(), y1.getZ() }, x);
		double utcTimeMillis = Interpolator.linear(X, new double[] { y0.getUtcTimeMillis(), y1.getUtcTimeMillis() }, x);
		IMUsensor y = new IMUsensor(y0.getType(), xVal, yVal, zVal, x, (long) utcTimeMillis);
		return y;
	}

}
