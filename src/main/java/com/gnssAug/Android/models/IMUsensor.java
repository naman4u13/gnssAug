package com.gnssAug.Android.models;

import com.gnssAug.Android.constants.AndroidSensor;

public class IMUsensor {

	private AndroidSensor type;
	private long utcTimeMillis;
	private long elapsedRealtimeNanos;
	private double[] val = new double[3];
	private double[] bias = new double[3];
	private long tRx;

	public IMUsensor(String[] data) {
		super();
		this.type = data[0].equals("UncalAccel") ? AndroidSensor.Accelerometer
				: data[0].equals("UncalGyro") ? AndroidSensor.Gyroscope : AndroidSensor.Magnetometer;
		this.utcTimeMillis = Long.parseLong(data[1]);
		this.elapsedRealtimeNanos = Long.parseLong(data[2]);
		this.val[0] = Double.parseDouble(data[3]);
		this.val[1] = Double.parseDouble(data[4]);
		this.val[2] = Double.parseDouble(data[5]);
		if (data.length > 6) {
			this.bias[0] = data[6].isBlank() ? 0 : Double.parseDouble(data[6]);
			this.bias[1] = data[7].isBlank() ? 0 : Double.parseDouble(data[7]);
			this.bias[2] = data[8].isBlank() ? 0 : Double.parseDouble(data[8]);
		}

	}

	// Constructor used for interpolated IMU sensor values, don't require
	// utcTimeMillis, elapsedRealtimeNanos and biases
	public IMUsensor(AndroidSensor type, double[] val, double[] bias, long tRx, long utcTimeMillis) {
		super();
		this.type = type;
		this.val = val;
		this.bias = bias;
		this.tRx = tRx;
		this.utcTimeMillis = utcTimeMillis;
	}

	public AndroidSensor getType() {
		return type;
	}

	public long getUtcTimeMillis() {
		return utcTimeMillis;
	}

	public long getElapsedRealtimeNanos() {
		return elapsedRealtimeNanos;
	}

	public double[] getVal() {
		return val;
	}

	public double[] getBias() {
		return bias;
	}

	public long gettRx() {
		return tRx;
	}

	// In millisec, used for timestamp, IMU tRx
	// doesn't need to be as accurate GNSS
	public void settRx(long bootGPStime) {
		tRx = (long) ((bootGPStime + elapsedRealtimeNanos) / 1e6);
	}

}
