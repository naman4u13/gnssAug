package com.gnssAug.Android.models;

import com.gnssAug.Android.constants.AndroidSensor;

public class IMUsensor {

	private AndroidSensor type;
	private long utcTimeMillis;
	private long elapsedRealtimeNanos;
	private double xVal;
	private double yVal;
	private double zVal;
	private double biasX;
	private double biasY;
	private double biasZ;
	private long tRx;

	public IMUsensor(String[] data) {
		super();
		this.type = data[0].equals("UncalAccel") ? AndroidSensor.Accelerometer
				: data[0].equals("UncalGyro") ? AndroidSensor.Gyroscope : AndroidSensor.Magnetometer;
		this.utcTimeMillis = Long.parseLong(data[1]);
		this.elapsedRealtimeNanos = Long.parseLong(data[2]);
		this.xVal = Double.parseDouble(data[3]);
		this.yVal = Double.parseDouble(data[4]);
		this.zVal = Double.parseDouble(data[5]);
		if (data.length > 6) {
			this.biasX = data[6].isBlank() ? 0 : Double.parseDouble(data[6]);
			this.biasY = data[7].isBlank() ? 0 : Double.parseDouble(data[7]);
			this.biasZ = data[8].isBlank() ? 0 : Double.parseDouble(data[8]);
		}
	}

	// Constructor used for interpolated IMU sensor values, don't require
	// utcTimeMillis, elapsedRealtimeNanos and biases
	public IMUsensor(AndroidSensor type, double xVal, double yVal, double zVal, long tRx, long utcTimeMillis) {
		super();
		this.type = type;
		this.xVal = xVal;
		this.yVal = yVal;
		this.zVal = zVal;
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

	public double getX() {
		return xVal;
	}

	public double getY() {
		return yVal;
	}

	public double getZ() {
		return zVal;
	}

	public double getBiasX() {
		return biasX;
	}

	public double getBiasY() {
		return biasY;
	}

	public double getBiasZ() {
		return biasZ;
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
