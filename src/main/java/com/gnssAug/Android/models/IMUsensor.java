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
	private double tRx;

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

	public double gettRx() {
		return tRx;
	}

	// Accurate upto nanosec range, as bootGPStime is chosen as Long, IMU tRx
	// doesn't need to be as accurate GNSS
	public void settRx(long bootGPStime) {
		tRx = (bootGPStime + elapsedRealtimeNanos) / 1e9;
	}

}
