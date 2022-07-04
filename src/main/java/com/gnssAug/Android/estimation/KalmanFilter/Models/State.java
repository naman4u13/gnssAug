package com.gnssAug.Android.estimation.KalmanFilter.Models;

public class State {

	private double[] p;
	private double[] v;
	private double[][] dcm;
	private double[] accBias;
	private double[] gyroBias;

	public State(double pX, double pY, double pZ, double vX, double vY, double vZ, double[][] dcm, double biasAccX,
			double biasAccY, double biasAccZ, double biasGyroX, double biasGyroY, double biasGyroZ) {
		super();
		this.p[0] = pX;
		this.p[1] = pY;
		this.p[2] = pZ;
		this.v[0] = vX;
		this.v[1] = vY;
		this.v[2] = vZ;
		this.dcm = dcm;
		this.accBias[0] = biasAccX;
		this.accBias[1] = biasAccY;
		this.accBias[2] = biasAccZ;
		this.gyroBias[0] = biasGyroX;
		this.gyroBias[1] = biasGyroY;
		this.gyroBias[2] = biasGyroZ;
	}

	public double[] getP() {
		return p;
	}

	public void setP(double[] p) {
		this.p = p;
	}

	public double[] getV() {
		return v;
	}

	public void setV(double[] v) {
		this.v = v;
	}

	public double[][] getDcm() {
		return dcm;
	}

	public void setDcm(double[][] dcm) {
		this.dcm = dcm;
	}

	public double[] getAccBias() {
		return accBias;
	}

	public void setAccBias(double[] accBias) {
		this.accBias = accBias;
	}

	public double[] getGyroBias() {
		return gyroBias;
	}

	public void setGyroBias(double[] gyroBias) {
		this.gyroBias = gyroBias;
	}

}
