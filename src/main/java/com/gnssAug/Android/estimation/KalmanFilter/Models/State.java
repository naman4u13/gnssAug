package com.gnssAug.Android.estimation.KalmanFilter.Models;

import org.ejml.simple.SimpleMatrix;

public class State {

	private double[] p = new double[3];
	private double[] v = new double[3];
	private SimpleMatrix dcm;
	private double[] accBias = new double[3];
	private double[] gyroBias = new double[3];
	private double[] rxClk = new double[2];

	public State(double pX, double pY, double pZ, double vX, double vY, double vZ, SimpleMatrix dcm, double biasAccX,
			double biasAccY, double biasAccZ, double biasGyroX, double biasGyroY, double biasGyroZ, double clkOff,
			double clkDrift) {
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
		this.rxClk[0] = clkOff;
		this.rxClk[1] = clkDrift;
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

	public SimpleMatrix getDcm() {
		return dcm;
	}

	public void setDcm(SimpleMatrix dcm) {
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

	public double[] getRxClk() {
		return rxClk;
	}

	public void setRxClk(double[] rxClk) {
		this.rxClk = rxClk;
	}

}
