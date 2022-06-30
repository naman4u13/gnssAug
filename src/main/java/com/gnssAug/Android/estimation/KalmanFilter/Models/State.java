package com.gnssAug.Android.estimation.KalmanFilter.Models;

public class State {

	private double[] p;
	private double[] v;
	private double[][] dcm;
	private double[] biasAcc;
	private double[] biasGyro;

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
		this.biasAcc[0] = biasAccX;
		this.biasAcc[1] = biasAccY;
		this.biasAcc[2] = biasAccZ;
		this.biasGyro[0] = biasGyroX;
		this.biasGyro[1] = biasGyroY;
		this.biasGyro[2] = biasGyroZ;
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

	public double[] getBiasAcc() {
		return biasAcc;
	}

	public void setBiasAcc(double[] biasAcc) {
		this.biasAcc = biasAcc;
	}

	public double[] getBiasGyro() {
		return biasGyro;
	}

	public void setBiasGyro(double[] biasGyro) {
		this.biasGyro = biasGyro;
	}

}
