package com.gnssAug.Android.estimation.KalmanFilter.Models;

public class Covariance {

	private double[][] cov;

	public Covariance(double pX, double pY, double pZ, double vX, double vY, double vZ, double attX, double attY,
			double attZ, double biasAccX, double biasAccY, double biasAccZ, double biasGyroX, double biasGyroY,
			double biasGyroZ) {
		super();
		this.cov = new double[15][15];
		cov[0][0] = pX;
		cov[1][1] = pY;
		cov[2][2] = pZ;
		cov[3][3] = vX;
		cov[4][4] = vY;
		cov[5][5] = vZ;
		cov[6][6] = attX;
		cov[7][7] = attY;
		cov[8][8] = attZ;
		cov[9][9] = biasAccX;
		cov[10][10] = biasAccY;
		cov[11][11] = biasAccZ;
		cov[12][12] = biasGyroX;
		cov[13][13] = biasGyroY;
		cov[14][14] = biasGyroZ;
	}

	public double getpX() {
		return cov[0][0];
	}

	public void setpX(double pX) {
		cov[0][0] = pX;
	}

	public double getpY() {
		return cov[1][1];
	}

	public void setpY(double pY) {
		cov[1][1] = pY;
	}

	public double getpZ() {
		return cov[2][2];
	}

	public void setpZ(double pZ) {
		this.cov[2][2] = pZ;
	}

	public double getvX() {
		return cov[3][3];
	}

	public void setvX(double vX) {
		cov[3][3] = vX;
	}

	public double getvY() {
		return cov[4][4];
	}

	public void setvY(double vY) {
		cov[4][4] = vY;
	}

	public double getvZ() {
		return cov[5][5];
	}

	public void setvZ(double vZ) {
		cov[5][5] = vZ;
	}

	public double getAttX() {
		return cov[6][6];
	}

	public void setAttX(double attX) {
		cov[6][6] = attX;
	}

	public double getAttY() {
		return cov[7][7];
	}

	public void setAttY(double attY) {
		cov[7][7] = attY;
	}

	public double getAttZ() {
		return cov[8][8];
	}

	public void setAttZ(double attZ) {
		cov[8][8] = attZ;
	}

	public double getBiasAccX() {
		return cov[9][9];
	}

	public void setBiasAccX(double biasAccX) {
		cov[9][9] = biasAccX;
	}

	public double getBiasAccY() {
		return cov[10][10];
	}

	public void setBiasAccY(double biasAccY) {
		cov[10][10] = biasAccY;
	}

	public double getBiasAccZ() {
		return cov[11][11];
	}

	public void setBiasAccZ(double biasAccZ) {
		cov[11][11] = biasAccZ;
	}

	public double getBiasGyroX() {
		return cov[12][12];
	}

	public void setBiasGyroX(double biasGyroX) {
		cov[12][12] = biasGyroX;
	}

	public double getBiasGyroY() {
		return cov[13][13];
	}

	public void setBiasGyroY(double biasGyroY) {
		cov[13][13] = biasGyroY;
	}

	public double getBiasGyroZ() {
		return cov[14][14];
	}

	public void setBiasGyroZ(double biasGyroZ) {
		cov[14][14] = biasGyroZ;
	}

	public double get(int i, int j) {
		return cov[i][j];
	}

	public void set(int i, int j, double val) {
		cov[i][j] = val;
	}
}
