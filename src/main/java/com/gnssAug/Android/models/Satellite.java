package com.gnssAug.Android.models;

public class Satellite extends GNSSLog {

	// Note this is GPS System time at time of Transmission
	private double t;
	// Corrected pseudorange
	private double pseudorange;
	private double[] satEcef;
	private double[] satVel;
	private double[] satEci;
	private double[] elevAzm;
	// Corrected range rate
	private double rangeRate;
	// Experimental param, used for outlier testing
	private boolean isOutlier;
	// Experimental param, true range of satellite
	private double trueRange;
	

	public Satellite(GNSSLog log, double t, double pseudorange, double[] satEcef, double[] satVel, double rangeRate) {
		super(log);
		this.t = t;
		this.pseudorange = pseudorange;
		this.satEcef = satEcef;
		this.satVel = satVel;
		this.rangeRate = rangeRate;
		compECI();

	}

	private void compECI() {
		satEci = new double[3];
		final double OMEGA_E_DOT = 7.2921151467E-5;// WGS-84 value of the Earth's rotation rate

		// eciArg = Earth_Rotation_Rate *(Propgation_Time)
		double eciArg = OMEGA_E_DOT * (gettRx() - t);
		satEci[0] = (satEcef[0] * Math.cos(eciArg)) + (satEcef[1] * Math.sin(eciArg));
		satEci[1] = -(satEcef[0] * Math.sin(eciArg)) + (satEcef[1] * Math.cos(eciArg));
		satEci[2] = satEcef[2];

	}

	public double getT() {
		return t;
	}

	public double getPseudorange() {
		return pseudorange;
	}

	public double[] getSatEcef() {
		return satEcef;
	}

	public double[] getSatVel() {
		return satVel;
	}

	public double[] getSatEci() {
		return satEci;
	}

	public double[] getElevAzm() {
		return elevAzm;
	}

	public void setElevAzm(double[] elevAzm) {
		this.elevAzm = elevAzm;
	}

	public double getRangeRate() {
		return rangeRate;
	}

	public void setPseudorange(double pseudorange) {
		this.pseudorange = pseudorange;
	}

	public boolean isOutlier() {
		return isOutlier;
	}

	public void setOutlier(boolean isOutlier) {
		this.isOutlier = isOutlier;
	}
	
	public double getTrueRange() {
		return trueRange;
	}

	public void setTrueRange(double trueRange) {
		this.trueRange = trueRange;
	}

}
