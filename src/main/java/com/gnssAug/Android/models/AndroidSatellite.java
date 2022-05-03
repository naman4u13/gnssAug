package com.gnssAug.Android.models;

public class AndroidSatellite extends AndroidGNSSLog {

	// Note this is GPS System time at time of Transmission
	private double t;
	// Corrected pseudorange
	private double pseudorange;
	private double[] satEcef;
	private double[] satVel;
	private double[] satEci;
	private double[] elevAzm;

	public AndroidSatellite(AndroidGNSSLog log, double t, double pseudorange, double[] satEcef, double[] satVel) {
		super(log);
		this.t = t;
		this.pseudorange = pseudorange;
		this.satEcef = satEcef;
		this.satVel = satVel;
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

}
