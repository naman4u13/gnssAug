package com.gnssAug.Rinex.models;

import java.util.Arrays;
import java.util.Calendar;

public class Satellite extends Observable {

	private double[] satEcef;
	private double satClkOff;
	// Note this is GPS System time at time of Transmission
	private double t;
	private double[] satVel;
	// Note this Clock Drift is derived, its not what we get from Ephemeris
	private double satClkDrift;
	private double[] satEci;
	// Note this is GPS System time at time of Reception + Receiver clock offset
	private double tRX;
	// time
	private Calendar time;
	// Elevation and Azimuth Angle
	private double[] elevAzm;
	private double phaseWindUp;
	private double ionoErr;
	private double tropoErr;

	public double[] getSatEcef() {
		return satEcef;
	}

	public void setSatEcef(double[] eCEF) {
		satEcef = eCEF;
	}

	public double[] getSatEci() {
		return satEci;
	}

	public void setSatEci(double[] eCI) {
		satEci = eCI;
	}

	public void updateECI(double rcvrClkOff) {
		compECI(rcvrClkOff);

	}

	public void compECI() {
		compECI(0);
	}

	public void compECI(double time) {
		satEci = new double[3];
		final double OMEGA_E_DOT = 7.2921151467E-5;// WGS-84 value of the Earth's rotation rate

		// eciArg = Earth_Rotation_Rate *(Propgation_Time)
		double eciArg = OMEGA_E_DOT * (tRX - time - t);
		satEci[0] = (satEcef[0] * Math.cos(eciArg)) + (satEcef[1] * Math.sin(eciArg));
		satEci[1] = -(satEcef[0] * Math.sin(eciArg)) + (satEcef[1] * Math.cos(eciArg));
		satEci[2] = satEcef[2];

	}

	public Satellite(Observable satModel, double[] eCEF, double satClkOff, double t, double tRX, double[] satVel,
			double satClkDrift, double[] ECI, Calendar time) {
		super(satModel);
		satEcef = eCEF;
		this.satClkOff = satClkOff;
		this.t = t;
		this.tRX = tRX;
		this.satVel = satVel;
		this.satClkDrift = satClkDrift;
		this.satEci = ECI;
		this.time = time;

	}

	public double getSatClkOff() {
		return satClkOff;
	}

	public void setSatClkOff(double satClkOff) {
		this.satClkOff = satClkOff;
	}

	public double getT() {
		return t;
	}

	public void setT(double t) {
		this.t = t;
	}

	@Override
	public String toString() {
		return super.toString() + "Satellite [satEcef=" + Arrays.toString(satEcef) + ", satClkOff=" + satClkOff + ", t=" + t
				+ ", tRX=" + tRX + ", satVel=" + Arrays.toString(satVel) + ", satClkDrift=" + satClkDrift + ", satEci="
				+ Arrays.toString(satEci) + "]";
	}

	public double[] getSatVel() {
		return satVel;
	}

	public void setSatVel(double[] satVel) {
		this.satVel = satVel;
	}

	public double getSatClkDrift() {
		return satClkDrift;
	}

	public void setSatClkDrift(double satClkDrift) {
		this.satClkDrift = satClkDrift;
	}

	public double gettRX() {
		return tRX;
	}

	public void settRX(double tRX) {
		this.tRX = tRX;
	}

	public Calendar getTime() {
		return time;
	}

	public void setTime(Calendar time) {
		this.time = time;
	}

	public double[] getElevAzm() {
		return elevAzm;
	}

	public void setElevAzm(double[] elevAzm) {
		this.elevAzm = elevAzm;
	}

	public double getPhaseWindUp() {
		return phaseWindUp;
	}

	public void setPhaseWindUp(double phaseWindUp) {
		this.phaseWindUp = phaseWindUp;
	}

	public double getIonoErr() {
		return ionoErr;
	}

	public void setIonoErr(double ionoErr) {
		this.ionoErr = ionoErr;

	}

	public double getTropoErr() {
		return tropoErr;
	}

	public void setTropoErr(double tropoErr) {
		this.tropoErr = tropoErr;

	}

}
