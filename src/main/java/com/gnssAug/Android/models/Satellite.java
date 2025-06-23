package com.gnssAug.Android.models;

public class Satellite extends GNSSLog implements Cloneable {

	// Note this is GPS System time at time of Transmission
	private double t;
	// Corrected pseudorange
	private double pseudorange;
	// Corrected phase
	private double phase;
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
	// Iono param required for PPP & CSDR
	private double ionoErr;
	// Tropo param required to be removed from Doppler Measurements
	private double tropoErr;
	// Phase lock indicator
	private boolean isPhaseLocked = false;
	// Experimental param
	private double clkRate;
	// Adaptive Var
	private double prVar;
	// Tropo Wet Mapping Function
	private double wetMF;
	
	public Satellite(GNSSLog log, double t, double pseudorange, double[] satEcef, double[] satVel, double rangeRate, double phase) {
		this(log,t,pseudorange,satEcef,satVel,rangeRate);
		this.phase = phase;

	}
	public Satellite(GNSSLog log, double t, double pseudorange, double[] satEcef, double[] satVel, double rangeRate) {
		super(log);
		this.t = t;
		this.pseudorange = pseudorange;
		this.satEcef = satEcef;
		this.satVel = satVel;
		this.rangeRate = rangeRate;
		compECI();

	} 
	
	@Override
	public Satellite clone() throws CloneNotSupportedException {
	    Satellite cloned = (Satellite) super.clone();
	    cloned.satEcef = this.satEcef.clone();
	    cloned.satVel = this.satVel.clone();
	    cloned.satEci = this.satEci.clone();
	    cloned.elevAzm = this.elevAzm.clone();
	    return cloned;
	}
	private void compECI() {
		satEci = new double[3];
		final double OMEGA_E_DOT = 7.2921151467E-5;// WGS-84 value of the Earth's rotation rate

		// eciArg = Earth_Rotation_Rate *(Propgation_Time)
		double eciArg = OMEGA_E_DOT * (gettRx() - t);
		satEci[0] = (satEcef[0] * Math.cos(eciArg)) + (satEcef[1] * Math.sin(eciArg));
		satEci[1] = -(satEcef[0] * Math.sin(eciArg)) + (satEcef[1] * Math.cos(eciArg));
		satEci[2] = satEcef[2];
		
//		double[] satVel_new = new double[3];
//		satVel_new[0] = (satVel[0] * Math.cos(eciArg)) + (satVel[1] * Math.sin(eciArg));
//		satVel_new[1] = -(satVel[0] * Math.sin(eciArg)) + (satVel[1] * Math.cos(eciArg));
//		satVel_new[2] = satVel[2];
//		satVel = new double[] {satVel_new[0],satVel_new[1],satVel_new[2]};
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

	public double getPhase() {
		return phase;
	}
	public void setPhase(double phase) {
		this.phase = phase;
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
	public double getIonoErr() {
		return ionoErr;
	}
	public void setIonoErr(double ionoErr) {
		this.ionoErr = ionoErr;
	}
	public boolean isPhaseLocked() {
		return isPhaseLocked;
	}
	public void setPhaseLocked(boolean isPhaseLocked) {
		this.isPhaseLocked = isPhaseLocked;
	}
	public double getTropoErr() {
		return tropoErr;
	}
	public void setTropoErr(double tropoErr) {
		this.tropoErr = tropoErr;
	}
	public void setClkRate(double clkRate)
	{
		this.clkRate = clkRate;
	}
	public double getClkRate()
	{
		return clkRate;
	}
	public double getPrVar() {
		return prVar;
	}
	public void setPrVar(double prVar) {
		this.prVar = prVar;
	}
	public double getWetMF() {
		return wetMF;
	}
	public void setWetMF(double wetMF) {
		this.wetMF = wetMF;
	}
	
}
