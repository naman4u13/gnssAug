package com.gnssAug.Android.models;

public class CycleSlipDetect {
	
	private Satellite sat;
	private double dopplerDR;
	private double carrierPhaseDR;
	private boolean isCS;
	private double ionoRate;
	private double wavelength;
	private double satVelCorr;
	private int index;
	private double approxCS;
	public CycleSlipDetect(Satellite sat, double dopplerDR, double carrierPhaseDR, double ionoRate,boolean isCS,double wavelength,double satVelCorr,int index,double approxCS) {
		super();
		this.sat = sat;
		this.dopplerDR = dopplerDR;
		this.carrierPhaseDR = carrierPhaseDR;
		this.ionoRate = ionoRate;
		this.isCS = isCS;
		this.wavelength = wavelength;
		this.satVelCorr = satVelCorr;
		this.index = index;
		this.approxCS = approxCS;
	}
	public Satellite getSat() {
		return sat;
	}
	public void setSat(Satellite sat) {
		this.sat = sat;
	}
	public double getDopplerDR() {
		return dopplerDR;
	}
	public void setDopplerDR(double dopplerDR) {
		this.dopplerDR = dopplerDR;
	}
	public double getCarrierPhaseDR() {
		return carrierPhaseDR;
	}
	public void setCarrierPhaseDR(double carrierPhaseDR) {
		this.carrierPhaseDR = carrierPhaseDR;
	}
	public boolean isCS() {
		return isCS;
	}
	public void setCS(boolean isCS) {
		this.isCS = isCS;
	}
	public double getIonoRate() {
		return ionoRate;
	}
	public void setIonoRate(double ionoRate) {
		this.ionoRate = ionoRate;
	}
	public double getWavelength() {
		return wavelength;
	}
	public void setWavelength(double wavelength) {
		this.wavelength = wavelength;
	}
	public double getSatVelCorr() {
		return satVelCorr;
	}
	public void setSatVelCorr(double satVelCorr) {
		this.satVelCorr = satVelCorr;
	}
	public int getIndex() {
		return index;
	}
	public void setIndex(int index) {
		this.index = index;
	}
	public double getApproxCS() {
		return approxCS;
	}
	public void setApproxCS(double approxCS) {
		this.approxCS = approxCS;
	}
	

}
