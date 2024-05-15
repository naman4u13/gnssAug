package com.gnssAug.Android.models;

import org.ejml.simple.SimpleMatrix;

public class CycleSlipDetect {

	private Satellite sat;
	private double dopplerDR;
	private double carrierPhaseDR;
	private boolean isCS = false;
	private boolean isRepaired = false;
	private double ionoRate;
	private double wavelength;
	private double satVelCorr;
	private SimpleMatrix unitLOS;
	private double time;
	private double clkDrift;
	private double trueDR;
	private double floatAmb;
	private double intAmb;
	public CycleSlipDetect(Satellite sat, double dopplerDR, double carrierPhaseDR, double ionoRate, boolean isCS,
			double wavelength, double satVelCorr, SimpleMatrix unitLOS, double time,double trueDR) {
		this( sat,  dopplerDR,  carrierPhaseDR,  ionoRate,  isCS,
				 wavelength,  satVelCorr,  unitLOS,  time);
		this.trueDR = trueDR;
		
	}
	public CycleSlipDetect(Satellite sat, double dopplerDR, double carrierPhaseDR, double ionoRate, boolean isCS,
			double wavelength, double satVelCorr, SimpleMatrix unitLOS, double time) {
		super();
		this.sat = sat;
		this.dopplerDR = dopplerDR;
		this.carrierPhaseDR = carrierPhaseDR;
		this.ionoRate = ionoRate;
		this.isCS = isCS;
		this.wavelength = wavelength;
		this.satVelCorr = satVelCorr;
		this.unitLOS = unitLOS;
		this.time = time;
	}

	public Satellite getSat() {
		return sat;
	}

	public double getDopplerDR() {
		return dopplerDR;
	}

	public double getCarrierPhaseDR() {
		return carrierPhaseDR;
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

	public double getWavelength() {
		return wavelength;
	}

	public double getSatVelCorr() {
		return satVelCorr;
	}

	public SimpleMatrix getUnitLOS() {
		return unitLOS;
	}

	public boolean isRepaired() {
		return isRepaired;
	}

	public void setRepaired(boolean isRepaired) {
		this.isRepaired = isRepaired;
	}

	public double getTime() {
		return time;
	}

	public double getClkDrift() {
		return clkDrift;
	}

	public void setClkDrift(double clkDrift) {
		this.clkDrift = clkDrift;
	}
	public double getTrueDR() {
		return trueDR;
	}
	public double getFloatAmb() {
		return floatAmb;
	}
	public void setFloatAmb(double floatAmb) {
		this.floatAmb = floatAmb;
	}
	public double getIntAmb() {
		return intAmb;
	}
	public void setIntAmb(double intAmb) {
		this.intAmb = intAmb;
	}

}
