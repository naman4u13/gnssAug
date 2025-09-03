package com.gnssAug.Android.models;

import org.ejml.simple.SimpleMatrix;

public class CycleSlipDetect {

	private Satellite sat;
	private com.gnssAug.Rinex.models.Satellite igs_sat;
	private double dopplerDR;
	private double carrierPhaseDR;
	private double prDR;
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
	private double floatAmbCov;
	private double intAmbCov;
	private double successRate;
	private double failureRate;
	private boolean isWLcomb = false;
	private boolean exclude = false;
	public CycleSlipDetect(CycleSlipDetect csdObj) {
		this( csdObj.getSat(),  csdObj.getDopplerDR(), csdObj.getCarrierPhaseDR(),  csdObj.getIonoRate(),  csdObj.isCS(),
				csdObj.getWavelength(),  csdObj.getSatVelCorr(),  csdObj.getUnitLOS(),  csdObj.getTime(),csdObj.getTrueDR());
		
		
	}
	
	public CycleSlipDetect(Satellite sat, double dopplerDR, double carrierPhaseDR, double ionoRate, boolean isCS,
			double wavelength, double satVelCorr, SimpleMatrix unitLOS, double time,double trueDR) {
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
		this.trueDR = trueDR;
		
	}
	
	public CycleSlipDetect(Satellite sat, double dopplerDR, double carrierPhaseDR, double ionoRate, boolean isCS,
			double wavelength, double satVelCorr, SimpleMatrix unitLOS, double time,double trueDR,double prDR) {
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
		this.trueDR = trueDR;
		this.prDR = prDR;
		
	}
	public CycleSlipDetect(com.gnssAug.Rinex.models.Satellite igs_sat, double dopplerDR, double carrierPhaseDR, double ionoRate, boolean isCS,
			double wavelength, double satVelCorr, SimpleMatrix unitLOS, double time,double trueDR,double prDR) {
		super();
		this.igs_sat = igs_sat;
		this.dopplerDR = dopplerDR;
		this.carrierPhaseDR = carrierPhaseDR;
		this.ionoRate = ionoRate;
		this.isCS = isCS;
		this.wavelength = wavelength;
		this.satVelCorr = satVelCorr;
		this.unitLOS = unitLOS;
		this.time = time;
		this.trueDR = trueDR;
		this.prDR = prDR;
		
	}
	public CycleSlipDetect(Satellite sat, double dopplerDR, double carrierPhaseDR, double ionoRate, boolean isCS,
			double wavelength, double satVelCorr, SimpleMatrix unitLOS, double time) {
		
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

	public double getFloatAmbCov() {
		return floatAmbCov;
	}

	public void setFloatAmbCov(double floatAmbCov) {
		this.floatAmbCov = floatAmbCov;
	}

	public double getIntAmbCov() {
		return intAmbCov;
	}

	public void setIntAmbCov(double intAmbCov) {
		this.intAmbCov = intAmbCov;
	}

	public double getSuccessRate() {
		return successRate;
	}

	public void setSuccessRate(double successRate) {
		this.successRate = successRate;
	}

	public double getFailureRate() {
		return failureRate;
	}

	public void setFailureRate(double failureRate) {
		this.failureRate = failureRate;
	}

	public boolean isWLcomb() {
		return isWLcomb;
	}

	public void setWLcomb(boolean isWLcomb) {
		this.isWLcomb = isWLcomb;
	}

	public double getPrDR() {
		return prDR;
	}

	public void setPrDR(double prDR) {
		this.prDR = prDR;
	}

	public com.gnssAug.Rinex.models.Satellite getIgs_sat() {
		return igs_sat;
	}

	public boolean isExclude() {
		return exclude;
	}

	public void setExclude(boolean exclude) {
		this.exclude = exclude;
	}
	
	
	
}
