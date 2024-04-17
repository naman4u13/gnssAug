package com.gnssAug.Android.models;

import org.ejml.simple.SimpleMatrix;

public class TDCP {

	private Satellite sat;
	private double deltaRange;
	private double satVelCorr;
	private SimpleMatrix unitLOS;
	private boolean isOutlier;
	
	public TDCP(Satellite sat, double deltaRange,double satVelCorr,SimpleMatrix unitLOS) {
		super();
		this.sat = sat;
		this.deltaRange = deltaRange;
		this.satVelCorr = satVelCorr;
		this.unitLOS = unitLOS;
	}
	public Satellite getSat() {
		return sat;
	}
	public void setSat(Satellite sat) {
		this.sat = sat;
	}
	public double getDeltaRange() {
		return deltaRange;
	}
	public void setDeltaRange(double deltaRange) {
		this.deltaRange = deltaRange;
	}
	public double getSatVelCorr() {
		return satVelCorr;
	}
	public void setSatVelCorr(double satVelCorr) {
		this.satVelCorr = satVelCorr;
	}
	public SimpleMatrix getUnitLOS() {
		return unitLOS;
	}
	public void setUnitLOS(SimpleMatrix unitLOS) {
		this.unitLOS = unitLOS;
	}
	public boolean isOutlier() {
		return isOutlier;
	}
	public void setOutlier(boolean isOutlier) {
		this.isOutlier = isOutlier;
	}
	
}
