package com.gnssAug.Rinex.models;

public class SatResidual {
	// GPS Time
	private double t;
	// Elevation Angle(Radians)
	private double elevAngle;
	// Residual(meter)
	private double residual;
	// Measurement Noise Std Dev(meter)
	private double noiseStdDev;
	// derived from Outlier detection
	private boolean isOutlier;
	
	public SatResidual(double t, double elevAngle, double residual, boolean isOutlier) {
		super();
		this.t = t;
		this.elevAngle = elevAngle;
		this.residual = residual;
		this.isOutlier = isOutlier;
	}

	public SatResidual(double t, double elevAngle, double residual, boolean isOutlier,double noiseStdDev) {
		super();
		this.t = t;
		this.elevAngle = elevAngle;
		this.residual = residual;
		this.isOutlier = isOutlier;
		this.noiseStdDev = noiseStdDev;
	}

	public double getT() {
		return t;
	}

	public double getElevAngle() {
		return elevAngle;
	}

	public double getResidual() {
		return residual;
	}

	public double getNoiseStdDev() {
		return noiseStdDev;
	}

	public boolean isOutlier() {
		return isOutlier;
	}
}
