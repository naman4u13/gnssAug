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
	private double CN0;
	
	public SatResidual(double t, double elevAngle, double residual, boolean isOutlier,double CN0) {
		super();
		this.t = t;
		this.elevAngle = elevAngle;
		this.residual = residual;
		this.isOutlier = isOutlier;
		this.CN0 = CN0;
	}

	public SatResidual(double t, double elevAngle, double residual, boolean isOutlier,double noiseStdDev,double CN0) {
		super();
		this.t = t;
		this.elevAngle = elevAngle;
		this.residual = residual;
		this.isOutlier = isOutlier;
		this.noiseStdDev = noiseStdDev;
		this.CN0 = CN0;
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

	public double getCN0() {
		return CN0;
	}
	
}
