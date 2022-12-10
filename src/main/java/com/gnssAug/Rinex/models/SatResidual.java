package com.gnssAug.Rinex.models;

public class SatResidual {
	// GPS Time
	private double t;
	// Elevation Angle(Radians)
	private double elevAngle;
	// Residual(meter)
	private double residual;
	// Measurement Noise Std Dev(meter)
	private double noise;
	// derived from Outlier detection
	private boolean isOutlier;

	public SatResidual(double t, double elevAngle, double residual) {
		super();
		this.t = t;
		this.elevAngle = elevAngle;
		this.residual = residual;
	}

	public SatResidual(double t, double elevAngle, double residual, double noise) {
		super();
		this.t = t;
		this.elevAngle = elevAngle;
		this.residual = residual;
		this.noise = noise;
	}

	public SatResidual(double t, double elevAngle, double residual, boolean isOutlier) {
		super();
		this.t = t;
		this.elevAngle = elevAngle;
		this.residual = residual;
		this.isOutlier = isOutlier;
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

	public double getNoise() {
		return noise;
	}

	public boolean isOutlier() {
		return isOutlier;
	}
}
