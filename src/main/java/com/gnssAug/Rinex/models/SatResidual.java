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
}
