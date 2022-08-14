package com.gnssAug.Rinex.models;

public class SatResidual {
	// GPS Time
	private double t;
	// Elevation Angle(Radians)
	private double elevAngle;
	// Residual(meter)
	private double residual;

	public SatResidual(double t, double elevAngle, double residual) {
		super();
		this.t = t;
		this.elevAngle = elevAngle;
		this.residual = residual;
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
}
