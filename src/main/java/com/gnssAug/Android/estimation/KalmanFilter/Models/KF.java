package com.gnssAug.Android.estimation.KalmanFilter.Models;

import org.ejml.simple.SimpleMatrix;

public class KF {

	private final double SpeedofLight = 299792458;

	// kinematics description
	private SimpleMatrix F, Q;

	// sytem state estimate
	private SimpleMatrix x, P;

	public void configure(double[][] F, double[][] Q) {
		this.F = new SimpleMatrix(F);
		this.Q = new SimpleMatrix(Q);

	}

	public void setState_ProcessCov(double[][] x, double[][] P) {
		this.x = new SimpleMatrix(x);
		this.P = new SimpleMatrix(P);
	}

	// Prediction step
	public void predict() {
		// x = F x
		x = F.mult(x);

		// P = F P F' + Q
		P = F.mult(P).mult(F.transpose()).plus(Q);

	}

	// Update Step
	public void update(double[][] _z, double[][] _R, double[][] _ze, double[][] _H) {

		SimpleMatrix z = new SimpleMatrix(_z);
		SimpleMatrix R = new SimpleMatrix(_R);
		SimpleMatrix ze = new SimpleMatrix(_ze);
		SimpleMatrix H = new SimpleMatrix(_H);
		SimpleMatrix Ht = H.transpose();
		SimpleMatrix K = null;

		// Kalman Gain
		K = P.mult(Ht).mult(((H.mult(P).mult(Ht)).plus(R)).invert());

		// Posterior State Estimate
		x = x.plus((K.mult(z.minus(ze))));
		SimpleMatrix KH = K.mult(H);
		SimpleMatrix I = SimpleMatrix.identity(KH.numRows());
		/*
		 * Posterior Estimate Error Joseph Form to ensure Positive Definiteness P =
		 * (I-KH)P(I-KH)' + KRK'
		 */
		P = ((I.minus(KH)).mult(P).mult((I.minus(KH)).transpose())).plus(K.mult(R).mult(K.transpose()));

	}

	public SimpleMatrix getState() {
		return x;
	}

	public SimpleMatrix getCovariance() {
		return P;
	}

}
