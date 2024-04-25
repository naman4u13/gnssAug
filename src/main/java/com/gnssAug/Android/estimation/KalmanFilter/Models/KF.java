package com.gnssAug.Android.estimation.KalmanFilter.Models;

import org.ejml.simple.SimpleMatrix;

public class KF {

	private final double SpeedofLight = 299792458;

	// kinematics description
	private SimpleMatrix phi, Q;

	// sytem state estimate
	private SimpleMatrix x, P;

	// Kalman Gain
	private SimpleMatrix K;

	// Covariance of prefit residual or innovation
	private SimpleMatrix Cvv;

	public void configure(double[][] phi, double[][] Q) {
		this.phi = new SimpleMatrix(phi);
		this.Q = new SimpleMatrix(Q);

	}

	public void configure(double[][] phi, SimpleMatrix Q) {
		this.phi = new SimpleMatrix(phi);
		this.Q = new SimpleMatrix(Q);

	}

	public void setState_ProcessCov(double[][] x, double[][] P) {
		this.x = new SimpleMatrix(x);
		this.P = new SimpleMatrix(P);
	}
	
	public void setState_ProcessCov(SimpleMatrix x, SimpleMatrix P) {
		this.x = x;
		this.P = P;
	}
	
	public void setProcessCov(SimpleMatrix P) {
		this.P = P;
	}

	// Prediction step
	public void predict() {
		// x = phi x
		x = phi.mult(x);

		// P = phi P phi' + Q
		P = phi.mult(P).mult(phi.transpose()).plus(Q);

	}

	// Update Step
	public void update(double[][] _z, double[][] _R, double[][] _ze, double[][] _H) {

		SimpleMatrix R = new SimpleMatrix(_R);
		SimpleMatrix H = new SimpleMatrix(_H);
		update(_z, R, _ze, H);
	}

	// Update Step
	public void update(double[][] _z, SimpleMatrix R, double[][] _ze, SimpleMatrix H) {

		SimpleMatrix z = new SimpleMatrix(_z);
		SimpleMatrix ze = new SimpleMatrix(_ze);
		update(z, R, ze, H);

	}
	
	// Update Step
		public void update(SimpleMatrix z, SimpleMatrix R, SimpleMatrix ze, SimpleMatrix H) {

			
			SimpleMatrix Ht = H.transpose();

			Cvv = ((H.mult(P).mult(Ht)).plus(R));
			// Kalman Gain
			K = P.mult(Ht).mult(Cvv.invert());

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

	public SimpleMatrix getKalmanGain() {
		return K;
	}

	public SimpleMatrix getPhi() {
		return phi;
	}

	public SimpleMatrix getCvv() {
		return Cvv;
	}

	public SimpleMatrix getQ() {
		return Q;
	}

}
