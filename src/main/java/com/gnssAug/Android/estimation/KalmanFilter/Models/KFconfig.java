package com.gnssAug.Android.estimation.KalmanFilter.Models;

import java.util.stream.IntStream;

import org.ejml.dense.row.MatrixFeatures_DDRM;
import org.ejml.simple.SimpleMatrix;

import com.gnssAug.Android.constants.ClockAllanVar;
import com.gnssAug.utility.LatLonUtil;

public class KFconfig extends KF {

	private final double SpeedofLight = 299792458;
	private final double c2 = SpeedofLight * SpeedofLight;
	// Typical Allan Variance Coefficients for TCXO (low quality)
	private final double h0 = ClockAllanVar.TCXO_low_quality.h0;
	private final double h_2 = ClockAllanVar.TCXO_low_quality.h_2;
	private final double sf = ClockAllanVar.TCXO_low_quality.sf;
	private final double sg = ClockAllanVar.TCXO_low_quality.sg;

	public void config(double deltaT, Flag flag) throws Exception {

		/*
		 * The process noise for position vector will be initialized in ENU frame and
		 * will then be changed to ECEF frame. Rotation matrix 'R' will be computed to
		 * perform the coordinate transform.
		 */
		double[] ecef = new double[] { getState().get(0), getState().get(1), getState().get(2) };

		if (flag == Flag.POSITION) {

			double[][] phi = new double[5][5];
			double[][] _Q = new double[5][5];
			IntStream.range(0, 5).forEach(i -> phi[i][i] = 1);
			phi[3][4] = deltaT;

//			double[] qENU = new double[] { 12, 12, 0.2 };
			double[] qENU = new double[] { 16, 16, 4 };
			// qECEF_std can have negative element
			IntStream.range(0, 3).forEach(i -> _Q[i][i] = qENU[i]);
			_Q[3][3] = ((sf * deltaT) + ((sg * Math.pow(deltaT, 3)) / 3));
			_Q[3][4] = (sg * Math.pow(deltaT, 2)) / 2;
			_Q[4][3] = (sg * Math.pow(deltaT, 2)) / 2;
			_Q[4][4] = sg * deltaT;
			SimpleMatrix _R = LatLonUtil.getEnu2EcefRotMat(ecef);
			SimpleMatrix R = new SimpleMatrix(5, 5);
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					R.set(i, j, _R.get(i, j));
				}
			}
			R.set(3, 3, 1);
			R.set(4, 4, 1);
			SimpleMatrix Q = new SimpleMatrix(_Q);
			Q = R.mult(Q).mult(R.transpose());
			if (!MatrixFeatures_DDRM.isPositiveDefinite(Q.getMatrix())) {

				throw new Exception("PositiveDefinite test Failed");
			}
			super.configure(phi, Q);

		} else if (flag == Flag.VELOCITY) {
			double[][] F = new double[8][8];
			double[][] _Q = new double[8][8];
			IntStream.range(0, 8).forEach(i -> F[i][i] = 1);
			IntStream.range(0, 4).forEach(i -> F[i][i + 4] = deltaT);
			// double[] qENU_std = new double[] { 8, 12, 2 };

			double[] qENU = new double[] { 0.25, 0.25, 0.1 };
			double[] q = new double[] { qENU[0], qENU[1], qENU[2], sg };
			for (int i = 0; i < 4; i++) {
				_Q[i][i] = q[i] * Math.pow(deltaT, 3) / 3;
				_Q[i][i + 4] = q[i] * Math.pow(deltaT, 2) / 2;
				_Q[i + 4][i] = q[i] * Math.pow(deltaT, 2) / 2;
				_Q[i + 4][i + 4] = q[i] * deltaT;
			}
			_Q[3][3] += (sf * deltaT);
			SimpleMatrix _R = LatLonUtil.getEnu2EcefRotMat(ecef);
			SimpleMatrix R = new SimpleMatrix(8, 8);
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					R.set(i, j, _R.get(i, j));
					R.set(i + 4, j + 4, _R.get(i, j));
				}
			}
			R.set(3, 3, 1);
			R.set(7, 7, 1);
			SimpleMatrix Q = new SimpleMatrix(_Q);
			Q = R.mult(Q).mult(R.transpose());
			if (!MatrixFeatures_DDRM.isPositiveDefinite(Q.getMatrix())) {

				throw new Exception("PositiveDefinite test Failed");
			}
			super.configure(F, Q);

		}
	}
	
	public void configIGS(double deltaT) throws Exception {
		double[][] phi = new double[5][5];
		double[][] Q = new double[5][5];
		Q[3][3] = (sf * deltaT) + ((sg * Math.pow(deltaT, 3)) / 3);
		Q[3][4] = (sg * Math.pow(deltaT, 2))/ 2;
		Q[4][3] = (sg * Math.pow(deltaT, 2)) / 2;
		Q[4][4] = sg * deltaT ;

		IntStream.range(0, 5).forEach(x -> phi[x][x] = 1);
		phi[3][4] = deltaT;
		super.configure(phi, Q);
	}

}
