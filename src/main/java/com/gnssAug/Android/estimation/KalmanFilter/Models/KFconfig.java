package com.gnssAug.Android.estimation.KalmanFilter.Models;

import java.util.stream.IntStream;

import com.gnssAug.Android.constants.ClockAllanVar;
import com.gnssAug.Android.utility.LatLonUtil;

public class KFconfig extends KF {

	private final double SpeedofLight = 299792458;
	private final double c2 = SpeedofLight * SpeedofLight;
	// Typical Allan Variance Coefficients for TCXO (low quality)
	private final double h0 = ClockAllanVar.TCXO_low_quality.h0;
	private final double h_2 = ClockAllanVar.TCXO_low_quality.h_2;
	private final double sf = ClockAllanVar.TCXO_low_quality.sf;
	private final double sg = ClockAllanVar.TCXO_low_quality.sg;

	public void config(double deltaT, Flag flag) {

		/*
		 * The process noise for position vector will be initialized in ENU frame and
		 * will then be changed to ECEF frame. Rotation matrix 'R' will be computed to
		 * perform the coordinate transform.
		 */
		double[] ecef = new double[] { getState().get(0), getState().get(1), getState().get(2) };

		if (flag == Flag.POSITION) {

			double[][] F = new double[5][5];
			double[][] Q = new double[5][5];
			IntStream.range(0, 5).forEach(i -> F[i][i] = 1);
			F[3][4] = deltaT;

			double[] qENU_std = new double[] { 5, 5, .1 };
			// qECEF_std can have negative element
			double[] qECEF_std = LatLonUtil.enu2ecef(qENU_std, ecef, false);
			double[] q = IntStream.range(0, 3).mapToDouble(i -> Math.pow(qECEF_std[i], 2)).toArray();
			IntStream.range(0, 3).forEach(i -> Q[i][i] = q[i]);
			Q[3][3] = ((sf * deltaT) + ((sg * Math.pow(deltaT, 3)) / 3));
			Q[3][4] = (sg * Math.pow(deltaT, 2)) / 2;
			Q[4][3] = (sg * Math.pow(deltaT, 2)) / 2;
			Q[4][4] = sg * deltaT;
			super.configure(F, Q);

		} else if (flag == Flag.VELOCITY) {
			double[][] F = new double[8][8];
			double[][] Q = new double[8][8];
			IntStream.range(0, 8).forEach(i -> F[i][i] = 1);
			IntStream.range(0, 4).forEach(i -> F[i][i + 4] = deltaT);
			double[] qENU_std = new double[] { 5, 5, 2 };
			// qECEF_std can have negative element
			double[] qECEF_std = LatLonUtil.enu2ecef(qENU_std, ecef, false);
			double[] qECEF = IntStream.range(0, 3).mapToDouble(i -> Math.pow(qECEF_std[i], 2)).toArray();
			double[] q = new double[] { qECEF[0], qECEF[1], qECEF[2], sg };
			for (int i = 0; i < 4; i++) {
				Q[i][i] = q[i] * Math.pow(deltaT, 3) / 3;
				Q[i][i + 4] = q[i] * Math.pow(deltaT, 2) / 2;
				Q[i + 4][i] = q[i] * Math.pow(deltaT, 2) / 2;
				Q[i + 4][i + 4] = q[i] * deltaT;
			}
			Q[3][3] += (sf * deltaT);
			super.configure(F, Q);
		}
	}

}
