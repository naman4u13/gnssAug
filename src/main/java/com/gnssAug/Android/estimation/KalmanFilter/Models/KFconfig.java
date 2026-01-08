package com.gnssAug.Android.estimation.KalmanFilter.Models;

import java.util.ArrayList;
import java.util.stream.IntStream;

import org.ejml.dense.row.MatrixFeatures_DDRM;
import org.ejml.simple.SimpleMatrix;

import com.gnssAug.Android.constants.ClockAllanVar;
import com.gnssAug.Android.constants.GnssDataConfig;
import com.gnssAug.utility.LatLonUtil;
import com.gnssAug.utility.Matrix;

public class KFconfig extends KF {

	private final double SpeedofLight = 299792458;
	private final double c2 = SpeedofLight * SpeedofLight;
	// Typical Allan Variance Coefficients for TCXO (low quality)
	private final double sf = 1;//ClockAllanVar.TCXO_low_quality.sf;
	private final double sg = 0.1;//ClockAllanVar.TCXO_low_quality.sg;

	public void config(double deltaT, Flag flag, int m, boolean useDoppler, boolean complementary, boolean useEstVel)
			throws Exception {

		/*
		 * The process noise for position vector will be initialized in ENU frame and
		 * will then be changed to ECEF frame. Rotation matrix 'R' will be computed to
		 * perform the coordinate transform.
		 */
		double[] ecef = new double[] { getState().get(0), getState().get(1), getState().get(2) };

		if (flag == Flag.POSITION) {
			int n = 3 + (2 * m);
			double[][] phi = new double[n][n];
			double[][] _Q = new double[n][n];
			IntStream.range(0, n).forEach(i -> phi[i][i] = 1);

			double[] qENU = GnssDataConfig.qENU_posRandWalk;

			// qECEF_std can have negative element
			IntStream.range(0, 3).forEach(i -> _Q[i][i] = qENU[i] * deltaT);
			SimpleMatrix R = new SimpleMatrix(n, n);
			R.insertIntoThis(0, 0, LatLonUtil.getEnu2EcefRotMat(ecef));
			for (int i = 3; i < 3 + m; i++) {
				_Q[i][i] = (sf * deltaT) + ((sg * Math.pow(deltaT, 3)) / 3);
				_Q[i][i + m] = (sg * Math.pow(deltaT, 2)) / 2;
				_Q[i + m][i] = (sg * Math.pow(deltaT, 2)) / 2;
				_Q[i + m][i + m] = sg * deltaT;
				phi[i][i + m] = deltaT;
				R.set(i, i, 1);
				R.set(i + m, i + m, 1);
			}
			SimpleMatrix Q = new SimpleMatrix(_Q);
			Q = R.mult(Q).mult(R.transpose());
			if (!MatrixFeatures_DDRM.isPositiveDefinite(Q.getMatrix())) {

				throw new Exception("PositiveDefinite test Failed");
			}
			super.configure(phi, Q);

		} else if (flag == Flag.VELOCITY) {
			if ((useDoppler || useEstVel) && complementary) {
				int n = 6 + (2 * m);
				double[][] phi = new double[n][n];
				SimpleMatrix Q = new SimpleMatrix(n, n);
				IntStream.range(0, n).forEach(i -> phi[i][i] = 1);
				for (int i = 0; i < 3 + m; i++) {
					phi[i][i + 3 + m] = deltaT;
					Q.set(i + 3 + m, i + 3 + m, 1e12);
				}
				if (!useEstVel) {
					for (int i = m; i < 3 + m; i++) {
						Q.set(i, i, 1e5);
					}
				}
//				if (!MatrixFeatures_DDRM.isPositiveDefinite(Q.getMatrix())) {
//
//					throw new Exception("PositiveDefinite test Failed");
//				}
				super.configure(phi, Q);

			} else {
				int n = 6 + (2 * m);
				double[][] phi = new double[n][n];
				double[][] _Q = new double[n][n];
				IntStream.range(0, n).forEach(i -> phi[i][i] = 1);

				double[] qENU = GnssDataConfig.qENU_velRandWalk;
				// Samsung 29th double[] qENU = new double[] { 0.05, 0.03, 0.0001 };
				double[] q = new double[3 + m];
				IntStream.range(0, 3).forEach(i -> q[i] = qENU[i]);
				double _sf = useDoppler ? 25 : sf;
				double _sg = useDoppler ? GnssDataConfig.clkDriftVar : sg;
				IntStream.range(3, 3 + m).forEach(i -> q[i] = _sg);
				for (int i = 0; i < 3 + m; i++) {
					_Q[i][i] = q[i] * Math.pow(deltaT, 3) / 3;
					_Q[i][i + 3 + m] = q[i] * Math.pow(deltaT, 2) / 2;
					_Q[i + 3 + m][i] = q[i] * Math.pow(deltaT, 2) / 2;
					_Q[i + 3 + m][i + 3 + m] = q[i] * deltaT;
					phi[i][i + 3 + m] = deltaT;
				}
				IntStream.range(3, 3 + m).forEach(i -> _Q[i][i] += (_sf * deltaT));

				SimpleMatrix _R = LatLonUtil.getEnu2EcefRotMat(ecef);
				SimpleMatrix R = new SimpleMatrix(n, n);
				for (int i = 0; i < 3; i++) {
					for (int j = 0; j < 3; j++) {
						R.set(i, j, _R.get(i, j));
						R.set(i + 3 + m, j + 3 + m, _R.get(i, j));
					}
				}
				for (int i = 0; i < m; i++) {
					R.set(3 + i, 3 + i, 1);
					R.set(6 + m + i, 6 + m + i, 1);
				}
				SimpleMatrix Q = new SimpleMatrix(_Q);
				Q = R.mult(Q).mult(R.transpose());
				if (!MatrixFeatures_DDRM.isPositiveDefinite(Q.getMatrix())) {

					throw new Exception("PositiveDefinite test Failed");
				}
				super.configure(phi, Q);
			}
		}
	}

	public void configIGS(double deltaT, int m) throws Exception {
		int n = 3 + (2 * m);
		double[][] phi = new double[n][n];
		double[][] Q = new double[n][n];
		for (int i = 3; i < 3 + m; i++) {
			Q[i][i] = (sf * deltaT) + ((sg * Math.pow(deltaT, 3)) / 3);
			Q[i][i + m] = (sg * Math.pow(deltaT, 2)) / 2;
			Q[i + m][i] = (sg * Math.pow(deltaT, 2)) / 2;
			Q[i + m][i + m] = sg * deltaT;
			phi[i][i + m] = deltaT;
		}
		IntStream.range(0, n).forEach(x -> phi[x][x] = 1);
		super.configure(phi, Q);
	}

	public void configDoppler(double deltaT, SimpleMatrix Cxx_dot_hat, int m, SimpleMatrix X) {
		int n = 3 + m;

		// Extra noise
		double[] qENU = new double[3 + m];
//		for(int i=0;i<3;i++)
//		{
//			qENU[i] = 0.1;
//		}
		for (int i = 0; i < m; i++) {
			qENU[i + 3] = GnssDataConfig.clkDriftVar;
		}
		SimpleMatrix omega = new SimpleMatrix(n, n);
		for (int i = 0; i < n; i++) {
			omega.set(i, i, qENU[i]);
		}
		SimpleMatrix R = new SimpleMatrix(n, n);
		R.insertIntoThis(0, 0, LatLonUtil.getEnu2EcefRotMat(new double[] { X.get(0), X.get(1), X.get(2) }));
		for (int i = 3; i < n; i++) {
			R.set(i, i, 1);
		}
		omega = R.mult(omega).mult(R.transpose());

		double[][] phi = new double[n][n];
		SimpleMatrix _Q = new SimpleMatrix(n, n);
		_Q.insertIntoThis(0, 0, (Cxx_dot_hat.plus(omega)).scale(deltaT));
		IntStream.range(0, n).forEach(i -> phi[i][i] = 1);
		double[][] Q = Matrix.matrix2Array(_Q);
		super.configure(phi, Q);
	}

	public void configAKFStatic(double deltaT, int m) throws Exception {
		int n = 3 + m;
		double[][] phi = new double[n][n];
		double[][] _Q = new double[n][n];
		IntStream.range(0, n).forEach(i -> phi[i][i] = 1);

		double[] qENU = GnssDataConfig.qENU_posRandWalk;

		// qECEF_std can have negative element
		IntStream.range(0, 3).forEach(i -> _Q[i][i] = qENU[i] * deltaT);

		for (int i = 3; i < 3 + m; i++) {
			_Q[i][i] = (sf * deltaT) + ((sg * Math.pow(deltaT, 3)) / 3);

		}
		SimpleMatrix Q = new SimpleMatrix(_Q);

		if (!MatrixFeatures_DDRM.isPositiveDefinite(Q.getMatrix())) {

			throw new Exception("PositiveDefinite test Failed");
		}
		super.configure(phi, Q);
	}
	
	public void configTDCP(double deltaT, int m, double[] refPos) throws Exception {

		int n = 3 + (2 * m);
		double[][] phi = new double[n][n];

		// Identity transition for all states (Velocity and Drifts are Random Walks)
		IntStream.range(0, n).forEach(i -> phi[i][i] = 1);

		double[] qENU = GnssDataConfig.qENU_velRandWalk;

		SimpleMatrix _Q = new SimpleMatrix(n, n);

		// 1. Velocity Process Noise (Indices 0, 1, 2)
		IntStream.range(0, 3).forEach(i -> _Q.set(i, i, qENU[i] * deltaT));

		// 2. Clock Drift Process Noise (Doppler & TDCP)
		// Structure: [Vel(3) | Doppler_Drifts(m) | TDCP_Drifts(m)]
		
		// Base Drift Noise (Oscillator): High variance for Android TCXO
		double baseDriftVar = 1.0; 
		
		// Bias Drift Noise (Inter-System Bias Rate): Low variance for thermal drift
		double biasDriftVar = 1.0e-4; 

		// Doppler Drifts (Indices 3 to 3+m-1)
		for (int i = 0; i < m; i++) {
			int stateIdx = 3 + i;
			if (i == 0) {
				// Base State (GPS/Ref): High Noise
				_Q.set(stateIdx, stateIdx, baseDriftVar * deltaT);
			} else {
				// Bias State (GAL/BDS relative to GPS): Low Noise
				_Q.set(stateIdx, stateIdx, biasDriftVar * deltaT);
			}
		}

		// TDCP Drifts (Indices 3+m to 3+2m-1)
		for (int i = 0; i < m; i++) {
			int stateIdx = 3 + m + i;
			if (i == 0) {
				// Base State (GPS/Ref): High Noise
				_Q.set(stateIdx, stateIdx, baseDriftVar * deltaT);
			} else {
				// Bias State (GAL/BDS relative to GPS): Low Noise
				_Q.set(stateIdx, stateIdx, biasDriftVar * deltaT);
			}
		}

		// Rotation Matrix to convert Velocity ENU -> ECEF
		SimpleMatrix _R = LatLonUtil.getEnu2EcefRotMat(refPos);
		SimpleMatrix R = new SimpleMatrix(n, n);
		
		// Insert Rotation for Velocity
		R.insertIntoThis(0, 0, _R);
		
		// Identity Rotation for Clock States (Drifts are already scalar)
		for (int i = 0; i < (2 * m); i++) {
			R.set(3 + i, 3 + i, 1);
		}

		SimpleMatrix Q = R.mult(_Q).mult(R.transpose());
		if (!MatrixFeatures_DDRM.isPositiveDefinite(Q.getMatrix())) {
			throw new Exception("PositiveDefinite test Failed");
		}
		super.configure(phi, Q);
	}
	
	public void configPPP(double deltaT, int codeClkOffNum, int clkDriftNum, int totalStateNum,
			ArrayList<double[]> ionoParams, boolean isAndroid) throws Exception 
	{
		configPPP( deltaT,  codeClkOffNum,  clkDriftNum,  totalStateNum,
			 ionoParams, isAndroid, false,false);
	}
	public void configPPP(double deltaT, int codeClkOffNum, int clkDriftNum, int totalStateNum,
			ArrayList<double[]> ionoParams, boolean isAndroid,boolean predictPhaseClock,boolean singlePhaseClock) throws Exception 
	{
		int phaseClKOffNum = 0;
		if(isAndroid)
		{
			if(singlePhaseClock)
			{
				phaseClKOffNum = 1;
				
			}
			else {			
				phaseClKOffNum = codeClkOffNum;
			}
		}
		configPPP( deltaT,  codeClkOffNum,  clkDriftNum,  totalStateNum,
			 ionoParams, isAndroid, phaseClKOffNum, predictPhaseClock);
	}
	
	public void configPPP(double deltaT, int codeClkOffNum, int clkDriftNum, int totalStateNum,
			ArrayList<double[]> ionoParams,boolean isAndroid,int phaseClKOffNum,boolean predictPhaseClock) throws Exception {

		double[] refPos = new double[] { getState().get(0), getState().get(1), getState().get(2) };
		int ionoParamNum = ionoParams.size();
		int clkOffNum = codeClkOffNum + phaseClKOffNum;
		int driftStartIndex = 6 + clkOffNum;
		// Its 16 cm^2/s in TECU^2/s, assuming L1 freq
		final double TECU_var = 0.001;
		double[][] phi = new double[totalStateNum][totalStateNum];
		IntStream.range(0, totalStateNum).forEach(i -> phi[i][i] = 1);
		IntStream.range(0, 3).forEach(i -> phi[i][i + 3 + clkOffNum] = deltaT);
		
		double[] qENU = GnssDataConfig.qENU_velRandWalk;
		SimpleMatrix _Q = new SimpleMatrix(totalStateNum, totalStateNum);

		// Position and Velocity
		for (int i = 0; i < 3; i++) {
			_Q.set(i, i, (qENU[i] * Math.pow(deltaT, 3) / 3) + (1e-16));
			_Q.set(i, i + 3 + clkOffNum, qENU[i] * Math.pow(deltaT, 2) / 2);
			_Q.set(i + 3 + clkOffNum, i, qENU[i] * Math.pow(deltaT, 2) / 2);
			_Q.set(i + 3 + clkOffNum, i + 3 + clkOffNum, qENU[i] * deltaT);
		}

		double clkOffVar = 1;
		
		double diffCodeBiasVar = 1.0e-8; 
		double diffPhaseBiasVar = 1.0e-8; 
		// Clock Drift
		double clkDriftVar = 0.1;
		double diffClkDriftVar = 1.0e-4; 
		if (isAndroid) {
			clkOffVar = 1e2;
			clkDriftVar = 1;
			diffCodeBiasVar = 1.0e-4; 
			diffPhaseBiasVar = 1.0e-4; 
		}
		double phaseClkOffVar = 1e4;
		// Receiver Code Clock Offset
		_Q.set(3, 3, clkOffVar*deltaT);
		// Rx DCB
		for (int i = 1; i < codeClkOffNum; i++) {
			_Q.set(i + 3, i + 3, diffCodeBiasVar * deltaT);

		}
		if (isAndroid) {
			// Receiver Phase Clock Offset
			int phaseStartIndex = 3 + codeClkOffNum;
			_Q.set(phaseStartIndex, phaseStartIndex, phaseClkOffVar * deltaT);
			// Rx Differential Phase Bias
			for (int i = 1; i < phaseClKOffNum; i++) {
				_Q.set(phaseStartIndex + i, phaseStartIndex + i, diffPhaseBiasVar * deltaT);

			}
			if(predictPhaseClock)
			{
				// Only couple the Base Phase (Index 0) to the Base Drift (Index 0)
				// We do NOT couple ISB phases to ISB drifts; 1e-2 noise handles them.
				
				
				
				// Apply only for i = 0 (Base System)
				int phaseIdx = phaseStartIndex; // + 0
				int driftIdx = driftStartIndex; // + 0
				
				// 1. Transition: Base Phase(t+1) = Base Phase(t) + Base Drift(t) * dt
				phi[phaseIdx][driftIdx] = deltaT;
				
				// 2. Integrated Process Noise for Base State
				// (+1 is the Android safety floor for phase jitter)
				_Q.set(phaseIdx, phaseIdx, (clkDriftVar * Math.pow(deltaT, 3) / 3) + 1);
				
				// 3. Covariance coupling
				_Q.set(phaseIdx, driftIdx, clkDriftVar * Math.pow(deltaT, 2) / 2);
				_Q.set(driftIdx, phaseIdx, clkDriftVar * Math.pow(deltaT, 2) / 2);
			}
		}
		
		for (int i = 0; i < clkDriftNum; i++) {
			double varToUse;
			if (i == 0) {
				// Base Oscillator Drift (High Dynamics)
				varToUse = clkDriftVar;
			} else {
				// Inter-System Bias Drift (Slow Thermal Variation)
				varToUse = diffClkDriftVar;
			}
			
			_Q.set(driftStartIndex + i, driftStartIndex + i, varToUse * deltaT);
		}
		// Tropo: More than 1cm/sqrt(hr)
		_Q.set(6 + clkOffNum + clkDriftNum, 6 + clkOffNum + clkDriftNum, (1e-8) * deltaT);

		// Ambiguities
		for (int i = 6 + clkOffNum + clkDriftNum + 1; i < totalStateNum - ionoParamNum; i++) {
			_Q.set(i, i, 1e-16);
		}

		// Ionosphere: 4 cm/sqrt(s)*sin(elevation)
		for (int i = 0; i < ionoParamNum; i++) {
			_Q.set(totalStateNum - ionoParamNum + i, totalStateNum - ionoParamNum + i,
					(TECU_var * deltaT) / Math.pow(Math.sin(ionoParams.get(i)[0]), 2));
		}

		SimpleMatrix _R = LatLonUtil.getEnu2EcefRotMat(refPos);
		SimpleMatrix R = new SimpleMatrix(totalStateNum, totalStateNum);
		for (int i = 0; i < totalStateNum; i++) {
			R.set(i, i, 1);
		}
		R.insertIntoThis(0, 0, _R);
		R.insertIntoThis(3 + clkOffNum, 3 + clkOffNum, _R);

		SimpleMatrix Q = R.mult(_Q).mult(R.transpose());
		if (!MatrixFeatures_DDRM.isPositiveDefinite(Q.getMatrix())) {
			throw new Exception("PositiveDefinite test Failed");
		}
		Q = Q.plus(Q.transpose()).scale(0.5);
		super.configure(phi, Q);
	}
}
