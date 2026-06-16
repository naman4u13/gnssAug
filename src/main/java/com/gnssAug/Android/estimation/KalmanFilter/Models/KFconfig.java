package com.gnssAug.Android.estimation.KalmanFilter.Models;

import java.util.ArrayList;
import java.util.stream.IntStream;

import org.ejml.dense.row.MatrixFeatures_DDRM;
import org.ejml.simple.SimpleMatrix;

import com.gnssAug.Android.constants.ClockAllanVar;
import com.gnssAug.Android.constants.GnssDataConfig;
import com.gnssAug.utility.LatLonUtil;
import com.gnssAug.utility.Matrix;

/**
 * Builds and applies the Kalman filter process model (phi, Q) for each predict step.
 *
 * <p>Each {@code configXxx} method assembles a state-transition matrix (phi) and
 * process-noise covariance (Q) for a specific filter mode, then delegates to
 * {@link KF#configure} to push them into the parent filter.</p>
 *
 * <p>Process noise is always constructed in ENU and rotated to ECEF via
 * {@link LatLonUtil#getEnu2EcefRotMat} before being handed to the filter.
 * Two PPP variants exist — {@link #configPPP_IGS} for geodetic receivers and
 * {@link #configPPP_Android} for smartphones — sharing a private core with
 * receiver-class-specific noise constants passed in as parameters.</p>
 */
public class KFconfig extends KF {

	private final double SpeedofLight = 299792458;
	private final double c2 = SpeedofLight * SpeedofLight;
	// Typical Allan Variance Coefficients for TCXO (low quality)
	private final double sf = 1;//ClockAllanVar.TCXO_low_quality.sf;
	private final double sg = 0.1;//ClockAllanVar.TCXO_low_quality.sg;

	/**
	 * Configures the process model for the Android SPP position / velocity filter.
	 *
	 * <p>State layout (POSITION): [pos(3) | clk(m) | clkDrift(m)].
	 * State layout (VELOCITY): [pos(3) | clk(m) | vel(3) | clkDrift(m)].
	 * Clock noise is derived from the two-state Allan-variance model (sf, sg) for a TCXO.
	 * The complementary VELOCITY branch treats the incoming velocity as near-perfect and
	 * injects a high Q on the velocity sub-states rather than integrating kinematics.</p>
	 *
	 * @param deltaT        time step (s)
	 * @param flag          POSITION or VELOCITY mode
	 * @param m             number of clock offset states (one per constellation)
	 * @param useDoppler    true when Doppler is the primary velocity source
	 * @param complementary true for complementary (position + velocity) mode
	 * @param useEstVel     true when an estimated velocity drives the complementary branch
	 */
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

	/**
	 * Configures the process model for an IGS-grade SPP position filter.
	 *
	 * <p>State layout: [pos(3) | clk(m) | clkDrift(m)].  Clock process noise is
	 * built from the two-state Allan-variance model (sf = white FM, sg = random-walk
	 * FM) suited to a stable geodetic oscillator.</p>
	 *
	 * @param deltaT time step (s)
	 * @param m      number of clock states (one per signal / constellation)
	 */
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

	/**
	 * Configures the process model for the Android Doppler velocity filter.
	 *
	 * <p>Q is built from the posterior velocity covariance of the previous epoch
	 * (Cxx_dot_hat) plus a diagonal extra-noise floor (omega) for clock drift
	 * dynamics.  The combined matrix is scaled by deltaT and inserted into the
	 * velocity block; the result is rotated ENU→ECEF.</p>
	 *
	 * @param deltaT       time step (s)
	 * @param Cxx_dot_hat  posterior velocity covariance from the previous epoch
	 * @param m            number of clock drift states
	 * @param X            current state vector (provides ECEF position for rotation)
	 */
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

	/**
	 * Configures the process model for the Android AKF static position filter.
	 *
	 * <p>State layout: [pos(3) | clk(m)].  No velocity states are included;
	 * position is modelled as a near-static random walk (tiny qENU).
	 * Clock noise uses Allan-variance coefficients without the cross-coupling
	 * drift term used in {@link #config}.</p>
	 *
	 * @param deltaT time step (s)
	 * @param m      number of clock states
	 */
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

	/**
	 * Configures the process model for the Android TDCP delta-position filter.
	 *
	 * <p>State layout: [vel(3) | Doppler_drifts(m) | TDCP_drifts(m)].
	 * The base oscillator drift (index 0 of each group) carries high process noise
	 * to track TCXO frequency jumps; inter-system bias drifts (index 1+) carry
	 * low process noise reflecting slow thermal variation.</p>
	 *
	 * @param deltaT  time step (s)
	 * @param m       number of clock drift states per measurement type
	 * @param refPos  current ECEF position (used to build the ENU rotation matrix)
	 */
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

	/**
	 * Configures the PPP process model for a geodetic (IGS) receiver.
	 *
	 * <p>State layout:
	 * [pos(3) | code_clk(codeClkOffNum) | vel(3) | clkDrift(clkDriftNum) | tropo(1) | amb(n) | iono(n)]
	 *
	 * <p>Uses tight noise constants appropriate for a stable geodetic oscillator and
	 * choke-ring / high-grade antenna: base clock σ = 1 m, code DCB σ = √(1e-8) m/√s.
	 * No separate phase-clock states are included; phase ISBs are absorbed into
	 * float carrier-phase ambiguities.</p>
	 *
	 * @param deltaT        time step (s)
	 * @param codeClkOffNum number of code clock offset states (1 base + n−1 ISBs)
	 * @param clkDriftNum   number of clock drift states
	 * @param totalStateNum total length of the state vector
	 * @param ionoParams    per-satellite list; each entry is {elevation_rad, ...},
	 *                      elevation used for the sin²(el) iono noise scaling
	 */
	public void configPPP_IGS(double deltaT, int codeClkOffNum, int clkDriftNum,
			int totalStateNum, ArrayList<double[]> ionoParams) throws Exception {
		configPPP_core(deltaT, codeClkOffNum, 0, clkDriftNum, totalStateNum, ionoParams,
				/*clkOffVar=*/        1.0,
				/*diffCodeBiasVar=*/  1.0e-8,
				/*diffPhaseBiasVar=*/ 1.0e-8,
				/*clkDriftVar=*/      0.1,
				/*diffClkDriftVar=*/  1.0e-4,
				/*predictPhaseClock=*/false);
	}

	/**
	 * Configures the PPP process model for an Android smartphone receiver.
	 *
	 * <p>State layout:
	 * [pos(3) | code_clk(codeClkOffNum) | phase_clk(phaseClKOffNum) | vel(3) |
	 *  clkDrift(clkDriftNum) | tropo(1) | amb(n) | iono(n)]
	 *
	 * <p>Uses looser noise constants reflecting a TCXO: base clock σ = 10 m,
	 * code/phase DCB σ = √(1e-4) m/√s.  A separate phase-clock block is appended
	 * after the code clocks to track the independent phase behaviour seen on Android
	 * hardware.</p>
	 *
	 * @param deltaT           time step (s)
	 * @param codeClkOffNum    number of code clock offset states (1 base + n−1 ISBs)
	 * @param clkDriftNum      number of clock drift states
	 * @param totalStateNum    total length of the state vector
	 * @param ionoParams       per-satellite list; each entry is {elevation_rad, ...}
	 * @param predictPhaseClock if true, couples the base phase state to the base drift
	 *                          state via phi[phase][drift] = dt (Wiener-process clock)
	 * @param singlePhaseClock  if true, only one phase clock state (base GPS) is added;
	 *                          if false, one phase clock per code system is added
	 */
	public void configPPP_Android(double deltaT, int codeClkOffNum, int clkDriftNum,
			int totalStateNum, ArrayList<double[]> ionoParams,
			boolean predictPhaseClock, boolean singlePhaseClock) throws Exception {
		int phaseClKOffNum = singlePhaseClock ? 1 : codeClkOffNum;
		configPPP_core(deltaT, codeClkOffNum, phaseClKOffNum, clkDriftNum, totalStateNum, ionoParams,
				/*clkOffVar=*/        1.0e2,
				/*diffCodeBiasVar=*/  1.0e-4,
				/*diffPhaseBiasVar=*/ 1.0e-4,
				/*clkDriftVar=*/      1.0,
				/*diffClkDriftVar=*/  1.0e-4,
				/*predictPhaseClock=*/predictPhaseClock);
	}

	// Shared PPP process model; receiver-class noise constants injected by callers above.
	private void configPPP_core(double deltaT, int codeClkOffNum, int phaseClKOffNum,
			int clkDriftNum, int totalStateNum, ArrayList<double[]> ionoParams,
			double clkOffVar, double diffCodeBiasVar, double diffPhaseBiasVar,
			double clkDriftVar, double diffClkDriftVar,
			boolean predictPhaseClock) throws Exception {

		double[] refPos = new double[] { getState().get(0), getState().get(1), getState().get(2) };
		int ionoParamNum = ionoParams.size();
		int clkOffNum = codeClkOffNum + phaseClKOffNum;
		int driftStartIndex = 6 + clkOffNum;
		// 16 cm²/s iono rate in TECU²/s at L1 frequency
		final double TECU_var = 0.001;
		final double phaseClkOffVar = 1e4;

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

		// Base receiver code clock
		_Q.set(3, 3, clkOffVar * deltaT);
		// Receiver code DCBs (inter-system biases relative to base clock)
		for (int i = 1; i < codeClkOffNum; i++) {
			_Q.set(i + 3, i + 3, diffCodeBiasVar * deltaT);
		}

		// Phase clock states (Android only, when phaseClKOffNum > 0)
		if (phaseClKOffNum > 0) {
			int phaseStartIndex = 3 + codeClkOffNum;
			_Q.set(phaseStartIndex, phaseStartIndex, phaseClkOffVar * deltaT);
			// Phase inter-system bias drifts
			for (int i = 1; i < phaseClKOffNum; i++) {
				_Q.set(phaseStartIndex + i, phaseStartIndex + i, diffPhaseBiasVar * deltaT);
			}
			if (predictPhaseClock) {
				// Couple base phase to base drift; ISB phases handled by noise floor
				int phaseIdx = phaseStartIndex;
				int driftIdx = driftStartIndex;
				phi[phaseIdx][driftIdx] = deltaT;
				// +1 is the Android safety floor for phase jitter
				_Q.set(phaseIdx, phaseIdx, (clkDriftVar * Math.pow(deltaT, 3) / 3) + 1);
				_Q.set(phaseIdx, driftIdx, clkDriftVar * Math.pow(deltaT, 2) / 2);
				_Q.set(driftIdx, phaseIdx, clkDriftVar * Math.pow(deltaT, 2) / 2);
			}
		}

		// Clock drifts: i=0 is the base oscillator (high dynamics); i>0 are inter-system bias drifts (slow thermal)
		for (int i = 0; i < clkDriftNum; i++) {
			double varToUse = (i == 0) ? clkDriftVar : diffClkDriftVar;
			_Q.set(driftStartIndex + i, driftStartIndex + i, varToUse * deltaT);
		}

		// Troposphere: >1 cm/√hr random walk on zenith wet delay
		_Q.set(6 + clkOffNum + clkDriftNum, 6 + clkOffNum + clkDriftNum, (1e-8) * deltaT);

		// Ambiguities: near-constant between cycle slips; tiny noise prevents singularity
		for (int i = 6 + clkOffNum + clkDriftNum + 1; i < totalStateNum - ionoParamNum; i++) {
			_Q.set(i, i, 1e-16);
		}

		// Ionosphere: ~4 cm/√s · sin(el) random walk on STEC; scaled by sin²(el) to give slant noise
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
