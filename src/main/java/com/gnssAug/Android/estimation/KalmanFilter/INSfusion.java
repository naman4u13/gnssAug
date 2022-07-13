package com.gnssAug.Android.estimation.KalmanFilter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.stream.IntStream;

import org.ejml.simple.SimpleMatrix;

import com.gnssAug.Android.constants.AndroidSensor;
import com.gnssAug.Android.constants.ClockAllanVar;
import com.gnssAug.Android.constants.ImuDataSheets;
import com.gnssAug.Android.estimation.LinearLeastSquare;
import com.gnssAug.Android.estimation.KalmanFilter.Models.State;
import com.gnssAug.Android.helper.GeomagneticField;
import com.gnssAug.Android.helper.Rotation;
import com.gnssAug.Android.models.IMUsensor;
import com.gnssAug.Android.models.Satellite;
import com.gnssAug.Android.utility.LatLonUtil;
import com.gnssAug.Android.utility.Matrix;

public class INSfusion {
	private final static double SpeedofLight = 299792458;

	public static TreeMap<Long, double[]> process(TreeMap<Long, HashMap<AndroidSensor, IMUsensor>> imuMap,
			TreeMap<Long, ArrayList<Satellite>> SatMap, ArrayList<Long> timeList, double[][] _dcm) throws Exception {
		ArrayList<ArrayList<Satellite>> SVlist = new ArrayList<ArrayList<Satellite>>(SatMap.values());
		// Convert DCM from Body Frame to ENU to Body frame to NED
		SimpleMatrix dcm = new SimpleMatrix(LatLonUtil.enu_ned_convert(_dcm));
		dcm = Rotation.reorthonormDcm(dcm);
		double[] ecef0 = LinearLeastSquare.process(SVlist.get(0), true);
		double[] llh0 = LatLonUtil.ecef2lla(ecef0);
		IntStream.range(0, 2).forEach(i -> llh0[i] = Math.toRadians(llh0[i]));
		// Velocity in ENU frame, zero initially
		double[] vel0 = new double[3];
		double rxClkOff = SpeedofLight * ecef0[3];
		State X = new State(llh0[0], llh0[1], llh0[2], vel0[0], vel0[1], vel0[2], dcm, 0, 0, 0, 0, 0, 0, rxClkOff, 0);
		// Attitude std deviation is 20 degree, values mentioned below in covariance
		// matrix is in radians
		double accBiasCov = Math.pow(ImuDataSheets.Pixel4.accTurnOnBias, 2);
		double gyroBiasCov = Math.pow(ImuDataSheets.Pixel4.gyroTurnOnBias, 2);
		double[] p0 = new double[] { 100, 100, 100, 0.1, 0.1, 0.1, 0.1156, 0.1156, 0.1156, accBiasCov, accBiasCov,
				accBiasCov, gyroBiasCov, gyroBiasCov, gyroBiasCov, 5, 0.005 };
		SimpleMatrix P = new SimpleMatrix(17, 17);
		for (int i = 0; i < 17; i++) {
			P.set(i, i, p0[i]);
		}

		// PSD of random walk - N
		double acc_SN = Math.pow(ImuDataSheets.Pixel4.vrw, 2);
		double gyro_SN = Math.pow(ImuDataSheets.Pixel4.arw, 2);
		// PSD of Bias Instability - B
		double acc_SB = 2 * Math.pow(ImuDataSheets.Pixel4.accInRunBias, 2) * Math.log(2)
				/ (Math.PI * Math.pow(0.4365, 2) * ImuDataSheets.Pixel4.accBiasCorrTime);
		double gyro_SB = 2 * Math.pow(ImuDataSheets.Pixel4.gyroInRunBias, 2) * Math.log(2)
				/ (Math.PI * Math.pow(0.4365, 2) * ImuDataSheets.Pixel4.gyroBiasCorrTime);
		// Receiver Clock phase and frequency PSD
		double sf = ClockAllanVar.TCXO_low_quality.sf;
		double sg = ClockAllanVar.TCXO_low_quality.sg;
		double[] q = new double[] { acc_SN, acc_SN, acc_SN, gyro_SN, gyro_SN, gyro_SN, acc_SB, acc_SB, acc_SB, gyro_SB,
				gyro_SB, gyro_SB, sf, sg };
		int n = SVlist.size();
		if (n != timeList.size()) {
			System.err.println("FATAL ERROR in INSfusion.java");
			throw new Exception("size of SVlist and timelist does not match");
		}

		int i = 0;
		Iterator<Entry<Long, HashMap<AndroidSensor, IMUsensor>>> iterator = imuMap.entrySet().iterator();
		double prevTime = iterator.next().getKey();
		while (timeList.get(i) <= prevTime) {
			i++;
		}
		long utcTimeMilli = SVlist.get(i).get(0).getUtcTimeMillis();
		TreeMap<Long, double[]> ecefMap = new TreeMap<Long, double[]>();
		while ((i < n) && iterator.hasNext()) {

			Entry<Long, HashMap<AndroidSensor, IMUsensor>> entry = iterator.next();
			long time = entry.getKey();
			HashMap<AndroidSensor, IMUsensor> imuSensor = entry.getValue();
			double tau = (time - prevTime) * 1e-3;
			// Update Total State and get estimated INS observables - {acc,gyro}
			double[][] estInsObs = predictTotalState(X, imuSensor, tau);
			prevTime = time;
			/*
			 * Get transition matrix Phi and covariance matrix Qk, for discrete state space
			 * filtering
			 */
			SimpleMatrix[] discParam = getDiscreteParams(X, estInsObs[0], estInsObs[1], tau, q);
			P = predictErrorState(X, P, discParam[0], discParam[1]);
			boolean onlyMagneto = true;
			ArrayList<Satellite> SV = null;
			if (time == timeList.get(i)) {
				SV = SVlist.get(i);
				utcTimeMilli = SV.get(0).getUtcTimeMillis();
				onlyMagneto = false;
				i++;

			}
			P = update(X, P, utcTimeMilli, imuSensor.get(AndroidSensor.Magnetometer), SV, onlyMagneto);
			ecefMap.put(time, LatLonUtil.lla2ecef(X.getP(), false));
		}
		return ecefMap;

	}

	private static double[][] predictTotalState(State X, HashMap<AndroidSensor, IMUsensor> imuSensor, double tau) {

		// Update Bias state
		double[] accBias = Arrays.stream(X.getAccBias())
				.map(i -> i * Math.exp(-tau / ImuDataSheets.Pixel4.accBiasCorrTime)).toArray();
		double[] gyroBias = Arrays.stream(X.getGyroBias())
				.map(i -> i * Math.exp(-tau / ImuDataSheets.Pixel4.gyroBiasCorrTime)).toArray();

		double[] obsAcc = imuSensor.get(AndroidSensor.Accelerometer).getVal();
		double[] estAcc = IntStream.range(0, 3).mapToDouble(j -> obsAcc[j] - accBias[j]).toArray();
		double[] obsGyro = imuSensor.get(AndroidSensor.Gyroscope).getVal();
		double[] estGyro = IntStream.range(0, 3).mapToDouble(j -> obsGyro[j] - gyroBias[j]).toArray();

		double lat = X.getP()[0];
		double lon = X.getP()[1];
		double alt = X.getP()[2];
		double[] vel = X.getV();
		double Rn = LatLonUtil.getNormalEarthRadius(lat);
		double Rm = LatLonUtil.getMeridianEarthRadius(lat);
		double earthAngularRate = LatLonUtil.omega_ie;

		// Attitude update
		SimpleMatrix oldDcm = new SimpleMatrix(X.getDcm());
		SimpleMatrix omega_b_ib = new SimpleMatrix(Matrix.getSkewSymMat(estGyro));
		double[] _omega_n_ie = new double[] { earthAngularRate * Math.cos(lat), 0, -earthAngularRate * Math.sin(lat) };
		SimpleMatrix omega_n_ie = new SimpleMatrix(Matrix.getSkewSymMat(_omega_n_ie));
		double[] _omega_n_en = new double[] { vel[1] / (Rn + alt), -vel[0] / (Rm + alt),
				-vel[1] * Math.tan(lat) / (Rn + alt) };
		SimpleMatrix omega_n_en = new SimpleMatrix(Matrix.getSkewSymMat(_omega_n_en));
		SimpleMatrix I = SimpleMatrix.identity(3);
		SimpleMatrix newDcm = (oldDcm.mult(I.plus(omega_b_ib.scale(tau))))
				.minus((omega_n_ie.plus(omega_n_en)).mult(oldDcm.scale(tau)));
		newDcm = Rotation.reorthonormDcm(newDcm);

		// Specific-Force Frame Transformation
		SimpleMatrix f_b_ib = new SimpleMatrix(3, 1, true, estAcc);
		SimpleMatrix f_n_ib = (oldDcm.plus(newDcm)).scale(0.5).mult(f_b_ib);

		// Velocity Update
		SimpleMatrix oldVel = new SimpleMatrix(3, 1, true, vel);
		SimpleMatrix g_n_b = new SimpleMatrix(3, 1, true, new double[] { 0, 0, LatLonUtil.getGravity(lat, alt) });
		SimpleMatrix newVel = oldVel
				.plus((f_n_ib.plus(g_n_b).minus((omega_n_en.plus((omega_n_ie.scale(2)))).mult(oldVel))).scale(tau));

		// Position Update
		double newAlt = alt - ((tau / 2) * (oldVel.get(2) + newVel.get(2)));
		double newLat = lat + ((tau / 2) * ((oldVel.get(0) / (Rm + alt)) + (newVel.get(0) / (Rm + newAlt))));
		double newRe = LatLonUtil.getMeridianEarthRadius(newLat);
		double newLon = lon + ((tau / 2) * ((oldVel.get(1) / ((Rn + alt) * Math.cos(lat)))
				+ (newVel.get(1) / ((newRe + newAlt) * Math.cos(newLat)))));

		// Receiver Clock update
		double oldRxClkDrift = X.getRxClk()[1];
		double oldRxClkOff = X.getRxClk()[0];
		double newRxClkDrift = oldRxClkDrift;
		double newRxClkOff = oldRxClkOff + (newRxClkDrift * tau);

		X.setAccBias(accBias);
		X.setGyroBias(gyroBias);
		X.setDcm(newDcm);
		X.setV(new double[] { newVel.get(0), newVel.get(1), newVel.get(2) });
		X.setP(new double[] { newLat, newLon, newAlt });
		X.setRxClk(new double[] { newRxClkOff, newRxClkDrift });

		return new double[][] { estAcc, estGyro };
	}

	private static SimpleMatrix predictErrorState(State X, SimpleMatrix P, SimpleMatrix phi, SimpleMatrix Qk) {
		P = (phi.mult(P).mult(phi.transpose())).plus(Qk);
		return P;
	}

	private static SimpleMatrix[] getDiscreteParams(State X, double[] estAcc, double[] estGyro, double tau,
			double[] q) {

		double lat = X.getP()[0];
		double lon = X.getP()[1];
		double alt = X.getP()[2];
		double vn = X.getV()[0];
		double ve = X.getV()[1];
		double vd = X.getV()[2];
		SimpleMatrix dcm = new SimpleMatrix(X.getDcm());
		double Rn = LatLonUtil.getNormalEarthRadius(lat);
		double Rm = LatLonUtil.getMeridianEarthRadius(lat);
		double omega_ie = LatLonUtil.omega_ie;
		double slat = Math.sin(lat);
		double clat = Math.cos(lat);
		double tlat = Math.tan(lat);
		// Instead of Re has used original denominator in rhos
		double rhoE = -vn / (Rm + alt);
		double rhoN = ve / (Rn + alt);
		double rhoD = -ve * tlat / (Rn + alt);
		double OmegaN = omega_ie * clat;
		double OmegaD = -omega_ie * slat;
		double omega_N = OmegaN + rhoN;
		double omega_E = rhoE;
		double omega_D = OmegaD + rhoD;
		// Unsure whether to use Rn+alt or Rm+alt for Re
		double Re = Rn + alt;
		double kD = vd / Re;
		double g = LatLonUtil.getGravity(lat, alt);
		double F63 = Math.pow(rhoN, 2) + Math.pow(rhoE, 2) - (2 * g / Re);

		SimpleMatrix Fpp = new SimpleMatrix(new double[][] { { 0, 0, -vn / Math.pow(Rm + alt, 2) },
				{ ve * slat / ((Rn + alt) * Math.pow(clat, 2)), 0, -ve / (Math.pow(Rn + alt, 2) * clat) },
				{ 0, 0, 0 } });
		SimpleMatrix Fpv = new SimpleMatrix(
				new double[][] { { 1 / (Rm + alt), 0, 0 }, { 0, 1 / ((Rn + alt) * clat), 0 }, { 0, 0, -1 } });
		SimpleMatrix Fpj = new SimpleMatrix(new double[3][3]);
		SimpleMatrix Fjp = new SimpleMatrix(new double[][] { { omega_ie * slat, 0, ve / Math.pow(Rn + alt, 2) },
				{ 0, 0, -vn / Math.pow(Rm + alt, 2) }, { (omega_ie * clat) + (ve / ((Rn + alt) * Math.pow(clat, 2))), 0,
						-(ve * tlat) / Math.pow(Rn + alt, 2) } });
		SimpleMatrix Fjv = new SimpleMatrix(
				new double[][] { { 0, -1 / (Rn + alt), 0 }, { 1 / (Rm + alt), 0, 0 }, { 0, tlat / (Rn + alt), 0 } });
		SimpleMatrix Fjj = new SimpleMatrix(Matrix.getSkewSymMat(new double[] { -omega_N, -omega_E, -omega_D }));
		SimpleMatrix Fvp = new SimpleMatrix(new double[][] {
				{ (-2 * OmegaN * ve) - (rhoN * ve / Math.pow(clat, 2)), 0, (rhoE * kD) - (rhoN * rhoD) },
				{ (2 * ((OmegaN * vn) + (OmegaD * vd))) + (rhoN * vn / Math.pow(clat, 2)), 0,
						(-rhoE * rhoD) - (kD * rhoN) },
				{ -2 * ve * OmegaD, 0, F63 } });
		SimpleMatrix Fvv = new SimpleMatrix(new double[][] { { kD, 2 * omega_D, -rhoE },
				{ -(omega_D + OmegaD), kD - (rhoE * tlat), omega_N + OmegaN }, { 2 * rhoE, -2 * omega_N, 0 } });
		SimpleMatrix Fvj = new SimpleMatrix(Matrix.getSkewSymMat(estAcc, true));
		SimpleMatrix zero = new SimpleMatrix(3, 3);
		SimpleMatrix identity = SimpleMatrix.identity(3);
		SimpleMatrix Fva = identity;
		SimpleMatrix Fjg = identity;
		SimpleMatrix Faa = identity.scale(-1 / ImuDataSheets.Pixel4.accBiasCorrTime);
		SimpleMatrix Fgg = identity.scale(-1 / ImuDataSheets.Pixel4.gyroBiasCorrTime);
		SimpleMatrix F = (Fpp.concatColumns(Fpv, Fpj, zero, zero))
				.concatRows(Fvp.concatColumns(Fvv, Fvj, dcm.scale(-1).mult(Fva), zero))
				.concatRows(Fjp.concatColumns(Fjv, Fjj, zero, dcm.mult(Fjg)))
				.concatRows(zero.concatColumns(zero, zero, Faa, zero))
				.concatRows(zero.concatColumns(zero, zero, zero, Fgg));
		F = F.concatColumns(new SimpleMatrix(15, 2));
		F = F.concatRows(new SimpleMatrix(2, 17));
		F.set(15, 16, 1);
		// First order approx of phi
		SimpleMatrix phi = SimpleMatrix.identity(17).plus(F.scale(tau));
		SimpleMatrix G = (zero.concatColumns(zero, zero, zero))
				.concatRows(dcm.scale(-1).concatColumns(zero, zero, zero))
				.concatRows(zero.concatColumns(dcm, zero, zero)).concatRows(zero.concatColumns(zero, identity, zero))
				.concatRows(zero.concatColumns(zero, zero, identity));
		G = G.concatColumns(new SimpleMatrix(15, 2));
		G = G.concatRows(new SimpleMatrix(2, 14));
		G.set(15, 12, 1);
		G.set(16, 13, 1);
		SimpleMatrix Q = new SimpleMatrix(14, 14);
		IntStream.range(0, 14).forEach(i -> Q.set(i, i, q[i]));
		// First order approx of Qk
		SimpleMatrix Qk = G.mult(Q).mult(G.transpose()).scale(tau);
		return new SimpleMatrix[] { phi, Qk };

	}

	public static SimpleMatrix update(State X, SimpleMatrix P, long utcTimeMilli, IMUsensor magData,
			ArrayList<Satellite> SV, boolean onlyMagneto) {

		double[] llh = X.getP();
		double[] vel = X.getV();
		double rxClkOff = X.getRxClk()[0];
		double rxClkDrift = X.getRxClk()[1];
		SimpleMatrix dcm = X.getDcm();
		double[] pos_ecef = LatLonUtil.lla2ecef(llh, false);
		double[] vel_ecef = LatLonUtil.ned2ecef(vel, pos_ecef, false);

		GeomagneticField gmf = new GeomagneticField((float) Math.toDegrees(llh[0]), (float) Math.toDegrees(llh[1]),
				(float) llh[2], utcTimeMilli);
		// True Magnetic Strength
		final double[] mag = new double[] { gmf.getX(), gmf.getY(), gmf.getZ() };
		// Convert from nanotesla to microtesla
		for (int i = 0; i < 3; i++) {
			mag[i] *= 1e-3;
		}

		// Subtract hard bias from magnetometer measurements, estimated magnetic
		// strength and convert from body frame to local-nav frame
		final SimpleMatrix mag_hat = dcm.mult(new SimpleMatrix(3, 1, true,
				IntStream.range(0, 3).mapToDouble(i -> magData.getVal()[i] - magData.getBias()[i]).toArray()));
		double[] zm = IntStream.range(0, 3).mapToDouble(i -> mag[i] - mag_hat.get(i)).toArray();
		double[][] Mag = Matrix.getSkewSymMat(mag, true);
		double[][] magNoise = new double[3][3];
		IntStream.range(0, 3).forEach(i -> magNoise[i][i] = ImuDataSheets.Pixel4.magnetoNoise);
		// Magnetometer noise covariance
		SimpleMatrix Rm = dcm.mult(new SimpleMatrix(magNoise)).mult(dcm.transpose());
		SimpleMatrix H = null;
		SimpleMatrix Z = null;
		SimpleMatrix R = null;
		if (onlyMagneto) {
			// Delta Magnetic
			Z = new SimpleMatrix(3, 1, true, zm);
			H = new SimpleMatrix(3, 17);
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					H.set(i, 6 + j, Mag[i][j]);
				}
			}
			R = Rm;
		} else {

			int n = SV.size();
			double[] z = new double[(2 * n) + 3];
			double[][] a = new double[n][3];
			R = new SimpleMatrix((2 * n) + 3, (2 * n) + 3);
			for (int i = 0; i < n; i++) {
				Satellite sat = SV.get(i);
				double pseudorange = sat.getPseudorange();
				double rangeRate = sat.getRangeRate();
				// Its not really a ECI, therefore don't get confused
				double[] satEcef = sat.getSatEci();
				// Approx Geometric Range
				double approxGR = Math.sqrt(IntStream.range(0, 3).mapToDouble(j -> satEcef[j] - pos_ecef[j])
						.map(j -> Math.pow(j, 2)).reduce((j, k) -> j + k).getAsDouble());
				final int index = i;
				IntStream.range(0, 3).forEach(j -> a[index][j] = (satEcef[j] - pos_ecef[j]) / approxGR);
				SimpleMatrix satVel = new SimpleMatrix(3, 1, true, sat.getSatVel());
				SimpleMatrix A = new SimpleMatrix(1, 3, true, a[index]);
				/*
				 * Observable derived from doppler and satellite velocity, refer Kaplan and
				 * personal notes
				 */
				double dopplerDerivedObs = rangeRate - A.mult(satVel).get(0);
				// Approx Pseudorange Range
				double approxPR = approxGR + rxClkOff;
				// Approx doppler derived observable
				double approxDopplerDerivedObs = -A.mult(new SimpleMatrix(3, 1, true, vel_ecef)).get(0) + rxClkDrift;
				// Delta PR
				z[i] = pseudorange - approxPR;
				// Delta doppler derived observable
				z[i + n] = dopplerDerivedObs - approxDopplerDerivedObs;

				R.set(i, i, Math.pow(sat.getReceivedSvTimeUncertaintyNanos() * SpeedofLight * 1e-9, 2));
				R.set(i + n, i + n, Math.pow(sat.getPseudorangeRateUncertaintyMetersPerSecond(), 2));
			}
			IntStream.range(0, 3).forEach(i -> z[(2 * n) + i] = zm[i]);

			H = new SimpleMatrix((2 * n) + 3, 17);
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < 3; j++) {
					H.set(i, j, -a[i][j]);
					H.set(i + n, j + 3, -a[i][j]);
				}
				H.set(i, 15, 1);
				H.set(i + n, 16, 1);
			}
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					H.set(i + (2 * n), 6 + j, Mag[i][j]);
					R.set((2 * n) + i, (2 * n) + j, Rm.get(i, j));

				}

			}

			Z = new SimpleMatrix(z.length, 1, true, z);
		}
		SimpleMatrix Ht = H.transpose();
		SimpleMatrix K = null;
		// Kalman Gain

		K = P.mult(Ht).mult(((H.mult(P).mult(Ht)).plus(R)).invert());

		// Posterior State Estimate
		// As prior deltaX is zero, there is no 'ze'
		SimpleMatrix deltaX = K.mult(Z);
		SimpleMatrix KH = K.mult(H);
		SimpleMatrix I = SimpleMatrix.identity(KH.numRows());
		// Posterior Estimate Error Joseph Form to ensure Positive Definiteness P =
		// (I-KH)P(I-KH)' + KRK'
		P = ((I.minus(KH)).mult(P).mult((I.minus(KH)).transpose())).plus(K.mult(R).mult(K.transpose()));

		// Update Total State
		pos_ecef[0] += deltaX.get(0);
		pos_ecef[1] += deltaX.get(1);
		pos_ecef[2] += deltaX.get(2);
		double[] _llh = LatLonUtil.ecef2lla(pos_ecef);
		IntStream.range(0, 2).forEach(i -> _llh[i] = Math.toRadians(_llh[i]));
		X.setP(_llh);
		vel_ecef[0] += deltaX.get(3);
		vel_ecef[1] += deltaX.get(4);
		vel_ecef[2] += deltaX.get(5);
		X.setV(LatLonUtil.ecef2ned(vel_ecef, pos_ecef, false));
		SimpleMatrix updateDcm = new SimpleMatrix(
				Matrix.getSkewSymMat(new double[] { deltaX.get(6), deltaX.get(7), deltaX.get(8) }));
		dcm = (SimpleMatrix.identity(3).minus(updateDcm)).mult(dcm);
		dcm = Rotation.reorthonormDcm(dcm);
		X.setDcm(dcm);
		X.setAccBias(IntStream.range(0, 3).mapToDouble(i -> X.getAccBias()[i] + deltaX.get(9 + i)).toArray());
		X.setGyroBias(IntStream.range(0, 3).mapToDouble(i -> X.getGyroBias()[i] + deltaX.get(12 + i)).toArray());
		X.setRxClk(IntStream.range(0, 2).mapToDouble(i -> X.getRxClk()[i] + deltaX.get(15 + i)).toArray());
		return P;
	}
}
