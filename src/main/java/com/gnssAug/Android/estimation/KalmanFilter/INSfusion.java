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
import com.gnssAug.Android.helper.Rotation;
import com.gnssAug.Android.models.IMUsensor;
import com.gnssAug.Android.models.Satellite;
import com.gnssAug.Android.utility.LatLonUtil;
import com.gnssAug.Android.utility.Matrix;

public class INSfusion {

	public static void process(TreeMap<Long, HashMap<AndroidSensor, IMUsensor>> imuMap,
			ArrayList<ArrayList<Satellite>> SVlist, ArrayList<Long> timeList, double[][] dcm) throws Exception {

		// Convert DCM from Body Frame to ENU to Body frame to NED
		dcm = LatLonUtil.enu_ned_convert(dcm);
		double[] ecef0 = LinearLeastSquare.process(SVlist.get(0), true);
		double[] llh0 = LatLonUtil.ecef2lla(ecef0);
		IntStream.range(0, 2).forEach(i -> llh0[i] = Math.toRadians(llh0[i]));
		// Velocity in ENU frame, zero initially
		double[] vel0 = new double[3];
		double rxClkOff = ecef0[3];
		State x = new State(llh0[0], llh0[1], llh0[2], vel0[0], vel0[1], vel0[2], dcm, 0, 0, 0, 0, 0, 0, rxClkOff, 0);
		// Attitude std deviation is 20 degree, values mentioned below in covariance
		// matrix is in radians
		double accBiasCov = Math.pow(ImuDataSheets.Pixel4.accTurnOnBias, 2);
		double gyroBiasCov = Math.pow(ImuDataSheets.Pixel4.gyroTurnOnBias, 2);
		double[] p0 = new double[] { 100, 100, 100, 0.1, 0.1, 0.1, 0.1156, 0.1156, 0.1156, accBiasCov, accBiasCov,
				accBiasCov, gyroBiasCov, gyroBiasCov, gyroBiasCov, 5, 0.005 };
		SimpleMatrix p = new SimpleMatrix(17, 17);
		IntStream.range(0, 17).forEach(i -> p.set(i, i, p0[i]));
		// PSD of random walk - N
		double acc_SN = Math.pow(ImuDataSheets.Pixel4.vrw, 2);
		double gyro_SN = Math.pow(ImuDataSheets.Pixel4.arw, 2);
		// PSD of Bias Instability - B
		double acc_SB = 2 * Math.pow(ImuDataSheets.Pixel4.accInRunBias, 2) * Math.log(2)
				/ (Math.PI * Math.pow(0.4365, 2) * ImuDataSheets.Pixel4.accBiasCorrTime);
		double gyro_SB = 2 * Math.pow(ImuDataSheets.Pixel4.gyroInRunBias, 2) * Math.log(2)
				/ (Math.PI * Math.pow(0.4365, 2) * ImuDataSheets.Pixel4.gyroBiasCorrTime);
		// Receiver Clock phase and frequency PSD
		double sf = ClockAllanVar.TCXO_low_qaulity.sf;
		double sg = ClockAllanVar.TCXO_low_qaulity.sg;
		double[] q = new double[] { acc_SN, gyro_SN, acc_SB, gyro_SB, sf, sg };
		int n = SVlist.size();
		if (n != timeList.size()) {
			System.err.println("FATAL ERROR in INSfusion.java");
			throw new Exception("size of SVlist and timelist does not match");
		}

		int i = 0;
		Iterator<Entry<Long, HashMap<AndroidSensor, IMUsensor>>> iterator = imuMap.entrySet().iterator();
		double prevTime = iterator.next().getKey();
		do {

			while (iterator.hasNext()) {
				Entry<Long, HashMap<AndroidSensor, IMUsensor>> entry = iterator.next();
				long time = entry.getKey();
				HashMap<AndroidSensor, IMUsensor> imuSensor = entry.getValue();
				double tau = time - prevTime;
				// Update Total State and get estimated INS observables - {acc,gyro}
				double[][] estInsObs = updateTotalState(x, imuSensor, tau);
				prevTime = time;
				/*
				 * Get transition matrix Phi and covariance matrix Qk, for discrete state space
				 * filtering
				 */
				SimpleMatrix[] discParam = getDiscreteParams(x, estInsObs[0], estInsObs[1], tau, q);
				updateErrorState(x, p, discParam[0], discParam[1]);
				if (time == timeList.get(i)) {

					i++;
				}

			}
		} while (i < n);

	}

	private static double[][] updateTotalState(State x, HashMap<AndroidSensor, IMUsensor> imuSensor, double tau) {

		// Update Bias state
		double[] accBias = Arrays.stream(x.getAccBias())
				.map(i -> i * Math.exp(-tau / ImuDataSheets.Pixel4.accBiasCorrTime)).toArray();
		double[] gyroBias = Arrays.stream(x.getGyroBias())
				.map(i -> i * Math.exp(-tau / ImuDataSheets.Pixel4.gyroBiasCorrTime)).toArray();

		double[] obsAcc = imuSensor.get(AndroidSensor.Accelerometer).getVal();
		double[] estAcc = IntStream.range(0, 3).mapToDouble(j -> obsAcc[j] - accBias[j]).toArray();
		double[] obsGyro = imuSensor.get(AndroidSensor.Gyroscope).getVal();
		double[] estGyro = IntStream.range(0, 3).mapToDouble(j -> obsGyro[j] - gyroBias[j]).toArray();

		double lat = x.getP()[0];
		double lon = x.getP()[1];
		double alt = x.getP()[2];
		double[] vel = x.getV();
		double Rn = LatLonUtil.getNormalEarthRadius(lat);
		double Rm = LatLonUtil.getMeridianEarthRadius(lat);
		double earthAngularRate = LatLonUtil.omega_ie;

		// Attitude update
		SimpleMatrix oldDcm = new SimpleMatrix(x.getDcm());
		SimpleMatrix omega_b_ib = new SimpleMatrix(Matrix.getSkewSymMat(estGyro));
		double[] _omega_n_ie = new double[] { earthAngularRate * Math.cos(lat), 0, -earthAngularRate * Math.sin(lat) };
		SimpleMatrix omega_n_ie = new SimpleMatrix(Matrix.getSkewSymMat(_omega_n_ie));
		double[] _omega_n_en = new double[] { vel[0], -vel[1], -vel[0] * Math.tan(lat) };
		Arrays.stream(_omega_n_en).forEach(i -> i = i / (Rn + alt));
		SimpleMatrix omega_n_en = new SimpleMatrix(Matrix.getSkewSymMat(_omega_n_en));
		SimpleMatrix I = SimpleMatrix.identity(3);
		SimpleMatrix newDcm = (oldDcm.mult(I.plus(omega_b_ib.scale(tau))))
				.minus((omega_n_ie.plus(omega_n_en)).mult(oldDcm.scale(tau)));
		newDcm = Rotation.reorthonormDcm(newDcm);

		// Specific-Force Frame Transformation
		SimpleMatrix f_b_ib = new SimpleMatrix(3, 1, true, estAcc);
		SimpleMatrix f_n_ib = f_b_ib.scale(0.5).mult((oldDcm.plus(newDcm)));

		// Velocity Update
		SimpleMatrix oldVel = new SimpleMatrix(3, 1, true, vel);
		SimpleMatrix g_n_b = new SimpleMatrix(3, 1, true, new double[] { 0, 0, LatLonUtil.getGravity(lat, alt) });
		SimpleMatrix newVel = oldVel
				.plus((f_n_ib.plus(g_n_b).minus(oldVel.mult(omega_n_en.plus((omega_n_ie.scale(2)))))).scale(tau));

		// Position Update
		double newAlt = alt - ((tau / 2) * (oldVel.get(2) + newVel.get(2)));
		double newLat = lat + ((tau / 2) * ((oldVel.get(0) / (Rm + alt)) + (newVel.get(0) / (Rm + newAlt))));
		double newRe = LatLonUtil.getMeridianEarthRadius(newLat);
		double newLon = lon + ((tau / 2) * ((oldVel.get(1) / ((Rn + alt) * Math.cos(lat)))
				+ (newVel.get(1) / ((newRe + newAlt) * Math.cos(newLat)))));

		// Receiver Clock update
		double oldRxClkDrift = x.getRxClk()[1];
		double oldRxClkOff = x.getRxClk()[0];
		double newRxClkDrift = oldRxClkDrift;
		double newRxClkOff = oldRxClkOff + (newRxClkDrift * tau);

		x.setAccBias(accBias);
		x.setGyroBias(gyroBias);
		x.setDcm(Matrix.matrix2Array(newDcm));
		x.setV(new double[] { newVel.get(0), newVel.get(1), newVel.get(2) });
		x.setP(new double[] { newLat, newLon, newAlt });
		x.setRxClk(new double[] { newRxClkOff, newRxClkDrift });

		return new double[][] { estAcc, estGyro };
	}

	private static void updateErrorState(State x, SimpleMatrix p, SimpleMatrix phi, SimpleMatrix Qk) {
		p = phi.mult(p).mult(phi.transpose()).plus(Qk);
	}

	private static SimpleMatrix[] getDiscreteParams(State x, double[] estAcc, double[] estGyro, double tau,
			double[] q) {

		double lat = x.getP()[0];
		double lon = x.getP()[1];
		double alt = x.getP()[2];
		double vn = x.getV()[0];
		double ve = x.getV()[1];
		double vd = x.getV()[2];
		SimpleMatrix dcm = new SimpleMatrix(x.getDcm());
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
		F.concatColumns(new SimpleMatrix(15, 2));
		F.concatRows(new SimpleMatrix(2, 17));
		F.set(15, 16, 1);
		// First order approx of phi
		SimpleMatrix phi = SimpleMatrix.identity(17).plus(F.scale(tau));
		SimpleMatrix G = (zero.concatColumns(zero, zero, zero))
				.concatRows(dcm.scale(-1).concatColumns(zero, zero, zero))
				.concatRows(zero.concatColumns(dcm, zero, zero)).concatRows(zero.concatColumns(zero, identity, zero))
				.concatRows(zero.concatColumns(zero, zero, identity));
		G.concatColumns(new SimpleMatrix(15, 2));
		G.concatRows(new SimpleMatrix(2, 14));
		G.set(15, 12, 1);
		G.set(16, 13, 1);
		SimpleMatrix Q = new SimpleMatrix(6, 1, true, q);
		// First order approx of Qk
		SimpleMatrix Qk = G.mult(Q).mult(G.transpose()).scale(tau);
		return new SimpleMatrix[] { phi, Qk };

	}
}
