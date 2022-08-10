package com.gnssAug.Android;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.TimeZone;
import java.util.TreeMap;
import java.util.stream.IntStream;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.orekit.data.DataProvidersManager;
import org.orekit.data.DirectoryCrawler;
import org.orekit.forces.gravity.potential.GravityFieldFactory;
import org.orekit.forces.gravity.potential.NormalizedSphericalHarmonicsProvider;
import org.orekit.frames.Frame;
import org.orekit.frames.FramesFactory;
import org.orekit.frames.ITRFVersion;
import org.orekit.models.earth.Geoid;
import org.orekit.models.earth.ReferenceEllipsoid;
import org.orekit.utils.IERSConventions;

import com.gnssAug.Android.constants.AndroidSensor;
import com.gnssAug.Android.estimation.LinearLeastSquare;
import com.gnssAug.Android.estimation.KalmanFilter.EKF;
import com.gnssAug.Android.estimation.KalmanFilter.INSfusion;
import com.gnssAug.Android.estimation.KalmanFilter.Models.Flag;
import com.gnssAug.Android.fileParser.DerivedCSV;
import com.gnssAug.Android.fileParser.GNSS_Log;
import com.gnssAug.Android.fileParser.GroundTruth;
import com.gnssAug.Android.helper.ComputeEleAzm;
import com.gnssAug.Android.helper.INS.IMUconfigure;
import com.gnssAug.Android.helper.INS.StateInitialization;
import com.gnssAug.Android.models.Derived;
import com.gnssAug.Android.models.GNSSLog;
import com.gnssAug.Android.models.IMUsensor;
import com.gnssAug.Android.models.Satellite;
import com.gnssAug.fileParser.Bias;
import com.gnssAug.fileParser.Clock;
import com.gnssAug.fileParser.Orbit;
import com.gnssAug.utility.Analyzer;
import com.gnssAug.utility.GraphPlotter;
import com.gnssAug.utility.LatLonUtil;
import com.gnssAug.utility.MathUtil;
import com.gnssAug.utility.Time;

public class Android {
	public static void posEstimate(boolean doPosErrPlot, double cutOffAng, double snrMask, int estimatorType,
			String[] obsvCodeList, String derived_csv_path, String gnss_log_path, String GTcsv, String bias_path,
			String clock_path, String orbit_path, boolean useIGS, boolean checkOutlier) {
		try {

			TimeZone.setDefault(TimeZone.getTimeZone("UTC"));
			HashMap<String, ArrayList<HashMap<String, Double>>> ErrMap = new HashMap<String, ArrayList<HashMap<String, Double>>>();

			ArrayList<Long> timeList = new ArrayList<Long>();
			ArrayList<double[]> trueLLHlist = new ArrayList<double[]>();
			ArrayList<double[]> trueEcefList = new ArrayList<double[]>();

			HashMap<Long, double[]> err = new HashMap<Long, double[]>();
			TreeMap<Long, ArrayList<Satellite>> SatMap = new TreeMap<Long, ArrayList<Satellite>>();
			HashMap<String, ArrayList<double[]>> estPosMap = new HashMap<String, ArrayList<double[]>>();
			HashMap<String, ArrayList<double[]>> estVelMap = new HashMap<String, ArrayList<double[]>>();
			ArrayList<Double> postUnitWeightList = new ArrayList<Double>();
			Bias bias = null;
			Orbit orbit = null;
			Clock clock = null;

			String path = "C:\\Users\\Naman\\Desktop\\rinex_parse_files\\google2\\test2";
			File output = new File(path + ".txt");
			PrintStream stream;
			stream = new PrintStream(output);
			System.setOut(stream);

			ArrayList<double[]> rxLLH = GroundTruth.processCSV(GTcsv);
			HashMap<Long, HashMap<String, HashMap<Integer, Derived>>> derivedMap = null;

			derivedMap = DerivedCSV.processCSV(derived_csv_path);
			GNSS_Log.process(gnss_log_path);
			TreeMap<Long, HashMap<String, ArrayList<GNSSLog>>> gnssLogMaps = GNSS_Log.getGnssLogMaps();
			ArrayList<IMUsensor> imuList = GNSS_Log.getImuList();

			if (useIGS) {

				orbit = new Orbit(orbit_path);
				bias = new Bias(bias_path);
				clock = new Clock(clock_path, bias);

			}

			int gtIndex = 0;
			for (long tRxMilli : gnssLogMaps.keySet()) {
				if (gtIndex >= rxLLH.size()) {
					break;
				}
				HashMap<String, ArrayList<GNSSLog>> gnssLogMap = gnssLogMaps.get(tRxMilli);
				GNSSLog entry = ((ArrayList<GNSSLog>) gnssLogMap.values().toArray()[0]).get(0);
				double tRx = entry.gettRx();
				int weekNo = entry.getWeekNo();
				if (tRxMilli != (rxLLH.get(gtIndex)[0] * 1000) || weekNo != rxLLH.get(gtIndex)[1]) {
					System.err.println("FATAL ERROR - GT timestamp does not match");
					continue;
				}
				double[] trueUserLLH = new double[] { rxLLH.get(gtIndex)[2], rxLLH.get(gtIndex)[3],
						rxLLH.get(gtIndex)[4] };
				double trueVelRms = rxLLH.get(gtIndex)[5];
				gtIndex++;

				Calendar time = Time.getDate(tRx, weekNo, 0);
				ArrayList<Satellite> satList = SingleFreq.process(tRx, derivedMap, gnssLogMap, time, obsvCodeList,
						weekNo, clock, orbit, useIGS);
				if (satList.size() < 4) {
					System.err.println("Less than 4 satellites");
					continue;
				}
				double[] userEcef = LinearLeastSquare.process(satList, false);
				satList.stream().forEach(i -> i.setElevAzm(ComputeEleAzm.computeEleAzm(userEcef, i.getSatEci())));
				filterSat(satList, cutOffAng, snrMask);
				double[] truePosEcef = LatLonUtil
						.lla2ecef(new double[] { trueUserLLH[0], trueUserLLH[1], trueUserLLH[2] - 61 }, true);
				if (checkOutlier) {
					checkOutlier(satList, truePosEcef);
				}
				if (satList.size() < 4) {
					System.err.println("Less than 4 satellites");
					continue;
				}
				double[] estEcefClk = null;
				switch (estimatorType) {
				case 1:
					// Implement LS method
					estEcefClk = LinearLeastSquare.process(satList, false);
					estPosMap.computeIfAbsent("LS", k -> new ArrayList<double[]>()).add(estEcefClk);
					break;
				case 2:
					// Implement WLS method
					estEcefClk = LinearLeastSquare.process(satList, true);
					estPosMap.computeIfAbsent("WLS", k -> new ArrayList<double[]>()).add(estEcefClk);
					double[] estVel = LinearLeastSquare.getEstVel(satList, estEcefClk);
					estVelMap.computeIfAbsent("WLS", k -> new ArrayList<double[]>()).add(estVel);

					// Implement WLS method
//					estEcefClk = LinearLeastSquare.process(satList, true, true);
//					estPosMap.computeIfAbsent("WLS_QualityControl", k -> new ArrayList<double[]>()).add(estEcefClk);
//					double postUnitW = LinearLeastSquare.getPostUnitW();
//					postUnitWeightList.add(Math.sqrt(postUnitW));
					break;
				case 3:
//					// Implement LS method
//					estEcefClk = LinearLeastSquare.process(satList, false);
//					estPosMap.computeIfAbsent("LS", k -> new ArrayList<double[]>()).add(estEcefClk);
					// Implement WLS method
					estEcefClk = LinearLeastSquare.process(satList, true);
					estPosMap.computeIfAbsent("WLS", k -> new ArrayList<double[]>()).add(estEcefClk);
					estVel = LinearLeastSquare.getEstVel(satList, estEcefClk);
					estVelMap.computeIfAbsent("WLS", k -> new ArrayList<double[]>()).add(estVel);

				}
				SatMap.put(tRxMilli, satList);
				trueLLHlist.add(trueUserLLH);
				trueEcefList.add(truePosEcef);
				timeList.add(tRxMilli);
			}
			// Get True Velocity
			TreeMap<Long, double[]> trueVelEcef = Analyzer.getVel(trueEcefList, timeList);
			if (estimatorType == 4) {
				TreeMap<Long, HashMap<AndroidSensor, IMUsensor>> imuMap = IMUconfigure.configure(timeList.get(0), 100,
						imuList);
				for (Map.Entry<Long, HashMap<AndroidSensor, IMUsensor>> entry : imuMap.entrySet()) {
					if (entry.getValue().size() != 3) {
						System.err.println("FATAL ERROR in IMU map");
						throw new Exception("Erroneous imu sampling rate");
					}
				}

				// Body Frame to ENU
				double[][] dcm = StateInitialization.initialize(imuMap, SatMap);
				TreeMap<Long, double[]> ecefMap = INSfusion.process(imuMap, SatMap, timeList, dcm);
				// GraphPlotter.graphGnssIns(ecefMap, trueECEFlist, timeList);
				int n = timeList.size();
				for (int i = 0; i < n; i++) {
					long time = timeList.get(i);
					double[] estEcef = ecefMap.get(time);
					estPosMap.computeIfAbsent("GNSS/INS fusion", k -> new ArrayList<double[]>()).add(estEcef);

				}

			}
			if (estimatorType == 5 || estimatorType == 3) {

				EKF ekf = new EKF();
//				 Implement EKF based on receiver’s position and clock offset errors as a
//				 random walk process
				TreeMap<Long, double[]> estStateMap_pos = ekf.process(SatMap, timeList, Flag.POSITION, false);
//				 Implement EKF based on receiver’s velocity and clock drift errors as a random
//				 walk process
				TreeMap<Long, double[]> estStateMap_vel = ekf.process(SatMap, timeList, Flag.VELOCITY, false);
//				 Implement EKF based on receiver’s velocity and clock drift errors as a random
//				 walk process along with doppler updates
				TreeMap<Long, double[]> estStateMap_vel_doppler = ekf.process(SatMap, timeList, Flag.VELOCITY, true);
				int n = timeList.size();
				for (int i = 0; i < n; i++) {
					long time = timeList.get(i);
					double[] estPos = estStateMap_pos.get(time);
					estPosMap.computeIfAbsent("EKF - pos. random walk", k -> new ArrayList<double[]>()).add(estPos);
					double[] estState = estStateMap_vel.get(time);
					estPos = null;
					double[] estVel = null;
					if (estState != null) {
						estPos = new double[] { estState[0], estState[1], estState[2] };
						estVel = new double[] { estState[3], estState[4], estState[5] };
					}
					estPosMap.computeIfAbsent("EKF - vel. random walk", k -> new ArrayList<double[]>()).add(estPos);
					estVelMap.computeIfAbsent("EKF - vel. random walk", k -> new ArrayList<double[]>()).add(estVel);
					estState = estStateMap_vel_doppler.get(time);
					estPos = null;
					estVel = null;
					if (estState != null) {
						estPos = new double[] { estState[0], estState[1], estState[2] };
						estVel = new double[] { estState[3], estState[4], estState[5] };
					}
					estPosMap.computeIfAbsent("EKF - vel. random walk + doppler", k -> new ArrayList<double[]>())
							.add(estPos);
					estVelMap.computeIfAbsent("EKF - vel. random walk + doppler", k -> new ArrayList<double[]>())
							.add(estVel);

				}

			}

			if (estimatorType == 6) {
				TreeMap<Long, HashMap<AndroidSensor, IMUsensor>> imuMap = null;
//				TreeMap<Long, HashMap<AndroidSensor, IMUsensor>> imuMap = IMUconfigure.configure(timeList.get(0), 100,
//						imuList);
				Analyzer.process(SatMap, imuMap, trueEcefList, trueVelEcef);

			}

			// Calculate Accuracy Metrics
			HashMap<String, ArrayList<double[]>> GraphPosMap = new HashMap<String, ArrayList<double[]>>();
			HashMap<String, ArrayList<double[]>> GraphVelMap = new HashMap<String, ArrayList<double[]>>();

			for (String key : estPosMap.keySet()) {
				ArrayList<Double>[] posErrList = new ArrayList[6];

				IntStream.range(0, 6).forEach(i -> posErrList[i] = new ArrayList<Double>());

				ArrayList<double[]> estPosList = estPosMap.get(key);

				int n = estPosList.size();
				if (n != trueEcefList.size()) {
					System.err.println("FATAL ERROR: EST and TRUE ecef list size does not match ");
					throw new Exception("FATAL ERROR: EST and TRUE ecef list size does not match ");
				}
				ArrayList<double[]> enuPosList = new ArrayList<double[]>();

				for (int i = 0; i < n; i++) {
					double[] estEcef = estPosList.get(i);
					if (estEcef == null) {
						continue;
					}
					double[] enu = LatLonUtil.ecef2enu(estEcef, trueEcefList.get(i), true);
					double[] estLLH = LatLonUtil.ecef2lla(estEcef);
					// Great Circle Distance
					double gcErr = LatLonUtil.getHaversineDistance(estLLH, trueLLHlist.get(i));
					enuPosList.add(enu);
					// error in East direction
					posErrList[0].add(Math.sqrt(enu[0] * enu[0]));
					// error in North direction
					posErrList[1].add(Math.sqrt(enu[1] * enu[1]));
					// error in Up direction
					posErrList[2].add(Math.sqrt(enu[2] * enu[2]));
					// 3d error
					posErrList[3].add(Math.sqrt(Arrays.stream(enu).map(j -> j * j).sum()));
					// 2d error
					posErrList[4].add(Math.sqrt((enu[0] * enu[0]) + (enu[1] * enu[1])));
					// Haversine Distance
					posErrList[5].add(gcErr);
				}

				GraphPosMap.put(key, enuPosList);

				// RMSE
				System.out.println("\n" + key);
				System.out.println("Position RMS - ");
				System.out.println(" E - " + RMS(posErrList[0]));
				System.out.println(" N - " + RMS(posErrList[1]));
				System.out.println(" U - " + RMS(posErrList[2]));
				System.out.println(" 3d Error - " + RMS(posErrList[3]));
				System.out.println(" 2d Error - " + RMS(posErrList[4]));
				System.out.println(" Haversine Distance - " + RMS(posErrList[5]));

				// 95th Percentile

				IntStream.range(0, 6).forEach(i -> Collections.sort(posErrList[i]));
				int q95 = (int) (n * 0.95);

				System.out.println("\n" + key + " 95%");
				System.out.println("RMS - ");
				System.out.println(" E - " + posErrList[0].get(q95));
				System.out.println(" N - " + posErrList[1].get(q95));
				System.out.println(" U - " + posErrList[2].get(q95));
				System.out.println(" 3d Error - " + posErrList[3].get(q95));
				System.out.println(" 2d Error - " + posErrList[4].get(q95));
				System.out.println(" Haversine distance - " + posErrList[5].get(q95));

			}
			for (String key : estVelMap.keySet()) {
				ArrayList<Double>[] velErrList = new ArrayList[6];
				ArrayList<double[]> estVelList = null;
				ArrayList<double[]> enuVelList = null;
				IntStream.range(0, 5).forEach(i -> velErrList[i] = new ArrayList<Double>());
				enuVelList = new ArrayList<double[]>();
				estVelList = estVelMap.get(key);
				int n = trueEcefList.size();
				for (int i = 0; i < n; i++) {
					double[] estVel = estVelList.get(i);
					long time = timeList.get(i);
					if (estVel == null || !trueVelEcef.containsKey(time)) {
						continue;
					}
					double[] trueVel = trueVelEcef.get(time);
					double[] velErr = IntStream.range(0, 3).mapToDouble(j -> estVel[j] - trueVel[j]).toArray();
					double[] enu = LatLonUtil.ecef2enu(velErr, trueEcefList.get(i), false);

					enuVelList.add(enu);
					// error in East direction
					velErrList[0].add(Math.sqrt(enu[0] * enu[0]));
					// error in North direction
					velErrList[1].add(Math.sqrt(enu[1] * enu[1]));
					// error in Up direction
					velErrList[2].add(Math.sqrt(enu[2] * enu[2]));
					// 3d error
					velErrList[3].add(Math.sqrt(Arrays.stream(enu).map(j -> j * j).sum()));
					// 2d error
					velErrList[4].add(Math.sqrt((enu[0] * enu[0]) + (enu[1] * enu[1])));

				}
				GraphVelMap.put(key, enuVelList);

				// RMSE
				System.out.println("\n" + key);
				System.out.println("Velocity RMS - ");
				System.out.println(" E - " + RMS(velErrList[0]));
				System.out.println(" N - " + RMS(velErrList[1]));
				System.out.println(" U - " + RMS(velErrList[2]));
				System.out.println(" 3d Error - " + RMS(velErrList[3]));
				System.out.println(" 2d Error - " + RMS(velErrList[4]));

				// 95th Percentile

//				IntStream.range(0, 5).forEach(i -> Collections.sort(velErrList[i]));
//				int q95 = (int) (n * 0.95);
//
//				System.out.println("\n" + key + " 95%");
//				System.out.println("RMS - ");
//				System.out.println(" E - " + velErrList[0].get(q95));
//				System.out.println(" N - " + velErrList[1].get(q95));
//				System.out.println(" U - " + velErrList[2].get(q95));
//				System.out.println(" 3d Error - " + velErrList[3].get(q95));
//				System.out.println(" 2d Error - " + velErrList[4].get(q95));

			}
			for (int i = 0; i < timeList.size(); i++) {
				timeList.set(i, (long) (timeList.get(i) * 1e-3));
			}
			// Plot Error Graphs
			GraphPlotter.graphENU(GraphPosMap, timeList, true);
			// Plot Error Graphs
			// GraphPlotter.graphENU(GraphVelMap, timeList, false);
			// Plot Posteriori Variance of Unit Weight
			GraphPlotter.graphPostUnitW(postUnitWeightList, timeList);

		} catch (Exception e) {
			// TODO: handle exception
			System.out.println(e.getMessage());
			e.printStackTrace();
		}
	}

	public static void filterSat(ArrayList<Satellite> satList, double cutOffAng, double snrMask) {
		if (cutOffAng >= 0) {
			satList.removeIf(i -> i.getElevAzm()[0] < Math.toRadians(cutOffAng));
		}
		if (snrMask >= 0) {
			satList.removeIf(i -> i.getCn0DbHz() < snrMask);
		}
	}

	public static void checkOutlier(ArrayList<Satellite> satList, double[] truePos) {
		final double SpeedofLight = 299792458;
		int n = satList.size();
		double max = Double.MIN_VALUE;
		int index = -1;
		HashSet<Satellite> indexSet = new HashSet<Satellite>();
		NormalDistribution norm = new NormalDistribution();
		for (int i = 0; i < n; i++) {
			Satellite sat = satList.get(i);
			double[] satPos = sat.getSatEci();
			double trueRange = MathUtil.getEuclidean(truePos, satPos);
			double pseudoRange = sat.getPseudorange();
			double err = Math.abs(pseudoRange - trueRange);
			double sigma = satList.get(i).getReceivedSvTimeUncertaintyNanos() * SpeedofLight * 1e-9;
			double test = err / sigma;
			if (1 - norm.cumulativeProbability(Math.abs(test)) < 0.001) {
				indexSet.add(sat);
				if (test > max) {
					max = test;
					index = i;
					sat.setFailsTest(true);
				}
			}
			sat.setTestStat(test);

		}
		if (index != -1) {
			satList.get(index).setOutlier(true);
			// satList.remove(index);
		}
		if (!indexSet.isEmpty()) {
			satList.removeAll(indexSet);
		}
	}

	public static double RMS(ArrayList<Double> list) {
		return Math.sqrt(list.stream().mapToDouble(x -> x * x).average().orElse(Double.NaN));
	}

	public static double MAE(ArrayList<Double> list) {
		return list.stream().mapToDouble(x -> x).average().orElse(Double.NaN);
	}

	public static Geoid buildGeoid() {
		// Semi-major axis or Equatorial radius
		final double ae = 6378137;
		// flattening
		final double f = 1 / 298.257223563;

		// Earth's rotation rate
		final double spin = 7.2921151467E-5;
		// Earth's universal gravitational parameter
		final double GM = 3.986004418E14;

		File orekitData = new File(
				"C:\\Users\\Naman\\Desktop\\rinex_parse_files\\orekit\\orekit-data-master\\orekit-data-master");
		DataProvidersManager manager = DataProvidersManager.getInstance();
		manager.addProvider(new DirectoryCrawler(orekitData));
		NormalizedSphericalHarmonicsProvider nhsp = GravityFieldFactory.getNormalizedProvider(50, 50);
		Frame frame = FramesFactory.getITRF(ITRFVersion.ITRF_2014, IERSConventions.IERS_2010, true);

		// ReferenceEllipsoid refElp = new ReferenceEllipsoid(ae, f, frame, GM, spin);
		Geoid geoid = new Geoid(nhsp, ReferenceEllipsoid.getWgs84(frame));
		return geoid;

	}

}
