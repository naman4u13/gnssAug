package com.gnssAug.IGS;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Optional;
import java.util.TimeZone;
import java.util.TreeMap;
import java.util.stream.IntStream;

import org.ejml.simple.SimpleMatrix;
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

import com.gnssAug.Android.constants.Measurement;
import com.gnssAug.Android.constants.State;
import com.gnssAug.Rinex.estimation.LinearLeastSquare;
import com.gnssAug.Rinex.fileParser.Antenna;
import com.gnssAug.Rinex.fileParser.Bias;
import com.gnssAug.Rinex.fileParser.Clock;
import com.gnssAug.Rinex.fileParser.IONEX;
import com.gnssAug.Rinex.fileParser.NavigationRNX;
import com.gnssAug.Rinex.fileParser.ObservationRNX;
import com.gnssAug.Rinex.fileParser.Orbit;
import com.gnssAug.Rinex.models.IonoCoeff;
import com.gnssAug.Rinex.models.NavigationMsg;
import com.gnssAug.Rinex.models.ObservationMsg;
import com.gnssAug.Rinex.models.SatResidual;
import com.gnssAug.Rinex.models.Satellite;
import com.gnssAug.Rinex.models.TimeCorrection;
import com.gnssAug.helper.ComputeEleAzm;
import com.gnssAug.helper.ComputeIonoCorr;
import com.gnssAug.helper.ComputeTropoCorr;
import com.gnssAug.utility.Analyzer;
import com.gnssAug.utility.GraphPlotter;
import com.gnssAug.utility.LatLonUtil;
import com.gnssAug.utility.MathUtil;
import com.gnssAug.utility.Time;

public class IGS {

	private static Geoid geoid = null;

	public static void posEstimate(String bias_path, String clock_path, String orbit_path, String ionex_path,
			String sinex_path, boolean useBias, boolean useGIM, boolean useIGS, boolean useSNX, String[] obsvCodeList,
			int minSat, double cutOffAng, double snrMask, boolean corrIono, boolean corrTropo, int estimatorType,
			boolean doAnalyze, boolean doTest, boolean outlierAnalyze) {
		try {
			TimeZone.setDefault(TimeZone.getTimeZone("UTC"));

			HashMap<String, ArrayList<double[]>> estPosMap = new HashMap<String, ArrayList<double[]>>();
			HashMap<String, ArrayList<double[]>> estVelMap = new HashMap<String, ArrayList<double[]>>();
			TreeMap<Long, ArrayList<Satellite>> satMap = new TreeMap<Long, ArrayList<Satellite>>();
			HashMap<String, HashMap<String, ArrayList<Double>>> satMeasNoiseMap = new HashMap<String, HashMap<String, ArrayList<Double>>>();

			HashMap<Measurement, HashMap<String, HashMap<String, ArrayList<SatResidual>>>> satResMap = new HashMap<Measurement, HashMap<String, HashMap<String, ArrayList<SatResidual>>>>();

			HashMap<Measurement, HashMap<String, ArrayList<Double>>> postVarOfUnitWeightMap = new HashMap<Measurement, HashMap<String, ArrayList<Double>>>();
			HashMap<State, HashMap<String, ArrayList<SimpleMatrix>>> Cxx_hat_map = new HashMap<State, HashMap<String, ArrayList<SimpleMatrix>>>();
			HashMap<String, ArrayList<double[]>> dopMap = new HashMap<String, ArrayList<double[]>>();
			HashMap<Measurement, TreeMap<String, ArrayList<Long>>> satCountMap = new HashMap<Measurement, TreeMap<String, ArrayList<Long>>>();
			ArrayList<Long> timeList = new ArrayList<Long>();
			Bias bias = null;
			Orbit orbit = null;
			Clock clock = null;
			Antenna antenna = null;
			IONEX ionex = null;
			
			String base_path = "C:\\Users\\naman.agarwal\\OneDrive - University of Calgary\\input_files";// "C:\\Users\\naman.agarwal\\Downloads\\GNSS\\input_files";

			String nav_path = base_path + "\\BRDC00IGS_R_20201000000_01D_MN.rnx\\BRDC00IGS_R_20201000000_01D_MN.rnx";

			String obs_path = base_path
					+ "\\LROC00FRA_R_20211190000_01D_30S_MO.crx\\LROC00FRA_R_20211190000_01D_30S_MO.rnx";

			String antenna_path = base_path + "\\complementary\\igs14.atx\\igs14.atx";

			String antenna_csv_path = base_path + "\\complementary\\antenna.csv";
			String path = "C:\\\\Users\\\\naman.agarwal\\\\Documents\\\\GNSS\\\\gnss_output\\\\IGS_rinex_output\\test2";
			// String path = "C:\\Users\\naman.agarwal\\Documents\\gnss_output\\test";
			File output = new File(path + ".txt");
			PrintStream stream;

			try {
				stream = new PrintStream(output);
				System.setOut(stream);
			} catch (FileNotFoundException e) { // TODO Auto-generated catch block
				e.printStackTrace();
			}

			geoid = buildGeoid();
			Map<String, Object> NavMsgComp = NavigationRNX.rinex_nav_process(nav_path, useIGS);
			@SuppressWarnings("unchecked")
			HashMap<Integer, ArrayList<NavigationMsg>> NavMsgs = (HashMap<Integer, ArrayList<NavigationMsg>>) NavMsgComp
					.getOrDefault("NavMsgs", null);
			IonoCoeff ionoCoeff = (IonoCoeff) NavMsgComp.get("ionoCoeff");
			TimeCorrection timeCorr = (TimeCorrection) NavMsgComp.getOrDefault("timeCorr", null);
			HashMap<String, Object> ObsvMsgComp = ObservationRNX.rinex_obsv_process(obs_path, useSNX, sinex_path,
					obsvCodeList, false, false);
			@SuppressWarnings("unchecked")
			ArrayList<ObservationMsg> ObsvMsgs = (ArrayList<ObservationMsg>) ObsvMsgComp.get("ObsvMsgs");
			// Note PVT algos will compute for Antenna Reference Point as it is independent
			// of frequency
			double[] rxARP = (double[]) ObsvMsgComp.get("ARP");
			HashMap<String, double[]> rxPCO = (HashMap<String, double[]>) ObsvMsgComp.get("PCO");

			int interval = (int) ObsvMsgComp.get("interval");

			if (useBias) {
				bias = new Bias(bias_path);

			}
			if (useIGS) {

				orbit = new Orbit(orbit_path);

				clock = new Clock(clock_path, bias);
				antenna = new Antenna(antenna_csv_path);

			}
			if (useGIM) {
				ionex = new IONEX(ionex_path);

			}
			double tRx0 = ObsvMsgs.get(0).getTRX();
			for (ObservationMsg obsvMsg : ObsvMsgs) {
				double tRx = obsvMsg.getTRX();
				double dayTime = tRx % 86400;
				long weekNo = obsvMsg.getWeekNo();
				Calendar time = Time.getDate(tRx, weekNo, 0);
				if (Time.getGPSTime(time)[0] != tRx) {
					System.err.println("FATAL ERROR TIME calendar");
					return;
				}
				ArrayList<Satellite> satList = null;
				satList = SingleFreq.process(obsvMsg, NavMsgs, obsvCodeList, useIGS, useBias, ionoCoeff, bias, orbit,
						clock, antenna, tRx, weekNo, time);
				if (satList.size() < minSat) {
					System.err.println("Less than " + minSat + " satellites");
					continue;
				}
				double[] refEcef = LinearLeastSquare.getEstPos(satList, rxPCO, false);
				satList.stream().forEach(i -> i.setElevAzm(ComputeEleAzm.computeEleAzm(rxARP, i.getSatEci())));
				filterSat(satList, rxARP, cutOffAng, snrMask, corrIono, corrTropo, ionex, ionoCoeff, time);
				if (satList.size() < minSat) {
					System.err.println("Less than " + minSat + " satellites");
					continue;
				}
				long tRxMilli = (long) (tRx * 1000);

				if (estimatorType == 1 || estimatorType == 2 || estimatorType == 3) {
					int[] arr = new int[] { estimatorType };
					if (estimatorType == 3) {
						arr = new int[] { 1, 2 };
					}
					for (int i : arr) {
						boolean isWLS = false;
						String estType = "LS";
						if (i == 2) {
							isWLS = true;
							estType = "WLS";
						}
						// Implement WLS method
						double[] estEcefClk = null;
						double[] estVel = null;
						for (Measurement type : new Measurement[] { Measurement.Pseudorange, Measurement.Doppler }) {
							if (type == Measurement.Pseudorange) {
								estEcefClk = LinearLeastSquare.getEstPos(satList, rxPCO, isWLS, doAnalyze, doTest,
										outlierAnalyze);
								estPosMap.computeIfAbsent(estType, k -> new ArrayList<double[]>()).add(estEcefClk);
							} else {
								estVel = LinearLeastSquare.getEstVel(satList, rxPCO, isWLS, doAnalyze, doTest,
										outlierAnalyze, estEcefClk);
								estVelMap.computeIfAbsent(estType, k -> new ArrayList<double[]>()).add(estVel);
							}
							if (doAnalyze) {

								double[] residual = LinearLeastSquare.getResidual(type);
								satResMap
										.computeIfAbsent(type,
												k -> new HashMap<String, HashMap<String, ArrayList<SatResidual>>>())
										.computeIfAbsent(estType, k -> new HashMap<String, ArrayList<SatResidual>>());
								ArrayList<Satellite> testedSatList = LinearLeastSquare.getTestedSatList(type);
								int n = testedSatList.size();
								for (int j = 0; j < n; j++) {
									Satellite sat = testedSatList.get(j);
									satResMap.get(type).get(estType)
											.computeIfAbsent(sat.getObsvCode().charAt(0) + "" + sat.getSVID(),
													k -> new ArrayList<SatResidual>())
											.add(new SatResidual(tRx - tRx0, sat.getElevAzm()[0], residual[j],
													sat.isOutlier(),sat.getCNo()));

								}
								if (doTest) {
									n = satList.size() - n;
								}
								satCountMap.computeIfAbsent(type, k -> new TreeMap<String, ArrayList<Long>>())
										.computeIfAbsent(estType, k -> new ArrayList<Long>()).add((long) n);
								postVarOfUnitWeightMap
										.computeIfAbsent(type, k -> new HashMap<String, ArrayList<Double>>())
										.computeIfAbsent(estType, k -> new ArrayList<Double>())
										.add(LinearLeastSquare.getPostVarOfUnitW(type));
								State state = (type != Measurement.Pseudorange) ? State.Velocity : State.Position;
								Cxx_hat_map.computeIfAbsent(state, k -> new HashMap<String, ArrayList<SimpleMatrix>>())
										.computeIfAbsent(estType, k -> new ArrayList<SimpleMatrix>())
										.add(LinearLeastSquare.getCxx_hat(type));
							}
							dopMap.computeIfAbsent(estType, k -> new ArrayList<double[]>())
									.add(LinearLeastSquare.getDop());
						}
					}
				}

				satMap.put(tRxMilli, satList);

				timeList.add(tRxMilli);
			}
			if (estimatorType == 3 || estimatorType == 4) {
				HashMap<Measurement, HashMap<String, HashMap<String, ArrayList<SatResidual>>>> satInnMap = new HashMap<Measurement, HashMap<String, HashMap<String, ArrayList<SatResidual>>>>();
				com.gnssAug.Rinex.estimation.EKF ekf = new com.gnssAug.Rinex.estimation.EKF();
				TreeMap<Long, double[]> estStateMap_pos = ekf.process(satMap, rxPCO, timeList, doAnalyze, doTest,
						outlierAnalyze, obsvCodeList);
				int n = timeList.size();
				if (doAnalyze) {
					satResMap.put(Measurement.Pseudorange,
							new HashMap<String, HashMap<String, ArrayList<SatResidual>>>());
					satResMap.get(Measurement.Pseudorange).put("EKF", new HashMap<String, ArrayList<SatResidual>>());
					satInnMap.put(Measurement.Pseudorange,
							new HashMap<String, HashMap<String, ArrayList<SatResidual>>>());
					satInnMap.get(Measurement.Pseudorange).put("EKF", new HashMap<String, ArrayList<SatResidual>>());
					satCountMap.put(Measurement.Pseudorange, new TreeMap<String, ArrayList<Long>>());
					Cxx_hat_map.put(State.Position, new HashMap<String, ArrayList<SimpleMatrix>>());
					postVarOfUnitWeightMap.put(Measurement.Pseudorange, new HashMap<String, ArrayList<Double>>());
					satMeasNoiseMap.put("EKF", new HashMap<String, ArrayList<Double>>());
					ArrayList<double[]> redundancyList = ekf.getRedundancyList();
					GraphPlotter.graphRedundancy(redundancyList);
				}
				for (int i = 1; i < n; i++) {
					long time = timeList.get(i);
					double[] estPos = estStateMap_pos.get(time);
					estPosMap.computeIfAbsent("EKF", k -> new ArrayList<double[]>()).add(estPos);
					if (doAnalyze) {
						ArrayList<Satellite> satList = ekf.getSatListMap().get(time);
						double[] residual = ekf.getResidualMap().get(time);
						int m = satList.size();
						long tRx = time / 1000;
						// double[] measNoise = ekf.getMeasNoiseMap().get(time);
						for (int j = 0; j < m; j++) {
							Satellite sat = satList.get(j);
							satResMap.get(Measurement.Pseudorange).get("EKF")
									.computeIfAbsent(sat.getSSI() + "" + sat.getSVID(),
											k -> new ArrayList<SatResidual>())
									.add(new SatResidual(tRx - tRx0, sat.getElevAzm()[0], residual[j],
											sat.isOutlier(),sat.getCNo()));

						}
						satCountMap.get(Measurement.Pseudorange).computeIfAbsent("EKF", k -> new ArrayList<Long>())
								.add(ekf.getSatCountMap().get(time));
						Cxx_hat_map.get(State.Position).computeIfAbsent("EKF", k -> new ArrayList<SimpleMatrix>())
								.add(ekf.getErrCovMap().get(time));
						postVarOfUnitWeightMap.get(Measurement.Pseudorange)
								.computeIfAbsent("EKF", k -> new ArrayList<Double>())
								.add(ekf.getPostVarOfUnitWMap().get(time));

						// For innovation vector
						satList = satMap.get(time);
						double[] innovation = ekf.getInnovationMap().get(time);
						m = satList.size();
						if (m != innovation.length) {
							throw new Exception("Fatal Error while mapping innovation sequence");
						}
						for (int j = 0; j < m; j++) {
							Satellite sat = satList.get(j);
							satInnMap.get(Measurement.Pseudorange).get("EKF")
									.computeIfAbsent(sat.getSSI() + "" + sat.getSVID(),
											k -> new ArrayList<SatResidual>())
									.add(new SatResidual(tRx - tRx0, sat.getElevAzm()[0], innovation[j],
											sat.isOutlier(),sat.getCNo()));

						}

					}
				}
				if (doAnalyze) {
					GraphPlotter.graphSatRes(satInnMap, outlierAnalyze, true);
				}
			}
			if (estimatorType == 1) {
				Analyzer.processIGS(satMap, rxARP, rxPCO, estPosMap, estVelMap, satResMap, outlierAnalyze);
			}
			// Calculate Accuracy Metrics
			HashMap<String, ArrayList<double[]>> GraphPosMap = new HashMap<String, ArrayList<double[]>>();
			HashMap<String, ArrayList<double[]>> GraphVelMap = new HashMap<String, ArrayList<double[]>>();

			for (String key : estPosMap.keySet()) {
				ArrayList<Double>[] posErrList = new ArrayList[5];

				IntStream.range(0, 5).forEach(i -> posErrList[i] = new ArrayList<Double>());

				ArrayList<double[]> estPosList = estPosMap.get(key);

				int n = estPosList.size();

				ArrayList<double[]> enuPosList = new ArrayList<double[]>();

				for (int i = 0; i < n; i++) {
					double[] estEcef = estPosList.get(i);
					if (estEcef == null) {
						continue;
					}
					double[] enu = LatLonUtil.ecef2enu(estEcef, rxARP, true);

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

				}

				GraphPosMap.put(key, enuPosList);

				// RMSE
				System.out.println("\n" + key);
				System.out.println("Position RMS - ");
				System.out.println(" E - " + MathUtil.RMS(posErrList[0]));
				System.out.println(" N - " + MathUtil.RMS(posErrList[1]));
				System.out.println(" U - " + MathUtil.RMS(posErrList[2]));
				System.out.println(" 3d Error - " + MathUtil.RMS(posErrList[3]));
				System.out.println(" 2d Error - " + MathUtil.RMS(posErrList[4]));

				// 95th Percentile

				IntStream.range(0, 5).forEach(i -> Collections.sort(posErrList[i]));
				int q95 = (int) (n * 0.95);

				System.out.println("\n" + key + " 95%");
				System.out.println("RMS - ");
				System.out.println(" E - " + posErrList[0].get(q95));
				System.out.println(" N - " + posErrList[1].get(q95));
				System.out.println(" U - " + posErrList[2].get(q95));
				System.out.println(" 3d Error - " + posErrList[3].get(q95));
				System.out.println(" 2d Error - " + posErrList[4].get(q95));

			}

			for (String key : estVelMap.keySet()) {
				ArrayList<Double>[] velErrList = new ArrayList[5];
				ArrayList<double[]> estVelList = null;
				ArrayList<double[]> enuVelList = null;
				IntStream.range(0, 5).forEach(i -> velErrList[i] = new ArrayList<Double>());
				enuVelList = new ArrayList<double[]>();
				estVelList = estVelMap.get(key);
				int n = estVelList.size();
				final double[] trueVel = new double[3];
				for (int i = 0; i < n; i++) {
					double[] estVel = estVelList.get(i);
					long time = timeList.get(i);
					if (estVel == null) {
						continue;
					}

					double[] velErr = IntStream.range(0, 3).mapToDouble(j -> estVel[j] - trueVel[j]).toArray();
					double[] enu = LatLonUtil.ecef2enu(velErr, rxARP, false);

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
//				for (int i = removalList.size() - 1; i >= 0; i--) {
//					Cxx_hat_map.get(State.Velocity).get(key).remove((int) removalList.get(i));
//				}

				GraphVelMap.put(key, enuVelList);

				// RMSE
				System.out.println("\n" + key);
				System.out.println("Velocity RMS - ");
				System.out.println(" E - " + MathUtil.RMS(velErrList[0]));
				System.out.println(" N - " + MathUtil.RMS(velErrList[1]));
				System.out.println(" U - " + MathUtil.RMS(velErrList[2]));
				System.out.println(" 3d Error - " + MathUtil.RMS(velErrList[3]));
				System.out.println(" 2d Error - " + MathUtil.RMS(velErrList[4]));

				// 95th Percentile

				IntStream.range(0, 5).forEach(i -> Collections.sort(velErrList[i]));
				int q95 = (int) (n * 0.95);

				System.out.println("\n" + key + " 95%");
				System.out.println("RMS - ");
				System.out.println(" E - " + velErrList[0].get(q95));
				System.out.println(" N - " + velErrList[1].get(q95));
				System.out.println(" U - " + velErrList[2].get(q95));
				System.out.println(" 3d Error - " + velErrList[3].get(q95));
				System.out.println(" 2d Error - " + velErrList[4].get(q95));

			}
			long t0 = timeList.get(0);
			for (int i = 0; i < timeList.size(); i++) {
				timeList.set(i, (long) ((timeList.get(i) - t0) * 1e-3));
			}

			// Plot Error Graphs
			if (Cxx_hat_map.isEmpty()) {
				GraphPlotter.graphENU(GraphPosMap, timeList, true);
				GraphPlotter.graphENU(GraphVelMap, timeList, false);
			} else {
				GraphPlotter.graphENU(GraphPosMap, timeList, true, Cxx_hat_map.get(State.Position));
				GraphPlotter.graphENU(GraphVelMap, timeList, false, Cxx_hat_map.get(State.Velocity));
			}

			if (doAnalyze) {

				GraphPlotter.graphSatRes(satResMap, outlierAnalyze);
				GraphPlotter.graphPostUnitW(postVarOfUnitWeightMap, timeList);
				// GraphPlotter.graphDOP(dopMap, satCountMap.get("EKF"), timeList);
				GraphPlotter.graphSatCount(satCountMap, timeList,30);
			}
		} catch (

		Exception e) {
			// TODO: handle exception
			System.out.println(e.getMessage());
			e.printStackTrace();
		}

	}

	public static void filterSat(ArrayList<Satellite> satList, double[] refEcef, double cutOffAng, double snrMask,
			boolean corrIono, boolean corrTropo, IONEX ionex, IonoCoeff ionoCoeff, Calendar time) {
		if (cutOffAng >= 0) {
			satList.removeIf(i -> i.getElevAzm()[0] < Math.toRadians(cutOffAng));
		}
		if (snrMask >= 0) {
			satList.removeIf(i -> i.getCNo() < snrMask);
		}
		if (corrIono || corrTropo) {
			double[] refLatLon = LatLonUtil.ecef2lla(refEcef);
			// Geocentric Latitude
			double gcLat = LatLonUtil.gd2gc(refLatLon[0], refLatLon[2]);
			int n = satList.size();
			for (int i = 0; i < n; i++) {
				Satellite sat = satList.get(i);
				double[] eleAzm = sat.getElevAzm();
				double ionoErr = 0;
				double tropoErr = 0;
				if (corrIono) {
					if (ionex != null) {
						ionoErr = ionex.computeIonoCorr(eleAzm[0], eleAzm[1], gcLat, refLatLon[1], sat.gettRX(),
								sat.getCarrier_frequency(), time);
					} else {
						if (Optional.ofNullable(ionoCoeff).isEmpty()) {
							System.out.println("You have not provided IonoCoeff");
							return;
						}
						ionoErr = ComputeIonoCorr.computeIonoCorr(eleAzm[0], eleAzm[1], refLatLon[0], refLatLon[1],
								sat.gettRX(), ionoCoeff, sat.getCarrier_frequency(), time);
					}
				}
				if (corrTropo) {
					ComputeTropoCorr tropo = new ComputeTropoCorr(refLatLon, time, geoid);
					tropoErr = tropo.getSlantDelay(eleAzm[0]);
				}
				sat.setPseudorange(sat.getPseudorange() - ionoErr - tropoErr);
			}
		}
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
				"C:\\Users\\naman.agarwal\\Documents\\orekit\\");
		/*
		 * File orekitData = new File(
		 * "C:\\Users\\naman.agarwal\\Downloads\\orekit\\orekit\\orekit-data-master\\orekit-data-master"
		 * );
		 */

		DataProvidersManager manager = DataProvidersManager.getInstance();
		manager.addProvider(new DirectoryCrawler(orekitData));
		NormalizedSphericalHarmonicsProvider nhsp = GravityFieldFactory.getNormalizedProvider(50, 50);
		Frame frame = FramesFactory.getITRF(ITRFVersion.ITRF_2014, IERSConventions.IERS_2010, true);

		// ReferenceEllipsoid refElp = new ReferenceEllipsoid(ae, f, frame, GM, spin);
		Geoid geoid = new Geoid(nhsp, ReferenceEllipsoid.getWgs84(frame));
		return geoid;

	}

}
