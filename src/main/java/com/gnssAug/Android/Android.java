package com.gnssAug.Android;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TimeZone;
import java.util.TreeMap;
import java.util.stream.IntStream;

import org.apache.commons.collections.set.ListOrderedSet;
import org.ejml.simple.SimpleMatrix;
import org.orekit.data.DataContext;
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
import com.gnssAug.Android.constants.GnssDataConfig;
import com.gnssAug.Android.constants.Measurement;
import com.gnssAug.Android.constants.State;
import com.gnssAug.Android.estimation.LLS_TDCP;
import com.gnssAug.Android.estimation.LLS_TDCP_ambFix;
import com.gnssAug.Android.estimation.LinearLeastSquare;
import com.gnssAug.Android.estimation.KalmanFilter.AKFDoppler;
import com.gnssAug.Android.estimation.KalmanFilter.EKF;
import com.gnssAug.Android.estimation.KalmanFilter.EKFDoppler;

import com.gnssAug.Android.estimation.KalmanFilter.EKF_TDCP_ambFix;
import com.gnssAug.Android.estimation.KalmanFilter.EKF_TDCP_ambFix2;
import com.gnssAug.Android.estimation.KalmanFilter.INSfusion;
import com.gnssAug.Android.estimation.KalmanFilter.EKFParent;
import com.gnssAug.Android.estimation.KalmanFilter.EKF_PPP;
import com.gnssAug.Android.estimation.KalmanFilter.Models.Flag;
import com.gnssAug.Android.fileParser.DerivedCSV;
import com.gnssAug.Android.fileParser.GNSS_Log;
import com.gnssAug.Android.fileParser.GroundTruth;
import com.gnssAug.Android.fileParser.GroundTruth_GSA;
import com.gnssAug.Android.models.CycleSlipDetect;
import com.gnssAug.Android.models.Derived;
import com.gnssAug.Android.models.GNSSLog;
import com.gnssAug.Android.models.IMUsensor;
import com.gnssAug.Android.models.Satellite;
import com.gnssAug.Android.models.TDCP;
import com.gnssAug.Rinex.fileParser.DCB_Bias;
import com.gnssAug.Rinex.fileParser.Antenna;
import com.gnssAug.Rinex.fileParser.Clock;
import com.gnssAug.Rinex.fileParser.IONEX;
import com.gnssAug.Rinex.fileParser.OSB_Bias;
import com.gnssAug.Rinex.fileParser.Orbit;
import com.gnssAug.Rinex.models.SatResidual;
import com.gnssAug.helper.ComputeEleAzm;
import com.gnssAug.helper.ComputeTropoCorr;
import com.gnssAug.helper.INS.IMUconfigure;
import com.gnssAug.helper.INS.StateInitialization;
import com.gnssAug.utility.Analyzer;
import com.gnssAug.utility.GraphPlotter;
import com.gnssAug.utility.LatLonUtil;
import com.gnssAug.utility.MathUtil;
import com.gnssAug.utility.SatUtil;
import com.gnssAug.utility.Time;

public class Android {
	private static Geoid geoid = null;

	public static void posEstimate(boolean doPosErrPlot, double cutOffAng, double snrMask, int estimatorType,
			String[] obsvCodeList, String derived_csv_path, String gnss_log_path, String GTcsv, String bias_path,
			String clock_path, String orbit_path, String ionex_path, boolean useIGS, boolean doAnalyze, boolean doTest,
			boolean outlierAnalyze, boolean mapDeltaRanges, Set<String> discardSet, boolean isGSA,boolean repairCS) {
		try {
			TimeZone.setDefault(TimeZone.getTimeZone("UTC"));
			HashMap<String, ArrayList<HashMap<String, Double>>> ErrMap = new HashMap<String, ArrayList<HashMap<String, Double>>>();

			ArrayList<Long> timeList = new ArrayList<Long>();
			ArrayList<double[]> trueLLHlist = new ArrayList<double[]>();
			ArrayList<double[]> trueEcefList = new ArrayList<double[]>();

			TreeMap<Long, ArrayList<Satellite>> satMap = new TreeMap<Long, ArrayList<Satellite>>();
			TreeMap<String, ArrayList<double[]>> estPosMap = new TreeMap<String, ArrayList<double[]>>();
			TreeMap<String, ArrayList<double[]>> estVelMap = new TreeMap<String, ArrayList<double[]>>();
			HashMap<Measurement, HashMap<String, HashMap<String, ArrayList<SatResidual>>>> satResMap = new HashMap<Measurement, HashMap<String, HashMap<String, ArrayList<SatResidual>>>>();
			HashMap<Measurement, HashMap<String, ArrayList<Double>>> postVarOfUnitWeightMap = new HashMap<Measurement, HashMap<String, ArrayList<Double>>>();
			HashMap<State, HashMap<String, ArrayList<SimpleMatrix>>> Cxx_hat_map = new HashMap<State, HashMap<String, ArrayList<SimpleMatrix>>>();
			HashMap<String, ArrayList<double[]>> dopMap = new HashMap<String, ArrayList<double[]>>();
			HashMap<Measurement, TreeMap<String, ArrayList<Long>>> satCountMap = new HashMap<Measurement, TreeMap<String, ArrayList<Long>>>();

			DCB_Bias bias = null;
			Orbit orbit = null;
			Clock clock = null;
			IONEX ionex = null;
			Antenna antenna = null;
			OSB_Bias osb_bias = null;
			String base_path = "/Users/naman.agarwal/Library/CloudStorage/OneDrive-UniversityofCalgary/input_files";
			String antenna_csv_path = base_path + "/complementary/antenna.csv";
			String path = "/Users/naman.agarwal/Library/CloudStorage/OneDrive-UniversityofCalgary/gnss_output/T-A-SIS-02_open_sky_pedestrian/ION_GNSS_2025/Samsung_Galaxy_S20+_5G_test";
			// "C:\\Users\\Naman\\Desktop\\rinex_parse_files\\google2\\2021-04-28-US-MTV-1\\test2";
			File output = new File(path + ".txt");
			PrintStream stream;
			stream = new PrintStream(output);
			System.setOut(stream);
			System.out.println("Discarded Satellites: " + discardSet.toString());
			System.out.println("Elevation Mask in degrees = " + cutOffAng);
			System.out
					.println("pseudorange_priorStdOfUnitW = " + Math.sqrt(GnssDataConfig.pseudorange_priorVarOfUnitW));
			System.out.println("doppler_priorStdOfUnitW = " + Math.sqrt(GnssDataConfig.doppler_priorVarOfUnitW));
			System.out.println("TDCP_priorStdOfUnitW = " + Math.sqrt(GnssDataConfig.tdcp_priorVarOfUnitW));
			System.out.println("Phase_priorStdOfUnitW = " + Math.sqrt(GnssDataConfig.phase_priorVarOfUnitW));
			System.out.println("GIM_TECU_priorStdOfUnitW = " + Math.sqrt(GnssDataConfig.GIM_TECU_variance));
			System.out.println("Q matrix for pos_rand_walk = " + Arrays.toString(GnssDataConfig.qENU_posRandWalk));
			System.out.println("Q matrix for vel_rand_walk = " + Arrays.toString(GnssDataConfig.qENU_velRandWalk));
			ArrayList<double[]> rxGT = new ArrayList<double[]>();
			if (isGSA) {
				// Its ECEF not LLH
				rxGT = GroundTruth_GSA.processCSV(GTcsv);
			} else {
				rxGT = GroundTruth.processCSV(GTcsv);
			}

			HashMap<Long, HashMap<String, HashMap<Integer, Derived>>> derivedMap = null;

			if (derived_csv_path != null) {
				derivedMap = DerivedCSV.processCSV(derived_csv_path);
			}
			GNSS_Log.process(gnss_log_path);
			TreeMap<Long, HashMap<String, ArrayList<GNSSLog>>> gnssLogMaps = GNSS_Log.getGnssLogMaps();
			ArrayList<IMUsensor> imuList = GNSS_Log.getImuList();

			if (useIGS) {

				orbit = new Orbit(orbit_path);
				bias = new DCB_Bias(bias_path);
				clock = new Clock(clock_path, bias);
				ionex = new IONEX(ionex_path);
				geoid = buildGeoid();
				antenna = new Antenna(antenna_csv_path);

			}
			ListOrderedSet ssiSet = new ListOrderedSet();
			for (int i = 0; i < obsvCodeList.length; i++) {
				ssiSet.add(obsvCodeList[i].charAt(0));
			}
			double[] refUserEcef = null;
			int gtIndex = 0;
			double tRx0 = ((ArrayList<GNSSLog>) gnssLogMaps.firstEntry().getValue().values().toArray()[0]).get(0)
					.gettRx();
			for (long tRxMilli : gnssLogMaps.keySet()) {

				if (gtIndex >= rxGT.size()) {
					break;
				}
				HashMap<String, ArrayList<GNSSLog>> gnssLogMap = gnssLogMaps.get(tRxMilli);
				GNSSLog entry = ((ArrayList<GNSSLog>) gnssLogMap.values().toArray()[0]).get(0);
				double tRx = entry.gettRx();
				int weekNo = entry.getWeekNo();

				if ((Math.abs(tRxMilli - (rxGT.get(gtIndex)[0] * 1000))) > 1 || weekNo != rxGT.get(gtIndex)[1]) {
					System.err.println("FATAL ERROR - GT timestamp does not match");
					if(isGSA)
					{
						while(((tRxMilli - (rxGT.get(gtIndex)[0] * 1000)))>0)
						{
							gtIndex++;
						}
						if(((tRxMilli - (rxGT.get(gtIndex)[0] * 1000)))<0)
						{
							continue;
						}
					}
					else
					{
						continue;
					}
				}

				double[] trueUserLLH = null;

				double[] truePosEcef = null;
				if (isGSA) {
					truePosEcef = new double[] { rxGT.get(gtIndex)[2], rxGT.get(gtIndex)[3], rxGT.get(gtIndex)[4] };
					trueUserLLH = LatLonUtil.ecef2lla(truePosEcef);

				} else {
					trueUserLLH = new double[] { rxGT.get(gtIndex)[2], rxGT.get(gtIndex)[3], rxGT.get(gtIndex)[4] };
					truePosEcef = LatLonUtil
							.lla2ecef(new double[] { trueUserLLH[0], trueUserLLH[1], trueUserLLH[2] - 61 }, true);
				}
				gtIndex++;
				Calendar time = Time.getDate(tRx, weekNo, 0);
				
				ArrayList<Satellite> satList = SingleFreq.process(tRx, derivedMap,osb_bias, antenna,  gnssLogMap, time, obsvCodeList,weekNo, clock, orbit, useIGS, discardSet,refUserEcef);
				int m = obsvCodeList.length;
				if (satList.size() < 3 + m) {
					System.err.println("Less than " + (3 + m) + " satellites");

					continue;

				}
				refUserEcef = new double[3];
				try {
					refUserEcef = LinearLeastSquare.getEstPos(satList, false, useIGS);
				} catch (org.ejml.data.SingularMatrixException e) {
					// TODO: handle exception
					e.printStackTrace();
					continue;
				}

				for (Satellite sat : satList) {
					sat.setElevAzm(ComputeEleAzm.computeEleAzm(truePosEcef, sat.getSatEci()));

				}
				filterSat(satList, cutOffAng, snrMask, truePosEcef, useIGS, ionex, time, estimatorType);
				if (satList.size() < 3 + m) {
					System.err.println("Less than " + (3 + m) + " satellites");

					continue;

				}
				double[] estEcefClk = null;

				if (estimatorType == 1 || estimatorType == 2 || estimatorType == 3 || estimatorType == 11) {
					int[] arr = new int[] { estimatorType };
					if (estimatorType == 3 || estimatorType == 11) {
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
						estEcefClk = null;
						double[] estVel = null;
						for (Measurement type : new Measurement[] { Measurement.Pseudorange, Measurement.Doppler }) {
							if (type == Measurement.Pseudorange) {
								estEcefClk = LinearLeastSquare.getEstPos(satList, isWLS, doAnalyze, doTest,
										outlierAnalyze, useIGS);
								estPosMap.computeIfAbsent(estType, k -> new ArrayList<double[]>()).add(estEcefClk);
							} else {
								estVel = LinearLeastSquare.getEstVel(satList, isWLS, doAnalyze, doTest, outlierAnalyze,
										estEcefClk, useIGS);
//								final double data = estVel[4];
//								satList.forEach(j->j.setClkRate(data));
								estVelMap.computeIfAbsent(estType, k -> new ArrayList<double[]>()).add(estVel);
							}
							if (doAnalyze && estimatorType != 11) {

								double[] residual = LinearLeastSquare.getResidual(type);
								SimpleMatrix Cyy = LinearLeastSquare.getCyy(type);
								satResMap
										.computeIfAbsent(type,
												k -> new HashMap<String, HashMap<String, ArrayList<SatResidual>>>())
										.computeIfAbsent(estType, k -> new HashMap<String, ArrayList<SatResidual>>());
								ArrayList<Satellite> testedSatList = LinearLeastSquare.getTestedSatList(type);
								int n = testedSatList.size();
								for (int j = 0; j < n; j++) {
									Satellite sat = testedSatList.get(j);
									satResMap.get(type).get(estType)
											.computeIfAbsent(sat.getObsvCode().charAt(0) + "" + sat.getSvid(),
													k -> new ArrayList<SatResidual>())
											.add(new SatResidual(tRx - tRx0, sat.getElevAzm()[0], residual[j],
													sat.isOutlier(), Math.sqrt(Cyy.get(j, j)), sat.getCn0DbHz()));

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
										.add(LinearLeastSquare.getCxx_hat(type, "ENU"));
							}
							if (type == Measurement.Pseudorange) {
								dopMap.computeIfAbsent(estType, k -> new ArrayList<double[]>())
										.add(LinearLeastSquare.getDop());
							}
						}
					}
				}

				satMap.put(tRxMilli, satList);

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
				double[][] dcm = StateInitialization.initialize(imuMap, satMap, useIGS);
				TreeMap<Long, double[]> ecefMap = INSfusion.process(imuMap, satMap, timeList, dcm, useIGS);
				// GraphPlotter.graphGnssIns(ecefMap, trueECEFlist, timeList);
				int n = timeList.size();
				for (int i = 0; i < n; i++) {
					long time = timeList.get(i);
					double[] estEcef = ecefMap.get(time);
					estPosMap.computeIfAbsent("GNSS/INS fusion", k -> new ArrayList<double[]>()).add(estEcef);

				}

			}
			boolean useDoppler = true;
			if (estimatorType == 5 || estimatorType == 6 || estimatorType == 7 || estimatorType == 8
					|| estimatorType == 9 || estimatorType == 11 || estimatorType == 12 || estimatorType == 13) {
				int m = obsvCodeList.length;
				EKF ekf = new EKF();

				TreeMap<Long, double[]> estStateMap_pos = null;
				TreeMap<Long, double[]> estStateMap_vel = null;
				TreeMap<Long, double[]> estStateMap_vel_doppler = null;
				TreeMap<Long, double[]> estStateMap_vel_doppler_complementary = null;
				TreeMap<Long, double[]> estStateMap_vel_estVel = null;
				TreeMap<Long, double[]> estStateMap_vel_estVel_complementary = null;
				int n = timeList.size();
				int[] estTypes = new int[] { estimatorType };
				String estName = "";
				if (((estimatorType == 9 && (!doAnalyze)) || (estimatorType == 11))) {
					estTypes = new int[] { 5, 6, 12 };
				}
				for (int type : estTypes) {
					switch (type) {
					case 5:
						// Implement EKF based on receiver’s position and clock offset errors as a
						// random walk process
						useDoppler = false;
						estStateMap_pos = ekf.process(satMap, timeList, Flag.POSITION, false, useIGS, obsvCodeList,
								doAnalyze, doTest, outlierAnalyze);
						estName = "PRW";
						for (int i = 0; i < n; i++) {
							long time = timeList.get(i);
							double[] estPos = estStateMap_pos.get(time);
							estPosMap.computeIfAbsent(estName, k -> new ArrayList<double[]>()).add(estPos);
						}

						break;
					case 6:
						// Implement EKF based on receiver’s velocity and clock drift errors as a random
						// walk process
						useDoppler = false;
						estStateMap_vel = ekf.process(satMap, timeList, Flag.VELOCITY, false, useIGS, obsvCodeList,
								doAnalyze, doTest, outlierAnalyze);
						estName = "VRW";
						for (int i = 0; i < n; i++) {
							long time = timeList.get(i);
							double[] estState = estStateMap_vel.get(time);
							double[] estPos = null;
							double[] estVel = null;
							if (estState != null) {
								estPos = new double[3 + m];
								estVel = new double[3 + m];
								for (int j = 0; j < 3 + m; j++) {
									estPos[j] = estState[j];
									estVel[j] = estState[j + 3 + m];
								}

							}
							estPosMap.computeIfAbsent(estName, k -> new ArrayList<double[]>()).add(estPos);
							estVelMap.computeIfAbsent(estName, k -> new ArrayList<double[]>()).add(estVel);

						}

						break;

					case 7:
						// Implement EKF based on receiver’s velocity and clock drift errors as a random
						// walk process along with doppler updates
						estStateMap_vel_doppler = ekf.process(satMap, timeList, Flag.VELOCITY, true, useIGS,
								obsvCodeList, doAnalyze, doTest, outlierAnalyze);

						estName = "EKF - vel. random walk + doppler";
						for (int i = 0; i < n; i++) {
							long time = timeList.get(i);
							double[] estState = estStateMap_vel_doppler.get(time);
							double[] estPos = null;
							double[] estVel = null;
							if (estState != null) {
								estPos = new double[3 + m];
								estVel = new double[3 + m];
								for (int j = 0; j < 3 + m; j++) {
									estPos[j] = estState[j];
									estVel[j] = estState[j + 3 + m];
								}
							}
							estPosMap.computeIfAbsent(estName, k -> new ArrayList<double[]>()).add(estPos);
							estVelMap.computeIfAbsent(estName, k -> new ArrayList<double[]>()).add(estVel);

						}

						break;

					// case 8:
					//
					// estStateMap_vel_doppler_complementary = ekf.process(satMap, timeList,
					// Flag.VELOCITY, true,
					// useIGS, obsvCodeList, doAnalyze, doTest, outlierAnalyze, true);
					//
					// estName = "EKF - vel. random walk + doppler + complementary equivalent";
					// for (int i = 0; i < n; i++) {
					// long time = timeList.get(i);
					// double[] estState = estStateMap_vel_doppler_complementary.get(time);
					// double[] estPos = null;
					// double[] estVel = null;
					// if (estState != null) {
					// estPos = new double[3+m];
					// estVel = new double[3+m];
					// for(int j=0;j<3+m;j++)
					// {
					// estPos[j] = estState[j];
					// estVel[j] = estState[j+3+m];
					// }
					// }
					// estPosMap.computeIfAbsent(estName, k -> new
					// ArrayList<double[]>()).add(estPos);
					// estVelMap.computeIfAbsent(estName, k -> new
					// ArrayList<double[]>()).add(estVel);
					//
					// }
					//
					// break;
					case 12:
						useDoppler = false;
						// Implement EKF based on receiver’s velocity and clock drift errors as a random
						// walk process along with estimated velocity updates
						estStateMap_vel_estVel = ekf.process(satMap, timeList, Flag.VELOCITY, false, useIGS,
								obsvCodeList, doAnalyze, doTest, outlierAnalyze, false, true);

						estName = "VRWD";
						for (int i = 0; i < n; i++) {
							long time = timeList.get(i);
							double[] estState = estStateMap_vel_estVel.get(time);
							double[] estPos = null;
							double[] estVel = null;
							if (estState != null) {
								estPos = new double[3 + m];
								estVel = new double[3 + m];
								for (int j = 0; j < 3 + m; j++) {
									estPos[j] = estState[j];
									estVel[j] = estState[j + 3 + m];
								}
							}
							estPosMap.computeIfAbsent(estName, k -> new ArrayList<double[]>()).add(estPos);
							estVelMap.computeIfAbsent(estName, k -> new ArrayList<double[]>()).add(estVel);

						}

						break;
					case 13:
						useDoppler = false;
						// Implement EKF based on receiver’s velocity and clock drift errors as a random
						// walk process along with estimated velocity updates
						estStateMap_vel_estVel_complementary = ekf.process(satMap, timeList, Flag.VELOCITY, false,
								useIGS, obsvCodeList, doAnalyze, doTest, outlierAnalyze, true, true);

						estName = "EKF - vel. random walk + estVel + complementary equivalent";
						for (int i = 0; i < n; i++) {
							long time = timeList.get(i);
							double[] estState = estStateMap_vel_estVel_complementary.get(time);
							double[] estPos = null;
							double[] estVel = null;
							if (estState != null) {
								estPos = new double[3 + m];
								estVel = new double[3 + m];
								for (int j = 0; j < 3 + m; j++) {
									estPos[j] = estState[j];
									estVel[j] = estState[j + 3 + m];
								}
							}
							estPosMap.computeIfAbsent(estName, k -> new ArrayList<double[]>()).add(estPos);
							estVelMap.computeIfAbsent(estName, k -> new ArrayList<double[]>()).add(estVel);

						}

					}
				}

				if (doAnalyze && estimatorType != 11) {
					Measurement[] measArr = useDoppler
							? new Measurement[] { Measurement.Pseudorange, Measurement.Doppler }
							: new Measurement[] { Measurement.Pseudorange };
					HashMap<Measurement, HashMap<String, HashMap<String, ArrayList<SatResidual>>>> satInnMap = new HashMap<Measurement, HashMap<String, HashMap<String, ArrayList<SatResidual>>>>();
					for (Measurement meas : measArr) {
						satResMap.put(meas, new HashMap<String, HashMap<String, ArrayList<SatResidual>>>());
						satInnMap.put(meas, new HashMap<String, HashMap<String, ArrayList<SatResidual>>>());
						postVarOfUnitWeightMap.put(meas, new HashMap<String, ArrayList<Double>>());
						satResMap.get(meas).put(estName, new HashMap<String, ArrayList<SatResidual>>());
						satInnMap.get(meas).put(estName, new HashMap<String, ArrayList<SatResidual>>());
						satCountMap.put(meas, new TreeMap<String, ArrayList<Long>>());
					}
					Cxx_hat_map.put(State.Position, new HashMap<String, ArrayList<SimpleMatrix>>());
					Cxx_hat_map.put(State.Velocity, new HashMap<String, ArrayList<SimpleMatrix>>());
					ArrayList<double[]> redundancyList = ekf.getRedundancyList(Measurement.Pseudorange);
					GraphPlotter.graphRedundancy(redundancyList, "Pseudorange ");
					if (useDoppler) {
						redundancyList = ekf.getRedundancyList(Measurement.Doppler);
						GraphPlotter.graphRedundancy(redundancyList, "Doppler ");
					}
					for (int i = 0; i < n; i++) {
						long time = timeList.get(i);
						HashMap<Measurement, ArrayList<Satellite>> satListMap = ekf.getSatListMap().get(time);

						if (satListMap == null) {
							continue;
						}
						HashMap<Measurement, double[]> residualMap = ekf.getResidualMap().get(time);
						HashMap<Measurement, Long> satCount = ekf.getSatCountMap().get(time);
						HashMap<Measurement, Double> postVarOfUnitW = ekf.getPostVarOfUnitWMap().get(time);
						HashMap<Measurement, double[]> innovation = ekf.getInnovationMap().get(time);
						ArrayList<Satellite> originalSatList = satMap.get(time);
						double tRx = time / 1000.0;
						for (Measurement meas : satListMap.keySet()) {
							ArrayList<Satellite> satList = satListMap.get(meas);
							int l = satList.size();
							for (int j = 0; j < l; j++) {
								Satellite sat = satList.get(j);
								satResMap.get(meas).get(estName)
										.computeIfAbsent(sat.getObsvCode().charAt(0) + "" + sat.getSvid() + "",
												k -> new ArrayList<SatResidual>())
										.add(new SatResidual(tRx - tRx0, sat.getElevAzm()[0], residualMap.get(meas)[j],
												sat.isOutlier(), sat.getCn0DbHz()));
								if (!outlierAnalyze) {
									sat = originalSatList.get(j);
								}
								satInnMap.get(meas).get(estName)
										.computeIfAbsent(sat.getObsvCode().charAt(0) + sat.getSvid() + "",
												k -> new ArrayList<SatResidual>())
										.add(new SatResidual(tRx - tRx0, sat.getElevAzm()[0], innovation.get(meas)[j],
												sat.isOutlier(), sat.getCn0DbHz()));

							}
							satCountMap.get(meas).computeIfAbsent(estName, k -> new ArrayList<Long>())
									.add(satCount.get(meas));
							postVarOfUnitWeightMap.get(meas).computeIfAbsent(estName, k -> new ArrayList<Double>())
									.add(postVarOfUnitW.get(meas));
						}
						for (State state : State.values()) {
							Cxx_hat_map.get(state).computeIfAbsent(estName, k -> new ArrayList<SimpleMatrix>())
									.add(ekf.getErrCovMap().get(time).get(state));

						}

					}

					GraphPlotter.graphSatRes(satInnMap, outlierAnalyze, true);

				}

			}

			if (estimatorType == 10 || estimatorType == 14 || estimatorType == 11) {
				int[] estArray = new int[] { estimatorType };
				if (estimatorType == 11) {
					estArray = new int[] { 10, 14 };
				}
				for (int estType : estArray) {
					String estName;
					TreeMap<Long, double[]> estStateMap;
					EKFParent ekf;
					if (estType == 10) {
						estName = "DBP Filter";
						ekf = new EKFDoppler();
						estStateMap = ((EKFDoppler) ekf).process(satMap, timeList, useIGS, obsvCodeList, doAnalyze,
								doTest, outlierAnalyze);
					} else {
						estName = "Proposed AKF";
						ekf = new AKFDoppler();
						estStateMap = ((AKFDoppler) ekf).process(satMap, timeList, useIGS, obsvCodeList, doAnalyze,
								doTest, true,outlierAnalyze);
					}

					// Doppler is used but no velocity is computed therefore useDoppler is set as
					// false
					useDoppler = false;

					int n = timeList.size();
					estPosMap.put(estName, new ArrayList<double[]>());
					HashMap<Measurement, HashMap<String, HashMap<String, ArrayList<SatResidual>>>> satInnMap = new HashMap<Measurement, HashMap<String, HashMap<String, ArrayList<SatResidual>>>>();
					if (doAnalyze && estimatorType != 11) {
						satResMap.put(Measurement.Pseudorange,
								new HashMap<String, HashMap<String, ArrayList<SatResidual>>>());
						satResMap.get(Measurement.Pseudorange).put(estName,
								new HashMap<String, ArrayList<SatResidual>>());
						satInnMap.put(Measurement.Pseudorange,
								new HashMap<String, HashMap<String, ArrayList<SatResidual>>>());
						satInnMap.get(Measurement.Pseudorange).put(estName,
								new HashMap<String, ArrayList<SatResidual>>());
						satCountMap.put(Measurement.Pseudorange, new TreeMap<String, ArrayList<Long>>());
						Cxx_hat_map.put(State.Position, new HashMap<String, ArrayList<SimpleMatrix>>());
						postVarOfUnitWeightMap.put(Measurement.Pseudorange, new HashMap<String, ArrayList<Double>>());
						ArrayList<double[]> redundancyList = ekf.getRedundancyList();
						GraphPlotter.graphRedundancy(redundancyList);
					}
					for (int i = 0; i < n; i++) {
						long time = timeList.get(i);
						double[] estPos = estStateMap.get(time);
						estPosMap.get(estName).add(estPos);
						if (estPos == null) {
							continue;
						}
						if (doAnalyze && estimatorType != 11) {
							ArrayList<Satellite> satList = ekf.getSatListMap().get(time);
							double[] residual = ekf.getResidualMap().get(time);
							int m = satList.size();
							double tRx = time / 1000.0;
							// double[] measNoise = ekf.getMeasNoiseMap().get(time);
							for (int j = 0; j < m; j++) {
								Satellite sat = satList.get(j);
								satResMap.get(Measurement.Pseudorange).get(estName)
										.computeIfAbsent(sat.getObsvCode().charAt(0) + "" + sat.getSvid(),
												k -> new ArrayList<SatResidual>())
										.add(new SatResidual(tRx - tRx0, sat.getElevAzm()[0], residual[j],
												sat.isOutlier(), sat.getCn0DbHz()));

							}
							satCountMap.get(Measurement.Pseudorange)
									.computeIfAbsent(estName, k -> new ArrayList<Long>())
									.add(ekf.getSatCountMap().get(time));
							Cxx_hat_map.get(State.Position).computeIfAbsent(estName, k -> new ArrayList<SimpleMatrix>())
									.add(ekf.getErrCovMap().get(time));
							postVarOfUnitWeightMap.get(Measurement.Pseudorange)
									.computeIfAbsent(estName, k -> new ArrayList<Double>())
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
								satInnMap.get(Measurement.Pseudorange).get(estName)
										.computeIfAbsent(sat.getObsvCode().charAt(0) + sat.getSvid() + "",
												k -> new ArrayList<SatResidual>())
										.add(new SatResidual(tRx - tRx0, sat.getElevAzm()[0], innovation[j],
												sat.isOutlier(), sat.getCn0DbHz()));

							}

						}
					}
					if (doAnalyze && estimatorType != 11) {
						GraphPlotter.graphSatRes(satInnMap, outlierAnalyze, true);
					}
				}

			}

//			if (estimatorType == 15) {
//				TreeMap<Long, double[]> estStateMap;
//				EKFParent ekf;
//				ekf = new EKF_PPP();
//				estStateMap = ((EKF_PPP) ekf).process(satMap, timeList, useIGS, obsvCodeList, doAnalyze, doTest,
//						outlierAnalyze);
//			}

			if (estimatorType == 16) {
				String estType = "LS TDCP";
				// TreeMap<String, ArrayList<double[]>> estTdcpVelMap = new TreeMap<String,
				// ArrayList<double[]>>();
				long prevTime = timeList.get(0);
				for (int i = 1; i < timeList.size(); i++) {
					long currentTime = timeList.get(i);
					ArrayList<Satellite> currentSatList = satMap.get(currentTime);
					ArrayList<Satellite> prevSatList = satMap.get(prevTime);
					double[] refPos = LinearLeastSquare.getEstPos(currentSatList, true, useIGS);
					double[] estVel = LLS_TDCP.getEstVel(currentSatList, prevSatList, false, doAnalyze, doTest,
							outlierAnalyze, refPos, useIGS, true);
					estVelMap.computeIfAbsent(estType, k -> new ArrayList<double[]>()).add(estVel);
					prevTime = currentTime;
					if (doAnalyze && estimatorType != 11) {
						double tRx = currentTime / 1000.0;
						double[] residual = LLS_TDCP.getResidual();
						SimpleMatrix Cyy = LLS_TDCP.getCyy();
						satResMap
								.computeIfAbsent(Measurement.TDCP,
										k -> new HashMap<String, HashMap<String, ArrayList<SatResidual>>>())
								.computeIfAbsent(estType, k -> new HashMap<String, ArrayList<SatResidual>>());
						ArrayList<TDCP> testedTdcpList = LLS_TDCP.getTestedTdcpList();
						int n = testedTdcpList.size();
						for (int j = 0; j < n; j++) {
							TDCP tdcp = testedTdcpList.get(j);
							Satellite sat = tdcp.getSat();
							satResMap.get(Measurement.TDCP).get(estType)
									.computeIfAbsent(sat.getObsvCode().charAt(0) + "" + sat.getSvid(),
											k -> new ArrayList<SatResidual>())
									.add(new SatResidual(tRx - tRx0, sat.getElevAzm()[0], residual[j], sat.isOutlier(),
											Math.sqrt(Cyy.get(j, j)), sat.getCn0DbHz()));

						}
						if (doTest) {
							n = LLS_TDCP.getCommSatCount() - n;
						}
						satCountMap.computeIfAbsent(Measurement.TDCP, k -> new TreeMap<String, ArrayList<Long>>())
								.computeIfAbsent(estType, k -> new ArrayList<Long>()).add((long) n);
						postVarOfUnitWeightMap
								.computeIfAbsent(Measurement.TDCP, k -> new HashMap<String, ArrayList<Double>>())
								.computeIfAbsent(estType, k -> new ArrayList<Double>())
								.add(LLS_TDCP.getPostVarOfUnitW());
						State state = State.Velocity;
						Cxx_hat_map.computeIfAbsent(state, k -> new HashMap<String, ArrayList<SimpleMatrix>>())
								.computeIfAbsent(estType, k -> new ArrayList<SimpleMatrix>())
								.add(LLS_TDCP.getCxx_hat("ENU"));
					}
				}

			}

			if (estimatorType == 17) {
				String estType = "TDCP-CSDR";
				// TreeMap<String, ArrayList<double[]>> estTdcpVelMap = new TreeMap<String,
				// ArrayList<double[]>>();
				long prevTime = timeList.get(0);
				for (int i = 1; i < timeList.size(); i++) {
					System.out.println("\n\n Epoch : " + i);
					long currentTime = timeList.get(i);
					ArrayList<Satellite> currentSatList = satMap.get(currentTime);
					ArrayList<Satellite> prevSatList = satMap.get(prevTime);
					double[] refPos = LinearLeastSquare.getEstPos(currentSatList, true, useIGS);
					double[] estVel = LLS_TDCP_ambFix.getEstVel(currentSatList, prevSatList, true, doAnalyze, refPos,
							useIGS, currentTime);
					estVelMap.computeIfAbsent(estType, k -> new ArrayList<double[]>()).add(estVel);
					prevTime = currentTime;
					if (doAnalyze && estimatorType != 11) {
						double tRx = currentTime / 1000.0;
						double[] residual = LLS_TDCP_ambFix.getResidual();
						SimpleMatrix Cyy = LLS_TDCP_ambFix.getCyy();
						satResMap
								.computeIfAbsent(Measurement.TDCP,
										k -> new HashMap<String, HashMap<String, ArrayList<SatResidual>>>())
								.computeIfAbsent(estType, k -> new HashMap<String, ArrayList<SatResidual>>());
						ArrayList<CycleSlipDetect> csdList = LLS_TDCP_ambFix.getCsdList();
						int n = csdList.size();
						for (int j = 0; j < n; j++) {
							CycleSlipDetect csd = csdList.get(j);
							Satellite sat = csd.getSat();
							satResMap.get(Measurement.TDCP).get(estType)
									.computeIfAbsent(sat.getObsvCode().charAt(0) + "" + sat.getSvid(),
											k -> new ArrayList<SatResidual>())
									.add(new SatResidual(tRx - tRx0, sat.getElevAzm()[0], residual[j], sat.isOutlier(),
											Math.sqrt(Cyy.get(j, j)), sat.getCn0DbHz()));

						}

						satCountMap.computeIfAbsent(Measurement.TDCP, k -> new TreeMap<String, ArrayList<Long>>())
								.computeIfAbsent(estType, k -> new ArrayList<Long>()).add((long) n);
						postVarOfUnitWeightMap
								.computeIfAbsent(Measurement.TDCP, k -> new HashMap<String, ArrayList<Double>>())
								.computeIfAbsent(estType, k -> new ArrayList<Double>())
								.add(LLS_TDCP_ambFix.getPostVarOfUnitW());
						State state = State.Velocity;
						Cxx_hat_map.computeIfAbsent(state, k -> new HashMap<String, ArrayList<SimpleMatrix>>())
								.computeIfAbsent(estType, k -> new ArrayList<SimpleMatrix>())
								.add(LLS_TDCP_ambFix.getCxx_hat("ENU"));
					}
				}
				GraphPlotter.graphAmbiguityCount(LLS_TDCP_ambFix.getAmbDetectedCountMap(),
						LLS_TDCP_ambFix.getAmbRepairedCountMap(), timeList);
				System.out.println("Ambiguity Detected Count: " + LLS_TDCP_ambFix.getAmbDetectedCount());
				System.out.println("Ambiguity Repaired Count: " + LLS_TDCP_ambFix.getAmbRepairedCount());
				System.out.println("Ambiguity Repair Percentage: "
						+ ((LLS_TDCP_ambFix.getAmbRepairedCount() * 100.0) / LLS_TDCP_ambFix.getAmbDetectedCount()));
			}
			if (estimatorType == 18 || estimatorType == 19 || estimatorType == 20) {
				HashMap<String, ArrayList<CycleSlipDetect>> satCSmap = new HashMap<String, ArrayList<CycleSlipDetect>>();
				String estName = "EKF TDCP";
				boolean innPhaseRate = false;
				boolean onlyDoppler = false;
				if (estimatorType == 19) {
					estName += " Innov Phase-Rate";
					innPhaseRate = true;
				} else if (estimatorType == 20) {
					estName = " Doppler EKF";
					onlyDoppler = true;
				}
				EKF_TDCP_ambFix2 ekf = new EKF_TDCP_ambFix2();
				TreeMap<Long, double[]> estStateMap = ekf.process(satMap, timeList, useIGS, obsvCodeList, doAnalyze,
						doTest, outlierAnalyze, innPhaseRate, onlyDoppler, trueEcefList);

				useDoppler = true;
				int n = timeList.size();
				estVelMap.put(estName, new ArrayList<double[]>());
				HashMap<Measurement, HashMap<String, HashMap<String, ArrayList<SatResidual>>>> satInnMap = new HashMap<Measurement, HashMap<String, HashMap<String, ArrayList<SatResidual>>>>();
				if (doAnalyze && estimatorType != 11) {
					satResMap.put(Measurement.TDCP, new HashMap<String, HashMap<String, ArrayList<SatResidual>>>());
					satResMap.get(Measurement.TDCP).put(estName, new HashMap<String, ArrayList<SatResidual>>());
					satInnMap.put(Measurement.TDCP, new HashMap<String, HashMap<String, ArrayList<SatResidual>>>());
					satInnMap.get(Measurement.TDCP).put(estName, new HashMap<String, ArrayList<SatResidual>>());
					satCountMap.put(Measurement.TDCP, new TreeMap<String, ArrayList<Long>>());
					Cxx_hat_map.put(State.Velocity, new HashMap<String, ArrayList<SimpleMatrix>>());
					postVarOfUnitWeightMap.put(Measurement.TDCP, new HashMap<String, ArrayList<Double>>());

				}
				for (int i = 0; i < n; i++) {
					long time = timeList.get(i);
					double[] estVel = estStateMap.get(time);

					if (estVel != null) {
						estVelMap.get(estName).add(estVel);
					} else {
						continue;
					}
					if (doAnalyze && estimatorType != 11) {
						ArrayList<Satellite> satList = ekf.getSatListMap().get(time);
						double[] residual = ekf.getResidualMap().get(time);
						double[] innovation = ekf.getInnovationMap().get(time);
						ArrayList<CycleSlipDetect> csdList = ekf.getCsdListMap().get(time);
						int m = satList.size();
						double tRx = time / 1000.0;
						// double[] measNoise = ekf.getMeasNoiseMap().get(time);
						for (int j = 0; j < m; j++) {
							Satellite sat = satList.get(j);
							satResMap.get(Measurement.TDCP).get(estName)
									.computeIfAbsent(sat.getObsvCode().charAt(0) + "" + sat.getSvid(),
											k -> new ArrayList<SatResidual>())
									.add(new SatResidual(tRx - tRx0, sat.getElevAzm()[0], residual[j], sat.isOutlier(),
											sat.getCn0DbHz()));

							satInnMap.get(Measurement.TDCP).get(estName)
									.computeIfAbsent(sat.getObsvCode().charAt(0) + "" + sat.getSvid(),
											k -> new ArrayList<SatResidual>())
									.add(new SatResidual(tRx - tRx0, sat.getElevAzm()[0], innovation[j],
											sat.isOutlier(), sat.getCn0DbHz()));

						}
						for (int j = 0; j < csdList.size(); j++) {
							CycleSlipDetect csdObj = csdList.get(j);
							Satellite sat = csdObj.getSat();
							String obsvCode = sat.getObsvCode();
							for (int k = 0; k < obsvCodeList.length; k++) {
								if (obsvCodeList[k].equals(obsvCode)) {
									csdObj.setClkDrift(estVel[3 + k]);
								}
							}
							satCSmap.computeIfAbsent(sat.getObsvCode().charAt(0) + "" + sat.getSvid(),
									k -> new ArrayList<CycleSlipDetect>()).add(csdObj);
						}
						satCountMap.get(Measurement.TDCP).computeIfAbsent(estName, k -> new ArrayList<Long>())
								.add(ekf.getSatCountMap().get(time));
						Cxx_hat_map.get(State.Velocity).computeIfAbsent(estName, k -> new ArrayList<SimpleMatrix>())
								.add(ekf.getErrCovMap().get(time));
						postVarOfUnitWeightMap.get(Measurement.TDCP)
								.computeIfAbsent(estName, k -> new ArrayList<Double>())
								.add(ekf.getPostVarOfUnitWMap().get(time));

						// For innovation vector

					}
				}
				if (doAnalyze && estimatorType != 11) {
					GraphPlotter.graphSatRes(satInnMap, outlierAnalyze, true);
					if (estimatorType == 18 || estimatorType == 19) {
						GraphPlotter.graphCycleSlip(satCSmap);
					}
				}
				if (!onlyDoppler) {
					GraphPlotter.graphAmbiguityCount(ekf.getAmbDetectedCountMap(), ekf.getAmbRepairedCountMap(),
							timeList);
					System.out.println("Ambiguity Detected Count: " + ekf.getAmbDetectedCount());
					System.out.println("Ambiguity Repaired Count: " + ekf.getAmbRepairedCount());
					System.out.println("Ambiguity Repair Percentage: "
							+ ((ekf.getAmbRepairedCount() * 100.0) / ekf.getAmbDetectedCount()));
				}
			}
			if (estimatorType == 22) {
				HashMap<Measurement, HashMap<String, HashMap<String, ArrayList<SatResidual>>>> satInnMap = new HashMap<Measurement, HashMap<String, HashMap<String, ArrayList<SatResidual>>>>();
				EKF_PPP ekf = new EKF_PPP();
				TreeMap<Long, double[]> estStateMap = ekf.process(satMap, timeList, obsvCodeList, doAnalyze, doTest,
						trueEcefList, false, repairCS, false);
				int n = timeList.size();
				estPosMap.put("PPP", new ArrayList<double[]>());
				estVelMap.put("PPP", new ArrayList<double[]>());
				for (int i = 0; i < n; i++) {
					long time = timeList.get(i);
					double[] estState = estStateMap.get(time);
					if (estState == null) {
						estPosMap.get("PPP").add(null);
						continue;
					}
					double[] estPos = new double[] {estState[0],estState[1],estState[2]};
					double[] estVel = new double[] {estState[3],estState[4],estState[5]};
					estPosMap.get("PPP").add(estPos);
					estVelMap.get("PPP").add(estVel);
				}
				HashMap<String, int[]> csCountMap = ekf.getCycleSlipCount();
				System.out.println("These satellite have more than 40% phase data with Cycle Slips");
				for (String satID : csCountMap.keySet()) {
					int[] csCount = csCountMap.get(satID);
					double percentage = (csCount[0] * 1.0) / csCount[1];
					if (percentage > 0.4) {
						System.out.print(satID + ", ");
					}
				}
				System.out.println();
				for (String satID : csCountMap.keySet()) {
					int[] csCount = csCountMap.get(satID);
					System.out.print(satID + " : " + csCount[0] + "/" + csCount[1] + " , ");

				}
				System.out.println();
				HashMap<Measurement, HashMap<String, ArrayList<Double>>> RedundancyNoMap = new HashMap<Measurement, HashMap<String, ArrayList<Double>>>();
				if (doAnalyze) {
					for (Measurement meas : List.of(Measurement.Pseudorange, Measurement.CarrierPhase,
							Measurement.Doppler, Measurement.GIM_Iono)) {
						satResMap.put(meas, new HashMap<String, HashMap<String, ArrayList<SatResidual>>>());
						satResMap.get(meas).put("PPP", new HashMap<String, ArrayList<SatResidual>>());
						satInnMap.put(meas, new HashMap<String, HashMap<String, ArrayList<SatResidual>>>());
						satInnMap.get(meas).put("PPP", new HashMap<String, ArrayList<SatResidual>>());
						satCountMap.put(meas, new TreeMap<String, ArrayList<Long>>());
						postVarOfUnitWeightMap.put(meas, new HashMap<String, ArrayList<Double>>());
						RedundancyNoMap.put(meas, new HashMap<String, ArrayList<Double>>());

					}
					dopMap.put("PPP", new ArrayList<double[]>());

				}
				for (int i = 1; i < n; i++) {
					long time = timeList.get(i);

					if (doAnalyze) {
						ArrayList<Satellite> satList = (ArrayList<Satellite>) ekf.getSatListMap().get(time);
						Map<Measurement, double[]> residualMap = (Map<Measurement, double[]>) ekf.getResidualMap()
								.get(time);
						Map<Measurement, double[]> innovationMap = (Map<Measurement, double[]>) ekf.getInnovationMap()
								.get(time);
						Map<Measurement, Double> postVarOfUnitWMap = (Map<Measurement, Double>) ekf
								.getPostVarOfUnitWMap().get(time);
						Map<Measurement, Double> redunMap = ekf.getRedundancyNoMap().get(time);
						int m = satList.size();
						long tRx = time / 1000;
						for (Measurement meas : List.of(Measurement.Pseudorange, Measurement.CarrierPhase,
								Measurement.Doppler, Measurement.GIM_Iono)) {
							int size = residualMap.get(meas).length;
							for (int j = 0; j < size; j++) {
								Satellite sat = satList.get(j);

								satResMap.get(meas).get("PPP")
										.computeIfAbsent(sat.getObsvCode() + "" + sat.getSvid(),
												k -> new ArrayList<SatResidual>())
										.add(new SatResidual(tRx - tRx0, sat.getElevAzm()[0], residualMap.get(meas)[j],
												sat.isOutlier(), sat.getCn0DbHz()));
								satInnMap.get(meas).get("PPP")
										.computeIfAbsent(sat.getObsvCode() + "" + sat.getSvid(),
												k -> new ArrayList<SatResidual>())
										.add(new SatResidual(tRx - tRx0, sat.getElevAzm()[0],
												innovationMap.get(meas)[j], sat.isOutlier(), sat.getCn0DbHz()));

							}
						}
						for (Measurement meas : List.of(Measurement.Pseudorange, Measurement.CarrierPhase,
								Measurement.Doppler, Measurement.GIM_Iono)) {

							postVarOfUnitWeightMap.get(meas).computeIfAbsent("PPP", k -> new ArrayList<Double>())
									.add(postVarOfUnitWMap.get(meas));

							RedundancyNoMap.get(meas).computeIfAbsent("PPP", k -> new ArrayList<Double>())
									.add(redunMap.get(meas));

						}
						satCountMap.get(Measurement.Pseudorange).computeIfAbsent("PPP", k -> new ArrayList<Long>())
								.add(ekf.getSatCountMap().get(time));
						dopMap.get("PPP").add(ekf.getDopMap().get(time));

					}
				}
				if (doAnalyze) {

					String[] ssiLabel = (String[]) ssiSet.stream().map(String::valueOf).toArray(String[]::new);
					GraphPlotter.graphSatRes(satInnMap, outlierAnalyze, true);
					GraphPlotter.graphRedundancyPPP(RedundancyNoMap, timeList);
					//GraphPlotter.createPPPplots(ekf, obsvCodeList, ssiLabel, timeList.get(0));
				}
			}
			TreeMap<Long, HashMap<AndroidSensor, IMUsensor>> imuMap = null;
			// TreeMap<Long, HashMap<AndroidSensor, IMUsensor>> imuMap =
			// IMUconfigure.configure(timeList.get(0), 100,
			// imuList);

			// Calculate Accuracy Metrics
			HashMap<String, ArrayList<double[]>> GraphPosMap = new HashMap<String, ArrayList<double[]>>();
			HashMap<String, ArrayList<double[]>> GraphVelMap = new HashMap<String, ArrayList<double[]>>();
			TreeMap<String, ArrayList<double[]>> trajectoryPosMap = new TreeMap<String, ArrayList<double[]>>();
			TreeMap<String, ArrayList<double[]>> trajectoryVelMap = new TreeMap<String, ArrayList<double[]>>();
			trajectoryPosMap.put("True", new ArrayList<double[]>());
			trajectoryVelMap.put("True", new ArrayList<double[]>());
			for (String key : estPosMap.keySet()) {
				ArrayList<Double>[] posErrList = new ArrayList[6];
				trajectoryPosMap.put(key, new ArrayList<double[]>());
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

					if (i == n - 1) {

						System.out.println("Converged Position RMS:");
						// error in East direction
						System.out.println("E  - " + Math.sqrt(enu[0] * enu[0]));
						// error in North direction
						System.out.println("N  - " + Math.sqrt(enu[1] * enu[1]));
						// error in Up direction
						System.out.println("U  - " + Math.sqrt(enu[2] * enu[2]));
						// 3d error
						System.out.println("3d Error - " + Math.sqrt(Arrays.stream(enu).map(j -> j * j).sum()));
						// 2d error
						System.out.println("2d Error - " + Math.sqrt((enu[0] * enu[0]) + (enu[1] * enu[1])));
					}


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
				System.out.println(" Haversine Distance - " + MathUtil.RMS(posErrList[5]));

				// 95th Percentile

				// IntStream.range(0, 6).forEach(i -> Collections.sort(posErrList[i]));
				// int q95 = (int) (n * 0.95);
				//
				// System.out.println("\n" + key + " 95%");
				// System.out.println("RMS - ");
				// System.out.println(" E - " + posErrList[0].get(q95));
				// System.out.println(" N - " + posErrList[1].get(q95));
				// System.out.println(" U - " + posErrList[2].get(q95));
				// System.out.println(" 3d Error - " + posErrList[3].get(q95));
				// System.out.println(" 2d Error - " + posErrList[4].get(q95));
				// System.out.println(" Haversine distance - " + posErrList[5].get(q95));

			}
			Set<Integer> tdcpEstSet = new HashSet<Integer>(Arrays.asList(16, 17, 18, 19, 20,22));
			if (tdcpEstSet.contains(estimatorType)) {
				TreeMap<Long, double[]> trueOrgVelEcef = Analyzer.getOriginalVel(trueEcefList, timeList);
				for (String key : estVelMap.keySet()) {
					trajectoryVelMap.put(key, new ArrayList<double[]>());
					ArrayList<Double>[] velErrList = new ArrayList[6];
					ArrayList<double[]> estVelList = null;
					ArrayList<double[]> enuVelList = null;
					IntStream.range(0, 5).forEach(i -> velErrList[i] = new ArrayList<Double>());
					enuVelList = new ArrayList<double[]>();
					estVelList = estVelMap.get(key);
					int n = trueOrgVelEcef.size();
					if (estVelList.size() != n) {
						System.err.println("FATAL ERROR: EST and TRUE ecef list size does not match ");
						throw new Exception("FATAL ERROR: EST and TRUE ecef list size does not match ");
					}

					for (int i = 0; i < n; i++) {
						double[] estVel = estVelList.get(i);
						long time = timeList.get(i + 1);
						double[] trueVel = trueOrgVelEcef.get(time);
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
					System.out.println(" E - " + MathUtil.RMS(velErrList[0]));
					System.out.println(" N - " + MathUtil.RMS(velErrList[1]));
					System.out.println(" U - " + MathUtil.RMS(velErrList[2]));
					System.out.println(" 3d Error - " + MathUtil.RMS(velErrList[3]));
					System.out.println(" 2d Error - " + MathUtil.RMS(velErrList[4]));

				}

			} else {
				for (String key : estVelMap.keySet()) {
					trajectoryVelMap.put(key, new ArrayList<double[]>());
					ArrayList<Double>[] velErrList = new ArrayList[6];
					ArrayList<double[]> estVelList = null;
					ArrayList<double[]> enuVelList = null;
					IntStream.range(0, 5).forEach(i -> velErrList[i] = new ArrayList<Double>());
					enuVelList = new ArrayList<double[]>();
					estVelList = estVelMap.get(key);
					int n = trueEcefList.size();
					ArrayList<Integer> removalList = new ArrayList<Integer>();
					for (int i = 0; i < n; i++) {
						double[] estVel = estVelList.get(i);
						long time = timeList.get(i);
						if (estVel == null || !trueVelEcef.containsKey(time)) {
							if (!trueVelEcef.containsKey(time) && (estVel != null)) {
								enuVelList.add(new double[3]);
								removalList.add(i);
							}
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
					// for (int i = removalList.size() - 1; i >= 0; i--) {
					// Cxx_hat_map.get(State.Velocity).get(key).remove((int) removalList.get(i));
					// }

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

					// IntStream.range(0, 5).forEach(i -> Collections.sort(velErrList[i]));
					// int q95 = (int) (n * 0.95);
					//
					// System.out.println("\n" + key + " 95%");
					// System.out.println("RMS - ");
					// System.out.println(" E - " + velErrList[0].get(q95));
					// System.out.println(" N - " + velErrList[1].get(q95));
					// System.out.println(" U - " + velErrList[2].get(q95));
					// System.out.println(" 3d Error - " + velErrList[3].get(q95));
					// System.out.println(" 2d Error - " + velErrList[4].get(q95));

				}
			}
			for (int i = 0; i < trueEcefList.size(); i++) {
				long time = timeList.get(i);
				double[] trueVel = new double[] { -999, -999, -999 };
				if (trueVelEcef.containsKey(time)) {
					trueVel = LatLonUtil.ecef2enu(trueVelEcef.get(time), trueEcefList.get(i), false);
				}
				trajectoryVelMap.get("True").add(trueVel);
				double[] trueTrajPos = LatLonUtil.ecef2enu(trueEcefList.get(i), trueEcefList.get(0), true);
				trajectoryPosMap.get("True").add(trueTrajPos);
				for (String key : estPosMap.keySet()) {
					double[] estPos = estPosMap.get(key).get(i);
					if (estPos != null) {
						trajectoryPosMap.get(key).add(LatLonUtil.ecef2enu(estPos, trueEcefList.get(0), true));
					} else {
						trajectoryPosMap.get(key).add(new double[] { -999, -999, -999 });
					}

				}

//				for (String key : estVelMap.keySet()) {
//					double[] estVel = estVelMap.get(key).get(i);
//					if (estVel != null) {
//
//						trajectoryVelMap.get(key).add(LatLonUtil.ecef2enu(estVel, trueEcefList.get(i), false));
//					} else {
//						trajectoryVelMap.get(key).add(new double[] { -999, -999, -999 });
//					}
//
//				}
			}
			// Trajectory.createCSV(trajectoryPosMap, trajectoryVelMap, path,
			// trueEcefList.size());

			System.out.println("\n\nPost Variance of Unit Weight Calculations");
			for (Measurement meas : postVarOfUnitWeightMap.keySet()) {
				System.out.println(meas.toString());
				for (String est_type : postVarOfUnitWeightMap.get(meas).keySet()) {
					System.out.println(est_type);
					ArrayList<Double> data = new ArrayList<Double>(postVarOfUnitWeightMap.get(meas).get(est_type));
					double sum = 0;
					int count = 0;
					for (int i = 0; i < data.size(); i++) {
						double val = data.get(i);
						if (val == 0 || val == -1) {
							continue;
						}
						sum += val;
						count++;

					}
					Collections.sort(data);
					double avg = sum / count;
					int q50 = (int) (count * 0.50);
					double median = data.get(q50);
					int _q75 = (int) (count * 0.75);
					double q75 = data.get(_q75);
					System.out.println("MEAN : " + avg);
					System.out.println("MEDIAN : " + median);
					System.out.println("Q75 : " + q75);

				}
			}

			long t0 = (long) (timeList.get(0) * 1e-3);
			long ctr = 0;
			System.out.println();
			for (int i = 0; i < timeList.size(); i++) {
				long t = (long) (timeList.get(i) * 1e-3);
				timeList.set(i, t - t0);

				if (i != (t - t0) - ctr) {
					ctr = (t - t0) - i;
					System.out.println("t-t0:" + (t - t0) + " i:" + i + " ctr:" + ctr);
				}
			}
			if (mapDeltaRanges) {
				GraphPlotter.graphDeltaRange(satMap, trueEcefList);
				GraphPlotter.graphTrajectory(trajectoryPosMap, trajectoryVelMap, trueEcefList.size());
			} else {
				GraphPlotter.graphTrajectory(trajectoryPosMap, trajectoryVelMap, trueEcefList.size());
				// Plot Error Graphs
				if (Cxx_hat_map.isEmpty()) {
					GraphPlotter.graphENU(GraphPosMap, timeList, true);
					GraphPlotter.graphENU(GraphVelMap, timeList, false);
				} else {
					GraphPlotter.graphENU(GraphPosMap, timeList, true, Cxx_hat_map.get(State.Position));
					GraphPlotter.graphENU(GraphVelMap, timeList, false, Cxx_hat_map.get(State.Velocity));
				}
				if (doAnalyze && estimatorType != 11) {
					GraphPlotter.graphSatRes(satResMap, outlierAnalyze);
					GraphPlotter.graphPostUnitW(postVarOfUnitWeightMap, timeList);
					// GraphPlotter.graphDOP(dopMap,satCountMap.get(Measurement.Pseudorange).get("WLS"),
					// timeList,1);
					GraphPlotter.graphSatCount(satCountMap, timeList, 1);

				}
			}
			if (doAnalyze && !tdcpEstSet.contains(estimatorType)) {
				Analyzer.processAndroid(satMap, imuMap, trueEcefList, trueVelEcef, estPosMap, estVelMap, satResMap,
						outlierAnalyze, useDoppler);
			}

		} catch (

		Exception e) {
			// TODO: handle exception
			System.out.println(e.getMessage());
			e.printStackTrace();
		}
	}

	public static void filterSat(ArrayList<Satellite> satList, double cutOffAng, double snrMask, double[] refEcef,
			boolean useIGS, IONEX ionex, Calendar time, int estimatorType) {
		if (cutOffAng >= 0) {
			satList.removeIf(i -> i.getElevAzm()[0] < Math.toRadians(cutOffAng));
		}
		if (snrMask >= 0) {
			satList.removeIf(i -> i.getCn0DbHz() < snrMask);
		}
		if (useIGS) {
			double[] refLatLon = LatLonUtil.ecef2lla(refEcef);
			// Geocentric Latitude
			double gcLat = LatLonUtil.gd2gc(refLatLon[0], refLatLon[2]);
			int n = satList.size();
			ComputeTropoCorr tropo = new ComputeTropoCorr(refLatLon, time, geoid);
			for (int i = 0; i < n; i++) {
				Satellite sat = satList.get(i);
				double[] eleAzm = sat.getElevAzm();
				double ionoErr = 0;
				double tropoErr = 0;

				ionoErr = ionex.computeIonoCorr(eleAzm[0], eleAzm[1], gcLat, refLatLon[1], sat.gettRx(),
						sat.getCarrierFrequencyHz(), time);

				double[] tropoParam = tropo.getSlantDelay(eleAzm[0]);
				tropoErr = tropoParam[0];
				double wetMF = tropoParam[1];

				Set<Integer> tdcpEstSet = new HashSet<Integer>(Arrays.asList(15, 16, 17, 18, 19, 20,22));
				if (tdcpEstSet.contains(estimatorType)) {
					sat.setIonoErr(ionoErr);
					ionoErr = 0;
				}
				sat.setPseudorange(sat.getPseudorange() - ionoErr - tropoErr);
				sat.setPhase(sat.getPhase() + ionoErr - tropoErr);
				sat.setWetMF(wetMF);
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
		File orekitData = new File("/Users/naman.agarwal/Documents/orekit/orekit-data-master/orekit-data-master");
		DataProvidersManager manager = DataContext.getDefault().getDataProvidersManager();
		manager.addProvider(new DirectoryCrawler(orekitData));
		
		NormalizedSphericalHarmonicsProvider nhsp = GravityFieldFactory.getNormalizedProvider(50, 50);
		Frame frame = FramesFactory.getITRF(ITRFVersion.ITRF_2014, IERSConventions.IERS_2010, true);

		// ReferenceEllipsoid refElp = new ReferenceEllipsoid(ae, f, frame, GM, spin);
		Geoid geoid = new Geoid(nhsp, ReferenceEllipsoid.getWgs84(frame));
		return geoid;

	}

}
