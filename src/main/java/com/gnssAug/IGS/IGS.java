package com.gnssAug.IGS;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.time.ZonedDateTime;
import java.time.ZoneOffset;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
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

import com.gnssAug.Rinex.constants.GnssDataConfig;
import com.gnssAug.Android.constants.EstimatorMode;
import com.gnssAug.Android.constants.Measurement;
import com.gnssAug.Android.constants.State;
import com.gnssAug.Rinex.estimation.EKF_PPP;
import com.gnssAug.Rinex.estimation.EKF_PPP_DF;
import com.gnssAug.Rinex.estimation.LinearLeastSquare;
import com.gnssAug.Rinex.fileParser.Antenna;
import com.gnssAug.Rinex.fileParser.DCB_Bias;
import com.gnssAug.Rinex.fileParser.Clock;
import com.gnssAug.Rinex.fileParser.IONEX;
import com.gnssAug.Rinex.fileParser.NavigationRNX;
import com.gnssAug.Rinex.fileParser.OSB_Bias;
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
import com.gnssAug.utility.MakeCSV;

/**
 * Top-level pipeline for IGS/geodetic receiver positioning using RINEX observation files
 * and MGEX precise products.
 *
 * <p>Supported estimator modes (set via {@link EstimatorMode}):
 * <ul>
 *   <li>{@code LLS_CODE} / {@code WLS_CODE} / {@code LLS_WLS_BOTH} — single-epoch
 *       least-squares (unweighted or weighted) on pseudorange only.</li>
 *   <li>{@code RINEX_EKF_CODE} — pseudorange EKF for smoothed code-only positioning.</li>
 *   <li>{@code IGS_PPP_FLOAT} — full UU-PPP float solution via {@link EKF_PPP}.</li>
 *   <li>{@code IGS_PPP_AR} — UU-PPP with BSD integer ambiguity resolution via
 *       {@link EKF_PPP} (sets {@code fixAmb = true}).</li>
 * </ul>
 *
 * <p><b>Data flow:</b>
 * <ol>
 *   <li>Parse RINEX NAV, RINEX OBS (optionally SINEX for ARP), MGEX orbit/clock/bias/IONEX.</li>
 *   <li>For each epoch: call {@link SingleFreq#process} to compute satellite positions,
 *       apply clock/bias/wind-up/relativistic corrections → raw {@link Satellite} list.</li>
 *   <li>Call {@link #filterSat} to apply elevation/SNR masks and iono/tropo corrections.</li>
 *   <li>Pass the epoch-keyed satellite map to the selected estimator.</li>
 *   <li>Compute ENU position errors relative to ARP; print RMS, 95th-percentile, and
 *       post-variance of unit weight per observable group.</li>
 *   <li>Plot ENU time series, residuals, DOP, redundancy, and PPP diagnostic plots
 *       via {@link GraphPlotter}.</li>
 * </ol>
 *
 * <p><b>Output:</b> results are redirected to a {@code .txt} file under
 * {@code gnss_output/IGS_rinex_output/PhD thesis/<obsName>}.
 */
public class IGS {

	/**
	 * Singleton Orekit {@link Geoid} instance shared across all filterSat calls.
	 * Initialised once by {@link #buildGeoid()} during {@link #posEstimate}.
	 * Uses EGM96 (50×50) gravity coefficients over WGS-84.
	 */
	private static Geoid geoid = null;

	/**
	 * Full positioning pipeline for an IGS/geodetic receiver session.
	 *
	 * <p>Parses all input files, iterates over RINEX observation epochs, and
	 * runs the selected estimator. Results (position errors, residuals, DOP,
	 * post-variance of unit weight) are printed to stdout (redirected to file)
	 * and plotted via {@link GraphPlotter}.
	 *
	 * @param obs_path          path to RINEX 3 observation file ({@code .rnx})
	 * @param nav_path          path to RINEX 3 mixed navigation file ({@code .rnx})
	 * @param osb_bias_path     path to MGEX OSB bias file ({@code .BIA}) — satellite
	 *                          code and phase biases in SINEX BIAS format
	 * @param dcb_bias_path     path to MGEX DCB file ({@code .BSX / .BIA}) — differential
	 *                          code biases used for clock product compatibility
	 * @param clock_path        path to MGEX precise clock file ({@code .CLK})
	 * @param orbit_path        path to MGEX precise orbit file ({@code .SP3})
	 * @param ionex_path        path to IONEX GIM file ({@code .INX / .I}) — used when
	 *                          {@code useGIM = true} for ionospheric corrections
	 * @param sinex_path        path to IGS SINEX solution file ({@code .SNX}) — used when
	 *                          {@code useSNX = true} to read the precise ARP ground truth
	 * @param useBias           if true, load and apply OSB and DCB satellite biases
	 * @param useGIM            if true, load the IONEX GIM for iono corrections
	 * @param useIGS            if true, use SP3/CLK precise orbit and clock products;
	 *                          if false, fall back to broadcast navigation
	 * @param useSNX            if true, read ARP from the SINEX file as ground truth
	 * @param obsvCodeList      signals to process (e.g. {@code "G1C", "G5Q"})
	 * @param minSat            minimum number of satellites required to process an epoch
	 * @param cutOffAng         elevation cut-off angle in degrees (satellites below are removed)
	 * @param snrMask           minimum C/N0 in dB-Hz (satellites below are removed); 0 = no mask
	 * @param corrIono          if true, apply ionospheric corrections via GIM or Klobuchar model
	 * @param corrTropo         if true, apply full slant tropospheric delay (dry + wet mapping)
	 * @param estimatorMode     positioning algorithm to use (see {@link EstimatorMode})
	 * @param doAnalyze         if true, collect and plot residuals, DOP, ambiguities, and state
	 * @param doTest            if true, apply chi-squared global model test + w-test to outliers
	 * @param outlierAnalyze    if true, flag outlier satellites in residual plots
	 * @param discardSet        set of satellite IDs to permanently exclude (e.g. {@code "G10"})
	 * @param repairCS          if true, attempt integer cycle-slip repair with LAMBDA
	 * @param predictPhaseClock if true, enable phase-clock prediction in the PPP process model
	 */
	public static void posEstimate(String obs_path, String nav_path, String osb_bias_path, String dcb_bias_path,
			String clock_path, String orbit_path, String ionex_path, String sinex_path,
			boolean useBias, boolean useGIM, boolean useIGS, boolean useSNX,
			String[] obsvCodeList, int minSat, double cutOffAng, double snrMask, boolean corrIono, boolean corrTropo,
			EstimatorMode estimatorMode, boolean doAnalyze, boolean doTest, boolean outlierAnalyze,
			Set<String> discardSet, boolean repairCS, boolean predictPhaseClock) {
		try {
			TimeZone.setDefault(TimeZone.getTimeZone("UTC"));

			HashMap<String, ArrayList<double[]>> estPosMap = new HashMap<String, ArrayList<double[]>>();
			HashMap<String, ArrayList<double[]>> estVelMap = new HashMap<String, ArrayList<double[]>>();
			TreeMap<Long, ArrayList<Satellite>> satMap = new TreeMap<Long, ArrayList<Satellite>>();
			HashMap<String, HashMap<String, ArrayList<Double>>> satMeasNoiseMap = new HashMap<String, HashMap<String, ArrayList<Double>>>();

			HashMap<Measurement, HashMap<String, HashMap<String, ArrayList<SatResidual>>>> satResMap = new HashMap<Measurement, HashMap<String, HashMap<String, ArrayList<SatResidual>>>>();
			HashMap<Measurement, HashMap<String, HashMap<String, ArrayList<SatResidual>>>> satInnMap = new HashMap<Measurement, HashMap<String, HashMap<String, ArrayList<SatResidual>>>>();
			HashMap<Measurement, HashMap<String, ArrayList<Double>>> postVarOfUnitWeightMap = new HashMap<Measurement, HashMap<String, ArrayList<Double>>>();
			HashMap<State, HashMap<String, ArrayList<SimpleMatrix>>> Cxx_hat_map = new HashMap<State, HashMap<String, ArrayList<SimpleMatrix>>>();
			HashMap<String, ArrayList<double[]>> dopMap = new HashMap<String, ArrayList<double[]>>();
			HashMap<Measurement, TreeMap<String, ArrayList<Long>>> satCountMap = new HashMap<Measurement, TreeMap<String, ArrayList<Long>>>();
			ArrayList<Long> timeList = new ArrayList<Long>();
			DCB_Bias dcb_bias = null;
			OSB_Bias osb_bias = null;
			Orbit orbit = null;
			Clock clock = null;
			Antenna antenna = null;
			IONEX ionex = null;

			String base_path = "/Users/naman.agarwal/Library/CloudStorage/OneDrive-UniversityofCalgary/input_files";

			String antenna_path = base_path + "/complementary/igs20.atx";

			String antenna_csv_path = base_path + "/complementary/antenna_igs20.csv";
			String obsName = fileNameOf(obs_path).replaceAll("\\.[^.]+$", "");
			String[] _obsNameParts = obsName.split("_");
			String _siteID   = _obsNameParts.length > 0 ? _obsNameParts[0] : obsName;
			String _epochTag = _obsNameParts.length > 2 ? _obsNameParts[2] : obsName;
			String _signals  = String.join("-", obsvCodeList);
			String _mode     = estimatorMode.name();
			String runLabel  = _siteID + "_" + _epochTag + "_" + _signals + "_" + _mode;
			String path = "/Users/naman.agarwal/Library/CloudStorage/OneDrive-UniversityofCalgary/gnss_output/IGS_rinex_output/PPP_AR/" + runLabel+"_discarded";
			File output = new File(path + ".txt");
			PrintStream stream;

			try {
				stream = new PrintStream(output);
				System.setOut(stream);
			} catch (FileNotFoundException e) { // TODO Auto-generated catch block
				e.printStackTrace();
			}
			// header printed after rxARP is available (see below)
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
			printRunHeader(obs_path, obsvCodeList, discardSet, rxARP, estimatorMode, repairCS, predictPhaseClock,
					doTest, doAnalyze, cutOffAng, snrMask, dcb_bias_path, clock_path, orbit_path, ionex_path, osb_bias_path);

			if (useBias) {
				osb_bias = new OSB_Bias(osb_bias_path);
				dcb_bias = new DCB_Bias(dcb_bias_path);

			}
			if (useIGS) {

				orbit = new Orbit(orbit_path);

				clock = new Clock(clock_path, dcb_bias);
				antenna = new Antenna(antenna_csv_path);

			}
			if (useSNX && antenna != null) {
				String antennaTypeKey = (String) ObsvMsgComp.getOrDefault("antennaType", null);
				if (antennaTypeKey != null) {
					for (String code : obsvCodeList) {
						double[] existing = rxPCO.get(code);
						if (existing == null || (existing[0] == 0.0 && existing[1] == 0.0 && existing[2] == 0.0)) {
							char ssi = code.charAt(0);
							int freq = Integer.parseInt(String.valueOf(code.charAt(1)));
							double[] enu = antenna.getRxPCO_ENU(antennaTypeKey, ssi, freq);
							if (enu != null) {
								double[] apc = LatLonUtil.enu2ecef(enu, rxARP, true);
								rxPCO.put(code, new double[]{ apc[0] - rxARP[0], apc[1] - rxARP[1], apc[2] - rxARP[2] });
							} else {
								System.err.println("Rx PCO unavailable in both SINEX and ANTEX for " + code + " (antenna: " + antennaTypeKey + ")");
							}
						}
					}
				}
			}
			if (useGIM) {
				ionex = new IONEX(ionex_path);

			}
			double tRx0 = ObsvMsgs.get(0).getTRX();
			int gtIndex = 0;
			double[] refEcef = null;
			for (ObservationMsg obsvMsg : ObsvMsgs) {
				gtIndex++;
				double tRx = obsvMsg.getTRX();
				double dayTime = tRx % 86400;
				long weekNo = obsvMsg.getWeekNo();
				Calendar time = Time.getDate(tRx, weekNo, 0);
				if (Time.getGPSTime(time)[0] != tRx) {
					System.err.println("FATAL ERROR TIME calendar");
					return;
				}
				ArrayList<Satellite> satList = null;
				satList = SingleFreq.process(obsvMsg, NavMsgs, obsvCodeList, useIGS, useBias, ionoCoeff, dcb_bias,
						osb_bias, orbit, clock, antenna, tRx, weekNo, time, discardSet, refEcef);

				if (satList.size() < minSat) {
					System.err.println("Less than " + minSat + " satellites");
					continue;
				}
				refEcef = LinearLeastSquare.getEstPos(satList, rxPCO, false);

				satList.stream().forEach(i -> i.setElevAzm(ComputeEleAzm.computeEleAzm(rxARP, i.getSatEci())));
				filterSat(satList, rxARP, cutOffAng, snrMask, corrIono, corrTropo, ionex, ionoCoeff, time,
						estimatorMode);
				if (satList.size() < minSat) {
					System.err.println("Less than " + minSat + " satellites");
					continue;
				}
				long tRxMilli = (long) (tRx * 1000);

				if (estimatorMode == EstimatorMode.LLS_CODE || estimatorMode == EstimatorMode.WLS_CODE || estimatorMode == EstimatorMode.LLS_WLS_BOTH) {
					int[] arr;
					if (estimatorMode == EstimatorMode.LLS_CODE)       arr = new int[] { 1 };
					else if (estimatorMode == EstimatorMode.WLS_CODE)  arr = new int[] { 2 };
					else                                                arr = new int[] { 1, 2 };
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
													sat.isOutlier(), sat.getCNo()));

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
							if (type == Measurement.Pseudorange) {
							dopMap.computeIfAbsent(estType, k -> new ArrayList<double[]>())
									.add(LinearLeastSquare.getDop());
							}
						}
					}
				}

				satMap.put(tRxMilli, satList);

				timeList.add(tRxMilli);
			}
			if (estimatorMode == EstimatorMode.LLS_WLS_BOTH || estimatorMode == EstimatorMode.RINEX_EKF_CODE) {
				satInnMap = new HashMap<Measurement, HashMap<String, HashMap<String, ArrayList<SatResidual>>>>();
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
									.computeIfAbsent(sat.getObsvCode() + "" + sat.getSVID(),
											k -> new ArrayList<SatResidual>())
									.add(new SatResidual(tRx - tRx0, sat.getElevAzm()[0], residual[j], sat.isOutlier(),
											sat.getCNo()));

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
									.computeIfAbsent(sat.getObsvCode() + "" + sat.getSVID(),
											k -> new ArrayList<SatResidual>())
									.add(new SatResidual(tRx - tRx0, sat.getElevAzm()[0], innovation[j],
											sat.isOutlier(), sat.getCNo()));

						}

					}
					
				}
				if (doAnalyze) {
					GraphPlotter.graphSatRes(satInnMap, outlierAnalyze, true);
				}
			}
			EKF_PPP ekf = null;
			if (estimatorMode == EstimatorMode.IGS_PPP_FLOAT || estimatorMode == EstimatorMode.IGS_PPP_AR) {
				boolean fixAmb = (estimatorMode == EstimatorMode.IGS_PPP_AR);
				ekf = new com.gnssAug.Rinex.estimation.EKF_PPP();
				TreeMap<Long, double[]> estStateMap_pos = ekf.process(satMap, rxPCO, timeList, doAnalyze, doTest, obsvCodeList,
						rxARP, true, repairCS, fixAmb, predictPhaseClock);

				int n = timeList.size();
				HashMap<String, int[]> csCountMap = ekf.getCycleSlipCount();
				System.out.println("SATELLITES WITH >20% CYCLE SLIPS");
				int _col = 0;
				for (String satID : csCountMap.keySet()) {
					int[] csCount = csCountMap.get(satID);
					if ((csCount[0] * 1.0) / csCount[1] > 0.2) {
						System.out.print("  " + satID);
						if (++_col % 8 == 0) System.out.println();
					}
				}
				System.out.println();
				System.out.println("CS DETECTED PER SATELLITE  (detected / total)");
				_col = 0;
				for (String satID : csCountMap.keySet()) {
					int[] csCount = csCountMap.get(satID);
					System.out.printf("  %-6s: %3d/%-4d", satID, csCount[0], csCount[1]);
					if (++_col % 4 == 0) System.out.println();
				}
				System.out.println();
				HashMap<Measurement, HashMap<String, ArrayList<Double>>> RedundancyNoMap = new HashMap<Measurement, HashMap<String, ArrayList<Double>>>();
				if (doAnalyze) {
					for (Measurement meas : List.of(Measurement.Pseudorange, Measurement.CarrierPhase,
							Measurement.Doppler,Measurement.GIM_Iono)) {
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
					double[] estPos = estStateMap_pos.get(time);
					estPosMap.computeIfAbsent("PPP", k -> new ArrayList<double[]>()).add(estPos);
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
								Measurement.Doppler,Measurement.GIM_Iono)) {
						int size = residualMap.get(meas).length;
						for (int j = 0; j < size; j++) {
							Satellite sat = satList.get(j);
							
								satResMap.get(meas).get("PPP")
								.computeIfAbsent(sat.getObsvCode() + "" + sat.getSVID(),
										k -> new ArrayList<SatResidual>())
								.add(new SatResidual(tRx - tRx0, sat.getElevAzm()[0], residualMap.get(meas)[j], sat.isOutlier(),
										sat.getCNo()));
								satInnMap.get(meas).get("PPP")
								.computeIfAbsent(sat.getObsvCode() + "" + sat.getSVID(),
										k -> new ArrayList<SatResidual>())
								.add(new SatResidual(tRx - tRx0, sat.getElevAzm()[0], innovationMap.get(meas)[j], sat.isOutlier(),
										sat.getCNo()));
								
							}
						}
						for (Measurement meas : List.of(Measurement.Pseudorange, Measurement.CarrierPhase,
								Measurement.Doppler,Measurement.GIM_Iono)) {
							
							postVarOfUnitWeightMap.get(meas)
							.computeIfAbsent("PPP", k -> new ArrayList<Double>()).add(postVarOfUnitWMap.get(meas));
							
							
							RedundancyNoMap.get(meas)
							.computeIfAbsent("PPP", k -> new ArrayList<Double>()).add(redunMap.get(meas));
							
						}
						satCountMap.get(Measurement.Pseudorange).computeIfAbsent("PPP", k -> new ArrayList<Long>())
								.add(ekf.getSatCountMap().get(time));
						dopMap.get("PPP").add(ekf.getDopMap().get(time));					
						
					}
				}
				if(doAnalyze)
				{
					long _dur = timeList.get(n - 1) - timeList.get(0);
					System.out.println("DATA SUMMARY");
					System.out.printf("  Epochs processed : %d%n", n - 1);
					System.out.printf("  Duration         : %.1f s  (%d min %02d s)%n",
							_dur / 1000.0, _dur / 60000, (_dur % 60000) / 1000);
					System.out.println("CYCLE SLIP STATISTICS");
					System.out.println("  Soft CS Detected : " + ekf.getCsDetectedCount());
					System.out.println("  CS Repaired      : " + ekf.getCsRepairedCount());
					System.out.println("  Hard Resets      : " + ekf.getHardResetCount());
					System.out.println("    GF discrepancy   : " + ekf.getResetCount_gfTest());
					System.out.println("    Consecutive CS   : " + ekf.getResetCount_consecutive());
					ListOrderedSet ssiSet = new ListOrderedSet();
					for (int i = 0; i < obsvCodeList.length;i++) {
						ssiSet.add(obsvCodeList[i].charAt(0)+"");
					}
					String[] ssiLabel = (String[]) ssiSet.toArray(new String[0]);
					GraphPlotter.graphSatRes(satInnMap, outlierAnalyze, true);
					GraphPlotter.graphRedundancyPPP(RedundancyNoMap, timeList);
					GraphPlotter.createPPPplots(ekf, obsvCodeList,ssiLabel, timeList.get(0));
					long _t0s = (timeList.get(0)      / 1000) % 86400;
					long _tNs = (timeList.get(n - 1)  / 1000) % 86400;
					String _utcTag = String.format("%02d%02d-%02d%02d",
							_t0s / 3600, (_t0s % 3600) / 60,
							_tNs / 3600, (_tNs % 3600) / 60);
					String taggedPath = path + "_" + _utcTag;
					output.renameTo(new File(taggedPath + ".txt"));
					GraphPlotter.writeEpochJSONL(taggedPath + "_epochs.jsonl", ekf, timeList,
							obsvCodeList, useSNX ? rxARP : null, estStateMap_pos);
				}
			}
		
			
//			if (estimatorType == 1) {
//				Analyzer.processIGS(satMap, rxARP, rxPCO, estPosMap, estVelMap, satResMap, outlierAnalyze);
//			}
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
					
					if (i == n - 1) {

						System.out.println("LAST EPOCH ERROR [m]");
						System.out.println("  E  : " + Math.sqrt(enu[0] * enu[0]));
						System.out.println("  N  : " + Math.sqrt(enu[1] * enu[1]));
						System.out.println("  U  : " + Math.sqrt(enu[2] * enu[2]));
						System.out.println("  3D : " + Math.sqrt(Arrays.stream(enu).map(j -> j * j).sum()));
						System.out.println("  2D : " + Math.sqrt((enu[0] * enu[0]) + (enu[1] * enu[1])));
					}

				}

				GraphPosMap.put(key, enuPosList);

				System.out.println("\n" + key + " POSITION RMS [m]");
				System.out.println("  E  : " + MathUtil.RMS(posErrList[0]));
				System.out.println("  N  : " + MathUtil.RMS(posErrList[1]));
				System.out.println("  U  : " + MathUtil.RMS(posErrList[2]));
				System.out.println("  3D : " + MathUtil.RMS(posErrList[3]));
				System.out.println("  2D : " + MathUtil.RMS(posErrList[4]));

				IntStream.range(0, 5).forEach(i -> Collections.sort(posErrList[i]));
				int q95 = (int) (n * 0.95);
				System.out.println("\n" + key + " 95th PERCENTILE [m]");
				System.out.println("  E  : " + posErrList[0].get(q95));
				System.out.println("  N  : " + posErrList[1].get(q95));
				System.out.println("  U  : " + posErrList[2].get(q95));
				System.out.println("  3D : " + posErrList[3].get(q95));
				System.out.println("  2D : " + posErrList[4].get(q95));

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
			System.out.println("\nPOST VARIANCE OF UNIT WEIGHT");
			for (Measurement meas : postVarOfUnitWeightMap.keySet()) {
				for (String est_type : postVarOfUnitWeightMap.get(meas).keySet()) {
					ArrayList<Double> data = new ArrayList<Double>(postVarOfUnitWeightMap.get(meas).get(est_type));
					double sum = 0;
					int count = 0;
					for (int i = 0; i < data.size(); i++) {
						double val = data.get(i);
						if (val <= 0) continue;
						sum += val;
						count++;
					}
					Collections.sort(data);
					double avg = sum / count;
					int q50 = (int) (count * 0.50);
					double median = data.get(q50);
					int _q75 = (int) (count * 0.75);
					double q75 = data.get(_q75);
					System.out.println("  " + meas.toString() + " [" + est_type + "]");
					System.out.println("    MEAN   : " + avg);
					System.out.println("    MEDIAN : " + median);
					System.out.println("    Q75    : " + q75);
				}
			}
			long t0 = timeList.get(0);
			for (int i = 0; i < timeList.size(); i++) {
				timeList.set(i, (long) ((timeList.get(i) - t0) * 1e-3));
			}

			// Plot Error Graphs
			if (Cxx_hat_map.isEmpty()) {
				GraphPlotter.graphENU(GraphPosMap, timeList, true);
//				GraphPlotter.graphENU(GraphVelMap, timeList, false);
			} else {
				GraphPlotter.graphENU(GraphPosMap, timeList, true, Cxx_hat_map.get(State.Position));
				GraphPlotter.graphENU(GraphVelMap, timeList, false, Cxx_hat_map.get(State.Velocity));
			}

			if (doAnalyze) {

				GraphPlotter.graphSatRes(satResMap, outlierAnalyze);
				GraphPlotter.graphPostUnitW(postVarOfUnitWeightMap, timeList);
				GraphPlotter.graphSatCount(satCountMap, timeList, 1);
				GraphPlotter.graphDOP(dopMap, satCountMap.get(Measurement.Pseudorange).get("PPP"),timeList);
				
				boolean makeCSV = false;
				if(makeCSV)
				{
					
					MakeCSV.exportPPPresultsToCSV(satMap, ekf,GraphPosMap,
							 satInnMap,
							timeList,  obsvCodeList,
							rxARP);
				}
			}
		} catch (

		Exception e) {
			// TODO: handle exception
			System.out.println(e.getMessage());
			e.printStackTrace();
		}

	}

	/**
	 * Applies satellite quality masking and measurement corrections in place.
	 *
	 * <p><b>Step 1 — Quality masking:</b>
	 * <ul>
	 *   <li>Elevation mask: removes satellites below {@code cutOffAng} degrees.</li>
	 *   <li>C/N0 mask: removes satellites with carrier-to-noise density below {@code snrMask} dB-Hz.</li>
	 * </ul>
	 *
	 * <p><b>Step 2 — Ionospheric correction:</b> if {@code corrIono} is true:
	 * <ul>
	 *   <li>If a loaded {@link IONEX} GIM is provided, use it for a precise slant iono delay
	 *       (GIM-interpolated VTEC mapped via thin-shell model).</li>
	 *   <li>Otherwise fall back to the Klobuchar broadcast model ({@link ComputeIonoCorr})
	 *       using coefficients from the RINEX NAV header.</li>
	 * </ul>
	 *
	 * <p><b>Step 3 — Tropospheric correction:</b> if {@code corrTropo} is true, computes
	 * the full slant delay (dry + wet) via {@link ComputeTropoCorr} (Neil mapping function,
	 * pressure/temperature from the Orekit geoid model). The wet mapping function is stored
	 * on each satellite for later use as a partial derivative in the EKF wet-tropo residual.
	 *
	 * <p><b>PPP mode special handling:</b> when the estimator is in PPP mode
	 * ({@link EstimatorMode#isPPPMode()}), the ionospheric delay is stored on the satellite
	 * object ({@link Satellite#setIonoErr}) but NOT subtracted from the measurements.
	 * The EKF estimates the iono as a state variable; the pre-stored value serves only
	 * as a pseudo-observation initialisation. For non-PPP modes the iono is pre-corrected.
	 *
	 * <p><b>Sign conventions applied to measurements:</b>
	 * <ul>
	 *   <li>Pseudorange: {@code PR -= ionoErr + tropoErr}</li>
	 *   <li>Carrier phase: {@code CP += ionoErr - tropoErr}
	 *       (iono advances phase, tropo delays it)</li>
	 * </ul>
	 *
	 * @param satList       satellite list to filter and correct (modified in place)
	 * @param refEcef       approximate receiver ECEF position used for tropo/iono computation
	 * @param cutOffAng     elevation cut-off angle in degrees; negative value disables the mask
	 * @param snrMask       minimum C/N0 in dB-Hz; negative or zero disables the mask
	 * @param corrIono      if true, compute and apply (or store) ionospheric delay
	 * @param corrTropo     if true, compute and apply full slant tropospheric delay
	 * @param ionex         loaded IONEX GIM object, or {@code null} to use Klobuchar
	 * @param ionoCoeff     Klobuchar coefficients from the RINEX NAV header (fallback)
	 * @param time          UTC epoch as {@link Calendar}, used for iono/tropo model input
	 * @param estimatorMode current estimator mode, determines PPP vs. non-PPP iono handling
	 */
	public static void filterSat(ArrayList<Satellite> satList, double[] refEcef, double cutOffAng, double snrMask,
			boolean corrIono, boolean corrTropo, IONEX ionex, IonoCoeff ionoCoeff, Calendar time, EstimatorMode estimatorMode) {
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
			ComputeTropoCorr tropo = new ComputeTropoCorr(refLatLon, time, geoid);
			for (int i = 0; i < n; i++) {
				Satellite sat = satList.get(i);
				double[] eleAzm = sat.getElevAzm();
				double ionoErr = 0;
				double tropoErr = 0;
				double wetMF = 1;
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

					double[] tropoParam = tropo.getSlantDelay(eleAzm[0]);
					tropoErr = tropoParam[0];
					wetMF = tropoParam[1];
				}
				if (estimatorMode.isPPPMode()) {
					sat.setIonoErr(ionoErr);
					ionoErr = 0;
				}
				sat.setPseudorange(sat.getPseudorange() - ionoErr - tropoErr);
				sat.setPhase(sat.getPhase() + ionoErr - tropoErr);
				sat.setWetMF(wetMF);

			}
		}
	}

	/**
	 * Prints a structured run header to stdout, capturing all configuration that
	 * affects the result. Output is redirected to the session's {@code .txt} log file,
	 * so every saved result is self-documenting.
	 *
	 * <p>The header includes: UTC timestamp, git commit hash, dataset file names,
	 * signal list, discard set, SINEX ground truth (if available), estimator settings,
	 * MGEX product file names, and all {@link GnssDataConfig} filter prior variances.
	 */
	private static void printRunHeader(String obs_path, String[] obsvCodeList, Set<String> discardSet,
			double[] rxARP, EstimatorMode estimatorMode, boolean repairCS, boolean predictPhaseClock,
			boolean doTest, boolean doAnalyze, double cutOffAng, double snrMask, String dcb_bias_path,
			String clock_path, String orbit_path, String ionex_path, String osb_bias_path) {
		String sep = "================================================";
		System.out.println(sep);
		System.out.println("Timestamp  : " + ZonedDateTime.now(ZoneOffset.UTC)
				.format(DateTimeFormatter.ofPattern("yyyy-MM-dd'T'HH:mm:ss'Z'")));
		System.out.println("Git Commit : " + getGitHash());
		System.out.println(sep);
		// Parse receiver / antenna from RINEX header (first occurrence only)
		String siteMarker = "", recType = "", antType = "";
		try (java.util.Scanner hSc = new java.util.Scanner(new java.io.File(obs_path))) {
			while (hSc.hasNextLine()) {
				String hl = hSc.nextLine();
				if (hl.contains("END OF HEADER")) break;
				if (hl.contains("MARKER NAME")     && siteMarker.isEmpty())
					siteMarker = hl.substring(0, Math.min(60, hl.length())).trim();
				if (hl.contains("REC # / TYPE / VERS") && recType.isEmpty())
					recType = hl.length() >= 40 ? hl.substring(20, 40).trim() : hl.substring(20).trim();
				if (hl.contains("ANT # / TYPE")    && antType.isEmpty())
					antType = hl.length() >= 40 ? hl.substring(20, 40).trim() : hl.substring(20).trim();
			}
		} catch (Exception ignored) {}
		System.out.println("DATASET");
		System.out.println("  File     : " + fileNameOf(obs_path));
		System.out.println("  Site     : " + siteMarker);
		System.out.println("  Receiver : " + recType);
		System.out.println("  Antenna  : " + antType);
		System.out.println("  Signals  : " + String.join("  ", obsvCodeList));
		System.out.println("  Discard  : " + discardSet);
		if (rxARP != null) {
			double[] llh = LatLonUtil.ecef2lla(rxARP);
			System.out.println("GROUND TRUTH (SINEX ARP)");
			System.out.printf("  Lat / Lon / H : %.8f deg  %.8f deg  %.4f m%n", llh[0], llh[1], llh[2]);
			System.out.printf("  ECEF (m)      : %.4f  %.4f  %.4f%n", rxARP[0], rxARP[1], rxARP[2]);
		}
		System.out.println("ESTIMATOR");
		System.out.println("  Mode      : " + estimatorMode.name());
		System.out.println("  repairCS          : " + repairCS);
		System.out.println("  predictPhaseClock : " + predictPhaseClock);
		System.out.println("  doTest    : " + doTest);
		System.out.println("  doAnalyze : " + doAnalyze);
		System.out.println("  elMask    : " + cutOffAng + " deg   snrMask : " + snrMask + " dB-Hz");
		System.out.println("PRODUCTS");
		System.out.println("  DCB   : " + fileNameOf(dcb_bias_path));
		System.out.println("  CLK   : " + fileNameOf(clock_path));
		System.out.println("  ORB   : " + fileNameOf(orbit_path));
		System.out.println("  IONEX : " + fileNameOf(ionex_path));
		System.out.println("  OSB   : " + fileNameOf(osb_bias_path));
		System.out.println("FILTER PARAMETERS");
		System.out.println("  pseudorange_priorStdOfUnitW = " + Math.sqrt(GnssDataConfig.pseudorange_priorVarOfUnitW));
		System.out.println("  doppler_priorStdOfUnitW     = " + Math.sqrt(GnssDataConfig.doppler_priorVarOfUnitW));
		System.out.println("  TDCP_priorStdOfUnitW        = " + Math.sqrt(GnssDataConfig.tdcp_priorVarOfUnitW));
		System.out.println("  phase_priorStdOfUnitW       = " + Math.sqrt(GnssDataConfig.phase_priorVarOfUnitW));
		System.out.println("  GIM_TECU_priorStdOfUnitW    = " + Math.sqrt(GnssDataConfig.GIM_TECU_variance));
		System.out.println("  Q_pos_randWalk              = " + Arrays.toString(GnssDataConfig.qENU_posRandWalk));
		System.out.println("  Q_vel_randWalk              = " + Arrays.toString(GnssDataConfig.qENU_velRandWalk));
		System.out.println("  nSamplesMC                  = " + (long) GnssDataConfig.nSamplesMC);
		System.out.println(sep);
	}

	/** Extracts the filename (everything after the last {@code /}) from an absolute path. */
	private static String fileNameOf(String path) {
		if (path == null) return "N/A";
		return path.substring(path.lastIndexOf('/') + 1);
	}

	/**
	 * Returns the current git commit short hash, with a {@code " (dirty)"} suffix if
	 * there are uncommitted changes. Falls back to {@code "unknown"} if git is not
	 * available or the command fails. Used in the run header to tie each result log
	 * to an exact code version.
	 */
	private static String getGitHash() {
		try {
			Process p = new ProcessBuilder("git", "rev-parse", "--short", "HEAD")
					.redirectErrorStream(true).start();
			String hash = new BufferedReader(new InputStreamReader(p.getInputStream())).readLine();
			Process p2 = new ProcessBuilder("git", "status", "--porcelain")
					.redirectErrorStream(true).start();
			boolean dirty = p2.getInputStream().read() != -1;
			return (hash != null ? hash : "unknown") + (dirty ? " (dirty)" : "");
		} catch (Exception e) {
			return "unknown";
		}
	}

	/**
	 * Initialises and returns an Orekit {@link Geoid} object backed by a 50×50 EGM96
	 * gravity field over the WGS-84 ellipsoid, referenced to ITRF 2014 / IERS 2010.
	 *
	 * <p>The geoid is used by {@link ComputeTropoCorr} to obtain the orthometric height
	 * (mean sea level altitude) of the station, which enters the pressure/temperature
	 * model for tropospheric delay computation.
	 *
	 * <p>Requires the Orekit data master directory at
	 * {@code ~/Documents/orekit/orekit-data-master/orekit-data-master}.
	 * This method is called once during {@link #posEstimate} and the result is cached
	 * in the {@link #geoid} static field.
	 *
	 * @return a configured {@link Geoid} instance
	 */
	public static Geoid buildGeoid() {
		// WGS-84 defining parameters
		final double ae = 6378137;          // semi-major axis (m)
		final double f = 1 / 298.257223563; // flattening

		final double spin = 7.2921151467E-5;  // Earth's rotation rate (rad/s)
		final double GM   = 3.986004418E14;   // Earth's gravitational parameter (m³/s²)
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
