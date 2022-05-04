package com.gnssAug.Android;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collections;
import java.util.HashMap;
import java.util.TimeZone;
import java.util.TreeMap;
import java.util.stream.IntStream;

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

import com.gnssAug.Android.estimation.LinearLeastSquare;
import com.gnssAug.Android.fileParser.DerivedCSV;
import com.gnssAug.Android.fileParser.GNSS_Log;
import com.gnssAug.Android.fileParser.GroundTruth;
import com.gnssAug.Android.helper.ComputeEleAzm;
import com.gnssAug.Android.models.Derived;
import com.gnssAug.Android.models.GNSSLog;
import com.gnssAug.Android.models.Satellite;
import com.gnssAug.Android.utility.GraphPlotter;
import com.gnssAug.Android.utility.LatLonUtil;
import com.gnssAug.Android.utility.Time;

public class Android {
	public static void posEstimate(boolean doPosErrPlot, double cutOffAng, double snrMask, int estimatorType,
			String[] obsvCodeList, String derived_csv_path, String gnss_log_path, String GTcsv) {
		try {
			TimeZone.setDefault(TimeZone.getTimeZone("UTC"));
			HashMap<String, ArrayList<HashMap<String, Double>>> ErrMap = new HashMap<String, ArrayList<HashMap<String, Double>>>();

			ArrayList<Long> timeList = new ArrayList<Long>();
			ArrayList<double[]> trueLLHlist = new ArrayList<double[]>();
			ArrayList<double[]> trueECEFlist = new ArrayList<double[]>();

			HashMap<String, ArrayList<double[]>> estPosMap = new HashMap<String, ArrayList<double[]>>();

			String path = "C:\\Users\\Naman\\Desktop\\rinex_parse_files\\google2\\test";
			File output = new File(path + ".txt");
			PrintStream stream;
			stream = new PrintStream(output);
			System.setOut(stream);

			ArrayList<double[]> rxLLH = GroundTruth.processCSV(GTcsv);
			HashMap<Long, HashMap<String, HashMap<Integer, Derived>>> derivedMap = null;
			TreeMap<Long, HashMap<String, ArrayList<GNSSLog>>> gnssLogMaps = null;
			derivedMap = DerivedCSV.processCSV(derived_csv_path);
			gnssLogMaps = GNSS_Log.process(gnss_log_path);

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
				gtIndex++;

				Calendar time = Time.getDate(tRx, weekNo, 0);
				ArrayList<Satellite> satList = SingleFreq.process(tRx, derivedMap, gnssLogMap, time, obsvCodeList,
						weekNo);
				if (satList.size() < 4) {
					continue;
				}
				double[] userEcef = LinearLeastSquare.process(satList, false);
				satList.stream().forEach(i -> i.setElevAzm(ComputeEleAzm.computeEleAzm(userEcef, i.getSatEci())));
				filterSat(satList, cutOffAng, snrMask);
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

					break;
				case 3:
					// Implement LS method
					estEcefClk = LinearLeastSquare.process(satList, false);
					estPosMap.computeIfAbsent("LS", k -> new ArrayList<double[]>()).add(estEcefClk);
					// Implement WLS method
					estEcefClk = LinearLeastSquare.process(satList, true);
					estPosMap.computeIfAbsent("WLS", k -> new ArrayList<double[]>()).add(estEcefClk);

				}
				trueLLHlist.add(trueUserLLH);
				trueECEFlist
						.add(LatLonUtil.lla2ecef(new double[] { trueUserLLH[0], trueUserLLH[1], trueUserLLH[2] - 61 }));
				timeList.add(tRxMilli);
			}
			// Calculate Accuracy Metrics
			HashMap<String, ArrayList<double[]>> GraphEnuMap = new HashMap<String, ArrayList<double[]>>();
			for (String key : estPosMap.keySet()) {
				ArrayList<Double>[] errList = new ArrayList[6];
				IntStream.range(0, 6).forEach(i -> errList[i] = new ArrayList<Double>());
				ArrayList<double[]> estPosList = estPosMap.get(key);
				int n = estPosList.size();
				ArrayList<double[]> enuList = new ArrayList<double[]>();
				for (int i = 0; i < n; i++) {
					double[] estEcef = estPosList.get(i);
					double[] enu = LatLonUtil.ecef2enu(estEcef, trueECEFlist.get(i));
					double[] estLLH = LatLonUtil.ecef2lla(estEcef);
					// Great Circle Distance
					double gcErr = LatLonUtil.getHaversineDistance(estLLH, trueLLHlist.get(i));
					enuList.add(enu);
					// error in East direction
					errList[0].add(Math.sqrt(enu[0] * enu[0]));
					// error in North direction
					errList[1].add(Math.sqrt(enu[1] * enu[1]));
					// error in Up direction
					errList[2].add(Math.sqrt(enu[2] * enu[2]));
					// 3d error
					errList[3].add(Math.sqrt(Arrays.stream(enu).map(j -> j * j).sum()));
					// 2d error
					errList[4].add(Math.sqrt((enu[0] * enu[0]) + (enu[1] * enu[1])));
					// Haversine Distance
					errList[5].add(gcErr);
				}

				GraphEnuMap.put(key, enuList);

				// RMSE
				System.out.println("\n" + key);
				System.out.println("RMS - ");
				System.out.println(" E - " + RMS(errList[0]));
				System.out.println(" N - " + RMS(errList[1]));
				System.out.println(" U - " + RMS(errList[2]));
				System.out.println(" 3d Error - " + RMS(errList[3]));
				System.out.println(" 2d Error - " + RMS(errList[4]));
				System.out.println(" Haversine Distance - " + RMS(errList[5]));

				// 95th Percentile
				IntStream.range(0, 6).forEach(i -> Collections.sort(errList[i]));
				int q95 = (int) (n * 0.95);

				System.out.println("\n" + key + " 95%");
				System.out.println("RMS - ");
				System.out.println(" E - " + errList[0].get(q95));
				System.out.println(" N - " + errList[1].get(q95));
				System.out.println(" U - " + errList[2].get(q95));
				System.out.println(" 3d Error - " + errList[3].get(q95));
				System.out.println(" 2d Error - " + errList[4].get(q95));
				System.out.println(" Haversine distance - " + errList[5].get(q95));

			}

			// Plot Error Graphs
			GraphPlotter.graphENU(GraphEnuMap, timeList);

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
