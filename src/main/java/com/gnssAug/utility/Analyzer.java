package com.gnssAug.utility;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.stream.IntStream;

import org.jfree.ui.RefineryUtilities;

import com.gnssAug.Android.constants.AndroidSensor;
import com.gnssAug.Android.models.IMUsensor;
import com.gnssAug.Android.models.Satellite;

public class Analyzer {

	public static void process(TreeMap<Long, ArrayList<Satellite>> SatMap,
			TreeMap<Long, HashMap<AndroidSensor, IMUsensor>> imuMap, HashMap<Long, double[]> err) throws IOException {
		ArrayList<ArrayList<Satellite>> SVlist = new ArrayList<ArrayList<Satellite>>(SatMap.values());
		HashMap<String, TreeMap<Integer, Double>> doppplerMap = new HashMap<String, TreeMap<Integer, Double>>();
		HashMap<String, TreeMap<Integer, Double>> rangeMap = new HashMap<String, TreeMap<Integer, Double>>();

		HashMap<String, double[]> firstVal = new HashMap<String, double[]>();
		long time0 = SatMap.firstKey();
		for (Long time : SatMap.keySet()) {
			int timeDiff = (int) ((time - time0) / 1e3);
			ArrayList<Satellite> satList = SatMap.get(time);
			for (Satellite sat : satList) {
				double rangeRate = sat.getPseudorangeRateMetersPerSecond();
				double range = sat.getPseudorange() / 1000;
				int svid = sat.getSvid();
				String code = sat.getObsvCode().charAt(0) + "";
				double[] first = null;
				if (firstVal.containsKey(code + svid)) {
					first = firstVal.get(code + svid);
				} else {
					first = new double[] { range, rangeRate };
					firstVal.put(code + svid, first);
				}
				rangeMap.computeIfAbsent(code + svid, k -> new TreeMap<Integer, Double>()).put(timeDiff,
						range - first[0]);

				doppplerMap.computeIfAbsent(code + svid, k -> new TreeMap<Integer, Double>()).put(timeDiff,
						rangeRate - first[1]);
			}
			if (err.containsKey(time)) {
				doppplerMap.computeIfAbsent("speed", k -> new TreeMap<Integer, Double>()).put(timeDiff,
						err.get(time)[1]);
				rangeMap.computeIfAbsent("position", k -> new TreeMap<Integer, Double>()).put(timeDiff,
						err.get(time)[0]);
			}
		}
		HashMap<String, Double> dopplerFirst = new HashMap<String, Double>();
		HashMap<String, Double> rangeFirst = new HashMap<String, Double>();
		for (String key : firstVal.keySet()) {
			rangeFirst.put(key, firstVal.get(key)[0]);
			dopplerFirst.put(key, firstVal.get(key)[1]);
		}

		GraphPlotter chart = new GraphPlotter("Range-Rate(in m/s)", dopplerFirst, doppplerMap);
		chart.pack();
		RefineryUtilities.positionFrameRandomly(chart);
		chart.setVisible(true);

		chart = new GraphPlotter("Range(in Km)", rangeFirst, rangeMap);
		chart.pack();
		RefineryUtilities.positionFrameRandomly(chart);
		chart.setVisible(true);

		GraphPlotter.graphIMU(imuMap);

	}

	public static Object[] plotVelAcc(ArrayList<double[]> ecefList, ArrayList<Long> time,
			TreeMap<Long, double[]> estVelMap) throws Exception {
		int n = ecefList.size();
		if (n != time.size()) {
			throw new Exception("FATAL ERROR while analyzing");
		}
		ArrayList<double[]> velList = new ArrayList<double[]>();
		ArrayList<double[]> accList = new ArrayList<double[]>();
		TreeMap<Long, double[]> velMap = new TreeMap<Long, double[]>();
		TreeMap<Long, double[]> errVelMap = new TreeMap<Long, double[]>();
		TreeMap<Long, double[]> accMap = new TreeMap<Long, double[]>();
		for (int i = 1; i < n; i++) {
			double[] ecef1 = ecefList.get(i - 1);
			double[] ecef2 = ecefList.get(i);
			double[] vel = IntStream.range(0, 3).mapToDouble(j -> ecef2[j] - ecef1[j]).toArray();

			velList.add(vel);
			if (i > 1) {
				double[] vel0 = velList.get(i - 2);
				double[] acc = IntStream.range(0, 3).mapToDouble(j -> vel[j] - vel0[j]).toArray();

				accList.add(acc);

			}

		}
		for (int i = 1; i < n - 1; i++) {
			long t = (long) (time.get(i) * 1e-3);
			final int _i = i;
			double[] vel = IntStream.range(0, 3).mapToDouble(j -> (velList.get(_i)[j] + velList.get(_i - 1)[j]) / 2)
					.toArray();
			velMap.put(t, vel);
			if (i > 1) {
				double[] acc = IntStream.range(0, 3)
						.mapToDouble(j -> (accList.get(_i - 1)[j] + accList.get(_i - 2)[j]) / 2).toArray();
				accMap.put(t, acc);
			}
			if (estVelMap.containsKey(t)) {
				double[] ecef = ecefList.get(i);
				double[] err = IntStream.range(0, 3).mapToDouble(j -> estVelMap.get(t)[j] - vel[j]).toArray();
				err = LatLonUtil.ecef2enu(err, ecef, false);
				errVelMap.put(t, err);

			}

		}

		GraphPlotter chart = new GraphPlotter(velMap, "True Vel (in m/s)");
		chart.pack();
		RefineryUtilities.positionFrameRandomly(chart);
		chart.setVisible(true);
		chart = new GraphPlotter(estVelMap, "Est Vel (in m/s)");
		chart.pack();
		RefineryUtilities.positionFrameRandomly(chart);
		chart.setVisible(true);
		chart = new GraphPlotter(errVelMap, "Err Vel ENU frame (in m/s)");
		chart.pack();
		RefineryUtilities.positionFrameRandomly(chart);
		chart.setVisible(true);
		chart = new GraphPlotter(accMap, "True Acc (in m/s)");
		chart.pack();
		RefineryUtilities.positionFrameRandomly(chart);
		chart.setVisible(true);

		return new Object[] { velMap, accMap };
	}
}
