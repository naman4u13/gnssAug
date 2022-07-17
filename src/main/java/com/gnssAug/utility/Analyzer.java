package com.gnssAug.utility;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeMap;

import org.jfree.ui.RefineryUtilities;

import com.gnssAug.Android.models.Satellite;

public class Analyzer {

	public static void process(TreeMap<Long, ArrayList<Satellite>> SatMap, HashMap<Long, double[]> err)
			throws IOException {
		ArrayList<ArrayList<Satellite>> SVlist = new ArrayList<ArrayList<Satellite>>(SatMap.values());
		HashMap<String, TreeMap<Integer, Double>> obsvMap = new HashMap<String, TreeMap<Integer, Double>>();
		HashMap<String, Double> firstVal = new HashMap<String, Double>();
		long time0 = SatMap.firstKey();
		for (Long time : SatMap.keySet()) {
			int timeDiff = (int) ((time - time0) / 1e3);
			ArrayList<Satellite> satList = SatMap.get(time);
			for (Satellite sat : satList) {
				double rangeRate = sat.getPseudorangeRateMetersPerSecond();
				int svid = sat.getSvid();
				String code = sat.getObsvCode().charAt(0) + "";
				double first = rangeRate;
				if (firstVal.containsKey(code + svid)) {
					first = firstVal.get(code + svid);
				} else {
					firstVal.put(code + svid, first);
				}
//				obsvMap.computeIfAbsent(code + svid, k -> new TreeMap<Integer, Double>()).put(timeDiff,
//						rangeRate - first);
			}
			if (err.containsKey(time)) {
				obsvMap.computeIfAbsent("speed", k -> new TreeMap<Integer, Double>()).put(timeDiff, err.get(time)[1]);
				obsvMap.computeIfAbsent("position", k -> new TreeMap<Integer, Double>()).put(timeDiff,
						err.get(time)[0]);
			}
		}

		GraphPlotter chart = new GraphPlotter(obsvMap);
		chart.pack();
		RefineryUtilities.positionFrameRandomly(chart);
		chart.setVisible(true);

	}
}
