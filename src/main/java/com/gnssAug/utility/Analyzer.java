package com.gnssAug.utility;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.TreeMap;
import java.util.stream.IntStream;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.jfree.ui.RefineryUtilities;

import com.gnssAug.Android.constants.AndroidSensor;
import com.gnssAug.Android.constants.Measurement;
import com.gnssAug.Android.models.IMUsensor;
import com.gnssAug.Android.models.Satellite;
import com.gnssAug.Rinex.models.SatResidual;

public class Analyzer {
	private final static double SpeedofLight = 299792458;

	public static void processAndroid(TreeMap<Long, ArrayList<Satellite>> satMap,
			TreeMap<Long, HashMap<AndroidSensor, IMUsensor>> imuMap, ArrayList<double[]> truePosEcef,
			TreeMap<Long, double[]> trueVelEcef, TreeMap<String, ArrayList<double[]>> estPosMap,
			TreeMap<String, ArrayList<double[]>> estVelMap,
			HashMap<Measurement, HashMap<String, HashMap<String, ArrayList<SatResidual>>>> satResMap,
			boolean outlierAnalyze) throws Exception {

		HashMap<String, TreeMap<Integer, Double>> dopplerMap = new HashMap<String, TreeMap<Integer, Double>>();
		HashMap<String, TreeMap<Integer, Double>> rangeMap = new HashMap<String, TreeMap<Integer, Double>>();

		if (truePosEcef.size() != satMap.size()) {
			throw new Exception("Error in Analyzer processing");
		}
		String estType = estPosMap.firstKey();
		if (estType.equals("EKF - Doppler")) {
			estPosMap.get(estType).remove(0);
			satMap.remove(satMap.firstKey());
			truePosEcef.remove(0);
		}
		long time0 = satMap.firstKey();
		int i = 0;
		for (Long time : satMap.keySet()) {
			double[] truePos = truePosEcef.get(i);
			i++;
			// true velocity list does not contain value for first and last epoch
			if (!trueVelEcef.containsKey(time)) {
				continue;
			}
			ArrayList<Satellite> satList = satMap.get(time);
			String[] obsvCodeList = findObsvCodeSetAndroid(satList);
			int m = obsvCodeList.length;
			double[] rxClkOff = new double[m];
			double[] rxClkDrift = new double[m];

			double[] trueVel = trueVelEcef.get(time);
			int timeDiff = (int) Math.round((time - time0) / 1e3);
			for (int j = 0; j < m; j++) {
				double[] arr = estPosMap.get(estType).get(i);
				if ((j + 3) < arr.length) {
					rxClkOff[j] = arr[j + 3];
//					rangeMap.computeIfAbsent("RxClkOff(offset of 100 added) " + obsvCodeList[j],
//							k -> new TreeMap<Integer, Double>()).put(timeDiff, 100 + rxClkOff[j]);
				}
				if (!estType.equals("EKF - Doppler")) {
					arr = estVelMap.get(estType).get(i);
					if ((j + 3) < arr.length) {
						rxClkDrift[j] = arr[j + 3];
//						dopplerMap.computeIfAbsent("RxClkDrift(offset of 10 added) " + obsvCodeList[j],
//								k -> new TreeMap<Integer, Double>()).put(timeDiff, 10 + rxClkDrift[j]);
					}
				}

			}

			for (Satellite sat : satList) {
				String obsvCode = sat.getObsvCode();
				double[] satPos = sat.getSatEci();
				double[] satVel = sat.getSatVel();
				double trueRange = MathUtil.getEuclidean(truePos, satPos);
				double[] unitLos = IntStream.range(0, 3).mapToDouble(j -> (satPos[j] - truePos[j]) / trueRange)
						.toArray();
				double[] relVel = IntStream.range(0, 3).mapToDouble(j -> satVel[j] - trueVel[j]).toArray();
				double trueRangeRate = IntStream.range(0, 3).mapToDouble(j -> unitLos[j] * relVel[j]).sum();
				double range = sat.getPseudorange();
				double rangeRate = sat.getRangeRate();
				for (int k = 0; k < m; k++) {
					if (obsvCode.equals(obsvCodeList[k])) {
						range -= rxClkOff[k];
						rangeRate -= rxClkDrift[k];
					}
				}
				int svid = sat.getSvid();
				String code = sat.getObsvCode().charAt(0) + "";
				rangeMap.computeIfAbsent(code + svid, k -> new TreeMap<Integer, Double>()).put(timeDiff,
						(range - trueRange));
				dopplerMap.computeIfAbsent(code + svid, k -> new TreeMap<Integer, Double>()).put(timeDiff,
						(rangeRate - trueRangeRate));
			}

		}

		GraphPlotter chart = new GraphPlotter("Error in Range-Rate(in m/s)", dopplerMap);
		chart.pack();
		RefineryUtilities.positionFrameRandomly(chart);
		chart.setVisible(true);

		chart = new GraphPlotter("Error in Range(in m)", rangeMap);
		chart.pack();
		RefineryUtilities.positionFrameRandomly(chart);
		chart.setVisible(true);

		if (outlierAnalyze) {
			// Creating a temp doppler sat res because true velocity list does not contain
			// value for first and last epoch
			/*
			 * HashMap<Measurement, HashMap<String, HashMap<String,
			 * ArrayList<SatResidual>>>> tempSatResMap = new HashMap<Measurement,
			 * HashMap<String, HashMap<String, ArrayList<SatResidual>>>>(); long end =
			 * (long) ((SatMap.lastKey()-time0)/1e3);
			 * 
			 * for(Measurement type:satResMap.keySet()) { tempSatResMap.put(type, new
			 * HashMap<String, HashMap<String, ArrayList<SatResidual>>>()); for(String
			 * est_type:satResMap.get(type).keySet()) {
			 * tempSatResMap.get(type).put(est_type, new HashMap<String,
			 * ArrayList<SatResidual>>()); for(String
			 * svid:satResMap.get(type).get(est_type).keySet()) { ArrayList<SatResidual>
			 * data = satResMap.get(type).get(est_type).get(svid); ArrayList<SatResidual>
			 * new_data = new ArrayList<SatResidual>(); for(int j=0;j<data.size();j++) {
			 * if(data.get(j).getT()>0&&data.get(j).getT()<end) { new_data.add(data.get(j));
			 * } } tempSatResMap.get(type).get(est_type).put(svid, new_data);
			 * 
			 * } } }
			 */
			chart = new GraphPlotter("Outlier in Doppler, based on DIA method(in m/s)", dopplerMap,
					satResMap.get(Measurement.Doppler).get(estType));
			chart.pack();
			RefineryUtilities.positionFrameRandomly(chart);
			chart.setVisible(true);
			chart = new GraphPlotter("Outlier in Range, based on DIA method(in m)", rangeMap,
					satResMap.get(Measurement.Pseudorange).get(estType));
			chart.pack();
			RefineryUtilities.positionFrameRandomly(chart);
			chart.setVisible(true);
		}

		// GraphPlotter.graphIMU(imuMap);

	}

	public static void processIGS(TreeMap<Long, ArrayList<com.gnssAug.Rinex.models.Satellite>> satMap, double[] rxARP,
			HashMap<String, double[]> rxPCO, HashMap<String, ArrayList<double[]>> estPosMap,
			HashMap<String, ArrayList<double[]>> estVelMap,
			HashMap<Measurement, HashMap<String, HashMap<String, ArrayList<SatResidual>>>> satResMap,
			boolean outlierAnalyze) throws Exception {
//		final double pseudorange_priorSigmaOfUnitW = Math.sqrt(0.182);
//		final double doppler_priorSimgaOfUnitW = Math.sqrt(2.397e-4);
		HashMap<String, TreeMap<Integer, Double>> dopplerMap = new HashMap<String, TreeMap<Integer, Double>>();
		HashMap<String, TreeMap<Integer, Double>> rangeMap = new HashMap<String, TreeMap<Integer, Double>>();
		long time0 = satMap.firstKey();
		double alpha = 0.01;
		String estType = "LS";
		int i = 0;
		if (estType.equals("EKF")) {
			satMap.remove(satMap.firstKey());
		}
		for (Long time : satMap.keySet()) {

			ArrayList<com.gnssAug.Rinex.models.Satellite> satList = satMap.get(time);
			String[] obsvCodeList = findObsvCodeSet(satList);
			int m = obsvCodeList.length;
			int n = satList.size();
			double[] rxClkOff = new double[m];
			double[] rxClkDrift = new double[m];
			double[] trueVel = new double[3];
			int timeDiff = (int) ((time - time0) / 1e3);
			for (int j = 0; j < m; j++) {
				rxClkOff[j] = estPosMap.get(estType).get(i)[j + 3];
				rangeMap.computeIfAbsent("RxClkOff(offset of 10 added) " + obsvCodeList[j],
						k -> new TreeMap<Integer, Double>()).put(timeDiff, 10 + rxClkOff[j]);
				if (estType.equals("WLS") || estType.equals("LS")) {
					rxClkDrift[j] = estVelMap.get(estType).get(i)[j + 3];
					dopplerMap.computeIfAbsent("RxClkDrift(offset of 10 added) " + obsvCodeList[j],
							k -> new TreeMap<Integer, Double>()).put(timeDiff, 10 + rxClkDrift[j]);
				}
			}

			for (int j = 0; j < n; j++) {
				com.gnssAug.Rinex.models.Satellite sat = satList.get(j);
				String obsvCode = sat.getObsvCode();
				double[] pco = rxPCO.get(obsvCode);
				double[] rxAPC = IntStream.range(0, 3).mapToDouble(x -> rxARP[x] + pco[x]).toArray();
				double[] satPos = sat.getSatEci();
				double trueRange = MathUtil.getEuclidean(rxAPC, satPos);
				double range = sat.getPseudorange();
				for (int k = 0; k < m; k++) {
					if (obsvCode.equals(obsvCodeList[k])) {
						range -= rxClkOff[k];
					}
				}
				int svid = sat.getSVID();
				String code = sat.getObsvCode().charAt(0) + "";

				rangeMap.computeIfAbsent(code + svid, k -> new TreeMap<Integer, Double>()).put(timeDiff,
						(range - trueRange));
				if (estType.equals("WLS") || estType.equals("LS")) {
					double[] satVel = sat.getSatVel();
					double[] relVel = IntStream.range(0, 3).mapToDouble(k -> satVel[k] - trueVel[k]).toArray();
					double[] unitLos = IntStream.range(0, 3).mapToDouble(k -> (satPos[k] - rxAPC[k]) / trueRange)
							.toArray();
					double trueRangeRate = IntStream.range(0, 3).mapToDouble(k -> unitLos[k] * relVel[k]).sum();
					double rangeRate = sat.getPseudoRangeRate();
					for (int k = 0; k < m; k++) {
						if (obsvCode.equals(obsvCodeList[k])) {
							rangeRate -= rxClkDrift[k];
						}
					}

					dopplerMap.computeIfAbsent(code + svid, k -> new TreeMap<Integer, Double>()).put(timeDiff,
							(rangeRate - trueRangeRate));
				}

			}
			i++;
		}

		// findOutliers(satMap, rangeMap, alpha);

		GraphPlotter chart = new GraphPlotter("Error in Range(in m)", rangeMap);
		chart.pack();
		RefineryUtilities.positionFrameRandomly(chart);
		chart.setVisible(true);

//		chart = new GraphPlotter("Range Outlier and Inliers(in m)", rangeMap, alpha);
//		chart.pack();
//		RefineryUtilities.positionFrameRandomly(chart);
//		chart.setVisible(true);
		if (estType.equals("WLS") || estType.equals("LS")) {
			chart = new GraphPlotter("Error in Range-Rate(in m/s)", dopplerMap);
			chart.pack();
			RefineryUtilities.positionFrameRandomly(chart);
			chart.setVisible(true);

			if (outlierAnalyze) {
				chart = new GraphPlotter("Outlier in Doppler, based on DIA method(in m/s)", dopplerMap,
						satResMap.get(Measurement.Doppler).get(estType));
				chart.pack();
				RefineryUtilities.positionFrameRandomly(chart);
				chart.setVisible(true);
			}

		}

//			chart = new GraphPlotter("Outlier and Inliers Range-Rate(in m/s)", dopplerMap, alpha);
//			chart.pack();
//			RefineryUtilities.positionFrameRandomly(chart);
//			chart.setVisible(true);
		if (outlierAnalyze) {
			chart = new GraphPlotter("Outlier in Range, based on DIA method(in m)", rangeMap,
					satResMap.get(Measurement.Pseudorange).get(estType));
			chart.pack();
			RefineryUtilities.positionFrameRandomly(chart);
			chart.setVisible(true);
		}

		// GraphPlotter.graphIMU(imuMap);

	}

	public static TreeMap<Long, double[]> getVel(ArrayList<double[]> ecefList, ArrayList<Long> time) throws Exception {
		int n = ecefList.size();
		if (n != time.size()) {
			throw new Exception("FATAL ERROR while analyzing");
		}
		ArrayList<double[]> velList = new ArrayList<double[]>();
		ArrayList<double[]> accList = new ArrayList<double[]>();
		TreeMap<Long, double[]> velMap = new TreeMap<Long, double[]>();

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
			long t = time.get(i);
			final int _i = i;
			double[] vel = IntStream.range(0, 3).mapToDouble(j -> (velList.get(_i)[j] + velList.get(_i - 1)[j]) / 2)
					.toArray();
			velMap.put(t, vel);
			if (i > 1) {
				double[] acc = IntStream.range(0, 3)
						.mapToDouble(j -> (accList.get(_i - 1)[j] + accList.get(_i - 2)[j]) / 2).toArray();
				accMap.put(t, acc);
			}

		}

//		GraphPlotter chart = new GraphPlotter(velMap, "True Vel (in m/s)");
//		chart.pack();
//		RefineryUtilities.positionFrameRandomly(chart);
//		chart.setVisible(true);
//		chart = new GraphPlotter(estVelMap, "Est Vel (in m/s)");
//		chart.pack();
//		RefineryUtilities.positionFrameRandomly(chart);
//		chart.setVisible(true);
//		chart = new GraphPlotter(errVelMap, "Err Vel ENU frame (in m/s)");
//		chart.pack();
//		RefineryUtilities.positionFrameRandomly(chart);
//		chart.setVisible(true);
//		chart = new GraphPlotter(accMap, "True Acc (in m/s)");
//		chart.pack();
//		RefineryUtilities.positionFrameRandomly(chart);
//		chart.setVisible(true);

		return velMap;
	}

	private static String[] findObsvCodeSet(ArrayList<com.gnssAug.Rinex.models.Satellite> satList) {
		LinkedHashSet<String> obsvCodeSet = new LinkedHashSet<String>();
		for (int i = 0; i < satList.size(); i++) {
			obsvCodeSet.add(satList.get(i).getObsvCode());
		}
		return obsvCodeSet.toArray(new String[0]);

	}

	private static String[] findObsvCodeSetAndroid(ArrayList<Satellite> satList) {
		LinkedHashSet<String> obsvCodeSet = new LinkedHashSet<String>();
		for (int i = 0; i < satList.size(); i++) {
			obsvCodeSet.add(satList.get(i).getObsvCode());
		}
		return obsvCodeSet.toArray(new String[0]);

	}

	private static void findOutliers(TreeMap<Long, ArrayList<com.gnssAug.Rinex.models.Satellite>> satMap,
			HashMap<String, TreeMap<Integer, Double>> rangeMap, double alpha) {
		long time0 = satMap.firstKey();
		NormalDistribution normal = new NormalDistribution();
		for (Long time : satMap.keySet()) {

			ArrayList<com.gnssAug.Rinex.models.Satellite> satList = satMap.get(time);
			int timeDiff = (int) ((time - time0) / 1e3);
			for (com.gnssAug.Rinex.models.Satellite sat : satList) {

				String svid = sat.getSSI() + "" + sat.getSVID();
				if (rangeMap.get(svid).containsKey(timeDiff)) {
					double err = rangeMap.get(svid).get(timeDiff);
					double pval = 1 - normal.cumulativeProbability(err);
					if (pval < alpha) {
						sat.setTrueOutlier(true);
					}
				}
			}
		}

	}
}
