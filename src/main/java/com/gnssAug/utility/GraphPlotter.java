package com.gnssAug.utility;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.TreeMap;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.ejml.simple.SimpleMatrix;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.ui.ApplicationFrame;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RefineryUtilities;

import com.gnssAug.Android.constants.AndroidSensor;
import com.gnssAug.Android.constants.Measurement;
import com.gnssAug.Android.models.IMUsensor;
import com.gnssAug.Android.models.Satellite;
import com.gnssAug.Rinex.models.SatResidual;

public class GraphPlotter extends ApplicationFrame {

	public GraphPlotter(HashMap<String, ArrayList<SatResidual>> satResMap, String name, boolean isSatRes, boolean flag,
			boolean outlierAnalyze) {

		super(name);
		name = name.split(":")[1];
		ArrayList<JFreeChart> charts = new ArrayList<JFreeChart>();
		if (isSatRes) {

			if (flag) {
				charts.add(ChartFactory.createScatterPlot(name, "GPS-time", name + "(in m or m/s)",
						createDatasetSatRes(satResMap, isSatRes, flag, outlierAnalyze)));

			} else {
				charts.add(ChartFactory.createScatterPlot(name + " vs Elevation Angle", "Elevation-Angle(in degrees)",
						name + "(in m or m/s)", createDatasetSatRes(satResMap, isSatRes, flag, outlierAnalyze)));

			}
		} else {

			if (flag) {
				charts.add(ChartFactory.createScatterPlot("Satellite-Measurement Noise Std Dev", "GPS-time",
						"Noise Std-Dev(in m or m/s)", createDatasetSatRes(satResMap, isSatRes, flag, outlierAnalyze)));
			} else {
				charts.add(ChartFactory.createScatterPlot("Satellite-Measurement Noise Std Dev vs Elevation Angle",
						"Elevation-Angle(in degrees)", "Noise Std-Dev(in m or m/s)",
						createDatasetSatRes(satResMap, isSatRes, flag, outlierAnalyze)));
			}
		}
		for (JFreeChart chart : charts) {
			final ChartPanel chartPanel = new ChartPanel(chart);
			chartPanel.setPreferredSize(new java.awt.Dimension(560, 370));
			chartPanel.setMouseZoomable(true, false);
			// ChartUtils.saveChartAsJPEG(new File(path + chartTitle + ".jpeg"), chart,
			// 1000, 600);
			setContentPane(chartPanel);
		}

	}

	public GraphPlotter(TreeMap<String, ArrayList<double[]>> trajectoryMap, int n, String name, boolean isHorizontal) {
		super(name + " Trajectory");
		JFreeChart chart;

		if (isHorizontal) {
			chart = ChartFactory.createScatterPlot("Horizontal " + name, "East", "North",
					createDataset2DTraj(trajectoryMap, n));
		} else {
			chart = ChartFactory.createXYLineChart("Vertical", "Time(in sec)", "Up",
					createDatasetVerticalTraj(trajectoryMap, n));
		}

		final ChartPanel chartPanel = new ChartPanel(chart);
		chartPanel.setPreferredSize(new java.awt.Dimension(560, 370));
		chartPanel.setMouseZoomable(true, false);
		// ChartUtils.saveChartAsJPEG(new File(path + chartTitle + ".jpeg"), chart,
		// 1000, 600);
		setContentPane(chartPanel);

	}

	public GraphPlotter(String title, String type, TreeMap<Long, Satellite> satMap) throws Exception {
		super(title + " " + type);
		// TODO Auto-generated constructor stub
		final JFreeChart chart;
		if (type.equals("Observables")) {
			chart = ChartFactory.createXYLineChart(title + " " + type, "GPS-time", type,
					createDatasetObservable(satMap));
		} else if(type.equals("Delta-Range")){
			chart = ChartFactory.createXYLineChart(title + " " + type, "GPS-time", type,
					createDatasetDeltaRange(satMap));
		}
		else
		{
			chart = ChartFactory.createXYLineChart(title + " " + type, "GPS-time", type,
					createDatasetSatellitePlot(satMap));
		}
		final ChartPanel chartPanel = new ChartPanel(chart);
		chartPanel.setPreferredSize(new java.awt.Dimension(560, 370));
		chartPanel.setMouseZoomable(true, false);
		// ChartUtils.saveChartAsJPEG(new File(path + chartTitle + ".jpeg"), chart,
		// 1000, 600);
		setContentPane(chartPanel);

	}

	public GraphPlotter(String title, ArrayList<Long> dataList, ArrayList<Long> timeList,int freq) throws IOException {
		super(title + " Satellite Count");
		// TODO Auto-generated constructor stub

		final JFreeChart chart = ChartFactory.createXYLineChart(title + " Satellite Count", "GPS-time",
				title + " Satellite Count", createDatasetSatCount(dataList, timeList,freq));
		final ChartPanel chartPanel = new ChartPanel(chart);
		chartPanel.setPreferredSize(new java.awt.Dimension(560, 370));
		chartPanel.setMouseZoomable(true, false);
		// ChartUtils.saveChartAsJPEG(new File(path + chartTitle + ".jpeg"), chart,
		// 1000, 600);
		setContentPane(chartPanel);

	}

	public GraphPlotter(ArrayList<double[]> dataList,String name) throws IOException {
		super(name+"Z measurement redundancy");
		// TODO Auto-generated constructor stub

		final JFreeChart chart = ChartFactory.createXYLineChart("Z measurement redundancy", "GPS-time", "Count",
				createDatasetRedundancy(dataList));
		final ChartPanel chartPanel = new ChartPanel(chart);
		chartPanel.setPreferredSize(new java.awt.Dimension(560, 370));
		chartPanel.setMouseZoomable(true, false);
		// ChartUtils.saveChartAsJPEG(new File(path + chartTitle + ".jpeg"), chart,
		// 1000, 600);
		setContentPane(chartPanel);

	}

	public GraphPlotter(HashMap<String, HashMap<String, ArrayList<Double>>> dataMap, ArrayList<Long> xList,
			boolean flag) throws IOException {
		super("DOP");
		// TODO Auto-generated constructor stub
		JFreeChart chart = null;
		if (flag) {
			chart = ChartFactory.createXYLineChart("DOP", "GPS-time", "DOP", createDatasetDOP(dataMap, xList));
		} else {
			chart = ChartFactory.createScatterPlot("DOP", "Satellite-Count", "DOP", createDatasetDOP(dataMap, xList));
		}
		final ChartPanel chartPanel = new ChartPanel(chart);
		chartPanel.setPreferredSize(new java.awt.Dimension(560, 370));
		chartPanel.setMouseZoomable(true, false);
		// ChartUtils.saveChartAsJPEG(new File(path + chartTitle + ".jpeg"), chart,
		// 1000, 600);
		setContentPane(chartPanel);

	}

	public GraphPlotter(TreeMap<Long, double[]> ecefMap) throws IOException {
		super("GNSS/INS");
		// TODO Auto-generated constructor stub

		final JFreeChart chart = ChartFactory.createXYLineChart("GNSS/INS", "GPS-time", "GNSS/INS",
				createDatasetConsecutiveINS(ecefMap));
		final ChartPanel chartPanel = new ChartPanel(chart);
		chartPanel.setPreferredSize(new java.awt.Dimension(560, 370));
		chartPanel.setMouseZoomable(true, false);
		// ChartUtils.saveChartAsJPEG(new File(path + chartTitle + ".jpeg"), chart,
		// 1000, 600);
		setContentPane(chartPanel);

	}

	public GraphPlotter(String name, IMUsensor[] imu, long[] time) throws IOException {
		super(name);
		// TODO Auto-generated constructor stub

		final JFreeChart chart = ChartFactory.createXYLineChart(name, "GPS-time", name, createImuDataset(time, imu));
		final ChartPanel chartPanel = new ChartPanel(chart);
		chartPanel.setPreferredSize(new java.awt.Dimension(560, 370));
		chartPanel.setMouseZoomable(true, false);
		// ChartUtils.saveChartAsJPEG(new File(path + chartTitle + ".jpeg"), chart,
		// 1000, 600);
		setContentPane(chartPanel);

	}

	public GraphPlotter(String name, HashMap<String, ArrayList<Double>> data, ArrayList<Long> timeList)
			throws IOException {
		super(name);
		// TODO Auto-generated constructor stub

		final JFreeChart chart = ChartFactory.createXYLineChart(name, "GPS-time", name,
				createDatasetPostUnitW(data, timeList));
		final ChartPanel chartPanel = new ChartPanel(chart);
		chartPanel.setPreferredSize(new java.awt.Dimension(560, 370));
		chartPanel.setMouseZoomable(true, false);
		// ChartUtils.saveChartAsJPEG(new File(path + chartTitle + ".jpeg"), chart,
		// 1000, 600);
		setContentPane(chartPanel);

	}

	public GraphPlotter(String name, HashMap<String, TreeMap<Integer, double[]>> map, boolean isElevAngle)
			throws IOException {

		super(name);
		// TODO Auto-generated constructor stub
		String xAxis = isElevAngle ? "Elevation-Angle(degrees)" : "GPS-Time(seconds)";

		final JFreeChart chart = ChartFactory.createScatterPlot(name, xAxis, name,
				createAnalyseDataset(map, isElevAngle));
		final ChartPanel chartPanel = new ChartPanel(chart);
		chartPanel.setPreferredSize(new java.awt.Dimension(560, 370));
		chartPanel.setMouseZoomable(true, false);
		// ChartUtils.saveChartAsJPEG(new File(path + chartTitle + ".jpeg"), chart,
		// 1000, 600);
		setContentPane(chartPanel);

	}

	public GraphPlotter(String name, HashMap<String, TreeMap<Integer, double[]>> map,
			HashMap<String, ArrayList<SatResidual>> outlierMap) throws Exception {
		super(name);
		// TODO Auto-generated constructor stub

		final JFreeChart chart = ChartFactory.createScatterPlot(name, "GPS-time", name,
				createAnalyseDataset3(map, outlierMap));
		final ChartPanel chartPanel = new ChartPanel(chart);
		chartPanel.setPreferredSize(new java.awt.Dimension(560, 370));
		chartPanel.setMouseZoomable(true, false);
		// ChartUtils.saveChartAsJPEG(new File(path + chartTitle + ".jpeg"), chart,
		// 1000, 600);
		setContentPane(chartPanel);

	}

	public GraphPlotter(String name, TreeMap<Long, double[]> ecefMap) throws IOException {
		super("GNSS/INS  " + name);
		// TODO Auto-generated constructor stub

		final JFreeChart chart = ChartFactory.createXYLineChart(name, "GPS-time", name, createDatasetGnssIns(ecefMap));
		final ChartPanel chartPanel = new ChartPanel(chart);
		chartPanel.setPreferredSize(new java.awt.Dimension(560, 370));
		chartPanel.setMouseZoomable(true, false);
		// ChartUtils.saveChartAsJPEG(new File(path + chartTitle + ".jpeg"), chart,
		// 1000, 600);
		setContentPane(chartPanel);

	}

	public GraphPlotter(String applicationTitle, String chartTitle, HashMap<String, ArrayList<double[]>> dataMap,
			ArrayList<Long> timeList) throws IOException {
		super(applicationTitle);
		// TODO Auto-generated constructor stub

		final JFreeChart chart = ChartFactory.createXYLineChart(applicationTitle, "GPS-time", chartTitle,
				createDataset3dErr(dataMap, timeList));
		final ChartPanel chartPanel = new ChartPanel(chart);
		chartPanel.setPreferredSize(new java.awt.Dimension(560, 370));
		chartPanel.setMouseZoomable(true, false);
		// ChartUtils.saveChartAsJPEG(new File(path + chartTitle + ".jpeg"), chart,
		// 1000, 600);
		setContentPane(chartPanel);

	}

	public GraphPlotter(String applicationTitle, String chartTitle, HashMap<String, double[]> dataMap,
			ArrayList<Long> timeList, boolean flag) throws IOException {
		super(applicationTitle + " Error");
		// TODO Auto-generated constructor stub

		final JFreeChart chart = ChartFactory.createXYLineChart(applicationTitle + " Error", "GPS-time", chartTitle,
				createDatasetENU(dataMap, timeList));
		final ChartPanel chartPanel = new ChartPanel(chart);
		chartPanel.setPreferredSize(new java.awt.Dimension(560, 370));
		chartPanel.setMouseZoomable(true, false);
		// ChartUtils.saveChartAsJPEG(new File(path + chartTitle + ".jpeg"), chart,
		// 1000, 600);
		setContentPane(chartPanel);

	}

	public GraphPlotter(String applicationTitle, String chartTitle, ArrayList<Long> timeList,
			HashMap<String, double[][]> dataMap) throws IOException {
		super(applicationTitle + " Error");
		// TODO Auto-generated constructor stub

		final JFreeChart chart = ChartFactory.createXYLineChart(applicationTitle + " Error", "GPS-time", chartTitle,
				createDatasetENU2(dataMap, timeList));
		final ChartPanel chartPanel = new ChartPanel(chart);
		chartPanel.setPreferredSize(new java.awt.Dimension(560, 370));
		chartPanel.setMouseZoomable(true, false);
		// ChartUtils.saveChartAsJPEG(new File(path + chartTitle + ".jpeg"), chart,
		// 1000, 600);
		setContentPane(chartPanel);

	}

	public GraphPlotter(String applicationTitle, HashMap<String, ArrayList<double[]>> dataMap, String unit)
			throws IOException {
		super(applicationTitle);
		// TODO Auto-generated constructor stub

		final JFreeChart chart = ChartFactory.createScatterPlot(applicationTitle, "East" + unit, "North" + unit,
				createDataset2dErr(dataMap));
		final ChartPanel chartPanel = new ChartPanel(chart);
		chartPanel.setPreferredSize(new java.awt.Dimension(560, 370));
		chartPanel.setMouseZoomable(true, false);
		// ChartUtils.saveChartAsJPEG(new File(path + "2d Error" + ".jpeg"), chart,
		// 1000, 600);
		setContentPane(chartPanel);

	}

	public static void graphPostUnitW(HashMap<Measurement, HashMap<String, ArrayList<Double>>> data,
			ArrayList<Long> timeList) throws IOException {
		for (Measurement key : data.keySet()) {
			String type;
			if (key == Measurement.Pseudorange) {
				type = "PseudoRange";
			} else {
				type = "Doppler";
			}
			HashMap<String, ArrayList<Double>> subData = data.get(key);
			GraphPlotter chart = new GraphPlotter(type + " Posteriori Variance of Unit Weight", subData, timeList);
			chart.pack();
			RefineryUtilities.positionFrameRandomly(chart);
			chart.setVisible(true);
		}
	}

	public static void graphENU(HashMap<String, ArrayList<double[]>> dataMap, ArrayList<Long> timeList, boolean isPos)
			throws IOException {

		String name = null;
		String unit = null;
		if (isPos) {
			name = "GNSS Position";
			unit = "(m)";
		} else {
			name = "GNSS Velocity";
			unit = "(m/s)";
		}
		String[] chartNames = new String[] { "E", "N", "U" };
		for (int i = 0; i < 3; i++) {
			final int index = i;
			HashMap<String, double[]> subDataMap = new HashMap<String, double[]>();
			for (String key : dataMap.keySet()) {
				ArrayList<double[]> data = dataMap.get(key);
				double[] arr = data.stream().mapToDouble(j -> j[index]).toArray();
				subDataMap.put(key, arr);
			}
			GraphPlotter chart = new GraphPlotter(name + "(" + chartNames[i] + ")", chartNames[i] + unit, subDataMap,
					timeList, true);
			chart.pack();
			RefineryUtilities.positionFrameRandomly(chart);
			chart.setVisible(true);
		}
		GraphPlotter chart = new GraphPlotter(name + " 2D-Error", dataMap, unit);
		chart.pack();
		RefineryUtilities.positionFrameRandomly(chart);
		chart.setVisible(true);
		chart = new GraphPlotter(name + " 3d-Error", "3d-Error" + unit, dataMap, timeList);
		chart.pack();
		RefineryUtilities.positionFrameRandomly(chart);
		chart.setVisible(true);

	}

	public static void graphENU(HashMap<String, ArrayList<double[]>> dataMap, ArrayList<Long> timeList, boolean isPos,
			HashMap<String, ArrayList<SimpleMatrix>> Cxx_hat_map) throws Exception {

		String name = null;
		String unit = null;
		if (isPos) {
			name = "GNSS Position";
			unit = "(m)";
		} else {
			name = "GNSS Velocity";
			unit = "(m/s)";
		}
		String[] chartNames = new String[] { "E", "N", "U" };
		for (int i = 0; i < 3; i++) {
			final int index = i;
			HashMap<String, double[][]> subDataMap = new HashMap<String, double[][]>();
			for (String key : dataMap.keySet()) {
				ArrayList<double[]> data = dataMap.get(key);
				ArrayList<SimpleMatrix> Cxx_hat_list = Cxx_hat_map.get(key);
				if (Cxx_hat_list.size() != data.size()) {
					throw new Exception("Size of position estimates and covariance does not match");
				}
				int n = data.size();
				double[][] arr = new double[n][3];
				for (int j = 0; j < n; j++) {
					double stdDev = Math.sqrt(Cxx_hat_list.get(j).get(i, i));
					arr[j][0] = -stdDev;
					arr[j][1] = data.get(j)[i];
					arr[j][2] = stdDev;
				}

				subDataMap.put(key, arr);
			}
			GraphPlotter chart = new GraphPlotter(name + "(" + chartNames[i] + ")", chartNames[i] + unit, timeList,
					subDataMap);
			chart.pack();
			RefineryUtilities.positionFrameRandomly(chart);
			chart.setVisible(true);
		}
		GraphPlotter chart = new GraphPlotter(name + " 2D-Error", dataMap, unit);
		chart.pack();
		RefineryUtilities.positionFrameRandomly(chart);
		chart.setVisible(true);
		chart = new GraphPlotter(name + " 3d-Error", "3d-Error" + unit, dataMap, timeList);
		chart.pack();
		RefineryUtilities.positionFrameRandomly(chart);
		chart.setVisible(true);

	}

	public static void graphGnssIns(TreeMap<Long, double[]> ecefMap, ArrayList<double[]> trueECEFlist,
			ArrayList<Long> timeList) throws IOException {

		GraphPlotter chart = new GraphPlotter(ecefMap);
		chart.pack();
		RefineryUtilities.positionFrameRandomly(chart);
		chart.setVisible(true);

		Iterator<Entry<Long, double[]>> iterator = ecefMap.entrySet().iterator();
		HashMap<String, TreeMap<Long, double[]>> map = new HashMap<String, TreeMap<Long, double[]>>();
		map.put("x-axis", new TreeMap<Long, double[]>());
		map.put("y-axis", new TreeMap<Long, double[]>());
		map.put("z-axis", new TreeMap<Long, double[]>());
		while (iterator.hasNext()) {
			Entry<Long, double[]> entry = iterator.next();
			long time = entry.getKey();
			double[] ecef = entry.getValue();

			map.get("x-axis").computeIfAbsent(time, k -> new double[2])[0] = ecef[0];
			map.get("y-axis").computeIfAbsent(time, k -> new double[2])[0] = ecef[1];
			map.get("z-axis").computeIfAbsent(time, k -> new double[2])[0] = ecef[2];
			int index = timeList.indexOf(time);
			if (index != -1) {
				map.get("x-axis").get(time)[1] = trueECEFlist.get(index)[0];
				map.get("y-axis").get(time)[1] = trueECEFlist.get(index)[1];
				map.get("z-axis").get(time)[1] = trueECEFlist.get(index)[2];
			}
		}
		Iterator<Entry<String, TreeMap<Long, double[]>>> iterator2 = map.entrySet().iterator();
		while (iterator2.hasNext()) {
			Entry<String, TreeMap<Long, double[]>> entry = iterator2.next();
			String name = entry.getKey();
			TreeMap<Long, double[]> data = entry.getValue();
			chart = new GraphPlotter(name, data);
			chart.pack();
			RefineryUtilities.positionFrameRandomly(chart);
			chart.setVisible(true);
		}

	}

	public static void graphSatRes(
			HashMap<Measurement, HashMap<String, HashMap<String, ArrayList<SatResidual>>>> satResMap) {
		graphSatRes(satResMap, false);
	}

	public static void graphSatRes(
			HashMap<Measurement, HashMap<String, HashMap<String, ArrayList<SatResidual>>>> satResMap,
			boolean outlierAnaylze) {
		graphSatRes(satResMap, outlierAnaylze, false);

	}

	public static void graphSatRes(
			HashMap<Measurement, HashMap<String, HashMap<String, ArrayList<SatResidual>>>> satResMap,
			boolean outlierAnaylze, boolean isInnov) {

		for (Measurement key : satResMap.keySet()) {
			HashMap<String, HashMap<String, ArrayList<SatResidual>>> subSatResMap = satResMap.get(key);
			String type;
			if (key == Measurement.Pseudorange) {
				type = "PseudoRange";
			} else {
				type = "Doppler";
			}
			String name = "Satellite-Residual";
			if (isInnov) {
				name = "Satellite-Innovation";
			}
			for (String subKey : subSatResMap.keySet()) {

				// For Satellite Residuals
				GraphPlotter chart = new GraphPlotter(subSatResMap.get(subKey), type + " " + subKey + ": " + name, true,
						true, outlierAnaylze);
				chart.pack();
				RefineryUtilities.positionFrameRandomly(chart);
				chart.setVisible(true);

				chart = new GraphPlotter(subSatResMap.get(subKey), type + " " + subKey + ": " + name, true, false,
						outlierAnaylze);
				chart.pack();
				RefineryUtilities.positionFrameRandomly(chart);
				chart.setVisible(true);

				// For Satellite measurement noise std dev

				chart = new GraphPlotter(subSatResMap.get(subKey),
						type + " " + subKey + ": Satellite-Measurement Noise Std Dev", false, true, outlierAnaylze);
				chart.pack();
				RefineryUtilities.positionFrameRandomly(chart);
				chart.setVisible(true);

				chart = new GraphPlotter(subSatResMap.get(subKey),
						type + " " + subKey + ": Satellite-Measurement Noise Std Dev", false, false, outlierAnaylze);
				chart.pack();
				RefineryUtilities.positionFrameRandomly(chart);
				chart.setVisible(true);

			}
		}

	}

	public static void graphTrueError(String name, HashMap<String, TreeMap<Integer, double[]>> map) throws IOException {

		GraphPlotter chart = new GraphPlotter(name, map, false);
		chart.pack();
		RefineryUtilities.positionFrameRandomly(chart);
		chart.setVisible(true);
		
		chart = new GraphPlotter(name, map, true);
		chart.pack();
		RefineryUtilities.positionFrameRandomly(chart);
		chart.setVisible(true);

	}

	public static void graphIMU(TreeMap<Long, HashMap<AndroidSensor, IMUsensor>> imuMap) throws IOException {
		int n = imuMap.size();
		IMUsensor[] acc = new IMUsensor[n];
		IMUsensor[] gyro = new IMUsensor[n];
		IMUsensor[] mag = new IMUsensor[n];
		long[] time = new long[n];
		int i = 0;
		for (long t : imuMap.keySet()) {
			time[i] = t;
			HashMap<AndroidSensor, IMUsensor> imu = imuMap.get(t);
			acc[i] = imu.get(AndroidSensor.Accelerometer);
			gyro[i] = imu.get(AndroidSensor.Gyroscope);
			mag[i] = imu.get(AndroidSensor.Magnetometer);
			i++;
		}

		GraphPlotter chart = new GraphPlotter("Accelerometer (in m/s2)", acc, time);
		chart.pack();
		RefineryUtilities.positionFrameRandomly(chart);
		chart.setVisible(true);
		chart = new GraphPlotter("Gyroscope (in rad/s)", gyro, time);
		chart.pack();
		RefineryUtilities.positionFrameRandomly(chart);
		chart.setVisible(true);
		chart = new GraphPlotter("Magnetometer (in microtesla)", mag, time);
		chart.pack();
		RefineryUtilities.positionFrameRandomly(chart);
		chart.setVisible(true);

	}

	public static void graphRedundancy(ArrayList<double[]> redundancyList) throws IOException {
		GraphPlotter chart = new GraphPlotter(redundancyList,"");
		chart.pack();
		RefineryUtilities.positionFrameRandomly(chart);
		chart.setVisible(true);
	}
	
	public static void graphRedundancy(ArrayList<double[]> redundancyList,String name) throws IOException {
		GraphPlotter chart = new GraphPlotter(redundancyList,name);
		chart.pack();
		RefineryUtilities.positionFrameRandomly(chart);
		chart.setVisible(true);
	}

	public static void graphTrajectory(TreeMap<String, ArrayList<double[]>> trajectoryPosMap,
			TreeMap<String, ArrayList<double[]>> trajectoryVelMap, int n) throws IOException {

		GraphPlotter chart = new GraphPlotter(trajectoryPosMap, n, "Position(in m)", true);
		chart.pack();
		RefineryUtilities.positionFrameRandomly(chart);
		chart.setVisible(true);
		chart = new GraphPlotter(trajectoryPosMap, n, "Position(in m)", false);
		chart.pack();
		RefineryUtilities.positionFrameRandomly(chart);
		chart.setVisible(true);
		chart = new GraphPlotter(trajectoryVelMap, n, "Velocity(in m/s)", true);
		chart.pack();
		RefineryUtilities.positionFrameRandomly(chart);
		chart.setVisible(true);
		chart = new GraphPlotter(trajectoryVelMap, n, "Velocity(in m/s)", false);
		chart.pack();
		RefineryUtilities.positionFrameRandomly(chart);
		chart.setVisible(true);

	}

	public static void graphSatCount(HashMap<Measurement, TreeMap<String, ArrayList<Long>>> satCountMap,
			ArrayList<Long> timeList,int freq) throws IOException {

		for (Measurement key : satCountMap.keySet()) {
			String type;
			if (key == Measurement.Pseudorange) {
				type = "PseudoRange ";
			} else {
				type = "Doppler ";
			}
			for (String subKey : satCountMap.get(key).keySet()) {
				type += subKey;
				GraphPlotter chart = new GraphPlotter(type, satCountMap.get(key).get(subKey), timeList,freq);
				chart.pack();
				RefineryUtilities.positionFrameRandomly(chart);
				chart.setVisible(true);
			}
		}
	}

	public static void graphDOP(HashMap<String, ArrayList<double[]>> dataMap, ArrayList<Long> satCountList,
			ArrayList<Long> timeList,int freq) throws Exception {

		HashMap<String, HashMap<String, ArrayList<Double>>> dopMap = new HashMap<String, HashMap<String, ArrayList<Double>>>();
		for (String key : dataMap.keySet()) {
			ArrayList<double[]> dopList = dataMap.get(key);
			ArrayList<Double> gdopList = new ArrayList<Double>();
			ArrayList<Double> pdopList = new ArrayList<Double>();
			ArrayList<Double> hdopList = new ArrayList<Double>();
			ArrayList<Double> vdopList = new ArrayList<Double>();
			ArrayList<Double> tdopList = new ArrayList<Double>();
			int n = dopList.size();
			if (n != timeList.size()) {
				throw new Exception("DOP list size does not match timeList size");
			}
			for (int i = 0; i < n; i++) {
				double[] dopDiag = dopList.get(i);
				gdopList.add(Math.sqrt(dopDiag[0] + dopDiag[1] + dopDiag[2] + dopDiag[3]));
				pdopList.add(Math.sqrt(dopDiag[0] + dopDiag[1] + dopDiag[2]));
				hdopList.add(Math.sqrt(dopDiag[0] + dopDiag[1]));
				vdopList.add(Math.sqrt(dopDiag[2]));
				tdopList.add(Math.sqrt(dopDiag[3]));
			}
			dopMap.computeIfAbsent("GDOP", k -> new HashMap<String, ArrayList<Double>>()).put(key, gdopList);
			dopMap.computeIfAbsent("PDOP", k -> new HashMap<String, ArrayList<Double>>()).put(key, pdopList);
			dopMap.computeIfAbsent("HDOP", k -> new HashMap<String, ArrayList<Double>>()).put(key, hdopList);
			dopMap.computeIfAbsent("VDOP", k -> new HashMap<String, ArrayList<Double>>()).put(key, vdopList);
			dopMap.computeIfAbsent("TDOP", k -> new HashMap<String, ArrayList<Double>>()).put(key, tdopList);

		}

		GraphPlotter chart = new GraphPlotter(dopMap, timeList, true);
		chart.pack();
		RefineryUtilities.positionFrameRandomly(chart);
		chart.setVisible(true);
		chart = new GraphPlotter(dopMap, satCountList, false);
		chart.pack();
		RefineryUtilities.positionFrameRandomly(chart);
		chart.setVisible(true);
		chart = new GraphPlotter("", satCountList, timeList,freq);
		chart.pack();
		RefineryUtilities.positionFrameRandomly(chart);
		chart.setVisible(true);

	}

	public static void graphDeltaRange(TreeMap<Long, ArrayList<Satellite>> satMap, ArrayList<double[]> truePosEcef)
			throws Exception {
		if (truePosEcef.size() != satMap.size()) {
			throw new Exception("Error in DeltaRange Plotting");
		}
		HashMap<String, TreeMap<Long, Satellite>> satListMap = new HashMap<String, TreeMap<Long, Satellite>>();
		long time0 = satMap.firstKey();
		int i = 0;
		for (long time : satMap.keySet()) {
			double[] truePos = truePosEcef.get(i);
			long t = (long) ((time - time0) * (1e-3));
			ArrayList<Satellite> satList = satMap.get(time);
			for (Satellite sat : satList) {
				double[] satPos = sat.getSatEci();
				double trueRange = MathUtil.getEuclidean(truePos, satPos);
				sat.setTrueRange(trueRange);
				String id = sat.getObsvCode().charAt(0) + "" + sat.getSvid();
				satListMap.computeIfAbsent(id, k -> new TreeMap<Long, Satellite>()).put(t, sat);
			}
			i++;
		}

		for (String id : satListMap.keySet()) {
			GraphPlotter chart = new GraphPlotter(id, "Delta-Range", satListMap.get(id));
			chart.pack();
			RefineryUtilities.positionFrameRandomly(chart);
			chart.setVisible(true);
			chart = new GraphPlotter(id, "Observables", satListMap.get(id));
			chart.pack();
			RefineryUtilities.positionFrameRandomly(chart);
			chart.setVisible(true);
			chart = new GraphPlotter(id, "Satellite-Plot", satListMap.get(id));
			chart.pack();
			RefineryUtilities.positionFrameRandomly(chart);
			chart.setVisible(true);
		}
	}

	private XYDataset createDataset3dErr(HashMap<String, ArrayList<double[]>> dataMap, ArrayList<Long> timeList) {
		XYSeriesCollection dataset = new XYSeriesCollection();
		for (String key : dataMap.keySet()) {
			final XYSeries series = new XYSeries(key);
			ArrayList<double[]> list = dataMap.get(key);
			for (int i = 0; i < list.size(); i++) {
				double[] data = list.get(i);
				double err = Math.sqrt(Arrays.stream(data).map(j -> j * j).sum());
				series.add(timeList.get(i), Double.valueOf(err));
			}
			dataset.addSeries(series);
		}

		return dataset;

	}

	private XYDataset createDatasetConsecutiveINS(TreeMap<Long, double[]> ecefMap) {
		XYSeriesCollection dataset = new XYSeriesCollection();
		final XYSeries series = new XYSeries("GNSS/INS");
		Iterator<Entry<Long, double[]>> iterator = ecefMap.entrySet().iterator();
		double[] prev = iterator.next().getValue();
		while (iterator.hasNext()) {
			Entry<Long, double[]> entry = iterator.next();
			long time = entry.getKey();
			double[] ecef = entry.getValue();
			double data = 0;
			for (int i = 0; i < 3; i++) {
				data += Math.pow(ecef[i] - prev[i], 2);
			}
			data = Math.sqrt(data);
			series.add(time, data);
			prev = ecef;
		}
		dataset.addSeries(series);

		return dataset;

	}

	private XYDataset createImuDataset(long[] time, IMUsensor[] imu) {
		XYSeriesCollection dataset = new XYSeriesCollection();

		final XYSeries valX = new XYSeries("Uncal Val X");
		final XYSeries valY = new XYSeries("Uncal Val Y");
		final XYSeries valZ = new XYSeries("Uncal Val Z");
		final XYSeries biasX = new XYSeries("Bias value X");
		final XYSeries biasY = new XYSeries("Bias value Y");
		final XYSeries biasZ = new XYSeries("Bias value Z");
		int n = imu.length;
		for (int i = 0; i < n; i++) {
			long t = (long) ((time[i] - time[0]) * 1e-3);
			valX.add(t, imu[i].getVal()[0]);
			valY.add(t, imu[i].getVal()[1]);
			valZ.add(t, imu[i].getVal()[2]);
			biasX.add(t, imu[i].getBias()[0]);
			biasY.add(t, imu[i].getBias()[1]);
			biasZ.add(t, imu[i].getBias()[2]);

		}
		dataset.addSeries(valX);
		dataset.addSeries(valY);
		dataset.addSeries(valZ);
		dataset.addSeries(biasX);
		dataset.addSeries(biasY);
		dataset.addSeries(biasZ);

		return dataset;

	}

	private XYDataset createAnalyseDataset(HashMap<String, TreeMap<Integer, double[]>> map, boolean isElevAngle) {
		XYSeriesCollection dataset = new XYSeriesCollection();
		ArrayList<Double> dataSeries = new ArrayList<Double>();
		for (String key : map.keySet()) {
			String name = key;
			final XYSeries series = new XYSeries(name);
			TreeMap<Integer, double[]> data = map.get(key);
			for (int x : data.keySet()) {
				double y = data.get(x)[0];
				if (isElevAngle) {
					series.add(data.get(x)[1], y);
				} else {
					series.add(x, y);
				}
				dataSeries.add(Math.abs(y));
			}
			dataset.addSeries(series);
		}

		double avg = dataSeries.stream().mapToDouble(i -> Math.abs(i)).average().orElse(0);
		Collections.sort(dataSeries);
		int n = dataSeries.size();
		double q50 = dataSeries.get((int) (n * 0.5));
		double q75 = dataSeries.get((int) (n * 0.75));
		double q95 = dataSeries.get((int) (n * 0.95));
		dataset.addSeries(new XYSeries("Average: " + avg));
		dataset.addSeries(new XYSeries("Q50: " + q50));
		dataset.addSeries(new XYSeries("Q75: " + q75));
		dataset.addSeries(new XYSeries("Q95: " + q95));
		return dataset;
	}

	private XYDataset createAnalyseDataset3(HashMap<String, TreeMap<Integer, double[]>> map,
			HashMap<String, ArrayList<SatResidual>> outlierMap) throws Exception {
		XYSeriesCollection dataset = new XYSeriesCollection();
		final XYSeries outliers = new XYSeries("outliers");
		final XYSeries inliers = new XYSeries("inliers");

		int in = 0;
		int out = 0;
		for (String key : map.keySet()) {

			TreeMap<Integer, double[]> data = map.get(key);
			if (outlierMap.containsKey(key)) {
				ArrayList<SatResidual> outlierList = outlierMap.get(key);

//				if (data.size() != outlierList.size()) {
//					throw new Exception("Error while Analysing data in GraphPlotter");
//				}
				int i = 0;
				for (int x : data.keySet()) {
					while (Math.abs(outlierList.get(i).getT() - x) > 0.5) {
						i++;

					}
					double y = data.get(x)[0];
					SatResidual satRes = outlierList.get(i);
					i++;

					if (satRes.isOutlier()) {
						outliers.add(x, y);
						out++;
					} else {
						inliers.add(x, y);
						in++;
					}

				}
			}
		}
		inliers.setKey("Inliers (Count: " + in + ")");
		outliers.setKey("Outliers (Count: " + out + ")");
		dataset.addSeries(outliers);
		dataset.addSeries(inliers);

		return dataset;
	}

	private XYDataset createDatasetGnssIns(TreeMap<Long, double[]> ecefMap) {
		XYSeriesCollection dataset = new XYSeriesCollection();

		final XYSeries est = new XYSeries("Estimated");
		final XYSeries truth = new XYSeries("True");
		Iterator<Entry<Long, double[]>> iterator = ecefMap.entrySet().iterator();
		while (iterator.hasNext()) {
			Entry<Long, double[]> entry = iterator.next();
			long time = entry.getKey();
			double[] ecef = entry.getValue();

			est.add(time, ecef[0]);
			if (ecef[1] != 0) {
				truth.add(time, ecef[1]);
			}
		}
		dataset.addSeries(est);
		dataset.addSeries(truth);

		return dataset;
	}

	private XYDataset createDatasetENU(HashMap<String, double[]> dataMap, ArrayList<Long> timeList) {
		XYSeriesCollection dataset = new XYSeriesCollection();
		for (String key : dataMap.keySet()) {
			final XYSeries series = new XYSeries(key);
			double[] data = dataMap.get(key);
			for (int i = 0; i < data.length; i++) {
				series.add(timeList.get(i), Double.valueOf(data[i]));
			}
			dataset.addSeries(series);
		}

		return dataset;

	}

	private XYDataset createDatasetENU2(HashMap<String, double[][]> dataMap, ArrayList<Long> timeList) {
		XYSeriesCollection dataset = new XYSeriesCollection();
		for (String key : dataMap.keySet()) {
			final XYSeries series = new XYSeries(key);
			final XYSeries std_minus = new XYSeries(key + ": Std-Dev(-)");
			final XYSeries std_plus = new XYSeries(key + ": Std-Dev(+)");
			double[][] data = dataMap.get(key);
			for (int i = 0; i < data.length; i++) {
				long t = timeList.get(i);
				std_minus.add(t, Double.valueOf(data[i][0]));
				series.add(t, Double.valueOf(data[i][1]));
				std_plus.add(t, Double.valueOf(data[i][2]));
			}
			dataset.addSeries(std_minus);
			dataset.addSeries(std_plus);
			dataset.addSeries(series);
		}

		return dataset;

	}

	private XYDataset createDatasetPostUnitW(HashMap<String, ArrayList<Double>> dataMap, ArrayList<Long> timeList) {
		XYSeriesCollection dataset = new XYSeriesCollection();

		for (String key : dataMap.keySet()) {
			ArrayList<Double> data = dataMap.get(key);
			final XYSeries series = new XYSeries(key + " Posteriori Variance of Unit Weight");
			double sum = 0;
			int count = 0;

			for (int i = 0; i < data.size(); i++) {
				double val = data.get(i);
				if (val == 0 || val == -1) {
					continue;
				}
				sum += val;
				count++;
				series.add(timeList.get(i), (Double) val);
			}
			Collections.sort(data);

			double _mean = ((int) ((sum / count) * 1e4)) / 1e4;
			int q50 = (int) (count * 0.5);
			double _Q50 = ((int) ((data.get(q50) * 1e4))) / 1e4;
			int q75 = (int) (count * 0.75);
			double _Q75 = ((int) ((data.get(q75) * 1e4))) / 1e4;
			// avg = Math.round(avg * 1000) / 1000;
			final XYSeries mean = new XYSeries(key + " Mean Post Var of Unit W: " + _mean);
			final XYSeries Q50 = new XYSeries(key + " Median Post Var of Unit W: " + _Q50);
			final XYSeries Q75 = new XYSeries(key + " Q75 Post Var of Unit W: " + _Q75);
			for (int i = 0; i < data.size(); i++) {
				long time = timeList.get(i);
				mean.add(time, _mean);
				Q50.add(time, _Q50);
				Q75.add(time, _Q75);
			}
			dataset.addSeries(series);
			dataset.addSeries(mean);
			dataset.addSeries(Q50);
			dataset.addSeries(Q75);

		}
		return dataset;

	}

	private XYDataset createDataset2dErr(HashMap<String, ArrayList<double[]>> dataMap) {
		XYSeriesCollection dataset = new XYSeriesCollection();
		for (String key : dataMap.keySet()) {
			final XYSeries series = new XYSeries(key);
			ArrayList<double[]> list = dataMap.get(key);
			for (int i = 0; i < list.size(); i++) {
				double[] data = list.get(i);
				series.add(data[0], data[1]);
			}
			dataset.addSeries(series);
		}

		return dataset;

	}

	private XYDataset createDatasetRedundancy(ArrayList<double[]> dataList) {
		XYSeriesCollection dataset = new XYSeriesCollection();

		final XYSeries diff = new XYSeries("Difference");
		final XYSeries satCount = new XYSeries("SatCount");
		final XYSeries rX = new XYSeries("rX");
		final XYSeries rW = new XYSeries("rW");
		final XYSeries rZ = new XYSeries("rZ");
		// final XYSeries rSum = new XYSeries("rSum");
		for (int i = 0; i < dataList.size(); i++) {
			double[] data = dataList.get(i);
			satCount.add(i, data[0]);
			// rSum.add(i, data[1]);
			diff.add(i, data[0] - data[4]);
			rX.add(i, data[2]);
			rW.add(i, data[3]);
			rZ.add(i, data[4]);
		}
		dataset.addSeries(diff);
		dataset.addSeries(satCount);
		dataset.addSeries(rX);
		dataset.addSeries(rW);
		dataset.addSeries(rZ);
		// dataset.addSeries(rSum);

		return dataset;

	}

	private XYDataset createDatasetDOP(HashMap<String, HashMap<String, ArrayList<Double>>> dataMap,
			ArrayList<Long> timeList) {
		XYSeriesCollection dataset = new XYSeriesCollection();
		for (String key : dataMap.keySet()) {
			for (String subKey : dataMap.get(key).keySet()) {
				final XYSeries series = new XYSeries(key + " " + subKey);
				ArrayList<Double> data = dataMap.get(key).get(subKey);
				for (int i = 0; i < data.size(); i++) {
					series.add(timeList.get(i), data.get(i));
				}
				dataset.addSeries(series);
			}
		}

		return dataset;

	}

	private XYDataset createDatasetSatCount(ArrayList<Long> data, ArrayList<Long> timeList,int freq) {
		XYSeriesCollection dataset = new XYSeriesCollection();
		final XYSeries series = new XYSeries("SatCount");
		double avg = data.stream().mapToDouble(i -> i).average().orElseThrow();
		avg = Math.round(avg * 100) / 100.0;
		long t0 = timeList.get(0);
		for (int i = 0; i < data.size(); i++) {
			long t = timeList.get(i);
			if(t-t0>freq)
			{
				series.add(timeList.get(i-1)+freq, null);
			}
			series.add(timeList.get(i), data.get(i));
			t0 = t;

		}
		final XYSeries avgSeries = new XYSeries("Avg:" + avg);
		dataset.addSeries(series);
		dataset.addSeries(avgSeries);

		return dataset;

	}

	private XYDataset createDatasetDeltaRange(TreeMap<Long, Satellite> satMap) throws Exception {

		XYSeriesCollection dataset = new XYSeriesCollection();
		final XYSeries pr_dr_series = new XYSeries("PR Derived DeltaRange");
		final XYSeries prRate_dr_series = new XYSeries("Doppler Derived");
		final XYSeries phase_dr_series = new XYSeries("Carrier-Phase Derived");
		final XYSeries true_dr_series = new XYSeries("True DeltaRange");
		
		long t_prev = -999;
		double pr_prev = 0;
		double prRate_prev = 0;
		double phase_prev = 0;
		double trueRange_prev = 0;
		for (Long t : satMap.keySet()) {

			Satellite sat = satMap.get(t);
			double trueRange = sat.getTrueRange();
			double pr = sat.getPseudorange();
			double prRate = sat.getPseudorangeRateMetersPerSecond();
			double phase = sat.getAccumulatedDeltaRangeMeters();
			if ((t - t_prev) > 1.1) {
				pr_dr_series.add(t, null);
				prRate_dr_series.add(t, null);
				phase_dr_series.add(t, null);
				true_dr_series.add(t, null);
			} else {
				double pr_dr = pr - pr_prev;
				double prRate_dr = (prRate + prRate_prev) * (t - t_prev) / 2;
				double phase_dr = phase - phase_prev;
				double true_dr = trueRange - trueRange_prev;
				pr_dr_series.add(t, (Double) pr_dr);
				prRate_dr_series.add(t, (Double) prRate_dr);
				phase_dr_series.add(t, (Double) phase_dr);
				true_dr_series.add(t, (Double) true_dr);

			}
			t_prev = t;
			pr_prev = pr;
			prRate_prev = prRate;
			phase_prev = phase;
			trueRange_prev = trueRange;
		
		}
		dataset.addSeries(true_dr_series);
		dataset.addSeries(pr_dr_series);
		dataset.addSeries(prRate_dr_series);
		//dataset.addSeries(phase_dr_series);

		return dataset;

	}

	private XYDataset createDatasetObservable(TreeMap<Long, Satellite> satMap) {
		XYSeriesCollection dataset = new XYSeriesCollection();
		final XYSeries pr_series = new XYSeries("Pseudorange");
		final XYSeries prRate_series = new XYSeries("Pseudorange Rate(Doppler)");
		double t0 = satMap.firstEntry().getKey();
		double tr0 = satMap.firstEntry().getValue().getTrueRange();
		double pr0 = satMap.firstEntry().getValue().getPseudorange();
		double prRate0 = satMap.firstEntry().getValue().getPseudorangeRateMetersPerSecond();
		for (Long t : satMap.keySet()) {
			Satellite sat = satMap.get(t);
			double pr = sat.getPseudorange();
			double prRate = sat.getPseudorangeRateMetersPerSecond();
			double tr = sat.getTrueRange();
			if ((t - t0) > 1.1) {
				pr_series.add(t0 + 1, null);
				prRate_series.add(t0 + 1, null);
			}
			pr_series.add(t, (Double) (pr));
			prRate_series.add(t, (Double) (prRate));
			t0 = t;
		}
		dataset.addSeries(pr_series);
		dataset.addSeries(prRate_series);
		return dataset;
	}
	
	private XYDataset createDatasetSatellitePlot(TreeMap<Long, Satellite> satMap) {
		XYSeriesCollection dataset = new XYSeriesCollection();
		final XYSeries x = new XYSeries("X");
		final XYSeries y = new XYSeries("Y");
		final XYSeries z = new XYSeries("Z");
		Long t0 = satMap.firstEntry().getKey();
		Satellite sat0 = satMap.get(t0);
		double[] ecef0 = sat0.getSatEci();
		
		for (Long t : satMap.keySet()) {
			Satellite sat = satMap.get(t);
			double[] ecef = sat.getSatEci();
		
			if ((t - t0) > 1.1) {
				x.add(t0 + 1, null);
				y.add(t0 + 1, null);
				z.add(t0 + 1, null);
			}
			x.add(t0 + 1, ecef[0]-ecef0[0]);
			y.add(t0 + 1, ecef[1]-ecef0[1]);
			z.add(t0 + 1, ecef[2]-ecef0[2]);
			t0 = t;
		}
		dataset.addSeries(x);
		dataset.addSeries(y);
		dataset.addSeries(z);
	
		return dataset;
	}

	private XYDataset createDatasetSatRes(HashMap<String, ArrayList<SatResidual>> dataMap, boolean isSatRes,
			boolean flag, boolean outlierAnalyze) {
		XYSeriesCollection dataset = new XYSeriesCollection();

		if (outlierAnalyze) {
			final XYSeries inlier = new XYSeries("inlier");
			final XYSeries outlier = new XYSeries("outlier");
			for (String key : dataMap.keySet()) {
				ArrayList<SatResidual> dataList = dataMap.get(key);

				double t0 = 0;
				for (int i = 0; i < dataList.size(); i++) {
					SatResidual satData = dataList.get(i);
					double x = 0;
					if (flag) {
						double t = satData.getT();

						x = t;
					} else {
						double elevAngle = Math.toDegrees(satData.getElevAngle());
						x = elevAngle;
					}
					double data;
					if (isSatRes) {
						data = satData.getResidual();
					} else {
						data = satData.getNoiseStdDev();
					}
					if (satData.isOutlier()) {
						outlier.add(x, data);
					} else {
						inlier.add(x, data);
					}

				}

			}
			dataset.addSeries(outlier);
			dataset.addSeries(inlier);
		} else {

			ArrayList<Double> dataSeries = new ArrayList<Double>();
			for (String key : dataMap.keySet()) {
				final XYSeries series = new XYSeries(key);
				ArrayList<SatResidual> dataList = dataMap.get(key);

				double t0 = 0;
				for (int i = 0; i < dataList.size(); i++) {
					SatResidual satData = dataList.get(i);
					double x = 0;
					if (flag) {
						double t = satData.getT();
						if (t - t0 > 1) {
							series.add(t0, null);
						}
						t0 = t;
						x = t;
					} else {
						double elevAngle = Math.toDegrees(satData.getElevAngle());
						x = elevAngle;
					}

					double data;
					if (isSatRes) {
						data = satData.getResidual();
					} else {
						data = satData.getNoiseStdDev();
					}
					series.add(x, data);
					dataSeries.add(data);
				}

				dataset.addSeries(series);
			}
			dataSeries.replaceAll(i->Math.abs(i));
			double avg = dataSeries.stream().mapToDouble(i -> i).average().orElse(0);
			Collections.sort(dataSeries);
			int n = dataSeries.size();
			double q50 = dataSeries.get((int) (n * 0.5));
			double q75 = dataSeries.get((int) (n * 0.75));
			double q95 = dataSeries.get((int) (n * 0.95));
			dataset.addSeries(new XYSeries("Average: " + avg));
			dataset.addSeries(new XYSeries("Q50: " + q50));
			dataset.addSeries(new XYSeries("Q75: " + q75));
			dataset.addSeries(new XYSeries("Q95: " + q95));
		}

		return dataset;
	}

	private static XYDataset createDataset2DTraj(TreeMap<String, ArrayList<double[]>> trajectoryMap, int n) {
		XYSeriesCollection dataset = new XYSeriesCollection();
		int j = 0;
		while (trajectoryMap.get("True").get(j)[0] == -999) {
			j++;
		}
		double[] start = trajectoryMap.get("True").get(j);
		j = n - 1;
		while (trajectoryMap.get("True").get(j)[0] == -999) {
			j--;
		}
		double[] end = trajectoryMap.get("True").get(j);
		final XYSeries s = new XYSeries("START");
		final XYSeries e = new XYSeries("END");
		s.add(start[0], start[1]);
		e.add(end[0], end[1]);
		dataset.addSeries(s);
		dataset.addSeries(e);
		for (String key : trajectoryMap.keySet()) {
			final XYSeries series = new XYSeries(key);
			ArrayList<double[]> trajectory = trajectoryMap.get(key);
			for (int i = 0; i < n; i++) {
				double[] data = trajectory.get(i);
				if ((data[0] != -999) && (data[1] != -999) && (data[2] != -999)) {
					series.add(data[0], data[1]);
				}
			}

			dataset.addSeries(series);
		}

		return dataset;

	}

	private static XYDataset createDatasetVerticalTraj(TreeMap<String, ArrayList<double[]>> trajectoryMap, int n) {
		XYSeriesCollection dataset = new XYSeriesCollection();

		for (String key : trajectoryMap.keySet()) {
			final XYSeries series = new XYSeries(key);
			ArrayList<double[]> trajectory = trajectoryMap.get(key);
			for (int i = 0; i < n; i++) {
				double[] data = trajectory.get(i);
				if ((data[0] != -999) && (data[1] != -999) && (data[2] != -999)) {
					series.add(i, data[2]);
				} else {
					series.add(i, null);
				}
			}

			dataset.addSeries(series);
		}
		return dataset;

	}
}