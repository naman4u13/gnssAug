package com.gnssAug.utility;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.TreeMap;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.ui.ApplicationFrame;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RefineryUtilities;

import com.gnssAug.Android.constants.AndroidSensor;
import com.gnssAug.Android.models.IMUsensor;
import com.gnssAug.Rinex.models.SatResidual;

public class GraphPlotter extends ApplicationFrame {

	public GraphPlotter(HashMap<Integer, ArrayList<SatResidual>> satResMap, String name, boolean flag) {
		super(name + ": Satellite-Residual");
		JFreeChart chart;
		if (flag) {
			chart = ChartFactory.createScatterPlot("Satellite-Residual", "GPS-time", "Satellite-Residual(in m)",
					createDatasetSatRes(satResMap, flag));
		} else {
			chart = ChartFactory.createScatterPlot("Satellite-Residual vs Elevation Angle",
					"Elevation-Angle(in degrees)", "Satellite-Residual(in m)", createDatasetSatRes(satResMap, flag));
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

	public GraphPlotter(TreeMap<Long, double[]> map, String name) throws IOException {
		super(name);
		// TODO Auto-generated constructor stub

		final JFreeChart chart = ChartFactory.createXYLineChart(name, "GPS-time", name, createDatasetTrueObsv(map));
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

	public GraphPlotter(String name, ArrayList<Double> data, ArrayList<Long> timeList) throws IOException {
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

	public GraphPlotter(String name, HashMap<String, TreeMap<Integer, Double>> map) throws IOException {
		super(name);
		// TODO Auto-generated constructor stub

		final JFreeChart chart = ChartFactory.createXYLineChart(name, "GPS-time", name, createAnalyseDataset(map));
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

	public static void graphPostUnitW(ArrayList<Double> data, ArrayList<Long> timeList) throws IOException {
		GraphPlotter chart = new GraphPlotter("Posteriori Variance of Unit Weight", data, timeList);
		chart.pack();
		RefineryUtilities.positionFrameRandomly(chart);
		chart.setVisible(true);
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

	public static void graphSatRes(HashMap<String, HashMap<Integer, ArrayList<SatResidual>>> satResMap) {
		for (String key : satResMap.keySet()) {
			GraphPlotter chart = new GraphPlotter(satResMap.get(key), key, true);
			chart.pack();
			RefineryUtilities.positionFrameRandomly(chart);
			chart.setVisible(true);

			chart = new GraphPlotter(satResMap.get(key), key, false);
			chart.pack();
			RefineryUtilities.positionFrameRandomly(chart);
			chart.setVisible(true);
		}

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

	private XYDataset createDataset3dErr(HashMap<String, ArrayList<double[]>> dataMap, ArrayList<Long> timeList) {
		XYSeriesCollection dataset = new XYSeriesCollection();
		for (String key : dataMap.keySet()) {
			final XYSeries series = new XYSeries(key);
			ArrayList<double[]> list = dataMap.get(key);
			for (int i = 0; i < list.size(); i++) {
				double[] data = list.get(i);
				double err = Math.sqrt(Arrays.stream(data).map(j -> j * j).sum());
				series.add(timeList.get(i) - timeList.get(0), Double.valueOf(err));
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

	private XYDataset createDatasetTrueObsv(TreeMap<Long, double[]> map) {
		XYSeriesCollection dataset = new XYSeriesCollection();
		final XYSeries x = new XYSeries("X");
		final XYSeries y = new XYSeries("Y");
		final XYSeries z = new XYSeries("Z");
		Iterator<Entry<Long, double[]>> iterator = map.entrySet().iterator();

		while (iterator.hasNext()) {
			Entry<Long, double[]> entry = iterator.next();
			long time = entry.getKey();
			double[] data = entry.getValue();

			x.add(time, data[0]);
			y.add(time, data[1]);
			z.add(time, data[2]);

		}
		dataset.addSeries(x);
		dataset.addSeries(y);
		dataset.addSeries(z);

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

	private XYDataset createAnalyseDataset(HashMap<String, TreeMap<Integer, Double>> map) {
		XYSeriesCollection dataset = new XYSeriesCollection();
		for (String key : map.keySet()) {
			String name = key;
			final XYSeries series = new XYSeries(name);
			TreeMap<Integer, Double> data = map.get(key);
			for (int x : data.keySet()) {
				double y = data.get(x);
				series.add(x, y);
			}
			dataset.addSeries(series);
		}
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
				series.add(timeList.get(i) - timeList.get(0), Double.valueOf(data[i]));
			}
			dataset.addSeries(series);
		}

		return dataset;

	}

	private XYDataset createDatasetPostUnitW(ArrayList<Double> data, ArrayList<Long> timeList) {
		XYSeriesCollection dataset = new XYSeriesCollection();

		final XYSeries series = new XYSeries("Posteriori Variance of Unit Weight");
		double sum = 0;
		int count = 0;
		for (int i = 0; i < data.size(); i++) {
			if (data.get(i) == 0) {
				continue;
			}
			sum += data.get(i);
			count++;
			series.add(timeList.get(i) - timeList.get(0), data.get(i));
		}
		double avg = sum / count;
		final XYSeries mean = new XYSeries("Mean Posteriori Variance of Unit Weight: " + avg);

		for (int i = 0; i < data.size(); i++) {
			mean.add(timeList.get(i) - timeList.get(0), avg);
		}
		dataset.addSeries(series);
		dataset.addSeries(mean);
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

	private XYDataset createDatasetSatRes(HashMap<Integer, ArrayList<SatResidual>> dataMap, boolean flag) {
		XYSeriesCollection dataset = new XYSeriesCollection();

		for (int key : dataMap.keySet()) {
			final XYSeries series = new XYSeries(key);
			ArrayList<SatResidual> dataList = dataMap.get(key);
			if (flag) {
				double t0 = 0;
				for (int i = 0; i < dataList.size(); i++) {
					SatResidual data = dataList.get(i);
					double t = data.getT();
					double satRes = data.getResidual();
					if (t - t0 > 1) {
						series.add(t0, null);
					}
					series.add(t, satRes);
					t0 = t;
				}
			} else {
				for (int i = 0; i < dataList.size(); i++) {
					SatResidual data = dataList.get(i);
					double elevAngle = Math.toDegrees(data.getElevAngle());
					double satRes = data.getResidual();

					series.add(elevAngle, satRes);

				}
			}

			dataset.addSeries(series);
		}

		return dataset;
	}

}