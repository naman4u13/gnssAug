package com.gnssAug.Android.utility;

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

public class GraphPlotter extends ApplicationFrame {

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

		final JFreeChart chart = ChartFactory.createXYLineChart(chartTitle, "GPS-time", chartTitle,
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
		super(applicationTitle);
		// TODO Auto-generated constructor stub

		final JFreeChart chart = ChartFactory.createXYLineChart(chartTitle, "GPS-time", chartTitle,
				createDatasetENU(dataMap, timeList));
		final ChartPanel chartPanel = new ChartPanel(chart);
		chartPanel.setPreferredSize(new java.awt.Dimension(560, 370));
		chartPanel.setMouseZoomable(true, false);
		// ChartUtils.saveChartAsJPEG(new File(path + chartTitle + ".jpeg"), chart,
		// 1000, 600);
		setContentPane(chartPanel);

	}

	public GraphPlotter(String applicationTitle, HashMap<String, ArrayList<double[]>> dataMap) throws IOException {
		super(applicationTitle);
		// TODO Auto-generated constructor stub

		final JFreeChart chart = ChartFactory.createScatterPlot("2D Error", "East(m)", "North(m)",
				createDataset2dErr(dataMap));
		final ChartPanel chartPanel = new ChartPanel(chart);
		chartPanel.setPreferredSize(new java.awt.Dimension(560, 370));
		chartPanel.setMouseZoomable(true, false);
		// ChartUtils.saveChartAsJPEG(new File(path + "2d Error" + ".jpeg"), chart,
		// 1000, 600);
		setContentPane(chartPanel);

	}

	public static void graphENU(HashMap<String, ArrayList<double[]>> dataMap, ArrayList<Long> timeList)
			throws IOException {

		String[] chartNames = new String[] { "E", "N", "U" };
		for (int i = 0; i < 3; i++) {
			final int index = i;
			HashMap<String, double[]> subDataMap = new HashMap<String, double[]>();
			for (String key : dataMap.keySet()) {
				ArrayList<double[]> data = dataMap.get(key);
				double[] arr = data.stream().mapToDouble(j -> j[index]).toArray();
				subDataMap.put(key, arr);
			}
			GraphPlotter chart = new GraphPlotter("GPS PVT Error - ", chartNames[i] + "(m)", subDataMap, timeList,
					true);
			chart.pack();
			RefineryUtilities.positionFrameRandomly(chart);
			chart.setVisible(true);
		}
		GraphPlotter chart = new GraphPlotter("2D-Error", dataMap);
		chart.pack();
		RefineryUtilities.positionFrameRandomly(chart);
		chart.setVisible(true);
		chart = new GraphPlotter("3d-Error", "3d-Error", dataMap, timeList);
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
}