package com.gnssAug;

import java.time.Duration;
import java.time.Instant;

import com.gnssAug.Android.Android;
import com.gnssAug.Android.utility.LatLonUtil;

public class MainApp {

	public static void main(String[] args) {

		Instant start = Instant.now();

		switch (1) {
		case 1:
			String[] obsvCodeList = new String[] { "G1C" };
			String basePath = "C:\\D drive\\Study\\Google Decimeter Challenge\\decimeter\\train\\2021-04-29-US-MTV-1\\Pixel5";
			String[] strList = basePath.split("\\\\");
			String MobName = strList[strList.length - 1];
			String obs_path = basePath + "\\supplemental\\" + MobName + "_GnssLog.20o";
			String derived_csv_path = basePath + "\\" + MobName + "_derived.csv";
			String gnss_log_path = basePath + "\\" + MobName + "_GnssLog.txt";
			String GTcsv = basePath + "\\" + "ground_truth.csv";
			Android.posEstimate(true, 0, 0, 3, obsvCodeList, derived_csv_path, gnss_log_path, GTcsv);
			break;

		case 2:
			double[] ecefr = LatLonUtil.lla2ecef(new double[] { 45.9132, 36.7484, 1877753.2 }, true);
			double[] enu = new double[] { 355601.3, -923083.2, 1041016.4 };
			double[] ecef = LatLonUtil.enu2ecef(enu, ecefr, true);
			System.out.println(ecef);
			break;

		}
		Instant end = Instant.now();
		System.out.println("EXECUTION TIME -  " + Duration.between(start, end));

	}
}
