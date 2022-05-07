package com.gnssAug;

import java.time.Duration;
import java.time.Instant;

import com.gnssAug.Android.Android;
import com.gnssAug.Android.helper.RotationMatrix;

public class MainApp {

	public static void main(String[] args) {

		Instant start = Instant.now();

		switch (1) {
		case 1:
			String[] obsvCodeList = new String[] { "G1C" };
			String basePath = "C:\\D drive\\Study\\Google Decimeter Challenge\\decimeter\\train\\2021-04-29-US-SJC-2\\Pixel4";
			String[] strList = basePath.split("\\\\");
			String MobName = strList[strList.length - 1];
			String obs_path = basePath + "\\supplemental\\" + MobName + "_GnssLog.20o";
			String derived_csv_path = basePath + "\\" + MobName + "_derived.csv";
			String gnss_log_path = basePath + "\\" + MobName + "_GnssLog.txt";
			String GTcsv = basePath + "\\" + "ground_truth.csv";
			Android.posEstimate(true, 5, -1, 4, obsvCodeList, derived_csv_path, gnss_log_path, GTcsv);
			break;

		case 2:

			double[][] test = new double[][] { { 2, 3 }, { -3, 5 }, { -4, -1 }, { 5, -2 } };
			for (int i = 0; i < 4; i++) {
				RotationMatrix.check(test[i][0], test[i][1]);
			}
			break;

		}
		Instant end = Instant.now();
		System.out.println("EXECUTION TIME -  " + Duration.between(start, end));

	}
}
