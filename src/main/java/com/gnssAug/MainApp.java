package com.gnssAug;

import java.time.Duration;
import java.time.Instant;
import java.util.Calendar;
import java.util.TimeZone;

import com.gnssAug.Android.Android;
import com.gnssAug.utility.LatLonUtil;

public class MainApp {

	public static void main(String[] args) {

		Instant start = Instant.now();

		switch (1) {
		case 1:
			String[] obsvCodeList = new String[] { "G1C" };
			String basePath = "C:\\D drive\\Study\\Google Decimeter Challenge\\decimeter\\train\\2021-04-29-US-MTV-1\\SamsungS20Ultra";
			String[] strList = basePath.split("\\\\");

			String[] date = strList[strList.length - 2].split("-");
			int year = Integer.parseInt(date[0]);
			int month = Integer.parseInt(date[1]);
			int day = Integer.parseInt(date[2]);
			Calendar cal = Calendar.getInstance(TimeZone.getTimeZone("UTC"));
			cal.set(year, month - 1, day, 0, 0, 0);
			int doy = cal.get(Calendar.DAY_OF_YEAR);
			String base_url = "C:\\D drive\\Study\\Google Decimeter Challenge\\input_files\\";
			String mobName = strList[strList.length - 1];
			String obs_path = basePath + "\\supplemental\\" + mobName + "_GnssLog.20o";
			String derived_csv_path = basePath + "\\" + mobName + "_derived.csv";
			String gnss_log_path = basePath + "\\" + mobName + "_GnssLog.txt";
			String GTcsv = basePath + "\\" + "ground_truth.csv";
			String bias_path = base_url + year + "_" + doy + "\\CAS0MGXRAP_" + year + doy + "0000_01D_01D_DCB.BSX";
			String clock_path = base_url + year + "_" + doy + "\\COD0MGXFIN_" + year + doy + "0000_01D_30S_CLK.CLK";
			String orbit_path = base_url + year + "_" + doy + "\\COD0MGXFIN_" + year + doy + "0000_01D_05M_ORB.SP3";
			Android.posEstimate(true, 0, 0, 3, obsvCodeList, derived_csv_path, gnss_log_path, GTcsv, bias_path,
					clock_path, orbit_path, false);
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
