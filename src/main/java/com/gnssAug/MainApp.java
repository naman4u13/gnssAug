package com.gnssAug;

import java.time.Duration;
import java.time.Instant;
import java.util.Calendar;
import java.util.TimeZone;

import com.gnssAug.Android.Android;
import com.gnssAug.IGS.IGS;

public class MainApp {

	public static void main(String[] args) {

		Instant start = Instant.now();
		String base_url = "C:\\Users\\naman.agarwal\\Downloads\\GNSS\\Google Decimeter Challenge\\input_files\\";
		//"C:\\D drive\\Study\\Google Decimeter Challenge\\input_files\\";
		switch (2) {
		case 1:
			String[] obsvCodeList = new String[] { "G1C", "E1C" };

			String basePath = "C:\\Users\\naman.agarwal\\Downloads\\GNSS\\Google Decimeter Challenge\\decimeter\\train\\2021-04-28-US-MTV-1\\Pixel5";
			//"C:\\D drive\\Study\\Google Decimeter Challenge\\decimeter\\train\\2021-04-28-US-MTV-1\\Pixel5";
			String[] strList = basePath.split("\\\\");

			String[] date = strList[strList.length - 2].split("-");
			int year = Integer.parseInt(date[0]);
			int month = Integer.parseInt(date[1]);
			int day = Integer.parseInt(date[2]);
			Calendar cal = Calendar.getInstance(TimeZone.getTimeZone("UTC"));
			cal.set(year, month - 1, day, 0, 0, 0);

			String doy = String.format("%03d", cal.get(Calendar.DAY_OF_YEAR));

			String mobName = strList[strList.length - 1];
			String obs_path = basePath + "\\supplemental\\" + mobName + "_GnssLog.20o";
			String derived_csv_path = basePath + "\\" + mobName + "_derived.csv";
			String gnss_log_path = basePath + "\\" + mobName + "_GnssLog.txt";
			String GTcsv = basePath + "\\" + "ground_truth.csv";
			String bias_path = base_url + year + "_" + doy + "\\CAS0MGXRAP_" + year + doy + "0000_01D_01D_DCB.BSX";
			String clock_path = base_url + year + "_" + doy + "\\COD0MGXFIN_" + year + doy + "0000_01D_30S_CLK.CLK";
			String orbit_path = base_url + year + "_" + doy + "\\COD0MGXFIN_" + year + doy + "0000_01D_05M_ORB.SP3";
			String ionex_path = base_url + year + "_" + doy + "\\igsg" + doy + "0.21I";
			Android.posEstimate(true, 0, 0, 6, obsvCodeList, derived_csv_path, gnss_log_path, GTcsv, bias_path,
					clock_path, orbit_path, ionex_path, false, true, false);
			break;

		case 2:
			year = 2021;
			doy = "" + 119;
			bias_path = base_url + year + "_" + doy + "\\CAS0MGXRAP_" + year + doy + "0000_01D_01D_DCB.BSX";
			clock_path = base_url + year + "_" + doy + "\\COD0MGXFIN_" + year + doy + "0000_01D_30S_CLK.CLK";
			orbit_path = base_url + year + "_" + doy + "\\COD0MGXFIN_" + year + doy + "0000_01D_05M_ORB.SP3";
			ionex_path = base_url + year + "_" + doy + "\\igsg" + doy + "0.21I";
			String sinex_path = base_url + year + "_" + doy + "\\igs21P21554.SNX";
			IGS.posEstimate(bias_path, clock_path, orbit_path, ionex_path, sinex_path, true, true, true, true,
					new String[] { "G1C" }, 4, 1, 0, true, true, 4, true, true, false);
			break;

		}
		Instant end = Instant.now();
		System.out.println("EXECUTION TIME -  " + Duration.between(start, end));

	}
}
