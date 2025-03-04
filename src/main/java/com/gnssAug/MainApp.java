package com.gnssAug;

import java.time.Duration;
import java.time.Instant;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.HashMap;
import java.util.Set;
import java.util.TimeZone;
import java.util.TreeMap;

import org.apache.commons.math3.distribution.NormalDistribution;

import com.gnssAug.Android.Android;
import com.gnssAug.Android.Android_Static;
import com.gnssAug.Android.fileParser.GNSS_Log;
import com.gnssAug.Android.fileParser.GroundTruth_GSA;
import com.gnssAug.Android.models.GNSSLog;
import com.gnssAug.IGS.IGS;
import com.gnssAug.utility.LatLonUtil;

public class MainApp {

	public static void main(String[] args) throws Exception {

		Instant start = Instant.now();
		String base_url = "/Users/naman.agarwal/Library/CloudStorage/OneDrive-UniversityofCalgary/Google Decimeter Challenge/input_files/";// "C:\\Users\\naman.agarwal\\Documents\\GNSS\\Google
																																			// Decimeter
																																			// Challenge\\input_files\\";
		// "C:\\D drive\\Study\\Google Decimeter Challenge\\input_files\\";
		boolean isMac = true;
		String sep = isMac ? "/" : "\\";
		switch (3) {
		case 1:
			//String[] obsvCodeList = new String[] { "G5X", "E5X", "C2I" };
			String[] obsvCodeList = new String[] {"G1C","E1C","C2I"};
			String basePath = "/Users/naman.agarwal/Library/CloudStorage/OneDrive-UniversityofCalgary/Google Decimeter Challenge/decimeter/train/2021-03-10-US-SVL-1/SamsungS20Ultra";
			// "C:\\Users\\naman.agarwal\\Documents\\GNSS\\Google Decimeter
			// Challenge\\decimeter\\train\\2021-04-29-US-SJC-2\\SamsungS20Ultra";
			// "C:\\D drive\\Study\\Google Decimeter
			// Challenge\\decimeter\\train\\2021-04-29-US-MTV-1\\Pixel4";
			Set<String> discardSet = false ? Set.of("E11", "E25" ,"E7", "E8") :Set.of(""); //Set.of("C11", "G12", "G2", "G30");// C33
			String[] strList = basePath.split("/");
			String[] date = strList[strList.length - 2].split("-");
			int year = Integer.parseInt(date[0]);
			int month = Integer.parseInt(date[1]);
			int day = Integer.parseInt(date[2]);
			Calendar cal = Calendar.getInstance(TimeZone.getTimeZone("UTC"));
			cal.set(year, month - 1, day, 0, 0, 0);
			String doy = String.format("%03d", cal.get(Calendar.DAY_OF_YEAR));
			String mobName = strList[strList.length - 1];
			String obs_path = basePath + sep + "supplemental" + sep + mobName + "_GnssLog.20o";
			String derived_csv_path = basePath + sep + mobName + "_derived.csv";
			String gnss_log_path = basePath + sep + mobName + "_GnssLog.txt";
			String GTcsv = basePath + sep + "ground_truth.csv";
			String bias_path = base_url + year + "_" + doy + sep + "CAS0MGXRAP_" + year + doy + "0000_01D_01D_DCB.BSX";
			String clock_path = base_url + year + "_" + doy + sep + "COD0MGXFIN_" + year + doy + "0000_01D_30S_CLK.CLK";
			String orbit_path = base_url + year + "_" + doy + sep + "COD0MGXFIN_" + year + doy + "0000_01D_05M_ORB.SP3";
			String ionex_path = base_url + year + "_" + doy + sep + "igsg" + doy + "0.21I";
			Android.posEstimate(true, 0, 0, 18, obsvCodeList, derived_csv_path, gnss_log_path, GTcsv, bias_path,
					clock_path, orbit_path, ionex_path, true, true,true, false, false, discardSet,false);
			break;

		case 2:
			year = 2021;
			doy = "" + 119;
			bias_path = base_url + year + "_" + doy + "/CAS0MGXRAP_" + year + doy + "0000_01D_01D_DCB.BSX";
			clock_path = base_url + year + "_" + doy + "/COD0MGXFIN_" + year + doy + "0000_01D_30S_CLK.CLK";
			orbit_path = base_url + year + "_" + doy + "/COD0MGXFIN_" + year + doy + "0000_01D_05M_ORB.SP3";
			ionex_path = base_url + year + "_" + doy + "/igsg" + doy + "0.21I";
			String sinex_path = base_url + year + "_" + doy + "/igs21P21554.SNX";
			IGS.posEstimate(bias_path, clock_path, orbit_path, ionex_path, sinex_path, true, true, true, true,
					new String[] { "G1C", "E1C" }, 4, 1, 0, true, true, 2, true, true, false);
			break;
		case 3:
			//obsvCodeList = new String[] {"G5X","E5X" };
			obsvCodeList = new String[] {"G1C","E1C","C2I"};
			basePath = "/Users/naman.agarwal/Library/CloudStorage/OneDrive-UniversityofCalgary/Google Decimeter Challenge/Static/Test_20210713_T-A-SIS-01_open_sky_static_assisted_4h/RxX_Samsung_Galaxy_S20+_5G/RawData/gnss_log_2021_07_13_09_40_48.txt";
//			GNSS_Log.process(basePath);
//			TreeMap<Long, HashMap<String, ArrayList<GNSSLog>>> gnssLogMaps = GNSS_Log.getGnssLogMaps();

			discardSet = true ? Set.of("C22") : Set.of("C11", "G12", "G2", "G30");// C33
			strList = basePath.split("/");
			date = strList[strList.length - 1].split("_");
			year = Integer.parseInt(date[2]);
			month = Integer.parseInt(date[3]);
			day = Integer.parseInt(date[4]);
			cal = Calendar.getInstance(TimeZone.getTimeZone("UTC"));
			cal.set(year, month - 1, day, 0, 0, 0);
			doy = String.format("%03d", cal.get(Calendar.DAY_OF_YEAR));
			mobName = strList[strList.length - 3];

			bias_path = base_url + year + "_" + doy + sep + "CAS0MGXRAP_" + year + doy + "0000_01D_01D_DCB.BSX";
			clock_path = base_url + year + "_" + doy + sep + "COD0MGXFIN_" + year + doy + "0000_01D_30S_CLK.CLK";
			orbit_path = base_url + year + "_" + doy + sep + "COD0MGXFIN_" + year + doy + "0000_01D_05M_ORB.SP3";
			ionex_path = base_url + year + "_" + doy + sep + "igsg" + doy + "0.21I";
			// Urban Static
			double[] trueEcef = new double[] {4183777.518, 862857.121, 4721221.153};
			double[] llh = LatLonUtil.ecef2lla(trueEcef);
			
			
			Android_Static.posEstimate(true, 0, 0,2, obsvCodeList,basePath, trueEcef, bias_path,
					clock_path, orbit_path, ionex_path, true, true,false, false,false,discardSet,mobName);
			break;
		case 4:
			//obsvCodeList = new String[] { "G5X", "E5X", "C2I" };
			obsvCodeList = new String[] {"G1C","E1C","C2I"};
			basePath = "/Users/naman.agarwal/Library/CloudStorage/OneDrive-UniversityofCalgary/Google Decimeter Challenge/Dynamic/Test_20210716_T-A-SIS-09_suburban_bike_assisted_0.5h/RxX_Xiaomi_Mi_10T_Pro/RawData/gnss_log_2021_07_16_15_19_59.txt";

			discardSet = true ? Set.of("G4") : Set.of("C11", "G12", "G2", "G30");// C33
			strList = basePath.split("/");
			date = strList[strList.length - 1].split("_");
			year = Integer.parseInt(date[2]);
			month = Integer.parseInt(date[3]);
			day = Integer.parseInt(date[4]);
			cal = Calendar.getInstance(TimeZone.getTimeZone("UTC"));
			cal.set(year, month - 1, day, 0, 0, 0);
			doy = String.format("%03d", cal.get(Calendar.DAY_OF_YEAR));
			mobName = strList[strList.length - 3];
			GTcsv ="/Users/naman.agarwal/Library/CloudStorage/OneDrive-UniversityofCalgary/Google Decimeter Challenge/Dynamic/Test_20210716_T-A-SIS-09_suburban_bike_assisted_0.5h/Reference_trajectory/Test_20210716_T-A-SIS-09_suburban_bike_assisted_ECEF.txt";
			bias_path = base_url + year + "_" + doy + sep + "CAS0MGXRAP_" + year + doy + "0000_01D_01D_DCB.BSX";
			clock_path = base_url + year + "_" + doy + sep + "COD0MGXFIN_" + year + doy + "0000_01D_30S_CLK.CLK";
			orbit_path = base_url + year + "_" + doy + sep + "COD0MGXFIN_" + year + doy + "0000_01D_05M_ORB.SP3";
			ionex_path = base_url + year + "_" + doy + sep + "igsg" + doy + "0.21I";
			
			Android.posEstimate(true, 0, 0, 14, obsvCodeList, null, basePath, GTcsv, bias_path,
					clock_path, orbit_path, ionex_path, true, true, true, true, false, discardSet,true);
			break;
		case 5:
			obsvCodeList = new String[] { "G5X","E5X","C5X"};
			//obsvCodeList = new String[] {"G1C","E1C"};
			basePath = "/Users/naman.agarwal/Library/CloudStorage/OneDrive-UniversityofCalgary/Google Decimeter Challenge/Personal Data Collection/Static/ASCM419739/Pixel 7 pro/gnss_log_2025_01_28_22_05_38.txt";

			discardSet = true ? Set.of("E8") : Set.of("C11", "G12", "G2", "G30");// C33
			strList = basePath.split("/");
			date = strList[strList.length - 1].split("_");
			year = Integer.parseInt(date[2]);
			month = Integer.parseInt(date[3]);
			day = Integer.parseInt(date[4]);
			cal = Calendar.getInstance(TimeZone.getTimeZone("UTC"));
			cal.set(year, month - 1, day, 0, 0, 0);
			doy = String.format("%03d", cal.get(Calendar.DAY_OF_YEAR));
			mobName = strList[strList.length - 2];

			bias_path = base_url + year + "_" + doy + sep + "CAS0OPSRAP_" + year + doy + "0000_01D_01D_DCB.BIA";
			clock_path = base_url + year + "_" + doy + sep + "COD0MGXFIN_" + year + doy + "0000_01D_30S_CLK.CLK";
			orbit_path = base_url + year + "_" + doy + sep + "COD0MGXFIN_" + year + doy + "0000_01D_05M_ORB.SP3";
			ionex_path = base_url + year + "_" + doy + sep + "IGS0OPSFIN_"+ year + doy + "0000_01D_02H_GIM.INX";
			// Urban Static
			
			llh = new double[] {51.081628,-114.134081,1110.130-16.7243};
			trueEcef = LatLonUtil.lla2ecef(llh, true);
			
			Android_Static.posEstimate(true, 0, 0,21, obsvCodeList,basePath, trueEcef, bias_path,
					clock_path, orbit_path, ionex_path, true, true,true, false,false,discardSet,mobName);
			break;
		}
		Instant end = Instant.now();
		System.out.println("EXECUTION TIME -  " + Duration.between(start, end));

	}
}
