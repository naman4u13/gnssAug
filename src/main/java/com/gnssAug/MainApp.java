package com.gnssAug;

import java.time.Duration;
import java.time.Instant;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Set;
import java.util.TimeZone;
import com.gnssAug.Android.Android;
import com.gnssAug.Android.Android_Static;
import com.gnssAug.Android.Android_Static_Rinex;
import com.gnssAug.IGS.IGS;
import com.gnssAug.helper.ComputeSolidEarthTide;
import com.gnssAug.utility.LatLonUtil;
import com.gnssAug.utility.Time;

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
			// String[] obsvCodeList = new String[] { "G5X", "E5X", "C2I" };
			String[] obsvCodeList = new String[] { "G1C", "E1C", "C2I" };
			String basePath = "/Users/naman.agarwal/Library/CloudStorage/OneDrive-UniversityofCalgary/Google Decimeter Challenge/decimeter/train/2021-03-10-US-SVL-1/SamsungS20Ultra";
			// "C:\\Users\\naman.agarwal\\Documents\\GNSS\\Google Decimeter
			// Challenge\\decimeter\\train\\2021-04-29-US-SJC-2\\SamsungS20Ultra";
			// "C:\\D drive\\Study\\Google Decimeter
			// Challenge\\decimeter\\train\\2021-04-29-US-MTV-1\\Pixel4";
			Set<String> discardSet = false ? Set.of("E11", "E25", "E7", "E8") : Set.of(""); // Set.of("C11", "G12",
																							// "G2", "G30");// C33
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
			String dcb_bias_path = base_url + year + "_" + doy + sep + "CAS0MGXRAP_" + year + doy
					+ "0000_01D_01D_DCB.BSX";
			String osb_bias_path = base_url + year + "_" + doy + sep + "WUM0MGXFIN_" + year + doy
					+ "0000_01D_01D_OSB.BIA";
			String clock_path = base_url + year + "_" + doy + sep + "COD0MGXFIN_" + year + doy + "0000_01D_30S_CLK.CLK";
			String orbit_path = base_url + year + "_" + doy + sep + "COD0MGXFIN_" + year + doy + "0000_01D_05M_ORB.SP3";
			String ionex_path = base_url + year + "_" + doy + sep + "igsg" + doy + "0.21I";
			Android.posEstimate(true, 0, 0, 18, obsvCodeList, derived_csv_path, gnss_log_path, GTcsv, dcb_bias_path,
					clock_path, orbit_path, ionex_path, true, true, true, false, false, discardSet, false,false);
			break;

		case 2:
			year = 2024;
			doy = "" + 200;
//			bias_path = base_url + year + "_" + doy + "/CAS0MGXRAP_" + year + doy + "0000_01D_01D_DCB.BSX";
//			clock_path = base_url + year + "_" + doy + "/COD0MGXFIN_" + year + doy + "0000_01D_30S_CLK.CLK";
//			orbit_path = base_url + year + "_" + doy + "/COD0MGXFIN_" + year + doy + "0000_01D_05M_ORB.SP3";
//			ionex_path = base_url + year + "_" + doy + "/igsg" + doy + "0.21I";
			discardSet = Set.of("");// Set.of("E25","E33","C8", "C13", "C38");
			dcb_bias_path = base_url + year + "_" + doy + sep + "CAS0OPSRAP_" + year + doy + "0000_01D_01D_DCB.BIA";
			osb_bias_path = base_url + year + "_" + doy + sep + "WUM0MGXFIN_" + year + doy + "0000_01D_01D_OSB.BIA";
			clock_path = base_url + year + "_" + doy + sep + "WUM0MGXFIN_" + year + doy + "0000_01D_30S_CLK.CLK";
			orbit_path = base_url + year + "_" + doy + sep + "WUM0MGXFIN_" + year + doy + "0000_01D_05M_ORB.SP3";
//			osb_bias_path = base_url + year + "_" + doy + sep + "COD0MGXFIN_" + year + doy + "0000_01D_01D_OSB.BIA";
//			clock_path = base_url + year + "_" + doy + sep + "COD0MGXFIN_" + year + doy + "0000_01D_30S_CLK.CLK";
//			orbit_path = base_url + year + "_" + doy + sep + "COD0MGXFIN_" + year + doy + "0000_01D_05M_ORB.SP3";
			ionex_path = base_url + year + "_" + doy + sep + "IGS0OPSFIN_" + year + doy + "0000_01D_02H_GIM.INX";

			String sinex_path = base_url + year + "_" + doy + sep + "IGS0OPSSNX_" + year + doy + "0000_01D_01D_SOL.SNX";

//			String sinex_path = base_url + year + "_" + doy + "/igs21P21554.SNX";
			IGS.posEstimate(osb_bias_path, dcb_bias_path, clock_path, orbit_path, ionex_path, sinex_path, true, true,
					true, true, new String[] { "G1C","G5Q","E1C","E5Q"}, 4, 1, 0, true, true, 4, true, true, false,
					discardSet,false);
			break;
		case 3:
			// obsvCodeList = new String[] {"G5X","E5X" };
			obsvCodeList = new String[] { "G1C","E1C","C2I"};
			basePath = "/Users/naman.agarwal/Library/CloudStorage/OneDrive-UniversityofCalgary/Google Decimeter Challenge/Static/Test_20210728_T-A-SIS-10_urban_static_assisted_4h/RxX_Xiaomi_Redmi_Note_9T/RawData/gnss_log_2021_07_28_14_01_01.txt";
//			GNSS_Log.process(basePath);
//			TreeMap<Long, HashMap<String, ArrayList<GNSSLog>>> gnssLogMaps = GNSS_Log.getGnssLogMaps();

			discardSet = false ? Set.of("C2I22","C2I5", "G1C23", "E5X18", "E5X2", "E5X12", "E5X36", "E1C25", "C2I13", "E1C18") :Set.of("E1C20","E1C27","E1C14","E1C21","C2I39","C2I42","C2I45","C2I41","C2I24");//,"C2I6", "E5X31", "G5X1", "E5X7", "G5X8", "C2I23", "E5X25", "E5X24", "C2I34" );
			strList = basePath.split("/");
			date = strList[strList.length - 1].split("_");
			year = Integer.parseInt(date[2]);
			month = Integer.parseInt(date[3]);
			day = Integer.parseInt(date[4]);
			cal = Calendar.getInstance(TimeZone.getTimeZone("UTC"));
			cal.set(year, month - 1, day, 0, 0, 0);
			doy = String.format("%03d", cal.get(Calendar.DAY_OF_YEAR));
			mobName = strList[strList.length - 3];

			dcb_bias_path = base_url + year + "_" + doy + sep + "CAS0MGXRAP_" + year + doy + "0000_01D_01D_DCB.BSX";
			clock_path = base_url + year + "_" + doy + sep + "COD0MGXFIN_" + year + doy + "0000_01D_30S_CLK.CLK";
			orbit_path = base_url + year + "_" + doy + sep + "COD0MGXFIN_" + year + doy + "0000_01D_05M_ORB.SP3";
			ionex_path = base_url + year + "_" + doy + sep + "igsg" + doy + "0.21I";
			osb_bias_path = base_url + year + "_" + doy + sep + "WUM0MGXFIN_" + year + doy + "0000_01D_01D_OSB.BIA";
			// Urban Static
			double[] trueEcef = new double[] {4183748.339, 862806.185, 4721229.282};//{ 4183777.518, 862857.121, 4721221.153 };
			double[] llh = LatLonUtil.ecef2lla(trueEcef);
			Set<String> phaseDiscardSet = Set.of("C2I22","C2I5", "G1C23", "E5X18", "E5X2", "E5X12", "E5X36", "E1C25", "C2I13", "E1C18");
			Android_Static.posEstimate(true, 0, 0, 22, obsvCodeList, basePath, trueEcef, dcb_bias_path, clock_path,
					orbit_path, ionex_path, osb_bias_path, true, true, true, false, false, discardSet, mobName, false);
			break;
		case 4:
			// obsvCodeList = new String[] { "G5X", "E5X", "C2I" };
			obsvCodeList = new String[] { "G1C","G5X","E1C","E5X","C2I" };
			basePath = "/Users/naman.agarwal/Library/CloudStorage/OneDrive-UniversityofCalgary/Google Decimeter Challenge/Dynamic/Test_20210719_T-A-SIS-02_open_sky_pedestrian_assisted_0.5h/RxX_Samsung_Galaxy_S20+_5G/RawData/gnss_log_2021_07_19_16_16_18.txt";

			discardSet = true ? Set.of("") : Set.of("C11", "G12", "G2", "G30");// C33
			strList = basePath.split("/");
			date = strList[strList.length - 1].split("_");
			year = Integer.parseInt(date[2]);
			month = Integer.parseInt(date[3]);
			day = Integer.parseInt(date[4]);
			cal = Calendar.getInstance(TimeZone.getTimeZone("UTC"));
			cal.set(year, month - 1, day, 0, 0, 0);
			doy = String.format("%03d", cal.get(Calendar.DAY_OF_YEAR));
			mobName = strList[strList.length - 3];
			GTcsv = "/Users/naman.agarwal/Library/CloudStorage/OneDrive-UniversityofCalgary/Google Decimeter Challenge/Dynamic/Test_20210719_T-A-SIS-02_open_sky_pedestrian_assisted_0.5h/Reference_trajectory/Test_20210719_T-A-SIS-02_open_sky_pedestrian_assisted_ECEF.txt";
			dcb_bias_path = base_url + year + "_" + doy + sep + "CAS0MGXRAP_" + year + doy + "0000_01D_01D_DCB.BSX";
			clock_path = base_url + year + "_" + doy + sep + "COD0MGXFIN_" + year + doy + "0000_01D_30S_CLK.CLK";
			orbit_path = base_url + year + "_" + doy + sep + "COD0MGXFIN_" + year + doy + "0000_01D_05M_ORB.SP3";
			ionex_path = base_url + year + "_" + doy + sep + "igsg" + doy + "0.21I";

			Android.posEstimate(true, 0, 0, 22, obsvCodeList, null, basePath, GTcsv, dcb_bias_path, clock_path,
					orbit_path, ionex_path, true, true, true, false, false, discardSet, true,false);
			break;
		case 5:
			// obsvCodeList = new String[] { "G5X","E5X","C5X"};
			obsvCodeList = new String[] {"G1C","G5X","E1C","E5X","C2I"};
			basePath = "/Users/naman.agarwal/Library/CloudStorage/OneDrive-UniversityofCalgary/Google Decimeter Challenge/Personal Data Collection/Static/ASCM419739/Pixel 7 pro/gnss_log_2025_01_28_22_05_38.txt";

			discardSet = false ? Set.of("G5X10", "G1C30", "G5X23", "G5X27", "G1C7") : Set.of("C2I43", "C2I30", "E1C8","E5X5", "G5X27", "E1C27", "E5X27", "E5X26", "E1C5");
			strList = basePath.split("/");
			date = strList[strList.length - 1].split("_");
			year = Integer.parseInt(date[2]);
			month = Integer.parseInt(date[3]);
			day = Integer.parseInt(date[4]);
			cal = Calendar.getInstance(TimeZone.getTimeZone("UTC"));
			cal.set(year, month - 1, day, 0, 0, 0);
			doy = String.format("%03d", cal.get(Calendar.DAY_OF_YEAR));
			mobName = strList[strList.length - 2];

			dcb_bias_path = base_url + year + "_" + doy + sep + "CAS0OPSRAP_" + year + doy + "0000_01D_01D_DCB.BIA";
			clock_path = base_url + year + "_" + doy + sep + "WUM0MGXFIN_" + year + doy + "0000_01D_30S_CLK.CLK";
			orbit_path = base_url + year + "_" + doy + sep + "WUM0MGXFIN_" + year + doy + "0000_01D_05M_ORB.SP3";
			ionex_path = base_url + year + "_" + doy + sep + "IGS0OPSFIN_" + year + doy + "0000_01D_02H_GIM.INX";
			osb_bias_path = base_url + year + "_" + doy + sep + "WUM0MGXFIN_" + year + doy + "0000_01D_30S_OSB.BIA";

			llh = new double[] { 51.081628, -114.134081, 1110.130 - 16.7243 };
			trueEcef = LatLonUtil.lla2ecef(llh, true);
			Android_Static.posEstimate(true, 0, 0, 22, obsvCodeList, basePath, trueEcef, dcb_bias_path, clock_path,
					orbit_path, ionex_path, osb_bias_path, true, true, true, false, false, discardSet, mobName, true);
			break;
		case 6:
			// obsvCodeList = new String[] { "G5X","E5X","C5X"};
			obsvCodeList = new String[] { "G1C"};
			basePath = "/Users/naman.agarwal/Library/CloudStorage/OneDrive-UniversityofCalgary/Google Decimeter Challenge/Static/2024-10-15-Open/Pixel9XL.24o";

			discardSet = true ? Set.of("") : Set.of("E8", "C30", "C43");
			strList = basePath.split("/");
			date = strList[strList.length - 2].split("-");
			year = Integer.parseInt(date[0]);
			month = Integer.parseInt(date[1]);
			day = Integer.parseInt(date[2]);
			cal = Calendar.getInstance(TimeZone.getTimeZone("UTC"));
			cal.set(year, month - 1, day, 0, 0, 0);
			doy = String.format("%03d", cal.get(Calendar.DAY_OF_YEAR));
			
			dcb_bias_path = base_url + year + "_" + doy + sep + "CAS0OPSRAP_" + year + doy + "0000_01D_01D_DCB.BIA";
			clock_path = base_url + year + "_" + doy + sep + "WUM0MGXFIN_" + year + doy + "0000_01D_30S_CLK.CLK";
			orbit_path = base_url + year + "_" + doy + sep + "WUM0MGXFIN_" + year + doy + "0000_01D_05M_ORB.SP3";
			ionex_path = base_url + year + "_" + doy + sep + "IGS0OPSFIN_" + year + doy + "0000_01D_02H_GIM.INX";
			osb_bias_path = base_url + year + "_" + doy + sep + "WUM0MGXFIN_" + year + doy + "0000_01D_01D_OSB.BIA";

			llh = new double[] { 48.67878205,19.22444203,456.649 };
			trueEcef = LatLonUtil.lla2ecef(llh, true);

			Android_Static_Rinex.posEstimate(basePath,osb_bias_path, dcb_bias_path, clock_path, orbit_path, ionex_path, true,
					true, true, trueEcef, obsvCodeList, 4, 0, 0, true, true, 5, true, false, false, discardSet, false);
			break;

		}
		Instant end = Instant.now();
		System.out.println("EXECUTION TIME -  " + Duration.between(start, end));

	}
}
