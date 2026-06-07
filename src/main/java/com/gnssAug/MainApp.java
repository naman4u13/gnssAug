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
import com.gnssAug.Android.constants.EstimatorMode;
import com.gnssAug.IGS.IGS;
import com.gnssAug.helper.ComputeSolidEarthTide;
import com.gnssAug.utility.LatLonUtil;
import com.gnssAug.utility.Time;

public class MainApp {

	static final int GDC_DYNAMIC          = 1; // Google Decimeter Challenge, dynamic
	static final int IGS_RINEX            = 2; // IGS reference station, RINEX
	static final int ANDROID_STATIC       = 3; // T-A-SIS-01, open-sky static, Samsung S20+
	static final int ANDROID_PEDESTRIAN   = 4; // T-A-SIS-02, open-sky pedestrian, Samsung S20+
	static final int PERSONAL_STATIC      = 5; // personal data collection, static (Calgary)
	static final int RINEX_ANDROID_STATIC = 6; // RINEX-format Android log, static

	static final int RUN = PERSONAL_STATIC;

	public static void main(String[] args) throws Exception {

		Instant start = Instant.now();
		String base_url = "/Users/naman.agarwal/Library/CloudStorage/OneDrive-UniversityofCalgary/Google Decimeter Challenge/input_files/";
		String sep = "/";
		switch (RUN) {

		case GDC_DYNAMIC: {
			// Samsung S20 Ultra | GDC train 2021-03-10 US-SVL-1 | G1C+E1C+C2I | EKF_TDCP
			String[] obsvCodeList = new String[] { "G1C", "E1C", "C2I" };
			String basePath = "/Users/naman.agarwal/Library/CloudStorage/OneDrive-UniversityofCalgary/Google Decimeter Challenge/decimeter/train/2021-03-10-US-SVL-1/SamsungS20Ultra";
			boolean useDiscardSet = false;
			Set<String> discardSet = useDiscardSet ? Set.of("E11", "E25", "E7", "E8") : Set.of("");
			String[] yd    = dateFromFolder(basePath, 2); // folder: 2021-03-10-US-SVL-1
			String year    = yd[0], doy = yd[1];
			String mobName = basePath.split("/")[basePath.split("/").length - 1];
			String derived_csv_path = basePath + sep + mobName + "_derived.csv";
			String gnss_log_path    = basePath + sep + mobName + "_GnssLog.txt";
			String GTcsv            = basePath + sep + "ground_truth.csv";
			String dir = pDir(base_url, year, doy);
			String dcb_bias_path = dir + "CAS0MGXRAP_" + year + doy + "0000_01D_01D_DCB.BSX";
			String osb_bias_path = dir + "WUM0MGXFIN_" + year + doy + "0000_01D_01D_OSB.BIA";
			String clock_path    = dir + "COD0MGXFIN_" + year + doy + "0000_01D_30S_CLK.CLK";
			String orbit_path    = dir + "COD0MGXFIN_" + year + doy + "0000_01D_05M_ORB.SP3";
			String ionex_path    = dir + "igsg"         + doy + "0.21I";
			Android.posEstimate(
					/*doPosErrPlot*/ true, /*elMask*/ 0, /*snrMask*/ 0,
					EstimatorMode.EKF_TDCP, obsvCodeList, derived_csv_path, gnss_log_path, GTcsv,
					dcb_bias_path, clock_path, orbit_path, ionex_path,
					/*useIGS*/ true, /*doAnalyze*/ true, /*doTest*/ true,
					/*outlierAnalyze*/ false, /*mapDeltaRanges*/ false,
					discardSet, /*repairCS*/ false, /*doTimeSlice*/ false);
			break;
		}

		case IGS_RINEX: {
			// IGS station | 2024 DOY 200 | G1C+G5Q | WUM products
			String year = "2024", doy = "200";
			Set<String> discardSet = Set.of("");
			String dir = pDir(base_url, year, doy);
			String dcb_bias_path = dir + "CAS0OPSRAP_" + year + doy + "0000_01D_01D_DCB.BIA";
			String osb_bias_path = dir + "WUM0MGXFIN_" + year + doy + "0000_01D_01D_OSB.BIA";
			String clock_path    = dir + "WUM0MGXFIN_" + year + doy + "0000_01D_30S_CLK.CLK";
			String orbit_path    = dir + "WUM0MGXFIN_" + year + doy + "0000_01D_05M_ORB.SP3";
			String ionex_path    = dir + "IGS0OPSFIN_" + year + doy + "0000_01D_02H_GIM.INX";
			String sinex_path    = dir + "IGS0OPSSNX_" + year + doy + "0000_01D_01D_SOL.SNX";
			IGS.posEstimate(osb_bias_path, dcb_bias_path, clock_path, orbit_path, ionex_path, sinex_path,
					true, true, true, true, new String[] { "G1C", "G5Q" }, 4, 1, 0,
					true, true, 5, true, false, false, discardSet, true);
			break;
		}

		case ANDROID_STATIC: {
			// Samsung S20+ | T-A-SIS-01 open-sky static 4h, 2021-07-13 | G1C+E1C+G5X+E5X+C2I | PPP_FLOAT
			String[] obsvCodeList = new String[] { "G1C", "E1C", "G5X", "E5X", "C2I" };
			String basePath = "/Users/naman.agarwal/Library/CloudStorage/OneDrive-UniversityofCalgary/Google Decimeter Challenge/Static/Test_20210713_T-A-SIS-01_open_sky_static_assisted_4h/RxX_Samsung_Galaxy_S20+_5G/RawData/gnss_log_2021_07_13_09_40_48.txt";
			boolean useDiscardSet = true;
			Set<String> discardSet = useDiscardSet
					? Set.of("E5X2", "E5X18", "E5X36", "C2I22", "C2I13", "E5X12", "E1C18", "G1C23", "E1C25")
					: Set.of("");
			String[] yd    = dateFromGnssLog(basePath, 0);
			String year    = yd[0], doy = yd[1];
			String mobName = basePath.split("/")[basePath.split("/").length - 3];
			String dir = pDir(base_url, year, doy);
			String dcb_bias_path = dir + "CAS0MGXRAP_" + year + doy + "0000_01D_01D_DCB.BSX";
			String clock_path    = dir + "COD0MGXFIN_" + year + doy + "0000_01D_30S_CLK.CLK";
			String orbit_path    = dir + "COD0MGXFIN_" + year + doy + "0000_01D_05M_ORB.SP3";
			String ionex_path    = dir + "igsg"         + doy + "0.21I";
			String osb_bias_path = dir + "WUM0MGXFIN_"  + year + doy + "0000_01D_01D_OSB.BIA";
			double[] trueEcef    = new double[] { 4183777.518, 862857.121, 4721221.153 };
			Android_Static.posEstimate(
					/*doPosErrPlot*/ true, /*elMask*/ 0, /*snrMask*/ 0,
					EstimatorMode.PPP_FLOAT, obsvCodeList, basePath, trueEcef,
					dcb_bias_path, clock_path, orbit_path, ionex_path, osb_bias_path,
					/*useIGS*/ true, /*doAnalyze*/ true, /*doTest*/ true,
					/*outlierAnalyze*/ false, /*mapDeltaRanges*/ false,
					discardSet, mobName, /*repairCS*/ true, /*doTimeSlice*/ false);
			break;
		}

		case ANDROID_PEDESTRIAN: {
			// Samsung S20+ | T-A-SIS-02 open-sky pedestrian 0.5h, 2021-07-19 | G1C+G5X+E1C+E5X+C2I | PPP_FLOAT
			String[] obsvCodeList = new String[] { "G1C", "G5X", "E1C", "E5X", "C2I" };
			String basePath = "/Users/naman.agarwal/Library/CloudStorage/OneDrive-UniversityofCalgary/Google Decimeter Challenge/Dynamic/Test_20210719_T-A-SIS-02_open_sky_pedestrian_assisted_0.5h/RxX_Samsung_Galaxy_S20+_5G/RawData/gnss_log_2021_07_19_16_16_18.txt";
			String GTcsv    = "/Users/naman.agarwal/Library/CloudStorage/OneDrive-UniversityofCalgary/Google Decimeter Challenge/Dynamic/Test_20210719_T-A-SIS-02_open_sky_pedestrian_assisted_0.5h/Reference_trajectory/Test_20210719_T-A-SIS-02_open_sky_pedestrian_assisted_ECEF.txt";
			boolean useDiscardSet = false;
			Set<String> discardSet = useDiscardSet ? Set.of("C11", "G12", "G2", "G30") : Set.of("");
			String[] yd = dateFromGnssLog(basePath, 0);
			String year = yd[0], doy = yd[1];
			String dir = pDir(base_url, year, doy);
			String dcb_bias_path = dir + "CAS0MGXRAP_" + year + doy + "0000_01D_01D_DCB.BSX";
			String clock_path    = dir + "COD0MGXFIN_" + year + doy + "0000_01D_30S_CLK.CLK";
			String orbit_path    = dir + "COD0MGXFIN_" + year + doy + "0000_01D_05M_ORB.SP3";
			String ionex_path    = dir + "igsg"         + doy + "0.21I";
			Android.posEstimate(
					/*doPosErrPlot*/ true, /*elMask*/ 0, /*snrMask*/ 0,
					EstimatorMode.PPP_FLOAT, obsvCodeList, null, basePath, GTcsv,
					dcb_bias_path, clock_path, orbit_path, ionex_path,
					/*useIGS*/ true, /*doAnalyze*/ true, /*doTest*/ true,
					/*outlierAnalyze*/ false, /*mapDeltaRanges*/ false,
					discardSet, /*repairCS*/ true, /*doTimeSlice*/ false);
			break;
		}

		case PERSONAL_STATIC: {
			// Pixel 4 | ASCM Calgary, 2025-11-21 | G1C+G5X+E1C+E5X | PPP_FLOAT | repair OFF
			String[] obsvCodeList = new String[] { "G1C", "G5X", "E1C", "E5X" };
			String basePath = "/Users/naman.agarwal/Library/CloudStorage/OneDrive-UniversityofCalgary/Google Decimeter Challenge/Personal Data Collection/Static/ASCM419739/Nov212025/Pixel 4/gnss_log_2025_11_21_17_32_25.txt";
			boolean useDiscardSet = true;
			Set<String> discardSet = useDiscardSet ? Set.of("G1C10") : Set.of("");
			String[] yd = dateFromGnssLog(basePath, 7); // Calgary MST = UTC+7
			String year = yd[0], doy = yd[1];
			String mobName    = basePath.split("/")[basePath.split("/").length - 2];
			double[] trueEcef = LatLonUtil.lla2ecef(new double[] { 51.081628, -114.134081, 1110.130 - 16.7243+1.195 }, true);
			String dir = pDir(base_url, year, doy);
			String dcb_bias_path = dir + "CAS1OPSRAP_" + year + doy + "0000_01D_01D_DCB.BIA";
			String clock_path    = dir + "COD0MGXFIN_" + year + doy + "0000_01D_30S_CLK.CLK";
			String orbit_path    = dir + "COD0MGXFIN_" + year + doy + "0000_01D_05M_ORB.SP3";
			String ionex_path    = dir + "IGS0OPSFIN_" + year + doy + "0000_01D_02H_GIM.INX";
			String osb_bias_path = dir + "WUM0MGXNRT_" + year + doy + "0000_01D_01D_OSB.BIA";
			Android_Static.posEstimate(
					/*doPosErrPlot*/ true, /*elMask*/ 0, /*snrMask*/ 0,
					EstimatorMode.PPP_FLOAT, obsvCodeList, basePath, trueEcef,
					dcb_bias_path, clock_path, orbit_path, ionex_path, osb_bias_path,
					/*useIGS*/ true, /*doAnalyze*/ true, /*doTest*/ true,
					/*outlierAnalyze*/ false, /*mapDeltaRanges*/ false,
					discardSet, mobName, /*repairCS*/ false, /*doTimeSlice*/ false);
			break;
		}

		case RINEX_ANDROID_STATIC: {
			// Pixel 9 XL | 2024-10-15 open-sky | G1C | RINEX_PPP_LOWCOST | WUM products
			String[] obsvCodeList = new String[] { "G1C" };
			String basePath = "/Users/naman.agarwal/Library/CloudStorage/OneDrive-UniversityofCalgary/Google Decimeter Challenge/Static/2024-10-15-Open/Pixel9XL.24o";
			boolean useDiscardSet = false;
			Set<String> discardSet = useDiscardSet ? Set.of("E8", "C30", "C43") : Set.of("");
			String[] yd = dateFromFolder(basePath, 2); // folder: 2024-10-15-Open
			String year = yd[0], doy = yd[1];
			double[] trueEcef    = LatLonUtil.lla2ecef(new double[] { 48.67878205, 19.22444203, 456.649 }, true);
			String dir = pDir(base_url, year, doy);
			String dcb_bias_path = dir + "CAS0OPSRAP_" + year + doy + "0000_01D_01D_DCB.BIA";
			String clock_path    = dir + "WUM0MGXFIN_" + year + doy + "0000_01D_30S_CLK.CLK";
			String orbit_path    = dir + "WUM0MGXFIN_" + year + doy + "0000_01D_05M_ORB.SP3";
			String ionex_path    = dir + "IGS0OPSFIN_" + year + doy + "0000_01D_02H_GIM.INX";
			String osb_bias_path = dir + "WUM0MGXFIN_" + year + doy + "0000_01D_01D_OSB.BIA";
			Android_Static_Rinex.posEstimate(basePath, osb_bias_path, dcb_bias_path, clock_path, orbit_path, ionex_path,
					true, true, true, trueEcef, obsvCodeList, 4, 0, 0,
					true, true, EstimatorMode.RINEX_PPP_LOWCOST, true, false, false, discardSet, false);
			break;
		}

		}
		Instant end = Instant.now();
		Duration _d = Duration.between(start, end);
		long _mins = _d.toMinutes();
		double _secs = (_d.toMillis() - _mins * 60000) / 1000.0;
		String _timeStr = _mins > 0
				? _mins + " min " + String.format("%.1f", _secs) + " s"
				: String.format("%.1f", _d.toMillis() / 1000.0) + " s";
		System.out.println("Execution Time : " + _timeStr);

	}

	/** Returns {year, doy} as Strings from a gnss_log_YYYY_MM_DD_HH_mm_ss.txt path, shifted by utcOffset hours. */
	private static String[] dateFromGnssLog(String gnssLogPath, int utcOffset) {
		String[] tok  = gnssLogPath.split("/");
		String[] date = tok[tok.length - 1].split("_");
		Calendar c = Calendar.getInstance(TimeZone.getTimeZone("UTC"));
		c.set(Integer.parseInt(date[2]), Integer.parseInt(date[3]) - 1, Integer.parseInt(date[4]),
			  Integer.parseInt(date[5]) + utcOffset, 0, 0);
		return new String[] { String.valueOf(c.get(Calendar.YEAR)),
							  String.format("%03d", c.get(Calendar.DAY_OF_YEAR)) };
	}

	/** Returns {year, doy} as Strings from a YYYY-MM-DD folder segment in the path. segFromEnd=2 means the parent folder. */
	private static String[] dateFromFolder(String path, int segFromEnd) {
		String[] parts = path.split("/");
		String[] date  = parts[parts.length - segFromEnd].split("-");
		Calendar c = Calendar.getInstance(TimeZone.getTimeZone("UTC"));
		c.set(Integer.parseInt(date[0]), Integer.parseInt(date[1]) - 1, Integer.parseInt(date[2]), 0, 0, 0);
		return new String[] { String.valueOf(c.get(Calendar.YEAR)),
							  String.format("%03d", c.get(Calendar.DAY_OF_YEAR)) };
	}

	/** Returns the product directory prefix: baseUrl/YYYY_DOY/ */
	private static String pDir(String baseUrl, String year, String doy) {
		return baseUrl + year + "_" + doy + "/";
	}
}
