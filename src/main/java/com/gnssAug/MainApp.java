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

/**
 * Application entry point for the gnssAug GNSS positioning pipeline.
 *
 * <p>Select one of the six {@code RUN} modes by changing the {@code RUN} constant.
 * Each mode corresponds to a specific dataset, receiver type, and estimator:
 *
 * <table border="1" summary="Run mode descriptions">
 *   <tr><th>Constant</th><th>Dataset</th><th>Receiver</th><th>Estimator</th></tr>
 *   <tr><td>{@link #GDC_DYNAMIC}</td>
 *       <td>Google Decimeter Challenge, US-SVL-1, 2021-03-10</td>
 *       <td>Samsung S20 Ultra</td><td>EKF_TDCP (dynamic)</td></tr>
 *   <tr><td>{@link #IGS_RINEX}</td>
 *       <td>IGS station AJAC, 2024 DOY 200</td>
 *       <td>Geodetic (RINEX)</td><td>IGS_PPP_FLOAT (static)</td></tr>
 *   <tr><td>{@link #ANDROID_STATIC}</td>
 *       <td>T-A-SIS-01 open-sky static 4h, 2021-07-13</td>
 *       <td>Samsung S20+</td><td>PPP_FLOAT (static)</td></tr>
 *   <tr><td>{@link #ANDROID_PEDESTRIAN}</td>
 *       <td>T-A-SIS-02 open-sky pedestrian 0.5h, 2021-07-19</td>
 *       <td>Samsung S20+</td><td>PPP_FLOAT (kinematic)</td></tr>
 *   <tr><td>{@link #PERSONAL_STATIC}</td>
 *       <td>ASCM Calgary, 2025-11-21</td>
 *       <td>Pixel 4</td><td>PPP_FLOAT (static)</td></tr>
 *   <tr><td>{@link #RINEX_ANDROID_STATIC}</td>
 *       <td>Pixel 9 XL open-sky, 2024-10-15 (RINEX format)</td>
 *       <td>Pixel 9 XL</td><td>RINEX_PPP_LOWCOST (static)</td></tr>
 * </table>
 *
 * <p><b>Product file naming convention (MGEX):</b> all precise orbit, clock, and bias
 * files follow the IGS long filename format: {@code AGENCY_TYPE_YYYYDDD0000_01D_RES_TYPE.EXT}.
 * Product directories are resolved by {@link #pDir} from {@code base_url}, year, and DOY.
 */
public class MainApp {

	/** Google Decimeter Challenge dynamic dataset — Samsung S20 Ultra, EKF_TDCP. */
	static final int GDC_DYNAMIC          = 1;
	/** IGS reference station AJAC, RINEX input — geodetic receiver, IGS_PPP_FLOAT. */
	static final int IGS_RINEX            = 2;
	/** T-A-SIS-01 open-sky static 4h — Samsung S20+, PPP_FLOAT. */
	static final int ANDROID_STATIC       = 3;
	/** T-A-SIS-02 open-sky pedestrian 0.5h — Samsung S20+, PPP_FLOAT. */
	static final int ANDROID_PEDESTRIAN   = 4;
	/** Personal static data collection (Calgary ASCM) — Pixel 4, PPP_FLOAT. */
	static final int PERSONAL_STATIC      = 5;
	/** RINEX-format Android log, static — Pixel 9 XL, RINEX_PPP_LOWCOST. */
	static final int RINEX_ANDROID_STATIC = 6;

	/**
	 * Active run mode. Change this constant to switch datasets.
	 * All other configuration (file paths, signals, estimator flags) is embedded
	 * in the corresponding {@code case} block inside {@link #main}.
	 */
	static final int RUN = PERSONAL_STATIC;

	/**
	 * Application entry point. Selects the active dataset via {@link #RUN}, constructs
	 * all file paths and configuration flags, invokes the appropriate positioning
	 * pipeline ({@link IGS}, {@link com.gnssAug.Android.Android_Static}, etc.),
	 * and prints total wall-clock execution time on completion.
	 */
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
			String igs_base   = "/Users/naman.agarwal/Library/CloudStorage/OneDrive-UniversityofCalgary/input_files";
			String obs_path   = igs_base + "/Highrate/AJAC00FRA_S_20242000530_15M_01S_MO.rnx";
			String nav_path   = igs_base + "/BRDC00IGS_R_20201000000_01D_MN.rnx/BRDC00IGS_R_20201000000_01D_MN.rnx";
			String dcb_bias_path = dir + "CAS0OPSRAP_" + year + doy + "0000_01D_01D_DCB.BIA";
			String osb_bias_path = dir + "WUM0MGXFIN_" + year + doy + "0000_01D_01D_OSB.BIA";
			String clock_path    = dir + "WUM0MGXFIN_" + year + doy + "0000_01D_30S_CLK.CLK";
			String orbit_path    = dir + "WUM0MGXFIN_" + year + doy + "0000_01D_05M_ORB.SP3";
			String ionex_path    = dir + "IGS0OPSFIN_" + year + doy + "0000_01D_02H_GIM.INX";
			String sinex_path    = dir + "IGS0OPSSNX_" + year + doy + "0000_01D_01D_SOL.SNX";
			IGS.posEstimate(obs_path, nav_path, osb_bias_path, dcb_bias_path, clock_path, orbit_path,
					ionex_path, sinex_path,
					/*useBias*/ true, /*useGIM*/ true, /*useIGS*/ true, /*useSNX*/ true,
					new String[] { "G1C", "G5Q" }, /*minSat*/ 4, /*cutOffAng*/ 1, /*snrMask*/ 0,
					/*corrIono*/ true, /*corrTropo*/ true, EstimatorMode.IGS_PPP_FLOAT,
					/*doAnalyze*/ true, /*doTest*/ false, /*outlierAnalyze*/ false,
					discardSet, /*repairCS*/ true, /*predictPhaseClock*/ false);
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

	/**
	 * Extracts the GPS year and day-of-year (DOY) from an Android GnssLog filename.
	 *
	 * <p>Expected filename format: {@code gnss_log_YYYY_MM_DD_HH_mm_ss.txt}.
	 * The {@code utcOffset} is added to the hour before computing the DOY, allowing
	 * local-time log filenames to be correctly converted to UTC-based GPS day.
	 * For example, Calgary MST (UTC−7) passes {@code utcOffset = 7}.
	 *
	 * @param gnssLogPath full path to the GnssLog file
	 * @param utcOffset   hours to add to the log timestamp to obtain UTC (e.g. 7 for MST)
	 * @return {@code String[]{year, doy}} where {@code doy} is zero-padded to 3 digits
	 */
	private static String[] dateFromGnssLog(String gnssLogPath, int utcOffset) {
		String[] tok  = gnssLogPath.split("/");
		String[] date = tok[tok.length - 1].split("_");
		Calendar c = Calendar.getInstance(TimeZone.getTimeZone("UTC"));
		c.set(Integer.parseInt(date[2]), Integer.parseInt(date[3]) - 1, Integer.parseInt(date[4]),
			  Integer.parseInt(date[5]) + utcOffset, 0, 0);
		return new String[] { String.valueOf(c.get(Calendar.YEAR)),
							  String.format("%03d", c.get(Calendar.DAY_OF_YEAR)) };
	}

	/**
	 * Extracts the GPS year and day-of-year (DOY) from a {@code YYYY-MM-DD} folder
	 * segment embedded in the file path.
	 *
	 * <p>{@code segFromEnd} controls which path segment to parse, counted from the end.
	 * For example, {@code segFromEnd = 2} picks the parent directory of the file,
	 * which is the convention used by GDC dataset folders (e.g. {@code 2021-03-10-US-SVL-1}).
	 *
	 * @param path       full file path containing the date folder
	 * @param segFromEnd index from the end of the path split by {@code /}
	 * @return {@code String[]{year, doy}} where {@code doy} is zero-padded to 3 digits
	 */
	private static String[] dateFromFolder(String path, int segFromEnd) {
		String[] parts = path.split("/");
		String[] date  = parts[parts.length - segFromEnd].split("-");
		Calendar c = Calendar.getInstance(TimeZone.getTimeZone("UTC"));
		c.set(Integer.parseInt(date[0]), Integer.parseInt(date[1]) - 1, Integer.parseInt(date[2]), 0, 0, 0);
		return new String[] { String.valueOf(c.get(Calendar.YEAR)),
							  String.format("%03d", c.get(Calendar.DAY_OF_YEAR)) };
	}

	/**
	 * Constructs the MGEX precise-product directory path for a given year and DOY.
	 * All product files for a given day are expected to reside in
	 * {@code baseUrl/YYYY_DOY/}.
	 *
	 * @param baseUrl root directory containing per-day product folders
	 * @param year    four-digit year string (e.g. {@code "2024"})
	 * @param doy     three-digit day-of-year string (e.g. {@code "200"})
	 * @return directory path ending with {@code /}, ready for filename concatenation
	 */
	private static String pDir(String baseUrl, String year, String doy) {
		return baseUrl + year + "_" + doy + "/";
	}
}
