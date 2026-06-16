package com.gnssAug.Rinex.fileParser;

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.Scanner;
import java.util.Set;
import java.util.stream.IntStream;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.orekit.bodies.CelestialBody;
import org.orekit.bodies.CelestialBodyFactory;
import org.orekit.frames.FramesFactory;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.TimeScalesFactory;

import com.gnssAug.Android.constants.Constellation;
import com.gnssAug.IGS.models.IGSAntenna;
import com.gnssAug.utility.LatLonUtil;
import com.gnssAug.utility.SatUtil;
import com.gnssAug.utility.StringUtil;
import com.gnssAug.utility.Time;
import com.gnssAug.utility.Vector;
import com.opencsv.CSVReader;
import com.opencsv.CSVWriter;

/**
 * Loads IGS satellite antenna PCO (Phase Center Offset) data from a pre-processed CSV
 * (produced by {@link #buildCSV}) and provides two services per epoch per satellite:
 * <ol>
 *   <li>PCO-corrected satellite phase center position in ECEF.</li>
 *   <li>Phase wind-up correction (Wu et al. 1993) using the nominal yaw-steering body frame.</li>
 * </ol>
 *
 * <p>Satellites entering Earth's umbra or undergoing a noon/midnight yaw maneuver are
 * excluded ({@link #getSatPC_windup} returns {@code null}) for the eclipse duration
 * plus a 30-minute recovery window (Kouba 2009).</p>
 */
public class Antenna {

	private HashMap<Character, HashMap<Integer, HashMap<Integer, ArrayList<IGSAntenna>>>> satAntMap;
	// antennaTypeKey (e.g. "TRM115000.00 NONE") → SSI → freq band → [E, N, U] in metres
	private HashMap<String, HashMap<Character, HashMap<Integer, double[]>>> rxAntMap = new HashMap<>();
	private static final double SPEED_OF_LIGHT = 299792458;
	private CelestialBody sun;
	private static final double ECLIPSE_RECOVERY_SECONDS = 1800.0;
	private final HashMap<String, Double> eclipseExitMap = new HashMap<>();
	private final HashMap<String, Boolean> prevInEclipseMap = new HashMap<>();
	private static final java.util.Set<String> _pcoWarnedSet = new java.util.HashSet<>();

	/**
	 * @param path  Path to the pre-processed satellite antenna CSV (produced by {@link #buildCSV}).
	 *              Orekit's {@code DataProvidersManager} must be initialised before constructing
	 *              this object — call {@code MainApp.buildGeoid()} first.
	 */
	public Antenna(String path) throws Exception {
		satAntMap = new HashMap<Character, HashMap<Integer, HashMap<Integer, ArrayList<IGSAntenna>>>>();
		readCSV(path);
		sun = CelestialBodyFactory.getSun();
	}

	/**
	 * Populates {@code satAntMap} from the pre-processed antenna CSV.
	 *
	 * <p>The CSV is produced by {@link #buildCSV} from an IGS ANTEX file. Each row covers
	 * one (satellite, frequency) pair with columns: TYPE, SVID, DAZI, ZEN, VALID_FROM,
	 * VALID_UNTIL, FREQUENCY, NEU/XYZ offsets, PCV_NOAZI, PCV_AZI.</p>
	 *
	 * <p>Only rows with TYPE = {@code 'S'} (satellite) are loaded. The CSV is ordered so
	 * that receiver antenna rows ({@code 'R'}) follow all satellite rows; hitting the first
	 * non-{@code 'S'} row is therefore used as an early-exit sentinel.</p>
	 *
	 * <p>Entries are indexed as: SSI (constellation char) → PRN → frequency band → list of
	 * {@link IGSAntenna}, where the list holds multiple validity-period entries for the same
	 * (satellite, frequency) pair. {@link #getSatPC_windup} walks the list in reverse to
	 * find the most recent valid entry for a given epoch.</p>
	 */
	private void readCSV(String path) throws Exception {
		try {
			CSVReader reader = new CSVReader(new FileReader(path));
			String[] line;
			reader.readNext(); // skip header row
			while ((line = reader.readNext()) != null) {

				char antType = line[0].charAt(0);
				if (antType == 'R') {
					// Receiver PCO row: antennaTypeKey → SSI → freq → [E, N, U] in metres
					String antennaKey = line[1].replaceAll("\\s+", " ").trim();
					char ssi = line[6].charAt(0);
					int freq = Integer.parseInt(line[6].substring(1));
					double[] neu_mm = StringUtil.str2arr_D(line[7]); // ANTEX order: [N, E, U] in mm
					double[] enu_m = { neu_mm[1] / 1000.0, neu_mm[0] / 1000.0, neu_mm[2] / 1000.0 };
					rxAntMap.computeIfAbsent(antennaKey, k -> new HashMap<>())
							.computeIfAbsent(ssi, k -> new HashMap<>())
							.put(freq, enu_m);
					continue;
				} else if (antType != 'S') {
					continue;
				}
				String _SVID = line[1];
				char SSI = _SVID.charAt(0);
				int SVID = Integer.parseInt(_SVID.substring(1));
				String _freq = line[6];
				int freq = Integer.parseInt(_freq.substring(1));
				IGSAntenna satAnt = new IGSAntenna(line[2], line[3], line[4], line[5], line[7], line[8], line[9]);
				satAntMap.computeIfAbsent(SSI, k -> new HashMap<Integer, HashMap<Integer, ArrayList<IGSAntenna>>>())
						.computeIfAbsent(SVID, k -> new HashMap<Integer, ArrayList<IGSAntenna>>())
						.computeIfAbsent(freq, k -> new ArrayList<IGSAntenna>()).add(satAnt);

			}
			reader.close();

		} catch (Exception e) {
			e.printStackTrace();
			throw new Exception("Error occured during reading and parsing of Antenna(.csv) file \n" + e);
		}
	}

	/**
	 * Returns the receiver Phase Center Offset as [E, N, U] in metres for the given
	 * antenna type, satellite system, and frequency band, loaded from the ANTEX CSV.
	 *
	 * @param antennaTypeKey  Normalised antenna type string, e.g. {@code "TRM115000.00 NONE"}.
	 *                        Must match the normalised form stored during {@link #readCSV}.
	 * @param ssi             Satellite system identifier char, e.g. {@code 'G'} for GPS.
	 * @param freq            Frequency band number, e.g. {@code 1} for L1, {@code 5} for L5.
	 * @return [E, N, U] in metres, or {@code null} if not found in the CSV.
	 */
	public double[] getRxPCO_ENU(String antennaTypeKey, char ssi, int freq) {
		HashMap<Character, HashMap<Integer, double[]>> ssiMap = rxAntMap.get(antennaTypeKey);
		if (ssiMap == null) return null;
		HashMap<Integer, double[]> freqMap = ssiMap.get(ssi);
		if (freqMap == null) return null;
		return freqMap.get(freq);
	}

	/**
	 * Computes the PCO-corrected satellite phase center position and phase wind-up
	 * correction for a single frequency observation (Wu et al. 1993).
	 *
	 * <p>The satellite body frame is constructed under the nominal yaw-steering model:
	 * Z-axis toward Earth (nadir), Y-axis along the solar panel axis (nadir × toSun),
	 * X-axis completing the right-hand frame. The PCO vector from the ANTEX CSV is
	 * rotated into ECEF using this frame and added to the antenna mechanical center.</p>
	 *
	 * <p>Satellites in Earth's umbra or at a noon/midnight yaw singularity (magJ &lt; 1e-6)
	 * are excluded: returns {@code null} for the eclipse duration plus a
	 * {@value #ECLIPSE_RECOVERY_SECONDS}-second recovery window (Kouba 2009).
	 * Callers must remove the satellite's wind-up history entry on {@code null} return
	 * to avoid stale continuity tracking.</p>
	 *
	 * @param SVID           Satellite PRN number.
	 * @param obsvCode       Observation code, e.g. {@code "G1C"}.
	 * @param GPSTime        GPS seconds-of-week at signal transmission.
	 * @param weekNo         GPS week number.
	 * @param satMC          Satellite ECEF position (metres), mechanical phase center.
	 * @param userECEF       Receiver ECEF position (metres). Pass {@code null} to skip
	 *                       wind-up (result[3] = 0).
	 * @param previousWindUp Wind-up value [metres] stored from the previous tracked epoch,
	 *                       used for phase-continuity unwrapping. Pass 0 on first use.
	 * @return {@code double[4]}: [0..2] PCO-corrected satellite APC in ECEF (metres),
	 *         [3] phase wind-up in metres; or {@code null} if eclipsed/excluded.
	 */
	public double[] getSatPC_windup(int SVID, String obsvCode, double GPSTime, long weekNo, double[] satMC,
			double[] userECEF, double previousWindUp) {

		double[] eccXYZ = new double[3];
		double[] satPC_windUp = new double[4];

		char SSI = obsvCode.charAt(0);
		int freq = Integer.parseInt(obsvCode.charAt(1) + "");
		double wavelength = SPEED_OF_LIGHT / Constellation.frequency.get(SSI).get(freq);

		HashMap<Integer, ArrayList<IGSAntenna>> svidAntMap = satAntMap.get(SSI).get(SVID);
		ArrayList<IGSAntenna> satAntList = svidAntMap != null ? svidAntMap.get(freq) : null;

		if (satAntList == null) {
			// No PCO data: return uncorrected position so the satellite is still usable
			System.arraycopy(satMC, 0, satPC_windUp, 0, 3);
			String _warnKey = SSI + SVID + "_" + freq;
			if (_pcoWarnedSet.add(_warnKey)) {
				System.err.println("Sat " + SVID + " PCO info unavailable for frequency - " + freq
						+ " ! (not calibrated in ATX — using uncorrected satellite position)");
			}
			return satPC_windUp;
		}

		// Walk backwards to find the most recent valid antenna entry for this epoch
		boolean foundValid = false;
		int _n = satAntList.size();
		for (int j = _n - 1; j >= 0; j--) {
			IGSAntenna satAnt = satAntList.get(j);
			if (satAnt.checkValidity(new double[] { GPSTime, weekNo })) {
				eccXYZ = satAnt.getEccXYZ();
				foundValid = true;
				break;
			}
		}
		if (!foundValid) {
			// eccXYZ = zero: PCO correction silently dropped rather than crashing
			System.err.println("No valid antenna for Sat " + SVID + " at time " + GPSTime);
		}

		// Build nominal yaw-steering body frame
		double R_sat = 0.0;
		for (int i = 0; i < 3; i++) R_sat += satMC[i] * satMC[i];
		R_sat = Math.sqrt(R_sat);

		// Body Z-axis: nadir (toward Earth center)
		double[] unitNadir = new double[3];
		for (int i = 0; i < 3; i++) unitNadir[i] = -satMC[i] / R_sat;

		// Sun position in EME2000 — same frame approximation used throughout this class
		double[] sunXYZ = getSunCoord(GPSTime, weekNo);
		double[] toSun = new double[3];
		double R_sun_sat = 0.0;
		for (int i = 0; i < 3; i++) { toSun[i] = sunXYZ[i] - satMC[i]; R_sun_sat += toSun[i] * toSun[i]; }
		R_sun_sat = Math.sqrt(R_sun_sat);
		double[] unitToSun = new double[3];
		for (int i = 0; i < 3; i++) unitToSun[i] = toSun[i] / R_sun_sat;

		// Body Y-axis: solar panel axis = nadir × toSun
		double[] j = Vector.crossProd(unitNadir, unitToSun);
		double magJ = Vector.mod(j);

		// Eclipse / yaw-singularity exclusion (Kouba 2009)
		String satKey = obsvCode.charAt(0) + "" + SVID;
		boolean inUmbra = isInUmbra(satMC, sunXYZ);
		if (inUmbra || magJ < 1e-6) {
			if (!Boolean.TRUE.equals(prevInEclipseMap.get(satKey))) {
				System.err.println("Eclipse entry: " + satKey + " at GPS time " + GPSTime);
			}
			prevInEclipseMap.put(satKey, true);
			return null;
		}
		if (Boolean.TRUE.equals(prevInEclipseMap.get(satKey))) {
			eclipseExitMap.put(satKey, GPSTime);
			prevInEclipseMap.put(satKey, false);
			System.err.println("Eclipse exit: " + satKey + ", recovery started at GPS time " + GPSTime);
		}
		Double exitTime = eclipseExitMap.get(satKey);
		if (exitTime != null && (GPSTime - exitTime) < ECLIPSE_RECOVERY_SECONDS) {
			return null;
		}

		for (int i = 0; i < 3; i++) j[i] /= magJ;

		// Body X-axis: cross(j, nadir) — already unit since j ⊥ unitNadir
		double[] i = Vector.crossProd(j, unitNadir);

		// Apply PCO: rotate ANTEX offset (X, Y, Z_body) into ECEF and add to satMC
		for (int y = 0; y < 3; y++) {
			satPC_windUp[y] = satMC[y] + eccXYZ[0] * i[y] + eccXYZ[1] * j[y] + eccXYZ[2] * unitNadir[y];
		}

		if (userECEF != null) {
			double[] userLL = LatLonUtil.ecef2lla(userECEF);
			double[] unitLOS = SatUtil.getUnitLOS(Arrays.copyOfRange(satPC_windUp, 0, 3), userECEF);

			// Receiver effective dipole vector Dr (fixed leveled antenna: East/North basis)
			double[] ar = new double[] { -Math.sin(userLL[1]), Math.cos(userLL[1]), 0.0 };           // East
			double[] br = new double[] { -Math.sin(userLL[0]) * Math.cos(userLL[1]),
					-Math.sin(userLL[0]) * Math.sin(userLL[1]), Math.cos(userLL[0]) };                // North
			double[] Dr = Vector.add(Vector.subtract(ar, Vector.scale(unitLOS, Vector.dotProd(unitLOS, ar))),
					Vector.crossProd(unitLOS, br));

			// Satellite effective dipole vector Ds (body X-axis; minus sign per Wu et al. 1993)
			double[] as = i;
			double[] bs = j;
			double[] Ds = Vector.add(Vector.subtract(as, Vector.scale(unitLOS, Vector.dotProd(unitLOS, as))),
					Vector.scale(Vector.crossProd(unitLOS, bs), -1.0));

			double modDr = Vector.mod(Dr);
			double modDs = Vector.mod(Ds);
			if (modDr == 0 || modDs == 0) {
				// Degenerate geometry (satellite at zenith): wind-up undefined, set to 0
				System.err.println("Invalid dipole magnitude for wind-up computation.");
				satPC_windUp[3] = 0.0;
				return satPC_windUp;
			}

			double cosTheta = Math.max(-1.0, Math.min(1.0, Vector.dotProd(Ds, Dr) / (modDs * modDr)));
			double sign = Math.signum(Vector.dotProd(unitLOS, Vector.crossProd(Ds, Dr)));
			double deltaPhi = sign * Math.acos(cosTheta); // raw angle in [-π, π]

			// Unwrap to nearest 2π multiple of the previous epoch's value to ensure continuity
			double twoPi = 2 * Math.PI;
			double previousRadians = (previousWindUp / wavelength) * twoPi;
			double diff = deltaPhi - previousRadians;
			deltaPhi -= Math.round(diff / twoPi) * twoPi;

			satPC_windUp[3] = (deltaPhi / twoPi) * wavelength;
		}
		return satPC_windUp;
	}

	/**
	 * Returns true if the satellite is inside Earth's umbra (full shadow).
	 *
	 * <p>Uses a cylindrical shadow model (point-source Sun assumption): the satellite is in
	 * umbra when it is on the anti-Sun side of Earth <em>and</em> its perpendicular distance
	 * from the Earth-Sun axis is less than Earth's radius (6378136 m). This slightly
	 * over-estimates the true umbra cone edge by a few seconds of arc, which is conservative
	 * and acceptable for exclusion purposes. Penumbra is not detected.</p>
	 *
	 * <p>Note: {@code satECEF} and {@code sunECEF} are in different frames (ECEF vs EME2000),
	 * but because Earth is at the origin of both and the Sun is ~1 AU away, the directional
	 * error from this frame mismatch is negligible (~0.003°).</p>
	 *
	 * @param satECEF  Satellite position in ECEF (metres).
	 * @param sunECEF  Sun position in EME2000 from {@link #getSunCoord} (metres).
	 * @return {@code true} if the satellite is within Earth's shadow cylinder.
	 */
	private boolean isInUmbra(double[] satECEF, double[] sunECEF) {
		double sunDist = Math.sqrt(sunECEF[0]*sunECEF[0] + sunECEF[1]*sunECEF[1] + sunECEF[2]*sunECEF[2]);
		double[] uSun = new double[]{ sunECEF[0]/sunDist, sunECEF[1]/sunDist, sunECEF[2]/sunDist };
		double proj = satECEF[0]*uSun[0] + satECEF[1]*uSun[1] + satECEF[2]*uSun[2];
		if (proj >= 0) return false; // satellite on Sun-side of Earth, cannot be in shadow
		double satR2 = satECEF[0]*satECEF[0] + satECEF[1]*satECEF[1] + satECEF[2]*satECEF[2];
		double perpDist2 = satR2 - proj * proj;
		return perpDist2 < 6378136.0 * 6378136.0;
	}

	/**
	 * Returns the Sun's position in the EME2000 inertial frame (metres) at the given epoch,
	 * queried via Orekit's {@link CelestialBody} interface.
	 *
	 * <p>GPS time is converted to TAI by adding 19 seconds (the fixed GPS–TAI offset)
	 * before constructing the {@link AbsoluteDate}.</p>
	 *
	 * @param GPSTime  GPS seconds-of-week at the epoch of interest.
	 * @param weekNo   GPS week number.
	 * @return Sun position vector [x, y, z] in EME2000 (metres).
	 */
	private double[] getSunCoord(double GPSTime, long weekNo) {
		Date date = Time.getDate(GPSTime + 19, weekNo, 0).getTime();
		AbsoluteDate absDate = new AbsoluteDate(date, TimeScalesFactory.getTAI());
		Vector3D coords = sun.getPVCoordinates(absDate, FramesFactory.getEME2000()).getPosition();
		return new double[] { coords.getX(), coords.getY(), coords.getZ() };
	}

	/**
	 * One-time offline utility: converts an IGS ANTEX (.atx) file into the flat CSV format
	 * consumed by {@link #Antenna(String)}. Re-run whenever a new ANTEX file is available.
	 *
	 * @param in_path   Path to the IGS ANTEX (.atx) file.
	 * @param out_path  Destination path for the output CSV.
	 */
	public static void buildCSV(String in_path, String out_path) throws Exception {
		try {
			Set<Character> SSIset = Set.of('G', 'R', 'E', 'C', 'I', 'S', 'J');

			String[] header = { "TYPE", "SVID/SrNo", "DAZI", "ZEN", "VALID_FROM", "VALID_UNTIL", "FREQUENCY", "NEU/XYZ",
					"PCV_NOAZI", "PCV_AZI" };
			FileWriter outputfile = new FileWriter(new File(out_path));
			CSVWriter writer = new CSVWriter(outputfile);
			writer.writeNext(header);

			File file = new File(in_path);
			Scanner input = new Scanner(file);
			input.useDelimiter("END OF HEADER");
			input.next();
			input.useDelimiter("START OF ANTENNA|END OF ANTENNA");
			input.next();
			while (input.hasNext()) {
				String str = input.next().trim();
				if (str.isBlank()) {
					continue;
				}
				String[] lines = str.split("\n");

				String[] type_sno = StringUtil.splitter(lines[0], false, 20, 20, 10, 10);
				char type = 'R';
				String typeCode = type_sno[1];
				if (typeCode.length() == 3 && SSIset.contains(typeCode.charAt(0))) {
					type = 'S';
				} else {
					typeCode = type_sno[0];
					// break;
				}

				double dazi = Double.parseDouble(StringUtil.splitter(lines[2], false, 2, 6, 52)[1]);
				String[] _zen = StringUtil.splitter(lines[3], false, 2, 6, 6, 6, 40);
				double[] zen = IntStream.range(1, 4).mapToDouble(x -> Double.parseDouble(_zen[x])).toArray();
				int col = (int) ((zen[1] - zen[0]) / zen[2]) + 1;

				int freqCount = Integer.parseInt(StringUtil.splitter(lines[4], false, 6, 54)[0]);

				String[] temp = StringUtil.splitter(lines[5], 60, 20);
				double[] validFrom = null;
				int index = 6;
				if (temp[1].equalsIgnoreCase("VALID FROM")) {

					String[] _validFrom = temp[0].split("\\s+");
					validFrom = Time.getGPSTime(_validFrom);
					index++;
				}

				temp = StringUtil.splitter(lines[6], 60, 20);
				double[] validUntil = null;

				if (temp[1].equalsIgnoreCase("VALID UNTIL")) {

					String[] _validUntil = temp[0].split("\\s+");
					validUntil = Time.getGPSTime(_validUntil);
					index++;
				}
				while (true) {
					temp = StringUtil.splitter(lines[index], 60, 20);
					if (!temp[1].equalsIgnoreCase("COMMENT")) {
						break;
					}
					index++;

				}
				while (freqCount > 0) {
					String freq = StringUtil.splitter(lines[index++], false, 3, 3, 54)[1];

					String[] _NEU_XYZ = StringUtil.splitter(lines[index++], false, 10, 10, 10, 30);
					double[] NEU_XYZ = IntStream.range(0, 3).mapToDouble(x -> Double.parseDouble(_NEU_XYZ[x]))
							.toArray();

					String[] NOAZI = lines[index++].trim().split("\\s+");
					double[] PCVnoazi = IntStream.range(1, NOAZI.length).mapToDouble(x -> Double.parseDouble(NOAZI[x]))
							.toArray();
					double[][] PCVazi = null;
					if (dazi > 0) {
						int row = (int) (360 / dazi) + 1;
						PCVazi = new double[row][col];
						for (int j = 0; j < row; j++) {

							String[] strArr = lines[index++].trim().split("\\s+");

							PCVazi[j] = IntStream.range(1, col + 1)
									.mapToDouble(x -> Double.parseDouble(strArr[x].trim())).toArray();

						}
					}
					String[] line = new String[] { type + "", typeCode, dazi + "", Arrays.toString(zen),
							Arrays.toString(validFrom), Arrays.toString(validUntil), freq, Arrays.toString(NEU_XYZ),
							Arrays.toString(PCVnoazi), Arrays.deepToString(PCVazi) };
					writer.writeNext(line);
					index++;
					freqCount--;

				}

			}
			writer.close();

		} catch (Exception e) {
			e.printStackTrace();
			throw new Exception(
					"Error occured during parsing of Antenna(.atx) file or while generating CSV table \n" + e);

		}

	}

}
