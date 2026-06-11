package com.gnssAug.IGS;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collections;
import java.util.HashMap;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.gnssAug.Android.constants.Constellation;
import com.gnssAug.Rinex.fileParser.Antenna;
import com.gnssAug.Rinex.fileParser.DCB_Bias;
import com.gnssAug.Rinex.fileParser.OSB_Bias;
import com.gnssAug.Rinex.fileParser.Clock;
import com.gnssAug.Rinex.fileParser.Orbit;
import com.gnssAug.Rinex.models.IonoCoeff;
import com.gnssAug.Rinex.models.NavigationMsg;
import com.gnssAug.Rinex.models.Observable;
import com.gnssAug.Rinex.models.ObservationMsg;
import com.gnssAug.Rinex.models.Satellite;
import com.gnssAug.helper.ComputeSatPos;
import com.gnssAug.utility.Closest;
import com.gnssAug.utility.Vector;

/**
 * Pre-processes one RINEX observation epoch for an IGS/geodetic receiver, producing
 * a list of corrected {@link Satellite} objects ready for {@link IGS#filterSat} and
 * the positioning filter.
 *
 * <p><b>Two operating modes, selected by {@code useIGS}:</b>
 * <ul>
 *   <li><b>IGS mode ({@code true}):</b> uses MGEX SP3/CLK precise products. Applies
 *       satellite PCO (via {@link Antenna}), carrier-phase wind-up with eclipse
 *       exclusion, and OSB hardware biases. Returns phase-center-corrected positions.</li>
 *   <li><b>Broadcast mode ({@code false}):</b> uses broadcast navigation messages.
 *       No PCO, wind-up, or bias correction is applied.</li>
 * </ul>
 *
 * <p><b>Eclipse exclusion:</b> when {@link Antenna#getSatPC_windup} returns {@code null}
 * (satellite in Earth's umbra or within the 30-minute post-eclipse recovery window),
 * the satellite is dropped from this epoch and its entry in {@code phase_windup_map}
 * is cleared to prevent stale continuity tracking when the satellite re-emerges.
 */
public class SingleFreq {

	/** Speed of light in vacuum (m/s). */
	final static double SpeedofLight = 299792458;

	/**
	 * Persistent carrier-phase wind-up accumulator, keyed by satellite ID string
	 * (e.g. {@code "G1C5"} = GPS PRN 5, signal G1C).
	 * Wind-up is accumulated across epochs; the previous value is passed to
	 * {@link Antenna#getSatPC_windup} to maintain phase continuity.
	 * Cleared for a satellite on eclipse entry or signal loss.
	 */
	private static HashMap<String, Double> phase_windup_map = new HashMap<String, Double>();

	/**
	 * Processes all observable codes for one RINEX observation epoch.
	 *
	 * <p><b>IGS mode — per-satellite processing steps:</b>
	 * <ol>
	 *   <li>Estimate signal travel time: {@code tSV = tRX − PR/c}.</li>
	 *   <li>Retrieve precise clock offset and drift from {@link Clock}.
	 *       The DCB is already embedded in the clock product by the clock parser.</li>
	 *   <li>Look up OSB satellite code and phase biases from {@link OSB_Bias}
	 *       (converted to metres using the signal wavelength).</li>
	 *   <li>Compute GPS transmission time: {@code t = tSV − satClkOff + codeBias/c}.</li>
	 *   <li>Interpolate satellite ECEF position and velocity at {@code t} from SP3
	 *       via {@link Orbit} (10th-order Lagrange polynomial).</li>
	 *   <li>Apply relativistic clock correction (dot-product term):
	 *       {@code Δt_rel = −2(r·ṙ)/c²}, update {@code satClkOff}, and recompute
	 *       the transmission time.</li>
	 *   <li>Call {@link Antenna#getSatPC_windup} to obtain:
	 *       <ul>
	 *         <li>Indices [0..2]: PCO-corrected satellite ECEF position.</li>
	 *         <li>Index [3]: cumulative phase wind-up in metres.</li>
	 *         <li>{@code null}: satellite is eclipsing — skip and clear wind-up state.</li>
	 *       </ul></li>
	 *   <li>Form corrected observables:
	 *       <ul>
	 *         <li>Pseudorange: {@code PR + c·satClkOff − codeBias}</li>
	 *         <li>Carrier phase: {@code CP + c·satClkOff − phaseBias − windUp[m]}</li>
	 *         <li>Doppler / pseudorange rate: {@code dPR + c·satClkDrift}</li>
	 *       </ul></li>
	 * </ol>
	 *
	 * <p><b>Broadcast mode:</b> delegates entirely to {@link ComputeSatPos#computeSatPos};
	 * no PCO, wind-up, or bias corrections are applied.
	 *
	 * @param obsvMsg      one epoch of RINEX observations from {@link ObservationMsg}
	 * @param NavMsgs      broadcast navigation messages keyed by PRN (broadcast mode only)
	 * @param obsvCodeList signals to process, e.g. {@code {"G1C", "G5Q"}}
	 * @param useIGS       {@code true} = precise SP3/CLK products; {@code false} = broadcast nav
	 * @param useBias      if {@code true}, apply OSB code and phase hardware biases
	 * @param ionoCoeff    Klobuchar coefficients from NAV header (not applied here, passed through)
	 * @param dcb_bias     DCB object (not used directly; embedded in the Clock product)
	 * @param osb__bias    OSB bias object for per-satellite, per-signal biases
	 * @param orbit        SP3 precise orbit for satellite position interpolation
	 * @param clock        CLK precise clock for satellite clock offset and drift
	 * @param antenna      ATX antenna model for PCO correction and phase wind-up
	 * @param tRX          receiver signal reception time in GPS seconds of week
	 * @param weekNo       GPS week number of the current epoch
	 * @param time         epoch as {@link Calendar} in UTC
	 * @param discardSet   satellite IDs to permanently exclude (e.g. {@code "G10"})
	 * @param refEcef      approximate receiver ECEF position for wind-up geometry;
	 *                     {@code null} on the first epoch (those satellites are skipped)
	 * @return list of corrected {@link Satellite} objects for this epoch,
	 *         ready for {@link IGS#filterSat} and the downstream EKF
	 */
	public static ArrayList<Satellite> process(ObservationMsg obsvMsg,
			HashMap<Integer, ArrayList<NavigationMsg>> NavMsgs, String[] obsvCodeList, boolean useIGS, boolean useBias,
			IonoCoeff ionoCoeff, DCB_Bias dcb_bias,OSB_Bias osb__bias, Orbit orbit, Clock clock, Antenna antenna, double tRX, long weekNo,
			Calendar time, Set<String> discardSet, double[] refEcef) throws Exception {
		ArrayList<Satellite> SV = new ArrayList<Satellite>();
		if (obsvCodeList.length > 1 && !useIGS) {
			throw new Exception("Multi-Constellation is not supported without IGS");
		}

		if (useIGS) {
			
			for (String obsvCode : obsvCodeList) {
				ArrayList<Observable> observables = obsvMsg.getObsvSat(obsvCode);
				if (observables == null) {
					continue;
				}
				observables.removeAll(Collections.singleton(null));
				int satCount = observables.size();
				char SSI = obsvCode.charAt(0);

				int polyOrder = 10;
				orbit.findPts(tRX, polyOrder);
				clock.findPts(tRX);
				for (int i = 0; i < satCount; i++) {
					Observable sat = observables.get(i);
					// PRN
					int SVID = sat.getSVID();
					String code = SSI + "" + SVID;
					if (discardSet.contains(code)) {
						continue;
					}
					double tSV = tRX - (sat.getPseudorange() / SpeedofLight);

					double[] sat_ClkOff_Drift = clock.getBiasAndDrift(tSV, SVID, obsvCode, false);
					double satClkOff = sat_ClkOff_Drift[0];
					double satClkDrift = sat_ClkOff_Drift[1];
					int freqBand = Integer.parseInt(obsvCode.charAt(1) + "");
					double wavelength = SpeedofLight / Constellation.frequency.get(SSI).get(freqBand);
					double satHardCodeBiasMeters = osb__bias.getOSBMeters(SSI,"C"+obsvCode.substring(1), SVID, tSV, wavelength);
					double satHardPhaseBiasMeters = osb__bias.getOSBMeters(SSI,"L"+obsvCode.substring(1), SVID, tSV, wavelength);
					
					// GPS System transmission time
					double t = tSV - satClkOff + satHardCodeBiasMeters / SpeedofLight;
					double[][] satPV = orbit.getPV(t, SVID, polyOrder, SSI);
					if (satPV == null) {
						System.err.println(SSI + "" + SVID + " MGEX data absent");
						continue;
					}
					double[] satECEF = satPV[0];
					double[] satVel = satPV[1];

					double relativistic_error = -2 * (Vector.dotProd(satECEF, satVel)) / Math.pow(SpeedofLight, 2);
					// Correct sat clock offset for relativistic error and recompute the Sat coords
					satClkOff += relativistic_error;
					t = tSV - satClkOff + satHardCodeBiasMeters / SpeedofLight;
					String key = obsvCode + SVID;
					double prev_previousWindUpCycles = phase_windup_map.containsKey(key)?phase_windup_map.get(key):0;
					double[] satPC_windup = antenna.getSatPC_windup(SVID, obsvCode, tRX, weekNo, satECEF, refEcef,
							prev_previousWindUpCycles);
					if (satPC_windup == null) {
						phase_windup_map.remove(key); // reset continuity on eclipse exclusion
						continue;
					}
					double phase_windup = satPC_windup[3];
					phase_windup_map.put(key, phase_windup);
					IntStream.range(0, 3).forEach(j -> satECEF[j] = satPC_windup[j]);
					sat.setPseudorange(sat.getPseudorange() + (SpeedofLight * satClkOff) - satHardCodeBiasMeters);
					sat.setPseudoRangeRate(sat.getPseudoRangeRate() + (SpeedofLight * satClkDrift));
					sat.setPhase(sat.getPhase() + (SpeedofLight * satClkOff) - satHardPhaseBiasMeters - phase_windup);
					Satellite _sat = new Satellite(sat, satECEF, satClkOff, t, tRX, satVel, satClkDrift, null, time);
					_sat.compECI();
					SV.add(_sat);

				}
			}
			
		} else {
			ArrayList<Observable> observables = obsvMsg.getObsvSat(obsvCodeList[0]);
			if (observables == null) {
				return SV;
			}
			observables.removeAll(Collections.singleton(null));
			int satCount = observables.size();
			// find out index of nav-msg inside the nav-msg list which is most suitable for
			// each obs-msg based on time
			int order[] = observables.stream().map(i -> NavMsgs.get(i.getSVID()))
					.map(i -> (ArrayList<Double>) i.stream().map(j -> j.getTOC()).collect(Collectors.toList()))
					.mapToInt(i -> Closest.findClosest(tRX, i)).toArray();

			for (int i = 0; i < order.length; i++) {

				Observable sat = observables.get(i);
				// PRN
				int SVID = sat.getSVID();
				String code = sat.getSSI() + "" + SVID;
				if (discardSet.contains(code)) {
					continue;
				}
				// IGS .BSX file DCB
				double ISC = 0;
				NavigationMsg NavMsg = NavMsgs.get(SVID).get(order[i]);
				if (useBias) {
					ISC = dcb_bias.getISC(obsvCodeList[0], SVID);

				}

				double tSV = tRX - (sat.getPseudorange() / SpeedofLight);

				Object[] SatParams = ComputeSatPos.computeSatPos(NavMsg, tSV, tRX, ISC);
				double[] ECEF_SatClkOff = (double[]) SatParams[0];
				double[] SatVel = (double[]) SatParams[1];
				// Note this Clock Drift is derived, it not what we get from Ephemeris
				double SatClkDrift = (double) SatParams[2];
				// GPS System time at time of transmission
				double t = (double) SatParams[3];
				// ECI coordinates
				double[] ECI = (double[]) SatParams[4];
				// AbsoluteDate date = new AbsoluteDate(time.getTime(),
				// TimeScalesFactory.getGPS());
				// double ele = tpf.getElevation(new Vector3D(Arrays.copyOfRange(ECEF_SatClkOff,
				// 0, 3)), frame, date);
				// double az = tpf.getAzimuth(new Vector3D(Arrays.copyOfRange(ECEF_SatClkOff, 0,
				// 3)), frame, date);
				sat.setPseudoRangeRate(sat.getPseudoRangeRate() + (SpeedofLight * SatClkDrift));
				sat.setPseudorange(sat.getPseudorange() + (SpeedofLight * ECEF_SatClkOff[3]));
				SV.add(new Satellite(sat, Arrays.copyOfRange(ECEF_SatClkOff, 0, 3), ECEF_SatClkOff[3], t, tRX, SatVel,
						SatClkDrift, ECI, time));

			}
		}

		return SV;
	}

}
