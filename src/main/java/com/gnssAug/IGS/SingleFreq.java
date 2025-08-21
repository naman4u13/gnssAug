package com.gnssAug.IGS;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collections;
import java.util.HashMap;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

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

public class SingleFreq {
	final static double SpeedofLight = 299792458;
	private static HashMap<String, Double> phase_windup_map = new HashMap<String, Double>();

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
					double satHardCodeBias = osb__bias.getOSB(SSI,"C"+obsvCode.substring(1), SVID, tSV);
					double satHardPhaseBias = osb__bias.getOSB(SSI,"L"+obsvCode.substring(1), SVID, tSV);
					
					// GPS System transmission time
					double t = tSV - (satClkOff-satHardCodeBias);
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
					t = tSV - (satClkOff-satHardCodeBias);
					String key = obsvCode + SVID;
					double prev_previousWindUpCycles = phase_windup_map.containsKey(key)?phase_windup_map.get(key):0;
					double[] satPC_windup = antenna.getSatPC_windup_new(SVID, obsvCode, tRX, weekNo, satECEF, refEcef,
							prev_previousWindUpCycles);
					double phase_windup = satPC_windup[3];
					phase_windup_map.put(key, phase_windup);
					IntStream.range(0, 3).forEach(j -> satECEF[j] = satPC_windup[j]);
					// fractional Wind up in cycles, will require further processing to correct for
					// full cycles and then multiply by wavelength
					/* double windup = satPC_windup[3]; */
					sat.setPseudorange(sat.getPseudorange() + (SpeedofLight * (satClkOff-satHardCodeBias)));
					sat.setPseudoRangeRate(sat.getPseudoRangeRate() + (SpeedofLight * satClkDrift));
					sat.setPhase(sat.getPhase() + (SpeedofLight * (satClkOff-satHardPhaseBias))-phase_windup);
					Satellite _sat = new Satellite(sat, satECEF, satClkOff, t, tRX, satVel, satClkDrift, null, time);
					_sat.compECI();
					/* _sat.setPhaseWindUp(windup); */
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
