package com.gnssAug.Android;

import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collections;
import java.util.HashMap;
import java.util.Set;
import java.util.TimeZone;
import java.util.stream.IntStream;

import com.gnssAug.Android.models.Derived;
import com.gnssAug.Android.models.GNSSLog;
import com.gnssAug.Android.models.Satellite;
import com.gnssAug.Rinex.fileParser.Antenna;
import com.gnssAug.Rinex.fileParser.Clock;
import com.gnssAug.Rinex.fileParser.OSB_Bias;
import com.gnssAug.Rinex.fileParser.Orbit;
import com.gnssAug.utility.MathUtil;
import com.gnssAug.utility.Vector;

public class SingleFreq {
	private final static double SpeedofLight = 299792458;
	private static HashMap<String, Double> phase_windup_map = new HashMap<String, Double>();
	public static ArrayList<Satellite> process(double tRX,
			HashMap<Long, HashMap<String, HashMap<Integer, Derived>>> derivedMap,OSB_Bias osb__bias,Antenna antenna,
			HashMap<String, ArrayList<GNSSLog>> gnssLogMap, Calendar time, String[] obsvCodeList, int weekNo,
			Clock clock, Orbit orbit, boolean useIGS,Set<String> discardSet, double[] refEcef) throws Exception {

		SimpleDateFormat sdf = new SimpleDateFormat("MM/dd/YYYY kk:mm:ss.SSS");
		String errStr = sdf.format(time.getTime());
		sdf.setTimeZone(TimeZone.getTimeZone("UTC"));
		ArrayList<Satellite> SV = new ArrayList<Satellite>();
		long key = Math.round(tRX * 1e3);
		HashMap<String, HashMap<Integer, Derived>> derObsvMap = null;
		if (!useIGS) {
			if (derivedMap.containsKey(key)) {
				derObsvMap = derivedMap.get(key);
			} else if (derivedMap.containsKey(key + 1)) {
				derObsvMap = derivedMap.get(key + 1);
			} else if (derivedMap.containsKey(key - 1)) {
				derObsvMap = derivedMap.get(key - 1);
			} else {

				ArrayList<GNSSLog> errLog = null;
				for (int i = 0; i < obsvCodeList.length; i++) {
					if (gnssLogMap.get(obsvCodeList[i]) != null) {
						errLog = gnssLogMap.get(obsvCodeList[i]);
						break;
					}
				}

				if (errLog == null) {
					System.err.println("No PR info available or captured = " + errStr);
					return SV;
				}
				errLog.removeAll(Collections.singleton(null));
				if (errLog.size() == 0) {
					System.err.println("No PR info available or captured = " + errStr);
					return SV;
				}

				System.err.println("Missing data in derived.csv at time = " + errStr);
				return SV;
			}
		}

		for (String obsvCode : obsvCodeList) {
			HashMap<Integer, Derived> navMap = null;
			if (!useIGS) {
				if (!derObsvMap.containsKey(obsvCode)) {

					System.err.println("No data for obsvCode " + obsvCode + " in derived.csv at time = " + errStr);
					continue;

				}

				navMap = derObsvMap.get(obsvCode);
			}
			ArrayList<GNSSLog> gnssLog = gnssLogMap.get(obsvCode);

			if (gnssLog == null) {

				System.err.println("No data for obsvCode " + obsvCode + " in Rinex Obs at time = " + errStr);
				continue;
			}
			gnssLog.removeAll(Collections.singleton(null));
			char SSI = obsvCode.charAt(0);

			int satCount = gnssLog.size();
			int polyOrder = 10;
			if (useIGS) {

				orbit.findPts(tRX, polyOrder);
				clock.findPts(tRX);
			}
			for (int i = 0; i < satCount; i++) {

				GNSSLog logObs = gnssLog.get(i);
				int svid = logObs.getSvid();
				String code = logObs.getObsvCode()+""+svid;
				
				if(discardSet.contains(code))
				{
					continue;
				}
				double t = 0;

				Derived navData = null;
				double[] satECEF = null;
				double[] satVel = null;
				double corrPR = 0;
				double corrPRrate = 0;
				double corrPhase = 0;
				String state = Integer.toBinaryString(logObs.getState());

				if (useIGS) {
					
					double tSV = logObs.gettTx();
					double[] sat_ClkOff_Drift = clock.getBiasAndDrift(tSV, svid, obsvCode, false);
					double satClkOff = sat_ClkOff_Drift[0];
					double satClkDrift = sat_ClkOff_Drift[1];
					double satHardCodeBias = osb__bias.getOSB(SSI,"C"+obsvCode.substring(1), svid, tSV);
					double satHardPhaseBias = osb__bias.getOSB(SSI,"L"+obsvCode.substring(1), svid, tSV);
					// GPS System transmission time
					t = tSV - (satClkOff-satHardCodeBias);
					double[][] satPV = orbit.getPV(t, svid, polyOrder, SSI);
					if (satPV == null) {
						System.err.println(SSI + "" + svid + " MGEX data absent");
						continue;
					}
					
					satECEF = satPV[0];
					satVel = satPV[1];
					double relativistic_error = -2 * (Vector.dotProd(satECEF, satVel)) / Math.pow(SpeedofLight, 2);
					// Correct sat clock offset for relativistic error and recompute the Sat coords
					satClkOff += relativistic_error;
					t = tSV - (satClkOff-satHardCodeBias);
					String satKey = obsvCode + svid;
					double prev_previousWindUpCycles = phase_windup_map.containsKey(satKey)?phase_windup_map.get(satKey):0;
					double[] satPC_windup = antenna.getSatPC_windup_new(svid, obsvCode, tRX, weekNo, satECEF, refEcef,
							prev_previousWindUpCycles);
					double phase_windup = satPC_windup[3];
					phase_windup_map.put(satKey, phase_windup);
					for(int j=0;j<3;j++)
					{
						satECEF[j] = satPC_windup[j];
					}
				
					double rawPR = (logObs.gettRx() - t) * SpeedofLight;
					corrPR = rawPR-(SpeedofLight*satHardCodeBias);
					corrPRrate = logObs.getPseudorangeRateMetersPerSecond() + (SpeedofLight * satClkDrift);
					corrPhase = logObs.getAccumulatedDeltaRangeMeters()+(SpeedofLight*(satClkOff-satHardPhaseBias))-phase_windup;

				} else {
					if (!navMap.containsKey(svid)) {

						System.err.println("No data for svid " + svid + " belonging to obsvCode" + obsvCode
								+ " in derived.csv at time = " + errStr);
						continue;
					}
					navData = navMap.get(svid);
					satECEF = navData.getSatECEF();
					satVel = navData.getSatVel();
					// For non-GPS signal there will be additional biases
					t = (navData.gettSV() / 1e9) - navData.getSatClkBias();
					double rawPR = navData.getRawPrM() + (navData.getSatClkBias() * SpeedofLight);
					corrPR = rawPR - navData.getIsrbM() - navData.getIonoDelayM() - navData.getTropoDelayM();
					corrPRrate = logObs.getPseudorangeRateMetersPerSecond() + (SpeedofLight * navData.getSatClkDrift());
					corrPhase = logObs.getAccumulatedDeltaRangeMeters() - navData.getIsrbM() - navData.getIonoDelayM() - navData.getTropoDelayM()+ (navData.getSatClkBias() * SpeedofLight);
				}

				// double diff = MathUtil.getEuclidean(satVel, navData.getSatVel());
				// NOTE: satClkDrift require investigation
				Satellite sat = new Satellite(logObs, t, corrPR, satECEF, satVel, corrPRrate,corrPhase);
				SV.add(sat);

			}
		}

		return SV;

	}
}
