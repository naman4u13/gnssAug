package com.gnssAug.Android;

import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collections;
import java.util.HashMap;
import java.util.TimeZone;

import com.gnssAug.Android.models.Derived;
import com.gnssAug.Android.models.GNSSLog;
import com.gnssAug.Android.models.Satellite;

public class SingleFreq {
	private final static double SpeedofLight = 299792458;

	public static ArrayList<Satellite> process(double tRX,
			HashMap<Long, HashMap<String, HashMap<Integer, Derived>>> derivedMap,
			HashMap<String, ArrayList<GNSSLog>> gnssLogMap, Calendar time, String[] obsvCodeList, int weekNo) {

		SimpleDateFormat sdf = new SimpleDateFormat("MM/dd/YYYY kk:mm:ss.SSS");
		String errStr = sdf.format(time.getTime());
		sdf.setTimeZone(TimeZone.getTimeZone("UTC"));
		ArrayList<Satellite> SV = new ArrayList<Satellite>();
		long key = Math.round(tRX * 1e3);
		HashMap<String, HashMap<Integer, Derived>> derObsvMap = null;

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

		for (String obsvCode : obsvCodeList) {

			if (!derObsvMap.containsKey(obsvCode)) {

				System.err.println("No data for obsvCode " + obsvCode + " in derived.csv at time = " + errStr);
				continue;
			}

			HashMap<Integer, Derived> navMap = derObsvMap.get(obsvCode);
			ArrayList<GNSSLog> gnssLog = gnssLogMap.get(obsvCode);

			if (gnssLog == null) {

				System.err.println("No data for obsvCode " + obsvCode + " in Rinex Obs at time = " + errStr);
				continue;
			}
			gnssLog.removeAll(Collections.singleton(null));
			char SSI = obsvCode.charAt(0);

			int satCount = gnssLog.size();

			for (int i = 0; i < satCount; i++) {
				GNSSLog logObs = gnssLog.get(i);
				int svid = logObs.getSvid();
				if (!navMap.containsKey(svid)) {

					System.err.println("No data for svid " + svid + " belonging to obsvCode" + obsvCode
							+ " in derived.csv at time = " + errStr);
					continue;
				}

				Derived navData = navMap.get(svid);

				double[] satECEF = navData.getSatECEF();
				// For non-GPS signal there will be additional biases
				double t = logObs.gettTx() - navData.getSatClkBias();
				// Corrected for satClk Off
				double rawPR = (logObs.gettRx() - t) * SpeedofLight;
				double corrPR = rawPR - navData.getIsrbM() - navData.getIonoDelayM() - navData.getTropoDelayM();
				double corrPRrate = logObs.getPseudorangeRateMetersPerSecond() + navData.getSatClkDrift();
				// NOTE: satClkDrift require investigation
				Satellite sat = new Satellite(logObs, t, corrPR, navData.getSatECEF(), navData.getSatVel(), corrPRrate);
				SV.add(sat);

			}
		}

		return SV;

	}
}
