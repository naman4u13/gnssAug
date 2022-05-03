package com.gnssAug.Android;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.HashMap;
import java.util.TimeZone;
import java.util.TreeMap;

import org.orekit.data.DataProvidersManager;
import org.orekit.data.DirectoryCrawler;
import org.orekit.forces.gravity.potential.GravityFieldFactory;
import org.orekit.forces.gravity.potential.NormalizedSphericalHarmonicsProvider;
import org.orekit.frames.Frame;
import org.orekit.frames.FramesFactory;
import org.orekit.frames.ITRFVersion;
import org.orekit.models.earth.Geoid;
import org.orekit.models.earth.ReferenceEllipsoid;
import org.orekit.utils.IERSConventions;

import com.gnssAug.Android.fileParser.DerivedCSV;
import com.gnssAug.Android.fileParser.GNSSLog;
import com.gnssAug.Android.fileParser.GroundTruth;
import com.gnssAug.Android.models.AndroidGNSSLog;
import com.gnssAug.Android.models.AndroidSatellite;
import com.gnssAug.Android.models.Derived;

public class Android {
	public static void posEstimate(boolean doPosErrPlot, double cutOffAng, int estimatorType, String[] obsvCodeList,
			String derived_csv_path, String gnss_log_path, String GTcsv) {
		try {
			TimeZone.setDefault(TimeZone.getTimeZone("UTC"));
			HashMap<String, ArrayList<HashMap<String, Double>>> ErrMap = new HashMap<String, ArrayList<HashMap<String, Double>>>();

			ArrayList<Calendar> timeList = new ArrayList<Calendar>();
			ArrayList<double[]> trueLLHlist = new ArrayList<double[]>();
			ArrayList<ArrayList<AndroidSatellite>> SVlist = new ArrayList<ArrayList<AndroidSatellite>>();
			double[] trueUserLLH = null;

			String path = "C:\\Users\\Naman\\Desktop\\rinex_parse_files\\google2\\test2";
			File output = new File(path + ".txt");
			PrintStream stream;
			stream = new PrintStream(output);
			System.setOut(stream);

			ArrayList<double[]> rxLLH = GroundTruth.processCSV(GTcsv);
			HashMap<Long, HashMap<String, HashMap<Integer, Derived>>> derivedMap = null;
			TreeMap<Long, HashMap<String, ArrayList<AndroidGNSSLog>>> gnssLogMaps = null;
			derivedMap = DerivedCSV.processCSV(derived_csv_path);
			gnssLogMaps = GNSSLog.process(gnss_log_path);

			int gtIndex = 0;
			for (long tRxMilli : gnssLogMaps.keySet()) {
				HashMap<String, ArrayList<AndroidGNSSLog>> gnssLogMap = gnssLogMaps.get(tRxMilli);
				AndroidGNSSLog entry = ((ArrayList<AndroidGNSSLog>) gnssLogMap.values().toArray()[0]).get(0);
				double tRx = entry.gettRx();
				int weekNo = entry.getWeekNo();
				if (tRxMilli != (rxLLH.get(gtIndex)[0] * 1000) || weekNo != rxLLH.get(gtIndex)[1]) {

					System.err.println("FATAL ERROR - GT timestamp does not match");

					continue;

				}
				trueUserLLH = new double[] { rxLLH.get(gtIndex)[2], rxLLH.get(gtIndex)[3], rxLLH.get(gtIndex)[4] };
				gtIndex++;

//				Calendar time = Time.getDate(tRx, weekNo, 0);
//				ArrayList<AndroidSatellite> SV = SingleFreq.process(tRx, derivedMap, gnssLogMap, time, obsvCodeList,
//						weekNo);

			}

		} catch (Exception e) {
			// TODO: handle exception
			System.out.println(e.getMessage());
			e.printStackTrace();
		}
	}

	public static Geoid buildGeoid() {
		// Semi-major axis or Equatorial radius
		final double ae = 6378137;
		// flattening
		final double f = 1 / 298.257223563;

		// Earth's rotation rate
		final double spin = 7.2921151467E-5;
		// Earth's universal gravitational parameter
		final double GM = 3.986004418E14;

		File orekitData = new File(
				"C:\\Users\\Naman\\Desktop\\rinex_parse_files\\orekit\\orekit-data-master\\orekit-data-master");
		DataProvidersManager manager = DataProvidersManager.getInstance();
		manager.addProvider(new DirectoryCrawler(orekitData));
		NormalizedSphericalHarmonicsProvider nhsp = GravityFieldFactory.getNormalizedProvider(50, 50);
		Frame frame = FramesFactory.getITRF(ITRFVersion.ITRF_2014, IERSConventions.IERS_2010, true);

		// ReferenceEllipsoid refElp = new ReferenceEllipsoid(ae, f, frame, GM, spin);
		Geoid geoid = new Geoid(nhsp, ReferenceEllipsoid.getWgs84(frame));
		return geoid;

	}

}
