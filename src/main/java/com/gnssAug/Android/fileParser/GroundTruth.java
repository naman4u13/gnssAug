package com.gnssAug.Android.fileParser;

import java.io.FileReader;
import java.util.ArrayList;

import com.opencsv.CSVReader;

public class GroundTruth {
	private static final long NumberMilliSecondsWeek = 604800000;

	public static ArrayList<double[]> processCSV(String path) throws Exception {

		ArrayList<double[]> data = new ArrayList<double[]>();
		try {
			CSVReader reader = new CSVReader(new FileReader(path));
			String[] line;
			String[] header = reader.readNext();
			// reads one line at a time
			while ((line = reader.readNext()) != null) {
				double millisSinceGpsEpoch = Double.parseDouble(line[2]);
				double GPStime = ((millisSinceGpsEpoch % NumberMilliSecondsWeek) / 1000);
				double weekNo = Math.floor(millisSinceGpsEpoch / NumberMilliSecondsWeek);
				double lat = Double.parseDouble(line[3]);
				double lon = Double.parseDouble(line[4]);
				double alt = Double.parseDouble(line[5]);

				data.add(new double[] { GPStime, weekNo, lat, lon, alt });

			}
			reader.close();

		} catch (Exception e) {
			// TODO: handle exception
			throw new Exception("Error occured during parsing of Ground Truth CSV file \n" + e);

		}
		return data;

	}

}
