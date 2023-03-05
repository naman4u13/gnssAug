package com.gnssAug.utility;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Set;
import java.util.TreeMap;
import java.util.stream.IntStream;

import com.opencsv.CSVWriter;

public class Trajectory {

	public static void createCSV(TreeMap<String, ArrayList<double[]>> trajectoryPosMap,
			TreeMap<String, ArrayList<double[]>> trajectoryVelMap, String path, int n) throws IOException {
		Set<String> pos = trajectoryPosMap.keySet();
		Set<String> vel = trajectoryVelMap.keySet();
		String[] header1 = new String[pos.size() + 1];
		String[] header2 = new String[vel.size() + 1];
		header1[0] = "Position";
		int i = 1;
		for (String key : pos) {
			header1[i] = key;
			i++;
		}
		i = 1;
		header2[0] = "Velocity";
		for (String key : vel) {
			header2[i] = key;
			i++;
		}
		CSVWriter writer = null;
		File file = new File(path + "trajectory.csv");
		// create FileWriter object with file as parameter
		FileWriter outputfile = new FileWriter(file);
		// create CSVWriter object filewriter object as parameter
		writer = new CSVWriter(outputfile);
		writer.writeNext(header1);
		writer.writeNext(header2);
		for (int j = 0; j < n; j++) {
			ArrayList<String> line = new ArrayList<String>();
			for (String key : pos) {

				double[] p = trajectoryPosMap.get(key).get(j);
				if ((p[0] != -999) && (p[1] != -999) && (p[2] != -999)) {
					line.add(p[0] + "");
					line.add(p[1] + "");
					line.add(p[2] + "");

				} else {
					IntStream.range(0, 3).forEach(k -> line.add(""));
				}

			}
			for (String key : vel) {

				double[] v = trajectoryVelMap.get(key).get(j);
				if ((v[0] != -999) && (v[1] != -999) && (v[2] != -999)) {
					line.add(v[0] + "");
					line.add(v[1] + "");
					line.add(v[2] + "");

				} else {
					IntStream.range(0, 3).forEach(k -> line.add(""));
				}

			}
			String[] _line = new String[line.size()];
			IntStream.range(0, line.size()).forEach(k -> _line[k] = line.get(k));
			writer.writeNext(_line);
		}
		// closing writer connection
		writer.close();
	}
}
