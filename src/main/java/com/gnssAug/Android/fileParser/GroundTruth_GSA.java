package com.gnssAug.Android.fileParser;

import java.io.FileReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;

import com.opencsv.CSVReader;

public class GroundTruth_GSA {

	private static final long NumberMilliSecondsWeek = 604800000;

	public static ArrayList<double[]> processCSV(String path) throws Exception {

		ArrayList<double[]> data = new ArrayList<double[]>();
		
			Path fileName = Path.of(path);
			String[] input = Files.readString(fileName).trim().split("\r\n|\r|\n");
			for(int i =2;i<input.length;i++)
			{
				String[] arr = input[i].trim().split("\\s+");
				double GPStime = Double.parseDouble(arr[0]);
				double weekNo = Double.parseDouble(arr[4]);
				double x = Double.parseDouble(arr[1]);
				double y = Double.parseDouble(arr[2]);
				double z = Double.parseDouble(arr[3]);
				data.add(new double[] { GPStime, weekNo, x, y, z });
			}
			
		return data;

	}
}
