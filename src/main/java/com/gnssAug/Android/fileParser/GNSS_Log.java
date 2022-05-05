package com.gnssAug.Android.fileParser;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.TreeMap;

import com.gnssAug.Android.models.GNSSLog;
import com.gnssAug.Android.models.IMUsensor;
import com.opencsv.CSVWriter;

public class GNSS_Log {

	private static TreeMap<Long, HashMap<String, ArrayList<GNSSLog>>> gnssLogMaps = null;
	private static ArrayList<IMUsensor> imuList = null;

	public static void process(String path) throws Exception {

		gnssLogMaps = new TreeMap<Long, HashMap<String, ArrayList<GNSSLog>>>();
		imuList = new ArrayList<IMUsensor>();
		ArrayList<String[]> logs = new ArrayList<String[]>();
		try {
			Path fileName = Path.of(path);
			String[] input = Files.readString(fileName).split("#");

			for (int i = 0; i < input.length - 1; i++) {
				if (input[i].trim().isBlank()) {
					continue;
				}
				String[] fields = input[i].trim().split(",");
				if (fields[0].equals("UncalAccel")) {

					break;
				}
			}
			String[] lines = input[input.length - 1].trim().split("\r\n|\r|\n");
			long bootGPStime = 0;

			// Map<String, AndroidSensor> imuMap = Map.of("Unc",
			for (String line : lines) {
				String[] data = line.trim().split(",");
				if (data[0].equals("Raw")) {
					GNSSLog log = new GNSSLog(data);
					long tRx = Math.round(log.gettRx() * 1e3);
					String obsvCode = log.getObsvCode();
					int svid = log.getSvid();
					gnssLogMaps.computeIfAbsent(tRx, k -> new HashMap<String, ArrayList<GNSSLog>>())
							.computeIfAbsent(obsvCode, k -> new ArrayList<GNSSLog>()).add(log);
					logs.add(log.toString().split(","));
					if ((log.getBootGPStime() - bootGPStime) / 1e6 > 1) {
						if (bootGPStime == 0) {
							bootGPStime = log.getBootGPStime();
						} else {
							System.err.println("ERROR in Android GNSS LOG elapsedtime computation");
							throw new Exception("ERROR in Android GNSS LOG elapsedtime computation");
						}

					}

				}
				if (data[0].equals("UncalAccel") || data[0].equals("UncalGyro") || data[0].equals("UncalMag")) {
					IMUsensor imu = new IMUsensor(data);
					imuList.add(imu);
				}
			}
			for (IMUsensor imu : imuList) {
				imu.settRx(bootGPStime);
			}
			Collections.sort(imuList, (o1, o2) -> (int) (o1.gettRx() - o2.gettRx()));

			// recordGNSSLogCSV(logs);
		} catch (Exception e) {
			// TODO: handle exception
			e.printStackTrace();
			throw new Exception(e);
		}

	}

	public static TreeMap<Long, HashMap<String, ArrayList<GNSSLog>>> getGnssLogMaps() {
		return gnssLogMaps;
	}

	public static ArrayList<IMUsensor> getImuList() {
		return imuList;
	}

	public static void recordGNSSLogCSV(ArrayList<String[]> logs) {
		String headerLine = "utcTimeMillis,TimeNanos,LeapSecond,TimeUncertaintyNanos,FullBiasNanos,BiasNanos,BiasUncertaintyNanos,DriftNanosPerSecond,DriftUncertaintyNanosPerSecond,HardwareClockDiscontinuityCount,Svid,TimeOffsetNanos,State,ReceivedSvTimeNanos,ReceivedSvTimeUncertaintyNanos,Cn0DbHz,PseudorangeRateMetersPerSecond,PseudorangeRateUncertaintyMetersPerSecond,AccumulatedDeltaRangeState,AccumulatedDeltaRangeMeters,AccumulatedDeltaRangeUncertaintyMeters,CarrierFrequencyHz,CarrierCycles,CarrierPhase,CarrierPhaseUncertainty,MultipathIndicator,SnrInDb,ConstellationType,AgcDb,BasebandCn0DbHz,FullInterSignalBiasNanos,FullInterSignalBiasUncertaintyNanos,SatelliteInterSignalBiasNanos,SatelliteInterSignalBiasUncertaintyNanos,CodeType,ChipsetElapsedRealtimeNanos,obsvCode,tRx,bootGPStime";
		String[] header = headerLine.split(",");

		String filePath = "C:\\D drive\\Study\\Google Decimeter Challenge\\decimeter\\train\\2021-04-29-US-SJC-2\\Pixel4\\log2.csv";
		CSVWriter writer = null;

		File file = new File(filePath);
		// create FileWriter object with file as parameter
		FileWriter outputfile;
		try {
			outputfile = new FileWriter(file);
			// create CSVWriter object filewriter object as parameter
			writer = new CSVWriter(outputfile);
			writer.writeNext(header);
			writer.writeAll(logs);
			writer.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

}
