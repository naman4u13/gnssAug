package com.gnssAug.utility;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.collections.set.ListOrderedSet;
import org.ejml.simple.SimpleMatrix;

import com.gnssAug.Android.constants.Measurement;
import com.gnssAug.Android.estimation.KalmanFilter.EKF_PPP3;
import com.gnssAug.Rinex.estimation.EKF_PPP;
import com.gnssAug.Rinex.models.SatResidual;
import com.gnssAug.Rinex.models.Satellite;
import com.opencsv.CSVWriter;

public class MakeCSV {

	public static void exportPPPresultsToCSV(TreeMap<Long, ArrayList<Satellite>> satMap, EKF_PPP ekf,
			HashMap<String, ArrayList<double[]>> GraphPosMap,
			HashMap<Measurement, HashMap<String, HashMap<String, ArrayList<SatResidual>>>> satInnMap,
			ArrayList<Long> timeList, String[] obsvCodeList, double[] refEcef) {

		String basePath = "/Users/naman.agarwal/Library/CloudStorage/OneDrive-UniversityofCalgary/Phd Proposal/Plots_CSV/Pixel4_Nov/GPS_GAL_L1_L5_PPP_noRepair_noSoftwareCSD/";
		new File(basePath).mkdirs(); // Ensure folder exists
		ListOrderedSet ssiSet = new ListOrderedSet();
		for (int i = 0; i < obsvCodeList.length; i++) {
			ssiSet.add(obsvCodeList[i].charAt(0) + "");
		}
		String[] ssiLabels = (String[]) ssiSet.toArray(new String[0]);
		long t0 = timeList.get(0);
		// 1. Clock Offset
		exportDenseCSV(basePath + "PPP_ClockOffset.csv", ekf.getClkOffMap(), obsvCodeList);

		// 2. Clock Drift
		exportDenseCSV(basePath + "PPP_ClockDrift.csv", ekf.getClkDriftMap(), ssiLabels);

		// 3. Ionosphere
		exportSparseCSV(basePath + "PPP_Iono.csv", ekf.getIonoMap());

		// 4. Troposphere
		exportScalarCSV(basePath + "PPP_Tropo.csv", ekf.getTropoMap());

		// 5. Ambiguities
		exportSparseCSV(basePath + "PPP_Ambiguity.csv", ekf.getAmbMap());

		// 6. Innovations or Prefit Residuals
		exportSatResToCSV(basePath + "Innovation", satInnMap);

		// 7. ENU positioning Error
		exportENUToCSV(basePath + "PPP_ENU.csv", GraphPosMap, timeList);

		// 7. DOP and Satellite count
		exportQualityMetricsToCSV(basePath + "PPP_DOP_SatCount.csv", satMap, refEcef);

	}
	// Inside MakeCSV.java

	public static void exportAndroidPPPToCSV(TreeMap<Long, ArrayList<com.gnssAug.Android.models.Satellite>> satMap, EKF_PPP3 ekf,
			HashMap<String, ArrayList<double[]>> GraphPosMap,
			HashMap<Measurement, HashMap<String, HashMap<String, ArrayList<SatResidual>>>> satInnMap,
			ArrayList<Long> timeList, String[] obsvCodeList, double[] refEcef) {

		String basePath = "/Users/naman.agarwal/Library/CloudStorage/OneDrive-UniversityofCalgary/Phd Proposal/Plots_CSV/Hervanta_SamsungA52/GPS_GAL_BEI_L1_PPP_unRepaired/";
		new File(basePath).mkdirs();
		ListOrderedSet ssiSet = new ListOrderedSet();
		for (int i = 0; i < obsvCodeList.length; i++) {
			ssiSet.add(obsvCodeList[i].charAt(0) + "");
		}
		String[] ssiLabels = (String[]) ssiSet.toArray(new String[0]);
		// --- 1. CLOCK OFFSETS (Code & Phase) ---
		// Export Code Clock
		TreeMap<Long, double[]> codeClk = ekf.getClkOffMap().get(Measurement.Pseudorange);
		if (codeClk != null) {
			exportDenseCSV(basePath + "PPP_ClockOffset_Code.csv", codeClk, obsvCodeList);
		}

		// Export Phase Clock
		TreeMap<Long, double[]> phaseClk = ekf.getClkOffMap().get(Measurement.CarrierPhase);
		if (phaseClk != null) {
			

			exportDenseCSV(basePath + "PPP_ClockOffset_Phase.csv", phaseClk, obsvCodeList);
		}

		// --- 2. CLOCK DRIFTS (Doppler & TDCP) ---
		// Export Doppler Drift
		TreeMap<Long, double[]> dopplerDrift = ekf.getClkDriftMap().get(Measurement.Doppler);
		if (dopplerDrift != null && !dopplerDrift.isEmpty()) {
			
			exportDenseCSV(basePath + "PPP_ClockDrift_Doppler.csv", dopplerDrift, ssiLabels);
		}

		// Export TDCP Drift
		TreeMap<Long, double[]> tdcpDrift = ekf.getClkDriftMap().get(Measurement.TDCP);
		if (tdcpDrift != null && !tdcpDrift.isEmpty()) {
			
			exportDenseCSV(basePath + "PPP_ClockDrift_TDCP.csv", tdcpDrift, ssiLabels);
		}

		// --- 3. OTHER EXPORTS (Unchanged) ---
		exportScalarCSV(basePath + "PPP_Tropo.csv", ekf.getTropoMap());
		exportSparseCSV(basePath + "PPP_Iono.csv", ekf.getIonoMap());
		exportSparseCSV(basePath + "PPP_Ambiguity.csv", ekf.getAmbMap());

		// 6. Innovations or Prefit Residuals
		exportSatResToCSV(basePath + "Innovation", satInnMap);

		// 7. ENU positioning Error
		exportENUToCSV(basePath + "PPP_ENU.csv", GraphPosMap, timeList);

		// 8. DOP and Satellite count
		exportAndroidQualityMetricsToCSV(basePath + "PPP_DOP_SatCount.csv", satMap, refEcef);

		System.out.println("Android PPP3 Export Completed.");
	}

	public static void exportSatResToCSV(String baseDir,
			HashMap<Measurement, HashMap<String, HashMap<String, ArrayList<SatResidual>>>> satResMap) {

		// Ensure directory exists
		new File(baseDir).mkdirs();

		// Loop through Measurements (e.g., Pseudorange, Doppler)
		for (Measurement meas : satResMap.keySet()) {
			HashMap<String, HashMap<String, ArrayList<SatResidual>>> estMap = satResMap.get(meas);

			// Loop through Estimators (e.g., WLS, EKF, PPP)
			for (String estName : estMap.keySet()) {
				HashMap<String, ArrayList<SatResidual>> satDataMap = estMap.get(estName);

				// Sanitize filename
				String cleanEstName = estName.replaceAll("[^a-zA-Z0-9.-]", "_");
				String fileName = baseDir + "_" + meas + "_" + cleanEstName + ".csv";

				try (CSVWriter writer = new CSVWriter(new FileWriter(fileName))) {
					// Write Header
					writer.writeNext(new String[] { "Time", "SatID", "Value", "Elevation", "IsOutlier", "CN0" });

					// Iterate through all satellites and their residuals
					for (String satID : satDataMap.keySet()) {
						for (SatResidual res : satDataMap.get(satID)) {
							// Assumption: SatResidual has getters like getTime(), getValue(), etc.
							// You might need to adjust these method names to match your SatResidual class
							String[] line = new String[] { String.valueOf(res.getT()), // Already relative time?
									satID, String.valueOf(res.getResidual()), String.valueOf(res.getElevAngle()), // Radians
																													// or
																													// Degrees?
									res.isOutlier() ? "1" : "0", String.valueOf(res.getCN0()) };
							writer.writeNext(line);
						}
					}
					System.out.println("Exported: " + fileName);

				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}

	// Export "Dense" data (like Clocks) where every column exists for every epoch
	public static void exportDenseCSV(String filePath, TreeMap<Long, double[]> map, String[] headers) {
		try (CSVWriter writer = new CSVWriter(new FileWriter(filePath))) {
			// Create Header: Time, G, E, ...
			String[] fullHeader = new String[headers.length + 1];
			fullHeader[0] = "Time";
			System.arraycopy(headers, 0, fullHeader, 1, headers.length);
			writer.writeNext(fullHeader);
			long t0 = map.firstKey();
			for (Map.Entry<Long, double[]> entry : map.entrySet()) {
				String[] line = new String[headers.length + 1];
				// Normalize time to seconds here to match Matlab
				line[0] = String.valueOf((entry.getKey() - t0) / 1000);
				double[] values = entry.getValue();
				for (int i = 0; i < values.length && i < headers.length; i++) {
					line[i + 1] = String.valueOf(values[i]);
				}
				writer.writeNext(line);
			}
			System.out.println("Exported: " + filePath);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	// Export "Sparse" data (like Iono/Ambiguity) as a "Long Table" (Time, ID,
	// Value)
	public static void exportSparseCSV(String filePath, TreeMap<Long, HashMap<String, Double>> map) {
		try (CSVWriter writer = new CSVWriter(new FileWriter(filePath))) {
			writer.writeNext(new String[] { "Time", "ID", "Value" });
			long t0 = map.firstKey();
			for (Map.Entry<Long, HashMap<String, Double>> entry : map.entrySet()) {
				String time = String.valueOf((entry.getKey() - t0) / 1000);
				for (Map.Entry<String, Double> inner : entry.getValue().entrySet()) {
					writer.writeNext(new String[] { time, inner.getKey(), String.valueOf(inner.getValue()) });
				}
			}
			System.out.println("Exported: " + filePath);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	// Export Scalar data (like Troposphere)
	public static void exportScalarCSV(String filePath, TreeMap<Long, Double> map) {
		try (CSVWriter writer = new CSVWriter(new FileWriter(filePath))) {
			writer.writeNext(new String[] { "Time", "Value" });
			long t0 = map.firstKey();
			for (Map.Entry<Long, Double> entry : map.entrySet()) {
				String time = String.valueOf((entry.getKey() - t0) / 1000);
				writer.writeNext(new String[] { time, String.valueOf(entry.getValue()) });
			}
			System.out.println("Exported: " + filePath);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static void exportENUToCSV(String filePath, HashMap<String, ArrayList<double[]>> dataMap,
			ArrayList<Long> timeList) {
		try (CSVWriter writer = new CSVWriter(new FileWriter(filePath))) {
			// Write Header
			writer.writeNext(new String[] { "Time", "Estimator", "East", "North", "Up" });

			long t0 = timeList.get(0); // Start time for normalization

			// Iterate through time steps
			for (int i = 0; i < timeList.size(); i++) {
				// Calculate relative time in seconds
				String time = String.valueOf((timeList.get(i) - t0));

				// Iterate through each estimator (e.g., "WLS", "EKF")
				for (String estName : dataMap.keySet()) {
					ArrayList<double[]> errList = dataMap.get(estName);

					// Safety check: ensure index exists (handles cases where lists might be
					// shorter)
					if (i < errList.size()) {
						double[] enu = errList.get(i);
						// Only write if data is not null
						if (enu != null && enu.length >= 3) {
							String[] line = new String[] { time, estName, String.valueOf(enu[0]), // East
									String.valueOf(enu[1]), // North
									String.valueOf(enu[2]) // Up
							};
							writer.writeNext(line);
						}
					}
				}
			}
			System.out.println("ENU Error CSV exported to: " + filePath);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static void exportQualityMetricsToCSV(String filePath, TreeMap<Long, ArrayList<Satellite>> satMap,
			double[] refEcef) {
		// Removed 'timeList' argument because we will use satMap's keys directly

		try (CSVWriter writer = new CSVWriter(new FileWriter(filePath))) {
			// 1. Setup Headers
			Set<String> allSignalTypes = new TreeSet<>();
			for (ArrayList<Satellite> sats : satMap.values()) {
				for (Satellite sat : sats) {
					allSignalTypes.add(sat.getObsvCode());
				}
			}
			ArrayList<String> signalList = new ArrayList<>(allSignalTypes);

			ArrayList<String> header = new ArrayList<>();
			header.add("Time");
			header.add("Total_Satellites");
			header.addAll(signalList);
			header.add("PDOP");
			header.add("HDOP");
			header.add("VDOP");
			header.add("GDOP");
			writer.writeNext(header.toArray(new String[0]));

			// Get sorted keys from the map itself
			ArrayList<Long> validTimes = new ArrayList<>(satMap.keySet());
			if (validTimes.isEmpty())
				return;

			long t0 = validTimes.get(0);

			// Pre-compute Rotation Matrix (ECEF to ENU)
			double[] lla = LatLonUtil.ecef2lla(refEcef);
			double lat = Math.toRadians(lla[0]);
			double lon = Math.toRadians(lla[1]);

			SimpleMatrix R = new SimpleMatrix(3, 3);
			R.set(0, 0, -Math.sin(lon));
			R.set(0, 1, Math.cos(lon));
			R.set(0, 2, 0);
			R.set(1, 0, -Math.sin(lat) * Math.cos(lon));
			R.set(1, 1, -Math.sin(lat) * Math.sin(lon));
			R.set(1, 2, Math.cos(lat));
			R.set(2, 0, Math.cos(lat) * Math.cos(lon));
			R.set(2, 1, Math.cos(lat) * Math.sin(lon));
			R.set(2, 2, Math.sin(lat));

			// 2. Iterate Epochs
			for (Long t : validTimes) {
				ArrayList<Satellite> sats = satMap.get(t);
				ArrayList<String> line = new ArrayList<>();

				// Time (Convert ms to relative seconds)
				line.add(String.valueOf((t - t0) / 1000.0));

				HashSet<String> uniqueSatIds = new HashSet<>();
				HashMap<String, Integer> currentSigCounts = new HashMap<>();

				// List for unique geometry satellites (Filter duplicates for DOP)
				ArrayList<Satellite> geometrySats = new ArrayList<>();
				HashSet<String> processedGeometrySats = new HashSet<>();

				if (sats != null) {
					for (Satellite sat : sats) {
						// 1. Count Signals
						String code = sat.getObsvCode();
						currentSigCounts.put(code, currentSigCounts.getOrDefault(code, 0) + 1);

						// 2. Identify Unique Satellites
						char constellation = sat.getObsvCode().charAt(0);
						String uniqueID = constellation + String.valueOf(sat.getSVID()); // Make sure using getSvid()
																							// matching your class
						uniqueSatIds.add(uniqueID);

						// 3. Filter for Geometry (DOP)
						if (!processedGeometrySats.contains(uniqueID)) {
							processedGeometrySats.add(uniqueID);
							geometrySats.add(sat);
						}
					}
				}

				line.add(String.valueOf(uniqueSatIds.size()));
				for (String sig : signalList) {
					line.add(String.valueOf(currentSigCounts.getOrDefault(sig, 0)));
				}

				// DOP Calculation
				double pdop = 99.9, hdop = 99.9, vdop = 99.9, gdop = 99.9;

				if (geometrySats.size() >= 4) {
					SimpleMatrix H = new SimpleMatrix(geometrySats.size(), 4);
					int row = 0;
					for (Satellite sat : geometrySats) {
						double[] satPos = sat.getSatEci();
						double range = Math.sqrt(Math.pow(satPos[0] - refEcef[0], 2)
								+ Math.pow(satPos[1] - refEcef[1], 2) + Math.pow(satPos[2] - refEcef[2], 2));

						double ex = (satPos[0] - refEcef[0]) / range;
						double ey = (satPos[1] - refEcef[1]) / range;
						double ez = (satPos[2] - refEcef[2]) / range;

						H.set(row, 0, -ex);
						H.set(row, 1, -ey);
						H.set(row, 2, -ez);
						H.set(row, 3, 1.0);
						row++;
					}

					try {
						SimpleMatrix Q_xyz = H.transpose().mult(H).invert();
						SimpleMatrix Q3_xyz = Q_xyz.extractMatrix(0, 3, 0, 3);
						SimpleMatrix Q_enu = R.mult(Q3_xyz).mult(R.transpose());

						double qEE = Q_enu.get(0, 0);
						double qNN = Q_enu.get(1, 1);
						double qUU = Q_enu.get(2, 2);
						double qTT = Q_xyz.get(3, 3);

						pdop = Math.sqrt(qEE + qNN + qUU);
						hdop = Math.sqrt(qEE + qNN);
						vdop = Math.sqrt(qUU);
						gdop = Math.sqrt(qEE + qNN + qUU + qTT);

					} catch (Exception e) {
						// Singular matrix
					}
				}

				line.add(String.valueOf(pdop));
				line.add(String.valueOf(hdop));
				line.add(String.valueOf(vdop));
				line.add(String.valueOf(gdop));

				writer.writeNext(line.toArray(new String[0]));
			}
			System.out.println("Quality Metrics exported to: " + filePath);

		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	public static void exportAndroidQualityMetricsToCSV(String filePath, TreeMap<Long, ArrayList<com.gnssAug.Android.models.Satellite>> satMap,
			double[] refEcef) {
		// Removed 'timeList' argument because we will use satMap's keys directly

		try (CSVWriter writer = new CSVWriter(new FileWriter(filePath))) {
			// 1. Setup Headers
			Set<String> allSignalTypes = new TreeSet<>();
			for (ArrayList<com.gnssAug.Android.models.Satellite> sats : satMap.values()) {
				for (com.gnssAug.Android.models.Satellite sat : sats) {
					allSignalTypes.add(sat.getObsvCode());
				}
			}
			ArrayList<String> signalList = new ArrayList<>(allSignalTypes);

			ArrayList<String> header = new ArrayList<>();
			header.add("Time");
			header.add("Total_Satellites");
			header.addAll(signalList);
			header.add("PDOP");
			header.add("HDOP");
			header.add("VDOP");
			header.add("GDOP");
			writer.writeNext(header.toArray(new String[0]));

			// Get sorted keys from the map itself
			ArrayList<Long> validTimes = new ArrayList<>(satMap.keySet());
			if (validTimes.isEmpty())
				return;

			long t0 = validTimes.get(0);

			// Pre-compute Rotation Matrix (ECEF to ENU)
			double[] lla = LatLonUtil.ecef2lla(refEcef);
			double lat = Math.toRadians(lla[0]);
			double lon = Math.toRadians(lla[1]);

			SimpleMatrix R = new SimpleMatrix(3, 3);
			R.set(0, 0, -Math.sin(lon));
			R.set(0, 1, Math.cos(lon));
			R.set(0, 2, 0);
			R.set(1, 0, -Math.sin(lat) * Math.cos(lon));
			R.set(1, 1, -Math.sin(lat) * Math.sin(lon));
			R.set(1, 2, Math.cos(lat));
			R.set(2, 0, Math.cos(lat) * Math.cos(lon));
			R.set(2, 1, Math.cos(lat) * Math.sin(lon));
			R.set(2, 2, Math.sin(lat));

			// 2. Iterate Epochs
			for (Long t : validTimes) {
				ArrayList<com.gnssAug.Android.models.Satellite> sats = satMap.get(t);
				ArrayList<String> line = new ArrayList<>();

				// Time (Convert ms to relative seconds)
				line.add(String.valueOf((t - t0) / 1000.0));

				HashSet<String> uniqueSatIds = new HashSet<>();
				HashMap<String, Integer> currentSigCounts = new HashMap<>();

				// List for unique geometry satellites (Filter duplicates for DOP)
				ArrayList<com.gnssAug.Android.models.Satellite> geometrySats = new ArrayList<>();
				HashSet<String> processedGeometrySats = new HashSet<>();

				if (sats != null) {
					for (com.gnssAug.Android.models.Satellite sat : sats) {
						// 1. Count Signals
						String code = sat.getObsvCode();
						currentSigCounts.put(code, currentSigCounts.getOrDefault(code, 0) + 1);

						// 2. Identify Unique Satellites
						char constellation = sat.getObsvCode().charAt(0);
						String uniqueID = constellation + String.valueOf(sat.getSvid()); // Make sure using getSvid()
																							// matching your class
						uniqueSatIds.add(uniqueID);

						// 3. Filter for Geometry (DOP)
						if (!processedGeometrySats.contains(uniqueID)) {
							processedGeometrySats.add(uniqueID);
							geometrySats.add(sat);
						}
					}
				}

				line.add(String.valueOf(uniqueSatIds.size()));
				for (String sig : signalList) {
					line.add(String.valueOf(currentSigCounts.getOrDefault(sig, 0)));
				}

				// DOP Calculation
				double pdop = 99.9, hdop = 99.9, vdop = 99.9, gdop = 99.9;

				if (geometrySats.size() >= 4) {
					SimpleMatrix H = new SimpleMatrix(geometrySats.size(), 4);
					int row = 0;
					for (com.gnssAug.Android.models.Satellite sat : geometrySats) {
						double[] satPos = sat.getSatEci();
						double range = Math.sqrt(Math.pow(satPos[0] - refEcef[0], 2)
								+ Math.pow(satPos[1] - refEcef[1], 2) + Math.pow(satPos[2] - refEcef[2], 2));

						double ex = (satPos[0] - refEcef[0]) / range;
						double ey = (satPos[1] - refEcef[1]) / range;
						double ez = (satPos[2] - refEcef[2]) / range;

						H.set(row, 0, -ex);
						H.set(row, 1, -ey);
						H.set(row, 2, -ez);
						H.set(row, 3, 1.0);
						row++;
					}

					try {
						SimpleMatrix Q_xyz = H.transpose().mult(H).invert();
						SimpleMatrix Q3_xyz = Q_xyz.extractMatrix(0, 3, 0, 3);
						SimpleMatrix Q_enu = R.mult(Q3_xyz).mult(R.transpose());

						double qEE = Q_enu.get(0, 0);
						double qNN = Q_enu.get(1, 1);
						double qUU = Q_enu.get(2, 2);
						double qTT = Q_xyz.get(3, 3);

						pdop = Math.sqrt(qEE + qNN + qUU);
						hdop = Math.sqrt(qEE + qNN);
						vdop = Math.sqrt(qUU);
						gdop = Math.sqrt(qEE + qNN + qUU + qTT);

					} catch (Exception e) {
						// Singular matrix
					}
				}

				line.add(String.valueOf(pdop));
				line.add(String.valueOf(hdop));
				line.add(String.valueOf(vdop));
				line.add(String.valueOf(gdop));

				writer.writeNext(line.toArray(new String[0]));
			}
			System.out.println("Quality Metrics exported to: " + filePath);

		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
