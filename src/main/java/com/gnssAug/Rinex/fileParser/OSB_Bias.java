package com.gnssAug.Rinex.fileParser;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;

import com.gnssAug.utility.StringUtil;
import com.gnssAug.utility.Time;

public class OSB_Bias {
    private final static double SpeedofLight = 299792458;
    private String path;
    private HashMap<Character, HashMap<Integer, HashMap<String, ArrayList<BiasEntry>>>> osbMap;  // SSI -> PRN -> OBS -> List<BiasEntry>

    private static class BiasEntry {
        double start;
        double end;
        double value;
        double std;
        String unit;

        BiasEntry(double start, double end, double value, double std, String unit) {
            this.start = start;
            this.end = end;
            this.value = value;
            this.std = std;
            this.unit = unit;
        }
    }

    public OSB_Bias(String path) throws Exception {
        this.path = path;
        osbMap = new HashMap<>();
        bsx_process();
        System.out.println();
    }

    private void bsx_process() throws Exception {
        try {
            File file = new File(path);
            Scanner input = new Scanner(file);
            input.useDelimiter("\\+BIAS/SOLUTION|\\-BIAS/SOLUTION");
            input.next();  // Skip to solution block
            String solutionBlock = input.next();
            Scanner biasSoln = new Scanner(solutionBlock);
            biasSoln.useDelimiter("\n");
            biasSoln.next();  // Skip header line in block
            while (biasSoln.hasNext()) {
                String line = biasSoln.next().trim();
                if (line.startsWith("OSB")) {
                	if (line.startsWith("OSB")) {
                	    // Parse using fixed positions to handle both formats robustly
                	    if (line.length() < 90) continue; // Too short
                	    String svn = line.substring(4, 9).trim();  // SVN
                	    String prnStr = line.substring(9, 13).trim();  // PRN
                	    String station = line.substring(13, 23).trim();  // STATION
                	    String obs1 = line.substring(23, 28).trim();  // OBS1
                	    String obs2 = line.substring(28, 33).trim();  // OBS2
                	    String startStr = line.substring(33, 48).trim();  // START
                	    String endStr = line.substring(48, 63).trim();  // END
                	    String unit = line.substring(63, 68).trim();  // UNIT
                	    String valueStr = line.substring(69, 90).trim();  // VALUE (21 chars)
                	    String stdStr = (line.length() >= 102) ? line.substring(91, 102).trim() : "0.0";  // STD (11 chars)

                	    // Only parse satellite OSB: station empty, obs2 empty, prn not empty
                	    if (!station.isEmpty() || !obs2.isEmpty() || prnStr.isEmpty()) {
                	        continue;
                	    }

                	    char ssi = prnStr.charAt(0);
                	    int prn = Integer.parseInt(prnStr.substring(1));
                	    String obs = obs1;
                	    String[] start = startStr.split(":");
                	    String[] end = endStr.split(":");
                	    double start_GPStime = Time.getGPSTimeFromYDOY(start)[0];
                	    double end_GPStime = Time.getGPSTimeFromYDOY(end)[0];
                	    double biasValue = 0.0;
                	    if(!valueStr.equals("nan"))
                	    {
                	    	double raw = Double.parseDouble(valueStr);
                	    	biasValue = unit.equals("cyc") ? raw : raw * 1e-9; // cyc: keep as cycles; ns: convert to seconds
                	    }
                	    double std = Double.parseDouble(stdStr);

                	    osbMap.computeIfAbsent(ssi, k -> new HashMap<>())
                	          .computeIfAbsent(prn, k -> new HashMap<>())
                	          .computeIfAbsent(obs, k -> new ArrayList<>())
                	          .add(new BiasEntry(start_GPStime, end_GPStime, biasValue, std, unit));
                	}}
            }
        } catch (Exception e) {
        	e.printStackTrace();
        	throw new Exception("Error occurred during parsing of OSB (.BIA) file \n" + e);
        }
    }

    private BiasEntry findEntry(char ssi, String signalType, int PRN, double GPStime) {
        HashMap<Integer, HashMap<String, ArrayList<BiasEntry>>> sysMap = osbMap.get(ssi);
        if (sysMap == null) return null;
        HashMap<String, ArrayList<BiasEntry>> prnMap = sysMap.get(PRN);
        if (prnMap == null) return null;
        ArrayList<BiasEntry> entries = prnMap.get(signalType);
        if (entries == null) return null;
        for (BiasEntry entry : entries) {
            if (isEpochInRange(GPStime, entry.start, entry.end)) return entry;
        }
        return null;
    }

    /** Returns the OSB in meters. ns (stored as seconds) → ×c; cyc (stored as cycles) → ×wavelength. */
    public double getOSBMeters(char ssi, String signalType, int PRN, double GPStime, double wavelength) throws Exception {
        BiasEntry entry = findEntry(ssi, signalType, PRN, GPStime);
        if (entry == null) return 0.0;
        return entry.unit.equals("cyc") ? entry.value * wavelength : entry.value * SpeedofLight;
    }

    // Helper method to check if epochTime is between start and end (format YYYY:DDD:SSSSS)
    private boolean isEpochInRange(double GPStime, double start, double end) {
        // Simple string comparison assuming format is sortable
        return GPStime >= start && GPStime<=end;
    }
}