package com.gnssAug.Rinex.fileParser;

import java.io.File;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;

import org.ejml.simple.SimpleMatrix;

import com.gnssAug.Android.constants.Constellation;
import com.gnssAug.utility.StringUtil;

public class DCB_Bias {
    private final static double SpeedofLight = 299792458;
    private String path;
    private HashMap<Character, HashMap<Integer, HashMap<String, HashMap<String, Double>>>> biasMap;
    private Map<String, Integer> GPSindexMap;
    private HashMap<Integer, double[]> C1Wmap;

    public DCB_Bias(String path) throws Exception {
        this.path = path;
        // FIX 1: Expanded map to include new CAS signals (C1L, C1X, etc.)
        // This prevents the map.get() returning null for valid CAS files.
        GPSindexMap = Map.ofEntries(
            Map.entry("C1C", 0), Map.entry("C1W", 1), Map.entry("C2C", 2),
            Map.entry("C2W", 3), Map.entry("C2S", 4), Map.entry("C2L", 5),
            Map.entry("C2X", 6), Map.entry("C5Q", 7), Map.entry("C5X", 8),
            Map.entry("C1L", 9), Map.entry("C1X", 10)
        );
        C1Wmap = new HashMap<Integer, double[]>();
        bsx_process();
        System.out.println("DCB Bias Parsing Completed.");
    }

    private void bsx_process() throws Exception {
        try {
            File file = new File(path);
            // FIX 2: Increased Matrix size from 9x9 to 12x12 to accommodate new signals
            HashMap<Integer, double[][]> GPSBiasMat = new HashMap<Integer, double[][]>();
            biasMap = new HashMap<Character, HashMap<Integer, HashMap<String, HashMap<String, Double>>>>();
            
            Scanner input = new Scanner(file);
            input.useDelimiter("\\+BIAS/SOLUTION|\\-BIAS/SOLUTION");
            if(input.hasNext()) input.next(); // Skip header
            
            if(input.hasNext()) {
                Scanner biasSoln = new Scanner(input.next());
                biasSoln.useDelimiter("\n");
                if(biasSoln.hasNext()) biasSoln.next(); // Skip column headers
                
                while (biasSoln.hasNext()) {
                    String line = biasSoln.next();
                    if(line.startsWith("*") || line.startsWith("-")) continue; // Skip comments/footer
                    
                    int len = line.length();
                    String[] soln;
                    // Keep your existing splitter logic
                    if (len >= 137) { 
                        soln = StringUtil.splitter(line, 5, 5, 4, 10, 5, 5, 15, 15, 5, 22, 12, 22, 12);
                    } else {
                        soln = StringUtil.splitter(line, 5, 5, 4, 10, 5, 5, 15, 15, 5, 22, 12);
                    }

                    if (soln != null && soln.length > 9 && soln[0].trim().equals("DSB")) {
                        // FIX 3: Check isBlank on the STATION field (soln[3])
                        // Satellite DCBs usually have empty station fields.
                        if (soln[3].isBlank()) {
                            // FIX 4: Added .trim() to handle loose column formatting in CAS files
                            String PRNstr = soln[2].trim(); 
                            char SSI = PRNstr.charAt(0);
                            int prn = Integer.parseInt(PRNstr.substring(1));
                            
                            String obs1 = soln[4].trim();
                            String obs2 = soln[5].trim();
                            double biasValue = Double.parseDouble(soln[9]) * 1e-9;

                            biasMap.computeIfAbsent(SSI,
                                    k -> new HashMap<Integer, HashMap<String, HashMap<String, Double>>>())
                                    .computeIfAbsent(prn, k -> new HashMap<String, HashMap<String, Double>>())
                                    .computeIfAbsent(obs1, k -> new HashMap<String, Double>()).put(obs2, biasValue);

                            if (SSI == 'G') {
                                // FIX 5: Null safety check. If the file has a weird signal not in our map, skip it.
                                Integer i = GPSindexMap.get(obs1);
                                Integer j = GPSindexMap.get(obs2);

                                if (i != null && j != null) {
                                    double[][] mat = GPSBiasMat.computeIfAbsent(prn, k -> new double[12][12]);
                                    mat[i][j] = biasValue;
                                    mat[j][i] = -biasValue;
                                }
                            }
                        }
                    }
                }
                biasSoln.close();
            }
            input.close();
            computeAllGPSBias(GPSBiasMat);

        } catch (Exception e) {
            e.printStackTrace();
            throw new Exception("Error occured during parsing of Bias(.BSX) file \n" + e.toString());
        }
    }

    public void computeAllGPSBias(HashMap<Integer, double[][]> GPSBiasMat) {
        for (int prn : GPSBiasMat.keySet()) {
            // FIX 6: Array size 12 to match map
            double[] C1W = new double[12];
            double[][] mat = GPSBiasMat.get(prn);
            
            // Legacy Calculations (Unchanged)
            C1W[0] = mat[1][0];
            C1W[1] = 0;
            C1W[2] = mat[0][3] + mat[3][2];
            C1W[3] = mat[1][3];
            C1W[4] = mat[1][3] + mat[3][4];
            C1W[5] = mat[1][3] + mat[3][5];
            C1W[6] = mat[1][3] + mat[3][6];
            C1W[7] = mat[1][0] + mat[0][7];
            C1W[8] = mat[1][0] + mat[0][8];
            
            // New Signal Calculations (Safe for old files because mat values defaults to 0.0)
            C1W[9] = mat[1][0] + mat[0][9];  // C1L
            C1W[10] = mat[1][0] + mat[0][10]; // C1X

            C1Wmap.put(prn, C1W);
        }
    }

    public double getISC(String obsvCode, int PRN) throws Exception {
        char SSI = obsvCode.charAt(0);
        String _obsvCode = 'C' + obsvCode.substring(1);
        Double ISC = 0.0;
        
        if (SSI == 'G') {
            Integer index = GPSindexMap.get(_obsvCode);
            if(index != null && C1Wmap.containsKey(PRN)) {
                ISC = C1Wmap.get(PRN)[index];
            }
        } else if (SSI == 'E') {
            if(biasMap.containsKey(SSI) && biasMap.get(SSI).containsKey(PRN)) {
               HashMap<String, HashMap<String, Double>> galileoMap = biasMap.get(SSI).get(PRN);
               
               if (_obsvCode.equals("C5Q")) {
                   ISC = galileoMap.getOrDefault("C1C", new HashMap<>()).getOrDefault("C5Q", null);
               } else if (_obsvCode.equals("C5X")) {
                   ISC = galileoMap.getOrDefault("C1X", new HashMap<>()).getOrDefault("C5X", null);
               }
            }
        } else if (SSI == 'C') {
			
			HashMap<String, HashMap<String, Double>> beidouMap = biasMap.get(SSI).get(PRN);
			// Calculate Frequency Ratios
	        double f2 = Constellation.frequency.get('C').get(2); // B1I
	        double f7 = Constellation.frequency.get('C').get(7); // B2b/B2I
	        double f6 = Constellation.frequency.get('C').get(6); // B3I
	        
	        double bds2FreqRatio = Math.pow(f2, 2) / Math.pow(f7, 2);
	        double bds13FreqRatio = Math.pow(f2, 2) / Math.pow(f6, 2);
	        
	        // Common Bias Term: C2I -> C6I (B1I - B3I)
	        // Used to align with the Iono-Free standard
	        double bias_C2I_C6I = beidouMap.get("C2I").getOrDefault("C6I", null);
	        double ISC_ionoFree_C6I = -(bds13FreqRatio) * bias_C2I_C6I / (1 - bds13FreqRatio);
	        
			if(_obsvCode.equals("C2I"))
			{
				ISC = -bias_C2I_C6I / (1 - bds13FreqRatio);
			}
			else if(_obsvCode.equals("C5X"))
			{
				double ISC_C6I_C1X = -beidouMap.get("C1X").getOrDefault("C6I", null);
				double ISC_C1X_C5X =  beidouMap.get("C1X").getOrDefault("C5X", null);
				ISC = ISC_ionoFree_C6I+ISC_C6I_C1X+ISC_C1X_C5X;
			}
			else if (_obsvCode.equals("C5P")) {
	            // *** NEW BLOCK FOR C5P ***
	            // Assuming C5P follows the same path as C5X (via C1X)
	            
	            // 1. C6I -> C1P
	            double bias_C1P_C6I = beidouMap.get("C1P").getOrDefault("C6I", null);
	            double ISC_C6I_C1P = -bias_C1P_C6I;

	            // 2. C1P -> C5P (Try looking up C1X->C5P)
	            double bias_C1P_C5P = beidouMap.get("C1P").getOrDefault("C5P",null);
	           

	            ISC = ISC_ionoFree_C6I + ISC_C6I_C1P + bias_C1P_C5P;
	        }
			else if(_obsvCode.equals("C1P"))
			{
				double ISC_C6I_C1P = -beidouMap.get("C1P").getOrDefault("C6I", null);
				
				ISC = ISC_ionoFree_C6I+ISC_C6I_C1P;
			}
			else {
				throw new Exception("Error in Beidou DCB initialization");
			}
		}
        return ISC;
    }
}
