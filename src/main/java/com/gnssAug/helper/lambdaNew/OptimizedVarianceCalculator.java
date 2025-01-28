package com.gnssAug.helper.lambdaNew;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.ejml.simple.SimpleMatrix;

public class OptimizedVarianceCalculator {

    public static SimpleMatrix calculateVariance(List<SimpleMatrix> allCandidates, SimpleMatrix aFixAll, int nSamples) {
        int nn = aFixAll.numRows();
        SimpleMatrix variance = new SimpleMatrix(nn, nn);

        // Step 1: Preprocess all samples into a map for quick lookup
        Map<String, Integer> sampleCounts = new HashMap<>();
        for (int jj = 0; jj < nSamples; jj++) {
            SimpleMatrix sample = aFixAll.extractVector(false, jj); // Extract column vector
            String sampleKey = sampleToKey(sample); // Serialize the matrix into a key
            sampleCounts.put(sampleKey, sampleCounts.getOrDefault(sampleKey, 0) + 1);
        }

        // Step 2: Process each candidate and compute contributions
        for (SimpleMatrix candidate : allCandidates) {
            String candidateKey = sampleToKey(candidate); // Serialize candidate to a key
            int count = sampleCounts.getOrDefault(candidateKey, 0); // Get the count from the map

            if (count > 0) {
                variance = variance.plus(candidate.mult(candidate.transpose()).scale((count * 1.0) / nSamples));
            }
        }

        return variance;
    }

    /**
     * Converts a SimpleMatrix into a unique string key for hash map storage.
     * This key is used to compare matrices.
     */
    private static String sampleToKey(SimpleMatrix matrix) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < matrix.numRows(); i++) {
            sb.append(matrix.get(i)).append(",");
        }
        return sb.toString();
    }

    
}
