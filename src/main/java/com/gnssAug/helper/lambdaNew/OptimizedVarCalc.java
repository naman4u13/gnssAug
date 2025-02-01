package com.gnssAug.helper.lambdaNew;

import org.ejml.simple.SimpleMatrix;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

public class OptimizedVarCalc {

    /**
     * Calculates the variance matrix based on the provided candidates and samples.
     *
     * @param allCandidatesKey List of candidate vectors (column vectors as SimpleMatrix objects).
     * @param aFixAll       Matrix containing all fixed vectors (column-wise SimpleMatrix).
     * @param nSamples      Number of samples (columns in aFixAll).
     * @param tolerance     Numerical tolerance for comparing matrices.
     * @return Variance matrix (SimpleMatrix).
     */
    public static Object[] calculateVariance(HashSet<String> allCandidatesKey, SimpleMatrix aFixAll, int nSamples, double tolerance) {
        int nn = aFixAll.numRows();
        String zeroKey ="";
        for(int i=0;i<nn;i++)
        {
        	zeroKey += "0.0,";
        }
        SimpleMatrix variance = new SimpleMatrix(nn, nn);
        
        // Input validation to avoid edge cases
        if (allCandidatesKey == null || allCandidatesKey.isEmpty() || aFixAll == null || aFixAll.numCols() == 0 || nSamples <= 0) {
            throw new IllegalArgumentException("Invalid inputs: Candidates or samples cannot be empty, and nSamples must be > 0.");
        }
        double failRate = 0;
        // Step 1: Preprocess all samples into a map for quick lookup with rounding to handle floating-point precision
        Map<String, Integer> sampleCounts = new HashMap<>();
        for (int jj = 0; jj < nSamples; jj++) {
            SimpleMatrix sample = aFixAll.extractVector(false, jj); // Extract column vector
            
            String sampleKey = serializeMatrix(sample, tolerance); // Serialize the matrix into a key with rounding
            if(areAllElementsIntegers(sample)&&!sampleKey.equals(zeroKey))
            {
            	failRate++;
            }
            sampleCounts.put(sampleKey, sampleCounts.getOrDefault(sampleKey, 0) + 1);
        }
        
        int zeroSampleCount = sampleCounts.getOrDefault(zeroKey, 0);
        double successRate = (zeroSampleCount*1.0)/nSamples;
        // Step 2: Process each candidate and compute contributions
        for (String candidateKey : allCandidatesKey) {
             // Serialize candidate to a key
            int count = sampleCounts.getOrDefault(candidateKey, 0); // Get the count from the map
            SimpleMatrix candidate = deserializeMatrix(candidateKey, nn);
            if (count > 0) {
                // Contribution to the variance matrix
                variance = variance.plus(candidate.mult(candidate.transpose()).scale((count * 1.0) / nSamples));
            }
            
        }
        
        // Step 3: Enforce symmetry on the variance matrix (optional but useful for numerical stability)
        variance = enforceSymmetry(variance);
        failRate = failRate/nSamples;
        return new Object[] {variance,successRate,failRate};
    }

    /**
     * Serializes a SimpleMatrix to a key with rounding to avoid numerical precision issues.
     *
     * @param matrix    The matrix to serialize (assumes column vector).
     * @param tolerance Numerical tolerance for rounding.
     * @return A serialized string key representing the matrix.
     */
    private static String serializeMatrix(SimpleMatrix matrix, double tolerance) {
        StringBuilder keyBuilder = new StringBuilder();
        for (int i = 0; i < matrix.numRows(); i++) {
            // Round each element to the specified tolerance
        	if(Double.isNaN(matrix.get(i, 0)))
        	{
        		return "NaN";
        	}
            double roundedValue = Math.round(matrix.get(i, 0) / tolerance) * tolerance;
            keyBuilder.append(roundedValue).append(",");
        }
        return keyBuilder.toString();
    }
    
    private static SimpleMatrix deserializeMatrix(String serializedMatrix, int numRows) {
        // Split the serialized string into individual values using a comma as the delimiter
        String[] values = serializedMatrix.split(",");

        // Ensure the number of rows matches the number of values in the serialized matrix
        if (values.length != numRows) {
            throw new IllegalArgumentException("The number of rows does not match the serialized matrix length.");
        }

        // Create a new SimpleMatrix with the specified number of rows and 1 column
        SimpleMatrix matrix = new SimpleMatrix(numRows, 1);

        // Parse each value and set it into the matrix
        for (int i = 0; i < numRows; i++) {
            double value = Double.parseDouble(values[i]); // Convert the string value back to a double
            matrix.set(i, 0, value); // Set the value in the matrix
        }

        return matrix; // Return the reconstructed SimpleMatrix
    }

    /**
     * Enforces symmetry on a matrix by averaging it with its transpose.
     *
     * @param matrix The input matrix (SimpleMatrix).
     * @return A symmetric version of the matrix.
     */
    private static SimpleMatrix enforceSymmetry(SimpleMatrix matrix) {
        return matrix.plus(matrix.transpose()).scale(0.5);
    }
    
    public static boolean areAllElementsIntegers(SimpleMatrix matrix) {
        int rows = matrix.numRows();
        int cols = matrix.numCols();

        // Iterate over all elements and check if they are integers
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                double element = matrix.get(i, j);
                // Check if the element is an integer (i.e., it's equal to its integer cast)
                if (element != (int) element) {  // If element is not equal to its integer part
                    return false;
                }
            }
        }
        return true;  // All elements are integers
    }
}