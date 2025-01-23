package com.gnssAug.helper.lambdaNew.models;

import java.util.Arrays;

public class LambdaResult {
    private double[] aFix;         // Fixed ambiguity vector (column)
    private double sqNorm;         // Squared norm of the ambiguity residuals
    private int nFixed;            // Number of integer-fixed ambiguity components
    private double successRate;    // Success rate (bootstrapping) for full ambiguity resolution
    private double[][] ZMatrix;    // Admissible Z-transformation matrix (unimodular)
    private double[][] QzHat;      // Variance-covariance matrix of the decorrelated ambiguities

    /**
     * Constructor for the LambdaResult class.
     *
     * @param aFix        Fixed ambiguity vector
     * @param sqNorm      Squared norm of the ambiguity residuals
     * @param nFixed      Number of integer-fixed ambiguity components
     * @param successRate Success rate for full ambiguity resolution
     * @param ZMatrix     Z-transformation matrix
     * @param QzHat       Variance-covariance matrix of decorrelated ambiguities
     */
    public LambdaResult(double[] aFix, double sqNorm, int nFixed, double successRate, double[][] ZMatrix, double[][] QzHat) {
        this.aFix = aFix;
        this.sqNorm = sqNorm;
        this.nFixed = nFixed;
        this.successRate = successRate;
        this.ZMatrix = ZMatrix;
        this.QzHat = QzHat;
    }

    // Getters and Setters

    public double[] getAFix() {
        return aFix;
    }

    public void setAFix(double[] aFix) {
        this.aFix = aFix;
    }

    public double getSqNorm() {
        return sqNorm;
    }

    public void setSqNorm(double sqNorm) {
        this.sqNorm = sqNorm;
    }

    public int getNFixed() {
        return nFixed;
    }

    public void setNFixed(int nFixed) {
        this.nFixed = nFixed;
    }

    public double getSuccessRate() {
        return successRate;
    }

    public void setSuccessRate(double successRate) {
        this.successRate = successRate;
    }

    public double[][] getZMatrix() {
        return ZMatrix;
    }

    public void setZMatrix(double[][] ZMatrix) {
        this.ZMatrix = ZMatrix;
    }

    public double[][] getQzHat() {
        return QzHat;
    }

    public void setQzHat(double[][] QzHat) {
        this.QzHat = QzHat;
    }

    // Utility methods for debugging and visualization

    /**
     * Converts the result to a readable string format for debugging.
     *
     * @return String representation of the LambdaResult object
     */
    @Override
    public String toString() {
        return "LambdaResult{" +
                "aFix=" + Arrays.toString(aFix) +
                ", sqNorm=" + sqNorm +
                ", nFixed=" + nFixed +
                ", successRate=" + successRate +
                ", ZMatrix=" + matrixToString(ZMatrix) +
                ", QzHat=" + matrixToString(QzHat) +
                '}';
    }

    /**
     * Converts a 2D array (matrix) to a readable string.
     *
     * @param matrix The 2D array to be converted
     * @return String representation of the matrix
     */
    private String matrixToString(double[][] matrix) {
        if (matrix == null) return "null";
        StringBuilder sb = new StringBuilder("[");
        for (double[] row : matrix) {
            sb.append(Arrays.toString(row)).append("\n");
        }
        sb.append("]");
        return sb.toString();
    }
}