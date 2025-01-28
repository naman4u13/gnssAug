package com.gnssAug.helper.lambdaNew;

import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.decomposition.chol.CholeskyDecompositionLDL_DDRM;
import org.ejml.dense.row.factory.DecompositionFactory_DDRM;
import org.ejml.interfaces.decomposition.CholeskyDecomposition_F64;
import org.ejml.simple.SimpleMatrix;

import com.gnssAug.helper.lambdaNew.DecomposeLtDL.DecompositionResult;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorIA_FFRT;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorIB;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorILS;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorIR;

import java.util.Arrays;
import java.util.Random;

public class ComputeSRNumerical {

    /**
     * Computes the success rate based on numerical simulations.
     *
     * @param QMat       Variance-Covariance matrix of the ambiguities
     * @param decorr     Flag to decorrelate ambiguities (0 = No, 1 = Yes)
     * @param nSamples   Number of samples for the numerical simulation
     * @param estimator  Estimator used in the numerical simulation
     * @param config     Configuration parameters for certain estimators
     * @return Result object containing Success Rate (SR), Failure Rate (FR), and Computational Time (timeCPU)
     * @throws IllegalArgumentException if not enough inputs are provided or if an invalid estimator is selected
     */
    public static Result computeSRNumerical(SimpleMatrix QMat, int decorr, int nSamples, int estimator, Object config) {
        // Problem dimensionality
        int nn = QMat.numCols();

        // Handle default parameters based on number of provided arguments
        // Since Java does not support variable number of arguments in the same way as MATLAB,
        // we assume that the caller provides all necessary arguments, possibly using default values.

        // Parameters for specific estimators
        int[] dimBlocks = null;
        double maxFR = 0.0;
        double betaIAB = 0.0;

        if (estimator == 2 || estimator == 4 || estimator == 6 || estimator == 7) {
            if (config == null) {
                // Set default values based on the estimator type
                dimBlocks = new int[]{nn / 2, (int) Math.ceil((double) nn / 2)};
                maxFR = 0.5 / 100.0;
                betaIAB = 0.7;
            } else {
                // Assuming config is an instance of Config class with appropriate fields
                if (estimator == 2 || estimator == 4) {
                    dimBlocks = (int[]) config;
                } else if (estimator == 6) {
                    maxFR = (double) config;
                } else if (estimator == 7) {
                    betaIAB = (double) config;
                }
            }
        }

        SimpleMatrix LMat;
        double[] dVec;

        // Z-transformation (optional) for decorrelating the ambiguities
        if (decorr == 1) {
            // Decorrelate the ambiguity vc-matrix and return a LtDL-decomposition
        	DecorrelateVCResult decorResult = DecorrelateVC.decorrelateVC(QMat.copy(),null);
            QMat = decorResult.getQzHat();
            LMat = decorResult.getLzMat();
            dVec = decorResult.getDzVec();
        } else {
            // Retrieve LtDL-decomposition of the original ambiguity vc-matrix
            DecompositionResult ltldlResult = DecomposeLtDL.decomposeLtDL(QMat.copy());
            LMat = ltldlResult.getLMat();
            dVec = ltldlResult.getDVec();
        }

        // Generation of numerical samples
        if (nSamples == 0) {
            // Get an approximative value of success rate based on IB formulation
            double P0 = ComputeSR_IBexact.computeSR_IBexact(dVec).getSR();

            // Compute the number of samples to be used
            nSamples = Utilities.computeNumSamples(P0);
        }

        // Initialize random number generator
        Random rand = new Random();

        // Create a Cholesky decomposition instance
        CholeskyDecompositionLDL_DDRM chol = new CholeskyDecompositionLDL_DDRM();

        // Convert SimpleMatrix to DMatrixRMaj for decomposition
        DMatrixRMaj matrix = QMat.getMatrix().copy();

        // Perform the decomposition
        if (!chol.decompose(matrix)) {
            throw new IllegalArgumentException("Matrix is not positive definite.");
        }

        // Retrieve the lower triangular matrix
        SimpleMatrix cholQMat = new SimpleMatrix(chol.getL());

        // Now cholQMat contains the lower triangular matrix from the Cholesky decomposition

        // Now cholQMat contains the lower triangular matrix from the Cholesky decomposition
        SimpleMatrix aHatAll = cholQMat.transpose().mult(generateRandn(nn, nSamples, rand));

        // Initialize all the ambiguity fixed vectors
        SimpleMatrix aFixAll = new SimpleMatrix(nn, nSamples);
        for (int i = 0; i < aFixAll.getNumElements(); i++) {
            aFixAll.set(i, Double.NaN);
        }

        // ESTIMATORS: numerical simulations for computing the success rate
        long startTime = System.nanoTime();

        switch (estimator) {
            case 1: // Use ILS
                for (int ii = 0; ii < nSamples; ii++) {
                    SimpleMatrix aHat = aHatAll.extractVector(false, ii);
                    SimpleMatrix aFix = EstimatorILS.estimatorILS(aHat, LMat, dVec, 1).getAFix();
                    for (int row = 0; row < aFix.numRows(); row++) {
                        aFixAll.set(row, ii, aFix.get(row));
                    }
                }
                break;

//            case 2: // Use VIB-ILS
//                for (int ii = 0; ii < nSamples; ii++) {
//                    SimpleMatrix aHat = aHatAll.extractVector(false, ii);
//                    SimpleMatrix aFix = estimatorVIB(aHat, LMat, dVec, "ILS", dimBlocks);
//                    for (int row = 0; row < aFix.numRows(); row++) {
//                        aFixAll.set(row, ii, aFix.get(row));
//                    }
//                }
//                break;

            case 3: // Use IB
                for (int ii = 0; ii < nSamples; ii++) {
                    SimpleMatrix aHat = aHatAll.extractVector(false, ii);
                    SimpleMatrix aFix = EstimatorIB.estimatorIB(aHat, LMat);
                    for (int row = 0; row < aFix.numRows(); row++) {
                        aFixAll.set(row, ii, aFix.get(row));
                    }
                }
                break;

//            case 4: // Use VIB-IR
//                for (int ii = 0; ii < nSamples; ii++) {
//                    SimpleMatrix aHat = aHatAll.extractVector(false, ii);
//                    SimpleMatrix aFix = estimatorVIB(aHat, LMat, dVec, "IR", dimBlocks);
//                    for (int row = 0; row < aFix.numRows(); row++) {
//                        aFixAll.set(row, ii, aFix.get(row));
//                    }
//                }
//                break;

            case 5: // Use IR
                for (int ii = 0; ii < nSamples; ii++) {
                    SimpleMatrix aHat = aHatAll.extractVector(false, ii);
                    SimpleMatrix aFix = EstimatorIR.estimatorIR(aHat);
                    for (int row = 0; row < aFix.numRows(); row++) {
                        aFixAll.set(row, ii, aFix.get(row));
                    }
                }
                break;

            case 6: // Use IA-FFRT (ILS w/ Fixed Failure-rate Ratio Test)
                for (int ii = 0; ii < nSamples; ii++) {
                    SimpleMatrix aHat = aHatAll.extractVector(false, ii);
                    SimpleMatrix aFix = EstimatorIA_FFRT.estimatorIA_FFRT(aHat, LMat, dVec, maxFR,null).getaFix();
                    for (int row = 0; row < aFix.numRows(); row++) {
                        aFixAll.set(row, ii, aFix.get(row));
                    }
                }
                break;

//            case 7: // Compute IAB
//                for (int ii = 0; ii < nSamples; ii++) {
//                    SimpleMatrix aHat = aHatAll.extractVector(false, ii);
//                    SimpleMatrix aFix = estimatorIAB(aHat, LMat, dVec, betaIAB);
//                    for (int row = 0; row < aFix.numRows(); row++) {
//                        aFixAll.set(row, ii, aFix.get(row));
//                    }
//                }
//                break;

            default:
                throw new IllegalArgumentException("ATTENTION: the estimator selected is not available! Use 1-7.");
        }

        double timeCPU = (System.nanoTime() - startTime) / 1e9 / nSamples;

        // Success Rate (SR) from numerical simulations
        int successCount = 0;
        for (int ii = 0; ii < nSamples; ii++) {
            boolean allZero = true;
            for (int jj = 0; jj < nn; jj++) {
                if (aFixAll.get(jj, ii) != 0.0) {
                    allZero = false;
                    break;
                }
            }
            if (allZero) {
                successCount++;
            }
        }
        double SR = (double) successCount / nSamples;

        // Failure Rate (FR) from numerical simulations
        int failureCount = 0;
        for (int ii = 0; ii < nSamples; ii++) {
            boolean anyNonZero = false;
            boolean allInteger = true;
            for (int jj = 0; jj < nn; jj++) {
                double val = aFixAll.get(jj, ii);
                if (val != 0.0) {
                    anyNonZero = true;
                }
                if (val != Math.round(val)) {
                    allInteger = false;
                    break;
                }
            }
            if (anyNonZero && allInteger) {
                failureCount++;
            }
        }
        double FR = (double) failureCount / nSamples;

        return new Result(SR, FR, timeCPU);
    }

    // Helper method to generate a matrix of random Gaussian numbers
    private static SimpleMatrix generateRandn(int rows, int cols, Random rand) {
        double[] data = new double[rows * cols];
        for (int i = 0; i < data.length; i++) {
            data[i] = rand.nextGaussian();
        }
        return new SimpleMatrix(rows, cols, true, data);
    }

    

    // Result class to hold output values
    public static class Result {
        public final double SR;
        public final double FR;
        public final double timeCPU;

        public Result(double SR, double FR, double timeCPU) {
            this.SR = SR;
            this.FR = FR;
            this.timeCPU = timeCPU;
        }
    }

    
}
