package com.gnssAug.helper.lambdaNew;



import org.ejml.simple.SimpleMatrix;

import com.gnssAug.helper.lambdaNew.DecomposeLtDL.DecompositionResult;
import com.gnssAug.helper.lambdaNew.TransformZ.TransformResult;

public class DecorrelateVC {

    /**
     * Decorrelate the ambiguity vc-matrix by a Z-transformation.
     * This function computes a decorrelation of the ambiguity vc-matrix, which 
     * is firstly decomposed in its LtDL form, and the latter is updated based 
     * on an admissible Z-transformation (reduction and ordering of conditional 
     * variances). This transformation matrix (unimodular) is also provided in
     * output as inv(Z'), later used for a straightforward back-transformation.
     *
     * @param qaHat      Variance-covariance matrix of the original ambiguities
     * @param aHat       Ambiguity float vector (column)                 [OPTIONAL]
     * @return DecorrelateVCResult containing:
     *         QzHat      Variance-covariance matrix of the decorrelated ambiguities
     *         LzMat      New LtDL-decomposition matrix L (lower unitriangular)
     *         dzVec      New LtDL-decomposition matrix D (diagonal elements)
     *         iZtMat     Inverse transpose of Z-transformation matrix (unimodular)
     *         ZMat       Z-transformation matrix (unimodular)            [OPTIONAL]
     *         zHat       Decorrelated ambiguity float vector (column)    [OPTIONAL]
     */
    public static DecorrelateVCResult decorrelateVC(SimpleMatrix qaHat, SimpleMatrix aHat) {
        // Compute LtDL-decomposition    | Qa_hat = La_mat' * diag(da_vec) * La_mat
        DecompositionResult ltDl = DecomposeLtDL.decomposeLtDL(qaHat.copy());
        SimpleMatrix laMat = ltDl.getLMat();
        double[] daVec = ltDl.getDVec();

        // Compute Z-transformation      | Reduction & ordering of cond. variances
        TransformResult zTrans = TransformZ.transformZ(laMat, daVec,null);
        SimpleMatrix lzMat = zTrans.getLmat();
        double[] dzVec = zTrans.getdVec();
        SimpleMatrix iZtMat = zTrans.getiZtMat();

        // Apply Z-transformation        | Qa_hat = iZt_mat * Qz_hat * iZt_mat'
        // Qz_hat = Lz_mat' * (dz_vec' .* Lz_mat)
        SimpleMatrix dzVecDiag = new SimpleMatrix(dzVec.length, dzVec.length);
       
        for (int i = 0; i < dzVec.length; i++) {
            dzVecDiag.set(i, i, dzVec[i]);
        }

        SimpleMatrix qzHat = lzMat.transpose().mult(dzVecDiag).mult(lzMat);
        SimpleMatrix zMat = null;
        SimpleMatrix zHat = null;

        // If provided, Z-transform also the ambiguity float vector
        if (aHat != null) {	
            
            
            // Retrieve the Z-transformation matrix
            zMat = iZtMat.invert().transpose();
            // Assuming rounding to nearest integer for unimodular
            zMat = roundMatrix(new SimpleMatrix(zMat));

            // Transform the ambiguity float vector
            zHat = zMat.transpose().mult(aHat);
        }

        return new DecorrelateVCResult(qzHat, lzMat, dzVec, iZtMat, zMat, zHat);
    }

    /**
     * Rounds each element of the matrix to the nearest integer.
     *
     * @param matrix The matrix to be rounded.
     * @return A new matrix with each element rounded.
     */
    private static SimpleMatrix roundMatrix(SimpleMatrix matrix) {
        SimpleMatrix rounded = new SimpleMatrix(matrix.numRows(), matrix.numCols());
        for (int i = 0; i < matrix.numRows(); i++) {
            for (int j = 0; j < matrix.numCols(); j++) {
                rounded.set(i, j, Math.round(matrix.get(i, j)));
            }
        }
        return rounded;
    }


}

/**
 * Class to hold the result of the decorrelateVC method.
 */
class DecorrelateVCResult {
    private SimpleMatrix qzHat;
    private SimpleMatrix lzMat;
    private double[] dzVec;
    private SimpleMatrix iZtMat;
    private SimpleMatrix zMat;
    private SimpleMatrix zHat;

    public DecorrelateVCResult(SimpleMatrix qzHat, SimpleMatrix lzMat, double[] dzVec,
                               SimpleMatrix iZtMat, SimpleMatrix zMat, SimpleMatrix zHat) {
        this.qzHat = qzHat;
        this.lzMat = lzMat;
        this.dzVec = dzVec;
        this.iZtMat = iZtMat;
        this.zMat = zMat;
        this.zHat = zHat;
    }

    public SimpleMatrix getQzHat() {
        return qzHat;
    }

    public SimpleMatrix getLzMat() {
        return lzMat;
    }

    public double[] getDzVec() {
        return dzVec;
    }

    public SimpleMatrix getIZtMat() {
        return iZtMat;
    }

    public SimpleMatrix getZMat() {
        return zMat;
    }

    public SimpleMatrix getZHat() {
        return zHat;
    }
}


