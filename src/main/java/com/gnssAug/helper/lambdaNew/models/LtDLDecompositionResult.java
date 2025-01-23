package com.gnssAug.helper.lambdaNew.models;
import org.apache.commons.math3.linear.*;

public class LtDLDecompositionResult {
    private final RealMatrix L;
    private final RealVector D;

    public LtDLDecompositionResult(RealMatrix L, RealVector D) {
        this.L = L;
        this.D = D;
    }

    public RealMatrix getL() {
        return L;
    }

    public RealVector getD() {
        return D;
    }
}