package com.gnssAug.helper.lambdaNew.models;
import org.apache.commons.math3.linear.*;

public class TransformZResult {
    private final RealMatrix L;
    private final RealVector D;
    private final RealMatrix iZt;

    public TransformZResult(RealMatrix L, RealVector D, RealMatrix iZt) {
        this.L = L;
        this.D = D;
        this.iZt = iZt;
    }

    public RealMatrix getL() {
        return L;
    }

    public RealVector getD() {
        return D;
    }

    public RealMatrix getIZt() {
        return iZt;
    }
}