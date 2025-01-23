package com.gnssAug.helper.lambdaNew.models;

import org.apache.commons.math3.linear.*;

public class DecorrelateResult {
    private final RealMatrix Qz_hat;
    private final RealMatrix Lz_mat;
    private final RealVector dz_vec;
    private final RealMatrix iZt_mat;
    private final RealMatrix Z_mat;
    private final RealVector z_hat;

    public DecorrelateResult(RealMatrix Qz_hat, RealMatrix Lz_mat, RealVector dz_vec, RealMatrix iZt_mat, RealMatrix Z_mat, RealVector z_hat) {
        this.Qz_hat = Qz_hat;
        this.Lz_mat = Lz_mat;
        this.dz_vec = dz_vec;
        this.iZt_mat = iZt_mat;
        this.Z_mat = Z_mat;
        this.z_hat = z_hat;
    }

    public RealMatrix getQz_hat() {
        return Qz_hat;
    }

    public RealMatrix getLz_mat() {
        return Lz_mat;
    }

    public RealVector getDz_vec() {
        return dz_vec;
    }

    public RealMatrix getiZt_mat() {
        return iZt_mat;
    }

    public RealMatrix getZ_mat() {
        return Z_mat;
    }

    public RealVector getZ_hat() {
        return z_hat;
    }
}