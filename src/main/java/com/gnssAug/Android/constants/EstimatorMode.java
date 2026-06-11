package com.gnssAug.Android.constants;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public enum EstimatorMode {

    // --- Code-based positioning ---
    LLS_CODE(1),
    WLS_CODE(2),
    LLS_WLS_BOTH(3),

    // --- INS fusion ---
    GNSS_INS(4),

    // --- Basic EKF ---
    EKF_POS_RW(5),
    EKF_VEL_RW(6),
    EKF_VEL_RW_DOPPLER(7),
    EKF_VEL_RW_ESTVEL(12),
    EKF_VEL_RW_ESTVEL_COMP(13),
    EKF_BASIC_MULTI(9),

    // --- DBP / AKF - thesis Pillar 1 ---
    EKF_DBP(10),
    AKF_DBP(14),

    // --- TDCP velocity ---
    LLS_TDCP_VEL(16),
    TDCP_CSDR(17),

    // --- EKF TDCP - thesis Pillar 2 ---
    EKF_TDCP(18),
    EKF_TDCP_PHASE_RATE(19),
    EKF_TDCP_DOPPLER_ONLY(20),
    EKF_TDCP_ALL_ESTIMATORS(21),

    // --- PPP - thesis Pillar 3 ---
    PPP_FLOAT(22),

    // --- Multi-mode analysis ---
    ALL_ANALYSIS(11),

    // --- RINEX-specific modes (Android_Static_Rinex) ---
    RINEX_EKF_CODE(104),
    RINEX_PPP_LOWCOST(105),
    RINEX_PPP_DF(106),

    // --- IGS station PPP modes ---
    IGS_PPP_FLOAT(201),
    IGS_PPP_AR(202);

    public final int legacyCode;

    EstimatorMode(int legacyCode) {
        this.legacyCode = legacyCode;
    }

    public boolean isLLSMode() {
        return this == LLS_CODE || this == WLS_CODE || this == LLS_WLS_BOTH;
    }

    public boolean isBasicEKFMode() {
        return this == EKF_POS_RW || this == EKF_VEL_RW || this == EKF_VEL_RW_DOPPLER
            || this == EKF_VEL_RW_ESTVEL || this == EKF_VEL_RW_ESTVEL_COMP || this == EKF_BASIC_MULTI;
    }

    public boolean isDBPMode() {
        return this == EKF_DBP || this == AKF_DBP;
    }

    public boolean isTDCPMode() {
        return this == LLS_TDCP_VEL || this == TDCP_CSDR || this == EKF_TDCP
            || this == EKF_TDCP_PHASE_RATE || this == EKF_TDCP_DOPPLER_ONLY
            || this == EKF_TDCP_ALL_ESTIMATORS;
    }

    public boolean isPPPMode() {
        return this == PPP_FLOAT || this == RINEX_PPP_LOWCOST || this == RINEX_PPP_DF
            || this == IGS_PPP_FLOAT || this == IGS_PPP_AR;
    }

    public boolean isAnalysisMode() {
        return this == ALL_ANALYSIS;
    }

    public boolean isTDCPorPPPMode() {
        return isTDCPMode() || isPPPMode();
    }

    public List<EstimatorMode> getSubModes() {
        switch (this) {
            case LLS_WLS_BOTH:    return Arrays.asList(LLS_CODE, WLS_CODE);
            case EKF_BASIC_MULTI: return Arrays.asList(EKF_POS_RW, EKF_VEL_RW, EKF_VEL_RW_ESTVEL);
            case ALL_ANALYSIS:    return Arrays.asList(LLS_CODE, WLS_CODE,
                                      EKF_POS_RW, EKF_VEL_RW, EKF_VEL_RW_ESTVEL,
                                      EKF_DBP, AKF_DBP);
            default:              return Collections.singletonList(this);
        }
    }

    public static EstimatorMode fromLegacyCode(int code) {
        for (EstimatorMode m : values()) {
            if (m.legacyCode == code) return m;
        }
        throw new IllegalArgumentException("Unknown estimatorType: " + code);
    }
}
