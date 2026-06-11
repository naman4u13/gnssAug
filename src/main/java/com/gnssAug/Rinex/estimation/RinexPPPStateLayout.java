package com.gnssAug.Rinex.estimation;

/**
 * State-vector index boundaries for pos_kfObj in EKF_PPP.
 *
 * Layout: pos(3) | clkOff(clkOffNum) | vel(3) | clkDrift(clkDriftNum) | tropo(1) | amb(n) | iono(ionoParamNum)
 */
final class RinexPPPStateLayout {

    final int clkOffNum;
    final int clkDriftNum;
    final int n;
    final int ionoParamNum;

    final int clkOffStart;    // = 3
    final int velStart;       // = 3 + clkOffNum
    final int clkDriftStart;  // = 6 + clkOffNum
    final int tropoIdx;       // = 6 + clkOffNum + clkDriftNum
    final int coreStateNum;   // = tropoIdx + 1
    final int ambStart;       // = coreStateNum
    final int ionoStart;      // = coreStateNum + n
    final int totalStateNum;  // = ionoStart + ionoParamNum

    RinexPPPStateLayout(int clkOffNum, int clkDriftNum, int n, int ionoParamNum) {
        this.clkOffNum    = clkOffNum;
        this.clkDriftNum  = clkDriftNum;
        this.n            = n;
        this.ionoParamNum = ionoParamNum;

        clkOffStart   = 3;
        velStart      = 3 + clkOffNum;
        clkDriftStart = 6 + clkOffNum;
        tropoIdx      = 6 + clkOffNum + clkDriftNum;
        coreStateNum  = tropoIdx + 1;
        ambStart      = coreStateNum;
        ionoStart     = coreStateNum + n;
        totalStateNum = ionoStart + ionoParamNum;
    }
}
