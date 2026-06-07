package com.gnssAug.Android.estimation.KalmanFilter;

/**
 * Encodes the state-vector index boundaries for the pos_kfObj in EKF_PPP3,
 * eliminating scattered arithmetic like {@code 3 + codeClkOffNum} or
 * {@code coreStateNum - 1} throughout the filter methods.
 *
 * State vector (pos_kfObj):
 *   [0:2]                                  ECEF position (3)
 *   [codeClkStart : phaseClkStart-1]       Code clock offsets (codeClkOffNum)
 *   [phaseClkStart : velStart-1]           Phase clock offsets (phaseClkOffNum)
 *   [velStart : clkDriftStart-1]           Velocity (3)
 *   [clkDriftStart : tropoIdx-1]           Clock drifts (clkDriftNum)
 *   [tropoIdx]                             Tropospheric zenith wet delay (1)
 *   [ambStart : ionoStart-1]              Float ambiguities (n, per epoch)
 *   [ionoStart : totalStateNum-1]          Ionospheric delays (ionoParamNum, per epoch)
 */
final class PPPStateLayout {

    final int codeClkOffNum;
    final int phaseClkOffNum;
    final int clkDriftNum;
    final int n;
    final int ionoParamNum;

    final int codeClkStart;   // = 3
    final int phaseClkStart;  // = 3 + codeClkOffNum
    final int velStart;       // = 3 + codeClkOffNum + phaseClkOffNum
    final int clkDriftStart;  // = 6 + codeClkOffNum + phaseClkOffNum
    final int tropoIdx;       // = 6 + codeClkOffNum + phaseClkOffNum + clkDriftNum
    final int coreStateNum;   // = tropoIdx + 1
    final int ambStart;       // = coreStateNum
    final int ionoStart;      // = coreStateNum + n
    final int totalStateNum;  // = ionoStart + ionoParamNum

    PPPStateLayout(int codeClkOffNum, int phaseClkOffNum, int clkDriftNum, int n, int ionoParamNum) {
        this.codeClkOffNum  = codeClkOffNum;
        this.phaseClkOffNum = phaseClkOffNum;
        this.clkDriftNum    = clkDriftNum;
        this.n              = n;
        this.ionoParamNum   = ionoParamNum;

        codeClkStart  = 3;
        phaseClkStart = 3 + codeClkOffNum;
        velStart      = 3 + codeClkOffNum + phaseClkOffNum;
        clkDriftStart = 6 + codeClkOffNum + phaseClkOffNum;
        tropoIdx      = 6 + codeClkOffNum + phaseClkOffNum + clkDriftNum;
        coreStateNum  = tropoIdx + 1;
        ambStart      = coreStateNum;
        ionoStart     = coreStateNum + n;
        totalStateNum = ionoStart + ionoParamNum;
    }
}
