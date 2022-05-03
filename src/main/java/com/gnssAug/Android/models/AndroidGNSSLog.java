package com.gnssAug.Android.models;

import java.util.Map;

public class AndroidGNSSLog {

	private static final long nanosInWeek = (long) (604800 * 1e9);
	private static final Map<Integer, String> freqMap = Map.of(1176450050, "5", 1575420030, "1", 1602000000, "1",
			1561097980, "2");
	private static final Map<Integer, String> SSIMap = Map.of(1, "G", 2, "S", 3, "R", 4, "J", 5, "C", 6, "E", 7, "I");
	private static final Map<Integer, Integer> qzssMap = Map.of(193, 1, 194, 2, 199, 3, 195, 4);

	// Raw,utcTimeMillis,TimeNanos,LeapSecond,TimeUncertaintyNanos,FullBiasNanos,BiasNanos,BiasUncertaintyNanos,DriftNanosPerSecond,DriftUncertaintyNanosPerSecond,HardwareClockDiscontinuityCount,Svid,TimeOffsetNanos,State,ReceivedSvTimeNanos,ReceivedSvTimeUncertaintyNanos,Cn0DbHz,PseudorangeRateMetersPerSecond,PseudorangeRateUncertaintyMetersPerSecond,AccumulatedDeltaRangeState,AccumulatedDeltaRangeMeters,AccumulatedDeltaRangeUncertaintyMeters,CarrierFrequencyHz,CarrierCycles,CarrierPhase,CarrierPhaseUncertainty,MultipathIndicator,SnrInDb,ConstellationType,AgcDb,BasebandCn0DbHz,FullInterSignalBiasNanos,FullInterSignalBiasUncertaintyNanos,SatelliteInterSignalBiasNanos,SatelliteInterSignalBiasUncertaintyNanos,CodeType,ChipsetElapsedRealtimeNanos\n
	// Raw,utcTimeMillis,TimeNanos,LeapSecond,TimeUncertaintyNanos,FullBiasNanos,BiasNanos,BiasUncertaintyNanos,DriftNanosPerSecond,DriftUncertaintyNanosPerSecond,HardwareClockDiscontinuityCount,Svid,TimeOffsetNanos,State,ReceivedSvTimeNanos,ReceivedSvTimeUncertaintyNanos,Cn0DbHz,PseudorangeRateMetersPerSecond,PseudorangeRateUncertaintyMetersPerSecond,AccumulatedDeltaRangeState,AccumulatedDeltaRangeMeters,AccumulatedDeltaRangeUncertaintyMeters,CarrierFrequencyHz,CarrierCycles,CarrierPhase,CarrierPhaseUncertainty,MultipathIndicator,SnrInDb,ConstellationType,AgcDb\n
	private long utcTimeMillis;
	private long timeNanos;
	private int leapSecond;
	private double timeUncertaintyNanos;
	private long fullBiasNanos;
	private double biasNanos;
	private double biasUncertaintyNanos;
	private double driftNanosPerSecond;
	private double driftUncertaintyNanosPerSecond;
	private int hardwareClockDiscontinuityCount;
	private int svid;
	private double timeOffsetNanos;
	private int state;
	private double tTx;
	private double receivedSvTimeUncertaintyNanos;
	private double cn0DbHz;
	private double pseudorangeRateMetersPerSecond;
	private double pseudorangeRateUncertaintyMetersPerSecond;
	private int accumulatedDeltaRangeState;
	private double accumulatedDeltaRangeMeters;
	private double accumulatedDeltaRangeUncertaintyMeters;
	private double carrierFrequencyHz;
	private long carrierCycles;
	private double carrierPhase;
	private double carrierPhaseUncertainty;
	private int multipathIndicator;
	private double snrInDb;
	private int constellationType;
	private double agcDb;
	private double basebandCn0DbHz;
	private double fullInterSignalBiasNanos;
	private double fullInterSignalBiasUncertaintyNanos;
	private double satelliteInterSignalBiasNanos;
	private double satelliteInterSignalBiasUncertaintyNanos;
	private String codeType;
	private long chipsetElapsedRealtimeNanos;
	private String obsvCode;
	private double tRx;
	private long bootGPStime;
	// GPS week;
	private int weekNo;

	public AndroidGNSSLog(AndroidGNSSLog log) {
		super();
		this.utcTimeMillis = log.utcTimeMillis;
		this.timeNanos = log.timeNanos;
		this.leapSecond = log.leapSecond;
		this.timeUncertaintyNanos = log.timeUncertaintyNanos;
		this.fullBiasNanos = log.fullBiasNanos;
		this.biasNanos = log.biasNanos;
		this.biasUncertaintyNanos = log.biasUncertaintyNanos;
		this.driftNanosPerSecond = log.driftNanosPerSecond;
		this.driftUncertaintyNanosPerSecond = log.driftUncertaintyNanosPerSecond;
		this.hardwareClockDiscontinuityCount = log.hardwareClockDiscontinuityCount;
		this.svid = log.svid;
		this.timeOffsetNanos = log.timeOffsetNanos;
		this.state = log.state;
		this.tTx = log.tTx;
		this.receivedSvTimeUncertaintyNanos = log.receivedSvTimeUncertaintyNanos;
		this.cn0DbHz = log.cn0DbHz;
		this.pseudorangeRateMetersPerSecond = log.pseudorangeRateMetersPerSecond;
		this.pseudorangeRateUncertaintyMetersPerSecond = log.pseudorangeRateUncertaintyMetersPerSecond;
		this.accumulatedDeltaRangeState = log.accumulatedDeltaRangeState;
		this.accumulatedDeltaRangeMeters = log.accumulatedDeltaRangeMeters;
		this.accumulatedDeltaRangeUncertaintyMeters = log.accumulatedDeltaRangeUncertaintyMeters;
		this.carrierFrequencyHz = log.carrierFrequencyHz;
		this.carrierCycles = log.carrierCycles;
		this.carrierPhase = log.carrierPhase;
		this.carrierPhaseUncertainty = log.carrierPhaseUncertainty;
		this.multipathIndicator = log.multipathIndicator;
		this.snrInDb = log.snrInDb;
		this.constellationType = log.constellationType;
		this.agcDb = log.agcDb;
		this.basebandCn0DbHz = log.basebandCn0DbHz;
		this.fullInterSignalBiasNanos = log.fullInterSignalBiasNanos;
		this.fullInterSignalBiasUncertaintyNanos = log.fullInterSignalBiasUncertaintyNanos;
		this.satelliteInterSignalBiasNanos = log.satelliteInterSignalBiasNanos;
		this.satelliteInterSignalBiasUncertaintyNanos = log.satelliteInterSignalBiasUncertaintyNanos;
		this.codeType = log.codeType;
		this.chipsetElapsedRealtimeNanos = log.chipsetElapsedRealtimeNanos;
		this.svid = log.svid;
		this.obsvCode = log.obsvCode;
		this.tRx = log.tRx;
		this.weekNo = log.weekNo;
		this.bootGPStime = log.bootGPStime;

	}

	public AndroidGNSSLog(String[] data) {
		super();
		this.utcTimeMillis = Long.parseLong(data[1]);
		this.timeNanos = Long.parseLong(data[2]);
		this.leapSecond = data[3].isBlank() ? 0 : Integer.parseInt(data[3]);
		this.timeUncertaintyNanos = data[4].isBlank() ? 0 : Double.parseDouble(data[4]);
		this.fullBiasNanos = Long.parseLong(data[5]);
		this.biasNanos = Double.parseDouble(data[6]);
		this.biasUncertaintyNanos = Double.parseDouble(data[7]);
		this.driftNanosPerSecond = Double.parseDouble(data[8]);
		this.driftUncertaintyNanosPerSecond = Double.parseDouble(data[9]);
		this.hardwareClockDiscontinuityCount = Integer.parseInt(data[10]);
		this.svid = Integer.parseInt(data[11]);
		this.timeOffsetNanos = Double.parseDouble(data[12]);
		this.state = Integer.parseInt(data[13]);
		// tTx
		this.tTx = Long.parseLong(data[14]) / 1e9;
		this.receivedSvTimeUncertaintyNanos = Double.parseDouble(data[15]);
		this.cn0DbHz = Double.parseDouble(data[16]);
		this.pseudorangeRateMetersPerSecond = Double.parseDouble(data[17]);
		this.pseudorangeRateUncertaintyMetersPerSecond = Double.parseDouble(data[18]);
		this.accumulatedDeltaRangeState = Integer.parseInt(data[19]);
		this.accumulatedDeltaRangeMeters = Double.parseDouble(data[20]);
		this.accumulatedDeltaRangeUncertaintyMeters = Double.parseDouble(data[21]);
		this.carrierFrequencyHz = Double.parseDouble(data[22]);
		this.carrierCycles = data[23].isBlank() ? 0 : Long.parseLong(data[23]);
		this.carrierPhase = data[24].isBlank() ? 0 : Double.parseDouble(data[24]);
		this.carrierPhaseUncertainty = data[25].isBlank() ? 0 : Double.parseDouble(data[25]);
		this.multipathIndicator = Integer.parseInt(data[26]);
		this.snrInDb = data[27].isBlank() ? 0 : Double.parseDouble(data[27]);
		this.constellationType = Integer.parseInt(data[28]);
		this.agcDb = data[29].isBlank() ? 0 : Double.parseDouble(data[29]);
		if (data.length > 30) {
			this.basebandCn0DbHz = data[30].isBlank() ? 0 : Double.parseDouble(data[30]);
			this.fullInterSignalBiasNanos = data[31].isBlank() ? 0 : Double.parseDouble(data[31]);
			this.fullInterSignalBiasUncertaintyNanos = data[32].isBlank() ? 0 : Double.parseDouble(data[32]);
			this.satelliteInterSignalBiasNanos = data[33].isBlank() ? 0 : Double.parseDouble(data[33]);
			this.satelliteInterSignalBiasUncertaintyNanos = data[34].isBlank() ? 0 : Double.parseDouble(data[34]);
			this.codeType = data[35];
			this.chipsetElapsedRealtimeNanos = Long.parseLong(data[36]);
		}
		int freq = (int) carrierFrequencyHz;
		if (freq >= 1598062460 && freq <= 1608750000) {
			freq = (int) 1602e6;
		}

		String freqID = freqMap.get(freq);
		String channel = freqID.equals("1") ? "C" : freqID.equals("2") ? "I" : "X";

		String ssi = SSIMap.get(constellationType);
		if (ssi.equals("J")) {
			svid = qzssMap.get(svid);
		}
		obsvCode = ssi + freqID + channel;
		// Only integer value, haven't subtracted the fractional value
		long nanosSinceGpsEpoch = timeNanos - fullBiasNanos;
		weekNo = Math.round(nanosSinceGpsEpoch / nanosInWeek);
		tRx = (nanosSinceGpsEpoch % nanosInWeek) + timeOffsetNanos - biasNanos;
		bootGPStime = Math.round(tRx - chipsetElapsedRealtimeNanos);
//		long num = 40772809490L;
//		if (((long) (bootGPStime / 1e4)) != num) {
//			System.out.println();
//		}
		tRx /= 1e9;

	}

	public long getUtcTimeMillis() {
		return utcTimeMillis;
	}

	public long getTimeNanos() {
		return timeNanos;
	}

	public int getLeapSecond() {
		return leapSecond;
	}

	public double getTimeUncertaintyNanos() {
		return timeUncertaintyNanos;
	}

	public long getFullBiasNanos() {
		return fullBiasNanos;
	}

	public double getBiasNanos() {
		return biasNanos;
	}

	public double getBiasUncertaintyNanos() {
		return biasUncertaintyNanos;
	}

	public double getDriftNanosPerSecond() {
		return driftNanosPerSecond;
	}

	public double getDriftUncertaintyNanosPerSecond() {
		return driftUncertaintyNanosPerSecond;
	}

	public int getHardwareClockDiscontinuityCount() {
		return hardwareClockDiscontinuityCount;
	}

	public int getSvid() {
		return svid;
	}

	public double getTimeOffsetNanos() {
		return timeOffsetNanos;
	}

	public int getState() {
		return state;
	}

	public double gettTx() {
		return tTx;
	}

	public double getReceivedSvTimeUncertaintyNanos() {
		return receivedSvTimeUncertaintyNanos;
	}

	public double getCn0DbHz() {
		return cn0DbHz;
	}

	public double getPseudorangeRateMetersPerSecond() {
		return pseudorangeRateMetersPerSecond;
	}

	public double getPseudorangeRateUncertaintyMetersPerSecond() {
		return pseudorangeRateUncertaintyMetersPerSecond;
	}

	public int getAccumulatedDeltaRangeState() {
		return accumulatedDeltaRangeState;
	}

	public double getAccumulatedDeltaRangeMeters() {
		return accumulatedDeltaRangeMeters;
	}

	public double getAccumulatedDeltaRangeUncertaintyMeters() {
		return accumulatedDeltaRangeUncertaintyMeters;
	}

	public double getCarrierFrequencyHz() {
		return carrierFrequencyHz;
	}

	public long getCarrierCycles() {
		return carrierCycles;
	}

	public double getCarrierPhase() {
		return carrierPhase;
	}

	public double getCarrierPhaseUncertainty() {
		return carrierPhaseUncertainty;
	}

	public int getMultipathIndicator() {
		return multipathIndicator;
	}

	public double getSnrInDb() {
		return snrInDb;
	}

	public double getConstellationType() {
		return constellationType;
	}

	public double getAgcDb() {
		return agcDb;
	}

	public double getBasebandCn0DbHz() {
		return basebandCn0DbHz;
	}

	public double getFullInterSignalBiasNanos() {
		return fullInterSignalBiasNanos;
	}

	public double getFullInterSignalBiasUncertaintyNanos() {
		return fullInterSignalBiasUncertaintyNanos;
	}

	public double getSatelliteInterSignalBiasNanos() {
		return satelliteInterSignalBiasNanos;
	}

	public double getSatelliteInterSignalBiasUncertaintyNanos() {
		return satelliteInterSignalBiasUncertaintyNanos;
	}

	public String getCodeType() {
		return codeType;
	}

	public long getChipsetElapsedRealtimeNanos() {
		return chipsetElapsedRealtimeNanos;
	}

	public String getObsvCode() {
		return obsvCode;
	}

	public double gettRx() {
		return tRx;
	}

	@Override
	public String toString() {
		return utcTimeMillis + ", " + timeNanos + ", " + leapSecond + ", " + timeUncertaintyNanos + ", " + fullBiasNanos
				+ ", " + biasNanos + ", " + biasUncertaintyNanos + ", " + driftNanosPerSecond + ", "
				+ driftUncertaintyNanosPerSecond + ", " + hardwareClockDiscontinuityCount + ", " + svid + ", "
				+ timeOffsetNanos + ", " + state + ", " + tTx + ", " + receivedSvTimeUncertaintyNanos + ", " + cn0DbHz
				+ ", " + pseudorangeRateMetersPerSecond + ", " + pseudorangeRateUncertaintyMetersPerSecond + ", "
				+ accumulatedDeltaRangeState + ", " + accumulatedDeltaRangeMeters + ", "
				+ accumulatedDeltaRangeUncertaintyMeters + ", " + carrierFrequencyHz + ", " + carrierCycles + ", "
				+ carrierPhase + ", " + carrierPhaseUncertainty + ", " + multipathIndicator + ", " + snrInDb + ", "
				+ constellationType + ", " + agcDb + ", " + basebandCn0DbHz + ", " + fullInterSignalBiasNanos + ", "
				+ fullInterSignalBiasUncertaintyNanos + ", " + satelliteInterSignalBiasNanos + ", "
				+ satelliteInterSignalBiasUncertaintyNanos + ", " + codeType + ", " + chipsetElapsedRealtimeNanos + ", "
				+ obsvCode + ", " + tRx + ", " + bootGPStime;
	}

	public long getBootGPStime() {
		return bootGPStime;
	}

	public int getWeekNo() {
		return weekNo;
	}

}
