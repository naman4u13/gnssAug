package com.gnssAug.Android.constants;

public class GnssDataConfig {
	
	final public static double pseudorange_priorVarOfUnitW =Math.pow(1.0299, 2);
	final public static double doppler_priorVarOfUnitW =0.008;//Math.pow(0.0979, 2);
	final public static double tdcp_priorVarOfUnitW =0.007;//Math.pow(0.05, 2);
	final public static double clkDriftVar = 1e5;
	//final public static double[] qENU_posRandWalk = new double[] { 25e-2,25,1 };
	final public static double[] qENU_posRandWalk = new double[] { 25,25,1 };
	//final public static double[] qENU_velRandWalk = new double[] { 0.0005, 0.5, 0.01};
	final public static double[] qENU_velRandWalk = new double[] { 0.5, 0.5, 0.01 };
}
