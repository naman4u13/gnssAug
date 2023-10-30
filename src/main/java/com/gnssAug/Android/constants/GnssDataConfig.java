package com.gnssAug.Android.constants;

public class GnssDataConfig {
	
	final public static double pseudorange_priorVarOfUnitW =7.9817;//Math.pow(2.0248456731316584, 2);
	final public static double doppler_priorVarOfUnitW =0.0602;//Math.pow(0.16522711641858306, 2);
	final public static double clkDriftVar = 1e5;
	//final public static double[] qENU_posRandWalk = new double[] { 25e-2,25,1 };
	final public static double[] qENU_posRandWalk = new double[] { 25,25,1 };
	//final public static double[] qENU_velRandWalk = new double[] { 0.0005, 0.5, 0.01};
	final public static double[] qENU_velRandWalk = new double[] { 0.5, 0.5, 0.01 };
}
