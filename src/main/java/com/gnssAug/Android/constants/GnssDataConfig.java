package com.gnssAug.Android.constants;

public class GnssDataConfig {
	
	final public static double pseudorange_priorVarOfUnitW =1.5379;//Math.pow(1.0299, 2);
	final public static double doppler_priorVarOfUnitW =Math.pow(0.6*0.25, 2);
	final public static double tdcp_priorVarOfUnitW =Math.pow(0.3*0.25, 2);
	final public static double clkDriftVar = 1e5;
	//final public static double[] qENU_posRandWalk = new double[] { 25e-2,25,1 };
	//final public static double[] qENU_posRandWalk = new double[] { 25,25,1 };
	
	//Static
	final public static double[] qENU_posRandWalk = new double[] { 1e-10,1e-10,1e-10};
	
	//final public static double[] qENU_velRandWalk = new double[] { 0.0005, 0.5, 0.01};
	//final public static double[] qENU_velRandWalk = new double[] { 0.5, 0.5, 0.01 };
	
	//Static
	final public static double[] qENU_velRandWalk = new double[] { 1e-10, 1e-10, 1e-10 };
}
