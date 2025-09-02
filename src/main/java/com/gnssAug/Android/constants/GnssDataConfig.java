package com.gnssAug.Android.constants;

public class GnssDataConfig {
	
	final public static double pseudorange_priorVarOfUnitW =  100;//Math.pow(2.47, 2);
	final public static double doppler_priorVarOfUnitW =0.001;//Math.pow(0.0616, 2);
	final public static double tdcp_priorVarOfUnitW = 0.001;//0.0006;
	final public static double phase_priorVarOfUnitW = 1;//0.0006;
	final public static double GIM_TECU_variance= Math.pow(3,2);//0.0006;
	final public static double clkDriftVar = 1e5;
	//final public static double[] qENU_posRandWalk = new double[] { 25e-2,25,1 };
//	final public static double[] qENU_posRandWalk = new double[] { 25,25,1 };
	
	//Static
	final public static double[] qENU_posRandWalk = new double[] { 1e-10,1e-10,1e-10};
	
	//final public static double[] qENU_velRandWalk = new double[] { 0.0005, 0.5, 0.01};
//	final public static double[] qENU_velRandWalk = new double[] { 0.5, 0.5, 0.01 };
	
	//Static
	final public static double[] qENU_velRandWalk = new double[] {  1e-16,1e-16,1e-16 };
	
	final public static double nSamplesMC = 1e5;
}
