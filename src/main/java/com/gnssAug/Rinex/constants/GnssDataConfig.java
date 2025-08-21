package com.gnssAug.Rinex.constants;

public class GnssDataConfig {
	
	final public static double pseudorange_priorVarOfUnitW =3e-2;
	final public static double doppler_priorVarOfUnitW =5e-4;//Math.pow(0.0979, 2);
	final public static double tdcp_priorVarOfUnitW =5e-4;//Math.pow(0.05, 2);
	final public static double phase_priorVarOfUnitW =3e-4;//Math.pow(0.05, 2);
	final public static double GIM_TECU_variance= Math.pow(6,2);//0.0006;
	final public static double clkDriftVar = 1e5;
	//final public static double[] qENU_posRandWalk = new double[] { 25e-2,25,1 };
	final public static double[] qENU_posRandWalk = new double[] { 1e-10,1e-10,1e-10};
	//final public static double[] qENU_velRandWalk = new double[] { 0.0005, 0.5, 0.01};
	final public static double[] qENU_velRandWalk = new double[] {  1e-10,1e-10,1e-10 };
}
