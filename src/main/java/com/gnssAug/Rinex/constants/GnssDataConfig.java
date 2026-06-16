package com.gnssAug.Rinex.constants;

// Geodetic receiver noise levels — tighter than Android.constants.GnssDataConfig
// (phase σ = 1 mm vs smartphone 1 cm; pseudorange σ = 0.1 m vs 3.16 m)
public class GnssDataConfig {
	
	final public static double pseudorange_priorVarOfUnitW =1e-2;
	final public static double doppler_priorVarOfUnitW =3e-4;//Math.pow(0.0979, 2);
	final public static double tdcp_priorVarOfUnitW =3e-4;//Math.pow(0.05, 2);
	final public static double phase_priorVarOfUnitW =1e-6;//Math.pow(0.05, 2);
	final public static double GIM_TECU_variance= Math.pow(2,2); // 2 TECU — IGS final GIM accuracy over Europe
	final public static double clkDriftVar = 1e5;
	//final public static double[] qENU_posRandWalk = new double[] { 25e-2,25,1 };
	final public static double[] qENU_posRandWalk = new double[] { 1e-10,1e-10,1e-10};
	//final public static double[] qENU_velRandWalk = new double[] { 0.0005, 0.5, 0.01};
	final public static double[] qENU_velRandWalk = new double[] {  1e-16,1e-16,1e-16 };
	final public static double nSamplesMC = 1e5;
}
