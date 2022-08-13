package com.gnssAug.utility;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.stream.IntStream;

import com.gnssAug.Android.models.Satellite;

public class SatUtil {

	public static double[][] getUnitLOS(ArrayList<Satellite> SV, double[] userXYZ) {
		int SVcount = SV.size();
		double[][] unitLOS = new double[SVcount][3];
		for (int k = 0; k < SVcount; k++) {
			Satellite sat = SV.get(k);
			// Line of Sight vector
			double[] LOS = IntStream.range(0, 3).mapToDouble(i -> sat.getSatEci()[i] - userXYZ[i]).toArray();
			double GeometricRange = Math.sqrt(Arrays.stream(LOS).map(i -> i * i).reduce(0.0, (i, j) -> i + j));
			// Converting LOS to unit vector
			unitLOS[k] = Arrays.stream(LOS).map(i -> i / GeometricRange).toArray();
		}
		return unitLOS;

	}

	public static double[] getUnitLOS(double[] satXYZ, double[] userXYZ) {

		double[] unitLOS = new double[3];

		// Line of Sight vector
		double[] LOS = IntStream.range(0, 3).mapToDouble(i -> satXYZ[i] - userXYZ[i]).toArray();
		double GeometricRange = Math.sqrt(Arrays.stream(LOS).map(i -> i * i).reduce(0.0, (i, j) -> i + j));
		// Converting LOS to unit vector
		unitLOS = Arrays.stream(LOS).map(i -> i / GeometricRange).toArray();

		return unitLOS;

	}

}
