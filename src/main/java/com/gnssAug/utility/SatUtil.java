package com.gnssAug.utility;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.Set;
import java.util.stream.IntStream;

import org.apache.commons.collections.set.ListOrderedSet;
import org.ejml.simple.SimpleMatrix;

import com.gnssAug.Android.constants.Measurement;
import com.gnssAug.Android.estimation.LinearLeastSquare;
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
	
	public static ArrayList<Satellite> createCopy(ArrayList<Satellite> original) throws CloneNotSupportedException
	{
		int n = original.size();
		ArrayList<Satellite> copy = new ArrayList<Satellite>(n);
		for(int i=0;i<n;i++)
		{
			copy.add(original.get(i).clone());
		}
		return copy;
	}
	
	public static String[] findObsvCodeArray(ArrayList<Satellite> satList) {
		LinkedHashSet<String> obsvCodeSet = new LinkedHashSet<String>();
		for (int i = 0; i < satList.size(); i++) {
			obsvCodeSet.add(satList.get(i).getObsvCode());
		}
		return obsvCodeSet.toArray(new String[0]);

	}
	
	public static ListOrderedSet findSSIset(String[] obsvCodeList)
	{
		ListOrderedSet ssiSet = new ListOrderedSet();
		for(int i=0;i<obsvCodeList.length;i++)
		{
			ssiSet.add(obsvCodeList[i].charAt(0));
		}
		return ssiSet;
	}
	public static LinkedHashSet<String> findObsvCodeSet(ArrayList<Satellite> satList) {
		LinkedHashSet<String> obsvCodeSet = new LinkedHashSet<String>();
		for (int i = 0; i < satList.size(); i++) {
			obsvCodeSet.add(satList.get(i).getObsvCode());
		}
		return obsvCodeSet;

	}
	
	public static Object[] resetVar(Measurement meas, String[] obsvCodeList, double[] estX, SimpleMatrix estXcov)
	{
		int m = obsvCodeList.length;
		ArrayList<Satellite> temp_satList = LinearLeastSquare.getTestedSatList(Measurement.Doppler);
		Set<String> obsvCode_current = SatUtil.findObsvCodeSet(temp_satList);
		ArrayList<Integer> indices = new ArrayList<Integer>();
		if(obsvCode_current.size()!=obsvCodeList.length)
		{
			for(int i =0;i<m;i++)
			{
				if(obsvCode_current.contains(obsvCodeList[i]))
				{
					indices.add(i);
				}
				
			}
		}
		if(indices.size()>0)
		{
			SimpleMatrix _estXcov = new SimpleMatrix(3+m,3+m);
			_estXcov.fill(1e10);
			_estXcov.insertIntoThis(0, 0, estXcov.extractMatrix(0, 3, 0, 3));
			double[] _estX = new double[3+m];
			System.arraycopy(estX, 0, _estX, 
                    0, 3);
			for(int i=0;i<indices.size();i++)
			{
				int index = indices.get(i);
				_estX[3+index] = estX[3+i];
				_estXcov.insertIntoThis(3+index, 0, estXcov.extractMatrix(3+i, 3+1+i, 0, 3));
				_estXcov.insertIntoThis(0,3+index, estXcov.extractMatrix(0, 3, 3+i, 3+1+i));
			}
			int k = 0,l = 0;
			for(int i=0;i<m;i++)
			{
				if(indices.contains(i))
				{
					l=0;
					for(int j=0;j<m;j++)
					{
						if(indices.contains(j))
						{
							_estXcov.set(i+3,j+3,estXcov.get(k+3, l+3));
							l++;
						}
					}
					k++;
				}
			}
			estX= _estX;
			estXcov = _estXcov;
		}
		
		return new Object[] {estX,estXcov};
		
	
	}

}
