package com.gnssAug.helper;

import java.util.ArrayList;
import java.util.stream.IntStream;

import org.ejml.simple.SimpleMatrix;

import com.gnssAug.utility.Matrix;

public class IntegerLeastSquares {
	
	private SimpleMatrix floatAmb;
	private SimpleMatrix floatCov;
	private SimpleMatrix intAmb;
	private SimpleMatrix intAmb2;
	private double minVal = Double.MAX_VALUE;
	private double minVal2 = Double.MAX_VALUE;
	public IntegerLeastSquares(SimpleMatrix floatAmb, SimpleMatrix floatCov) {
		super();
		this.floatAmb = floatAmb;
		this.floatCov = floatCov;
	}

	public boolean process()
	{
		int n = floatAmb.getNumElements();
		ArrayList<int[]> searchSpace = new ArrayList<int[]>();
		for(int i=0;i<n;i++)
		{
			double floatSD = Math.sqrt(floatCov.get(i, i));
			int[] ss = IntStream.range((int)Math.floor((-3)+floatAmb.get(i)),(int)Math.ceil((3)+floatAmb.get(i))+1).toArray();
			searchSpace.add(ss);
		}
		bruteForceSearch(searchSpace, 0, n, new ArrayList<Integer>());
		double ratioTest = minVal2/minVal;
		if(ratioTest>3)
		{
			return true;
		}
		return false;
		
	}
	
	private void bruteForceSearch(ArrayList<int[]> searchSpace,int i,int n,ArrayList<Integer> intAmbList)
	{
		if(i<n)
		{
			int[] ss = searchSpace.get(i);
			for(int j=0;j<ss.length;j++)
			{
				int ele = ss[j];
				intAmbList.add(ele);
				bruteForceSearch(searchSpace,i+1,n,intAmbList);
				intAmbList.remove(intAmbList.size()-1);
			}
			
		}
		else
		{
			SimpleMatrix intAmbVec = Matrix.ArrayList2Vector(intAmbList);
			double val = ((floatAmb.minus(intAmbVec).transpose()).mult((floatCov.invert())).mult(floatAmb.minus(intAmbVec))).get(0);
			if(val<minVal)
			{
				minVal2 = minVal;
				minVal = val;
				intAmb2 = intAmb;
				intAmb = new SimpleMatrix(intAmbVec);
			}
			
		}
		
	}

	public SimpleMatrix getIntAmb() {
		return intAmb;
	}
	
	
}
