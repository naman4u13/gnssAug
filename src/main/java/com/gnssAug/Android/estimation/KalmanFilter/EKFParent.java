package com.gnssAug.Android.estimation.KalmanFilter;

import java.util.ArrayList;
import java.util.TreeMap;

import org.ejml.simple.SimpleMatrix;

import com.gnssAug.Android.estimation.KalmanFilter.Models.KFconfig;
import com.gnssAug.Android.models.Satellite;

public class EKFParent {
	protected final double SpeedofLight = 299792458;
	protected KFconfig kfObj;

	protected double[] innovation;
	protected double[] temp_innovation;
	protected TreeMap<Long, double[]> innovationMap;
	protected TreeMap<Long, double[]> residualMap;
	protected TreeMap<Long, double[]> measNoiseMap;
	protected TreeMap<Long, Double> postVarOfUnitWMap;
	// Posteriori Err Cov
	protected TreeMap<Long, SimpleMatrix> errCovMap;
	protected ArrayList<double[]> redundancyList;
	// Satellite Count
	protected TreeMap<Long, Long> satCountMap;
	protected TreeMap<Long, ArrayList<Satellite>> satListMap;

	
	
	
	public TreeMap<Long, double[]> getInnovationMap() {
		return innovationMap;
	}

	public TreeMap<Long, SimpleMatrix> getErrCovMap() {
		return errCovMap;
	}

	public TreeMap<Long, double[]> getResidualMap() {
		return residualMap;
	}

	public TreeMap<Long, Double> getPostVarOfUnitWMap() {
		return postVarOfUnitWMap;
	}

	public ArrayList<double[]> getRedundancyList() {
		return redundancyList;
	}

	public TreeMap<Long, Long> getSatCountMap() {
		return satCountMap;
	}

	public TreeMap<Long, ArrayList<Satellite>> getSatListMap() {
		return satListMap;
	}

	public TreeMap<Long, double[]> getMeasNoiseMap() {
		return measNoiseMap;
	}

}
