package com.gnssAug.Rinex.estimation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.TreeMap;
import java.util.stream.IntStream;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.ejml.dense.row.MatrixFeatures_DDRM;
import org.ejml.simple.SimpleMatrix;

import com.gnssAug.Android.estimation.KalmanFilter.Models.KFconfig;
import com.gnssAug.Rinex.models.Satellite;
import com.gnssAug.utility.LatLonUtil;
import com.gnssAug.utility.Matrix;
import com.gnssAug.utility.Weight;

public class EKF {

	private final double SpeedofLight = 299792458;
	private KFconfig kfObj;
	private double prObsNoiseVar;
	private double[] innovation;
	private TreeMap<Long, double[]> innovationMap;
	private TreeMap<Long, double[]> residualMap;
	private TreeMap<Long, double[]> measNoiseMap;
	private TreeMap<Long, Double> postVarOfUnitWMap;
	// Posteriori Err Cov
	private TreeMap<Long, SimpleMatrix> errCovMap;
	private ArrayList<double[]> redundancyList;
	// Satellite Count
	private TreeMap<Long, Long> satCountMap;
	private TreeMap<Long, ArrayList<Satellite>> satListMap;

	public EKF() {
		kfObj = new KFconfig();
	}

	public TreeMap<Long, double[]> process(TreeMap<Long, ArrayList<Satellite>> SatMap, HashMap<String, double[]> rxPCO,
			ArrayList<Long> timeList, boolean doAnalyze, boolean doTest,boolean outlierAnalyze, String[] obsvCodeList) throws Exception {

		int m = obsvCodeList.length;
		int n = 3 + (2 * m);

		/*
		 * state XYZ is intialized WLS generated Rx position estimated using first epoch
		 * data, a-priori estimate error covariance matrix for state XYZ is therefore
		 * assigned 25 m^2 value. Other state variables are assigned infinite(big)
		 * variance
		 */
		double[][] x = new double[n][1];
		double[][] P = new double[n][n];
		double[] intialECEF = LinearLeastSquare.getEstPos(SatMap.firstEntry().getValue(), rxPCO, true);
		IntStream.range(0, 3 + m).forEach(i -> x[i][0] = intialECEF[i]);

		IntStream.range(0, 3 + m).forEach(i -> P[i][i] = 10);
		IntStream.range(3 + m, 3 + (2 * m)).forEach(i -> P[i][i] = 1e5);

		kfObj.setState_ProcessCov(x, P);
		if (doAnalyze) {
			innovationMap = new TreeMap<Long, double[]>();
			errCovMap = new TreeMap<Long, SimpleMatrix>();
			residualMap = new TreeMap<Long, double[]>();
			postVarOfUnitWMap = new TreeMap<Long, Double>();
			redundancyList = new ArrayList<double[]>();
			satCountMap = new TreeMap<Long, Long>();
			satListMap = new TreeMap<Long, ArrayList<Satellite>>();
			measNoiseMap = new TreeMap<Long, double[]>();
		}
		// Begin iteration or recursion
		return iterate(SatMap, rxPCO, timeList, doAnalyze, doTest,outlierAnalyze, obsvCodeList);

	}

	private TreeMap<Long, double[]> iterate(TreeMap<Long, ArrayList<Satellite>> SatMap, HashMap<String, double[]> rxPCO,
			ArrayList<Long> timeList, boolean doAnalyze, boolean doTest,boolean outlierAnalyze, String[] obsvCodeList) throws Exception {
		TreeMap<Long, double[]> estStateMap = new TreeMap<Long, double[]>();
		int m = obsvCodeList.length;
		int x_size = 3 + (2 * m);
		long time = timeList.get(0);
		// Start from 2nd epoch
		for (int i = 1; i < timeList.size(); i++) {
			long currentTime = timeList.get(i);
			ArrayList<Satellite> satList = SatMap.get(currentTime);
			double deltaT = (currentTime - time) / 1e3;
			// Perform Predict and Update
			runFilter(deltaT, satList, rxPCO, currentTime, doAnalyze, doTest,outlierAnalyze, obsvCodeList );
			// Fetch Posteriori state estimate and estimate error covariance matrix
			SimpleMatrix x = kfObj.getState();
			SimpleMatrix P = kfObj.getCovariance();

			double[] estState = new double[x_size];
			IntStream.range(0, x_size).forEach(j -> estState[j] = x.get(j));

			// Add position estimate to the list
			estStateMap.put(currentTime, estState);
			if (doAnalyze) {
				// Convert to ENU frame
				SimpleMatrix R = new SimpleMatrix(x_size, x_size);
				R.insertIntoThis(0, 0, new SimpleMatrix(LatLonUtil.getEcef2EnuRotMat(estState)));
				IntStream.range(3, x_size).forEach(j -> R.set(j, j, 1));

				SimpleMatrix errCov = R.mult(P).mult(R.transpose());
				errCovMap.put(currentTime, errCov);
				innovationMap.put(currentTime, innovation);
			}
			/*
			 * Check whether estimate error covariance matrix is positive semidefinite
			 * before further proceeding
			 */
			if (!MatrixFeatures_DDRM.isPositiveSemidefinite(P.getMatrix())) {

				throw new Exception("PositiveDefinite test Failed");
			}
			time = currentTime;

		}

		return estStateMap;
	}

	private void runFilter(double deltaT, ArrayList<Satellite> satList, HashMap<String, double[]> rxPCO,
			long currentTime, boolean doAnalyze, boolean doTest,boolean outlierAnalyze, String[] obsvCodeList) throws Exception {

		// Satellite count
		int n = satList.size();
		int m = obsvCodeList.length;
		// Last update VC-matrix
		SimpleMatrix priorP = new SimpleMatrix(kfObj.getCovariance());
		SimpleMatrix priorX = new SimpleMatrix(kfObj.getState());
		// Assign Q and F matrix
		kfObj.configIGS(deltaT, m);
		kfObj.predict();
		
		boolean isWeighted = false;
		SimpleMatrix x = kfObj.getState();
		final double[] estPos = new double[] { x.get(0), x.get(1), x.get(2) };
		double[] rxClkOff = new double[m];// in meters
		for (int i = 0; i < m; i++) {
			rxClkOff[i] = x.get(i + 3);
		}

		/*
		 * H is the Jacobian matrix of partial derivatives Observation StateModel(h) of
		 * with respect to x
		 */
		SimpleMatrix H = getJacobian(satList, estPos, obsvCodeList);

		// Measurement vector
		double[][] z = new double[n][1];
		// Estimated Measurement vector
		double[][] ze = new double[n][1];
		// Measurement Noise
		double[][] _R = new double[n][n];
		innovation = new double[n];
		for (int i = 0; i < n; i++) {

			Satellite sat = satList.get(i);
			String obsvCode = sat.getObsvCode();
			double[] pco = rxPCO.get(obsvCode);
			double[] rxAPC = IntStream.range(0, 3).mapToDouble(j -> estPos[j] + pco[j]).toArray();
			z[i][0] = sat.getPseudorange();
			ze[i][0] = Math.sqrt(IntStream.range(0, 3).mapToDouble(j -> rxAPC[j] - sat.getSatEci()[j]).map(j -> j * j)
					.reduce(0, (j, k) -> j + k));
			for (int j = 0; j < m; j++) {
				if (obsvCode.equals(obsvCodeList[j])) {
					ze[i][0] += rxClkOff[j];
				}
			}

			innovation[i] = z[i][0] - ze[i][0];
		}
		double priorVarOfUnitW = 0.146;
		if (isWeighted) {
			double[][] weight = Weight.computeCovInvMat(satList);
			SimpleMatrix Cyy = null;
			priorVarOfUnitW = 0.03;
			double[][] cov = new double[n][n];
			double max = Double.MIN_VALUE;
			for (int i = 0; i < n; i++) {
				max = Math.max(weight[i][i], max);
			}
			for (int i = 0; i < n; i++) {
				weight[i][i] = weight[i][i] / max;
			}

			for (int i = 0; i < n; i++) {
				cov[i][i] = priorVarOfUnitW / weight[i][i];
			}
			Cyy = new SimpleMatrix(cov);
			for (int i = 0; i < n; i++) {
				_R[i][i] = Cyy.get(i, i);
			}
		} else {
			for (int i = 0; i < n; i++) {
				_R[i][i] = priorVarOfUnitW;
			}
		}
		SimpleMatrix R = new SimpleMatrix(_R);
		ArrayList<Satellite> testedSatList = new ArrayList<Satellite>(satList);	
		if(doTest&&!outlierAnalyze)
		{
			performTesting(R, H, n, m, satList, testedSatList, z, ze);
		}
		kfObj.update(z, R, ze, H);
		if(doAnalyze)
		{
			performAnalysis(testedSatList, rxPCO, R, H, priorP, currentTime, doTest, priorVarOfUnitW, n,obsvCodeList);
		}
		if(outlierAnalyze)
		{
			kfObj.setState_ProcessCov(priorX, priorP);
			kfObj.predict();
			performTesting(R, H, n, m, satList, testedSatList, z, ze);
			kfObj.update(z, R, ze, H);
		}
		
		
		
	}

	private void performAnalysis(ArrayList<Satellite> testedSatList,HashMap<String, double[]> rxPCO,SimpleMatrix R,SimpleMatrix H,SimpleMatrix priorP,long currentTime,boolean doTest,double priorVarOfUnitW,int n,String[] obsvCodeList)
	{
		
			int _n = testedSatList.size();
			int m = obsvCodeList.length;
			double[] residual = new double[_n];
			double[] measNoise = new double[_n];
			SimpleMatrix x = kfObj.getState();
			final double[] _estPos = new double[] { x.get(0), x.get(1), x.get(2) };
			double[] rxClkOff = new double[m];// in meters
			for (int i = 0; i < m; i++) {
				rxClkOff[i] = x.get(i + 3);
			}
			double[][] ze = new double[_n][1];
			double[][] z = new double[_n][1];
			for (int i = 0; i < _n; i++) {

				Satellite sat = testedSatList.get(i);
				String obsvCode = sat.getSSI() + "" + sat.getFreqID() + "C";
				double[] pco = rxPCO.get(obsvCode);
				double[] rxAPC = IntStream.range(0, 3).mapToDouble(j -> _estPos[j] + pco[j]).toArray();
				z[i][0] = sat.getPseudorange();
				ze[i][0] = Math.sqrt(IntStream.range(0, 3).mapToDouble(j -> rxAPC[j] - sat.getSatEci()[j])
						.map(j -> j * j).reduce(0, (j, k) -> j + k));
				for (int j = 0; j < m; j++) {
					if (obsvCode.equals(obsvCodeList[j])) {
						ze[i][0] += rxClkOff[j];
					}
				}
				residual[i] = z[i][0] - ze[i][0];
				
			}
			// Post-fit residual
			SimpleMatrix e_post_hat = new SimpleMatrix(_n, 1, true, residual);
			SimpleMatrix Cyy_inv = R.invert();

			// Compute Redundancies
			SimpleMatrix K = kfObj.getKalmanGain();

			SimpleMatrix HK = H.mult(K);
			SimpleMatrix phi = kfObj.getPhi();
			SimpleMatrix Cvv = kfObj.getCvv();
			SimpleMatrix HtCvvInvH = H.transpose().mult(Cvv.invert()).mult(H);
			SimpleMatrix Q = kfObj.getQ();
			double rX = phi.mult(priorP).mult(phi.transpose()).mult(HtCvvInvH).trace();
			double rW = Q.mult(HtCvvInvH).trace();
			double rZ = SimpleMatrix.identity(HK.numRows()).minus(HK).trace();
			double rSum = rX + rW + rZ;
			if (_n - rSum > 0.01) {
				System.err.println("FATAL ERROR: Redundancy sum is wrong");
			}
			double postVarOfUnitW = 0;
			try {
			postVarOfUnitW = e_post_hat.transpose().mult(Cyy_inv).mult(e_post_hat).get(0) * priorVarOfUnitW / rZ;
			}
			catch (Exception e) {
				// TODO: handle exception
				System.out.print("");
			}
			redundancyList.add(new double[] { _n, rSum, rX, rW, rZ });
			postVarOfUnitWMap.put(currentTime, postVarOfUnitW);
			residualMap.put(currentTime, residual);
			if (doTest == true) {
				_n = n - _n;
			}
			satCountMap.put(currentTime, (long) _n);
			satListMap.put(currentTime, testedSatList);
			measNoiseMap.put(currentTime, measNoise);
		
	}
	
	void performTesting(SimpleMatrix R,SimpleMatrix H,int n,int m,ArrayList<Satellite> satList,ArrayList<Satellite> testedSatList,double[][] z,double[][] ze) throws Exception
	{
			// Pre-fit residual/innovation
			SimpleMatrix v = new SimpleMatrix(n, 1, true, innovation);
			SimpleMatrix P = kfObj.getCovariance();
			SimpleMatrix Cvv = ((H.mult(P).mult(H.transpose())).plus(R));
			SimpleMatrix Cvv_inv = Cvv.invert();
			double globalTq = v.transpose().mult(Cvv_inv).mult(v).get(0);
			ChiSquaredDistribution csd = new ChiSquaredDistribution(n);
			double alpha = 0.01;
			if (globalTq == 0) {
				throw new Exception("Error: T stat is zero");
			}
			// Detection
			double globalPVal = 1 - csd.cumulativeProbability(globalTq);
			int _n = testedSatList.size();
			while (globalPVal < alpha&&_n>(n/2)) {
				
				double max_w = Double.MIN_VALUE;
				int index = -1;
				for(int i =0;i<_n;i++)
				{
					SimpleMatrix cv = new SimpleMatrix(_n,1);
					cv.set(i, 1);
					double w = Math.abs(cv.transpose().mult(Cvv_inv).mult(v).get(0))/Math.sqrt(cv.transpose().mult(Cvv_inv).mult(cv).get(0));
					if (w > max_w) {
						max_w = w;
						index = i;
					}
					
				}
				satList.get(satList.indexOf(testedSatList.remove(index))).setOutlier(true);
				_n = testedSatList.size();
				SimpleMatrix R_ = new SimpleMatrix(_n, _n);
				double[][] z_ = new double[_n][1];
				double[][] ze_ = new double[_n][1];
				SimpleMatrix H_ = new SimpleMatrix(_n, 3 + (2 * m));
				SimpleMatrix v_ = new SimpleMatrix(_n,1);
				int j =0;
				for (int i = 0; i < _n+1; i++) {
					if (i!=index) {
						R_.set(j, j, R.get(i, i));
						v_.set(j, v.get(i));
						z_[j][0] = z[i][0];
						ze_[j][0] = ze[i][0];
						for (int k = 0; k < 3 + (2 * m); k++) {
							H_.set(j, k, H.get(i, k));							
						}
						j++;
					}
				}
				R = new SimpleMatrix(R_);
				H = new SimpleMatrix(H_);
				z = new double[_n][1];
				ze = new double[_n][1];
				v = new SimpleMatrix(v_);
				for (int i = 0; i < _n; i++) {
					z[i] = z_[i];
					ze[i] = ze_[i];
				}
				Cvv = ((H.mult(P).mult(H.transpose())).plus(R));
				Cvv_inv = Cvv.invert();
				globalTq = v.transpose().mult(Cvv_inv).mult(v).get(0);
				if (globalTq == 0) {
					throw new Exception("Error: T stat is zero");
				}
				csd = new ChiSquaredDistribution(_n-1);
				globalPVal = 1 - csd.cumulativeProbability(globalTq);
			}
			
			
	}
	
	

	private SimpleMatrix getJacobian(ArrayList<Satellite> satList, double[] estECEF, String[] obsvCodeList) {
		int m = obsvCodeList.length;
		int n = satList.size();
		int rows = n;
		double[][] H = new double[rows][3 + (2 * m)];

		for (int i = 0; i < n; i++) {
			Satellite sat = satList.get(i);
			String obsvCode = sat.getObsvCode();
			// Line of Sight vector
			// Its not really a ECI, therefore don't get confused
			double[] LOS = IntStream.range(0, 3).mapToDouble(j -> sat.getSatEci()[j] - estECEF[j]).toArray();
			// Geometric Range
			double GR = Math.sqrt(Arrays.stream(LOS).map(j -> j * j).reduce(0.0, (j, k) -> j + k));
			// Converting LOS to unit vector
			final int _i = i;
			IntStream.range(0, 3).forEach(j -> H[_i][j] = -LOS[j] / GR);

			for (int j = 0; j < m; j++) {
				if (obsvCode.equals(obsvCodeList[j])) {
					H[i][3 + j] = 1;
				}
			}

		}

		return new SimpleMatrix(H);

	}

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
