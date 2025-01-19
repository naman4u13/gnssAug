package com.gnssAug.Android.estimation.KalmanFilter;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.TreeMap;
import java.util.stream.IntStream;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.ejml.dense.row.MatrixFeatures_DDRM;
import org.ejml.simple.SimpleMatrix;

import com.gnssAug.Android.constants.GnssDataConfig;
import com.gnssAug.Android.constants.Measurement;
import com.gnssAug.Android.estimation.LinearLeastSquare;
import com.gnssAug.Android.estimation.KalmanFilter.Models.KFconfig;
import com.gnssAug.Android.models.Satellite;
import com.gnssAug.utility.LatLonUtil;
import com.gnssAug.utility.Matrix;
import com.gnssAug.utility.SatUtil;
import com.gnssAug.utility.Vector;
import com.gnssAug.utility.Weight;
import com.opencsv.CSVWriter;

public class AKFDoppler_Static extends EKFParent {

	private HashMap<String, ArrayList<double[]>[][]> adaptVarMap = new HashMap<String, ArrayList<double[]>[][]>();
	private HashMap<String, ArrayList<double[]>[][]> adaptVarMapCont = new HashMap<String, ArrayList<double[]>[][]>();
	public AKFDoppler_Static() {
		kfObj = new KFconfig();
	}

	public TreeMap<Long, double[]> process(TreeMap<Long, ArrayList<Satellite>> SatMap, ArrayList<Long> timeList,
			boolean useIGS, String[] obsvCodeList, boolean doAnalyze, boolean doTest_vel, boolean isAdapt,boolean doTest_pos)
			throws Exception {

		for (String obsvCode : obsvCodeList) {
			adaptVarMap.put(obsvCode, new ArrayList[18][10]);
			adaptVarMapCont.put(obsvCode, new ArrayList[18][10]);
		}

		boolean isWeighted = true;
		int m = obsvCodeList.length;
		int n = 3 + m;
		double[][] _X = new double[n][1];
		double[][] P = new double[n][n];
		/*
		 * state XYZ is intialized WLS generated Rx position estimated using first epoch
		 * data, a-priori estimate error covariance matrix for state XYZ is therefore
		 * assigned 25 m^2 value. Other state variables are assigned infinite(big)
		 * variance
		 */
		double[] intialECEF = LinearLeastSquare.getEstPos(SatMap.firstEntry().getValue(), true, useIGS);
		IntStream.range(0, n).forEach(i -> _X[i][0] = intialECEF[i]);

		// Total State
		SimpleMatrix X = new SimpleMatrix(_X);
		IntStream.range(0, n).forEach(i -> P[i][i] = 100);
		// Error State intialized as zero
		double[][] x = new double[n][1];
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

		return iterate(X, SatMap, timeList, useIGS, obsvCodeList, doAnalyze, isWeighted, doTest_vel, isAdapt,doTest_pos);

	}

	TreeMap<Long, double[]> iterate(SimpleMatrix X, TreeMap<Long, ArrayList<Satellite>> SatMap,
			ArrayList<Long> timeList, boolean useIGS, String[] obsvCodeList, boolean doAnalyze, boolean isWeighted,
			boolean doTest_vel, boolean isAdapt,boolean doTest_pos) throws Exception {
		TreeMap<Long, double[]> estStateMap = new TreeMap<Long, double[]>();
		int m = obsvCodeList.length;
		int x_size = 3 + m;
		long time = timeList.get(0);
		
		// Start from 2nd epoch
		for (int i = 1; i < timeList.size(); i++) {
			long currentTime = timeList.get(i);
			double deltaT = (currentTime - time) / 1e3;
			ArrayList<Satellite> satList = SatMap.get(currentTime);
			
			runFilter(X, currentTime, deltaT, satList, obsvCodeList, doAnalyze, useIGS, i, isWeighted,
					isAdapt,doTest_pos);
			SimpleMatrix P = kfObj.getCovariance();
			double[] estState = new double[x_size];
			IntStream.range(0, x_size).forEach(j -> estState[j] = X.get(j));
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
			if (!MatrixFeatures_DDRM.isPositiveSemidefinite(P.getMatrix())) {
				throw new Exception("PositiveDefinite test Failed");
			}
			time = currentTime;
		}

		boolean makeCSV = false;
		if (makeCSV) {
			String filePath = "/Users/naman.agarwal/Library/CloudStorage/OneDrive-UniversityofCalgary/GPS/ION-GNSS-2024/Plots/T-A-SIS-10_urban_static/Xiaomi_Mi_10T_Pro_AKF.csv";
			File file = new File(filePath);
			try {
				// create FileWriter object with file as parameter
				FileWriter outputfile = new FileWriter(file);
				// create CSVWriter object filewriter object as parameter
				CSVWriter writer = new CSVWriter(outputfile);
				// create a List which contains String array
				List<String[]> entryList = new ArrayList<String[]>();
				String[] header = new String[] { "ObsvCode", "ElevAngle", "CN0", "ResidualSq", "Redundancy","SVID","AdaptVar" };
				writer.writeNext(header);
				for (String key : adaptVarMapCont.keySet()) {
					ArrayList<double[]>[][] dataArray = adaptVarMapCont.get(key);
					for (int i = 0; i < 18; i++) {
						for (int j = 0; j < 10; j++) {
							ArrayList<double[]> dataList = dataArray[i][j];
							
							if (dataList != null) {
								
								for (double[] data : dataList) {

									String[] entry = new String[7];
									entry[0] = key;
									entry[1] = i + "";
									entry[2] = j + "";
									entry[3] = data[0] + "";
									entry[4] = data[1] + "";
									entry[5] = data[2] + "";
									entry[6] = data[3] + "";
									entryList.add(entry);
								}
							}

						}
					}
				}
				writer.writeAll(entryList);
				writer.close();
			} catch (IOException err) {
				// TODO Auto-generated catch block
				err.printStackTrace();
			}
		}

		return estStateMap;
	}

	

	// Innovation Based Testing
	private void runFilter(SimpleMatrix X, long currentTime, double deltaT, ArrayList<Satellite> satList,
			String[] obsvCodeList, boolean doAnalyze, boolean useIGS, int ct,
			boolean isWeighted, boolean isAdapt,boolean doTest_pos) throws Exception {

		boolean useAndroidW = false;
		// Satellite count
		int n = satList.size();
		int m = obsvCodeList.length;

		// Last update VC-matrix
		SimpleMatrix priorP = new SimpleMatrix(kfObj.getCovariance());
		SimpleMatrix priorX = new SimpleMatrix(X);
		// Assign Q and F matrix
		kfObj.configAKFStatic(deltaT, m);
		kfObj.predict();

		double[] estPos = new double[] { X.get(0), X.get(1), X.get(2) };
		double[] rxClkOff = new double[m];// in meters
		for (int i = 0; i < m; i++) {
			rxClkOff[i] = X.get(i + 3);
		}

		/*
		 * H is the Jacobian matrix of partial derivatives Observation StateModel(h) of
		 * with respect to x
		 */
		double[][] _H = new double[n][3 + m];
		// Measurement vector
		double[][] z = new double[n][1];
		// Estimated Measurement vector
		double[][] ze = new double[n][1];
		// Measurement Noise
		double[][] _R = new double[n][n];
		innovation = new double[n];
		temp_innovation = new double[n];

		if (isWeighted) {
			if (useAndroidW) {
				for (int i = 0; i < n; i++) {
					_R[i][i] = Math.pow(satList.get(i).getReceivedSvTimeUncertaintyNanos() * SpeedofLight * 1e-9, 2);
				}
			} else {
//				LinearLeastSquare.getEstPos(satList, true, true, false, false, useIGS);
//				SimpleMatrix Cyy = LinearLeastSquare.getCyy_updated(Measurement.Pseudorange);
//				_R = Matrix.matrix2Array(Cyy);
				SimpleMatrix Cyy = Weight.getNormCyy(satList, GnssDataConfig.pseudorange_priorVarOfUnitW);
				_R = Matrix.matrix2Array(Cyy);
			}
		} else {
			for (int i = 0; i < n; i++) {
				_R[i][i] = GnssDataConfig.pseudorange_priorVarOfUnitW;
			}
		}
		/*
		 * Notice measurement 'z' is difference b/w PR and estimated PR, because it's a
		 * complimentary filter predicted error state is zero and therefore Hx = 0 or
		 * estimated measurement(ze) is zero
		 */
		for (int i = 0; i < n; i++) {
			Satellite sat = satList.get(i);
			// Its not really a ECI, therefore don't get confused
			String obsvCode = sat.getObsvCode();
			double[] satEcef = sat.getSatEci();
			double PR = sat.getPseudorange();
			double[] LOS = IntStream.range(0, 3).mapToDouble(j -> satEcef[j] - estPos[j]).toArray();
			// Approx Geometric Range
			double approxGR = Vector.mod(LOS);
			// Approx Pseudorange Range
			double approxPR = approxGR;
			for (int j = 0; j < m; j++) {
				if (obsvCode.equals(obsvCodeList[j])) {
					approxPR += rxClkOff[j];
					_H[i][3 + j] = 1;
				}
			}

			z[i][0] = PR - approxPR;
			final int _i = i;
			IntStream.range(0, 3).forEach(j -> _H[_i][j] = -LOS[j] / approxGR);
			innovation[i] = z[i][0] - ze[i][0];
		}
		SimpleMatrix R = new SimpleMatrix(_R);
		SimpleMatrix H = new SimpleMatrix(_H);
		kfObj.update(z, R, ze, H);
		SimpleMatrix x = kfObj.getState();
		for (int i = 0; i < 3 + m; i++) {
			X.set(i, X.get(i) + x.get(i));
			x.set(i, 0);
		}
		ArrayList<Satellite> testedSatList = new ArrayList<Satellite>(satList);
		if (isAdapt) {
			R = adapt(X, satList, obsvCodeList, H, R, n);
			kfObj.setState_ProcessCov(new SimpleMatrix(3 + m, 1), priorP);
			for (int i = 0; i < 3 + m; i++) {
				X.set(i, priorX.get(i));
			}
			kfObj.predict();
			if(doTest_pos)
			{
				
				performTesting(R, H, n, m, satList, testedSatList, z, ze);

			}
			kfObj.update(z, R, ze, H);
			x = kfObj.getState();
			for (int i = 0; i < 3 + m; i++) {
				X.set(i, X.get(i) + x.get(i));
				x.set(i, 0);
			}
		}

		if (doAnalyze) {
			performAnalysis(X, satList, R, H, priorP, currentTime, n, obsvCodeList, ct);
		}

	}

	private void performAnalysis(SimpleMatrix X, ArrayList<Satellite> satList, SimpleMatrix R, SimpleMatrix H,
			SimpleMatrix priorP, long currentTime, int n, String[] obsvCodeList, int ct) {

		double[] measNoise = new double[n];
		double[] residual = (double[]) get_z_ze_res(X, satList, obsvCodeList)[2];
		// Post-fit residual
		SimpleMatrix e_post_hat = new SimpleMatrix(n, 1, true, residual);
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
		if (n - rSum > 0.01) {
			System.err.println("FATAL ERROR: Redundancy sum is wrong");
		}
		double postVarOfUnitW = e_post_hat.transpose().mult(Cyy_inv).mult(e_post_hat).get(0) / rZ;
		redundancyList.add(new double[] { n, rSum, rX, rW, rZ });
		postVarOfUnitWMap.put(currentTime, postVarOfUnitW);
		residualMap.put(currentTime, residual);
		satCountMap.put(currentTime, (long) n);
		satListMap.put(currentTime, satList);
		measNoiseMap.put(currentTime, measNoise);
		for (int i = 0; i < n; i++) {
			Satellite sat = satList.get(i);
			sat.setPrVar(R.get(i,i));
			
		}

	}

	private SimpleMatrix adapt(SimpleMatrix X, ArrayList<Satellite> satList, String[] obsvCodeList, SimpleMatrix H,
			SimpleMatrix R, int n) {
		double[] residual = (double[]) get_z_ze_res(X, satList, obsvCodeList)[2];
		// Compute Redundancies
		SimpleMatrix K = kfObj.getKalmanGain();
		SimpleMatrix HK = H.mult(K);
		SimpleMatrix redunMat = SimpleMatrix.identity(HK.numRows()).minus(HK);
		SimpleMatrix adaptR = new SimpleMatrix(R);
		for (int i = 0; i < n; i++) {
			Satellite sat = satList.get(i);
			String obsvCode = sat.getObsvCode();
			double elevAng = Math.toDegrees(sat.getElevAzm()[0]);
			double cn0 = sat.getCn0DbHz();
			
			int i1 = (int) (elevAng / 5);
			i1 = i1 < 18 ? i1 : 17;
			int i2 = (int) (cn0 / 5);
			i2 = i2 < 10 ? i2 : 9;
			if (adaptVarMap.get(obsvCode)[i1][i2] == null) {
				adaptVarMap.get(obsvCode)[i1][i2] = new ArrayList<double[]>();
				adaptVarMapCont.get(obsvCode)[i1][i2] = new ArrayList<double[]>();
			}
			ArrayList<double[]> varList = adaptVarMap.get(obsvCode)[i1][i2];
			if(varList.size()>10)
			{
				varList.remove(0);
			}
			varList.add(new double[] { Math.pow(residual[i], 2), redunMat.get(i, i),sat.getSvid() });
			if (varList.size() < 5) {
				adaptVarMapCont.get(obsvCode)[i1][i2].add(new double[] { Math.pow(residual[i], 2), redunMat.get(i, i),sat.getSvid(),adaptR.get(i,i) });
				continue;
			}
			double num = varList.stream().mapToDouble(j -> j[0]).sum();
			double denom = varList.stream().mapToDouble(j -> j[1]).sum();
			double adaptVar = num / denom;
			adaptR.set(i, i, adaptVar);
			adaptVarMapCont.get(obsvCode)[i1][i2].add(new double[] { Math.pow(residual[i], 2), redunMat.get(i, i),sat.getSvid(),adaptR.get(i,i) });

		}
		return adaptR;

	}

	private Object[] get_z_ze_res(SimpleMatrix X, ArrayList<Satellite> satList, String[] obsvCodeList) {
		int n = satList.size();
		int m = obsvCodeList.length;
		double[][] z = new double[n][1];
		double[][] ze = new double[n][1];
		double[] residual = new double[n];
		double[] estPos = new double[] { X.get(0), X.get(1), X.get(2) };
		double[] rxClkOff = new double[m];// in meters
		for (int i = 0; i < m; i++) {
			rxClkOff[i] = X.get(i + 3);
		}
		for (int i = 0; i < n; i++) {
			Satellite sat = satList.get(i);
			String obsvCode = sat.getObsvCode();
			double approxPR = Math.sqrt(IntStream.range(0, 3).mapToDouble(j -> estPos[j] - sat.getSatEci()[j])
					.map(j -> j * j).reduce(0, (j, k) -> j + k));
			for (int j = 0; j < m; j++) {
				if (obsvCode.equals(obsvCodeList[j])) {
					approxPR += rxClkOff[j];
				}
			}
			z[i][0] = sat.getPseudorange() - approxPR;
			residual[i] = z[i][0] - ze[i][0];
		}

		return new Object[] { z, ze, residual };
	}
	
	private Object[] performTesting(SimpleMatrix R, SimpleMatrix H, int n, int m, ArrayList<Satellite> satList,
			ArrayList<Satellite> testedSatList, double[][] z, double[][] ze) throws Exception {
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
		while (globalPVal < alpha && _n > (n / 2)) {

			double max_w = Double.MIN_VALUE;
			int index = -1;
			for (int i = 0; i < _n; i++) {
				SimpleMatrix cv = new SimpleMatrix(_n, 1);
				cv.set(i, 1);
				double w = Math.abs(cv.transpose().mult(Cvv_inv).mult(v).get(0))
						/ Math.sqrt(cv.transpose().mult(Cvv_inv).mult(cv).get(0));
				if (w > max_w) {
					max_w = w;
					index = i;
				}
			}
			//System.out.print(" " + index + " ");
			satList.get(satList.indexOf(testedSatList.remove(index))).setOutlier(true);
			_n = testedSatList.size();
			SimpleMatrix R_ = new SimpleMatrix(_n, _n);
			double[][] z_ = new double[_n][1];
			double[][] ze_ = new double[_n][1];
			SimpleMatrix H_ = new SimpleMatrix(_n, 3 + m);
			SimpleMatrix v_ = new SimpleMatrix(_n, 1);
			int j = 0;
			for (int i = 0; i < _n + 1; i++) {
				if (i != index) {
					R_.set(j, j, R.get(i, i));
					v_.set(j, v.get(i));
					z_[j][0] = z[i][0];
					ze_[j][0] = ze[i][0];
					for (int k = 0; k < 3 + m; k++) {
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
			csd = new ChiSquaredDistribution(_n);
			globalPVal = 1 - csd.cumulativeProbability(globalTq);
		}
		return new Object[] { R, H, z, ze };

	}

}
