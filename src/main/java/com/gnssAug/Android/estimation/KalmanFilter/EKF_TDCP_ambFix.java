package com.gnssAug.Android.estimation.KalmanFilter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.stream.IntStream;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.MatrixFeatures_DDRM;
import org.ejml.dense.row.factory.DecompositionFactory_DDRM;
import org.ejml.interfaces.decomposition.CholeskyDecomposition_F64;
import org.ejml.interfaces.decomposition.CholeskyLDLDecomposition_F64;
import org.ejml.simple.SimpleMatrix;
import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.RealMatrix;
import org.orekit.estimation.measurements.gnss.IntegerLeastSquareSolution;
import org.orekit.estimation.measurements.gnss.LambdaMethod;
import org.orekit.estimation.measurements.gnss.SimpleRatioAmbiguityAcceptance;

import com.gnssAug.Android.constants.GnssDataConfig;
import com.gnssAug.Android.estimation.LinearLeastSquare;
import com.gnssAug.Android.estimation.KalmanFilter.Models.KFconfig;
import com.gnssAug.Android.models.CycleSlipDetect;
import com.gnssAug.Android.models.Satellite;
import com.gnssAug.helper.Decorrel;
import com.gnssAug.helper.FixingSolution;
import com.gnssAug.helper.ILS_LAMBDA;
import com.gnssAug.utility.Combination;
import com.gnssAug.utility.LatLonUtil;
import com.gnssAug.utility.Matrix;
import com.gnssAug.utility.SatUtil;
import com.gnssAug.utility.Weight;


public class EKF_TDCP_ambFix extends EKFParent {
	private long ambDetectedCount = 0;
	private long ambRepairedCount = 0;

	public EKF_TDCP_ambFix() {
		kfObj = new KFconfig();
	}

	protected TreeMap<Long, Integer> ambDetectedCountMap;
	protected TreeMap<Long, Integer> ambRepairedCountMap;

	public TreeMap<Long, double[]> process(TreeMap<Long, ArrayList<Satellite>> SatMap, ArrayList<Long> timeList,
			boolean useIGS, String[] obsvCodeList, boolean doAnalyze, boolean doTest, boolean outlierAnalyze,
			boolean innPhaseRate) throws Exception {
		boolean isWeighted = true;
		int m = obsvCodeList.length;
		int n = 3 + m;
		SimpleMatrix x = new SimpleMatrix(n, 1);
		SimpleMatrix P = new SimpleMatrix(n, n);
		/*
		 * state XYZ is intialized WLS generated Rx position estimated using first epoch
		 * data, a-priori estimate error covariance matrix for state XYZ is therefore
		 * assigned 25 m^2 value. Other state variables are assigned infinite(big)
		 * variance
		 */
		double[] refPos = LinearLeastSquare.getEstPos(SatMap.firstEntry().getValue(), true, useIGS);
		double[] intialVel = LinearLeastSquare.getEstVel(SatMap.firstEntry().getValue(), true, refPos, useIGS);
		IntStream.range(0, n).forEach(i -> x.set(i, intialVel[i]));

		// Total State
		IntStream.range(0, n).forEach(i -> P.set(i, i, 100));
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
		ambDetectedCountMap = new TreeMap<Long, Integer>();
		ambRepairedCountMap = new TreeMap<Long, Integer>();
		return iterate(SatMap, timeList, useIGS, obsvCodeList, doAnalyze, doTest, outlierAnalyze, isWeighted,
				innPhaseRate);
	}

	TreeMap<Long, double[]> iterate(TreeMap<Long, ArrayList<Satellite>> SatMap, ArrayList<Long> timeList,
			boolean useIGS, String[] obsvCodeList, boolean doAnalyze, boolean doTest, boolean outlierAnalyze,
			boolean isWeighted, boolean innPhaseRate) throws Exception {
		TreeMap<Long, double[]> estStateMap = new TreeMap<Long, double[]>();
		int m = obsvCodeList.length;
		int x_size = 3 + m;
		long prevTime = timeList.get(0);
		// Start from 2nd epoch
		for (int i = 1; i < timeList.size(); i++) {
			System.out.println("\n\n Epoch : " + i);
			long currentTime = timeList.get(i);
			double deltaT = (currentTime - prevTime) / 1e3;
			ArrayList<Satellite> currSatList = SatMap.get(currentTime);
			ArrayList<Satellite> prevSatList = SatMap.get(prevTime);
			ArrayList<CycleSlipDetect> csdList = new ArrayList<CycleSlipDetect>();
			double[] refPos = LinearLeastSquare.getEstPos(currSatList, true, useIGS);
			int n_curr = currSatList.size();
			int n_prev = prevSatList.size();
			for (int j = 0; j < n_curr; j++) {
				Satellite current_sat = currSatList.get(j);
				String satID = current_sat.getObsvCode() + current_sat.getSvid();
				for (int k = 0; k < n_prev; k++) {
					Satellite prev_sat = prevSatList.get(k);
					String prev_satID = prev_sat.getObsvCode() + prev_sat.getSvid();
					if (satID.equals(prev_satID)) {
						SimpleMatrix satEci = new SimpleMatrix(3, 1, true, current_sat.getSatEci());
						SimpleMatrix prev_satEci = new SimpleMatrix(3, 1, true, prev_sat.getSatEci());
						SimpleMatrix unitLOS = new SimpleMatrix(1, 3, true,
								SatUtil.getUnitLOS(current_sat.getSatEci(), refPos));
						double ionoRate = current_sat.getIonoErr() - prev_sat.getIonoErr();
						double tropoRate = current_sat.getTropoErr() - prev_sat.getTropoErr();
						double dopplerDR = ((current_sat.getRangeRate() + prev_sat.getRangeRate()) / 2) - tropoRate
								+ ionoRate;
						double phaseDR = current_sat.getPhase() - prev_sat.getPhase() + ionoRate;
						double satVelCorr = unitLOS.mult(satEci.minus(prev_satEci)).get(0);
						double wavelength = SpeedofLight / current_sat.getCarrierFrequencyHz();
						double approxCS = Math.abs(phaseDR - dopplerDR);
						if (innPhaseRate) {
							if (approxCS < 100 * wavelength) {
								csdList.add(new CycleSlipDetect(current_sat, dopplerDR, phaseDR, ionoRate, false,
										wavelength, satVelCorr, i, approxCS / wavelength, unitLOS));
							}
						} else {
							if (approxCS < 5 * wavelength) {
								csdList.add(new CycleSlipDetect(current_sat, dopplerDR, phaseDR, ionoRate, false,
										wavelength, satVelCorr, i, approxCS / wavelength, unitLOS));
							} else if (approxCS < 100 * wavelength) {
								csdList.add(new CycleSlipDetect(current_sat, dopplerDR, phaseDR, ionoRate, true,
										wavelength, satVelCorr, i, approxCS / wavelength, unitLOS));
							}
						}

					}
				}
			}

			runFilter(currentTime, deltaT, csdList, obsvCodeList, doAnalyze, doTest, outlierAnalyze, useIGS, i,
					isWeighted, refPos, innPhaseRate);
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
//			if (!MatrixFeatures_DDRM.isPositiveSemidefinite(P.getMatrix())) {
//				throw new Exception("PositiveDefinite test Failed");
//			}
			prevTime = currentTime;
		}
		return estStateMap;
	}

	// Innovation Based Testing
	private void runFilter(long currentTime, double deltaT, ArrayList<CycleSlipDetect> csdList, String[] obsvCodeList,
			boolean doAnalyze, boolean doTest, boolean outlierAnalyze, boolean useIGS, int ct, boolean isWeighted,
			double[] refPos, boolean innPhaseRate) throws Exception {

		// Satellite count
		int n = csdList.size();
		ArrayList<Satellite> satList = new ArrayList<Satellite>();
		for (int i = 0; i < n; i++) {
			satList.add(csdList.get(i).getSat());
		}
		int m = obsvCodeList.length;
		// Last update VC-matrix
		SimpleMatrix priorP = new SimpleMatrix(kfObj.getCovariance());

		// Assign Q and F matrix
		kfObj.configTDCP(deltaT, m, refPos);
		;
		kfObj.predict();
		SimpleMatrix x = kfObj.getState();
		SimpleMatrix z = new SimpleMatrix(n, 1);

		SimpleMatrix H = new SimpleMatrix(n, 3 + m);
		SimpleMatrix unitLOS = new SimpleMatrix(SatUtil.getUnitLOS(satList, refPos));
		H.insertIntoThis(0, 0, unitLOS.scale(-1));
		for (int i = 0; i < n; i++) {
			CycleSlipDetect csdObj = csdList.get(i);
			z.set(i, csdObj.getDopplerDR() - csdObj.getSatVelCorr());
			String obsvCode = satList.get(i).getObsvCode();
			for (int j = 0; j < m; j++) {
				if (obsvCodeList[j].equals(obsvCode)) {
					H.set(i, 3 + j, 1);
				}
			}
		}

		SimpleMatrix ze = H.mult(x);
		innovation = Matrix.matrix2ArrayVec(z.minus(ze));

		// Measurement Noise
		SimpleMatrix doppler_Cyy = Weight.getNormCyy(satList, GnssDataConfig.doppler_priorVarOfUnitW);
		SimpleMatrix R = new SimpleMatrix(doppler_Cyy);

		ArrayList<Satellite> testedSatList = new ArrayList<Satellite>(satList);
		if (doTest && !outlierAnalyze) {
			Object[] params = performTesting(R, H, n, m, satList, testedSatList, z, ze, csdList, false);
			R = (SimpleMatrix) params[0];
			H = (SimpleMatrix) params[1];
			z = (SimpleMatrix) params[2];
			ze = (SimpleMatrix) params[3];
		}
		// Perform Update Step
		kfObj.update(z, R, ze, H);

		// TDCP integration begins
		x = kfObj.getState();
		int ambCount = 0;
		if (innPhaseRate) {
			z = new SimpleMatrix(n, 1);
			H = new SimpleMatrix(n, 3 + m);
			unitLOS = new SimpleMatrix(SatUtil.getUnitLOS(satList, refPos));
			H.insertIntoThis(0, 0, unitLOS.scale(-1));
			for (int i = 0; i < n; i++) {
				CycleSlipDetect csdObj = csdList.get(i);
				z.set(i, csdObj.getCarrierPhaseDR() - csdObj.getSatVelCorr());
				String obsvCode = satList.get(i).getObsvCode();
				for (int j = 0; j < m; j++) {
					if (obsvCodeList[j].equals(obsvCode)) {
						H.set(i, 3 + j, 1);
					}
				}
			}
			ze = H.mult(x);
			innovation = Matrix.matrix2ArrayVec(z.minus(ze));
			for (int i = 0; i < n; i++) {
				CycleSlipDetect csdObj = csdList.get(i);
				double wavelength = csdObj.getWavelength();
				double approxCS = Math.abs(innovation[i]);
				if(approxCS>wavelength)
				{
					csdObj.setCS(true);
					ambCount++;
				}
					
			}
			
		} else {
			testedSatList = new ArrayList<Satellite>();
			ArrayList<CycleSlipDetect> testedCsdList = new ArrayList<CycleSlipDetect>();

			for (int i = 0; i < n; i++) {
				if (!csdList.get(i).isCS()) {
					testedSatList.add(csdList.get(i).getSat());
					testedCsdList.add(csdList.get(i));
				}
			}
			int tested_n = testedSatList.size();
			
			z = new SimpleMatrix(tested_n, 1);
			H = new SimpleMatrix(tested_n, 3 + m);
			SimpleMatrix testedUnitLOS = new SimpleMatrix(SatUtil.getUnitLOS(testedSatList, refPos));
			H.insertIntoThis(0, 0, testedUnitLOS.scale(-1));
			for (int i = 0; i < tested_n; i++) {
				CycleSlipDetect csdObj = testedCsdList.get(i);
				z.set(i, csdObj.getCarrierPhaseDR() - csdObj.getSatVelCorr());
				String obsvCode = testedSatList.get(i).getObsvCode();
				for (int j = 0; j < m; j++) {
					if (obsvCodeList[j].equals(obsvCode)) {
						H.set(i, 3 + j, 1);
					}
				}
			}

			ze = H.mult(x);
			innovation = Matrix.matrix2ArrayVec(z.minus(ze));
			// Measurement Noise
			SimpleMatrix tested_tdcp_Cyy = Weight.getNormCyy(testedSatList, GnssDataConfig.tdcp_priorVarOfUnitW);
			R = new SimpleMatrix(tested_tdcp_Cyy);

			// Testing for CS
			performTesting(R, H, tested_n, m, satList, testedSatList, z, ze, csdList, true);
			
			// Resume full SatList or CSDList
			for (int i = 0; i < n; i++) {
				if (csdList.get(i).isCS()) {
					ambCount++;
				}
			}
		}
		ambDetectedCountMap.put(currentTime, ambCount);
		ambDetectedCount += ambCount;
		SimpleMatrix x_new = new SimpleMatrix(3 + m + ambCount, 1);
		x_new.insertIntoThis(0, 0, x);
		SimpleMatrix P = kfObj.getCovariance();
		SimpleMatrix P_new = new SimpleMatrix(3 + m + ambCount, 3 + m + ambCount);
		P_new.insertIntoThis(0, 0, P);
		kfObj.setState_ProcessCov(x_new, P_new);
		for (int i = 0; i < ambCount; i++) {
			P_new.set(3 + m + i, 3 + m + i, 1e12);
		}
		H = new SimpleMatrix(n, 3 + m + ambCount);
		H.insertIntoThis(0, 0, unitLOS.scale(-1));
		z = new SimpleMatrix(n, 1);
		int ctr = 3 + m;
		for (int i = 0; i < n; i++) {
			CycleSlipDetect csdObj = csdList.get(i);
			z.set(i, csdObj.getCarrierPhaseDR() - csdObj.getSatVelCorr());
			String obsvCode = satList.get(i).getObsvCode();
			double wavelength = csdObj.getWavelength();
			for (int j = 0; j < m; j++) {
				if (obsvCodeList[j].equals(obsvCode)) {
					H.set(i, 3 + j, 1);
					if (csdObj.isCS()) {
						H.set(i, ctr, wavelength);
						ctr++;
					}
				}
			}
		}
		ze = H.mult(x_new);
		innovation = Matrix.matrix2ArrayVec(z.minus(ze));
		SimpleMatrix tdcp_Cyy = Weight.getNormCyy(satList, GnssDataConfig.tdcp_priorVarOfUnitW);
		R = new SimpleMatrix(tdcp_Cyy);
		// Perform Update Step
		kfObj.update(z, R, ze, H);
		x = kfObj.getState();
		P = kfObj.getCovariance();
		if (ambCount > 0) {
			SimpleMatrix floatAmb = x.extractMatrix(3 + m, 3 + m + ambCount, 0, 1);
			SimpleMatrix floatAmbCov = P.extractMatrix(3 + m, 3 + m + ambCount, 3 + m, 3 + m + ambCount);
			System.out.println("Float Ambiguity");
			System.out.println(floatAmb.toString());
			System.out.println("Float Ambiguity Covariance");
			System.out.println(floatAmbCov.toString());
			System.out.println("Fixed Ambiguity Sequence");
			RealMatrix Cxx_floatAmb_hat = new Array2DRowRealMatrix(Matrix.matrix2Array(floatAmbCov));
			
			CholeskyLDLDecomposition_F64<DMatrixRMaj> chol = DecompositionFactory_DDRM.cholLDL(ambCount);
			if( !chol.decompose(floatAmbCov.getMatrix()))
				   throw new RuntimeException("Cholesky failed!");
			double[] diagonal = chol.getDiagonal();
			
			LambdaMethod lm = new LambdaMethod();

			// Full ambiguity Resolution
			int[] indirection = IntStream.range(0, ambCount).toArray();
			IntegerLeastSquareSolution[] ILSsol = lm.solveILS(5, Matrix.matrix2ArrayVec(floatAmb), indirection,
					Cxx_floatAmb_hat);
			
//			Jama.Matrix Qahat = new Jama.Matrix(Matrix.matrix2Array(floatAmbCov));
//			Jama.Matrix ahat = new Jama.Matrix(Matrix.matrix2Array(floatAmb));
//			Decorrel decor = new Decorrel( Qahat, ahat);
//			SimpleMatrix Qzhat = new SimpleMatrix(decor.getQzhat().getArray());
//			SimpleMatrix zhat = new SimpleMatrix(decor.getzhat().getArray());
//			SimpleMatrix D = new SimpleMatrix(decor.getD().getArray());
//			double[] Ps = new double[D.numRows()];
//	        NormalDistribution normalDistribution = new NormalDistribution();
//			for (int i = 0;i < D.numRows();i++){
//	            double cdf = normalDistribution.probability(-Double.MAX_VALUE, 0.5/Math.sqrt(D.get(i,i)));
//	            if(i==0)
//	            {
//	            	Ps[i] =1*(2 * cdf - 1);
//	            }
//	            else
//	            {
//	            	Ps[i] =Ps[i-1]*(2 * cdf - 1);
//	            }
//	            
//	        }
//			if(ambCount>5)
//			{
//				System.out.println();
//			}
			if (ILSsol.length > 0) {
				SimpleRatioAmbiguityAcceptance ratioTest = new SimpleRatioAmbiguityAcceptance(1.0 / 3.0);
				IntegerLeastSquareSolution acceptedILS = ratioTest.accept(ILSsol);

				IntegerLeastSquareSolution finalSol = acceptedILS;
				int[] finalComb = indirection;
				// Partial Ambiguity Resolution
				if (finalSol == null && ambCount > 1) {
					int count = ambCount - 1;
					int[] arr = IntStream.range(0, ambCount).toArray();
					while (acceptedILS == null && count > 0) {
						ArrayList<int[]> combs = Combination.getCombination(arr, ambCount, count);

						double sqDist = Double.MAX_VALUE;
						for (int i = 0; i < combs.size(); i++) {
							int[] comb = combs.get(i);
							double[] newFloatAmb = IntStream.range(0, count).mapToDouble(j -> floatAmb.get(comb[j]))
									.toArray();

							ILSsol = lm.solveILS(5, newFloatAmb, comb, Cxx_floatAmb_hat);

							acceptedILS = ratioTest.accept(ILSsol);
							if (acceptedILS != null) {
								double newSqDist = acceptedILS.getSquaredDistance();
								if (newSqDist < sqDist) {
									sqDist = newSqDist;
									finalSol = new IntegerLeastSquareSolution(acceptedILS.getSolution(), newSqDist);
									finalComb = IntStream.range(0, comb.length).map(j -> comb[j]).toArray();
								}

							}
						}
						if (finalSol != null) {

							break;
						}
						count--;
					}
				}
				if (finalSol != null) {
					System.out.println(Arrays.toString(finalComb));
					FixingSolution.process(finalComb, x, P, finalSol, m);
					ambRepairedCountMap.put(currentTime, finalComb.length);
					ambRepairedCount += finalComb.length;
				}
			}
		}
		if (doAnalyze) {
			SimpleMatrix residual = z.minus((H.mult(x)));
			performAnalysis(residual, satList, R, H, currentTime, obsvCodeList);
		}
		x = x.extractMatrix(0, 3 + m, 0, 1);
		P = P.extractMatrix(0, 3 + m, 0, 3 + m);
		kfObj.setState_ProcessCov(x, P);

	}

	private Object[] performTesting(SimpleMatrix R, SimpleMatrix H, int n, int m, ArrayList<Satellite> satList,
			ArrayList<Satellite> testedSatList, SimpleMatrix z, SimpleMatrix ze, ArrayList<CycleSlipDetect> csdList,
			boolean isTDCP) throws Exception {
		// Pre-fit residual/innovation
		SimpleMatrix v = new SimpleMatrix(n, 1, true, innovation);
		SimpleMatrix P = kfObj.getCovariance();
		SimpleMatrix Cvv = ((H.mult(P).mult(H.transpose())).plus(R));
		SimpleMatrix Cvv_inv = Cvv.invert();
		double globalTq = v.transpose().mult(Cvv_inv).mult(v).get(0);
		ChiSquaredDistribution csd = new ChiSquaredDistribution(n);
		double alpha = 0.01;
		if (isTDCP) {
			alpha = 0.01;
		}
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
			int rem_index = satList.indexOf(testedSatList.remove(index));
			if (isTDCP) {

				satList.get(rem_index).setOutlier(true);
				csdList.get(rem_index).setCS(true);
			}
			_n = testedSatList.size();
			SimpleMatrix R_ = new SimpleMatrix(_n, _n);
			SimpleMatrix z_ = new SimpleMatrix(_n, 1);
			SimpleMatrix ze_ = new SimpleMatrix(_n, 1);
			SimpleMatrix H_ = new SimpleMatrix(_n, 3 + m);
			SimpleMatrix v_ = new SimpleMatrix(_n, 1);
			int j = 0;
			for (int i = 0; i < _n + 1; i++) {
				if (i != index) {
					R_.set(j, j, R.get(i, i));
					v_.set(j, v.get(i));
					z_.set(j, z.get(i));
					ze_.set(j, ze.get(i));
					for (int k = 0; k < 3 + m; k++) {
						H_.set(j, k, H.get(i, k));
					}
					j++;
				}
			}
			R = new SimpleMatrix(R_);
			H = new SimpleMatrix(H_);
			z = new SimpleMatrix(z_);
			ze = new SimpleMatrix(ze_);
			v = new SimpleMatrix(v_);

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

	private void performAnalysis(SimpleMatrix e_post_hat, ArrayList<Satellite> satList, SimpleMatrix R, SimpleMatrix H,
			long currentTime, String[] obsvCodeList) {

		int _n = satList.size();
		double[] measNoise = new double[_n];

		// Post-fit residual
		SimpleMatrix Cyy_inv = R.invert();

		// Compute Redundancies
		SimpleMatrix K = kfObj.getKalmanGain();
		SimpleMatrix HK = H.mult(K);

		double rZ = SimpleMatrix.identity(HK.numRows()).minus(HK).trace();

		double postVarOfUnitW = e_post_hat.transpose().mult(Cyy_inv).mult(e_post_hat).get(0) / rZ;

		for (int i = 0; i < satList.size(); i++) {
			double wavelength = SpeedofLight / satList.get(i).getCarrierFrequencyHz();
			innovation[i] = innovation[i] / wavelength;
			e_post_hat.set(i, e_post_hat.get(i) / wavelength);
		}
		satListMap.put(currentTime, satList);
		postVarOfUnitWMap.put(currentTime, postVarOfUnitW);
		residualMap.put(currentTime, Matrix.matrix2ArrayVec(e_post_hat));
		satCountMap.put(currentTime, (long) _n);
		measNoiseMap.put(currentTime, measNoise);

	}

	public long getAmbDetectedCount() {
		return ambDetectedCount;
	}

	public long getAmbRepairedCount() {
		return ambRepairedCount;
	}

	public TreeMap<Long, Integer> getAmbDetectedCountMap() {
		return ambDetectedCountMap;
	}

	public TreeMap<Long, Integer> getAmbRepairedCountMap() {
		return ambRepairedCountMap;
	}

}
