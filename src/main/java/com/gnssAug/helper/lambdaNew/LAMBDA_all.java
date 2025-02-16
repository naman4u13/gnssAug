package com.gnssAug.helper.lambdaNew;

import java.util.HashMap;

import org.ejml.data.DMatrixRMaj;

//LAMBDA 4.0 | Least-squares AMBiguity Decorrelation Adjustment toolbox
//The class represents the main routine of LAMBDA 4.0 implementation. By
//providing a float ambiguity vector and its associated variance-covariance
//(vc-)matrix, the class allows resolving for the ambiguity. By default,
//the Integer Least Square (ILS, search-and-shrink) solution is adopted as 
//the optimal Integer (I-)estimator, but other estimators are also available.
//
//More information can be found in the official LAMBDA 4.0 Documentation 
//[RD01], or in other references as [RD02][RD03][RD04][RD05].
//
//-------------------------------------------------------------------------
//INPUTS:
//aHat       Ambiguity float vector (column)
//qaHat      Variance-covariance matrix of the original ambiguities
//method     Estimator (0-9) adopted, see section "METHODS"
//varArgs    Optional input parameters, which replace default values
//
//OUTPUTS:
//aFix       Ambiguity fixed vector (column)
//sqNorm     Squared norm of the ambiguity residuals (aHat - aFix)
//nFixed     Number of integer-fixed ambiguity components
//sr         Success rate (bootstrapping) for Full Ambiguity Resolution
//zMat       Admissible Z-transformation matrix (unimodular)
//qzHat      Variance-covariance matrix of the decorrelated ambiguities
//
//DEPENDENCIES:
//Several functionalities from LAMBDA 4.0 toolbox.
//    > Use "addpath('LAMBDA_toolbox')"
//
//REFERENCES:
//[RD01] Massarweh, L., Verhagen, S., and Teunissen, P.J.G. (2024). New 
//    LAMBDA toolbox for mixed-integer models: Estimation and Evaluation. 
//    GPS Solut NN, XXX (2024), submitted. DOI: not yet available.
//[RD02] Teunissen, P.J.G. (1993, August). Least-squares estimation of 
//    the integer GPS ambiguities. In Invited lecture, section IV theory 
//    and methodology, IAG general meeting, Beijing, China (pp. 1-16). 
//[RD03] Teunissen, P.J.G. (1995, November). The least-squares ambiguity 
//    decorrelation adjustment: a method for fast GPS integer ambiguity 
//    estimation. Journal of Geodesy 70, 65–82. 
//    DOI: 10.1007/BF00863419
//[RD04] De Jonge, P., Tiberius, C.C.J.M. (1996). The LAMBDA method for 
//    integer ambiguity estimation: implementation aspects. Publications 
//    of the Delft Computing Centre, LGR-Series, 12(12), 1-47.
//[RD05] Verhagen, S. (2005) The GNSS integer ambiguities: Estimation 
//    and validation. PhD thesis, Delft University of Technology.
//
//-------------------------------------------------------------------------
//METHODS:
//0 - Float solution
//1 - Integer Rounding (IR)
//2 - Integer Bootstrapping (IB)
//3 - Integer Least-Squares (ILS) by search-and-shrink         [DEFAULT]
//4 - Integer Least-Squares (ILS) by enumeration
//5 - Partial Ambiguity Resolution (PAR), based on ILS estimator
//6 - Vectorial IB (VIB), based on IR or ILS estimator
//7 - Integer Aperture with Fixed Failure-rate Ratio Test (IA-FFRT)
//8 - Integer Aperture Bootstrapping (IAB)
//9 - Best Integer Equivariant (BIE)
//
//-------------------------------------------------------------------------
//OPTIONS:
//Customized inputs used in different METHODS, such as
//    nCands    = Number of integer candidates.
//    minSR     = Minimum success rate threshold for PAR.
//    typeEstim = Estimator used in VIB partitioned blocks.
//    dimBlocks = Dimensionality of each VIB block.
//    maxFR     = Maximum failure rate threshold, within [0.05%-1%].
//    betaIAB   = Aperture coefficient for IAB.
//    alphaBIE  = Probability for BIE approximation.
//
//-------------------------------------------------------------------------
//Copyright: Geoscience & Remote Sensing department @ TUDelft | 01/06/2024
//Contact email:    LAMBDAtoolbox-CITG-GRS@tudelft.nl
//-------------------------------------------------------------------------
//Created by
//01/06/2024  - Lotfi Massarweh
//    Implementation for LAMBDA 4.0 toolbox, based on LAMBDA 3.0
//
//Modified by
//dd/mm/yyyy  - Name Surname (author)
//    >> Changes made in this new version
//-------------------------------------------------------------------------
//START

import org.ejml.simple.SimpleMatrix;
import com.gnssAug.Android.constants.GnssDataConfig;
import com.gnssAug.helper.lambdaNew.ComputeSR_IBexact.SR_IB;
import com.gnssAug.helper.lambdaNew.Estimators.*;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorBIE.EstimatorBIEResult;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorIA_FFRT.IAFFRTResult;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorILS.ILSResult;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorPAR.PARResult;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorPAR_FFRT.PARResult_FFRT;

public class LAMBDA_all {

	/**
	 * Encapsulates the results of the LAMBDA computation.
	 */
	public static class LambdaAllResult {
		private SimpleMatrix aFix; // Ambiguity fixed vector (column)
		private SimpleMatrix qFix; // Variance of fixed ambiguity vector

		private int nFixed; // Number of integer-fixed ambiguity components
		private double sr; // Success rate (bootstrapping) for Full Ambiguity Resolution
		private double approxSR;
		private double approxFR;
		public LambdaAllResult(SimpleMatrix aFix, SimpleMatrix qFix, int nFixed, double sr) {
			this.aFix = aFix;
			this.qFix = qFix;
			this.nFixed = nFixed;
			this.sr = sr;

		}
		public LambdaAllResult(SimpleMatrix aFix, SimpleMatrix qFix, int nFixed, double sr,double approxSR,double approxFR) {
			this.aFix = aFix;
			this.qFix = qFix;
			this.nFixed = nFixed;
			this.sr = sr;
			this.approxSR = approxSR;
			this.approxFR = approxFR;

		}

		public SimpleMatrix getaFix() {
			return aFix;
		}

		public SimpleMatrix getqFix() {
			return qFix;
		}

		public int getnFixed() {
			return nFixed;
		}

		public double getSr() {
			return sr;
		}
		public double getApproxSR()
		{
			return approxSR;
		}
		public double getApproxFR()
		{
			return approxFR;
		}

	}

	public static HashMap<EstimatorType, LambdaAllResult> computeLambda(SimpleMatrix aHat, SimpleMatrix qaHat,
			boolean estimateVar) throws Exception {
		// Problem dimensionality
		int nn = aHat.numRows();

		// Check number of input arguments
		if (aHat == null || qaHat == null) {
			throw new IllegalArgumentException(
					"ATTENTION: float ambiguity vector and its variance-covariance matrix are both needed in input!");
		}

		// Set default method if not provided
		// In Java, method is a required parameter, so defaulting is handled by the
		// caller
		// If method is not provided, it should be set to 3 (ILS estimator)

		// Check main inputs: "qaHat" & "aHat".
		Utilities.checkMainInputs(qaHat, aHat);

		// Origin-translation of ambiguities | Only for "numerical" reasons
		// Round toward zero, so the new origin is within (-1,1)
		SimpleMatrix aOrigin = aHat.copy();
		for (int i = 0; i < aOrigin.numRows(); i++) {
			aOrigin.set(i, 0, Math.floor(aOrigin.get(i, 0)));
		}
		aHat = aHat.minus(aOrigin);

		// PRE-PROCESS: decorrelate ambiguities by an admissible Z-transformation
		DecorrelateVCResult decorrelationResult = DecorrelateVC.decorrelateVC(qaHat, aHat);
		SimpleMatrix qzHat = decorrelationResult.getQzHat();
		SimpleMatrix lzMat = decorrelationResult.getLzMat();
		double[] dzVec = decorrelationResult.getDzVec();
		SimpleMatrix iZtMat = decorrelationResult.getIZtMat();
		SimpleMatrix zMat = decorrelationResult.getZMat();
		SimpleMatrix zHat = decorrelationResult.getZHat();

		// ADDITIONAL: computation of success rate & number of fixed components
		SR_IB srResult = ComputeSR_IBexact.computeSR_IBexact(dzVec);
		HashMap<EstimatorType,Integer> nFixedMap = new HashMap<EstimatorType,Integer>();
		double sr = srResult.getSR();
		// OPTIONAL PARAMETERS: set default values or get additional inputs
		int nCands = 1;
		double minSR = 0.99;
		double maxFR = 1 / 100.0;
		double alphaBIE = 1e-6;

		ILSResult ilsResult = new EstimatorILS().estimatorILS(zHat, lzMat, dzVec, nCands);
		PARResult parResult = EstimatorPAR.estimatorPAR(zHat, lzMat, dzVec, nCands, minSR, null, estimateVar);
		IAFFRTResult iaFfrtResult = new EstimatorIA_FFRT().estimatorIA_FFRT(zHat, lzMat, dzVec, maxFR, null);
		EstimatorBIE estBIE = new EstimatorBIE();
		// Call BIE-estimator (recursive implementation)
		double chi2BIE = 2.0 * GammaIncompleteInverse.gammaincinv(1.0 - alphaBIE, nn / 2.0);
		EstimatorBIEResult BieResult = estBIE.estimatorBIE(zHat, lzMat, dzVec, chi2BIE, null, qzHat);
		PARResult_FFRT parResult_ffrt = EstimatorPAR_FFRT.estimatorPAR_FFRT(zHat, lzMat, dzVec, nCands, minSR,
				estimateVar);

		HashMap<EstimatorType, SimpleMatrix> aFixMap = new HashMap<EstimatorType, SimpleMatrix>();
		HashMap<EstimatorType, SimpleMatrix> qFixMap = new HashMap<EstimatorType, SimpleMatrix>();
		HashMap<EstimatorType, double[]> srfrMap = new HashMap<EstimatorType, double[]>();

		aFixMap.put(EstimatorType.ILS, ilsResult.getAFix());
		aFixMap.put(EstimatorType.PAR, parResult.getaPAR());
		aFixMap.put(EstimatorType.IA_FFRT, iaFfrtResult.getaFix());
		aFixMap.put(EstimatorType.BIE, BieResult.getaBIE());
		aFixMap.put(EstimatorType.PAR_FFRT, parResult_ffrt.getaPAR());
		
		nFixedMap.put(EstimatorType.ILS, nn);
		nFixedMap.put(EstimatorType.PAR, parResult.getnFixed());
		nFixedMap.put(EstimatorType.IA_FFRT, iaFfrtResult.getnFixed());
		nFixedMap.put(EstimatorType.BIE, nn);
		nFixedMap.put(EstimatorType.PAR_FFRT, parResult_ffrt.getnFixed());

		qFixMap.put(EstimatorType.PAR, (SimpleMatrix) parResult.getStats()[0]);
		qFixMap.put(EstimatorType.PAR_FFRT, (SimpleMatrix) parResult_ffrt.getStats()[0]);
		if (estimateVar) {
			HashMap<EstimatorType, Object[]> varCalResMap = ComputeVariance.computeVarianceAll(qzHat, 0, maxFR,
					(int) GnssDataConfig.nSamplesMC, iaFfrtResult.getMuRatio());
			qFixMap.put(EstimatorType.ILS, (SimpleMatrix) varCalResMap.get(EstimatorType.ILS)[0]);

			qFixMap.put(EstimatorType.IA_FFRT, (SimpleMatrix) varCalResMap.get(EstimatorType.IA_FFRT)[0]);
			qFixMap.put(EstimatorType.BIE, (SimpleMatrix) varCalResMap.get(EstimatorType.BIE)[0]);

			srfrMap.put(EstimatorType.ILS, new double[] { (double) varCalResMap.get(EstimatorType.ILS)[1],
					(double) varCalResMap.get(EstimatorType.ILS)[2] });
			srfrMap.put(EstimatorType.PAR,
					new double[] { (double) parResult.getStats()[1], (double) parResult.getStats()[2] });
			srfrMap.put(EstimatorType.IA_FFRT, new double[] { (double) varCalResMap.get(EstimatorType.IA_FFRT)[1],
					(double) varCalResMap.get(EstimatorType.IA_FFRT)[2] });
			srfrMap.put(EstimatorType.BIE, new double[] { (double) varCalResMap.get(EstimatorType.BIE)[1],
					(double) varCalResMap.get(EstimatorType.BIE)[2] });
			srfrMap.put(EstimatorType.PAR_FFRT,
					new double[] { (double) parResult_ffrt.getStats()[1], (double) parResult_ffrt.getStats()[2] });
		} else {
			qFixMap.put(EstimatorType.ILS, new SimpleMatrix(nn, nn));
			qFixMap.put(EstimatorType.BIE, new SimpleMatrix(nn, nn));
			if (iaFfrtResult.getnFixed() == 0) {
				qFixMap.put(EstimatorType.IA_FFRT, new SimpleMatrix(qzHat));
			} else {
				qFixMap.put(EstimatorType.IA_FFRT, new SimpleMatrix(nn, nn));
			}
			qFixMap.put(EstimatorType.ILS, new SimpleMatrix(nn, nn));
		}

		HashMap<EstimatorType, LambdaAllResult> result = new HashMap<EstimatorType, LambdaAllResult>();
		for (EstimatorType est : new EstimatorType[] { EstimatorType.ILS, EstimatorType.PAR, EstimatorType.IA_FFRT,
				EstimatorType.BIE, EstimatorType.PAR_FFRT }) {
			// Back Z-transformation with translation to the old origin
			SimpleMatrix aFix = iZtMat.mult(aFixMap.get(est));
			aFix = aFix.plus(aOrigin);
			aFixMap.put(est, aFix);

			SimpleMatrix qFix = iZtMat.mult(qFixMap.get(est)).mult(iZtMat.transpose());
			qFixMap.put(est, qFix);
			System.out.println("Fixed Ambiguity Variance : " + est.toString());
			System.out.println(qFix.toString());
			int nFixed = nFixedMap.get(est);
			if (estimateVar) {
				double approxSR = (double) srfrMap.get(est)[0];
				double approxFR = (double) srfrMap.get(est)[1];
				System.out.println("Approximate Success Rate : " + est.toString() + "  " + approxSR * 100);
				System.out.println("Approximate Failure Rate : " + est.toString() + "  " + approxFR * 100);
				result.put(est, new LambdaAllResult(aFix, qFix, nFixed, sr,approxSR,approxFR));
				
			}
			else {
			result.put(est, new LambdaAllResult(aFix, qFix, nFixed, sr));}
		}
		return result;

	}

}
