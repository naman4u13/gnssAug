package com.gnssAug.helper.lambdaNew;

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

import com.gnssAug.utility.Matrix;
import com.gnssAug.Android.constants.GnssDataConfig;
import com.gnssAug.helper.lambda.Decorrel;
import com.gnssAug.helper.lambdaNew.ComputeSR_IBexact.SR_IB;
import com.gnssAug.helper.lambdaNew.ComputeVariance.VarianceResult;
import com.gnssAug.helper.lambdaNew.Estimators.*;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorBIE.EstimatorBIEResult;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorIA_FFRT.IAFFRTResult;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorILS.ILSResult;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorPAR.PARResult;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorPAR_FFRT.PARResult_FFRT;

public class LAMBDA {

	/**
	 * Encapsulates the results of the LAMBDA computation.
	 */
	public static class LambdaResult {
		private SimpleMatrix aFix; // Ambiguity fixed vector (column)
		private SimpleMatrix qFix; // Variance of fixed ambiguity vector
		private double sqNorm; // Squared norm of the ambiguity residuals
		private int nFixed; // Number of integer-fixed ambiguity components
		private double sr; // Success rate (bootstrapping) for Full Ambiguity Resolution
		private SimpleMatrix zMat; // Admissible Z-transformation matrix (unimodular)
		private SimpleMatrix qzHat; // Variance-covariance matrix of the decorrelated ambiguities

		public LambdaResult(SimpleMatrix aFix, SimpleMatrix qFix, double sqNorm, int nFixed, double sr, SimpleMatrix zMat,
				SimpleMatrix qzHat) {
			this.aFix = aFix;
			this.qFix = qFix;
			this.sqNorm = sqNorm;
			this.nFixed = nFixed;
			this.sr = sr;
			this.zMat = zMat;
			this.qzHat = qzHat;
		}

		public SimpleMatrix getaFix() {
			return aFix;
		}
		public SimpleMatrix getqFix()
		{
			return qFix;
		}
		public double getSqNorm() {
			return sqNorm;
		}

		public int getnFixed() {
			return nFixed;
		}

		public double getSr() {
			return sr;
		}

		public SimpleMatrix getzMat() {
			return zMat;
		}

		public SimpleMatrix getQzHat() {
			return qzHat;
		}
		
	}

	/**
	 * Computes the LAMBDA adjustment.
	 *
	 * @param aHat    Ambiguity float vector (column)
	 * @param qaHat   Variance-covariance matrix of the original ambiguities
	 * @param method  Estimator (0-9) adopted, see METHODS section
	 * @param varArgs Optional input parameters, which replace default values
	 * @return LambdaResult containing aFix, sqNorm, nFixed, sr, zMat, and qzHat
	 * @throws IllegalArgumentException if input arguments are insufficient or
	 *                                  invalid
	 */
	public static LambdaResult computeLambda(SimpleMatrix aHat, SimpleMatrix qaHat, int method, Object... varArgs) {
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
		int nFixed = nn;
		double sr = srResult.getSR();
		// OPTIONAL PARAMETERS: set default values or get additional inputs
		int nCands = 1;
		double minSR = 0.99;
		String typeEstim = "ILS";
		int[] dimBlocks = new int[] { (int) Math.floor(nn / 2.0), (int) Math.ceil(nn / 2.0) };
		double maxFR = 1 / 100.0;
		double betaIAB = 0.5;
		double alphaBIE = 1e-6;
		
		if (varArgs != null) {
			if (varArgs.length > 0 && (method == 3 || method == 4 || method == 5|| method == 10)) {
				nCands = (int) varArgs[0];
				if (varArgs.length > 1 && (method == 5||method == 10)) {
					minSR = (double) varArgs[1];
				}
			}
			if (varArgs.length > 0 && method == 6) {
				typeEstim = (String) varArgs[0];
				if (varArgs.length > 1) {
					dimBlocks = (int[]) varArgs[1];
				}
			}
			if (varArgs.length > 0 && method == 7) {
				maxFR = (double) varArgs[0];
			}
			if (varArgs.length > 0 && method == 8) {
				betaIAB = (double) varArgs[0];
			}
			if (varArgs.length > 0 && method == 9) {
				alphaBIE = (double) varArgs[0];
			}
		}

		SimpleMatrix zFix = null;
		SimpleMatrix QzFix = null;
		double sqNorm = 0.0;

		// METHODS: define the estimator (see LAMBDA 4.0 toolbox Documentation)
		switch (method) {
		case 0: // Compute Float
			nFixed = 0;
			break;
		case 1: // Compute IR
			zFix = EstimatorIR.estimatorIR(zHat);
			break;
		case 2: // Compute IB
			zFix = EstimatorIB.estimatorIB(zHat, lzMat);
			break;
		case 3: // Compute ILS (shrink-and-search) [DEFAULT]
			ILSResult ilsResult = EstimatorILS.estimatorILS(zHat, lzMat, dzVec, nCands);
			zFix = ilsResult.getAFix().extractVector(false, 0);
			sqNorm = ilsResult.getSqNorm()[0];
			QzFix =  ComputeVariance.computeVariance(qzHat, 1, 0, null,(int) GnssDataConfig.nSamplesMC).getVariance();
			break;
//         case 4: // Compute ILS (enumeration) based on an initial ellipsoid
//             double chi2 = computeInitialEllipsoid(zHat, lzMat, dzVec, nCands);
//             EstimatorILSEnumResult ilsEnumResult = estimatorILS_enum(zHat, lzMat, dzVec, nCands, chi2);
//             zFix = ilsEnumResult.zFix;
//             sqNorm = ilsEnumResult.sqNorm;
//             break;
		case 5: // Compute PAR
			PARResult parResult = EstimatorPAR.estimatorPAR(zHat, lzMat, dzVec, nCands, minSR, null);
			zFix = parResult.getaPAR();
			nFixed = parResult.getnFixed();
			sr = parResult.getSR_PAR();
			QzFix = parResult.getQPAR();
			break;
//         case 6: // Compute VIB
//             EstimatorVIBResult vibResult = estimatorVIB(zHat, lzMat, dzVec, typeEstim, dimBlocks);
//             zFix = vibResult.zFix;
//             nFixed = vibResult.nFixed;
//             break;
		case 7: // Compute IA-FFRT (ILS with Fixed Failure-rate Ratio Test)
			IAFFRTResult iaFfrtResult = EstimatorIA_FFRT.estimatorIA_FFRT(zHat, lzMat, dzVec, maxFR, null);
			zFix = iaFfrtResult.getaFix();
			sqNorm = iaFfrtResult.getsqNorm();
			nFixed = iaFfrtResult.getnFixed();
			QzFix =  ComputeVariance.computeVariance(qzHat, 2, 0, maxFR,(int) GnssDataConfig.nSamplesMC).getVariance();
			break;
//         case 8: // Compute IAB (Integer Aperture Bootstrapping)
//             EstimatorIABResult iabResult = estimatorIAB(zHat, lzMat, dzVec, betaIAB);
//             zFix = iabResult.zFix;
//             nFixed = iabResult.nFixed;
//             break;
		case 9: // Compute BIE based on chi-squared inverse CDF
			double chi2BIE = 2.0 * GammaIncompleteInverse.gammaincinv(1.0 - alphaBIE, nn / 2.0);

			// Call BIE-estimator (recursive implementation)
			EstimatorBIEResult BieResult = EstimatorBIE.estimatorBIE(zHat, lzMat, dzVec, chi2BIE, null);
			zFix = BieResult.getaBIE();
			break;
		case 10: // Compute PAR
			PARResult_FFRT parResult_ffrt = EstimatorPAR_FFRT.estimatorPAR_FFRT(zHat, lzMat, dzVec, nCands, minSR, null);
			zFix = parResult_ffrt.getaPAR();
			nFixed = parResult_ffrt.getnFixed();
			sr = parResult_ffrt.getSR_PAR();
			QzFix = parResult_ffrt.getQPAR();
			break;
		default:
			throw new IllegalArgumentException("ATTENTION: the method selected is not available! Use 0-10.");
		}

		// Check if fixed solution is rejected, e.g. METHOD = 7 (IA-FFRT) or 8 (IAB)
		if (nFixed == 0) {
			SimpleMatrix aFix = aHat.plus(aOrigin); // Back-translation to the old origin
			double finalSqNorm = 0.0; // Squared norm of float vector is zero
			return new LambdaResult(aFix,null, finalSqNorm, nFixed, sr, zMat, qzHat);
		}

		// Back Z-transformation with translation to the old origin
		SimpleMatrix aFix = iZtMat.mult(zFix);
		aFix = aFix.plus(aOrigin);
		
		if(QzFix!=null)
		{
			SimpleMatrix QFix = iZtMat.mult(QzFix).mult(iZtMat.transpose());
			System.out.println("Fixed Ambiguity Variance");
			System.out.println(QFix.toString());
			
		}
		
		// Squared norm of ambiguity residuals (invariant to any Z-transformations)
		if (method == 1 || method == 2 || method == 5 || method == 6 || method == 8 || method == 9|| method == 10) {
			SimpleMatrix dzInverse = new SimpleMatrix(dzVec.length, dzVec.length);
			for (int i = 0; i < dzVec.length; i++) {
				dzInverse.set(i, 0, 1.0 / dzVec[i]);
			}
			SimpleMatrix dzDiv = dzInverse;

			double residual = 0.0;
			for (int i = 0; i < dzVec.length; i++) {
				residual += Math.pow(dzDiv.get(i, 0) * (zHat.get(i, 0) - zFix.get(i, 0)), 2);
			}
			sqNorm = residual;
		}

		return new LambdaResult(aFix, QzFix,sqNorm, nFixed, sr, zMat, qzHat);
	}

}
