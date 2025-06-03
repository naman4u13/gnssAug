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
import com.gnssAug.helper.lambda.Decorrel;
import com.gnssAug.helper.lambdaNew.ComputeSR_IBexact.SR_IB;
import com.gnssAug.helper.lambdaNew.Estimators.*;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorBIE.EstimatorBIEResult;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorIA_FFRT.IAFFRTResult;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorILS.ILSResult;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorPAR.PARResult;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorPAR_FFRT.PARResult_FFRT;
import com.gnssAug.utility.Matrix;

public class LAMBDA {

	/**
	 * Encapsulates the results of the LAMBDA computation.
	 */
	public static class LambdaResult {
		private SimpleMatrix aFix; // Ambiguity fixed vector (column)
		private SimpleMatrix QaFix; // Variance of fixed ambiguity vector
		private double sqNorm; // Squared norm of the ambiguity residuals
		private int nFixed; // Number of integer-fixed ambiguity components
		private double sr; // Success rate (bootstrapping) for Full Ambiguity Resolution
		private SimpleMatrix zMat; // Admissible Z-transformation matrix (unimodular)
		private SimpleMatrix qzHat; // Variance-covariance matrix of the decorrelated ambiguities

		public LambdaResult(SimpleMatrix aFix, SimpleMatrix QaFix, double sqNorm, int nFixed, double sr,
				SimpleMatrix zMat, SimpleMatrix qzHat) {
			this.aFix = aFix;
			this.QaFix = QaFix;
			this.sqNorm = sqNorm;
			this.nFixed = nFixed;
			this.sr = sr;
			this.zMat = zMat;
			this.qzHat = qzHat;
		}

		public SimpleMatrix getaFix() {
			return aFix;
		}

		public SimpleMatrix getQaFix() {
			return QaFix;
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
	 * @param QaHat   Variance-covariance matrix of the original ambiguities
	 * @param method  Estimator (0-9) adopted, see METHODS section
	 * @param varArgs Optional input parameters, which replace default values
	 * @return LambdaResult containing aFix, sqNorm, nFixed, sr, zMat, and qzHat
	 * @throws Exception 
	 * @throws IllegalArgumentException if input arguments are insufficient or
	 *                                  invalid
	 */
	public static LambdaResult computeLambda(SimpleMatrix aHat, SimpleMatrix QaHat, EstimatorType method,
			boolean estimateVar, Object... varArgs) throws Exception {
		// Problem dimensionality
		int nn = aHat.numRows();

		// Check number of input arguments
		if (aHat == null || QaHat == null) {
			throw new IllegalArgumentException(
					"ATTENTION: float ambiguity vector and its variance-covariance matrix are both needed in input!");
		}

		// Set default method if not provided
		// In Java, method is a required parameter, so defaulting is handled by the
		// caller
		// If method is not provided, it should be set to 3 (ILS estimator)

		// Check main inputs: "qaHat" & "aHat".
		Utilities.checkMainInputs(QaHat, aHat);

		// Origin-translation of ambiguities | Only for "numerical" reasons
		// Round toward zero, so the new origin is within (-1,1)
		SimpleMatrix aOrigin = aHat.copy();
		for (int i = 0; i < aOrigin.numRows(); i++) {
			aOrigin.set(i, 0, Math.floor(aOrigin.get(i, 0)));
		}
		aHat = aHat.minus(aOrigin);

		// PRE-PROCESS: decorrelate ambiguities by an admissible Z-transformation
		DecorrelateVCResult decorrelationResult = DecorrelateVC.decorrelateVC(QaHat, aHat);
		SimpleMatrix QzHat = decorrelationResult.getQzHat();
		SimpleMatrix lzMat = decorrelationResult.getLzMat();
		double[] dzVec = decorrelationResult.getDzVec();
		SimpleMatrix iZtMat = decorrelationResult.getIZtMat();
		SimpleMatrix zMat = decorrelationResult.getZMat();
		SimpleMatrix zHat = decorrelationResult.getZHat();
		
		
		
		  //Compute Z matrix based on the decomposition  Q=L^T*D*L
//        Decorrel decorrel = new Decorrel(new Jama.Matrix(Matrix.matrix2Array(qaHat)),new Jama.Matrix(Matrix.matrix2Array(aHat)));
//        SimpleMatrix qzHat = new SimpleMatrix(decorrel.getQzhat().getArray());
//        SimpleMatrix zMat = new SimpleMatrix(decorrel.getZ().getArray());
//        SimpleMatrix lzMat = new SimpleMatrix(decorrel.getL().getArray());
//        double[] dzVec = Matrix.matrix2ArrayVec(new SimpleMatrix(decorrel.getD().getArray()));
//        SimpleMatrix zHat = new SimpleMatrix(decorrel.getzhat().getArray());
//        SimpleMatrix iZtMat = new SimpleMatrix(decorrel.getiZt().getArray());
		
//		double orgAmbDcrNo = Matrix.computeDecorrelationNumber(qaHat);
//		double zTransAmbDcrNo = Matrix.computeDecorrelationNumber(qzHat);
		double orgAmbDcrNo = Math.sqrt(Matrix.computeCorrelationMatrix(QaHat).normF());
		double zTransAmbDcrNo = Math.sqrt(Matrix.computeCorrelationMatrix(QzHat).normF());
		System.out.println("Original Float Ambiguities Decorrelation no. :"+orgAmbDcrNo);
		System.out.println("z-tranformed Float Ambiguities Decorrelation no. :"+zTransAmbDcrNo);
		System.out.println("Increase in decorrelation: " + (orgAmbDcrNo-zTransAmbDcrNo)); 
		System.out.println("z-transformed float ambiguity variance: \n" + QzHat.toString()); 
		System.out.println("z-transformed float ambiguity: \n" + zHat.toString()); 
		System.out.println("z Matrix : \n" + zMat.toString()); 

		// ADDITIONAL: computation of success rate & number of fixed components
		SR_IB srResult = ComputeSR_IBexact.computeSR_IBexact(dzVec);
		int nFixed = nn;
		double sr = srResult.getSR();
		// OPTIONAL PARAMETERS: set default values or get additional inputs
		int nCands = 1;
		double minSR = 0.99;
		double maxFR = 1 / 100.0;
		double alphaBIE = 1e-6;
		
		if (varArgs != null) {
			if (varArgs.length > 0 && (method == EstimatorType.ILS || method == EstimatorType.PAR
					|| method == EstimatorType.PAR_FFRT)) {
				nCands = (int) varArgs[0];
				if (varArgs.length > 1 && (method == EstimatorType.PAR || method == EstimatorType.PAR_FFRT)) {
					minSR = (double) varArgs[1];
				}
			}

			if (varArgs.length > 0 && (method == EstimatorType.IA_FFRT)) {
				maxFR = (double) varArgs[0];
			}

			if (varArgs.length > 0 && (method == EstimatorType.BIE)) {
				alphaBIE = (double) varArgs[0];
			}
		}

		SimpleMatrix zFix = null;
		SimpleMatrix QzFix = null;
		double sqNorm = 0.0;
		Object[] stats = null;
		// METHODS: define the estimator (see LAMBDA 4.0 toolbox Documentation)
		switch (method) {

		case ILS: // Compute ILS (shrink-and-search) [DEFAULT]
			ILSResult ilsResult = new EstimatorILS().estimatorILS(zHat, lzMat, dzVec, nCands);
			zFix = ilsResult.getAFix().extractVector(false, 0);
			sqNorm = ilsResult.getSqNorm()[0];
			if (estimateVar) {
				stats =  ComputeVariance.computeVariance(QzHat, 1, 0, null, (int) GnssDataConfig.nSamplesMC, null);
				QzFix = (SimpleMatrix) stats[0];
						
			} else {
				QzFix = new SimpleMatrix(nn, nn);
			}
			break;

		case PAR: // Compute PAR
			PARResult parResult = EstimatorPAR.estimatorPAR(zHat, lzMat, dzVec, nCands, minSR, null, estimateVar);
			zFix = parResult.getaPAR();
			nFixed = parResult.getnFixed();
			sr = parResult.getSR_PAR();
			stats  = parResult.getStats();
			QzFix = (SimpleMatrix) stats[0];

			break;
		case IA_FFRT: // Compute IA-FFRT (ILS with Fixed Failure-rate Ratio Test)
			IAFFRTResult iaFfrtResult = new EstimatorIA_FFRT().estimatorIA_FFRT(zHat, lzMat, dzVec, maxFR, null);
			zFix = iaFfrtResult.getaFix();
			sqNorm = iaFfrtResult.getsqNorm();
			nFixed = iaFfrtResult.getnFixed();
			if (estimateVar&&nFixed!=0) {
				stats = ComputeVariance
						.computeVariance(QzHat, 2, 0, maxFR, (int) GnssDataConfig.nSamplesMC, iaFfrtResult.getMuRatio());
				QzFix = (SimpleMatrix) stats[0];
			} else {
				if(nFixed==0)
				{
					QzFix = new SimpleMatrix(QzHat);
				}
				else
				{
					QzFix = new SimpleMatrix(nn, nn);
				}
				
			}
			break;

		case BIE: // Compute BIE based on chi-squared inverse CDF
			double chi2BIE = 2.0 * GammaIncompleteInverse.gammaincinv(1.0 - alphaBIE, nn / 2.0);
			EstimatorBIE estBIE = new EstimatorBIE();
			// Call BIE-estimator (recursive implementation)
			EstimatorBIEResult BieResult = estBIE.estimatorBIE(zHat, lzMat, dzVec, chi2BIE, null,QzHat);
			zFix = BieResult.getaBIE();
			stats =ComputeVariance.computeVariance(QzHat, 3,0 , null, (int) GnssDataConfig.nSamplesMC, null);
			QzFix = (SimpleMatrix) stats[0];
			//BIEvariance.computeBIEvariance(estBIE, zFix,qzHat,zHat);
			break;
		case PAR_FFRT: // Compute PAR-FFRt
			PARResult_FFRT parResult_ffrt = EstimatorPAR_FFRT.estimatorPAR_FFRT(zHat, lzMat, dzVec, nCands, minSR,
					estimateVar);
			zFix = parResult_ffrt.getaPAR();
			nFixed = parResult_ffrt.getnFixed();
			sr = parResult_ffrt.getSR_PAR();
			stats  = parResult_ffrt.getStats();
			QzFix = (SimpleMatrix) stats[0];
			break;
		
			
		default:
			throw new IllegalArgumentException("ATTENTION: the method selected is not available! Use 0-10.");
		}

		// Check if fixed solution is rejected, e.g. METHOD = 7 (IA-FFRT) or 8 (IAB)
		if (nFixed == 0) {
			SimpleMatrix aFix = aHat.plus(aOrigin); // Back-translation to the old origin
			double finalSqNorm = 0.0; // Squared norm of float vector is zero
			
			return new LambdaResult(aFix, QaHat, finalSqNorm, nFixed, sr, zMat, QzHat);
		}

		// Back Z-transformation with translation to the old origin
		SimpleMatrix aFix = iZtMat.mult(zFix);
		aFix = aFix.plus(aOrigin);

		SimpleMatrix QaFix = iZtMat.mult(QzFix).mult(iZtMat.transpose());
		System.out.println("Fixed Ambiguity Variance");
		System.out.println(QaFix.toString());

		// Squared norm of ambiguity residuals (invariant to any Z-transformations)
		if (method == EstimatorType.PAR || method == EstimatorType.BIE || method == EstimatorType.PAR_FFRT) {
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
		if(estimateVar)
		{
			double approxSR = (double) stats[1];
			double approxFR = (double) stats[2];
			System.out.println("Approximate Success Rate : " + approxSR*100);
			System.out.println("Approximate Failure Rate : " + approxFR*100);
		}

		return new LambdaResult(aFix, QaFix, sqNorm, nFixed, sr, zMat, QzHat);
	}

}
