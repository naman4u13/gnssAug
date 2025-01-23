package com.gnssAug.helper.lambdaNew;

import java.util.Arrays;

import org.apache.commons.math3.linear.*;

import com.gnssAug.helper.lambdaNew.models.DecorrelateResult;
import com.gnssAug.helper.lambdaNew.models.LambdaResult;
import com.gnssAug.utility.Vector;


public class Lambda {


	// Main function: LAMBDA implementation
	public static LambdaResult LAMBDA(double[] a_hat, double[][] Qa_hat, int method, Object... optionalParams) {
		//--------------------------------------------------------------------------
		// Problem dimensionality
		int nn = a_hat.length;

		// Check # of input arguments
		if (a_hat == null || Qa_hat == null) {
			throw new IllegalArgumentException(
					"ATTENTION: float ambiguity vector and its variance-covariance matrix are both needed in input!");
		}
		if (method < 0 || method > 9) {
			method = 3; // By default, we use ILS estimator (search-and-shrink)
		}


		// Check main inputs
		Utilities.checkMainInputs(Qa_hat, a_hat);

		//--------------------------------------------------------------------------   
		// Origin-translation of ambiguities (for numerical reasons)
		double[] a_origin = Utilities.fix(a_hat);
		a_hat = Vector.subtract(a_hat, a_origin);

		//--------------------------------------------------------------------------   
		// PRE-PROCESS: decorrelate ambiguities by an admissible Z-transformation
		DecorrelateResult decorrelated = DecorrelateVC.decorrelateVC(Qa_hat, a_hat);
		RealVector z_hat = decorrelated.getZ_hat();

		//--------------------------------------------------------------------------   
		// ADDITIONAL: compute success rate & number of fixed components
		double SR = (double) SuccessRate.computeSR_IBexact(decorrelated.getDz_vec().toArray())[0];
		int nFixed = nn;

		//--------------------------------------------------------------------------   
		// OPTIONAL PARAMETERS: set default values or get additional inputs
		int nCands = 1;
		double minSR = 0.99;
		String typeEstim = "ILS";
		int[] dimBlocks = { nn / 2, nn - nn / 2 };
		double maxFR = 0.1 / 100.0;
		double betaIAB = 0.5;
		double alphaBIE = 1e-6;

		if (method == 3 || method == 4 || method == 5) {
			nCands = optionalParams.length > 0 ? (int) optionalParams[0] : 1;
			if (method == 5) {
				minSR = optionalParams.length > 1 ? (double) optionalParams[1] : 0.99;
			}
		} else if (method == 6) {
			typeEstim = optionalParams.length > 0 ? (String) optionalParams[0] : "ILS";
			dimBlocks = optionalParams.length > 1 ? (int[]) optionalParams[1] : dimBlocks;
		} else if (method == 7) {
			maxFR = optionalParams.length > 0 ? (double) optionalParams[0] : 0.1 / 100.0;
		} else if (method == 8) {
			betaIAB = optionalParams.length > 0 ? (double) optionalParams[0] : 0.5;
		} else if (method == 9) {
			alphaBIE = optionalParams.length > 0 ? (double) optionalParams[0] : 1e-6;
		}

		//--------------------------------------------------------------------------   
		// METHODS: define the estimator
		RealVector z_fix = null;
		double sqnorm = 0;

		switch (method) {
		case 0: // Compute Float
			nFixed = 0;
			break;

		case 1: // Compute IR
			z_fix = Estimators.estimatorIR(z_hat);
			break;

		case 2: // Compute IB
			z_fix = Estimators.estimatorIB(z_hat, decorrelated.getLz_mat());
			break;

		case 3: // Compute ILS (shrink-and-search)
			Object[] ILS_soln = Estimators.estimatorILS(z_hat, decorrelated.getLz_mat(), decorrelated.getDz_vec().toArray(), nCands);
			z_fix = ((RealVector[])ILS_soln[0])[0];
			sqnorm = ((double[])ILS_soln[1])[0];
			break;

		/*  case 4: // Compute ILS (enumeration)
            double chi2 = computeInitialEllipsoid(z_hat, decorrelated.Lz_mat, decorrelated.dz_vec, nCands);
            z_fix = estimatorILS_enum(z_hat, decorrelated.Lz_mat, decorrelated.dz_vec, nCands, chi2);
            break;*/

		case 5: // Compute PAR
			Object[] parResult = Estimators.estimatorPAR(z_hat, decorrelated.getLz_mat(), decorrelated.getDz_vec().toArray(), nCands, minSR,0);
			z_fix = (RealVector)parResult[0];
			nFixed = (int) parResult[1];
			break;

//		case 6: // Compute VIB
//			VectorialIBResult vibResult = estimatorVIB(z_hat, decorrelated.Lz_mat, decorrelated.dz_vec, typeEstim, dimBlocks);
//			z_fix = vibResult.z_fix;
//			nFixed = vibResult.nFixed;
//			break;

		case 7: // Compute IA-FFRT
			Object[] iaffrtResult = Estimators.estimatorIAFFRT(z_hat, decorrelated.getLz_mat(), decorrelated.getDz_vec(), maxFR,null);
			z_fix = (RealVector) iaffrtResult[0];
			sqnorm = (double) iaffrtResult[1];
			nFixed = (int) iaffrtResult[2];
			break;

//		case 8: // Compute IAB
//			ApertureBootstrappingResult iabResult = estimatorIAB(z_hat, decorrelated.Lz_mat, decorrelated.dz_vec, betaIAB);
//			z_fix = iabResult.z_fix;
//			nFixed = iabResult.nFixed;
//			break;

		case 9: // Compute BIE
			double chi2_BIE = 2 * GammaIncompleteInverse.gammaincinv(1 - alphaBIE, nn / 2.0);
			z_fix = (RealVector) Estimators.estimatorBIE(z_hat, decorrelated.getLz_mat(), decorrelated.getDz_vec(), chi2_BIE)[0];
			break;

		default:
			throw new IllegalArgumentException("ATTENTION: the method selected is not available! Use 0-9.");
		}

		//--------------------------------------------------------------------------   
		// Check if fixed solution is rejected
		if (nFixed == 0) {
			double[] a_fix = Vector.add(a_hat, a_origin);
			return new LambdaResult(a_fix, 0, nFixed, SR, decorrelated.getZ_mat().getData(), decorrelated.getQz_hat().getData());
		}

		//--------------------------------------------------------------------------   
		// Back Z-transformation with translation to the old origin
		double[] a_fix = decorrelated.getiZt_mat().operate(z_fix).toArray();
		a_fix = Vector.add(a_fix, a_origin);

		//--------------------------------------------------------------------------   
		// Squared norm of ambiguity residuals
		if (Arrays.asList(1, 2, 5, 6, 8, 9).contains(method)) {
			
			 // Compute the residual vector: z_hat - z_fix
	        RealVector residual = z_hat.subtract(z_fix);

	        // Compute Dz_vec = Lz_mat' \ residual
	        // Note: Transpose Lz_mat and solve for Dz_vec
	        DecompositionSolver solver = new LUDecomposition(decorrelated.getLz_mat().transpose()).getSolver();
	        RealVector Dz_vec = solver.solve(residual);

	        // Compute squared norm: diag(Dz_vec' * (dz_vec' .\ Dz_vec))
	        double[] dz_vec = decorrelated.getDz_vec().toArray();
	        double[] sqnorm_dim = new double[dz_vec.length];
	        for (int i = 0; i < dz_vec.length; i++) {
	        	sqnorm_dim[i] = Dz_vec.getEntry(i) * Dz_vec.getEntry(i) / dz_vec[i];
	        }
	        // Sum all elements of sqnorm to get the final scalar value
	       sqnorm = Arrays.stream(sqnorm_dim).sum();
			
		}

		return new LambdaResult(a_fix, sqnorm, nFixed, SR, decorrelated.getZ_mat().getData(), decorrelated.getQz_hat().getData());
	}


}