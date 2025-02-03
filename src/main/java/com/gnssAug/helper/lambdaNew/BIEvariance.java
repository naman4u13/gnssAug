package com.gnssAug.helper.lambdaNew;

import java.util.ArrayList;

import org.ejml.simple.SimpleMatrix;

import com.gnssAug.Android.constants.GnssDataConfig;
import com.gnssAug.helper.lambdaNew.Estimators.EstimatorBIE;
import com.gnssAug.utility.Matrix;

public class BIEvariance {

	public static void computeBIEvariance(EstimatorBIE estBIE,SimpleMatrix zFix, SimpleMatrix qHat, SimpleMatrix aHat) throws Exception
	{
		int n = aHat.numRows();
		SimpleMatrix allZcand = estBIE.getAllIntCandidates();
		int m = allZcand.numCols();
		SimpleMatrix qHat_inv = qHat.invert();
		SimpleMatrix TVector = estBIE.getExpTlist();
		double expTsum = TVector.elementSum();
		SimpleMatrix bie = allZcand.mult(TVector).scale(1/expTsum);
		
		if(Math.abs(bie.minus(zFix).elementSum())>1e-3)
		{
			throw new Exception("Something wrong with BIE esimtation");
		}
		SimpleMatrix bieRes = new SimpleMatrix(allZcand);
		for(int i=0;i<bieRes.numCols();i++)
		{
			bieRes.setColumn(i, 0, Matrix.matrix2ArrayVec(bieRes.extractVector(false,i).minus(bie)));
		}
		// conditional variance
		SimpleMatrix bieVar1 = bieRes.mult((TVector.diag()).scale(1/expTsum)).mult(bieRes.transpose());
		
		
		// Yu, Xianwen, Siqi Xia, and Wang Gao. "A practical method for calculating reliable integer float estimator in GNSS precise positioning." Survey Review 53.377 (2021): 97-107.
		SimpleMatrix dT_daHatT_list = new SimpleMatrix(m,n);
		SimpleMatrix dT_daHatT_sum = new SimpleMatrix(1,n);
		for(int i=0;i<m;i++)
		{
			SimpleMatrix zi =  allZcand.extractVector(false, i);
			double ti = TVector.get(i,0);
			SimpleMatrix dti_daHatT= (aHat.minus(zi)).transpose().mult(qHat_inv).scale(-ti);
			dT_daHatT_list.insertIntoThis(i, 0, dti_daHatT);
			dT_daHatT_sum = dT_daHatT_sum.plus(dti_daHatT);
			
		}
		
		SimpleMatrix dP_daHatT_sum = new SimpleMatrix(1,n);
		SimpleMatrix zdp_dahat_sum = new SimpleMatrix(n,n);
		
		for(int i=0;i<m;i++)
		{
			SimpleMatrix dti_daHatT =  dT_daHatT_list.extractVector(true, i);
			double ti = TVector.get(i,0);
			SimpleMatrix dpi_daHatT= (dti_daHatT.scale(expTsum).minus(dT_daHatT_sum.scale(ti))).scale(1/Math.pow(expTsum,2));
			dP_daHatT_sum = dP_daHatT_sum.plus(dpi_daHatT);
			SimpleMatrix zi =  allZcand.extractVector(false, i);
			SimpleMatrix zdp_dahat = zi.mult(dpi_daHatT);
			zdp_dahat_sum = zdp_dahat_sum.plus(zdp_dahat);
			
		}
		SimpleMatrix bieVar2 = zdp_dahat_sum.mult(qHat).mult(zdp_dahat_sum.transpose());
		
		SimpleMatrix bieVar3 = (SimpleMatrix) ComputeVariance.computeVariance(qHat, 3,0 , null, (int) GnssDataConfig.nSamplesMC, null)[0];
		System.out.println("BIE variance 1 : "+bieVar1);
		System.out.println("BIE variance 2 : "+bieVar2);
		System.out.println("BIE variance 3 : "+bieVar3);

	}
}
