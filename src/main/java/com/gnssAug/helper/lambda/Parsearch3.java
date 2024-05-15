package com.gnssAug.helper.lambda;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Locale;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.interpolation.UnivariateInterpolator;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.json.JSONArray;
import org.json.JSONObject;

import Jama.Matrix;

public class Parsearch3 {
	private Matrix zhat = null;
	private Matrix Qzhat = null;
	private Matrix Z = null;
	private Matrix L = null;
	private Matrix D = null;

	private double mu;
	private int ncands;
	private double P0;
	private Matrix zpar = null;
	private double[] sqnorm = null;
	private Matrix Qzpar = null;
	private Matrix Zpar = null;
	private double Ps = Double.NaN;
	private int nfixed = 0;
	private Matrix zfixed = null;

	public Parsearch3(Matrix zhat, Matrix Qzhat, Matrix Z, Matrix L, Matrix D, double P0, int ncands) {
		this.zhat = zhat.copy();
		this.Qzhat = Qzhat.copy();
		this.Z = Z.copy();
		this.L = L.copy();
		this.D = D.copy();
		this.P0 = P0;
		this.ncands = ncands;

		parsearch();
	}

	public Matrix getzpar() {
		if (zpar != null) {
			return zpar.copy();
		} else {
			return null;
		}
	}

	public double[] getSqnorm() {
		if (sqnorm != null) {
			return sqnorm.clone();
		} else {
			return null;
		}
	}

	public Matrix getQzpar() {
		if (Qzpar != null) {
			return Qzpar.copy();
		} else {
			return null;
		}
	}

	public Matrix getZpar() {
		if (Zpar != null) {
			return Zpar.copy();
		} else {
			return null;
		}
	}
	
	public double getPs(){
        return Ps;
    }

	
	public double getMu() {
		return mu;
	}

	public int getNfixed() {
		return nfixed;
	}

	public Matrix getzfixed() {
		if (zfixed != null) {
			return zfixed.copy();
		} else {
			return null;
		}
	}

	private void parsearch() {
		int n = Qzhat.getRowDimension();
		zpar = null;
		Qzpar = null;
		Zpar = null;
		sqnorm = null;
		zfixed = zhat;
		nfixed = 0;
		int k = 1;
		Ps = 1;
        NormalDistribution normalDistribution = new NormalDistribution();
		while ((k <= n) && zpar == null) {
			Ssearch ssearch = new Ssearch(zhat.getMatrix(k - 1, n - 1, new int[] { 0 }),
					L.getMatrix(k - 1, n - 1, k - 1, n - 1), D.getMatrix(k - 1, n - 1, k - 1, n - 1), ncands);
			zpar = ssearch.getafixed();
			Ps = 1;
			if (zpar != null) {
				for (int i = k - 1;i < D.getColumnDimension();i++){
	                double cdf = normalDistribution.probability(-Double.MAX_VALUE, 0.5/Math.sqrt(D.get(i,i)));
	                Ps *= (2 * cdf - 1);
	            }
				mu = 1;
				if(Ps < P0){
                    mu = ratioinv(1-P0,1-Ps,n);
                }
				sqnorm = ssearch.getSqnorm();
				System.out.println("Combination count: "+(n-k+1));
				System.out.println("failure rate :" + (1-Ps) );
				System.out.println("MU :" + mu );
				System.out.println("ratio :" + (sqnorm[0] / sqnorm[1]) );
				if ((sqnorm[0] / sqnorm[1]) < mu) {

					Qzpar = Qzhat.getMatrix(k - 1, n - 1, k - 1, n - 1).copy();
					Zpar = Z.getMatrix(0, n - 1, k - 1, n - 1);

					Matrix QP = Qzhat.getMatrix(0, k - 2, k - 1, n - 1)
							.times(Qzhat.getMatrix(k - 1, n - 1, k - 1, n - 1).inverse());
					if (k == 1) {
						zfixed = zpar.copy();
					} else {
						zfixed = new Matrix(n, ncands);
						for (int i = 1; i <= ncands; i++) {
							zfixed.setMatrix(0, k - 2, new int[] { i - 1 },
									zhat.getMatrix(0, k - 2, new int[] { 0 })
											.minus(QP.times(zhat.getMatrix(k - 1, n - 1, new int[] { 0 })
													.minus(zpar.getMatrix(0, n - k, new int[] { i - 1 })))));
						}
						zfixed.setMatrix(k - 1, n - 1, 0, ncands - 1, zpar);
					}
					nfixed = n - k + 1;
					
					return;
				} else {
					zpar = null;
				}

			}
			k += 1;
		}
	}
	
	 static public double ratioinv(double Pf_FIX, double PfILS, int n){
	        int kPf = (int) Math.round(Pf_FIX*1000);
	        if(kPf != 1 && kPf != 10) return Double.NaN;
	        if(n < 1) return Double.NaN;
	        double[] PfList = null;
	        double[] muList = null;
	        double[][] table_double = null;
	        try {
	        	File initialFile = new File("/Users/naman.agarwal/eclipse-workspace/gnssAug/src/main/res/ratiotab.txt");
	            InputStream ratioTableInputStream = new FileInputStream(initialFile);
	            String ratioTableStr = readTextInputStream(ratioTableInputStream);

	            JSONObject jsonObject = new JSONObject(ratioTableStr);
	            JSONArray jsonArray = jsonObject.getJSONArray(String.format(Locale.CHINESE,"table%d",kPf));

	            if(n > jsonArray.getJSONArray(0).length() - 1){
	                n = jsonArray.getJSONArray(0).length() - 1;
	            }

	            PfList = new double[jsonArray.length()];
	            muList = new double[jsonArray.length()];
	            for(int i = 0;i < jsonArray.length();i++){
	                PfList[i] = jsonArray.getJSONArray(i).getDouble(0);
	                muList[i] = jsonArray.getJSONArray(i).getDouble(n);
	            }
	        }catch (Exception e){
	            e.printStackTrace();
	            return Double.NaN;
	        }

	        UnivariateInterpolator interpolator = new LinearInterpolator();
	        UnivariateFunction function = interpolator.interpolate(PfList, muList);

	        return function.value(PfILS);
	    }
	 
	 static private String readTextInputStream(InputStream is) throws Exception {
	        InputStreamReader reader = new InputStreamReader(is);
	        BufferedReader bufferedReader = new BufferedReader(reader);
	        StringBuffer buffer = new StringBuffer("");
	        String str;
	        while ((str = bufferedReader.readLine()) != null) {
	            buffer.append(str);
	            buffer.append("\n");
	        }
	        return buffer.toString();
	    }

}
