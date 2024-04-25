package com.gnssAug.helper.lambda;

import org.apache.commons.math3.distribution.NormalDistribution;

import Jama.Matrix;

public class Parsearch2 {
	private Matrix zhat = null;
	private Matrix Qzhat = null;
	private Matrix Z = null;
	private Matrix L = null;
	private Matrix D = null;

	private int ncands;
	private double mu;
	private Matrix zpar = null;
	private double[] sqnorm = null;
	private Matrix Qzpar = null;
	private Matrix Zpar = null;

	private int nfixed = 0;
	private Matrix zfixed = null;

	public Parsearch2(Matrix zhat, Matrix Qzhat, Matrix Z, Matrix L, Matrix D, double mu, int ncands) {
		this.zhat = zhat.copy();
		this.Qzhat = Qzhat.copy();
		this.Z = Z.copy();
		this.L = L.copy();
		this.D = D.copy();
		this.mu = mu;
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
		while ((k <= n) && zpar == null) {
			Ssearch ssearch = new Ssearch(zhat.getMatrix(k - 1, n - 1, new int[] { 0 }),
					L.getMatrix(k - 1, n - 1, k - 1, n - 1), D.getMatrix(k - 1, n - 1, k - 1, n - 1), ncands);
			zpar = ssearch.getafixed();
			if (zpar != null) {
				sqnorm = ssearch.getSqnorm();
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

}
