package com.gnssAug.helper;

import java.util.HashSet;

import org.ejml.simple.SimpleMatrix;
import org.orekit.estimation.measurements.gnss.IntegerLeastSquareSolution;

public class FixingSolution {

	// Notation follow from ENGO 625 notes Fall 2023
		public static void process(int[] comb, SimpleMatrix x, SimpleMatrix Cxx_hat, IntegerLeastSquareSolution fixAmb,
				int m) {
			int n = x.getNumElements();

			int l = comb.length;

			SimpleMatrix a_hat = new SimpleMatrix(l, 1);
			SimpleMatrix a_inv_hat = new SimpleMatrix(l, 1);
			SimpleMatrix C_a = new SimpleMatrix(l, l);
			HashSet<Integer> exclude = new HashSet<Integer>();
			for (int i = 0; i < l; i++) {
				a_hat.set(i, x.get(3 + m + comb[i]));
				a_inv_hat.set(i, fixAmb.getSolution()[i]);
				exclude.add(3 + m + comb[i]);
				for (int j = 0; j < l; j++) {
					C_a.set(i, j, Cxx_hat.get(3 + m + comb[i], 3 + m + comb[j]));
				}
			}
			SimpleMatrix C_a_inv = C_a.invert();

			for (int i = 0; i < n; i++) {
				if (!exclude.contains(i)) {
					SimpleMatrix C_ba = new SimpleMatrix(1, l);
					SimpleMatrix b_hat = new SimpleMatrix(1, 1);
					SimpleMatrix C_b_hat = new SimpleMatrix(1, 1);
					b_hat.set(0, x.get(i));
					C_b_hat.set(0, Cxx_hat.get(i, i));
					for (int j = 0; j < l; j++) {
						C_ba.set(0, j, Cxx_hat.get(3 + m + comb[j], i));
					}
					SimpleMatrix b_inv_hat = b_hat.minus(C_ba.mult(C_a_inv).mult((a_hat.minus(a_inv_hat))));
					SimpleMatrix Cb_inv_hat = C_b_hat.minus(C_ba.mult(C_a_inv).mult(C_ba.transpose()));
					x.set(i, b_inv_hat.get(0));
					Cxx_hat.set(i, i, Cb_inv_hat.get(0));
				}

			}

		}
}
