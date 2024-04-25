package com.gnssAug.helper;

import org.hipparchus.linear.RealMatrix;
import org.orekit.estimation.measurements.gnss.LambdaMethod;

public class ILS_LAMBDA extends LambdaMethod {

	public ILS_LAMBDA() {
		super();
		// TODO Auto-generated constructor stub
	}

	public double[] getDiag() {
		return getDiagReference();
	}

	public double[] getLow() {
		return getLowReference();
	}

	public double[] getDecorrelated() {
		return getDecorrelatedReference();
	}
	
}
