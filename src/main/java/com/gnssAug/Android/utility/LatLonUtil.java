package com.gnssAug.Android.utility;

import org.ejml.simple.SimpleMatrix;

public class LatLonUtil {

	// All are WGS-84 params
	// Semi-major axis or Equatorial radius
	public static final double a = 6378137;
	// flattening
	public static final double f = 1 / 298.257223563;
	// Semi-minor axis or Polar radius
	public static final double b = 6356752.314245;
	private static final double e = Math.sqrt((Math.pow(a, 2) - Math.pow(b, 2)) / Math.pow(a, 2));
	public static final double e2 = Math.sqrt((Math.pow(a, 2) - Math.pow(b, 2)) / Math.pow(b, 2));
	// Earth Angular rate in rad/s
	public static final double omega_ie = 7.292115e-5;
	// Earth’s gravitational constant in m^3/s^2
	public static final double mu = 3.986004418e14;

	public static double[] ecef2lla(double[] ECEF) {

		double x = ECEF[0];
		double y = ECEF[1];
		double z = ECEF[2];
		double[] lla = { 0, 0, 0 };
		double lat = 0, lon, height, N, theta, p;

		p = Math.sqrt(Math.pow(x, 2) + Math.pow(y, 2));

		theta = Math.atan((z * a) / (p * b));

		lon = Math.atan(y / x);

		if (x < 0) {
			if (y > 0) {
				lon = Math.PI + lon;
			} else {
				lon = -Math.PI + lon;
			}
		}
		for (int i = 0; i < 3; i++) {
			lat = Math.atan(((z + Math.pow(e2, 2) * b * Math.pow(Math.sin(theta), 3))
					/ ((p - Math.pow(e, 2) * a * Math.pow(Math.cos(theta), 3)))));
			theta = Math.atan((Math.tan(lat) * b) / (a));

		}

		N = a / (Math.sqrt(1 - (Math.pow(e, 2) * Math.pow(Math.sin(lat), 2))));

		height = (p * Math.cos(lat)) + ((z + (Math.pow(e, 2) * N * Math.sin(lat))) * Math.sin(lat)) - N;

		lon = lon * 180 / Math.PI;

		lat = lat * 180 / Math.PI;
		lla[0] = lat;
		lla[1] = lon;
		lla[2] = height;
		return lla;
	}

	public static double[] lla2ecef(double[] lla, boolean isDegree) {

		double finv = 1 / f;
		double esq = (2 - 1 / finv) / finv;
		double dtr = 1;
		if (isDegree) {
			dtr = Math.PI / 180;
		}

		double lat = lla[0] * dtr; // rad
		double lon = lla[1] * dtr; // rad
		double alt = lla[2]; // m

		double N = a / Math.sqrt(1 - esq * Math.pow(Math.sin(lat), 2));

		double x = (N + alt) * Math.cos(lat) * Math.cos(lon);
		double y = (N + alt) * Math.cos(lat) * Math.sin(lon);
		double z = ((1 - esq) * N + alt) * Math.sin(lat);

		double[] ecef = { x, y, z }; // m
		return ecef;
	}

	public static double getHaversineDistance(double[] LatLon1, double[] LatLon2) {

		double lat1 = LatLon1[0];
		double lon1 = LatLon1[1];
		double lat2 = LatLon2[0];
		double lon2 = LatLon2[1];
		// The math module contains a function named toRadians which converts from
		// degrees to radians.
		lon1 = Math.toRadians(lon1);
		lon2 = Math.toRadians(lon2);
		lat1 = Math.toRadians(lat1);
		lat2 = Math.toRadians(lat2);

		// Haversine formula
		double dlon = lon2 - lon1;
		double dlat = lat2 - lat1;
		double a = Math.pow(Math.sin(dlat / 2), 2) + Math.cos(lat1) * Math.cos(lat2) * Math.pow(Math.sin(dlon / 2), 2);

		double c = 2 * Math.asin(Math.sqrt(a));

		// Used Mean Earth Radius
		double r = 6371000;

		// calculate the result
		return (c * r);
	}

	public static double getVincentyDistance(double[] LatLon1, double[] LatLon2) {
		double lat1 = LatLon1[0];
		double lon1 = LatLon1[1];
		double lat2 = LatLon2[0];
		double lon2 = LatLon2[1];

		double L = Math.toRadians(lon2 - lon1);

		double U1 = Math.atan((1 - f) * Math.tan(Math.toRadians(lat1)));

		double U2 = Math.atan((1 - f) * Math.tan(Math.toRadians(lat2)));

		double sinU1 = Math.sin(U1), cosU1 = Math.cos(U1);

		double sinU2 = Math.sin(U2), cosU2 = Math.cos(U2);

		double cosSqAlpha, sinSigma, cos2SigmaM, cosSigma, sigma;

		double lambda = L, lambdaP, iterLimit = 100;

		do {

			double sinLambda = Math.sin(lambda), cosLambda = Math.cos(lambda);

			sinSigma = Math.sqrt((cosU2 * sinLambda)

					* (cosU2 * sinLambda)

					+ (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda)

							* (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda)

			);

			if (sinSigma == 0)
				return 0;

			cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda;

			sigma = Math.atan2(sinSigma, cosSigma);

			double sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma;

			cosSqAlpha = 1 - sinAlpha * sinAlpha;

			cos2SigmaM = cosSigma - 2 * sinU1 * sinU2 / cosSqAlpha;

			double C = f / 16 * cosSqAlpha * (4 + f * (4 - 3 * cosSqAlpha));

			lambdaP = lambda;

			lambda = L + (1 - C) * f * sinAlpha
					* (sigma + C * sinSigma * (cos2SigmaM + C * cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM)));

		} while (Math.abs(lambda - lambdaP) > 1e-12 && --iterLimit > 0);

		if (iterLimit == 0)
			return 0;

		double uSq = cosSqAlpha * (a * a - b * b) / (b * b);

		double A = 1 + uSq / 16384

				* (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)));

		double B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)));

		double deltaSigma =

				B * sinSigma

						* (cos2SigmaM + B / 4

								* (cosSigma

										* (-1 + 2 * cos2SigmaM * cos2SigmaM)
										- B / 6 * cos2SigmaM

												* (-3 + 4 * sinSigma * sinSigma)

												* (-3 + 4 * cos2SigmaM * cos2SigmaM)));

		double s = b * A * (sigma - deltaSigma);
		return s;
	}

	public static double[] enu2ecef(double[] enu, double[] ECEFr) {
		double[] LLH = ecef2lla(ECEFr);
		double lat = LLH[0];
		double lon = LLH[1];
		double[] ECEF = new double[3];
		ECEF[0] = (-Math.sin(lon) * enu[0]) + (-Math.sin(lat) * Math.cos(lon) * enu[1])
				+ (Math.cos(lat) * Math.cos(lon) * enu[2]) + ECEFr[0];
		ECEF[1] = (Math.cos(lon) * enu[0]) + (-Math.sin(lat) * Math.sin(lon) * enu[1])
				+ (Math.cos(lat) * Math.sin(lon) * enu[2]) + ECEFr[1];
		ECEF[2] = (0 * enu[0]) + (Math.cos(lat) * enu[1]) + (Math.sin(lat) * enu[2]) + ECEFr[2];
		return ECEF;
	}

	public static double[] ned2ecef(double[] ned, double[] ECEFr) {
		double[] ECEF = enu2ecef(enu_ned_convert(ned), ECEFr);
		return ECEF;
	}

	public static double[] ecef2ned(double[] ecef, double[] refEcef) {
		double[] ned = enu_ned_convert(ecef2enu(ecef, refEcef));
		return ned;
	}

	// Geodetic Latitude to Geocentric Latitude
	public static double gd2gc(double gdLat, double gdAlt) {

		gdLat = Math.toRadians(gdLat);

		double N = a / Math.sqrt(1 - (e2 * Math.pow(Math.sin(gdLat), 2)));
		double gcLat = Math.atan((1 - (e2 * (N / (N + gdAlt)))) * Math.tan(gdLat));
		gcLat = Math.toDegrees(gcLat);
		return gcLat;

	}

	// Reference -
	// https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#From_ECEF_to_ENU
	public static double[] ecef2enu(double[] ecef, double[] refEcef) {
		double[] _diff = new double[] { ecef[0] - refEcef[0], ecef[1] - refEcef[1], ecef[2] - refEcef[2] };
		SimpleMatrix diff = new SimpleMatrix(3, 1, false, _diff);
		double[] llh = ecef2lla(refEcef);
		double lat = Math.toRadians(llh[0]);
		double lon = Math.toRadians(llh[1]);
		double[][] _R = new double[][] { { -Math.sin(lon), Math.cos(lon), 0 },
				{ -Math.sin(lat) * Math.cos(lon), -Math.sin(lat) * Math.sin(lon), Math.cos(lat) },
				{ Math.cos(lat) * Math.cos(lon), Math.cos(lat) * Math.sin(lon), Math.sin(lat) } };
		SimpleMatrix R = new SimpleMatrix(_R);
		SimpleMatrix _enu = R.mult(diff);
		double[] enu = new double[] { _enu.get(0), _enu.get(1), _enu.get(2) };
		return enu;
	}

	public static double[] enu_ned_convert(double[] x) {
		SimpleMatrix Y = enu_ned_convert(new SimpleMatrix(3, 1, true, x));
		double[] y = new double[] { Y.get(0), Y.get(1), Y.get(2) };
		return y;
	}

	public static double[][] enu_ned_convert(double[][] x) {
		SimpleMatrix Y = enu_ned_convert(new SimpleMatrix(x));
		double[][] y = Matrix.matrix2Array(Y);
		return y;
	}

	// DCM for ENU to NED and vice versa is same
	public static SimpleMatrix enu_ned_convert(SimpleMatrix x) {

		double[][] c = new double[][] { { 0, 1, 0 }, { 1, 0, 0 }, { 0, 0, -1 } };
		SimpleMatrix C = new SimpleMatrix(c);
		SimpleMatrix y = C.mult(x);
		return y;

	}

	/*
	 * The radius of curvature for east-west motion is known as the transverse
	 * radius of curvature, Re, it is also known as normal radius or prime vertical
	 * radius
	 */
	public static double getNormalEarthRadius(double lat) {
		double Rn = a / Math.sqrt(1 - Math.pow(e * Math.sin(lat), 2));
		return Rn;
	}

	/*
	 * The radius of curvature for north-south motion is known as the meridian
	 * radius
	 */
	public static double getMeridianEarthRadius(double lat) {
		double Rm = (a * (1 - Math.pow(e, 2))) / Math.pow(1 - Math.pow(e * Math.sin(lat), 2), 1.5);
		return Rm;
	}

	/*
	 * Reference: Principles of GNSS, Inertial, and Multi-Sensor Integrated
	 * Navigation Systems (GNSS Technology and Applications) by Paul D. Groves
	 * Section: 2.3.5 Specific Force, Gravitation, and Gravity
	 */
	public static double getGravity(double lat, double alt) {
		double g0 = (9.7803253359 * (1 + (0.001931853 * Math.pow(Math.sin(lat), 2))))
				/ Math.sqrt(1 - Math.pow(e * Math.sin(lat), 2));
		double g = g0
				* (1 - ((2 * alt / a) * (1 + f + (Math.pow(omega_ie * a, 2) * b / mu))) + (3 * Math.pow(alt / a, 2)));
		return g;
	}
}
