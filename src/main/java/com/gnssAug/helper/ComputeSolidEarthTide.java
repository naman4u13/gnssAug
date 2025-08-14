package com.gnssAug.helper;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.orekit.data.BodiesElements;
import org.orekit.data.DataContext;
import org.orekit.data.DirectoryCrawler;
import org.orekit.data.FundamentalNutationArguments;
import org.orekit.frames.Frame;
import org.orekit.frames.FramesFactory;
import org.orekit.models.earth.displacement.StationDisplacement;
import org.orekit.models.earth.displacement.TidalDisplacement;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.TimeScalesFactory;
import org.orekit.utils.Constants;
import org.orekit.utils.IERSConventions;
import org.orekit.utils.PVCoordinatesProvider;
import org.orekit.bodies.CelestialBodyFactory; // Added import for Sun/Moon providers

import com.gnssAug.utility.LatLonUtil;

import java.io.File;
import java.time.ZonedDateTime;

public class ComputeSolidEarthTide implements StationDisplacement {

	private static TidalDisplacement tideModel;
	private static FundamentalNutationArguments fna;
	static {
//        // Configure Orekit data (adjust path to your orekit-data directory)
//    	File orekitData = new File("/Users/naman.agarwal/Documents/orekit/orekit-data-master/orekit-data-master");
//        DataContext.getDefault().getDataProvidersManager().addProvider(new DirectoryCrawler(orekitData));

		PVCoordinatesProvider sun = CelestialBodyFactory.getSun();
		PVCoordinatesProvider moon = CelestialBodyFactory.getMoon();
		// Use IERS_2010 for high precision; false for tide-free system (add permanent
		// if needed)
		tideModel = new TidalDisplacement(Constants.EGM96_EARTH_EQUATORIAL_RADIUS,
				Constants.JPL_SSD_SUN_EARTH_PLUS_MOON_MASS_RATIO, Constants.JPL_SSD_EARTH_MOON_MASS_RATIO, sun, moon,
				IERSConventions.IERS_2010, false);
		fna = IERSConventions.IERS_2010.getNutationArguments(TimeScalesFactory.getTT());
	}

	/**
	 * Calculates the solid Earth tide displacement vector in the ECEF frame using
	 * Orekit for high precision. This yields the displacement to be ADDED to a
	 * station's position to get the instantaneous "tide free" position.
	 *
	 * @param stationCoords The geodetic coordinates of the station [lat (degrees),
	 *                      lon (degrees), alt (m)].
	 * @param dateTime      The time of observation.
	 * @return The displacement vector in ECEF coordinates (meters).
	 */
	public static double[] calculateTimeVaryingTides(double[] stationCoords, boolean isLLA, ZonedDateTime dateTime) {
		AbsoluteDate date = new AbsoluteDate(dateTime.getYear(), dateTime.getMonthValue(), dateTime.getDayOfMonth(),
				dateTime.getHour(), dateTime.getMinute(), (double) dateTime.getSecond(), TimeScalesFactory.getUTC());

		Frame earthFrame = FramesFactory.getITRF(IERSConventions.IERS_2010, true); // With EOP for precision

		double[] rStation = stationCoords;
		if (isLLA) {
			rStation = LatLonUtil.lla2ecef(stationCoords, true);
		}
		Vector3D referencePoint = new Vector3D(rStation[0], rStation[1], rStation[2]);
		BodiesElements elements = fna.evaluateAll(date);
		Vector3D displacement = tideModel.displacement(elements, earthFrame, referencePoint);

		return new double[] { displacement.getX(), displacement.getY(), displacement.getZ() };
	}

	/**
	 * Calculates the correction to transform a "tide free" position to a "mean
	 * tide" position.
	 * 
	 * @param stationCoords The geodetic coordinates of the station [lat (radians),
	 *                      lon (radians), alt (m)].
	 * @return The correction vector in ECEF coordinates (meters).
	 */
	public static double[] getMeanTideCorrection(double[] stationCoords) {
		// Fixed bug: Use geodetic directly (assume lat/lon in radians)

		double[] llh =  LatLonUtil.ecef2lla(stationCoords,false);
		double sinLat = Math.sin(llh[0]);

		double p2 = 0.5 * ((3 * Math.pow(sinLat, 2)) - 1);

		double dUp = (-0.1206 + (0.0001 * p2)) * p2;
		double dNorth = (-0.0252 - (0.0001 * p2)) * 0.5 * Math.sin(2 * llh[0]);

		double[] dEnu = new double[] { 0, dNorth, dUp }; // East, North, Up

		// ENU to ECEF (transpose of ECEF to ENU)
		return LatLonUtil.enu2ecef(dEnu, stationCoords, false);
	}

	@Override
	public Vector3D displacement(BodiesElements elements, Frame earthFrame, Vector3D referencePoint) {
		// TODO Auto-generated method stub
		return null;
	}
}