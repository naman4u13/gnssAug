package com.gnssAug.utility;

import java.time.Instant;
import java.time.ZoneId;
import java.time.ZonedDateTime;
import java.time.temporal.ChronoUnit;
import java.util.Calendar;
import java.util.TimeZone;
import java.util.stream.IntStream;

public class Time {
	private static final long NumberMilliSecondsWeek = 604800000;
	private static final long NumberSecondsWeek = 604800;
	// GPS epoch: 1980-01-06 00:00:00 UTC
    private static final Instant GPS_EPOCH = Instant.parse("1980-01-06T00:00:00Z");

    // Current difference between GPS time and UTC in seconds (as of July 2025)
    private static final int LEAP_SECONDS = 18;

    /**
     * Correctly converts GPS week and seconds of the week to a ZonedDateTime.
     * This function accounts for the leap second difference between GPS time and UTC.
     *
     * @param gpsWeek The GPS week number.
     * @param gpsSecondsOfWeek The number of seconds into the GPS week.
     * @return The correct ZonedDateTime in UTC.
     */
    public static ZonedDateTime convertToZonedDateTime(int gpsWeek, double gpsSecondsOfWeek) {
        // Total seconds since GPS epoch on the GPS timescale
        long totalGpsSeconds = (long) ((gpsWeek * 7 * 86400) + gpsSecondsOfWeek);

        // Compute the instant on the continuous GPS timescale
        Instant gpsInstant = GPS_EPOCH.plus(totalGpsSeconds, ChronoUnit.SECONDS);

        // Convert from GPS time to UTC by subtracting the leap seconds
        Instant utcInstant = gpsInstant.minusSeconds(LEAP_SECONDS);

        // Return as ZonedDateTime in UTC
        return ZonedDateTime.ofInstant(utcInstant, ZoneId.of("UTC"));
    }
	
	// Convert GPS timestamp to GPS TOW and Week No.
	// Note this func does not handle conversion of UTC timestamp to GPS time
	// Java does not support handling leap seconds, so Calendar or Date classes will
	// not account for leap seconds
	public static double[] getGPSTime(int year, int month, int day, int hour, int minute, double sec) {

		Calendar cal = Calendar.getInstance(TimeZone.getTimeZone("UTC"));
		cal.set(year, month, day, hour, minute, 0);
		// Nomenclature is unixTime because getTimeInMillis does not account for leap
		// seconds
		double unixTime = cal.getTimeInMillis() + (1000 * sec);
		cal.set(1980, 0, 6, 0, 0, 0);
		long GPSEpoch = cal.getTimeInMillis();
		double GPSTime = ((unixTime - GPSEpoch)) % NumberMilliSecondsWeek;
		double weekNo = Math.floor(((unixTime - GPSEpoch)) / NumberMilliSecondsWeek);

		GPSTime = GPSTime / 1000;
		return new double[] { GPSTime, weekNo };
	}

	// Remember this method will not work if you have milliseconds
	// So handling of millisec outside this func is must
	public static double[] getGPSTime(Calendar time) {

		Calendar cal = Calendar.getInstance(TimeZone.getTimeZone("UTC"));

		double unixTime = time.getTimeInMillis() / 1000;

		cal.set(1980, 0, 6, 0, 0, 0);
		long GPSEpoch = cal.getTimeInMillis() / 1000;
		double GPSTime = ((unixTime - GPSEpoch)) % NumberSecondsWeek;
		long weekNo = (long) (((unixTime - GPSEpoch)) / NumberSecondsWeek);

		return new double[] { GPSTime, weekNo };
	}
	
	public static double UnixTimeToGPSTimeMilli(double unixTime)
	{
		Calendar cal = Calendar.getInstance(TimeZone.getTimeZone("UTC"));
		cal.set(1980, 0, 6, 0, 0, 0);
		double GPSEpoch = cal.getTimeInMillis();
		double GPSTime = unixTime - GPSEpoch;
	

		return GPSTime;
	}

	public static Calendar getDate(double GPSTime, long weekNo, double longitude) {
		// There is no such thing as local GPS Time, the variable is a way of
		// representing
		// GPS time in local TIMEZONE

		long localGPSTime = (long) (4.32 * 1000 * (1E4) * (longitude / 180)
				+ (GPSTime * 1000 + (weekNo * NumberMilliSecondsWeek)));
		Calendar cal = Calendar.getInstance(TimeZone.getTimeZone("UTC"));
		cal.set(1980, 0, 6, 0, 0, 0);
		long GPSEpoch = (cal.getTimeInMillis() / 1000) * 1000;
		// GPSEpoch = Math.round((double) GPSEpoch / 1000) * 1000;
		long unixTime = (localGPSTime) + GPSEpoch;
		cal.setTimeInMillis(unixTime);
		return cal;

	}

	public static Calendar convertToUTC(Calendar localGPStime, double longitude) {

		long UTCtime = (localGPStime.getTimeInMillis()) - (long) (4.32 * (1E4) * (longitude / 180) * 1000);
		Calendar cal = Calendar.getInstance(TimeZone.getTimeZone("UTC"));
		cal.setTimeInMillis(UTCtime);
		return cal;

	}

	public static double[] getGPSTime(String[] strTime) {

		int[] tArr = IntStream.range(0, 5).map(x -> Integer.parseInt(strTime[x])).toArray();
		double sec = Double.parseDouble(strTime[5]);
		double[] GPStime = getGPSTime(tArr[0], tArr[1] - 1, tArr[2], tArr[3], tArr[4], sec);
		return GPStime;
	}
	
	public static ZonedDateTime convertUsingToInstant(Calendar calendar) {
        // Get the Calendar's timezone
        ZoneId zoneId = calendar.getTimeZone().toZoneId();
        
        // Convert to ZonedDateTime
        return calendar.toInstant().atZone(zoneId);
    }
	
	public static double[] getGPSTimeFromYDOY(String[] timeArr) {
		
		return getGPSTimeFromYDOY(Integer.parseInt(timeArr[0]), Integer.parseInt(timeArr[1]), Double.parseDouble(timeArr[2]));
	}
	
	public static double[] getGPSTimeFromYDOY(int year, int dayOfYear, double timeOfDaySeconds) {
		Calendar cal = Calendar.getInstance(TimeZone.getTimeZone("UTC"));
		cal.set(Calendar.YEAR, year);
		cal.set(Calendar.DAY_OF_YEAR, dayOfYear);
		cal.set(Calendar.HOUR_OF_DAY, 0);
		cal.set(Calendar.MINUTE, 0);
		cal.set(Calendar.SECOND, 0);
		cal.set(Calendar.MILLISECOND, 0);

		// Add timeOfDaySeconds to the calendar time
		long unixTime = cal.getTimeInMillis() + (long) (timeOfDaySeconds * 1000);

		cal.set(1980, 0, 6, 0, 0, 0);
		long GPSEpoch = cal.getTimeInMillis();
		double GPSTime = ((unixTime - GPSEpoch)) % NumberMilliSecondsWeek;
		double weekNo = Math.floor(((unixTime - GPSEpoch)) / NumberMilliSecondsWeek);

		GPSTime = GPSTime / 1000;
		return new double[] { GPSTime, weekNo };
	}

}
