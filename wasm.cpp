#ifndef EM_PORT_API
#if defined(__EMSCRIPTEN__)
#include <emscripten.h>
#if defined(__cplusplus)
#define EM_PORT_API(rettype) extern "C" rettype EMSCRIPTEN_KEEPALIVE
#else
#define EM_PORT_API(rettype) rettype EMSCRIPTEN_KEEPALIVE
#endif
#else
#if defined(__cplusplus)
#define EM_PORT_API(rettype) extern "C" rettype
#else
#define EM_PORT_API(rettype) rettype
#endif
#endif
#endif

#include "gpssim.h"
#include <time.h>
#include <iostream>

bool readEphemeris(gpstime_t &g0, bool timeoverwrite)
{
    ////////////////////////////////////////////////////////////
	// Read ephemeris
	////////////////////////////////////////////////////////////
    ephem_t eph[EPHEM_ARRAY_SIZE][MAX_SAT];
    ionoutc_t ionoutc;
    ionoutc.enable = TRUE;
	int neph = readRinexNavAll(eph, &ionoutc, "brdc");
    int verb = FALSE;
    gpstime_t gmin,gmax;
    datetime_t tmin,tmax;
    int sv, i;

	if (neph==0)
	{
		fprintf(stderr, "ERROR: No ephemeris available.\n");
		return false;
	}
	else if (neph==-1)
	{
		fprintf(stderr, "ERROR: ephemeris file not found.\n");
		return false;
	}

	if ((verb==TRUE)&&(ionoutc.vflg==TRUE))
	{
		fprintf(stderr, "  %12.3e %12.3e %12.3e %12.3e\n", 
			ionoutc.alpha0, ionoutc.alpha1, ionoutc.alpha2, ionoutc.alpha3);
		fprintf(stderr, "  %12.3e %12.3e %12.3e %12.3e\n", 
			ionoutc.beta0, ionoutc.beta1, ionoutc.beta2, ionoutc.beta3);
		fprintf(stderr, "   %19.11e %19.11e  %9d %9d\n",
			ionoutc.A0, ionoutc.A1, ionoutc.tot, ionoutc.wnt);
		fprintf(stderr, "%6d\n", ionoutc.dtls);
	}

	for (sv=0; sv<MAX_SAT; sv++) 
	{
		if (eph[0][sv].vflg==1)
		{
			gmin = eph[0][sv].toc;
			tmin = eph[0][sv].t;
			break;
		}
	}

	gmax.sec = 0;
	gmax.week = 0;
	tmax.sec = 0;
	tmax.mm = 0;
	tmax.hh = 0;
	tmax.d = 0;
	tmax.m = 0;
	tmax.y = 0;
	for (sv=0; sv<MAX_SAT; sv++)
	{
		if (eph[neph-1][sv].vflg == 1)
		{
			gmax = eph[neph-1][sv].toc;
			tmax = eph[neph-1][sv].t;
			break;
		}
	}

	if (g0.week>=0) // Scenario start time has been set.
	{
		if (timeoverwrite==true)
		{
			gpstime_t gtmp;
			datetime_t ttmp;
			double dsec;

			gtmp.week = g0.week;
			gtmp.sec = (double)(((int)(g0.sec))/7200)*7200.0;

			dsec = subGpsTime(gtmp,gmin);

			// Overwrite the UTC reference week number
			ionoutc.wnt = gtmp.week;
			ionoutc.tot = (int)gtmp.sec;

			// Iono/UTC parameters may no longer valid
			//ionoutc.vflg = FALSE;

			// Overwrite the TOC and TOE to the scenario start time
			for (sv=0; sv<MAX_SAT; sv++)
			{
				for (i=0; i<neph; i++)
				{
					if (eph[i][sv].vflg == 1)
					{
						gtmp = incGpsTime(eph[i][sv].toc, dsec);
						gps2date(&gtmp,&ttmp);
						eph[i][sv].toc = gtmp;
						eph[i][sv].t = ttmp;

						gtmp = incGpsTime(eph[i][sv].toe, dsec);
						eph[i][sv].toe = gtmp;
					}
				}
			}
		}
		else
		{
			if (subGpsTime(g0, gmin)<0.0 || subGpsTime(gmax, g0)<0.0)
			{
				fprintf(stderr, "ERROR: Invalid start time.\n");
				fprintf(stderr, "tmin = %4d/%02d/%02d,%02d:%02d:%02.0f (%d:%.0f)\n", 
					tmin.y, tmin.m, tmin.d, tmin.hh, tmin.mm, tmin.sec,
					gmin.week, gmin.sec);
				fprintf(stderr, "tmax = %4d/%02d/%02d,%02d:%02d:%02.0f (%d:%.0f)\n", 
					tmax.y, tmax.m, tmax.d, tmax.hh, tmax.mm, tmax.sec,
					gmax.week, gmax.sec);
				exit(1);
			}
		}
	}
	else
	{
		g0 = gmin;
		t0 = tmin;
	}

	fprintf(stderr, "Start time = %4d/%02d/%02d,%02d:%02d:%02.0f (%d:%.0f)\n", 
		t0.y, t0.m, t0.d, t0.hh, t0.mm, t0.sec, g0.week, g0.sec);
	fprintf(stderr, "Duration = %.1f [sec]\n", ((double)numd)/10.0);

	// Select the current set of ephemerides
	ieph = -1;

	for (i=0; i<neph; i++)
	{
		for (sv=0; sv<MAX_SAT; sv++)
		{
			if (eph[i][sv].vflg == 1)
			{
				dt = subGpsTime(g0, eph[i][sv].toc);
				if (dt>=-SECONDS_IN_HOUR && dt<SECONDS_IN_HOUR)
				{
					ieph = i;
					break;
				}
			}
		}

		if (ieph>=0) // ieph has been set
			break;
	}

	if (ieph == -1)
	{
		fprintf(stderr, "ERROR: No current set of ephemerides has been found.\n");
		exit(1);
	}
}

EM_PORT_API(bool)
staticMode(double lat, double lng, double height, int unix_time, double duration, int data_format, double samp_freq, bool timeoverwrite)
{
    gpstime_t g0;
    g0.week = -1;

    datetime_t t0;

    double llh[3];
    double xyz[USER_MOTION_SIZE][3];
    int iduration = USER_MOTION_SIZE;
    if (duration == 0.0)
    {
        duration = (double)iduration / 10.0;
    }
    if (samp_freq == 0.0)
    {
        samp_freq = 2.6e6;
    }
    int iq_buff_size = (int)samp_freq;

    llh[0] = lat;
    llh[1] = lng;
    llh[2] = height;

    if (data_format != SC01 && data_format != SC08 && data_format != SC16)
    {
        std::cout << "ERROR: Invalid I/Q data format.(1/8/16)" << std::endl;
        return false;
    }
    if (samp_freq < 1.0e6)
    {
        std::cout << "ERROR: Invalid sampling frequency." << std::endl;
        return false;
    }
    if (timeoverwrite)
    {
        time_t timer;
        struct tm *gmt;

        time(&timer);
        gmt = gmtime(&timer);

        t0.y = gmt->tm_year + 1900;
        t0.m = gmt->tm_mon + 1;
        t0.d = gmt->tm_mday;
        t0.hh = gmt->tm_hour;
        t0.mm = gmt->tm_min;
        t0.sec = (double)gmt->tm_sec;

        date2gps(&t0, &g0);
    }
    else
    {
        if (unix_time)
        {
            time_t time = static_cast<time_t>(unix_time);
            std::tm *ptm = gmtime(&time);

            t0.y = ptm->tm_year + 1900;
            t0.m = ptm->tm_mon + 1;
            t0.d = ptm->tm_mday;
            t0.hh = ptm->tm_hour;
            t0.mm = ptm->tm_min;
            t0.sec = (double)ptm->tm_sec;
            if (t0.y <= 1980 || t0.m < 1 || t0.m > 12 || t0.d < 1 || t0.d > 31 ||
                t0.hh < 0 || t0.hh > 23 || t0.mm < 0 || t0.mm > 59 || t0.sec < 0.0 || t0.sec >= 60.0)
            {
                std::cout << "ERROR: Invalid date and time." << std::endl;
                return false;
            }
            t0.sec = floor(t0.sec);
            date2gps(&t0, &g0);
            llh[0] = llh[0] / R2D; // convert to RAD
            llh[1] = llh[1] / R2D; // convert to RAD
            llh2xyz(llh, xyz[0]);  // Convert llh to xyz
        }
    }
    if (duration < 0.0 || (duration > STATIC_MAX_DURATION))
    {
        std::cout << "ERROR: Invalid duration." << std::endl;
        return false;
    }
    iduration = (int)(duration * 10.0 + 0.5);

    // Buffer size
    samp_freq = floor(samp_freq / 10.0);
    iq_buff_size = (int)samp_freq; // samples per 0.1sec
    samp_freq *= 10.0;

    double delt = 1.0 / samp_freq;

    int numd = iduration;

    // Set user initial position
    llh2xyz(llh, xyz[0]);

    fprintf(stderr, "xyz = %11.1f, %11.1f, %11.1f\n", xyz[0][0], xyz[0][1], xyz[0][2]);
	fprintf(stderr, "llh = %11.6f, %11.6f, %11.1f\n", llh[0]*R2D, llh[1]*R2D, llh[2]);
}