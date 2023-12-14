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

template <size_t R, size_t C>
bool readEphemeris(gpstime_t &g0, datetime_t &t0, bool timeoverwrite, int &numd, double &samp_freq, int &data_format, double (&xyz)[R][C], bool staticLocationMode)
{
	////////////////////////////////////////////////////////////
	// Read ephemeris
	////////////////////////////////////////////////////////////
	ephem_t eph[EPHEM_ARRAY_SIZE][MAX_SAT];
	ionoutc_t ionoutc;
	ionoutc.enable = TRUE;
	int neph = readRinexNavAll(eph, &ionoutc, "brdc");
	int verb = FALSE;
	gpstime_t gmin, gmax;
	datetime_t tmin, tmax;
	int sv, i;
	double dt;
	int iTable;
	double path_loss;
	int ibs; // boresight angle index
	double ant_gain;
	int path_loss_enable = TRUE;
	int isamp;
	int gain[MAX_CHAN];
	int ip,qp;
	clock_t tend;
	char outfile[MAX_CHAR];
	strcpy(outfile, "gpssim.bin");

	// Buffer size
	samp_freq = floor(samp_freq / 10.0);
	int iq_buff_size = (int)samp_freq; // samples per 0.1sec
	samp_freq *= 10.0;

	double delt = 1.0 / samp_freq;

	if (neph == 0)
	{
		fprintf(stderr, "ERROR: No ephemeris available.\n");
		return false;
	}
	else if (neph == -1)
	{
		fprintf(stderr, "ERROR: ephemeris file not found.\n");
		return false;
	}

	if ((verb == TRUE) && (ionoutc.vflg == TRUE))
	{
		fprintf(stderr, "  %12.3e %12.3e %12.3e %12.3e\n",
				ionoutc.alpha0, ionoutc.alpha1, ionoutc.alpha2, ionoutc.alpha3);
		fprintf(stderr, "  %12.3e %12.3e %12.3e %12.3e\n",
				ionoutc.beta0, ionoutc.beta1, ionoutc.beta2, ionoutc.beta3);
		fprintf(stderr, "   %19.11e %19.11e  %9d %9d\n",
				ionoutc.A0, ionoutc.A1, ionoutc.tot, ionoutc.wnt);
		fprintf(stderr, "%6d\n", ionoutc.dtls);
	}

	for (sv = 0; sv < MAX_SAT; sv++)
	{
		if (eph[0][sv].vflg == 1)
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
	for (sv = 0; sv < MAX_SAT; sv++)
	{
		if (eph[neph - 1][sv].vflg == 1)
		{
			gmax = eph[neph - 1][sv].toc;
			tmax = eph[neph - 1][sv].t;
			break;
		}
	}

	if (g0.week >= 0) // Scenario start time has been set.
	{
		if (timeoverwrite == true)
		{
			gpstime_t gtmp;
			datetime_t ttmp;
			double dsec;

			gtmp.week = g0.week;
			gtmp.sec = (double)(((int)(g0.sec)) / 7200) * 7200.0;

			dsec = subGpsTime(gtmp, gmin);

			// Overwrite the UTC reference week number
			ionoutc.wnt = gtmp.week;
			ionoutc.tot = (int)gtmp.sec;

			// Iono/UTC parameters may no longer valid
			// ionoutc.vflg = FALSE;

			// Overwrite the TOC and TOE to the scenario start time
			for (sv = 0; sv < MAX_SAT; sv++)
			{
				for (i = 0; i < neph; i++)
				{
					if (eph[i][sv].vflg == 1)
					{
						gtmp = incGpsTime(eph[i][sv].toc, dsec);
						gps2date(&gtmp, &ttmp);
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
			if (subGpsTime(g0, gmin) < 0.0 || subGpsTime(gmax, g0) < 0.0)
			{
				fprintf(stderr, "ERROR: Invalid start time.\n");
				fprintf(stderr, "tmin = %4d/%02d/%02d,%02d:%02d:%02.0f (%d:%.0f)\n",
						tmin.y, tmin.m, tmin.d, tmin.hh, tmin.mm, tmin.sec,
						gmin.week, gmin.sec);
				fprintf(stderr, "tmax = %4d/%02d/%02d,%02d:%02d:%02.0f (%d:%.0f)\n",
						tmax.y, tmax.m, tmax.d, tmax.hh, tmax.mm, tmax.sec,
						gmax.week, gmax.sec);
				return false;
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
	fprintf(stderr, "Duration = %.1f [sec]\n", ((double)numd) / 10.0);

	// Select the current set of ephemerides
	int ieph = -1;

	for (i = 0; i < neph; i++)
	{
		for (sv = 0; sv < MAX_SAT; sv++)
		{
			if (eph[i][sv].vflg == 1)
			{
				dt = subGpsTime(g0, eph[i][sv].toc);
				if (dt >= -SECONDS_IN_HOUR && dt < SECONDS_IN_HOUR)
				{
					ieph = i;
					break;
				}
			}
		}

		if (ieph >= 0) // ieph has been set
			break;
	}

	if (ieph == -1)
	{
		fprintf(stderr, "ERROR: No current set of ephemerides has been found.\n");
		return false;
	}

	////////////////////////////////////////////////////////////
	// Baseband signal buffer and output file
	////////////////////////////////////////////////////////////

	// Allocate I/Q buffer
	short *iq_buff = static_cast<short *>(calloc(2 * iq_buff_size, 2));
	signed char *iq8_buff = NULL;

	if (iq_buff == NULL)
	{
		fprintf(stderr, "ERROR: Failed to allocate 16-bit I/Q buffer.\n");
		return false;
	}

	if (data_format == SC08)
	{
		iq8_buff = static_cast<signed char *>(calloc(2 * iq_buff_size, 1));
		if (iq8_buff == NULL)
		{
			fprintf(stderr, "ERROR: Failed to allocate 8-bit I/Q buffer.\n");
			return false;
		}
	}
	else if (data_format == SC01)
	{
		iq8_buff = static_cast<signed char *>(calloc(iq_buff_size / 4, 1)); // byte = {I0, Q0, I1, Q1, I2, Q2, I3, Q3}
		if (iq8_buff == NULL)
		{
			fprintf(stderr, "ERROR: Failed to allocate compressed 1-bit I/Q buffer.\n");
			return false;
		}
	}

	// Open output file
	// "-" can be used as name for stdout
	FILE *fp;
	if (strcmp("-", outfile))
	{
		if (NULL == (fp = fopen(outfile, "wb")))
		{
			fprintf(stderr, "ERROR: Failed to open output file.\n");
			exit(1);
		}
	}
	else
	{
		fp = stdout;
	}

	////////////////////////////////////////////////////////////
	// Initialize channels
	////////////////////////////////////////////////////////////
	channel_t chan[MAX_CHAN];
	gpstime_t grx;
	double elvmask = 0.0; // in degree
	double ant_pat[37];
	clock_t tstart;
	int iumd;

	// Clear all channels
	for (i = 0; i < MAX_CHAN; i++)
		chan[i].prn = 0;

	// Clear satellite allocation flag
	for (sv = 0; sv < MAX_SAT; sv++)
		allocatedSat[sv] = -1;

	// Initial reception time
	grx = incGpsTime(g0, 0.0);

	// Allocate visible satellites
	allocateChannel(chan, eph[ieph], ionoutc, grx, xyz[0], elvmask);

	for (i = 0; i < MAX_CHAN; i++)
	{
		if (chan[i].prn > 0)
			fprintf(stderr, "%02d %6.1f %5.1f %11.1f %5.1f\n", chan[i].prn,
					chan[i].azel[0] * R2D, chan[i].azel[1] * R2D, chan[i].rho0.d, chan[i].rho0.iono_delay);
	}

	////////////////////////////////////////////////////////////
	// Receiver antenna gain pattern
	////////////////////////////////////////////////////////////

	for (i = 0; i < 37; i++)
		ant_pat[i] = pow(10.0, -ant_pat_db[i] / 20.0);

	////////////////////////////////////////////////////////////
	// Generate baseband signals
	////////////////////////////////////////////////////////////

	tstart = clock();

	// Update receiver time
	grx = incGpsTime(grx, 0.1);

	for (iumd = 1; iumd < numd; iumd++)
	{
		for (i = 0; i < MAX_CHAN; i++)
		{
			if (chan[i].prn > 0)
			{
				// Refresh code phase and data bit counters
				range_t rho;
				sv = chan[i].prn - 1;

				// Current pseudorange
				if (!staticLocationMode)
					computeRange(&rho, eph[ieph][sv], &ionoutc, grx, xyz[iumd]);
				else
					computeRange(&rho, eph[ieph][sv], &ionoutc, grx, xyz[0]);

				chan[i].azel[0] = rho.azel[0];
				chan[i].azel[1] = rho.azel[1];

				// Update code phase and data bit counters
				computeCodePhase(&chan[i], rho, 0.1);
#ifndef FLOAT_CARR_PHASE
				chan[i].carr_phasestep = (int)round(512.0 * 65536.0 * chan[i].f_carr * delt);
#endif
				// Path loss
				path_loss = 20200000.0 / rho.d;

				// Receiver antenna gain
				ibs = (int)((90.0 - rho.azel[1] * R2D) / 5.0); // covert elevation to boresight
				ant_gain = ant_pat[ibs];

				// Signal gain
				if (path_loss_enable == TRUE)
					gain[i] = (int)(path_loss * ant_gain * 128.0); // scaled by 2^7
				else
					gain[i] = 128; // hold the power level constant
			}
		}

		for (isamp = 0; isamp < iq_buff_size; isamp++)
		{
			int i_acc = 0;
			int q_acc = 0;

			for (i = 0; i < MAX_CHAN; i++)
			{
				if (chan[i].prn > 0)
				{
#ifdef FLOAT_CARR_PHASE
					iTable = (int)floor(chan[i].carr_phase * 512.0);
#else
					iTable = (chan[i].carr_phase >> 16) & 0x1ff; // 9-bit index
#endif
					ip = chan[i].dataBit * chan[i].codeCA * cosTable512[iTable] * gain[i];
					qp = chan[i].dataBit * chan[i].codeCA * sinTable512[iTable] * gain[i];

					// Accumulate for all visible satellites
					i_acc += ip;
					q_acc += qp;

					// Update code phase
					chan[i].code_phase += chan[i].f_code * delt;

					if (chan[i].code_phase >= CA_SEQ_LEN)
					{
						chan[i].code_phase -= CA_SEQ_LEN;

						chan[i].icode++;

						if (chan[i].icode >= 20) // 20 C/A codes = 1 navigation data bit
						{
							chan[i].icode = 0;
							chan[i].ibit++;

							if (chan[i].ibit >= 30) // 30 navigation data bits = 1 word
							{
								chan[i].ibit = 0;
								chan[i].iword++;
								/*
								if (chan[i].iword>=N_DWRD)
									fprintf(stderr, "\nWARNING: Subframe word buffer overflow.\n");
								*/
							}

							// Set new navigation data bit
							chan[i].dataBit = (int)((chan[i].dwrd[chan[i].iword] >> (29 - chan[i].ibit)) & 0x1UL) * 2 - 1;
						}
					}

					// Set current code chip
					chan[i].codeCA = chan[i].ca[(int)chan[i].code_phase] * 2 - 1;

					// Update carrier phase
#ifdef FLOAT_CARR_PHASE
					chan[i].carr_phase += chan[i].f_carr * delt;

					if (chan[i].carr_phase >= 1.0)
						chan[i].carr_phase -= 1.0;
					else if (chan[i].carr_phase < 0.0)
						chan[i].carr_phase += 1.0;
#else
					chan[i].carr_phase += chan[i].carr_phasestep;
#endif
				}
			}

			// Scaled by 2^7
			i_acc = (i_acc + 64) >> 7;
			q_acc = (q_acc + 64) >> 7;

			// Store I/Q samples into buffer
			iq_buff[isamp * 2] = (short)i_acc;
			iq_buff[isamp * 2 + 1] = (short)q_acc;
		}

		if (data_format == SC01)
		{
			for (isamp = 0; isamp < 2 * iq_buff_size; isamp++)
			{
				if (isamp % 8 == 0)
					iq8_buff[isamp / 8] = 0x00;

				iq8_buff[isamp / 8] |= (iq_buff[isamp] > 0 ? 0x01 : 0x00) << (7 - isamp % 8);
			}

			fwrite(iq8_buff, 1, iq_buff_size / 4, fp);
		}
		else if (data_format == SC08)
		{
			for (isamp = 0; isamp < 2 * iq_buff_size; isamp++)
				iq8_buff[isamp] = iq_buff[isamp] >> 4; // 12-bit bladeRF -> 8-bit HackRF

			fwrite(iq8_buff, 1, 2 * iq_buff_size, fp);
		}
		else // data_format==SC16
		{
			fwrite(iq_buff, 2, 2 * iq_buff_size, fp);
		}

		//
		// Update navigation message and channel allocation every 30 seconds
		//

		int igrx = (int)(grx.sec * 10.0 + 0.5);

		if (igrx % 300 == 0) // Every 30 seconds
		{
			// Update navigation message
			for (i = 0; i < MAX_CHAN; i++)
			{
				if (chan[i].prn > 0)
					generateNavMsg(grx, &chan[i], 0);
			}

			// Refresh ephemeris and subframes
			// Quick and dirty fix. Need more elegant way.
			for (sv = 0; sv < MAX_SAT; sv++)
			{
				if (eph[ieph + 1][sv].vflg == 1)
				{
					dt = subGpsTime(eph[ieph + 1][sv].toc, grx);
					if (dt < SECONDS_IN_HOUR)
					{
						ieph++;

						for (i = 0; i < MAX_CHAN; i++)
						{
							// Generate new subframes if allocated
							if (chan[i].prn != 0)
								eph2sbf(eph[ieph][chan[i].prn - 1], ionoutc, chan[i].sbf);
						}
					}

					break;
				}
			}

			// Update channel allocation
			if (!staticLocationMode)
				allocateChannel(chan, eph[ieph], ionoutc, grx, xyz[iumd], elvmask);
			else
				allocateChannel(chan, eph[ieph], ionoutc, grx, xyz[0], elvmask);

			// Show details about simulated channels
			if (verb == TRUE)
			{
				fprintf(stderr, "\n");
				for (i = 0; i < MAX_CHAN; i++)
				{
					if (chan[i].prn > 0)
						fprintf(stderr, "%02d %6.1f %5.1f %11.1f %5.1f\n", chan[i].prn,
								chan[i].azel[0] * R2D, chan[i].azel[1] * R2D, chan[i].rho0.d, chan[i].rho0.iono_delay);
				}
			}
		}

		// Update receiver time
		grx = incGpsTime(grx, 0.1);

		// Update time counter
		fprintf(stderr, "\rTime into run = %4.1f", subGpsTime(grx, g0));
		fflush(stdout);
	}

	tend = clock();

	fprintf(stderr, "\nDone!\n");

	// Free I/Q buffer
	free(iq_buff);

	// Close file
	fclose(fp);

	// Process time
	fprintf(stderr, "Process time = %.1f [sec]\n", (double)(tend - tstart) / CLOCKS_PER_SEC);

	return (0);
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

	int numd = iduration;

	// Set user initial position
	llh2xyz(llh, xyz[0]);

	fprintf(stderr, "xyz = %11.1f, %11.1f, %11.1f\n", xyz[0][0], xyz[0][1], xyz[0][2]);
	fprintf(stderr, "llh = %11.6f, %11.6f, %11.1f\n", llh[0] * R2D, llh[1] * R2D, llh[2]);

	readEphemeris(g0, t0, timeoverwrite, numd, samp_freq, data_format, xyz, true);
	return true;
}