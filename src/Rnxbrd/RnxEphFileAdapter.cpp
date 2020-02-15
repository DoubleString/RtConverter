/*
 * RnxEphFileAdapter.cpp
 *
 *  Created on: 2018/3/1
 *      Author: doublestring
 */
#include <iostream>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include "../../include/RtConverter/Deploy.h"
#include "../../include/Rnxbrd/OrbitClk.h"
#include "../../include/Rnxobs/Patterns.h"
#include "../../include/Rnxobs/Com.h"
using namespace std;
using namespace bamboo;
// read from file to acquire ephemeris clock & orbit
RnxEphFileAdapter::~RnxEphFileAdapter() {
	// the stl container will delete automaticly
	v_closeRnxEph();
}
RnxEphFileAdapter::RnxEphFileAdapter() {
	memset(neph, 0, sizeof(int) * MAXSYS);
	mjd0 = mjd1 = leap = 0;
	sod0 = sod1 = ver = 0.0;
	memset(ionc, 0, sizeof(char) * MAXSYS * 2 * LEN_EPHION);
	memset(ion, 0, sizeof(double) * MAXSYS * 2 * 4);
	memset(timc, 0, sizeof(char) * MAXSYS * 2 * LEN_EPHION);
	memset(tim, 0, sizeof(double) * MAXSYS * 2 * 4);
	isOpen = 0;
	orbType = ORB_BRD;
	clkType = CLK_BRD;
}
void RnxEphFileAdapter::v_openRnxEph(string ephcmd) {
	int curmjd;
	this->isOpen = 1;
	Deploy dly = Deploy::s_getConfigures();
	int iyear, imon, iday, ih, imin;
	double dsec, cursod;
	curmjd = atoi(ephcmd.substr(0, index_string(ephcmd.c_str(), ':')).c_str());
	cursod = atof(ephcmd.substr(index_string(ephcmd.c_str(), ':') + 1).c_str());

	mjd0 = curmjd;
	sod0 = 0.0;
	mjd1 = curmjd + 1;
	sod1 = 0.0;
	mjd2date(curmjd, cursod, &iyear, &imon, &iday, &ih, &imin, &dsec);
	memset(curFile, 0, sizeof(char) * 1024);
	strcpy(curFile, dly.ephdir);
	Patterns::s_getInstance()->m_getPatternName(true, "", "", iyear, imon, iday,
			ih, curFile);

	if (curFile[strlen(curFile) - 1] == 'p') {
		this->m_readRnxnav('M');
	} else {
		this->m_readRnxnav('G');
	}
}
void RnxEphFileAdapter::v_closeRnxEph() {
	int isys;
	for (isys = 0; isys < MAXSYS; isys++) {
		this->gnssEph[isys].clear();
	}
	this->glsEph.clear();
	memset(neph, 0, sizeof(int) * MAXSYS);
	this->isOpen = 0;
}
int RnxEphFileAdapter::v_readClk(const char* cprn, int mjd, double sod,
		double* sclk) {
	// using rinex navigation file to acquire satellite clock
	int wk, iode = -1;
	double sow, tgd, dtmin;
	mjd2wksow(mjd, sod, &wk, &sow);
	*sclk = 0.0;
	if (cprn[0] != 'R')
		if (cprn[0] == 'C')
			return this->m_brd2xyz("ynn", cprn, wk, sow - 14.0, NULL, sclk,
					&dtmin, &tgd, &iode);
		else
			return this->m_brd2xyz("ynn", cprn, wk, sow, NULL, sclk, &dtmin,
					&tgd, &iode);
	else
		return this->m_gls2xyz("ynn", cprn, wk, sow, NULL, sclk, &dtmin, &iode);
}
int RnxEphFileAdapter::v_readClk(const char* cprn, int mjd, double sod,
		double* sclk, int* iode) {
	int wk;
	double sow, tgd, dtmin;
	mjd2wksow(mjd, sod, &wk, &sow);
	*sclk = 0.0;
	if (cprn[0] != 'R')
		if (cprn[0] == 'C')
			return this->m_brd2xyz("ynn", cprn, wk, sow - 14.0, NULL, sclk,
					&dtmin, &tgd, iode);
		else
			return this->m_brd2xyz("ynn", cprn, wk, sow, NULL, sclk, &dtmin,
					&tgd, iode);
	else
		return this->m_gls2xyz("ynn", cprn, wk, sow, NULL, sclk, &dtmin, iode);
}
int RnxEphFileAdapter::v_readOrbit(const char* cprn, int mjd, double sod,
		double* xsat) {
	// using rinex navigation file to acquire orbit
	int wk, iode = -1;
	double sow, tgd, dtmin, sclk;
	mjd2wksow(mjd, sod, &wk, &sow);
	if (cprn[0] != 'R')
		if (cprn[0] == 'C')
			return this->m_brd2xyz("yyy", cprn, wk, sow - 14.0, xsat, &sclk,
					&dtmin, &tgd, &iode);
		else
			return this->m_brd2xyz("yyy", cprn, wk, sow, xsat, &sclk, &dtmin,
					&tgd, &iode);
	else
		return this->m_gls2xyz("yyy", cprn, wk, sow, xsat, &sclk, &dtmin, &iode);
}
int RnxEphFileAdapter::v_readOrbit(const char* cprn, int mjd, double sod,
		double* xsat, int* iode) {
	// using rinex navigation file to acquire orbit
	// using rinex navigation file to acquire orbit
	int wk;
	double sow, tgd, dtmin, sclk;
	mjd2wksow(mjd, sod, &wk, &sow);
	if (cprn[0] != 'R')
		if (cprn[0] == 'C')
			return this->m_brd2xyz("yyy", cprn, wk, sow - 14.0, xsat, &sclk,
					&dtmin, &tgd, iode);
		else
			return this->m_brd2xyz("yyy", cprn, wk, sow, xsat, &sclk, &dtmin,
					&tgd, iode);
	else
		return this->m_gls2xyz("yyy", cprn, wk, sow, xsat, &sclk, &dtmin, iode);
}
void RnxEphFileAdapter::m_readRnxnav(char csys) {
	// read rinex navigation file
	int i, k, isys, len, iy, im, id, ih, imin;
	double dsec, dt1, dt0,d_aode;
	char buf[2048], cline[1024], *ptr;
	list<GPSEPH>::iterator gpsItr;
	list<GLSEPH>::iterator glsItr;
	GPSEPH eph;
	GLSEPH ephr;
	ifstream in;
	string line;
	streampos pos;
	bool already;
	in.open(curFile, ios::in | ios::binary);
	if (!in.is_open()) {
		cout << "***ERROR(m_readRnxnav):can't open file " << curFile << endl;
		return;
	}
	while (getline(in, line)) {
		if (strstr(line.c_str(), "END OF HEADER"))
			break;
		if (strstr(line.c_str(), "RINEX VERSION / TYPE") != NULL) {
			sscanf(line.c_str(), "%lf", &ver);
			if (ver > 3.0) {
				if (csys == 'M' && line[40] != 'M') {
					printf(
							"###WARNING(read_rnxnav):there only broadcast for line[40]!\n");
				}
			}
		} else if (strstr(line.c_str(), "PGM / RUN BY / DATE") != NULL)
			continue;
		else if (strstr(line.c_str(), "COMMENT") != NULL)
			continue;
		else if (strstr(line.c_str(), "ION ALPHA") != NULL) {

		} else if (strstr(line.c_str(), "ION BETA") != NULL) {

		}
		//rinex 3.00 3.01 3.02
		else if (strstr(line.c_str(), "IONOSPHERIC CORR") != NULL) {
			if (!strncmp(line.c_str(), "GPS ", 4)
					|| !strncmp(line.c_str(), "GPSA", 4)) {
				i = index_string(SYS, 'G');
				k = 0;
			} else if (!strncmp(line.c_str(), "GPSB", 4)) {
				i = index_string(SYS, 'G');
				k = 1;
			} else if (!strncmp(line.c_str(), "GAL ", 4)) {
				i = index_string(SYS, 'E');
				k = 0;
			} else if (!strncmp(line.c_str(), "BDS ", 4)
					|| !strncmp(line.c_str(), "BDSA", 4)) {
				i = index_string(SYS, 'C');
				k = 0;
			} else if (!strncmp(line.c_str(), "BDSB", 4)) {
				i = index_string(SYS, 'C');
				k = 1;
			} else if (!strncmp(line.c_str(), "QZS ", 4)
					|| !strncmp(line.c_str(), "QZSA", 4)) {
				i = index_string(SYS, 'J');
				k = 0;
			} else if (!strncmp(line.c_str(), "QZSB", 4)) {
				i = index_string(SYS, 'J');
				k = 1;
			} else
				printf("###WARNING(read_rnxnav):unknown ionospheric corr!\n");
			strncpy(ionc[i][k], line.c_str(), 4);
			strcpy(buf, line.c_str());
			ptr = buf;
			while (*ptr != '\0') {
				if (*ptr == 'D')
					*ptr = 'e';
				ptr++;
			}
			sscanf(buf, "%*s%lf%lf%lf%lf", ion[i][k], ion[i][k] + 1,
					ion[i][k] + 2, ion[i][k] + 3);
		} else if (strstr(line.c_str(), "DELTA-UTC: A0,A1,T,W") != NULL) {
			i = index_string(SYS, csys);
			if (i == -1) {
				printf(
						"$$$MESSAGE(read_rnxnav):DELTA-UTC: A0,A1,T,W is only valid for single system!\n");
				continue;
			}
			strcpy(buf, line.c_str());
			ptr = buf;
			while (*ptr != '\0') {
				if (*ptr == 'D')
					*ptr = 'e';
				ptr++;
			}
			sscanf(buf, "%lf%lf%lf%lf", tim[i][0], tim[i][0] + 1, tim[i][0] + 2,
					tim[i][0] + 3);
		} else if (strstr(line.c_str(), "TIME SYSTEM CORR") != NULL) {
			if (!strncmp(line.c_str(), "GPUT", 4)) {
				i = index_string(SYS, 'G');
				k = 0;
			} else if (!strncmp(line.c_str(), "GLUT", 4)) {
				i = index_string(SYS, 'R');
				k = 0;
			} else if (!strncmp(line.c_str(), "GAUT", 4)) {
				i = index_string(SYS, 'E');
				k = 0;
			} else if (!strncmp(line.c_str(), "BDUT", 4)) {
				i = index_string(SYS, 'C');
				k = 0;
			} else if (!strncmp(line.c_str(), "SBUT", 4)) {
				i = index_string(SYS, 'S');
				k = 0;
			} else if (!strncmp(line.c_str(), "QZUT", 4)) {
				i = index_string(SYS, 'J');
				k = 0;
			} else if (!strncmp(line.c_str(), "GPGA", 4)) {
				i = index_string(SYS, 'G');
				k = 1;
			} else if (!strncmp(line.c_str(), "GLGP", 4)) {
				i = index_string(SYS, 'R');
				k = 1;
			} else if (!strncmp(line.c_str(), "QZGP", 4)) {
				i = index_string(SYS, 'J');
				k = 1;
			} else {
				printf(
						"$$$MESSAGE(read_rnxnav):unknown unknown TIME SYSTEM CORR!\n");
				continue;
			}
			strncpy(timc[i][k], line.c_str(), 4);
			strcpy(buf, line.c_str());
			ptr = buf;
			while (*ptr != '\0') {
				if (*ptr == 'D')
					*ptr = 'e';
				ptr++;
			}
			sscanf(buf, "%*5c%lf%lf%lf%lf", tim[i][k], tim[i][k] + 1,
					tim[i][k] + 2, tim[i][k] + 3);
		} else if (strstr(line.c_str(), "LEAP SECONDS") != NULL)
			sscanf(line.c_str(), "%d", &leap);
	}
	while (!in.eof()) {
		// acquire the ephemeris
		pos = in.tellg();
		getline(in, line);
		if (!len_trim(line.c_str()))
			continue;
		in.seekg(pos,ios::beg);
		if (ver < 3.0) {
			buf[0] = csys;
		} else
			buf[0] = line[0];
		switch (buf[0]) {
		case 'G':
		case 'E':
		case 'C':
		case 'J':
			isys = index_string(SYS, buf[0]);
			len = 0;
			memset(buf, 0, sizeof(char) * 2048);
			for (i = 0; i < 8; i++) {
				getline(in, line);
				strcpy(cline, line.c_str());
				if (i != 7) {
					filleph(cline, ver);
				}
				if (ver < 3) {
					if (buf[0] != csys)
						buf[0] = csys;
					strncpy(buf + 1 + len, cline, strlen(cline));
				} else
					strncpy(buf + len, cline, strlen(cline));
				len += strlen(cline);
			}
			ptr = buf;
			while (*ptr != '\0') {
				if (*ptr == '\n' || *ptr == '\r')
					*ptr = ' ';

				if (*ptr == 'D')
					*ptr = 'e';
				ptr++;
			}
			if (ver < 3.0) {
				buf[0] = csys;
				if (buf[1] == ' ')
					buf[1] = '0';
			}
			sscanf(buf,
					"%s%d%d%d%d%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
					eph.cprn, &iy, &im, &id, &ih, &imin, &dsec, &eph.a0,
					&eph.a1, &eph.a2, &d_aode, &eph.crs, &eph.dn, &eph.m0,
					&eph.cuc, &eph.e, &eph.cus, &eph.roota, &eph.toe, &eph.cic,
					&eph.omega0, &eph.cis, &eph.i0, &eph.crc, &eph.omega,
					&eph.omegadot, &eph.idot, &eph.resvd0, &eph.week,
					&eph.resvd1, &eph.accu, &eph.hlth, &eph.tgd, &eph.aodc);
			/////////////////////////// for BDS III test,must be removed //////////////////////////////////
			if (atoi(eph.cprn + 1) > 14) {
				eph.hlth = 0.0;
			}
			///////////////////////////////////////////////////////////////////////////////////////////////
			if (eph.hlth > 0 || isys == -1)
				continue;
			yr2year(iy);
			if (iy < 2000)
				continue;
			// time of clock
			eph.mjd = md_julday(iy, im, id); // TIME IN BDS TIME AND SHOULD NOT CHANGE IT INTO GPST BECAUSE THERE ARE OTHER TIME TAG IN THE EPHEMERIS
			eph.sod = ih * 3600.0 + imin * 60.0 + dsec;

			// adapter to bds
			if (eph.cprn[0] == 'C') {
				//the week in broadcast file generated by WHU is GPS week,
				//but that in IGS meraged file is BDS week
				mjd2wksow(eph.mjd, eph.sod, &k, &dt1);
				if (k != (int) eph.week)
					eph.week = 1356 + eph.week;
				eph.tgd1 = eph.aodc;
			}
			//check time
			dt0 = 0.0;
			dt1 = 0.0;
			if (mjd0 != 0) {
				dt0 = eph.mjd + eph.sod / 86400.0 - mjd0;
			}
			if (mjd1 != 0) {
				dt1 = eph.mjd + eph.sod / 86400.0 - mjd1;
			}
			if (dt0 < -1.0 / 24.0 || dt1 > 1.0 / 24.0)
				continue;
			already = false;
			for (gpsItr = this->gnssEph[isys].begin();
					gpsItr != this->gnssEph[isys].end(); gpsItr++) {
				if (strstr((*gpsItr).cprn, eph.cprn) && (*gpsItr).mjd == eph.mjd
						&& (*gpsItr).sod == eph.sod) {
					already = true;
					break;
				}
			}
			if (!already) {
				neph[isys] = neph[isys] + 1;
				eph.aode = genAode(buf[0], eph.mjd,eph.sod, eph.toe,static_cast<int>(d_aode),&eph);
				this->gnssEph[isys].push_back(eph);
			}
			break;
		case 'R':
			isys = index_string(SYS, buf[0]);
			len = 0;
			memset(buf, 0, sizeof(char) * 2048);
			for (i = 0; i < 4; i++) {
				getline(in, line);
				strcpy(cline, line.c_str());
				if (ver < 3) {
					if (buf[0] != csys)
						buf[0] = csys;
					strncpy(buf + 1 + len, cline, strlen(cline));
				} else
					strncpy(buf + len, cline, strlen(cline));
				len += strlen(cline);
			}
			ptr = buf;
			while (*ptr != '\0') {
				if (*ptr == '\n' || *ptr == '\r') {
					*ptr = ' ';
				}
				if (*ptr == 'D') {
					*ptr = 'e';
				}
				ptr++;
			}
			if (ver < 3.0) {
				buf[0] = csys;
				if (buf[1] == ' ')
					buf[1] = '0';
			}
			sscanf(buf,
					"%s%d%d%d%d%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
					ephr.cprn, &iy, &im, &id, &ih, &imin, &dsec, &ephr.tau,
					&ephr.gamma, &ephr.tk, &ephr.pos[0], &ephr.vel[0],
					&ephr.acc[0], &ephr.health, &ephr.pos[1], &ephr.vel[1],
					&ephr.acc[1], &ephr.frenum, &ephr.pos[2], &ephr.vel[2],
					&ephr.acc[2], &ephr.age);

			if (ephr.health > 0.0)
				continue;
			yr2year(iy);
			if (iy < 2000)
				continue;
			ephr.mjd = md_julday(iy, im, id);
			ephr.sod = ih * 3600.0 + imin * 60.0 + dsec;
			//check time
			dt0 = 0.0;
			dt1 = 0.0;
			if (mjd0 != 0) {
				dt0 = ephr.mjd + ephr.sod / 86400.0 - mjd0;
			}
			if (mjd1 != 0) {
				dt1 = ephr.mjd + ephr.sod / 86400.0 - mjd1;
			}
			if (dt0 < -1.0 / 24.0 || dt1 > 1.0 / 24.0)
				continue;
			already = false;
			for (glsItr = this->glsEph.begin(); glsItr != this->glsEph.end();
					glsItr++) {
				if (strstr((*glsItr).cprn, ephr.cprn)
						&& (*glsItr).mjd == ephr.mjd
						&& (*glsItr).sod == ephr.sod) {
					already = true;
					break;
				}
			}
			if (!already) {
				neph[isys] = neph[isys] + 1;
				ephr.aode = genAode(buf[0], ephr.mjd,ephr.sod,0.0, ephr.aode,NULL);
				this->glsEph.push_back(ephr);
			}
			break;
		default:
			getline(in, line);
			break;
		}
	}
	in.close();
}
list<GPSEPH> RnxEphFileAdapter::m_getCurrentEph_G(string cprn, int mjd, double sod) {
	int isys,wk;
	double dt,sow,dtmin = 2.0;
	list<GPSEPH> list_out;
	list<GPSEPH>::iterator gpsItr;
	mjd2wksow(mjd, sod, &wk, &sow);
	isys = index_string(SYS, cprn[0]);
	for (gpsItr = this->gnssEph[isys].begin(); gpsItr != this->gnssEph[isys].end(); ++gpsItr) {
		if (strstr((*gpsItr).cprn, cprn.c_str())) {
			dt = (wk - (*gpsItr).week) * 168.0 + (sow - (*gpsItr).toe) / 3600.0;
			if (fabs(dt) < dtmin) {
				list_out.push_back(*gpsItr);
			}
		}
	}
	return list_out;
}
list<GLSEPH> RnxEphFileAdapter::m_getCurrentEph_R(string cprn, int mjd, double sod) {
	double dt,dtmin = 1.0;
	list<GLSEPH> list_out;
	list<GLSEPH>::iterator glsItr;
	for (glsItr = this->glsEph.begin(); glsItr != this->glsEph.end(); ++glsItr) {
		if (strstr((*glsItr).cprn, cprn.c_str())) {
			dt = (mjd - (*glsItr).mjd) * 24.0 + (sod - (*glsItr).sod) / 3600.0;
			if (fabs(dt) < dtmin) {
				list_out.push_back(*glsItr);
			}
		}
	}
	return list_out;
}
