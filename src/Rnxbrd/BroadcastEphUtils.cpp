/*
 * BroadcastEphUtils.cpp
 *
 *  Created on: 2018/3/2
 *      Author: doublestring
 */

#include <iostream>
#include <cmath>
#include "../../include/Rnxbrd/OrbitClk.h"
#include "../../include/Rnxobs/Com.h"
using namespace std;
using namespace bamboo;
int BroadcastEphUtils::m_brd2xyz(const char* mode,const char* cprn,int wk,double sow,double* xsat,double* clk,double* dtmin,double* tgd,int* iode){
	int i,isys,lfind;
	double GMEARTH, EARTH_ROTATE;
	double dt, dtc, pos[3], svpos[3], geovel[3], detavel[3];
	double a, xn, xm, ex, e, v0, vs, vc, phi, ccc, sss, du, dr, di, r, u, xi,
			xx, yy,intv;
	double xnode, xpdot, ypdot, asc, xinc, xp, yp, asctrm, v;
	double edot, cosf2, cose2, nudot, rdot, udot, idot;
	list<GPSEPH>::iterator gpsItr;
	GPSEPH* eph;
	intv = 2.0;
	switch (cprn[0]) {
	case 'G':
		GMEARTH = 3.986005E14;
		EARTH_ROTATE = 7.2921151467E-5;
		break;
	case 'E':
		GMEARTH = 3.986004418E14;
		EARTH_ROTATE = 7.2921151467E-5;
		break;
	case 'C':
		GMEARTH = 3.986004418E14;
		EARTH_ROTATE = 7.2921150E-5;
		break;
	case 'J':
		GMEARTH = 3.986005E14;
		EARTH_ROTATE = 7.2921150E-5;
		break;
	default:
		GMEARTH = GME;
		EARTH_ROTATE = E_ROTATE;
		break;
	}
	lfind = 0;
	*clk = 0;
	*dtmin = 8.0;
	isys = index_string(SYS,cprn[0]);
	double tmax = 0.0,ttag;
	for(gpsItr = this->gnssEph[isys].begin();gpsItr != this->gnssEph[isys].end();++gpsItr){
		if(strstr((*gpsItr).cprn,cprn)){
			dt = (wk - (*gpsItr).week) * 168.0
							+ (sow - (*gpsItr).toe) / 3600.0;
			if (*iode >= 0) {
				if ((*gpsItr).aode == *iode) {
					*dtmin = fabs(dt);
					eph = &(*gpsItr);
					lfind = 1;
					break;
				}
			} else if(*iode == -1){
				if (fabs(dt) <= *dtmin) {
					*dtmin = fabs(dt);
					eph = &(*gpsItr);
					lfind = 1;
				}
			} else if(*iode == -2){
				ttag = (*gpsItr).mjd * 86400.0 + (*gpsItr).sod;
				// For the real-time process : always the newest
				if (ttag > tmax && fabs(dt) <= intv) {
					tmax = ttag;
					*dtmin = fabs(dt);
					lfind = 1;
					eph = &(*gpsItr);
				}
			}
		}
	}
	if(lfind == 0)
		return 0;
	*iode = eph->aode;
	dt = (wk - eph->week) * 604800.0 + sow - eph->toe;
	dtc = (wk * 7 + 44244 - eph->mjd) * 86400.0 + sow - eph->sod;

	*clk = eph->a0 + (eph->a1 + eph->a2 * dtc) * dtc;
	*tgd = eph->tgd;
//	if(cprn[0] == 'C'){
//		// change the broadcast sclock into L2/L5
//		r_bds=(BDS_B1*BDS_B1)/(BDS_B2*BDS_B2);
//		*clk=*clk-r_bds/(r_bds-1)*eph->tgd+1/(r_bds-1)*eph->tgd1;
//	}
	if (strstr(mode, "ynn") != NULL)
		return 1;
	a = eph->roota * eph->roota;
	xn = sqrt(GMEARTH / a / a / a);
	xn = xn + eph->dn;
	xm = eph->m0 + xn * dt;
	ex = xm;
	e = eph->e;
	for (i = 0; i < 12; i++)
		ex = xm + e * sin(ex);
	v0 = 1.0 - e * cos(ex);
	vs = sqrt(1.0 - e * e) * sin(ex) / v0;
	vc = (cos(ex) - e) / v0;
	v = fabs(asin(vs));
	// decide vs
	if (vc >= 0) {
		if (vs < 0)
			v = 2.0 * PI - v;
	} else {
		if (vs <= 0)
			v = PI + v;
		else
			v = PI - v;
	}

	phi = v + eph->omega;
	ccc = cos(2.0 * phi);
	sss = sin(2.0 * phi);
	du = eph->cuc * ccc + eph->cus * sss;
	dr = eph->crc * ccc + eph->crs * sss;
	di = eph->cic * ccc + eph->cis * sss;
	r = a * (1.0 - e * cos(ex)) + dr;
	u = phi + du;

	xi = eph->i0 + eph->idot * dt + di;
	xx = r * cos(u);
	yy = r * sin(u);

	if (strstr(cprn, "C01") != NULL || strstr(cprn, "C02") != NULL
			|| strstr(cprn, "C03") != NULL || strstr(cprn, "C04") != NULL
			|| strstr(cprn, "C05") != NULL) {
		xnode = eph->omega0  + eph->omegadot * dt;
	} else {
		xnode = eph->omega0 + (eph->omegadot - EARTH_ROTATE) * dt;
	}

	//IF THE TIME SYSTEM IS DIFFERENT,IT IS BETTER TO COMPUTE THE ORBIT IN THEIR OWN SYSTEM,
	//BECAUSE TOE WILL BE DIFFERENT IF YOU TRANSFORM IT TO OTHER SYSTEM

	xnode = xnode - EARTH_ROTATE * eph->toe;
	xsat[0] = xx * cos(xnode) - yy * cos(xi) * sin(xnode);
	xsat[1] = xx * sin(xnode) + yy * cos(xi) * cos(xnode);
	xsat[2] = yy * sin(xi);

	for (i = 0; i < 3; i++)
		svpos[i] = xsat[i];

	if (strstr(cprn, "C01") != NULL || strstr(cprn, "C02") != NULL
			|| strstr(cprn, "C03") != NULL || strstr(cprn, "C04") != NULL
			|| strstr(cprn, "C05") != NULL) {
		// rx(-5.0*rad) for GEO
		pos[0] = xsat[0];
		pos[1] = cos(-5.0 * DEG2RAD) * xsat[1] + sin(-5.0 * DEG2RAD) * xsat[2];
		pos[2] = sin(5.0 * DEG2RAD) * xsat[1] + cos(-5.0 * DEG2RAD) * xsat[2];

		// rz(wearth*dt)
		xsat[0] = pos[0] * cos(EARTH_ROTATE * dt)
				+ pos[1] * sin(EARTH_ROTATE * dt);
		xsat[1] = -1.0 * pos[0] * sin(EARTH_ROTATE * dt)
				+ pos[1] * cos(EARTH_ROTATE * dt);
		xsat[2] = pos[2];
	}
	if (strstr(mode, "yyn") != NULL)
		return 1;

	asc = xnode;
	xinc = xi;
	xp = xx;
	yp = yy;

	edot = xn / (1 - e * cos(ex));
	cosf2 = cos(v / 2) * cos(v / 2);
	cose2 = cos(ex / 2) * cos(ex / 2);
	nudot = sqrt((1 + e) / (1 - e)) * cosf2 * edot / cose2;

	rdot =a * e * edot * sin(ex)
	+ 2 * (eph->crs * ccc - eph->crc * sss) * nudot;

	udot = (1 + 2.0 * eph->cus * ccc - 2.0 * eph->cuc * sss)
			* nudot;

	idot = eph->idot
			+ 2.0 * (eph->cic * ccc - eph->cis * sss) * nudot;

	if (strstr(cprn, "C01") != NULL || strstr(cprn, "C02") != NULL
			|| strstr(cprn, "C03") != NULL || strstr(cprn, "C04") != NULL
			|| strstr(cprn, "C05") != NULL) {
		asctrm = eph->omegadot;
	} else {
		asctrm = eph->omegadot - EARTH_ROTATE;
	}

	xpdot = rdot * cos(u) - r * sin(u) * udot;
	ypdot = rdot * sin(u) + r * cos(u) * udot;

	xsat[3] = xpdot * cos(asc) - ypdot * cos(xinc) * sin(asc)
			+ yp * sin(asc) * sin(xinc) * idot - xp * sin(asc) * asctrm
			- yp * cos(xinc) * cos(asc) * asctrm;
	xsat[4] = xpdot * sin(asc) + ypdot * cos(xinc) * cos(asc)
			- yp * cos(asc) * sin(xinc) * idot + xp * cos(asc) * asctrm
			- yp * cos(xinc) * sin(asc) * asctrm;
	xsat[5] = ypdot * sin(xinc) + yp * cos(xinc) * idot;

	if (strstr(cprn, "C01") != NULL || strstr(cprn, "C02") != NULL
			|| strstr(cprn, "C03") != NULL || strstr(cprn, "C04") != NULL
			|| strstr(cprn, "C05") != NULL) {
		// rx(-5.0*rad) for GEO
		geovel[0] = xsat[3];
		geovel[1] = cos(-5.0 * DEG2RAD) * xsat[4]
				+ sin(-5.0 * DEG2RAD) * xsat[5];
		geovel[2] = sin(5.0 * DEG2RAD) * xsat[4]
				+ cos(-5.0 * DEG2RAD) * xsat[5];

		// rz(wearth*dt)
		detavel[0] = geovel[0] * cos(EARTH_ROTATE * dt)
				+ geovel[1] * sin(EARTH_ROTATE * dt);
		detavel[1] = -1.0 * geovel[0] * sin(EARTH_ROTATE * dt)
				+ geovel[1] * cos(EARTH_ROTATE * dt);
		detavel[2] = geovel[2];

		// dot(R(w*t))*R(i)*XPOS
		geovel[0] = svpos[0];
		geovel[1] = cos(-5.0 * DEG2RAD) * svpos[1]
				+ sin(-5.0 * DEG2RAD) * svpos[2];
		geovel[2] = sin(5.0 * DEG2RAD) * svpos[1]
				+ cos(-5.0 * DEG2RAD) * svpos[2];

		// rz(wearth*dt)
		xsat[3] = EARTH_ROTATE
				* (-1.0 * geovel[0] * sin(EARTH_ROTATE * dt)
						+ geovel[1] * cos(EARTH_ROTATE * dt));
		xsat[4] = EARTH_ROTATE
				* (-1.0 * geovel[0] * cos(EARTH_ROTATE * dt)
						- geovel[1] * sin(EARTH_ROTATE * dt));
		xsat[5] = 0.0;

		// summary
		xsat[3] += detavel[0];
		xsat[4] += detavel[1];
		xsat[5] += detavel[2];
	}
	return 1;
}
int BroadcastEphUtils::m_gls2xyz(const char* mode,const char* cprn,int wk,double sow,double* xsat,double* clk,double* dtmin,int* iode){
	double dt = 0.0, dt1 = 0.0, dt2 = 0.0, x[6],sod,intv;
	int i, nsign,lfind,mjd;
	list<GLSEPH>::iterator glsItr;
	GLSEPH* ephr;
	lfind = 0;
	*clk = 0.0;
	*dtmin = 12.0;
	wksow2mjd(wk,sow,&mjd,&sod);
	memset(x, 0, sizeof(double) * 6);
	double tmax = 0.0,ttag;
	intv = 2.0;
	for(glsItr = this->glsEph.begin();glsItr != this->glsEph.end();++glsItr){
		if(strstr((*glsItr).cprn,cprn)){
			dt = (mjd - (*glsItr).mjd) * 24.0
					+ (sod - (*glsItr).sod) / 3600.0;
			if (*iode >= 0) {
				if (*iode == (*glsItr).aode) {
					*dtmin = fabs(dt);
					ephr = &(*glsItr);
					lfind = 1;
					break;
				}
			} else if(*iode == -1){
				if (fabs(dt) <= *dtmin) {
					*dtmin = fabs(dt);
					ephr = &(*glsItr);
					lfind = 1;
				}
			} else if(*iode == -2){
				ttag = (*glsItr).mjd * 86400.0 + (*glsItr).sod;
				if (ttag > tmax && fabs(dt) <= intv) {
					tmax = ttag;
					*dtmin = fabs(dt);
					ephr = &(*glsItr);
					lfind = 1;
				}
			}
		}
	}
	if(!lfind)
		return 0;
	*iode = ephr->aode;
	dt = (mjd - ephr->mjd) * 86400.0 + sod - ephr->sod;
	*clk = ephr->tau + ephr->gamma * dt;
	if (strstr(mode, "ynn") != NULL)
		return 1;

	dt1 = dt - (int) (dt / 60.0) * 60.0;
	dt2 = dt1;
	nsign = 1;
	if (dt < 0.0)
		nsign = -1;

	glsinit(x, *ephr);

	if (dt1 != 0.0)
		glsrkf4(dt1, x, *ephr);

	while (ABS(dt-dt1) >= 60.0) {
		if (dt - dt1 != 0.0) {
			glsrkf4(nsign * 60.0, x,  *ephr);
			dt1 = dt1 + nsign * 60.0;
		}
	}
	pz902wgs84(mjd, sod, x, xsat, "MCC");
	for (i = 0; i < 3; i++)
		xsat[i] = xsat[i] * 1e3;
	for (i = 3; i < 6; i++)
		xsat[i] = x[i] * 1e3;
	return 1;
}

void BroadcastEphUtils::glsinit(double* x,GLSEPH& eph){
	int i;
	for (i = 0; i < 3; i++)
		x[i] = eph.pos[i];
	for (i = 3; i < 6; i++)
		x[i] = eph.vel[i - 3];
}
void BroadcastEphUtils::glsrkf4(double h, double* x, GLSEPH& eph){
	int i, j;
	double coef[4], xj[6], acc[6], funct[6];
	coef[0] = h / 2;
	coef[1] = h / 2;
	coef[2] = h;
	coef[3] = h;
	for (i = 0; i < 6; i++)
		xj[i] = x[i];
	glsfright(xj, acc, eph);
	for (i = 0; i < 6; i++)
		funct[i] = xj[i];
	for (j = 0; j < 3; j++) {
		for (i = 0; i < 6; i++) {
			xj[i] = x[i] + coef[j] * acc[i];
			funct[i] = funct[i] + coef[j + 1] * acc[i] / 3.0;
		}
		glsfright(xj, acc, eph);
	}
	for (i = 0; i < 6; i++)
		x[i] = funct[i] + h * acc[i] / 6.0;
}
void BroadcastEphUtils::glsfright(double* x, double* acc,GLSEPH& eph){
	const double c20 = -1.08262575e-3;
	const double GMEARTH = 3.986004418e5;
	const double EROT = 7.292115e-5;
	const double ER = 6.378136e3;
	int i;
	double factor1, factor2, r, vec1[3], vec2[3];
	r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
	factor1 = -GMEARTH / pow(r, 3);
	factor2 = 1.5 * GMEARTH * c20 * ER * ER / pow(r, 5);
	vec1[0] = EROT * EROT * x[0] + 2.0 * EROT * x[4];
	vec1[1] = EROT * EROT * x[1] - 2.0 * EROT * x[3];
	vec1[2] = 0.0;
	vec2[0] = x[0] * (1.0 - 5.0 * pow((x[2] / r), 2));
	vec2[1] = x[1] * (1.0 - 5.0 * pow((x[2] / r), 2));
	vec2[2] = x[2] * (3.0 - 5.0 * pow((x[2] / r), 2));
	for (i = 0; i < 3; i++) {
		acc[i] = x[i + 3];
		acc[i + 3] = eph.acc[i] + vec1[i] + factor2 * vec2[i] + factor1 * x[i];
	}
}
void BroadcastEphUtils::pz902wgs84(int mjd, double sod, double* pos, double *xsat, const char* trans) {
	double factor, biase[3], mat[3][3], vec[3];
	int i, j;
	if (timdif(mjd, sod, 54362, 86367.0) >= 0.0) {
		factor = 1.0;
//		biase[0] = -0.36e-3;
//		biase[1] = 0.08e-3;
//		biase[2] = 0.18e-3;

		biase[0] = 0.008e-3;
		biase[1] = 0.001e-3;
		biase[2] = 0.001e-3;
		for (i = 0; i < 3; i++)
			for (j = 0; j < 3; j++) {
				if (i == j) {
					mat[i][j] = 1.0;
				} else {
					mat[i][j] = 0.0;
				}
			}
	} else {
		if (strstr(trans, "MCC") != NULL) {
			factor = 1.0 + 22.0e-9;
			biase[0] = -0.47e-3;
			biase[1] = -0.51e-3;
			biase[2] = -1.56e-3;
			mat[0][0] = 1.0;
			mat[0][1] = -1.728e-6;
			mat[0][2] = -1.7e-8;
			mat[1][0] = 1.728e-6;
			mat[1][1] = 1.0;
			mat[1][2] = 7.6e-8;
			mat[2][0] = 1.7e-8;
			mat[2][1] = -7.6e-8;
			mat[2][2] = 1.0;
		} else if (strstr(trans, "RUS") != NULL) {
			factor = 1.0 - 1.2e-7;
			biase[0] = -1.1e-3;
			biase[1] = -0.3e-3;
			biase[2] = -0.9e-3;
			mat[0][0] = 1.0;
			mat[0][1] = -8.2e-7;
			mat[0][2] = 0.0;
			mat[1][0] = 8.2e-7;
			mat[1][1] = 1.0;
			mat[1][2] = 0.0;
			mat[2][0] = 0.0;
			mat[2][1] = 0.0;
			mat[2][2] = 1.0;
		} else if (strstr(trans, "LMU") != NULL) {
			factor = 1.0;
			biase[0] = 0.0;
			biase[1] = 0.0;
			biase[2] = 0.0;
			mat[0][0] = 1.0;
			mat[0][1] = -1.6e-6;
			mat[0][2] = 0.0;
			mat[1][0] = 1.6e-6;
			mat[1][1] = 1.0;
			mat[1][2] = 0.0;
			mat[2][0] = 0.0;
			mat[2][1] = 0.0;
			mat[2][2] = 1.0;
		} else if (strstr(trans, "MIT") != NULL) {
			factor = 1.0;
			biase[0] = 0.0;
			biase[1] = 2.5e-3;
			biase[2] = 0.0;
			mat[0][0] = 1.0;
			mat[0][1] = -1.9e-6;
			mat[0][2] = 0.0;
			mat[1][0] = 1.9e-6;
			mat[1][1] = 1.0;
			mat[1][2] = 0.0;
			mat[2][0] = 0.0;
			mat[2][1] = 0.0;
			mat[2][2] = 1.0;
		} else {
			printf("***ERROR(pz902wgs84): unknown transformation type %s!\n",
					trans);
			exit(1);
		}
	}
	matmpy((double*)mat, pos, vec, 3, 3, 1);
	for (i = 0; i < 3; i++)
		xsat[i] = biase[i] + factor * vec[i];
}
