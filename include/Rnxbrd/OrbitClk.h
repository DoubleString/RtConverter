/*
 * Orbitclk.h
 *
 *  Created on: 2018/2/24
 *      Author: doublestring
 */

#ifndef INCLUDE_ORBITCLK_ORBITCLK_H_
#define INCLUDE_ORBITCLK_ORBITCLK_H_
#include <vector>
#include <string>
#include <list>
#include "../RtConverter/Const.h"
#include "../Rtklib/rtklib.h"
using namespace std;
namespace bamboo{
class GPSEPH{
public:
	GPSEPH(){
		mjd = 0;
		sod = 0.0;
	}
	char cprn[LEN_PRN];
	int mjd;
	double sod;
	// a[0]: SV clock offset
	// a[1]: SV clock drift
	// a[2]: SV clock drift rate
	double a0, a1, a2;
	// aode: age of ephemeris upload
	// crs, crc: Ortital radius correction
	// dn: Mean motion difference
	// m0: Mean anomaly at reference epoch
	// e: Eccentricity
	// cuc, cus: Latitude argument correction
	// roota: Square root of semi-major axis
	unsigned int aode;
	double crs, dn;
	double m0, cuc, e;
	double cus, roota;
	// toe, week: Ephemerides reference epoch in seconds with the week
	// cis, cic: Inclination correction
	// omega0: Longtitude of ascending node at the begining of the week
	// i0: Inclination at reference epoch
	// omega: Argument of perigee
	// omegadot: Rate of node's right ascension
	double toe, cic, omega0;
	double cis, i0, crc;
	double omega, omegadot;
	// idot: Rate of inclination angle
	// sesvd0:
	// resvd1:
	// accu: SV accuracy
	// hlth: SV health
	// tgd: Time group delay
	// aodc: Age of clock parameter upload
	double idot, resvd0, week, resvd1;
	double accu, hlth, tgd, tgd1, aodc;
};
class GLSEPH{
public:
	GLSEPH(){
		mjd = 0;
		sod = 0.0;
	}
	char cprn[LEN_PRN];
	int mjd,aode;
	double sod;
	// tau: SV clock bias
	// gama: SV relative frequency bias
	// tk: message frame time (tk+nd*86400) in seconds of the UTC week
	// pos: coordinate at ephemerides reference epoch in PZ-90
	// vel: velocity at ephemerides reference epoch in PZ-90
	// acc: acceleration at ephemerides reference epoch in PZ-90
	// health: SV health
	// frenum: frequency number
	// age: age of operation information
	double tau;
	double gamma;
	double tk;
	double pos[3];
	double vel[3];
	double acc[3];
	double health;
	double frenum;
	double age;
};
class OrbitAdapter{
public:
	OrbitAdapter(){orbType = ORB_NONE;}
	virtual ~OrbitAdapter(){}
	virtual int v_readOrbit(const char* cprn,int mjd,double sod,double* xsat){ cout << "***(ERROR):Base Class is not implement!" << endl;return 0;}
	virtual int v_readOrbit(const char* cprn,int mjd,double sod,double* xsat,int* iode){cout << "***(ERROR):Base Class is not implement!" << endl;return 0;}
	int m_getOrtType() { return orbType; }
protected:
	int orbType;
};
class ClkAdapter{
public:
	ClkAdapter(){memset(lRead,0,sizeof(int) * MAXSAT);clkType = CLK_NONE;}
	virtual ~ClkAdapter(){}
	virtual int v_readClk(const char* cprn,int mjd,double sod,double* sclk){ cout << "***(ERROR):Base Class is not implement!" << endl;return 0;}
	virtual int v_readClk(const char* cprn,int mjd,double sod,double* sclk,int* iode){ cout << "***(ERROR):Base Class is not implement!" << endl;return 0;}
	virtual int v_readClkDrift(const char* cprn,int mjd,double sod,double* vclk){cout << "***(ERROR):Base class is not implement" << endl; return 0;}

	virtual int v_lRead(int psat){return lRead[psat];}
	virtual void v_sRead(int psat,int value){lRead[psat] = value;}
	int m_getClkType(){return clkType;}
protected:
	int clkType;
	int lRead[MAXSAT];
};

// Be careful,it is not thread safe
class OrbitClkAdapter: public ClkAdapter,public OrbitAdapter {
public:
	OrbitClkAdapter(){isOutputOpen = false;outVer = 3.1;}
	virtual ~OrbitClkAdapter(){
		v_closeRnxEph();
		v_closeOutFile();
	};
	virtual bool v_isOpen(){cout << "***(ERROR):Base Class is not implement!" << endl; return false;}
	virtual void v_openRnxEph(string){ cout << "***(ERROR):Base Class is not implement!" << endl;};
	virtual void v_closeRnxEph(){};

	virtual void v_openOutFile(string,double ver){isOutputOpen = true;}
	virtual void v_outFileHeader(){}
	virtual void v_outFileEph(){}
	virtual void v_closeOutFile(){isOutputOpen = false;}
	virtual bool v_isOutFileOpen(){return isOutputOpen;}
protected:
	bool isOutputOpen;
	double outVer;
	char outFile[1024];
};
class BroadcastEphUtils{
public:
	virtual ~BroadcastEphUtils(){}
	virtual int m_brd2xyz(const char* mode,const char* cprn,int wk,double sow,double* xsat,double* clk,double* dtmin,double* tgd,int* iode);
	virtual int m_gls2xyz(const char* mode,const char* cprn,int wk,double sow,double* xsat,double* clk,double* dtmin,int* iode);
	void glsinit(double* x, GLSEPH& eph);
	void glsrkf4(double h, double* x, GLSEPH& eph);
	void glsfright(double* x, double* acc, GLSEPH& eph);
	void pz902wgs84(int mjd, double sod, double* pos, double *xsat,const char* trans);
protected:
	int neph[MAXSYS];
	list<GPSEPH> gnssEph[MAXSYS];
	list<GLSEPH> glsEph;
};
class RnxEphFileAdapter : public OrbitClkAdapter,public BroadcastEphUtils{
public:
	RnxEphFileAdapter();
	virtual ~RnxEphFileAdapter();
	virtual void v_openRnxEph(string);
	virtual inline bool v_isOpen(){return this->isOpen;}
	virtual void v_closeRnxEph();
	virtual int v_readOrbit(const char* cprn,int mjd,double sod,double* xsat);
	virtual int v_readClk(const char* cprn,int mjd,double sod,double* sclk);
    // for real-time oc
	virtual int v_readOrbit(const char* cprn,int mjd,double sod,double* xsat,int* iode);
	virtual int v_readClk(const char* cprn,int mjd,double sod,double* sclk,int* iode);

	// for clock drift
	virtual int v_readClkDrift(const char* cprn,int mjd,double sod,double* vclk){return 0;}
	// ionospheric correction parameters
	char ionc[MAXSYS][2][LEN_EPHION];
	double ion[MAXSYS][2][4];

	// for search 
	list<GPSEPH> m_getCurrentEph_G(string cprn,int mjd,double sod);
	list<GLSEPH> m_getCurrentEph_R(string cprn,int mjd,double sod);
protected:
	void m_readRnxnav(char mode);
	char curFile[1024];
	int mjd0,mjd1;
	double sod0,sod1;
	int isOpen;
private:
	// header part
	double ver;
	// corrections to transform the system to UTC or other time systems
	char timc[MAXSYS][2][LEN_EPHION];
	double tim[MAXSYS][2][4];
	// number of leap second since 6-Jan-980
	int leap;
};
}
#endif /* INCLUDE_ORBITCLK_ORBITCLK_H_ */
