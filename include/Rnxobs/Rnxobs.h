/*
 * Rnxobs.h
 *
 *  Created on: 2018/2/3
 *      Author: doublestring
 */

#ifndef INCLUDE_BAMBOO_RNXOBS_H_
#define INCLUDE_BAMBOO_RNXOBS_H_
#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <string>
#include "../RtConverter/Const.h"
#include "../Rtklib/rtklib.h"
using namespace std;
namespace bamboo{
class Rnxobs{
public:
	Rnxobs(){
		// Reset Varibles
		iRead = 0;
		mjd = 0;
		sod = 0;
		isOutputOpen = false;
		outVer = 2.11;
		memset(flag,0,sizeof(int) * MAXSAT * MAXFREQ);
		memset(obs,0,sizeof(double) * MAXSAT * 2 * MAXFREQ);
		memset(fob,0,sizeof(char) * MAXSAT * MAXFREQ * 2 * LEN_OBSTYPE);
		memset(outFile,0,sizeof(outFile));
		memset(snr,0,sizeof(snr));
	}
	Rnxobs(string name){
		// Reset Varibles
		iRead = 0;
		staname = name;
		mjd = 0;
		sod = 0;
		isOutputOpen = false;
		outVer = 2.11;
		memset(flag,0,sizeof(int) * MAXSAT * MAXFREQ);
		memset(obs,0,sizeof(double) * MAXSAT * 2 * MAXFREQ);
		memset(fob,0,sizeof(char) * MAXSAT * MAXFREQ * 2 * LEN_OBSTYPE);
		memset(outFile,0,sizeof(outFile));
		memset(snr,0,sizeof(snr));
	}
	virtual ~Rnxobs(){
		v_closeRnx();
	};
	virtual void v_openRnx(string){}
	virtual void v_readEpoch(int& mjd,double& sod){}
	virtual void v_closeRnx(){}
	virtual int v_isOpen(){return 0;}

	void m_checkObs(int isat,double* obs,double* dop,double * snr);
	int getRead(){return iRead;}
	void setRead(int iset){iRead = iset;}
	/// Saving Part
	virtual void v_openOutFile(string,double ver);
	virtual void v_outFileHeader();
	virtual void v_outFileObs(int otmjd,double otsod);
	virtual void v_closeOutFile();
	virtual bool v_isOutFileOpen(){return isOutputOpen;}

	double obs[MAXSAT][2*MAXFREQ],dop[MAXSAT][2*MAXFREQ],sod;
	char fob[MAXSAT][2*MAXFREQ][LEN_OBSTYPE];
	double snr[MAXSAT][MAXFREQ];
	string staname,rectype, anttype;
	int mjd,flag[MAXSAT][MAXFREQ];
protected:
	int iRead;
	bool isOutputOpen;
	double outVer;
	char outFile[1024];
};
class RnxobsFile_sat:public Rnxobs{
public:
	RnxobsFile_sat(string);//
	~RnxobsFile_sat(){
		v_closeRnx();
	};
	virtual void v_openRnx(string); //
	virtual void v_readEpoch(int&,double&);
	virtual void v_closeRnx();
	inline int v_isOpen(){return in.is_open();}
protected:
	//RnxFile Head
	void m_readRnxHead();
	double ver;
	char sys[4];
	string  recnum, antnum,mark;
	double x, y, z, e, n, h, intv;
	int nprn, fact1, fact2, nobstype[MAXSYS];
	string obstype[MAXSYS][MAXOBSTYP], cprn_ob[MAXSAT];
	int t0[6], t1[6], nsys;


	char usetype[MAXSAT][MAXFREQ];
	int tstore[MAXSAT][MAXFREQ];
private:
	ifstream in;
	string dir;
};
}
#endif /* INCLUDE_BAMBOO_RNXOBS_H_ */
