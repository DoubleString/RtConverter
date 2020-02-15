#include "../../include/Rnxobs/RtRinexStream.h"
#include "../../include/Rnxobs/Com.h"
#include "../../include/Rtklib/rtklib_fun.h"
#include "../../include/RtConverter/RtConvItem.h"
#include "../../include/RtConverter/RtConverter.h"
#include "../../include/Rnxobs/Rnxobs.h"
#include "../../include/Rnxbrd/OrbitClk.h"
#include <list>

using namespace bamboo;
Deploy RtRinexStream::configs_sav;
RtRinexStream::~RtRinexStream() {
	lcont = false;
}
void RtRinexStream::openStream() {
	configs_sav = Deploy::s_getConfigures();
	def_thread_t pid_handle;
	lcont = true;
	time(&lastCheck);
#ifdef _WIN32
	CreateThread(NULL, 0, s_pthRniexStream_win, this, 0, NULL);
#else
	if (0 != pthread_create(&pid_handle, NULL, &s_pthRinexStream, this)) {
		cout
			<< "***ERROR(s_VRSMain):cant create thread to handle request!"
			<< endl;
		exit(1);
	}
#endif
}
void RtRinexStream::openStream(list<RtConverter*>& svrs_in) {
	configs_sav = Deploy::s_getConfigures();
	def_thread_t pid_handle;
	lcont = true;
	time(&lastCheck);
#ifdef _WIN32
	CreateThread(NULL, 0, s_pthRniexStream_win, this, 0, NULL);
#else
	if (0 != pthread_create(&pid_handle, NULL, &s_pthRinexStream, this)) {
		cout
			<< "***ERROR(s_VRSMain):cant create thread to handle request!"
			<< endl;
		exit(1);
	}
#endif
	this->m_svrs = svrs_in;
}
void RtRinexStream::closeStream() {
	lcont = false;
}
#ifdef _WIN32
DWORD WINAPI RtRinexStream::s_pthRniexStream_win(LPVOID args) {
	RtRinexStream* str = (RtRinexStream*)args;
	str->m_routinue();
	return NULL;
}
#endif
void RtRinexStream::m_routinue() {
	time_t tt;
	/// open first
	int curmjd,mjd,isat;
	double sod;
	list<string>::iterator sitr;
	list<RtConverter*>::iterator convItr;
	vector<RnxobsFile_sat*>::iterator rnxItr;
	RnxobsFile_sat* rnx;
	list<GPSEPH> eph_g;
	list<GLSEPH> eph_r;
	for (sitr = configs_sav.post_stalist.begin(); sitr != configs_sav.post_stalist.end(); sitr++) {
		rnx = new RnxobsFile_sat(*sitr);
		m_rnxs.push_back(rnx);
	}
	mOacs = new RnxEphFileAdapter();
	mjd = configs_sav.mjd0;
	sod = configs_sav.sod0;
	while (lcont) {
		time(&tt);
		if (tt - lastCheck > 5) {
			lastCheck = tt;
			m_adaptConfigures();
		}
		/// open observation here
		for (rnxItr = m_rnxs.begin(); rnxItr != m_rnxs.end(); rnxItr++) {
			if (!(*rnxItr)->v_isOpen()) {
				(*rnxItr)->v_openRnx(toString(mjd) + ":" + toString(sod));
			}
			(*rnxItr)->v_readEpoch(mjd, sod);
		}
		/// open brdm here 
		if(!mOacs->v_isOpen())
			mOacs->v_openRnxEph(toString(mjd) + ":" + toString(sod));
		/// will update to the main thread
		RtConvItem curitem = m_makeupItems(mjd,sod,m_rnxs);
		for (convItr = m_svrs.begin(); convItr != m_svrs.end(); convItr++) {
			if(((*convItr)->getType() == RtConverter::ConvType::obsbnc || (*convItr)->getType() == RtConverter::ConvType::obsbin) && curitem.obslist.size() > 0)
				(*convItr)->inputObs(curitem);
		}
		/// got brdm here
		for (isat = 0; isat < configs_sav.nprn; isat++) {
			if (configs_sav.cprn[isat][0] == 'R') {
				eph_r.clear();
				if (m_makeupEph_R(eph_r,configs_sav.cprn[isat], mjd, sod)) {
					for (convItr = m_svrs.begin(); convItr != m_svrs.end(); convItr++) {
						if ((*convItr)->getType() == RtConverter::ConvType::eph) 
							(*convItr)->inputBrdm_R(configs_sav.cprn[isat],eph_r);
					}
				}
			}
			else {
				eph_g.clear();
				if (m_makeupEph_G(eph_g,configs_sav.cprn[isat], mjd, sod)) {
					for (convItr = m_svrs.begin(); convItr != m_svrs.end(); convItr++) {
						if ((*convItr)->getType() == RtConverter::ConvType::eph)
							(*convItr)->inputBrdm_G(configs_sav.cprn[isat],eph_g);
					}
				}
			}
		}
		curmjd = mjd;
		timinc(mjd, sod, configs_sav.dintv, &mjd, &sod);
		if (curmjd != mjd) {
			for (rnxItr = m_rnxs.begin(); rnxItr != m_rnxs.end();rnxItr++) {
				(*rnxItr)->v_closeRnx();
				(*rnxItr)->v_closeOutFile();
			}
			mOacs->v_closeRnxEph();
		}
		if ((mjd - configs_sav.mjd1) * 86400.0 + sod - configs_sav.sod1 >= 0.0)
			break;
		sleepms(NINT(configs_sav.dintv * 1000/configs_sav.speed));
	}
	for (rnxItr = m_rnxs.begin(); rnxItr != m_rnxs.end(); rnxItr++) {
		delete (*rnxItr);
	}
	m_rnxs.clear();
}
RtConvItem RtRinexStream::m_makeupItems(int mjd,double sod,vector<RnxobsFile_sat*>& rnxs) {
	int isat,isys,ifreq,week;
	double sow;
	RtConvItem item_ep;
	vector<RnxobsFile_sat*>::iterator rnxItr;
	for (rnxItr = rnxs.begin(); rnxItr != rnxs.end(); rnxItr++) {
		map<string, RtSatObs> mapObs;
		for (isat = 0; isat < configs_sav.nprn; isat++) {
			isys = index_string(SYS, configs_sav.cprn[isat][0]);
			RtSatObs sat;
			for (ifreq = 0; ifreq < configs_sav.nfreq[isys]; ifreq++) {
				if ((*rnxItr)->obs[isat][ifreq] != 0.0) {
					sat.obs[ifreq] = (*rnxItr)->obs[isat][ifreq];
					sat.obs[MAXFREQ + ifreq] = (*rnxItr)->obs[isat][MAXFREQ + ifreq];

					sat.fob[ifreq] = (*rnxItr)->fob[isat][ifreq];
					sat.fob[MAXFREQ + ifreq] = (*rnxItr)->fob[isat][MAXFREQ + ifreq];

					sat.snr[ifreq] = (*rnxItr)->snr[isat][ifreq];
					mapObs[configs_sav.cprn[isat]] = sat;
				}
			}
		}
		if(mapObs.size() > 0)
			item_ep.obslist[(*rnxItr)->staname] = mapObs;
	}
	mjd2wksow(mjd,sod,&week,&sow);
	item_ep.curt = gpst2time(week, sow);
	return item_ep;
}
void RtRinexStream::m_adaptConfigures() {
	Deploy config_new = Deploy::s_getConfigures();
	/// diff the new configures and saved configures and update the memory
	list<string> add_, del_;
	vector<RnxobsFile_sat*>::iterator rnxItr;
	list<string>::iterator strItr, strItr_mem;
	for (strItr = config_new.post_stalist.begin(); strItr != config_new.post_stalist.end(); strItr++) {
		bool lfind = false;
		for (strItr_mem = configs_sav.post_stalist.begin(); strItr_mem != configs_sav.post_stalist.end(); ++strItr_mem) {
			if ((*strItr) == (*strItr_mem)) {
				lfind = true;
				break;
			}
		}
		if (!lfind) add_.push_back(*strItr);
	}
	for (strItr_mem = configs_sav.post_stalist.begin(); strItr_mem != configs_sav.post_stalist.end(); ++strItr_mem) {
		bool lfind = false;
		for (strItr = config_new.post_stalist.begin(); strItr != config_new.post_stalist.end(); strItr++) {
			if ((*strItr) == (*strItr_mem)) {
				lfind = true;
				break;
			}
		}
		if (!lfind) del_.push_back(*strItr_mem);
	}
	/// adapt the current process
	for (strItr = add_.begin(); strItr != add_.end(); strItr++) {
		RnxobsFile_sat* rnx = new RnxobsFile_sat(*strItr);
		this->m_rnxs.push_back(rnx);
	}
	for (strItr = del_.begin(); strItr != del_.end(); strItr++) {
		for (rnxItr = m_rnxs.begin(); rnxItr != m_rnxs.end(); rnxItr++) {
			if ((*rnxItr)->staname == (*strItr)) {
				(*rnxItr)->v_closeRnx();
				delete (*rnxItr);
				m_rnxs.erase(rnxItr);
				break;
			}
		}
	}
	configs_sav = config_new;
}
void* RtRinexStream::s_pthRinexStream(void* args) {
	RtRinexStream* str = (RtRinexStream*)args;
	str->m_routinue();
	return NULL;
}
bool RtRinexStream::m_makeupEph_G(list<GPSEPH>& list_in,string cprn,int mjd,double sod) {
	list_in = mOacs->m_getCurrentEph_G(cprn,mjd,sod);
	return true;
}
bool RtRinexStream::m_makeupEph_R(list<GLSEPH>& list_in,string cprn,int mjd, double sod) {
	list_in = mOacs->m_getCurrentEph_R(cprn,mjd, sod);
	return true;
}


