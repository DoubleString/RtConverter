#include "../../include/RtConverter/RtConverter.h"
#include "../../include/RtConverter/Const.h"
#include "../../include/Rnxobs/RtRinexStream.h"
#include "../../include/VrsObs/VrsNtripClient.h"
#include "../../include/Rnxobs/Com.h"
#include "../../include/Rtklib/rtklib_fun.h"
#include <algorithm>
#include <signal.h>
using namespace bamboo;
using namespace std;
RtConverter::RtConverter() {
	initlock(&_mutex);
	obsdelay = 5;
	port = -1;
	type = ConvType::obsbnc;
}
RtConverter::RtConverter(int delay,unsigned int port_in, ConvType type_in) {
	initlock(&_mutex);
	obsdelay = delay;
	port = port_in;
	type = type_in;
}
RtConverter::~RtConverter() {
}
/// which will be called in the other threads
int RtConverter::inputObs(RtConvItem& item) {
	def_lock(&_mutex);
	if (obsqueue.empty()) {
		obsqueue.push_back(item);
		def_unlock(&_mutex);
		return 1;
	}
	RtConvItem& front = this->obsqueue.front();
	RtConvItem& end = this->obsqueue.back();
	if (item.curt.time < front.curt.time) {
		def_unlock(&_mutex);
		return 0;
	}
	if (item.curt.time > end.curt.time) {
		this->obsqueue.push_back(item);
		def_unlock(&_mutex);
		return 1;
	}
	list<RtConvItem>::iterator itrlist;
	for (itrlist = obsqueue.begin(); itrlist != obsqueue.end(); ++itrlist) {
		if (item.curt.time <= (*itrlist).curt.time)
			break;
	}
	if ((*itrlist).curt.time != item.curt.time) {
		obsqueue.insert(itrlist, item);
		def_unlock(&_mutex);
		return 1;
	}
	/// queue already exist the same time,will add into the memory
	std::map<std::string, std::map<std::string, RtSatObs> >::iterator itr, obitr_in = item.obslist.begin();
	std::map<std::string, std::map<std::string, RtSatObs> >::iterator obitr_mem = (*itrlist).obslist.begin();
	while (obitr_in != item.obslist.end()) {
		itr = (*itrlist).obslist.find((*obitr_in).first);
		if (itr != (*itrlist).obslist.end()) {
			++obitr_in;
			continue;
		}
		(*itrlist).obslist[(*obitr_in).first] = (*obitr_in).second;
		++obitr_in;
	}
	def_unlock(&_mutex);
	return 1;
}
void RtConverter::beginProcess() {
	def_thread_t pid_handle;
#ifdef _WIN32
	CreateThread(NULL, 0, s_pthProcess_win, this, 0, NULL);
#else
	if (0 != pthread_create(&pid_handle, NULL, &s_pthProcess, this)) {
		cout
			<< "***ERROR(s_VRSMain):cant create thread to handle request!"
			<< endl;
		exit(1);
	}
#endif
}
#ifdef _WIN32
DWORD WINAPI RtConverter::s_pthProcess_win(LPVOID args) {
	RtConverter* svr = (RtConverter*)args;
	if(svr->type == ConvType::obsbnc || svr->type == ConvType::obsbin)
		svr->routineObs();
	else
		svr->routineEph();
	return NULL;
}
#endif
void* RtConverter::s_pthProcess(void* args) {
	RtConverter* svr = (RtConverter*)args;
	if(svr->type == ConvType::obsbnc || svr->type == ConvType::obsbin)
		svr->routineObs();
	else
		svr->routineEph();
	return NULL;
}
void RtConverter::routineObs() {
	time_t now;
	int nbyte;
	char cmd[256];
	unsigned char buff[1024];
	list<RtConvItem> item_sd;
	sprintf(cmd, ":%d", port);
	strinit(&svr);
	stropen(&svr, STR_TCPSVR, STR_MODE_RW, cmd);
	while (true) {
		/* test for re-read configures*/
		item_sd.clear();
		def_lock(&_mutex);
		time(&now);
		list<RtConvItem>::reverse_iterator itrlist = obsqueue.rbegin();
		while (itrlist != obsqueue.rend()) {
			if (now - (*itrlist).gent >= obsdelay) {
				break;
			}
			++itrlist;
		}
		if (itrlist != obsqueue.rend()) {
			gtime_t sd_t = (*itrlist).curt;
			list<RtConvItem>::iterator itr = obsqueue.begin();
			while (itr != obsqueue.end()) {
				if ((*itr).curt.time > sd_t.time)  break;
				item_sd.push_back(*itr);
				itr = obsqueue.erase(itr);
			}
		}
		def_unlock(&_mutex);

		/// here to sending the observation 
		list<RtConvItem>::iterator itr = item_sd.begin();
		while (itr != item_sd.end()) {
			char strbuf[1024] = { 0 };
			time2str((*itr).curt, strbuf, 2);
			printf("Thread %010d,begin to sending observation %s\n", def_pthread_id_self(), strbuf);
			/// generate the buff here
			if(this->type == ConvType::obsbnc){
				memset(buff, 0, sizeof(buff));
				makeupTimeBuffer((*itr).curt, (char*)buff, nbyte);
				strwrite(&svr, buff, nbyte);
				std::map<std::string, std::map<std::string, RtSatObs> >::iterator obitr_tosend = (*itr).obslist.begin();
				while (obitr_tosend != (*itr).obslist.end()) {
					map<string, RtSatObs>::iterator satitr = (*obitr_tosend).second.begin();
					while (satitr != (*obitr_tosend).second.end()) {
						memset(buff, 0, sizeof(buff));
						makeupObsBuffer((*obitr_tosend).first.c_str(), (*satitr).first.c_str(), (*satitr).second, (char*)buff, nbyte);
						strwrite(&svr, buff, nbyte);
						++satitr;
					}
					obitr_tosend++;
				}
			}
			else{
				std::map<std::string, std::map<std::string, RtSatObs> >::iterator obitr_tosend = (*itr).obslist.begin();
				while (obitr_tosend != (*itr).obslist.end()) {
					if(rtcmexEncoder.m_encodeRtcmExMsg((*obitr_tosend).first.c_str(),(*itr).curt,(*obitr_tosend).second)){
						strwrite(&svr,(unsigned char*)rtcmexEncoder.buff,rtcmexEncoder.nbyte);
					}
					obitr_tosend++;
				}
			}
			++itr;
		}
		sleepms(NINT(1000 / Deploy::s_getSpeed()));
	}
	strclose(&svr);
}
void RtConverter::makeupTimeBuffer(gtime_t tt, char* buffin, int& nbyte) {
	int week;
	double sow;
	sow = time2gpst(tt, &week);
	sprintf(buffin, "> %04d %14.7lf\r\n", week, sow);
	nbyte = strlen(buffin);
	return;
}
void RtConverter::makeupObsBuffer(const char* sitname, const char* cprn, RtSatObs& obs, char* buffin, int& nbyte) {
	int ifreq, hasobs = false;
	string sit_out = string(sitname);
	transform(sit_out.begin(), sit_out.end(), sit_out.begin(), ::toupper);
	sprintf(buffin, "%s %s", sit_out.c_str(), cprn);
	for (ifreq = 0; ifreq < MAXFREQ; ifreq++) {
		hasobs = false;
		if (obs.obs[ifreq] != 0.0) {
			hasobs = true;
			sprintf(buffin + strlen(buffin), " %s %14.3lf", obs.fob[ifreq].c_str(), obs.obs[ifreq]);
		}
		if (obs.obs[MAXFREQ + ifreq] != 0.0) {
			hasobs = true;
			sprintf(buffin + strlen(buffin), " %s %14.3lf", obs.fob[MAXFREQ + ifreq].c_str(), obs.obs[MAXFREQ + ifreq]);
		}
		if (hasobs) {
			sprintf(buffin + strlen(buffin), " S%s %7.3lf",obs.fob[MAXFREQ + ifreq].c_str() + 1, obs.snr[ifreq]);
		}	
	}
	strcat(buffin, "\r\n");
	nbyte = strlen(buffin);
	return;
}
void RtConverter::inputBrdm_G(string cprn, list<GPSEPH>& list_in) {
	if (list_in.size() == 0) return;
	list<GPSEPH>::iterator ephItr;
	def_lock(&_mutex);
	for (ephItr = list_in.begin(); ephItr != list_in.end(); ++ephItr) {
		map<string, list<GPSEPH> >::iterator mapItr = ephqueue_G.find(cprn);
		if (mapItr == ephqueue_G.end()) {
			list<GPSEPH> ltp;
			ltp.push_back(*ephItr);
			ephqueue_G[cprn] = ltp;
			continue;
		}
		/// only contains 1 value,adapt the position,eg: 0 pos for 2:00,1 pos for 0:00
		if ((*mapItr).second.size() < 2) {
			double dt = ((*ephItr).mjd - ephqueue_G[cprn].front().mjd) * 86400.0 +
				(*ephItr).sod - ephqueue_G[cprn].front().sod;
			if (dt > 0.01)
				ephqueue_G[cprn].push_front(*ephItr);
			else if(dt < -0.01)
				ephqueue_G[cprn].push_back(*ephItr);
			continue;
		}
		/// if contains two value,will check
		double dt = ((*ephItr).mjd - ephqueue_G[cprn].front().mjd) * 86400.0 +
			(*ephItr).sod - ephqueue_G[cprn].front().sod;
		if (dt > 0.01) {
			ephqueue_G[cprn].pop_back();
			ephqueue_G[cprn].push_front(*ephItr);
		}else if (dt < -0.01) {
			double dt = ((*ephItr).mjd - ephqueue_G[cprn].back().mjd) * 86400.0 +
				(*ephItr).sod - ephqueue_G[cprn].back().sod;
			if (dt > 0.01) {
				ephqueue_G[cprn].pop_back();
				ephqueue_G[cprn].push_back(*ephItr);
			}
		}
	}
	def_unlock(&_mutex);
}
void RtConverter::inputBrdm_R(string cprn,list<GLSEPH>& list_in) {
	if (list_in.size() == 0) return;
	list<GLSEPH>::iterator ephItr;
	def_lock(&_mutex);
	for (ephItr = list_in.begin(); ephItr != list_in.end(); ++ephItr) {
		map<string, list<GLSEPH> >::iterator mapItr = ephqueue_R.find(cprn);
		if (mapItr == ephqueue_R.end()) {
			list<GLSEPH> ltp;
			ltp.push_back(*ephItr);
			ephqueue_R[cprn] = ltp;
			continue;
		}
		/// only contains 1 value,adapt the position,eg: 0 pos for 2:00,1 pos for 0:00
		if ((*mapItr).second.size() < 2) {
			double dt = ((*ephItr).mjd - ephqueue_R[cprn].front().mjd) * 86400.0 +
				(*ephItr).sod - ephqueue_R[cprn].front().sod;
			if (dt > 0.01)
				ephqueue_R[cprn].push_front(*ephItr);
			else if (dt < -0.01)
				ephqueue_R[cprn].push_back(*ephItr);
			continue;
		}
		/// if contains two value,will check
		double dt = ((*ephItr).mjd - ephqueue_R[cprn].front().mjd) * 86400.0 +
			(*ephItr).sod - ephqueue_R[cprn].front().sod;
		if (dt > 0.01) {
			ephqueue_R[cprn].pop_back();
			ephqueue_R[cprn].push_front(*ephItr);
		}
		else if (dt < -0.01) {
			double dt = ((*ephItr).mjd - ephqueue_R[cprn].back().mjd) * 86400.0 +
				(*ephItr).sod - ephqueue_R[cprn].back().sod;
			if (dt > 0.01) {
				ephqueue_R[cprn].pop_back();
				ephqueue_R[cprn].push_back(*ephItr);
			}
		}
	}
	def_unlock(&_mutex);
}
void RtConverter::makeupEphBuffer_G(GPSEPH& gpsEph, char* buffin, int& nbyte) {
	int iy, im, id, ih, imin, isys;
	double sec;
	isys = index_string(SYS, gpsEph.cprn[0]);
	mjd2date(gpsEph.mjd, gpsEph.sod, &iy, &im, &id, &ih, &imin, &sec);
	sprintf(buffin + strlen(buffin),"%3s %04d %02d %02d %02d %02d %02d%19.12e%19.12e%19.12e\r\n", gpsEph.cprn, iy, im, id, ih, imin, (int)sec, gpsEph.a0,
		gpsEph.a1, gpsEph.a2);
	sprintf(buffin + strlen(buffin),"    %19.12e%19.12e%19.12e%19.12e\r\n",(double)gpsEph.aode, gpsEph.crs, gpsEph.dn, gpsEph.m0);
	sprintf(buffin + strlen(buffin),"    %19.12e%19.12e%19.12e%19.12e\r\n", gpsEph.cuc, gpsEph.e, gpsEph.cus, gpsEph.roota);
	sprintf(buffin + strlen(buffin),"    %19.12e%19.12e%19.12e%19.12e\r\n", gpsEph.toe, gpsEph.cic, gpsEph.omega0, gpsEph.cis);
	sprintf(buffin + strlen(buffin),"    %19.12e%19.12e%19.12e%19.12e\r\n", gpsEph.i0, gpsEph.crc, gpsEph.omega, gpsEph.omegadot);
	sprintf(buffin + strlen(buffin),"    %19.12e%19.12e%19.12e%19.12e\r\n", gpsEph.idot, gpsEph.resvd0, gpsEph.week, gpsEph.resvd1);
	sprintf(buffin + strlen(buffin),"    %19.12e%19.12e%19.12e%19.12e\r\n", gpsEph.accu, gpsEph.hlth, gpsEph.tgd, SYS[isys] == 'C' ? gpsEph.tgd1 : gpsEph.aodc);
	if (SYS[isys] == 'E')
		sprintf(buffin + strlen(buffin), "    %19.12e\r\n",9.999e8);
	else
		sprintf(buffin + strlen(buffin), "    %19.12e%19.12e\r\n", 0.0, 0.0);
	nbyte = strlen(buffin);
}
void RtConverter::makeupEphBuffer_R(GLSEPH& glsEph, char* buffin, int& nbyte) {
	int iy, im, id, ih, imin, isys;
	double sec;
	mjd2date(glsEph.mjd, glsEph.sod, &iy, &im, &id, &ih, &imin, &sec);
	sprintf(buffin + strlen(buffin), "%3s %04d %02d %02d %02d %02d %02d%19.12e%19.12e%19.12e\r\n", glsEph.cprn, iy, im, id, ih, imin, (int)sec, glsEph.tau, glsEph.gamma, glsEph.tk);
	sprintf(buffin + strlen(buffin), "    %19.12e%19.12e%19.12e%19.12e\r\n", glsEph.pos[0], glsEph.vel[0],glsEph.acc[0], glsEph.health);
	sprintf(buffin + strlen(buffin), "    %19.12e%19.12e%19.12e%19.12e\r\n", glsEph.pos[1], glsEph.vel[1],glsEph.acc[1], glsEph.frenum);
	sprintf(buffin + strlen(buffin), "    %19.12e%19.12e%19.12e%19.12e\r\n", glsEph.pos[2], glsEph.vel[2],glsEph.acc[2], glsEph.age);
	nbyte = strlen(buffin);
}
void RtConverter::routineEph() {
	time_t now;
	int nbyte;
	char cmd[256];
	unsigned char buff[1024];
	sprintf(cmd, ":%d", port);
	strinit(&svr);
	stropen(&svr, STR_TCPSVR, STR_MODE_RW, cmd);
	while (true) {
		/// here to sending the ephemeris
		def_lock(&_mutex);
		std::map<string, list<GPSEPH> >::iterator mapItr_G;
		std::map<string, list<GLSEPH> >::iterator mapItr_R;
		for (mapItr_G = ephqueue_G.begin(); mapItr_G != ephqueue_G.end(); ++mapItr_G) {
			list<GPSEPH>::iterator ephItr;
			for (ephItr = (*mapItr_G).second.begin(); ephItr != (*mapItr_G).second.end(); ++ephItr) {
				memset(buff, 0, sizeof(buff));
				makeupEphBuffer_G(*ephItr, (char*)buff, nbyte);
				strwrite(&svr, buff, nbyte);
			}
		}
		for (mapItr_R = ephqueue_R.begin(); mapItr_R != ephqueue_R.end(); ++mapItr_R) {
			list<GLSEPH>::iterator ephItr;
			for (ephItr = (*mapItr_R).second.begin(); ephItr != (*mapItr_R).second.end(); ++ephItr) {
				memset(buff, 0, sizeof(buff));
				makeupEphBuffer_R(*ephItr, (char*)buff, nbyte);
				strwrite(&svr, buff, nbyte);
			}
		}
		def_unlock(&_mutex);
		sleepms(5000);
	}
	strclose(&svr);
}
int main(int argc, char* args[]) {
	signal(SIGPIPE, SIG_IGN);
	Deploy::s_initInstance(argc,args);
	list<RtConverter*> rt_svrs;
	Deploy dly = Deploy::s_getConfigures();
	strinitcom();

	RtConverter rtconv_nrtk(5,dly.nrtkport,RtConverter::ConvType::obsbin);
	RtConverter rtconv_rtk(10,dly.rtkport, RtConverter::ConvType::obsbnc);
	RtConverter rtconv_eph(1, dly.ephport, RtConverter::ConvType::eph);
	RtRinexStream poststr;
	VrsNtripClient vrsstr;
	rt_svrs.push_back(&rtconv_nrtk);
	rt_svrs.push_back(&rtconv_rtk);
	rt_svrs.push_back(&rtconv_eph);

	/// post observation 
	poststr.openStream(rt_svrs);
	/// vrs observation 
	vrsstr.openStream(rt_svrs);

	/// server thread 
	rtconv_nrtk.beginProcess();
	rtconv_rtk.beginProcess();
	rtconv_eph.beginProcess();
	while (true) {
		Deploy::s_updateConfigures();
		sleepms(5000);
	}
}
