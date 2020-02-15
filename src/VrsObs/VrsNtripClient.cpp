#include "../../include/VrsObs/VrsNtripClient.h"
#include "../../include/Rnxobs/Com.h"
using namespace std;
using namespace bamboo;
VrsNtripClient::~VrsNtripClient() {
	/// will delete all the memory
	lcont = false;
}
void VrsNtripClient::openStream() {
	lcont = true;
	configs_sav = Deploy::s_getConfigures();
	time(&lastCheck);
	def_thread_t pid_handle;
#ifdef _WIN32
	CreateThread(NULL, 0, s_pthVrsStream_win, this, 0, NULL);
#else
	if (0 != pthread_create(&pid_handle, NULL, &s_pthVrsStream, this)) {
		cout
			<< "***ERROR(s_VRSMain):cant create thread to handle request!"
			<< endl;
		exit(1);
	}
#endif
}
void VrsNtripClient::openStream(list<RtConverter*> svrs_in) {
	lcont = true;
	configs_sav = Deploy::s_getConfigures();
	time(&lastCheck);
	def_thread_t pid_handle;
#ifdef _WIN32
	CreateThread(NULL, 0, s_pthVrsStream_win, this, 0, NULL);
#else
	if (0 != pthread_create(&pid_handle, NULL, &s_pthVrsStream, this)) {
		cout
			<< "***ERROR(s_VRSMain):cant create thread to handle request!"
			<< endl;
		exit(1);
	}
#endif
	this->m_svrs = svrs_in;
}
#ifdef _WIN32
DWORD WINAPI VrsNtripClient::s_pthVrsStream_win(LPVOID lp) {
	VrsNtripClient* svr = (VrsNtripClient*)lp;
	svr->m_routine();
	return NULL;
}
#endif
void* VrsNtripClient::s_pthVrsStream(void* args) {
	VrsNtripClient* svr = (VrsNtripClient*)args;
	svr->m_routine();
	return NULL;
}
void VrsNtripClient::m_routine() {
	time_t tt;
	bool lcheck;
	int nread = 0, i;
	char cmd[1024] = { 0 }, buff[1024] = { 0 }, strbuf[1024] = { 0 };
	list<VrsStaItem>::iterator staItr;
	list<RtConverter*>::iterator rtItr;
	map<string, list<VrsStaItem>>::iterator mapItr;
	map<string, stream_t*>::iterator streamItr;
	map<string, rtcm_t*>::iterator rtcmItr;
	map<string, VrsStaItem>::iterator vrsItr;
	for (mapItr = configs_sav.rt_mounts.begin(); mapItr != configs_sav.rt_mounts.end(); ++mapItr) {
		for (staItr = (*mapItr).second.begin(); staItr != (*mapItr).second.end(); ++staItr) {
			m_newConnection((*staItr));
		}
	}
	while (lcont) {
		lcheck = false;
		time(&tt);
		if (tt - lastCheck > 5) {
			lcheck = true;
			lastCheck = tt;
			m_adaptConfigures();
		}
		for (streamItr = m_streams.begin(); streamItr != m_streams.end(); ++streamItr) {
			if (m_stas[(*streamItr).first].type == VrsStaItem::stream_type::vrs && lcheck) {
				//// sending GGA here
				m_sendGGAReq((*streamItr).second, m_stas[(*streamItr).first]);
			}
			nread = strread((*streamItr).second, (unsigned char*)buff, 1024);
			for (i = 0; i < nread; i++) {
				if (input_rtcm3(m_rtcms[(*streamItr).first], buff[i]) == 1) {
					RtConvItem curitem = m_makeupItems((*streamItr).first, m_rtcms[(*streamItr).first]);
					if (curitem.obslist.size() > 0) {
						for (rtItr = m_svrs.begin(); rtItr != m_svrs.end(); ++rtItr) {
							(*rtItr)->inputObs(curitem);
						}
					}
				}
			}
		}
		sleepms(50);
	}
	/// will delete memory here
	for (streamItr = m_streams.begin(); streamItr != m_streams.end(); ++streamItr) {
		strclose((*streamItr).second);
		free((*streamItr).second);
	}
	for (rtcmItr = m_rtcms.begin(); rtcmItr != m_rtcms.end(); ++rtcmItr) {
		free_rtcm((*rtcmItr).second);
		free((*rtcmItr).second);
	}
	m_streams.clear();
	m_rtcms.clear();
	m_stas.clear();
}
void VrsNtripClient::m_sendGGAReq(stream_t* stream, VrsStaItem& item) {
	char buff[1024] = { 0 };
	gtime_t gtime = { 0 };
	double h, ep[6], pos[3], dms1[3], dms2[3], dop = 1.0;
	int solq = 1,nsend;
	char *p = (char *)buff, *q, sum;
	gtime.time = time(NULL);
	time2epoch(gtime, ep);
	ecef2pos(item.x, pos);
	deg2dms(fabs(pos[0])*R2D, dms1, 7);
	deg2dms(fabs(pos[1])*R2D, dms2, 7);
	h = pos[2];
	p += sprintf(p, "$GPGGA,%02.0f%02.0f%05.2f,%02.0f%010.7f,%s,%03.0f%010.7f,%s,%d,%02d,%.1f,%.3f,M,%.3f,M,%.1f,",
		ep[3], ep[4], ep[5], dms1[0], dms1[1] + dms1[2] / 60.0, pos[0] >= 0 ? "N" : "S",
		dms2[0], dms2[1] + dms2[2] / 60.0, pos[1] >= 0 ? "E" : "W", solq,
		0, dop, pos[2] - h, h, 1.0);

	for (q = (char *)buff + 1, sum = 0; *q; q++) sum ^= *q; /* check-sum */
	p += sprintf(p, "*%02X%c%c", sum, 0x0D, 0x0A);
	nsend =  p - (char *)buff;
	strwrite(stream, (unsigned char*)buff, nsend);
	return;
}
RtConvItem VrsNtripClient::m_makeupItems(string staname, rtcm_t* rtcm) {
	int iobs, isys, ifreq, gotobs = false;
	RtConvItem item;
	char cprn[16] = { 0 }, code[4] = { 0 };
	map<string, RtSatObs> mapObs;
	for (iobs = 0; iobs < rtcm->obs.n; iobs++) {
		RtSatObs sat;
		satno2id(rtcm->obs.data[iobs].sat, cprn);
		if (-1 == pointer_string(configs_sav.nprn, configs_sav.cprn, cprn))
			continue;
		isys = index_string(SYS, cprn[0]);
		for (ifreq = 0; ifreq < sizeof(rtcm->obs.data[iobs].L) / sizeof(double); ifreq++) {
			if (rtcm->obs.data[iobs].P[ifreq] != 0.0)
				gotobs = true;
			strcpy(code, code2obs(rtcm->obs.data[iobs].code[ifreq], NULL));
			sat.obs[ifreq] = rtcm->obs.data[iobs].L[ifreq];
			sat.fob[ifreq] = "L" + string(code);

			sat.obs[MAXFREQ + ifreq] = rtcm->obs.data[iobs].P[ifreq];
			sat.fob[MAXFREQ + ifreq] = "C" + string(code);

			sat.snr[ifreq] = rtcm->obs.data[iobs].SNR[ifreq] / 4;
		}
		if (gotobs) mapObs[cprn] = sat;
	}
	item.curt = rtcm->time;
	if (mapObs.size() > 0)
		item.obslist[staname] = mapObs;
	return item;
}
void VrsNtripClient::m_deleteConnection(string staname) {
	map<string, stream_t*>::iterator streamItr;
	map<string, rtcm_t*>::iterator rtcmItr;
	map<string, VrsStaItem>::iterator vrsItr;

	streamItr = m_streams.find(staname);
	rtcmItr = m_rtcms.find(staname);
	vrsItr = m_stas.find(staname);

	if (streamItr != m_streams.end()) {
		strclose((*streamItr).second);
		free((*streamItr).second);
		m_streams.erase(streamItr);
	}
	if (rtcmItr != m_rtcms.end()) {
		free_rtcm((*rtcmItr).second);
		free((*rtcmItr).second);
		m_rtcms.erase(rtcmItr);
	}
	if (vrsItr != m_stas.end())
		m_stas.erase(vrsItr);
}
void VrsNtripClient::m_newConnection(VrsStaItem& item) {
	stream_t* str_m;
	rtcm_t* rtcm_m;
	str_m = (stream_t*)calloc(1, sizeof(stream_t));
	strinit(str_m);
	rtcm_m = (rtcm_t*)calloc(1, sizeof(rtcm_t));
	init_rtcm(rtcm_m);
	stropen(str_m, STR_NTRIPCLI, STR_MODE_RW, item.strpath.c_str());
	strsettimeout(str_m, 60000, 10000); /// 60s for timeout 10s for reconnect
	m_streams[item.staname] = str_m;
	m_rtcms[item.staname] = rtcm_m;
	m_stas[item.staname] = item;
}
void VrsNtripClient::m_adaptConfigures() {
	Deploy config_new = Deploy::s_getConfigures();
	list<VrsStaItem>::iterator staItr;
	list<RtConverter*>::iterator rtItr;
	map<string, list<VrsStaItem>>::iterator mapItr;
	map<string, VrsStaItem>::iterator vrsItr;
	for (vrsItr = m_stas.begin(); vrsItr != m_stas.end(); ++vrsItr) {
		(*vrsItr).second.lcheck = false;
	}
	for (mapItr = config_new.rt_mounts.begin(); mapItr != config_new.rt_mounts.end(); ++mapItr) {
		for (staItr = (*mapItr).second.begin(); staItr != (*mapItr).second.end(); ++staItr) {
			if (m_stas.find((*staItr).staname) == m_stas.end()) {
				/// new station
				m_newConnection((*staItr));
			}
			else {
				/// already in memory,will check the configure
				m_stas[(*staItr).staname].lcheck = true;
				VrsStaItem& item_mem = m_stas[(*staItr).staname];
				if (item_mem == (*staItr))
					continue;
				// delete first 
				m_deleteConnection((*staItr).staname);
				// add new here
				m_newConnection((*staItr));
			}
		}
	}
	list<string> _dels;
	list<string>::iterator strItr;
	for (vrsItr = m_stas.begin(); vrsItr != m_stas.end(); ++vrsItr) {
		if ((*vrsItr).second.lcheck == false) {
			_dels.push_back((*vrsItr).first);
		}
	}
	for (strItr = _dels.begin(); strItr != _dels.end(); ++strItr) {
		m_deleteConnection(*strItr);
	}
	configs_sav = config_new;
}
