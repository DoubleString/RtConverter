/*
 * RtcmExEncoder.cpp
 *
 *  Created on: 2020/01/20 17:15:25
 *      Author: juntao
 */

#include "../../include/Encoder/RtcmExEncoder.h"
#include "../../include/Rtklib/rtklib_fun.h"
#include "../../include/Rnxobs/Com.h"
using namespace std;
using namespace bamboo;
RtcmExEncoder::RtcmExEncoder(){
	nbyte = len = nbit = 0;
	init_rtcm(&m_rtcm);
}
RtcmExEncoder::~RtcmExEncoder(){
	free_rtcm(&m_rtcm);
}
/*            +----------+--------+-----------+--------------------+----------+
 *            | preamble |   length   |    data message    |  parity  |
 *            +----------+--------+-----------+--------------------+----------+
 *            |<-- 8 --->|<--- 16 --->|<--- length x 8 --->|<-- 32 -->| */
int RtcmExEncoder::m_encodeRtcmExMsg(const char* sitname,gtime_t tt,map<string, RtSatObs>& obs){
	unsigned int crc;
	int i = 0;
	nbyte = len = nbit = 0;
	setbitu((unsigned char*)buff, i, 8, RTCMEX_PREAMB); i += 8;
	setbitu((unsigned char*)buff, i, 16, 0);        i += 16;
	if(0 == (nbit = m_encodeObsMsg(sitname,tt,obs))){
		nbyte = nbit = len = 0;
		return 0;
	}
	/* padding to align 8 bit boundary */
	for (i = nbit; i % 8; i++) {
		setbitu((unsigned char*)buff, i, 1, 0);
	}
	/* message length (header+data) (bytes) */
	if ((len = i / 8) > sizeof(buff) - 4) {
		nbyte = nbit = len = 0;
		return 0;
	}
	/* message length without header and parity */
	setbitu((unsigned char*)buff, 8, 16, len - 3);

	/* crc-24q */
	crc = rtk_crc24q((unsigned char*)buff, len);
	setbitu((unsigned char*)buff, i, 32, crc);
	/* length total (bytes) */
	nbyte = len + 4;
	return 1;
}
int RtcmExEncoder::m_encodeObsMsg(const char* sitname,gtime_t tt,map<string, RtSatObs>& obs_in){
	int isys,ifreq,gotobs[MAXSYS] = { 0 },lsync,wk;
	char code[8];
	time2gpst(tt,&wk);
	strncpy(buff + 3,sitname,4);
	setbitu((unsigned char*)buff, 7 * 8, 16, wk); /*set the gps week here*/
	m_rtcm.obs.n = 0;
	m_rtcm.time = tt;
	map<string,RtSatObs>::iterator mapItr;
	for(mapItr = obs_in.begin();mapItr != obs_in.end();++mapItr){
		string cprn = (*mapItr).first;
		RtSatObs& obs = (*mapItr).second;
		isys = index_string(SYS, cprn[0]);
		gotobs[isys] = 1;
		m_rtcm.obs.data[m_rtcm.obs.n].sat = satid2no(cprn.c_str());
		if(m_rtcm.obs.data[m_rtcm.obs.n].sat == 0) continue;
		for(ifreq = 0;ifreq < MAXFREQ;ifreq++){
			m_rtcm.obs.data[m_rtcm.obs.n].L[ifreq] = obs.obs[ifreq];
			m_rtcm.obs.data[m_rtcm.obs.n].P[ifreq] = obs.obs[MAXFREQ + ifreq];
			m_rtcm.obs.data[m_rtcm.obs.n].SNR[ifreq] = static_cast<unsigned char>(obs.snr[ifreq] / 0.25);
			strcpy(code, obs.fob[MAXFREQ + ifreq].c_str());
			if (cprn[0] == 'C' && code[1] == '2')
				code[1] = '1';
			if (cprn[0] == 'G' && code[1] == '2')
				strcpy(code, "L2W");
			m_rtcm.obs.data[m_rtcm.obs.n].code[ifreq] = obs2code(code + 1, NULL);
		}
		m_rtcm.obs.n++;
	}
	nbyte = 3 + 4 + 2; // pre,len,sitname
	/////////////////// GPS /////////////////
	if (gotobs[index_string(SYS, 'G')]) {
		lsync = gotobs[index_string(SYS, 'C')] || gotobs[index_string(SYS, 'R')] || gotobs[index_string(SYS, 'E')];
		gen_rtcm3(&m_rtcm, 1074, lsync);
		if(nbyte + m_rtcm.nbyte >  sizeof(buff) - 4) return 0;
		memcpy(buff + nbyte, m_rtcm.buff, m_rtcm.nbyte);
		nbyte = nbyte + m_rtcm.nbyte;
	}
	////////////////// BDS ////////////////
	if (gotobs[index_string(SYS, 'C')]) {
		lsync = gotobs[index_string(SYS, 'R')] || gotobs[index_string(SYS, 'E')];
		gen_rtcm3(&m_rtcm, 1124, lsync);
		if(nbyte + m_rtcm.nbyte >  sizeof(buff) - 4) return 0;
		memcpy(buff + nbyte, m_rtcm.buff, m_rtcm.nbyte);
		nbyte = nbyte + m_rtcm.nbyte;
	}
	//////////////////// GLO ////////////////
	if (gotobs[index_string(SYS, 'R')]) {
		lsync = gotobs[index_string(SYS, 'E')];
		gen_rtcm3(&m_rtcm, 1084, lsync);
		if(nbyte + m_rtcm.nbyte >  sizeof(buff) - 4) return 0;
		memcpy(buff + nbyte, m_rtcm.buff, m_rtcm.nbyte);
		nbyte = nbyte + m_rtcm.nbyte;
	}
	if (gotobs[index_string(SYS, 'E')]) {
		gen_rtcm3(&m_rtcm, 1094, 0);
		if(nbyte + m_rtcm.nbyte >  sizeof(buff) - 4) return 0;
		memcpy(buff + nbyte, m_rtcm.buff, m_rtcm.nbyte);
		nbyte = nbyte + m_rtcm.nbyte;
	}
	return nbyte * 8;
}






