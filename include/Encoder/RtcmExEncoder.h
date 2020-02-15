/*
 * Rtcm3ExEncoder.h
 *
 *  Created on: 2020/01/20 17:09:35
 *      Author: juntao
 */

#ifndef INCLUDE_ENCODER_RTCMEXENCODER_H_
#define INCLUDE_ENCODER_RTCMEXENCODER_H_
#include <string>
#include <map>
#include "../Rtklib/rtklib_fun.h"
#include "../RtConverter/RtConvItem.h"
using namespace std;
namespace bamboo{
#define RTCMEX_PREAMB 0xDA
class RtcmExEncoder{
public:
	RtcmExEncoder();
	~RtcmExEncoder();
	int m_encodeRtcmExMsg(const char* sitname,gtime_t tt,map<string, RtSatObs>& obs);
	char buff[2048];
	int nbyte;
protected:
	int m_encodeObsMsg(const char* sitname,gtime_t tt,map<string, RtSatObs>& obs);
	int len,nbit;
	rtcm_t m_rtcm;
};
}
#endif /* INCLUDE_ENCODER_RTCMEXENCODER_H_ */
