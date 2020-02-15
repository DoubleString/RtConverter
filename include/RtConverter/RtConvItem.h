#ifndef RTCONVITEM
#define RTCONVITEM
#include <map>
#include <string>
#include <time.h>
#include "../RtConverter/Const.h"
namespace bamboo {
	class RtSatObs {
	public:
		RtSatObs() { memset(obs, 0, sizeof(obs)); }
		double obs[MAXFREQ * 2];
		double snr[MAXFREQ];
		std::string fob[MAXFREQ * 2];
	};
	class RtConvItem {
	public:
		RtConvItem() { time(&gent); }
		time_t gent;
		gtime_t curt; /// current time
		std::map<std::string, std::map<std::string, RtSatObs> > obslist;
	};
}
#endif
