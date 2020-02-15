#ifndef RTRINEXSTREAM_H_
#define RTRINEXSTREAM_H_
#include "../RtConverter/Deploy.h"
#include "../Rnxobs/Rnxobs.h"
#include "../RtConverter/RtConvItem.h"
#include "../Rnxbrd/OrbitClk.h"
#include <list>
namespace bamboo {
	class RtConverter;
	class RtRinexStream {
	public:
		~RtRinexStream();
		void openStream();
		void openStream(list<RtConverter*>&);
		void closeStream();
		static Deploy* m_getCurConfigures() { return &configs_sav; }
	protected:
#ifdef _WIN32
		static DWORD WINAPI s_pthRniexStream_win(LPVOID lp);
#endif
		void m_adaptConfigures();
		void m_routinue();
		RtConvItem m_makeupItems(int mjd,double sod,vector<RnxobsFile_sat*>&);
		vector<RnxobsFile_sat*> m_rnxs;
		RnxEphFileAdapter* mOacs;
		bool lcont;
		list<RtConverter*> m_svrs;
		time_t lastCheck;
		bool m_makeupEph_G(list<GPSEPH>&,string cprn,int mjd,double sod);
		bool m_makeupEph_R(list<GLSEPH>&,string cprn,int mjd,double sod);
	protected:
		static void* s_pthRinexStream(void*);
		static Deploy configs_sav;
	};
}
#endif

