#ifndef VRSNTRIPCLIENT_H_
#define VRSNTRIPCLIENT_H_
#include "../RtConverter/Deploy.h"
#include "../RtConverter/RtConverter.h"
#include "../Rtklib/rtklib_fun.h"
#include <list>
#include <map>
using namespace std;
namespace bamboo {
	class VrsNtripClient {
	public:
		~VrsNtripClient();
		void openStream();
		void openStream(list<RtConverter*>);
	protected:
#ifdef _WIN32
		static DWORD WINAPI s_pthVrsStream_win(LPVOID lp);
#endif
		RtConvItem m_makeupItems(string,rtcm_t*);
		static void* s_pthVrsStream(void*);
		void m_routine();
		void m_sendGGAReq(stream_t*,VrsStaItem&);
		void m_adaptConfigures();


		void m_newConnection(VrsStaItem&);
		void m_deleteConnection(string staname);

		bool lcont;
		time_t lastCheck;
		Deploy configs_sav;
		list<RtConverter*> m_svrs;
		map<string, stream_t*> m_streams;
		map<string, rtcm_t*> m_rtcms;
		map<string, VrsStaItem> m_stas;
	};
}
#endif