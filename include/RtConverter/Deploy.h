#ifndef DEPLOY_H_
#define DEPLOY_H_
#include "Const.h"
#include <string>
#include <map>
#include <list>
#include "../Rnxobs/Log.h"
using namespace std;
namespace bamboo {
	class VrsStaItem{
	public:
		enum stream_type{
			vrs, stream
		};
		int type; // identify whether it is vrs or mountpoint
		double x[3];
		string staname;
		string strpath;
		int lcheck;
	};
	inline bool operator ==(const VrsStaItem & v1, const VrsStaItem & v2) {
		return v1.type == v2.type && fabs(v1.x[0] - v2.x[0]) < 0.01 &&
			fabs(v1.x[1] - v2.x[1]) < 0.01 && fabs(v1.x[2] - v2.x[2]) < 0.01 && v1.strpath == v2.strpath && v1.staname == v2.staname;
	}
	class Deploy {
	public:
		Deploy();
		inline static void s_initInstance(int argc,char* args[]) {
			if (sInstance == NULL)
				sInstance = new Deploy(argc,args);
		}
		inline static int s_getSpeed() {
			return sInstance->speed;
		}
		~Deploy();
		static Deploy s_getConfigures();
		static bool s_updateConfigures ();
		int mjd0,mjd1,speed;
		double sod0, sod1,seslen,dintv;
		int nfreq[MAXSYS],nprn,rtkport,nrtkport,ephport;
		char freq[MAXSYS][MAXFREQ][LEN_FREQ];
		char obsdir[256],outdir[256],ephdir[1024];
		string cprn[MAXSAT];
		time_t lastAct;
		list<string> post_stalist;
		map<string, list<VrsStaItem>> rt_mounts;
	protected:
		Deploy(int,char*[]);
	
		bool m_checkStation(string);
		void m_readConfiguresJson(bool);
		char f_jsonConfigures[1024];
		static Deploy* sInstance;
		static def_lock_t s_mutex;
		static int s_count;
	};
}
#endif


