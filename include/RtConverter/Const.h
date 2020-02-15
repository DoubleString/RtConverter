#ifndef INCLUDE_BAMBOO_CONST_H_
#define INCLUDE_BAMBOO_CONST_H_
#include <string>
#include <cmath>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <time.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <map>
#include <list>
#include <iostream>
#ifdef _WIN32
#include <winsock2.h>
#include <windows.h>
#include <condition_variable>
#else
#include <pthread.h>
#include <unistd.h>
#endif
#include "../Rtklib/rtklib.h"
using namespace std;
namespace bamboo {
	const unsigned char RTCM2PREAMB = 0x66;       /* rtcm ver.2 frame preamble */
	const unsigned char RTCM3PREAMB = 0xD3;      /* rtcm ver.3 frame preamble */
	const unsigned char ATOMPREAMB = 0xDB;      /* atom preamble */
	const int MAXFREQ = 3;
	const int MAXSYS = 7;
	const int MAXSAT = 170;
	const int MAXSIT = 500;
	const int MAXOC = 120;
	const int MAXOBSTYP = 50;
	const int MAXPORT = 10;
	const int MAXICS = 21;
	const int MAXECLIP = 120;
	const int MAXCOEFPAR = 3 + MAXSYS + 1 + 1 + 2 + 1; // XYZ RECCLK ZTD AMB ION DCB (ION maybe two)
	const int MAXDGR = 7;
	const int MAXPNT = 2 * (MAXDGR + 1);
	const int MAXVARS = 3 * (MAXICS + 1);

	const int ORB_NONE = 0;
	const int ORB_RTSP3 = 1;
	const int ORB_SP3 = 2;
	const int ORB_BRD = 3;
	const int ORB_RTBRD = 4;
	const int ORB_PREBRD = 5;

	const int CLK_NONE = 0;
	const int CLK_BRD = 6;
	const int CLK_RTBRD = 7;
	const int CLK_PREBRD = 8;
	const int CLK_SP3 = 9;



	const int MINCMNSIT = 4;
	const int MAX_GRID = 73;
	const int MeanEarthRadius = 6371000;
	const int IONO_SHELL_H = 450000;
	const int MAXTABLE = 20;
	const int MAXBASE = 20;

	const int MAXREFBLINE = 20;
	const int MAXBLINELEN = 1024;

	const char SYS[] = "GRECSJI";
	const char SYSNAME[][4] = { "gps","glo","gal","bds","","","" };
	const double MAXWND = 1E-2;
	const char OBSTYPE[] = "PWCIXSAQLDBYMZN "; //"WPCIXSAQLDBYMZN ";

	const double PI = 3.1415926535897932384626433832795028841971693993;
	const double VEL_LIGHT = 299792458.0;
	const double E_MAJAXIS = 6378136.60;

	const double DEG2RAD = PI / 180.0;
	const double RAD2DEG = 180.0 / PI;

	const double ARCSEC2RAD = PI / (3600.0*180.0);
	const double CONE = 2.7182818284590452353602874713526624977572470936999;

	const double EARTH_A = 6378137;
	const double EARTH_ALPHA = 1 / 298.257223563;
	const double ESQUARE = (EARTH_A * EARTH_A - (EARTH_A - EARTH_A * EARTH_ALPHA) * (EARTH_A - EARTH_A * EARTH_ALPHA)) / (EARTH_A * EARTH_A);

	const double OFF_GPS2TAI = 19.0;
	const double OFF_TAI2TT = 32.184;
	const double OFF_GPS2TT = OFF_GPS2TAI + OFF_TAI2TT;
	const double OFF_MJD2JD = 2400000.50;

	const double E_ROTATE = 7.2921151467E-5;
	const double GME = 3.986004415E14; //TT-compatible
	const double GPS_FREQ = 10230000.0;
	const double GPS_L1 = 154 * GPS_FREQ;
	const double GPS_L2 = 120 * GPS_FREQ;
	const double GPS_L5 = 115 * GPS_FREQ;

	const double GLS_FREQ = 178000000.0;
	const double GLS_L1 = 9 * GLS_FREQ;
	const double GLS_L2 = 7 * GLS_FREQ;
	const double GLS_dL1 = 562500.0;
	const double GLS_dL2 = 437500.0;

	const double GAL_E1 = GPS_L1;
	const double GAL_E5 = 1191795000.0;
	const double GAL_E5a = GPS_L5;
	const double GAL_E5b = 1207140000.0;
	const double GAL_E6 = 1278750000.0;

	const double BDS_B1 = 1561098000.0;
	const double BDS_B2 = 1207140000.0;
	const double BDS_B3 = 1268520000.0;
	const double QZS_L1 = GPS_L1;
	const double QZS_L2 = GPS_L2;
	const double QZS_L5 = GPS_L5;
	const double QZS_LEX = 1278750000.0;

	const int LEN_OBSTYPE = 4;
	const int LEN_FREQ = 4;
	const int LEN_ANTENNA = 21;
	const int LEN_PRN = 4;
	const int LEN_PORT = 256;
	const int LEN_STRING = 1024;
	const int LEN_SATTYPE = 20;
	const int LEN_PCVTYPE = 5;
	const int LEN_CLKTYPE = 5;
	const int LEN_COSPARID = 7;
	const int LEN_SITENAME = 5;
	const int LEN_ZTDMAP = 3;
	const int LEN_EPHION = 4;
#ifdef _WIN32
#define def_thread_t    HANDLE
#define def_dev_t       HANDLE 
#define def_cond_t      CONDITION_VARIABLE 
#define def_lock_t      CRITICAL_SECTION
#define def_initlock(f) InitializeCriticalSection(f)
#define def_lock(f)     EnterCriticalSection(f)
#define def_unlock(f)   LeaveCriticalSection(f)
#define def_destroylock(f)  DeleteCriticalSection(f)

#define def_initcond(c)  InitializeConditionVariable(c)
#define def_notify(c)    WakeConditionVariable(c)
#define def_wait(c,f)    SleepConditionVariableCS(c,f,INFINITE)
#define def_destroycond(c) 

#define FILEPATHSEP '\\'
#define def_pthread_id_self()  GetCurrentThreadId()
#else

#define def_thread_t    pthread_t
#define def_dev_t       int 
#define def_cond_t      pthread_cond_t
#define def_lock_t      pthread_mutex_t
#define def_initlock(f) pthread_mutex_init(f,NULL)
#define def_lock(f)     pthread_mutex_lock(f)
#define def_unlock(f)   pthread_mutex_unlock(f)
#define def_destroylock(f)  pthread_mutex_destroy(f)

#define def_initcond(c) pthread_cond_init(c,NULL)
#define def_notify(c) pthread_cond_signal(c)
#define def_wait(c,f) pthread_cond_wait(c,f)
#define def_destroycond(c) pthread_cond_destroy(c)

#define FILEPATHSEP '/'

#define def_pthread_id_self()  pthread_self()
#endif


#define BDSTOINT(type, value) (type)(round(value))

#define BDSADDBITS(a, b) {bitbuffer = (bitbuffer<<(a)) \
                       |(BDSTOINT(long long,b)&((1ULL<<a)-1)); \
                       numbits += (a); \
                       while(numbits >= 8) { \
                       buffer[size++] = bitbuffer>>(numbits-8);numbits -= 8;}}
#define BDSADDBITSFLOAT(a,b,c) {long long i = BDSTOINT(long long,(b)/(c)); \
                             BDSADDBITS(a,i)};
	template<typename T> string toString(const T& t) {
		ostringstream oss;  //
		oss << t;             //
		return oss.str();
	}
	template<class T>
	inline T MAX(T a, T b) {
		return a > b ? a : b;
	}
	template<class T>
	inline T MIN(T a, T b) {
		return a > b ? b : a;
	}

	template<class T, class P>
	inline T SIGN(T a, P b) {
		return b > 0 ? a : -a;
	}
	template<class T>
	inline int NINT(T a) {
		return (int)(a + SIGN(1, a) * 0.5);
	}
	template<class T>
	inline T ABS(T a) {
		return a > 0 ? a : -a;
	}
}
#endif /* INCLUDE_BAMBOO_CONST_H_ */
