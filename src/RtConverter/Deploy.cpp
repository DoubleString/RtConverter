#include "../../include/RtConverter/Deploy.h"
#include "../../include/Json/json.h"
#include "../../include/Rnxobs/Com.h"
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include <cstdio>
using namespace std;
using namespace bamboo;
Deploy* Deploy::sInstance = NULL;
int Deploy::s_count = 0;
def_lock_t Deploy::s_mutex;
Deploy::Deploy(int nargc,char* args[]) {
	/// using init value 
	int iargc;
	lastAct = 0;
	ephport = -1;
	if (s_count++ == 0)
		def_initlock(&s_mutex);
	for (iargc = 0; iargc < nargc; iargc++) {
		if (!strcmp(args[iargc], "-h") || !strcmp(args[iargc], "-help")) {
			cout << "Usage:RtConverter -conf configures.json" << endl;
			exit(1);
		}
		if (!strcmp(args[iargc], "-conf")) {
			strcpy(f_jsonConfigures, args[++iargc]);
		}
	}
	this->m_readConfiguresJson(true);
}
Deploy::Deploy() {
	lastAct = 0;
	ephport = -1;
	if (s_count++ == 0)
		def_initlock(&s_mutex);
}
Deploy::~Deploy() {
}
Deploy Deploy::s_getConfigures() {
	def_lock(&s_mutex);
	Deploy ret;
	ret = *sInstance;
	def_unlock(&s_mutex);
	return ret;
}
bool Deploy::s_updateConfigures() {
	struct stat st;
	int result,lupdate = false;
	result = stat(sInstance->f_jsonConfigures, &st);
	//��ʾ�ļ�״̬��Ϣ
	if (result != 0) {
		printf("file = %s,monitor file is not exist!", sInstance->f_jsonConfigures);
		return false;
	}
	if (sInstance->lastAct != st.st_mtime) {
		def_lock(&s_mutex);
		sInstance->m_readConfiguresJson(false);
		def_unlock(&s_mutex);
		lupdate = true;
	}
	return lupdate;
}
void Deploy::m_readConfiguresJson(bool lexit) {
	fstream f;
	struct stat st;
	int iy, im, id, ih, imin, i, isys, j, ivrs, nsys, result,isat;
	double dsec, seslen;
	string item;
	char value[1024] = { 0 }, freq[MAXSYS][LEN_STRING] = { 0 }, cmd[1024] = { 0 };

	Json::CharReaderBuilder rbuilder;
	rbuilder["collectComments"] = false;
	Json::Value root;
	JSONCPP_STRING errs;
	f.open(f_jsonConfigures, ios::in);
	if (!f.is_open()) {
		printf("file = %s,failed to open file", f_jsonConfigures);
		if (lexit) exit(1); else return;
	}
	bool parse_ok = Json::parseFromStream(rbuilder, f, &root, &errs);
	f.close();
	if (!parse_ok) {
		printf("file = %s,failed to parse file", f_jsonConfigures);
		if (lexit) exit(1); else return;
	}
	post_stalist.clear();
	rt_mounts.clear();
	try {
		item = "post-configures";
		bamboo::excludeAnnoValue(value, root[item]["time"].asCString());
		sscanf(value, "%d%d%d%d%d%lf", &iy, &im, &id, &ih, &imin, &dsec);
		yr2year(iy);
		bamboo::excludeAnnoValue(value, root[item]["seslen"].asCString());
		sscanf(value, "%lf", &seslen);
		this->mjd0 = md_julday(iy, im, id);
		this->sod0 = ih * 3600 + imin * 60 + dsec;
		this->mjd1 = (int)(this->mjd0 + (this->sod0 + seslen) / 86400);
		this->sod1 = this->sod0 + seslen
			- (this->mjd1 - this->mjd0) * 86400.0;
		bamboo::excludeAnnoValue(value, root[item]["interval"].asCString());
		sscanf(value, "%lf", &dintv);

		bamboo::excludeAnnoValue(value, root[item]["speed"].asCString());
		sscanf(value, "%d", &speed);

		bamboo::excludeAnnoValue(value, root[item]["freqused"].asCString());
		split_string(true, value, ' ', ' ', ' ', &nsys, (char*)freq,
			LEN_STRING);
		for (i = 0; i < nsys; i++) {
			isys = -1;
			j = index_string(freq[i], ':');
			if (strstr(freq[i], "GPS") != NULL)
				isys = index_string(SYS, 'G');
			else if (strstr(freq[i], "GLS") != NULL)
				isys = index_string(SYS, 'R');
			else if (strstr(freq[i], "GAL") != NULL)
				isys = index_string(SYS, 'E');
			else if (strstr(freq[i], "BDS") != NULL
				|| strstr(freq[i], "CPS") != NULL
				|| strstr(freq[i], "CMS") != NULL)
				isys = index_string(SYS, 'C');
			else if (strstr(freq[i], "QZS") != NULL)
				isys = index_string(SYS, 'J');
			else {
				printf("freq = %s: unknown GNSS system", freq[i]);
				exit(1);
			}
			if (isys != -1) {
				split_string(true, freq[i] + j + 1, ' ', ' ', '_',
					this->nfreq + isys, (char*)(this->freq[isys]),
					LEN_FREQ);
				for (j = 0; j < this->nfreq[isys]; j++) {
					if (strstr(this->freq[isys][j], "L1") != NULL
						|| strstr(this->freq[isys][j], "G1") != NULL
						|| strstr(this->freq[isys][j], "E1") != NULL)
						strcpy(this->freq[isys][j], "L1");
					else if (strstr(this->freq[isys][j], "L2") != NULL
						|| strstr(this->freq[isys][j], "G2") != NULL
						|| strstr(this->freq[isys][j], "E2") != NULL
						|| strstr(this->freq[isys][j], "B1") != NULL) {
						strcpy(this->freq[isys][j], "L2");
					}
					else if (strstr(this->freq[isys][j], "L5") != NULL
						|| strstr(this->freq[isys][j], "E5a") != NULL)
						strcpy(this->freq[isys][j], "L5");
					else if (strstr(this->freq[isys][j], "E6") != NULL
						|| strstr(this->freq[isys][j], "B3") != NULL
						|| strstr(this->freq[isys][j], "LEX") != NULL)
						strcpy(this->freq[isys][j], "L6");
					else if (strstr(this->freq[isys][j], "E5b") != NULL
						|| strstr(this->freq[isys][j], "B2") != NULL)
						strcpy(this->freq[isys][j], "L7");
					else if (strstr(this->freq[isys][j], "E5") != NULL)
						strcpy(this->freq[isys][j], "L8");
				}
			}
		}
		
		bamboo::excludeAnnoValue(obsdir, root[item]["obsdir"].asCString());
		bamboo::excludeAnnoValue(outdir, root[item]["outdir"].asCString());
		bamboo::excludeAnnoValue(ephdir, root[item]["brdmdir"].asCString());
		for (i = 0; i < root[item]["sta-list"].size(); i++) {
			bamboo::excludeAnnoValue(value, root[item]["sta-list"][i].asCString());
			if (!m_checkStation(value)) {
				cout << "Existing same station,will continue for " << value << endl;
				continue;
			}
			post_stalist.push_back(value);
		}
		item = "satellites";
		const char* system[] = { "GPS","BDS","GAL","GLO" };
		nprn = 0;
		cprn->clear();
		for (isys = 0; isys < sizeof(system) / sizeof(char*); isys++) {
			for (isat = 0; isat < root[item][system[isys]].size(); isat++) {
				cprn[nprn++] = string(root[item][system[isys]][isat].asCString());
			}
		}
		item = "rt-configures";
		for (i = 0; i < root[item]["vrs-list"].size(); i++) {
			string svr = root[item]["vrs-list"][i]["server"].asCString();
			list<VrsStaItem> list_items;
			bamboo::excludeAnnoValue(value, root[item]["vrs-list"][i]["type"].asCString());
			int type = VrsStaItem::stream_type::stream;
			if (strstr(value, "vrs"))  
				type = VrsStaItem::stream_type::vrs;
			for (ivrs = 0; ivrs < root[item]["vrs-list"][i]["sta-list"].size(); ivrs++) {
				VrsStaItem vrsitem;
				bamboo::excludeAnnoValue(value, root[item]["vrs-list"][i]["sta-list"][ivrs]["name"].asCString());
				vrsitem.staname = value;
				vrsitem.type = type;
				if (vrsitem.type == VrsStaItem::stream_type::vrs) {
					bamboo::excludeAnnoValue(value, root[item]["vrs-list"][i]["sta-list"][ivrs]["xyz"].asCString());
					sscanf(value, "%lf%lf%lf", vrsitem.x, vrsitem.x + 1, vrsitem.x + 2);
				}
				if (!m_checkStation(vrsitem.staname)) {
					cout << "Existing same station,will continue for " << vrsitem.staname << endl;
					continue;
				}
				memset(cmd, 0, sizeof(cmd));
				if (vrsitem.type == VrsStaItem::stream_type::stream) {
					if ('/' == svr.back())
						sprintf(cmd, "%s%s", svr.c_str(), vrsitem.staname.c_str());
					else
						sprintf(cmd, "%s/%s", svr.c_str(), vrsitem.staname.c_str());
				}
				else {
					sprintf(cmd, "%s", svr.c_str());
				}
				vrsitem.strpath = string(cmd);
				list_items.push_back(vrsitem);
			}
			rt_mounts[svr] = list_items;
		}
		bamboo::excludeAnnoValue(value, root[item]["nrtk-port"].asCString());
		sscanf(value, "%d", &nrtkport);

		bamboo::excludeAnnoValue(value, root[item]["rtk-port"].asCString());
		sscanf(value, "%d", &rtkport);

		bamboo::excludeAnnoValue(value, root[item]["eph-port"].asCString());
		sscanf(value, "%d", &ephport);
		
		result = stat(f_jsonConfigures, &st);
		lastAct = st.st_mtime;
	}
	catch (...) {
		printf("file = %s,failed to parse database item", f_jsonConfigures);
		if (lexit) exit(1); else return;
	}
}
bool Deploy::m_checkStation(string sta) {
	list<string>::iterator strItr;
	list<VrsStaItem>::iterator itemItr;
	map<string, list<VrsStaItem>>::iterator mapItr;
	/// should check more precisely
	transform(sta.begin(), sta.end(), sta.begin(), ::toupper);
	for (strItr = post_stalist.begin(); strItr != post_stalist.end(); strItr++) {
		string value = (*strItr);
		transform(value.begin(), value.end(), value.begin(), ::toupper);
		if (!strncmp(sta.c_str(), value.c_str(), 4)) return false;
	}
	for (mapItr = rt_mounts.begin(); mapItr != rt_mounts.end(); mapItr++) {
		for (itemItr = (*mapItr).second.begin(); itemItr != (*mapItr).second.end(); ++itemItr) {
			string value = (*itemItr).staname;
			transform(value.begin(), value.end(), value.begin(), ::toupper);
			if (!strncmp(sta.c_str(), value.c_str(), 4)) return false;
		}
	}
	return true;
}
