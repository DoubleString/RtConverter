/*
 * PatternName.cpp
 *
 *  Created on: 2018年2月6日
 *      Author: doublestring
 */


#include <iostream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include "../../include/Rnxobs/Patterns.h"
#include "../../include/Rnxobs/Com.h"
#include "../../include/RtConverter/Const.h"
using namespace std;
using namespace bamboo;
Patterns* Patterns::m_Instance = NULL;
Patterns::Patterns(){
	memset(this->fform,0,sizeof(this->fform));
	memset(this->fname,0,sizeof(this->fname));
	this->MAXFILE = 100;
	Patterns::readFile();
}
Patterns::~Patterns(){
}
void Patterns::readFile(){
	nfile = 0;
	strcpy(fname[nfile], "rnxo");
	strcpy(fform[nfile++], "-STANAM--DDD-0.-YY-o");
	strcpy(fname[nfile], "rnxn");
	strcpy(fform[nfile++], "-brd0--DDD-0.-YY-n");
	strcpy(fname[nfile], "rnxg");
	strcpy(fform[nfile++], "-brd0--DDD-0.-YY-g");
	strcpy(fname[nfile], "rnxl");
	strcpy(fform[nfile++], "-brd0--DDD-0.-YY-l");
	strcpy(fname[nfile], "rnxb");
	strcpy(fform[nfile++], "-brd0--DDD-0.-YY-b");
	strcpy(fname[nfile], "rnxc");
	strcpy(fform[nfile++], "-brd0--DDD-0.-YY-c");
	strcpy(fname[nfile], "rnxq");
	strcpy(fform[nfile++], "-brd0--DDD-0.-YY-q");
	strcpy(fname[nfile], "rnxp");
	strcpy(fform[nfile++], "-brd0--DDD-0.-YY-p");
}
void Patterns::m_getPatternName(int ldefined, const char* keyword, const char* param_list, int iyear,
		int imonth, int iday, int ihour, char* name) {
	char varword[1024];
	char form[1024] = { '\0' };
	int npar, nword, mjd, week, wd, idoy;
	bool lfound;
	char word[40][256], param[40][256];
	int i, j, k, index;
	if (ldefined)
		strcpy(form, name);
	else {
		for (i = 0; i < nfile; i++) {
			if (!strcmp(fname[i], keyword)) {
				strcpy(form, fform[i]);
				break;
			}
		}
		if (len_trim(form) == 0) {
			printf("keyword = %s,keyword not defined in panda_file_name",keyword);
			exit(1);
		}
	}
	split_string(false, form, ' ', ' ', '-', &nword, (char*)word, 256);
	strcpy(varword, param_list);
	split_string(true, varword, ' ', ' ', ':', &npar, (char*)param, 256);
	mjd = md_julday(iyear, imonth, iday);
	gpsweek(iyear, imonth, iday, &week, &wd);
	mjd2doy(mjd, &iyear, &idoy);
	sprintf(param[npar], "Y=%1.1d", iyear % 10);
	sprintf(param[npar + 1], "YY=%2.2d", iyear % 100);
	sprintf(param[npar + 2], "YYYY=%4.4d", iyear);
	sprintf(param[npar + 3], "DDD=%3.3d", idoy);
	sprintf(param[npar + 4], "HH=%2.2d", ihour);
	sprintf(param[npar + 5], "WWW=%4.4d", week);
	sprintf(param[npar + 6], "GPSW=%4.4d", week);
	sprintf(param[npar + 7], "WWWWD=%5.5d", week * 10 + wd);
	sprintf(param[npar + 8], "MJD=%5.5d", mjd);
	npar += 9;
	memset(name,0,strlen(name) * sizeof(char));
	for (i = 0; i < nword; i++) {
		if (len_trim(word[i]) == 0)
			continue;
		trim(word[i]);
		if (i % 2 != 0) {
			lfound = false;
			for (j = 0; j < npar; j++) {
				//if (strstr(param[j], word[i]) != NULL) {
				if (!strncmp(param[j], word[i],strlen(word[i]))) {
					lfound = true;
					k = len_trim(name);
					index = index_string(param[j], '=');
					substringEx(varword, param[j], index + 1,
							strlen(param[j]) - index - 1);
					strcpy(name + k, varword);
					break;
				}
			}
			if (!lfound) {
				printf("word = %s,keyword = %s,word is not defined in keyword",word[i],keyword);
				exit(1);
			}
		} else {
			k = len_trim(name);
			strcpy(name + k, word[i]);
		}
	}
}
