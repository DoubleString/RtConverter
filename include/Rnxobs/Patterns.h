/*
 * PatternName.h
 *
 *  Created on: 2018/2/6
 *      Author: doublestring
 */

#ifndef INCLUDE_COM_PATTERNS_H_
#define INCLUDE_COM_PATTERNS_H_
#include <iostream>
#include <cstdlib>
using namespace std;
namespace bamboo{
class Patterns{
public:
	Patterns();
	~Patterns();
	void m_getPatternName(int ldefined, const char* keyword, const char* param_list, int iyear,
			int imonth, int iday, int ihour, char* name);
	static inline Patterns* s_getInstance(){
		if(Patterns::m_Instance == NULL){
			Patterns::m_Instance = new Patterns();
		}
		return Patterns::m_Instance;
	}
	static inline void s_destroyInstance(){
		if(m_Instance != NULL)
			delete Patterns::m_Instance;
		Patterns::m_Instance = NULL;
	}
public:
	void readFile();
	static Patterns* m_Instance;
	char fname[100][256];
	char fform[100][256];
	int nfile;
	int MAXFILE;
};
}
#endif /* INCLUDE_COM_PATTERNS_H_ */
