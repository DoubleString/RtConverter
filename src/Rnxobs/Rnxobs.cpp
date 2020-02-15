/*
 * Rnxobs.cpp
 *
 *  Created on: 2018年5月5日
 *      Author: doublestring
 */
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdlib>

#include "../../include/Rnxobs/Rnxobs.h"
#include "../../include/Rnxobs/Com.h"
#include "../../include/RtConverter/Deploy.h"
using namespace std;
using namespace bamboo;
// Done at 8.26
void Rnxobs::v_openOutFile(string obscmd,double ver){

}
void Rnxobs::v_outFileHeader(){

}
/// exclude the glonass & galileo
void Rnxobs::v_outFileObs(int otmjd,double otsod){

}
void Rnxobs::v_closeOutFile(){

}
void Rnxobs::m_checkObs(int psat,double* obs,double *dop,double* snr){
}
