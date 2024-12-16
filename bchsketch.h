#ifndef __BCH_SKETCH_H
#define __BCH_SKETCH_H

#include <fstream>
#include <unordered_map>
#include <NTL/vec_ZZ.h>
#include <NTL/GF2XFactoring.h>
#include <NTL/GF2EXFactoring.h>

using namespace std; // may be needed to compile on some platforms; may need to be removed on others

NTL_CLIENT

/************************ I/O ************************************/
void ReadDFile(long &d, istream &infile);
void ReadBioInput(vec_GF2 & vect_i, istream &infile, unsigned long & m);

/************************ BCH Syndrome Encoding/Decoding ***************/
void BCHSyndromeCompute(vec_GF2E &answer, const vec_GF2 & bio_vec, long d, const GF2E & generator);
bool BCHSyndromeDecode(vec_GF2 & translated_errors, const vec_GF2E & syndrome, long d, const vec_GF2E & powToElement );

/************************ Helper Functions for Polynomial ***************/
void initializeGF2K(long m);
GF2E initializeGF2EforBCH(vec_GF2E & powToElement, long m );
GF2E evaluatePolyAtElement(const vec_GF2 & polynomial, const GF2E & point);
void translateErrors(vec_GF2 & errors, vec_GF2E & located_errors, const vec_GF2E& powToElement);

#endif
