
#include <bits/stdc++.h>
#include "bchsketch.h"
using namespace std;

/***************** Converting integer to elements of GF(2^e) and back **************************/


// Uses the bits of the polynomial representation of an element of GF2E to come up with
// an integer with the same bits
static 
void BinElemToNum(ZZ & a, const GF2E & e) {
	GF2X g = rep(e); // convert to polynomial
	long numBytes=NumBytes(g);
	unsigned char * buffer = new unsigned char[numBytes];
	BytesFromGF2X(buffer, g, numBytes);
	ZZFromBytes(a, buffer, numBytes);
	delete [] buffer;
}


// Uses binary representation of an integer to come up with an
// element of GF(2^m) by using the bits of the integer
// as bits of a polynomial and then reducing modulo the irreducible polynomial that
// generates the field.
static 
void NumToBinElem(GF2E & e, const ZZ & a, unsigned long m) {
    GF2X g;
	long numBytes=NumBytes(a);
	unsigned char * buffer = new unsigned char[numBytes];
	BytesFromZZ(buffer, a, numBytes);
	GF2XFromBytes(g, buffer, numBytes);
	delete [] buffer;
    conv (e, g); // Convert from polynomial to field element
}



/******************************************** I/O ROUTINES ******************************/

// Reads in a vector of integers from a file, converting them to elements of GF(2^m)
// Designed to read a .set file after the ReadSetParams routine 
void ReadBioInput(vec_GF2 & vect_i, istream &infile, unsigned long & m)
{
	char c1='\0';
	unsigned int r =0;
	while(c1 !='['){
		infile>>c1;
	}
	unsigned long length =0;
  	m=0;
	vector<long> bin_vec ;
	while(!infile.eof()){
		infile>>r;
		infile>>c1;
		bin_vec.push_back(r);
		length++;
	}
	unsigned long temp_n=0;
	while(temp_n<length){
		m++;
		temp_n = (2<<m-1) -1;
	}
	
	long n = temp_n;
	while(bin_vec.size()<n){
		bin_vec.push_back(0);
	}

	vect_i.SetLength(n);
	for(int i=0; i<n;i++){
		GF2 temp;
		temp = bin_vec[i];
		vect_i.put(i,temp);
	}
}



// Reads in the field size m and the desired error-tolerance t from
// the input file, where they are assumed to be present in the format
// t=<integer>(no spaces around '=' are allowed)
// Returns  the minimum distance of the code d=2t+1
void ReadDFile(long &d, istream &infile)
{
	long t = 0;  // t = max set diff tolerated
	int count = 1; // # of params to be read in
	char c1, c2;
	ZZ temp;
	
	while (count > 0)
	{
		c1 = '\0';
	  	c2 = '\0';
		while (c2 != '=' && !infile.eof())
			infile >> c1 >> c2;
		count--;
		switch (c1)
		{
			case 't': infile >> t; break;
		}
	}
	if (infile.eof() || t<=0)
	{
		cerr << "Bad input format!" << endl;
		exit(-1);
	}

	d = 2*t+1; // t = max set diff tolerated, d=2*t+1

	return;
}

