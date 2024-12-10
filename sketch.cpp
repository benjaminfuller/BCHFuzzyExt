// sketch.cpp
//
// C++ code by Kevin Harmon and Leonid Reyzin
//
// Finds the PinSketch (BCH-based secure sketch) of an input set.
//
// See pinsketch.txt for detailed documentation
// This code and explanatory notes
// are hosted at http://www.cs.bu.edu/~reyzin/code/fuzzy.html

//

#include "pinsketch.h"
#include <cassert>

void initializeGF2K(long m){
	GF2X irrP;
	BuildSparseIrred(irrP, m); // fix irreducible poly of deg m over GF(2)
	GF2E::init(irrP); // fix field as GF(2^m)
}

GF2E findPrimitiveElement(vec_GF2E & element_vector, long m ){
	element_vector.SetLength((2<<(m-1))-1);
	long found=0;
	GF2E xx;
	do
	{
		xx = NTL::random_GF2E();
		found=1;
		element_vector[(2<<(m-1))-2] = power(xx,0);
		if(IsZero(xx)){
			continue;
		}
		for(int i=0; i< (2<<m-1)-2;i++){
			element_vector[i] = power(xx,i);
			if(i>=1 && IsOne(element_vector[i])){
				found=0;
				break;
			}
		}
		cout<<"Value of found "<<found<<endl;
		cout<<"Element"<<element_vector<<endl;
		
	} while (found==0);	
	return xx;
}

GF2E evaluatePolyAtElement(const vec_GF2 & polynomial, const GF2E & point){
	GF2E running_sum = GF2E(0);
	for(long i=0;i<polynomial.length();i++){
		if(IsOne(polynomial[i])){
			running_sum= running_sum+power(point,i);
		}
	}
	return running_sum;
}

int main(int argc, char** argv)
{
	long d; // minimum distance of the code;
	// can handle set difference up to t=(d-1)/2 elements
	// sketch is (d-1)/2 elements long

	int len;// length of argv[1] if it exists


	if (argc != 2 || (len=strlen(argv[1]))<5 || strcmp(&argv[1][len-4], ".set"))
	{
		cerr << "Usage: sketch A.set" << endl;
		if (argc == 2)
			cerr << "(file must be named `*.set`)" << endl;
		return -1;
	}

// Fix field and error-tolerance
	ifstream infile(argv[1]);
        if (!infile) {
          cerr << "Could not open file for reading!\n";
          return -1;
        }
	ReadDFile(d, infile);
	unsigned int long m =0;
	vec_GF2 vect_i;
	ReadBioInput(vect_i, infile, m);
	infile.close();
	initializeGF2K(m);
	
	GF2E generator;
	vec_GF2E genToPow;
	generator = findPrimitiveElement(genToPow,m);
	assert(!IsZero(generator));
	cout<<"Gen to Pow"<<genToPow<<endl;
	cout<<"Generator"<<generator<<endl;
	for(long i =0;i<d;i++){
		GF2E gi = power(generator,i);
		GF2E syn_i =evaluatePolyAtElement(vect_i, power(generator,i));
		cout<<"Syndrome "<<i<<" "<<syn_i<<endl;
	}
	
	vec_GF2E syndrome;

// // read in set  
// 	vec_GF2 vector;
// 	// vec_GF2E ss;
// 	ReadBioInput(vector, infile, m);

// // compute secure sketch of the set
// 	BCHSyndromeCompute(ss, vector, d);

// // write the sketch to file with same name and .ss extension
// 	strcpy(&argv[1][strlen(argv[1])-4], ".ss");
// 	ofstream outfile (argv[1], ios::out | ios::trunc);
//         if (!outfile) {
//           cerr << "Could not open file for writing!\n";
//           return -1;
//         }
// 	OutputSS(outfile, ss, d, GF2E::modulus());
// 	outfile.close();

// 	cout << "Secure sketch written to file `" << argv[1] << "`." << endl;
	
	return 0;
}
