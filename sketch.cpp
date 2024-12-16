
#include "bchsketch.h"
#include <cassert>


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
	vec_GF2 original_bio;
	ReadBioInput(original_bio, infile, m);
	infile.close();
	initializeGF2K(m);
	GF2E generator;
	vec_GF2E genToPow;
	generator = initializeGF2EforBCH(genToPow,m);
	assert(!IsZero(generator));
	vec_GF2E syndrome;
	syndrome.SetLength(d);
	vec_GF2E syndrome2;
	syndrome2.SetLength(d);
	vec_GF2 error_vect;
	error_vect.SetLength(original_bio.length());
	error_vect.put(1,1);
	vec_GF2 noisy_bio =original_bio + error_vect;
	// cout<<vect_i_test<<endl;
	BCHSyndromeCompute(syndrome,original_bio,d, generator);
	BCHSyndromeCompute(syndrome2,noisy_bio,d, generator);
	vec_GF2E final_syndrome;
	final_syndrome.SetLength(d);
	final_syndrome = syndrome+syndrome2;
	vec_GF2 decoded_value;
	BCHSyndromeDecode(decoded_value,final_syndrome,d, genToPow);
	vec_GF2 corrected_value = decoded_value+ noisy_bio;
	assert(corrected_value== original_bio);

	cout<<"Test passed successfully"<<endl;
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
