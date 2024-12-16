
#include "bchsketch.h"
#include <cassert>
#include "sketch.h"

void syndromeTest(long d, vec_GF2 &original_bio, GF2E &generator, vec_GF2E &genToPow);
void codeOffsetTest(long d, vec_GF2 &original_bio, GF2E &generator, vec_GF2E &genToPow, vector<vec_GF2> & powToMinimum);

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
	vector<vec_GF2> powToMinimumPolynomial;
	generator = initializeGF2EforBCH(powToMinimumPolynomial, genToPow,m);
	assert(!IsZero(generator));
    syndromeTest(d, original_bio, generator, genToPow);
    codeOffsetTest(d, original_bio, generator, genToPow, powToMinimumPolynomial);

	return 0;
}

void syndromeTest(long d, vec_GF2 &original_bio, GF2E &generator, vec_GF2E &genToPow)
{
    vec_GF2E syndrome;
    syndrome.SetLength(d);
    vec_GF2E syndrome2;
    syndrome2.SetLength(d);
    vec_GF2 error_vect;
    error_vect.SetLength(original_bio.length());
    error_vect.put(1, 1);
    vec_GF2 noisy_bio = original_bio + error_vect;
    // cout<<vect_i_test<<endl;
    BCHSyndromeCompute(syndrome, original_bio, d, generator);
    BCHSyndromeCompute(syndrome2, noisy_bio, d, generator);
    vec_GF2E final_syndrome;
    final_syndrome.SetLength(d);
    final_syndrome = syndrome + syndrome2;
    vec_GF2 decoded_value;
    BCHSyndromeDecode(decoded_value, final_syndrome, d, genToPow);
    vec_GF2 corrected_value = decoded_value + noisy_bio;
    assert(corrected_value == original_bio);
	cout<<"Syndrome test passed"<<endl;
}


void codeOffsetTest(long d, vec_GF2 &original_bio, GF2E &generator, vec_GF2E &genToPow, vector<vec_GF2> & powToMinimum)
{
	vec_GF2 genPoly;
	findGeneratorPolynomial(genPoly, powToMinimum, d);
}