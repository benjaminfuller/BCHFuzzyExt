

#include "bchsketch.h"
#include <math.h>
#include <vector>
#include <assert.h>

void initializeGF2K(long m){
	GF2X irrP;
	BuildSparseIrred(irrP, m); // fix irreducible poly of deg m over GF(2)
	GF2E::init(irrP); // fix field as GF(2^m)
}

void multiplyVecGF2(vec_GF2 & result, const vec_GF2 & a, const vec_GF2 & b){
	clear(result);
	vec_GF2 temp;
	vec_GF2 temp2=b;
	cout<<a<<" "<<b<<endl;
	for(long i=0;i<=a.length();i++){
		if(a.get(i)==1){
			if(temp.length()<result.length()){
				temp.SetLength(result.length());
			}
			if(temp2.length()<result.length()){
			temp2.SetLength(result.length());
			}
			assert(result.length()>=temp.length());
			assert(result.length()>=temp2.length());
			shift(temp,temp2,i);
			cout<<temp.length()<<endl;
			cout<<"Temp to add "<<temp<<endl;
			result+=temp;
		}
	}
	cout<<"Result "<<result<<endl;
}
void vecGF2fromLong(vec_GF2 & polynomial, long value){
	if(value==0){
		polynomial.SetLength(1);
		clear(polynomial);
		return;
	}
	long digits = floor(log2(value))+1;
	polynomial.SetLength(digits);
	clear(polynomial);
	for(long i=0;i<=digits;i++)
	{
		if (((value>>i) %2)==1){
			polynomial.put(i, 1);
		}
		else{
			polynomial.put(i, 0);
		}
	}
}

void findGeneratorPolynomial(vec_GF2 & result, vector<vec_GF2> powerToMinPoly, long d){
	cout<<powerToMinPoly[0]<<" "<<powerToMinPoly[1]<<endl;
	vector<bool> isFresh = {0};
	cout<<"Find generator poly"<<endl;
	long lengthResult=0;
	for(long i=0; i<powerToMinPoly.size();i++){
		cout<<powerToMinPoly[i]<<endl;
		isFresh[i]=1;
		for(long j=0;j<i;j++){
			if (powerToMinPoly[j]==powerToMinPoly[i])
			{
				isFresh[i]=0;
			}
		}
	}
	cout<<"Min poly to include"<<endl;
	for(long i=0;i<d;i++){
		if(isFresh[i]==1){
			cout<<powerToMinPoly[i]<<endl;
			lengthResult+=powerToMinPoly[i].length();

		}
	}
	result=powerToMinPoly[0];
	result.SetLength(lengthResult);
	vec_GF2 tempResult = result;
	for(long i=1;i<d;i++){
		if(isFresh[i]==1){
			multiplyVecGF2(tempResult, result, powerToMinPoly[i]);
			result = tempResult;
		}
	}
	cout<<"Size of generating poly "<<result<<endl;
	
}

void findMinimumPolynomial(vec_GF2 & minimumPolynomial, const GF2E & generatorPower, long m){
	GF2E result;
	for(long i=1; i< (2<<m);i++){
		vecGF2fromLong(minimumPolynomial,i);
		result = evaluatePolyAtElement(minimumPolynomial, generatorPower);
		if(IsZero(result)){
			return;
		}
	}
}
//This takes the BCH code and finds a primitive element
// It also creates a vector mapping i -> g^i for the found generator
// This vector should cover all of the field except for 0
// TODO return a second vector that is the inverse mapping
// Had issues with equality when using an unordered set
GF2E initializeGF2EforBCH(vector<vec_GF2> & powerToMinPoly, vec_GF2E & powToElement, long m ){
	powToElement.SetLength((2<<(m-1))-1);
	long found=0;
	GF2E xx;
	do
	{
		xx = NTL::random_GF2E();
		found=1;
		// element_vector[(2<<(m-1))-1] = power(xx,0);
		if(IsZero(xx)){
			found=0;
			continue;
		}
		for(int i=0; i< (2<<m-1)-1;i++){
			GF2E temp_power = power(xx,i);
			powToElement[i] = temp_power;
			if(i>=1 && IsOne(powToElement[i])){
				found=0;
				break;
			}
		}
	} while (found==0);	
	//Now we're going to compute minimum polynomials for each power
	vec_GF2 minimumPolynomial;
	for(int i=0; i< (2<<m-1)-1;i++){
		findMinimumPolynomial(minimumPolynomial, powToElement[i], m);
		powerToMinPoly.push_back(minimumPolynomial);
	}
	cout<<powerToMinPoly.size()<<endl;
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

void BCHSyndromeCompute(vec_GF2E & ss, const vec_GF2 & bio_vector, long d, const GF2E & generator)
{
	for(long i =0;i<d;i++){
		GF2E gi = power(generator,i);
		ss[i]=evaluatePolyAtElement(bio_vector, power(generator,i));
		
	}
}

void translateErrors(vec_GF2 & errors, vec_GF2E & located_errors, const vec_GF2E& powToElement){
	for(long i=0;i<located_errors.length();i++){
		for(long j=0;j<powToElement.length();j++){
			if(powToElement[j]==located_errors[i]){
				errors.put(j, 1);
				continue;
			}
		}
	}
	return;
}
///////////////////////////////////////////////////////////////////////////
// PURPOSE:
// Returns true if f fully factors into distinct roots
// (i.e., if f is a product of distinct monic degree-1 polynomials
// times possibly a constant)
// and false otherwise.
// If f is zero, returns false.
//
// ALGORITHM:
// Let m=GF2E::degree() (i.e., the field is GF(2^m)).
// The check is accomplished by checking if f divides X^{2^m} - X, 
// or equivalently if X^{2^m}-X is 0 mod f. 
// X^{2^m} - X has 2^m distinct roots -- namely,
// every element of the field is a root.  Hence, f divides it if and only
// if f has all its roots and they are all distinct.
//
// RUNNING TIME:
// Depends on NTL's implementation of FrobeniusMap, but for inputs of degree
// e that is relatively small compared m, should take e^{\log_2 3} m
// operations in GF(2^m).  Note that \log_2 3 is about 1.585.
static
bool CheckIfDistinctRoots(const GF2EX & f)
{
	if (IsZero(f))
	  return false;
	// We hanlde degree 0 and degree 1 case separately, so that later
	// we can assume X mod f is the same as X
	if (deg(f) == 0 || deg(f) == 1)
	  return true;

	GF2EXModulus F;
	// F is the same as f, just more efficient modular operations
	build(F, f);  	

	GF2EX h;
	FrobeniusMap (h, F); // h = X^{2^m} mod F

	// If X^{2^m} - X = 0 mod F, then X^{2^m} mod F
	// should be just X (because degree of F > 1)
	return (IsX(h));
}

bool BCHSyndromeDecode(vec_GF2 & translated_errors, const vec_GF2E & syndrome, long d, const vec_GF2E & powToElement )
{
        long i;	

	GF2EX r1, r2, r3, v1, v2, v3,  q, temp;
	GF2EX *Rold, *Rcur, *Rnew, *Vold, *Vcur, *Vnew, *tempPointer;

        // Use pointers to avoid moving polynomials around
        // An assignment of polynomials requires copying the coefficient vector;
        // we will not assign polynomials, but will swap pointers instead
	Rold = &r1;
	Rcur = &r2;
	Rnew = &r3;

	Vold = &v1;
	Vcur = &v2;
	Vnew = &v3;

	vec_GF2E answer;

	SetCoeff(*Rold, d-1, 1); // Rold holds z^{d-1}

	// Rcur=S(z)/z where S is the syndrome poly, Rcur = \sum S_j z^{j-1}
        // Note that because we index arrays from 0, S_j is stored in ss[j-1]
	for (i=0; i<d-1; i++){ 
	  SetCoeff (*Rcur, i, syndrome[i]);
	}
	// Vold is already 0 -- no need to initialize
	// Initialize Vcur to 1
	SetCoeff(*Vcur, 0, 1); // Vcur = 1

	// Now run Euclid, but stop as soon as degree of Rcur drops below
	// (d-1)/2
	// This will take O(d^2) operations in GF(2^m)
	long t = (d-1)/2;

	while (deg(*Rcur) >=  t) {
	  // Rold = Rcur*q + Rnew
	  DivRem(q, *Rnew, *Rold, *Rcur);

	  // Vnew = Vold - qVcur)
	  mul(temp, q, *Vcur);
	  sub (*Vnew, *Vold, temp);


          // swap everything
	  tempPointer = Rold;	
	  Rold = Rcur;
	  Rcur = Rnew;
	  Rnew = tempPointer;

	  tempPointer = Vold;
	  Vold = Vcur;
	  Vcur = Vnew;
	  Vnew = tempPointer;
	}

	if (IsZero(ConstTerm(*Vcur)))
	  return false;

	MakeMonic(*Vcur);
	if (CheckIfDistinctRoots(*Vcur) == false){
	  return false;
	}
    
	answer = FindRoots(*Vcur);
	for (i = 0; i < answer.length(); ++i){
		answer[i] = inv(answer[i]);
	}
	
	translated_errors.SetLength(powToElement.length());
	translateErrors(translated_errors, answer, powToElement);
	vec_GF2E test;
	test.SetLength(d);
	BCHSyndromeCompute (test, translated_errors, d, powToElement[1]);
	return (test==syndrome);
	
}

