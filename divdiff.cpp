//
// divdiff version 1.0. These routines are introduced in the paper:
// 
// This program is licensed under a Creative Commons Attribution 4.0 International License:
// http://creativecommons.org/licenses/by/4.0/
//

#include<random>
#include<stdio.h>
#include<string.h>
#include<math.h>

double* invPowersOf2;

class ld{
private:
	double mantissa;
	int exponent;
public:
	ld(){	mantissa = 0.5; exponent = 1;}
	void normalize(){ int tmp; mantissa = frexp(mantissa,&tmp); exponent += tmp;}
	void init_expmu(double mu){ double e = mu*1.4426950408889634; exponent = ceil(e); mantissa = pow(2.,e - ceil(e)); }
	void print(){
		double exp10, m;
		exp10 = (exponent*0.30102999566398114);
		m = mantissa*pow(10,exp10 - floor(exp10)); 
		exp10 = floor(exp10); if(m<1){ exp10--; m*=10; }
		if((exp10<7)&&(exp10>-7)) printf("%.17f",m*pow(10,exp10)); else printf("%.17f*^%.0f",m,exp10);
	}
	ld operator =(ld const &obj){ mantissa = obj.mantissa; exponent = obj.exponent;	return *this;}
	ld operator =(double const &obj){ mantissa = obj; exponent = 0; normalize(); return *this;}
	ld operator +(ld const &obj){ // important restriction: it is assumed here that each of the summands is not equal to zero
		ld res;
		if(obj.exponent >= exponent){
			res.mantissa = obj.mantissa + mantissa*invPowersOf2[obj.exponent - exponent];
			res.exponent = obj.exponent; res.normalize();
		} else{
			res.mantissa = mantissa + obj.mantissa*invPowersOf2[exponent - obj.exponent];
			res.exponent = exponent; res.normalize();
		}
		return res;
	}
	ld operator -(ld const &obj){ // important restriction: it is assumed here that each of the summands is not equal to zero
		ld res;
		if(obj.exponent >= exponent){
			res.mantissa = mantissa*invPowersOf2[obj.exponent - exponent] - obj.mantissa;
			res.exponent = obj.exponent; res.normalize();
		} else{
			res.mantissa = mantissa - obj.mantissa*invPowersOf2[exponent - obj.exponent];
			res.exponent = exponent; res.normalize();
		}
		return res;
	}
	ld operator +=(ld const &obj){ // important restriction: it is assumed here that each of the summands is not equal to zero
		if(obj.exponent >= exponent){
			mantissa = obj.mantissa + mantissa*invPowersOf2[obj.exponent - exponent];
			exponent = obj.exponent; normalize();
		} else{
			mantissa = mantissa + obj.mantissa*invPowersOf2[exponent - obj.exponent];
			exponent = exponent; normalize();
		}
		return *this;
	}
	ld operator -=(ld const &obj){ // important restriction: it is assumed here that each of the summands is not equal to zero
		if(obj.exponent >= exponent){
			mantissa = mantissa*invPowersOf2[obj.exponent - exponent] - obj.mantissa;
			exponent = obj.exponent; normalize();
		} else{
			mantissa = mantissa - obj.mantissa*invPowersOf2[exponent - obj.exponent];
			exponent = exponent; normalize();
		}
		return *this;
	}
	ld operator *(ld const &obj){
		ld res;	res.mantissa = mantissa * obj.mantissa;	res.exponent = exponent + obj.exponent;
		res.normalize(); return res;
	}
	ld operator /(ld const &obj){
		ld res;	res.mantissa = mantissa / obj.mantissa;
		res.exponent = exponent - obj.exponent;	res.normalize(); return res;
	}
	ld operator *(double const &obj){ ld res; res.mantissa = mantissa * obj; res.exponent = exponent; res.normalize(); return res; }
	ld operator /(double const &obj){ ld res; res.mantissa = mantissa / obj; res.exponent = exponent; res.normalize(); return res; }
	ld operator *=(ld const &obj){ mantissa *= obj.mantissa; exponent += obj.exponent; normalize(); return *this;}
	ld operator /=(ld const &obj){ mantissa /= obj.mantissa; exponent -= obj.exponent; normalize(); return *this;}
	ld operator *=(double const &obj){ mantissa *= obj; normalize(); return *this;}
	ld operator /=(double const &obj){ mantissa /= obj; normalize(); return *this;}
	int operator >=(double const &r){ // important restriction: it is assumed here that both values of mantissa are not negative
		if(r == 0) return (mantissa >= 0);
		else{
			ld R; R = r;
			if(exponent > R.exponent ) return 1;
			else if((exponent == R.exponent)&&(mantissa >= R.mantissa)) return 1;
			else return 0;
		}
	}
	double get_double(){ return mantissa * pow(2.,exponent); }
};

int extralen = 30;

int maxexp = 100000;

double *z, *z2; ld *h, *ddd, *divdiff; 
int CurrentLength; int s; double mu; ld expmu;
int maxlen = 10001; int smax = 500;

void PrintList(ld* list, int len, const char* namelist){
	int i,flag=1;
	printf("%s={",namelist);
	for(i=0;i<len;i++){
		list[i].print();
		if(i<len-1) printf(","); else printf("};\n");
	}
}

void PrintList(double* list, int len, const char* namelist){
	int i,flag=1;
	printf("%s={",namelist);
	for(i=0;i<len;i++){
		printf("%.17g",list[i]);
		if(i<len-1) printf(","); else printf("};\n");
	}
}

void Backupz(int len){ z2 = new double[maxlen]; memcpy(z2,z,maxlen*sizeof(double));}
void Restorez(int len){ memcpy(z,z2,maxlen*sizeof(double)); delete[] z2;}

void AllocateMem(){
	z = new double[maxlen]; h = new ld[maxlen+extralen+1];
	divdiff = new ld[maxlen+extralen+1]; ddd = new ld[maxlen*smax];
	CurrentLength = 0; s = 1;
}

void FreeMem(){
	delete[] z; delete[] h; delete[] divdiff; delete[] ddd;
}

void init(){
	invPowersOf2 = new double[maxexp];
	double curr=1; for(int i=0;i<maxexp;i++){ invPowersOf2[i] = curr; curr/=2; }
	AllocateMem();
}

void clear_up(){ delete[] invPowersOf2; FreeMem();}

template <typename T> T min (T a, T b) { return (a<b?a:b);}
template <typename T> T max (T a, T b) { return (a<b?b:a);}

double mean(double* z, int n){
	double sum=0; int i;
	for(i=0;i<n;i++) sum+=z[i];
	return sum/n;
}

double maxAbsDiff(double* z, int len){
	double zmax = z[0], zmin = z[0]; int i;
	for(i=1;i<len;i++) { zmin = min(zmin,z[i]); zmax = max(zmax,z[i]);}
	return fabs(zmax-zmin);
}

int s_changed(double znew = z[CurrentLength-1]){
	return fabs(znew-mu)/3.5 > s;
}

void AddAll(int len, int force_s);

void AddElement(double znew, int force_s = 0, double force_mu = 0){
	int j,k,n,l,N; ld curr; n = CurrentLength; l = n + 1; N = maxlen+extralen; z[n] = znew; CurrentLength++;
	if(CurrentLength==1){
		s = (force_s == 0) ? 1 : force_s;
		mu = (force_mu == 0) ? z[0] : force_mu;
		expmu.init_expmu(mu);
		h[0] = 1; for(k=1;k<=N;k++) h[k] = h[k-1]/s;
		if(mu != z[0]) for(k=N;k>0;k--) h[k-1] += h[k]*(z[0]-mu)/k;
		curr = expmu*h[0]; for(k=0;k<s-1;k++) { ddd[k*maxlen] = curr; curr*=h[0];}
		divdiff[0].init_expmu(z[0]); // alternatively: divdiff[0] = curr;
	} else if(s_changed()||(CurrentLength>=maxlen)) AddAll(CurrentLength, force_s);
	else{
		for(k=N;k>n;k--) h[k-1] += h[k]*(z[n]-mu)/k; curr = expmu*h[n];
		for(k=n;k>=1;k--) h[k-1] = (h[k-1]*n + h[k]*(z[n]-z[n-k]))/(n-k+1);
		for(k=0;k<s-1;k++){
			ddd[k*maxlen+n] = curr;
			curr = ddd[k*maxlen]*h[n]; for(j=1;j<=n;j++) curr += ddd[k*maxlen+j]*h[n-j];
		}
		divdiff[n] = curr;
	}
}

void RemoveElement(){ 
	int k,n,N;
	if(CurrentLength>=1){ 
		n = CurrentLength - 1; N = maxlen+extralen;
		for(k=1;k<=n;k++) h[k-1] = (h[k-1]*(n-k+1) - h[k]*(z[n]-z[n-k]))/n;
		for(k=n+1;k<=N;k++) h[k-1] -= h[k]*(z[n]-mu)/k;
		CurrentLength--;
	}
}

void AddAll(int len, int force_s = 0){ // input is taken from z, output is written to divdiff, size of z should be not smaller than len
	int i,s; CurrentLength = 0;
	if(force_s == 0) s = (int)ceil(maxAbsDiff(z,len)/3.5); else s = force_s;
	if((s > smax)||(len >= maxlen)){
		i = maxlen; Backupz(i); FreeMem();
		if(s > smax) smax = max(smax*2,s);
		if(len >= maxlen) maxlen = max(maxlen*2,len);
		AllocateMem(); Restorez(i);
	}
	AddElement(z[0],s,mean(z,len));
	for(i=1;i<len;i++) AddElement(z[i]);
	// calculates the vector (d[z_0], 1! d[z_0,z_1], ..., n! d[z_0,z_1,...,z_n]).
}

int main(){
	double newz, sigma=40; int nn=10000;
	init();	std::mt19937 gen(100); std::normal_distribution<float> distr(0, sigma); 
	for(int i=0;i<nn;i++){
		newz = distr(gen); 
		AddElement(newz); 
	}
	printf("Length=%d\n",CurrentLength);
	PrintList(z,CurrentLength,"z");
        PrintList(divdiff,CurrentLength,"DivDiffsRel"); // the vector (d[z_0], 1! d[z_0,z_1], ..., n! d[z_0,z_1,...,z_n])
	clear_up();
	return 0;
}

