#include <mutex>								// g++47 needs it! clang++4.0/g++44 doesn't.
#include <thread>
#include "libfbj_fet.hpp"
#include "libfbj_math.hpp"
#include "libfbj_base.hpp"

bool ROC_VERBOSE=false;

#define ChkNaN1(x)     if (std::isnan(x))                               return std::numeric_limits<double>::signaling_NaN()
#define ChkNaN2(x,y)   if (std::isnan(x)||std::isnan(y))                return std::numeric_limits<double>::signaling_NaN()
#define ChkNaN3(x,y,z) if (std::isnan(x)||std::isnan(y)||std::isnan(z)) return std::numeric_limits<double>::signaling_NaN()

// http://www.johndcook.com/cpp_phi.html (public domain)
double cdf_norms_2sided_pv(double x)
{
	ChkNaN1(x);
	
	// constants
	static const double a1 =  0.254829592;
	static const double a2 = -0.284496736;
	static const double a3 =  1.421413741;
	static const double a4 = -1.453152027;
	static const double a5 =  1.061405429;
	static const double p  =  0.3275911;
	
	// A&S formula 7.1.26
	x = fabs(x)/sqrt(2.0);
	double t = 1.0/(1.0 + p*x);
	return (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
}

// chisq nu is double, because nu doesn't have to be an integer, which makes a diff in results
// function qxxx = quantile of xxx distribution; cdf_xxx_p is inv of qxxx.
#if defined(_STAT_USE_GSL)
	#include <gsl/gsl_cdf.h> // gsl P=left_tail Q=right_tail
// students t, nu = d.f = sample size - 1; 1st func same as boost::math::cdf(boost::math::complement(dist,std::fabs(t)))*2;
	double cdf_t_p(double t,double nu)	{	if (isinf(t)) return t>0?1:0; return gsl_cdf_tdist_P(t,nu); }
	double cdf_t_q(double t,double nu)	{	if (isinf(t)) return t>0?0:1; return gsl_cdf_tdist_Q(t,nu); }
	double cdf_t_2sided_pv(double t,double nu)	{	if (isinf(t)) return 0; return gsl_cdf_tdist_Q(std::fabs(t),nu)*2; }
	double qstudents_t(double p, double nu)		{							return gsl_cdf_tdist_Pinv(p,nu); }
// chi-square distribution
	double qchisq_1df(double p)					{								   return gsl_cdf_chisq_Pinv(p,1); }
	double qchisq(double p,double nu)			{								   return gsl_cdf_chisq_Pinv(p,nu); }
	double cdf_chisq_1df_q(double x)			{	if (x>0 && isinf(x)) return 0; return x<=0?1:gsl_cdf_chisq_Q(x,1);	}
	double cdf_chisq_1df_p(double x)			{	if (x>0 && isinf(x)) return 1; return x<=0?0:gsl_cdf_chisq_P(x,1);	}
	double cdf_chisq_q(double x,double nu)		{	if (x>0 && isinf(x)) return 0; return x<=0?1:gsl_cdf_chisq_Q(x,nu); }
	double cdf_chisq_p(double x,double nu)		{	if (x>0 && isinf(x)) return 1; return x<=0?0:gsl_cdf_chisq_P(x,nu); }
// F distribution
	double cdf_f_q(double x, double nu1, double nu2)	{ return gsl_cdf_fdist_Q(x,nu1,nu2); }
	double cdf_f_p(double x, double nu1, double nu2)	{ return gsl_cdf_fdist_P(x,nu1,nu2); }
// standard normal distribution
	double pdf_norms(double z)					{	return gsl_sf_erf_Z(z); }
	double cdf_norms(double z)					{	return gsl_cdf_ugaussian_P(z); }
	double cdf_norms_q(double z)				{	return gsl_cdf_ugaussian_Q(z); }
	double cdf_norms_p(double z)				{	return gsl_cdf_ugaussian_P(z); }
	//double cdf_norms_2sided_pv(double z)		{	return gsl_cdf_ugaussian_Q(std::fabs(z))*2; }
	double qnorms(double p)						{	return gsl_cdf_ugaussian_Pinv(p);	}
// non-standard normal distribution
	double cdf_norm_p(double mu, double sd, double x) {	return gsl_cdf_gaussian_P(x-mu,sd); }
	double cdf_norm_q(double mu, double sd, double x) {	return gsl_cdf_gaussian_Q(x-mu,sd); }
	double cdf_norm_2sided_pv(double mu, double sd, double x){	return gsl_cdf_gaussian_Q(x-mu,sd)*2; }
	double qnorm(double mu, double sd, double p){	return gsl_cdf_gaussian_Pinv(p,sd)+mu; }
// binomial distribution
	double pdf_binomial	  (double k, double p, double n){ if (k<0||n<=0) exit_error("binomial k<0||n<=0"); return gsl_ran_binomial_pdf(my_round(k),p,my_round(n)); }
	double cdf_binomial_le(double k, double p, double n){ if (k<0||n<=0) exit_error("binomial k<0||n<=0"); return gsl_cdf_binomial_P  (my_round(k),p,my_round(n)); } // not tested whether k is included
	double cdf_binomial_ge(double k, double p, double n){ if (k<0||n<=0) exit_error("binomial k<0||n<=0"); return gsl_cdf_binomial_Q  (my_round(k),p,my_round(n)); } // not tested whether k is included
#elif defined(_STAT_USE_BOOST)
	#include <boost/math/distributions.hpp>
// students t, nu = d.f = sample size - 1; 1st func same as boost::math::cdf(boost::math::complement(dist,std::fabs(t)))*2;
	double cdf_t_q(double t,double nu)		{	ChkNaN2(t,nu); if (isinf(t)) return t>0?0:1; boost::math::students_t dist(nu); return 1-boost::math::cdf(dist,t); }
	double cdf_t_p(double t,double nu)		{	ChkNaN2(t,nu); if (isinf(t)) return t>0?1:0; boost::math::students_t dist(nu); return   boost::math::cdf(dist,t); }
	double cdf_t_2sided_pv(double t, double nu){ChkNaN2(t,nu); if (isinf(t)) return 0; 	boost::math::students_t dist(nu);	return(1-boost::math::cdf(dist,std::fabs(t)))*2;}
	double qstudents_t(double p, double nu)	{	ChkNaN2(p,nu);
		if (p<0 || p>1) exit_error("p cannot be <0 or >1 for qnorms(p)");
		if (p==0) return -INFINITY;
		if (p==1) return INFINITY;
		boost::math::students_t dist(nu); return boost::math::quantile(dist,p);}
// chi-square distribution
	double qchisq_1df(double p)				{	ChkNaN1(p);
		if (p<0 || p>1) exit_error("p cannot be <0 or >1 for qnorms(p)");
		if (p==0) return 0;
		if (p==1) return INFINITY;
		static const	boost::math::chi_squared dist(1);	return boost::math::quantile(dist,p); }
	double qchisq(double p,double nu)		{	ChkNaN2(p,nu);
		if (p<0 || p>1) exit_error("p cannot be <0 or >1 for qnorms(p)");
		if (p==0) return 0;
		if (p==1) return INFINITY;
		static const	boost::math::chi_squared dist(nu);	return boost::math::quantile(dist,p); }
	double cdf_chisq_1df_q(double x)		{	ChkNaN1(x);    if (x>0 && isinf(x)) return 0; static const	boost::math::chi_squared dist(1);	return 1-(x<=0?0:boost::math::cdf(dist,x)); }
	double cdf_chisq_1df_p(double x)		{	ChkNaN1(x);    if (x>0 && isinf(x)) return 1; static const	boost::math::chi_squared dist(1);	return   (x<=0?0:boost::math::cdf(dist,x)); }
	double cdf_chisq_q(double x,double nu)	{	ChkNaN2(x,nu); if (x>0 && isinf(x)) return 0; 				boost::math::chi_squared dist(nu);	return 1-(x<=0?0:boost::math::cdf(dist,x)); }
	double cdf_chisq_p(double x,double nu)	{	ChkNaN2(x,nu); if (x>0 && isinf(x)) return 1; 				boost::math::chi_squared dist(nu);	return   (x<=0?0:boost::math::cdf(dist,x)); }
// F distribution
	double cdf_f_q(double x, double nu1, double nu2)	{ ChkNaN3(x,nu1,nu2); boost::math::fisher_f dist(nu1,nu2);	return 1-boost::math::cdf(dist,x); }
	double cdf_f_p(double x, double nu1, double nu2)	{ ChkNaN3(x,nu1,nu2); boost::math::fisher_f dist(nu1,nu2);	return   boost::math::cdf(dist,x); }
// standard normal distribution
	double pdf_norms(double z)				{	ChkNaN1(z); static const	boost::math::normal	dist;			return   boost::math::pdf(dist,z); }
	double cdf_norms(double z)				{	ChkNaN1(z); static const	boost::math::normal	dist;			return	 boost::math::cdf(dist,z); }
	double cdf_norms_q(double z)			{	ChkNaN1(z); static const	boost::math::normal	dist;			return 1-boost::math::cdf(dist,z); }
	double cdf_norms_p(double z)			{	ChkNaN1(z); static const	boost::math::normal	dist;			return	 boost::math::cdf(dist,z); }
	//double cdf_norms_2sided_pv(double z)	{	static const	boost::math::normal	dist;			return(1-boost::math::cdf(dist,std::fabs(z)))*2;}
	double qnorms(double p)					{	ChkNaN1(p);
		if (p<0 || p>1) exit_error("p cannot be <0 or >1 for qnorms(p)");
		if (p==0) return -INFINITY;
		if (p==1) return INFINITY;
		static const boost::math::normal dist; return boost::math::quantile(dist,p);}
// non-standard normal distribution
	double cdf_norm_q(double mu, double sd, double x)		{	boost::math::normal dist(mu,sd);	return 1-boost::math::cdf(dist,x); }
	double cdf_norm_p(double mu, double sd, double x)		{	boost::math::normal dist(mu,sd);	return	 boost::math::cdf(dist,x); }
	double cdf_norm_2sided_pv(double mu,double sd,double x)	{	boost::math::normal dist(mu,sd);	return(1-boost::math::cdf(dist,std::fabs(x)))*2;}
	double qnorm(double mu, double sd, double p)			{
		if (p<0 || p>1) exit_error("p cannot be <0 or >1 for qnorms(p)");
		if (p==0) return -INFINITY;
		if (p==1) return INFINITY;
		boost::math::normal dist(mu,sd); return boost::math::quantile(dist,p);}
// binomial distribution
	double pdf_binomial   (double k, double p, double n)	{	if (k<0||n<=0) exit_error("binomial k<0||n<=0"); boost::math::binomial flip(my_round(n),p);	return boost::math::pdf(flip,my_round(k)); }
	double cdf_binomial_le(double k, double p, double n)	{	if (k<0||n<=0) exit_error("binomial k<0||n<=0"); boost::math::binomial flip(my_round(n),p);	return boost::math::cdf(flip,my_round(k)); }
	double cdf_binomial_ge(double k, double p, double n)	{	if (k<0||n<=0) exit_error("binomial k<0||n<=0"); boost::math::binomial flip(my_round(n),p);
																double rk=my_round(k); if (rk==0) return 1; else return 1-boost::math::cdf(flip,rk-1); }
#else
	#error Must define _STAT_USE_GSL or _STAT_USE_BOOST
#endif
double cdf_binomial_tail(double k, double p, double n)	{	if (k<0||n<=0) exit_error("binomial k<0||n<=0");
	if (k<=n*p)	return cdf_binomial_le(k,p,n);
	else		return cdf_binomial_ge(k,p,n); }

//cout << cdf_norm_p(0,1, INFINITY) << endl;		// 1
//cout << cdf_norm_p(0,1,-INFINITY) << endl;		// 0
//cout << pdf_norms( INFINITY) << endl;				// 0
//cout << pdf_norms(-INFINITY) << endl;				// 0
//cout << cdf_chisq_1df_q( INFINITY) << endl;		// 0
//cout << cdf_chisq_1df_q(-INFINITY) << endl;		// 1
//cout << cdf_t_2sided_pv(100, INFINITY) << endl;	// 0
//cout << cdf_t_2sided_pv(100,-INFINITY) << endl;	// 0

double prior_strength(int N, int k)
{
	if (k<0) exit_error("For prior strength calculation, k cannot be <0");
	if (N<0) exit_error("For prior strength calculation, N cannot be <0");
	if (k<2) return 0;
	if (N<1) return 0;
	static const std::map<int, std::map<int,double,std::greater<int>>, std::greater<int> > table = {
		{ 1,    { {2,24},   {3,4.83}, {4,1.17}, {5,0.65}, {6,0.45}, {8,0.27}, {16,0.11 } } } ,
		{ 2,    { {2,24},   {3,2.19}, {4,0.91}, {5,0.56}, {6,0.40}, {8,0.26}, {16,0.10 } } } ,
		{ 4,    { {2,24},   {3,1.72}, {4,0.81}, {5,0.52}, {6,0.38}, {8,0.25}, {16,0.10 } } } ,
		{ 10,   { {2,8.57}, {3,1.52}, {4,0.77}, {5,0.50}, {6,0.37}, {8,0.24}, {16,0.10 } } } ,
		{ 100,  { {2,6.19}, {3,1.42}, {4,0.74}, {5,0.49}, {6,0.37}, {8,0.24}, {16,0.10 } } }
	};
	return table.lower_bound(N)->second.lower_bound(k)->second;
}

// --------------------- truncated normal distribution ---------------------

/*/ http://en.wikipedia.org/wiki/Truncated_normal_distribution
C++:
cout<<pdf_tnorm_UpperTail(0.12,0.23,0,0.45)<<std::endl;
cout<<cdf_tnorm_UpperTail(0.12,0.23,0,0.45)<<std::endl;
cout<<pdf_tnorm_LowerTail(0.12,0.23,1,0.45)<<std::endl;
cout<<cdf_tnorm_LowerTail(0.12,0.23,1,0.45)<<std::endl;
cout<<pdf_tnorm(0.12,0.23,0,1,0.45)<<std::endl;
cout<<cdf_tnorm(0.12,0.23,0,1,0.45)<<std::endl;

cout << cdf_tnorm(0,1,0,1, INFINITY) << ' ';
cout << cdf_tnorm(0,1,0,1,-INFINITY) << ' ';
cout << cdf_tnorm_UpperTail(0,1,0, INFINITY) << ' ';
cout << cdf_tnorm_UpperTail(0,1,0,-INFINITY) << ' ';
cout << cdf_tnorm_LowerTail(0,1,0, INFINITY) << ' ';
cout << cdf_tnorm_LowerTail(0,1,0,-INFINITY) << endl;
// output 1 0 1 0 1 0
 
cout << pdf_tnorm(0,1,0,1, INFINITY) << ' ';
cout << pdf_tnorm(0,1,0,1,-INFINITY) << ' ';
cout << pdf_tnorm_UpperTail(0,1,0, INFINITY) << ' ';
cout << pdf_tnorm_UpperTail(0,1,0,-INFINITY) << ' ';
cout << pdf_tnorm_LowerTail(0,1,0, INFINITY) << ' ';
cout << pdf_tnorm_LowerTail(0,1,0,-INFINITY) << endl;
// output 0 0 0 0 0 0
 
R:
library(msm)
dtnorm(0.45,0.12,0.23,0,Inf) # 0.8864199
ptnorm(0.45,0.12,0.23,0,Inf) # 0.8917503
dtnorm(0.45,0.12,0.23,-Inf,1)# 0.6197135
ptnorm(0.45,0.12,0.23,-Inf,1)# 0.9243856
dtnorm(0.45,0.12,0.23,0,1)   # 0.8865025
ptnorm(0.45,0.12,0.23,0,1)   # 0.8918334
*/
 
double pdf_tnorm(double mu, double sd, double a, double b, double x)
{
	if (x<a) return 0;
	if (x>b) return 0;
	double xi = (x-mu)/sd;
	double alpha = (a-mu)/sd;
	double beta = (b-mu)/sd;
	double Z = cdf_norms(beta) - cdf_norms(alpha);
	return pdf_norms(xi) / ( sd * Z );
}

double cdf_tnorm(double mu, double sd, double a, double b, double x)
{
	if (x<a) return 0;
	if (x>b) return 1;
	double xi = (x-mu)/sd;
	double alpha = (a-mu)/sd;
	double beta = (b-mu)/sd;
	double Z = cdf_norms(beta) - cdf_norms(alpha);
	return (cdf_norms(xi) - cdf_norms(alpha))/ Z ;
}

double pdf_tnorm_UpperTail(double mu, double sd, double a, double x)
{
	if (x<a) return 0;
	double xi = (x-mu)/sd;
	double alpha = (a-mu)/sd;
	double Z = 1 - cdf_norms(alpha);
	return pdf_norms(xi) / ( sd * Z );
}

double cdf_tnorm_UpperTail(double mu, double sd, double a, double x)
{
	if (x<a) return 0;
	double xi = (x-mu)/sd;
	double alpha = (a-mu)/sd;
	double Z = 1 - cdf_norms(alpha);
	return (cdf_norms(xi) - cdf_norms(alpha))/ Z ;
}

double pdf_tnorm_LowerTail(double mu, double sd, double b, double x)
{
	if (x>b) return 0;
	double xi = (x-mu)/sd;
	double beta = (b-mu)/sd;
	double Z = cdf_norms(beta);
	return pdf_norms(xi) / ( sd * Z );
}

double cdf_tnorm_LowerTail(double mu, double sd, double b, double x)
{
	if (x>b) return 1;
	double xi = (x-mu)/sd;
	double beta = (b-mu)/sd;
	double Z = cdf_norms(beta);
	return cdf_norms(xi)/Z ;
}

// --------------------- Kolmogorov distribution ---------------------

void _Kolmogorov_mMultiply(double *A,double *B,double *C,int m)
{
	int i,j,k; double s;
	for(i=0;i<m;i++)
		for(j=0; j<m; j++)
		{	s=0.; for(k=0;k<m;k++) s+=A[i*m+k]*B[k*m+j]; C[i*m+j]=s; }
}

void _Kolmogorov_mPower(double *A,int eA,double *V,int *eV,int m,int n)
{
	double *B;int eB,i;
	if (n==1) { for (i=0;i<m*m;i++) V[i]=A[i]; *eV=eA; return; }
	_Kolmogorov_mPower(A,eA,V,eV,m,n/2);
	B=(double*)malloc((m*m)*sizeof(double));
	_Kolmogorov_mMultiply(V,V,B,m);
	eB=2*(*eV);
	if(n%2==0){for(i=0;i<m*m;i++) V[i]=B[i]; *eV=eB;}
	else {_Kolmogorov_mMultiply(A,B,V,m); *eV=eA+eB;}
	if (V[(m/2)*m+(m/2)]>1e140) {for(i=0;i<m*m;i++) V[i]=V[i]*1e-140;*eV+=140;} free(B);
}

// George Marsaglia, Wai Wan Tsang, Jingbo Wang. Evaluating Kolmogorov's Distribution. J Stat Software (2003) 8:18.
double Kolmogorov_K(int n,double d)
{
	int k,m,i,j,g,eH,eQ;
	double h,s,*H,*Q;
	//OMIT NEXT LINE IF YOU REQUIRE >7 DIGIT ACCURACY IN THE RIGHT TAIL
	s=d*d*n; if(s>7.24||(s>3.76&&n>99)) return 1-2*exp(-(2.000071+.331/sqrt(n)+1.409/n)*s);
	k=(int)(n*d)+1; m=2*k-1; h=k-n*d; H=(double*)malloc((m*m)*sizeof(double)); Q=(double*)malloc((m*m)*sizeof(double)); for(i=0;i<m;i++) for(j=0;j<m;j++)
		if(i-j+1<0) H[i*m+j]=0; else H[i*m+j]=1;
	for(i=0;i<m;i++) {H[i*m]-=pow(h,i+1); H[(m-1)*m+i]-=pow(h,(m-i));} H[(m-1)*m]+=(2*h-1>0?pow(2*h-1,m):0);
	for(i=0;i<m;i++) for(j=0;j<m;j++)
		if(i-j+1>0) for(g=1;g<=i-j+1;g++) H[i*m+j]/=g;
	eH=0;
	_Kolmogorov_mPower(H,eH,Q,&eQ,m,n);
	s=Q[(k-1)*m+k-1];
	for (i=1;i<=n;i++) {s=s*i/n; if(s<1e-140) {s*=1e140; eQ-=140;}} s*=pow(10.,eQ); free(H); free(Q); return s;
}

double Kolmogorov_K(double n1, double n2, double d)
{
	return Kolmogorov_K( n1*n2/(n1+n2), d);
}

// George Marsaglia, Wai Wan Tsang, Jingbo Wang. Evaluating Kolmogorov's Distribution. J Stat Software (2003) 8:18.
// I changed n from int to double, and it does asymtotic only.
double Kolmogorov_PVal(double n, double d)
{
	double s = d*d*n;
	long double core = -(2.000071+.331/sqrt(n)+1.409/n)*s;
	//	return 1 - 2 * exp(core); // original Kolmogorov K(n, d) = Pr(Dn < d) = 1 - PValue
	return exp(log(2)+core);
}

// For two-sample test, n = n.x*n.y/(n.x + n.y). Ref: R::ks.test
double Kolmogorov_PVal(double n1, double n2, double d)
{
	double n = n1*n2 / (n1+n2);
	return Kolmogorov_PVal(n,d);
}

// George Marsaglia, Wai Wan Tsang, Jingbo Wang. Evaluating Kolmogorov's Distribution. J Stat Software (2003) 8:18.
// I changed n from int to double, and it does asymtotic only. return -log10(p), good for very small p.
double Kolmogorov_PVal_NegLog(double n, double d)
{
	double s = d*d*n;
	long double core = -(2.000071+.331/sqrt(n)+1.409/n)*s;
	return - (log(2)+core) * log10(M_E);
}

// For two-sample test, n = n.x*n.y/(n.x + n.y). Ref: R::ks.test
double Kolmogorov_PVal_NegLog(double n1, double n2, double d)
{
	double n = n1*n2 / (n1+n2);
	return Kolmogorov_PVal_NegLog(n,d);
}

/* ============== p=0 problem ============
 int tbl_gtp[3][2]={595,367,628,693,178,298};
 int genetic_model[3]={ 1,2,3 }; // additive
 double x2,p;
 Cochran_Armitage_test_for_trend(tbl_gtp,genetic_model,3,x2,p);
 std::cout << "x2=" << x2 << " p=" << p << '\n';
 x2=85.7132 p=0 if calculated by boost, but correct (p=2.08023e-20) by GSL !
 http://jean-pierre.moreau.pagesperso-orange.fr/Cplus/chi2_cpp.txt compute p=0
 http://jean-pierre.moreau.pagesperso-orange.fr/Cplus/distri_cpp.txt compute qx(prob[X>=x2])=1 (strange, and px=p=5.33278e-17)
  ======================================== */

// The following 4 functions were taken from
// http://www.mathsisfun.com/combinatorics/combinations-permutations.html
// http://www.mathsisfun.com/combinatorics/combinations-permutations-calculator.html
// A permutation is an ordered combination. Unless specified, repetition is not allowed.
// n is the number of items to choose from, r is number of items one can choose.
// n must >=r , but these function does not check for the range of r for faster computation.
// it's very important to use long double inside the functions, but the return could be double.

long double num_combinations(int n, int r) // num_combinations = n! / r! (n-r)! = n(n-1)..(r+1) / (n-r)(n-r-1)..1
{
	int n_r=n-r;
	if (r < n_r) std::swap(r, n_r);
	long double result=1;
	for (int i=n;   i>r; --i) result*=i;
	for (int i=n_r; i>1; --i) result/=i;
	return result;
}

long double num_permutations(int n, int r) // num_permutations = n! / (n-r)! = n(n-1)..(n-r+1)
{
	long double result=1;
	for (int i=n; i>n-r; --i) result*=i;
	return result;
}

long double num_combinations_rep(int n, int r) // num_combinations with repitition = (n+r-1)! / r! (n-1)!
{
	int n_1   = n-1;
	int n_1pr = n-1+r;
	if (r < n_1) std::swap(r, n_1);
	long double result=1;
	for (int i=n_1pr; i>r; --i) result*=i;
	for (int i=n_1;   i>1; --i) result/=i;
	return result;
}

long double num_permutations_rep(int n, int r) // num_permutations with repitition = n^r
{
	long double result = 1;
	for (int i=0; i<r; ++i) result*=n;
	return result;
}

// -------- log10 version, slower but will not go out of bound ---------

double log10_num_combinations(int n, int r) // num_combinations = n! / r! (n-r)! = n(n-1)..(r+1) / (n-r)(n-r-1)..1
{
	int n_r=n-r;
	if (r < n_r) std::swap(r, n_r);
	double result = 0;
	for (int i=n;   i>r; --i) result+=log10(i);
	for (int i=n_r; i>1; --i) result-=log10(i);
	return result;
}

double log10_num_permutations(int n, int r) // num_permutations = n! / (n-r)! = n(n-1)..(n-r+1)
{
	double result = 0;
	for (int i=n; i>n-r; --i) result+=log10(i);
	return result;
}

double log10_num_combinations_rep(int n, int r) // num_combinations with repitition = (n+r-1)! / r! (n-1)!
{
	int n_1   = n-1;
	int n_1pr = n-1+r;
	if (r < n_1) std::swap(r, n_1);
	double result = 0;
	for (int i=n_1pr; i>r; --i) result+=log10(i);
	for (int i=n_1;   i>1; --i) result-=log10(i);
	return result;
}

double log10_num_permutations_rep(int n, int r) // num_permutations with repitition = n^r
{
	double result = 0;
	for (int i=0; i<r; ++i) result+=log10(n);
	return result;
}

//  cout << num_combinations(5, 3) << '\n';		// 10
//  cout << num_permutations(5, 3) << '\n';		// 60
//  cout << num_combinations_rep(5, 3) << '\n';	// 35
//  cout << num_permutations_rep(5, 3) << '\n';	// 125
//  cout << num_combinations(5, 2) << '\n';		// 10
//  cout << num_permutations(5, 2) << '\n';		// 20
//  cout << num_combinations_rep(5, 2) << '\n';	// 15
//  cout << num_permutations_rep(5, 2) << '\n';	// 25

// the following 2 taken from
// Understanding the relationship between risks and odds ratios.
// Shrier I, Steele R. Clin J Sport Med. 2006 Mar;16(2):107-10.
// P0 is prevalence of disease in the control group
// (I think the author means unexposed group, because RR = P1 / P0)
double convert_OR_to_RR(double OR, double P0)
{
	return OR / ( (1-P0) + (P0*OR) );
}

double convert_RR_to_OR(double RR, double P0)
{
	return (RR - P0*RR) / (1 - P0*RR);
}

// ------------------------------------------- sample size and power -------------------------------------------

// The following 2 functions are based on the article
// A simple method of sample size calculation for linear and logistic regression.
// Hsieh FY, Bloch DA & Larsen MD. Statist. Med. 17, 1623Ð1634 (1998) formula 1 and 2
// none of them were tested yet
double sample_size_logistic_regression_cont(double P1, double beta_star, double alpha, double power)
{
	double b = qnorms(alpha/2)+qnorms(1-power);
	double d = P1 * (1-P1) * beta_star * beta_star;
	return b*b / d;
}

// P1 and P2 are the event rates at X=0 and 1, respectively. B is the proportion of the sample with X=1.
double sample_size_logistic_regression_bin(double P1, double P2, double B, double alpha, double power)
{
	double P = (1-B)*P1 + B*P2;
	double b=qnorms(alpha/2)*sqrt(P*(1-P)/B) + qnorms(1-power)*sqrt(P1*(1-P1)+P2*(1-P2)*(1-B)/B);
	double d=(P1-P2)*(P1-P2)*(1-B);
	return b*b / d;
}

// Kelsey et al. Methods in Observational Epidemiology (1986) Table 10-15
// online http://www.sph.emory.edu/~cdckms/sample%20size%202%20grps%20case%20control.html
//     r+1    pb(1-pb)(Zbeta+Zalpha/2)^2
// n = --- x ---------------------------
//      r            (p1-p2)^2
// r is ratio of controls to cases for CC / ratio of unexposed to exposed samlpes for cohort
// p1 is % of cs exposed for CC / % of exposed samples develop the dz for cohort
// p0 is % of ct exposed for CC / % of unexposed samples develop the dz for cohort
// OR = odds ratio
// return n, number of cases for CC / number of exposed individuals studied for cohort
double sample_size_Kelsey_bin(double p0, double OR, double r, double alpha, double power)
{
	double p1 = p0*OR / (p0*(OR-1)+1);
	double pCs = 1/(r+1);
	double pCt = r/(r+1);
	double pb = pCs*p1 + pCt*p0;	// weighted average of p0 and p1
	double v = pb*(1-pb);			// a measure of variability (similar to square of standard deviation)
	double di = p1-p0;				// d*
	double b = qnorms(alpha/2)+qnorms(1-power); // beta=1-power
	return ((r+1)/r) * v * b*b / (di*di);
}

// sd = standard deviation of the studied variable in the population, di = d* (diff in % or means)
double sample_size_Kelsey_cont(double sd, double di, double r, double alpha, double power)
{
	double b = qnorms(alpha/2)+qnorms(1-power); // beta=1-power
	return ((r+1)/r) * sd * sd * b*b / (di*di);
}

// the following 2 functions: give n, return power
double power_Kelsey_bin(double p0, double OR, int n, double r, double alpha)
{
	double p1 = p0*OR / (p0*(OR-1)+1);
	double pCs = 1/(r+1);
	double pCt = r/(r+1);
	double pb = pCs*p1 + pCt*p0;	// weighted average of p0 and p1
	double v = pb*(1-pb);			// a measure of variability (similar to square of standard deviation)
	double di = p1-p0;				// d*
	return 1 - cdf_norms_p( -di * sqrt(n*r/(r+1)/v) - qnorms(alpha/2) );
}

double power_Kelsey_cont(double sd, double di, int n, double r, double alpha)
{
	return 1 - cdf_norms_p( -di * sqrt(n*r/(r+1)/sd/sd) - qnorms(alpha/2) );
}

// taken from the javascript from http://www.stat.ubc.ca/~rollin/stats/ssize/b2.html
// Reference: The cal are the customary ones based on the normal approximation to the binomial distribution.
// Ref: Categorical Data - Estimation of Sample Size and Power for Comparing Two Binomial
// Proportions in Bernard Rosner's Fundamentals of Biostatistics.
int sample_size_proportion_2samples(double p1, double p2, int side, double alpha, double power)
{
	double q1 = 1-p1;
	double q2 = 1-p2;
	double pb = (p1+p2)/2;
	double qb = 1-pb;
	double del = std::fabs(p1-p2);
	double za = qnorms(alpha/side);
	double zb = qnorms(1-power);
	double sqrt2pq = sqrt(pb * qb * 2);
	double sqrtsumpqs = sqrt(p1 * q1 + p2 * q2);
	double sss=((sqrt2pq * za + sqrtsumpqs * zb)/del);
	return my_round(sss*sss); // prv int(sss*sss+.5), should not be a problem since sss*sss>=0
}

/* removed because they are in libfbj_math.hpp
struct genetic_power_analysis_return_type {
	double power_dom;
	double power_rec;
	double power_add;
}; */

// power of caes-control SNP association study on binary trait
// Menashe I, Rosenberg PS, Chen BE. PGA: power cal for CC genet assoc analyses. BMC Genet (2008) 13;9:36.
// results not right, need Debug.
genetic_power_analysis_return_type power_CC_SNP_assoc(double pd, double p1, double Dp, double r2,
															 double pr, double rr1, double rr2,
															 int n1, int n0,
															 double alpha)
{
	if (Dp==-9 && r2==-9) exit_error("Need either D' or r2.");
	if (Dp!=-9 && r2!=-9) exit_error("Please provide either D' or r2, not both.");
	double D;
	double Dmax = fmin( p1*(1-pd) , pd*(1-p1) );
	if (Dp!=-9)
	{
		D = Dp * Dmax;
	}
	if (r2!=-9)
	{
		double denominator = p1*(1-pd)*pd*(1-p1);
		double r2max = Dmax * Dmax / denominator;
		D = sqrt( r2 * r2max * denominator);
	}
	double h[2][2];		// freq for hap AB Ab aB ab
	h[1][1] = pd*p1 + D;
	h[1][0] = pd*(1-p1) - D;
	h[0][1] = (1-pd)*p1 - D;
	h[0][0] = (1-pd)*(1-p1) + D;
	double PrXaXb[3][3];// prob of gtp Xa=i Xb=j (# of copies)
	PrXaXb[0][0] = h[0][0]*h[0][0];
	PrXaXb[0][1] = h[0][0]*h[0][1]*2;
	PrXaXb[0][2] = h[0][1]*h[0][1];
	PrXaXb[1][0] = h[0][0]*h[1][0]*2;
	PrXaXb[1][1] = h[1][1]*h[0][0]*2 + h[0][1]*h[1][0]*2;
	PrXaXb[1][2] = h[1][1]*h[0][1]*2;
	PrXaXb[2][0] = h[1][0]*h[1][0];
	PrXaXb[2][1] = h[1][1]*h[1][0]*2;
	PrXaXb[2][2] = h[1][1]*h[1][1];
	double r[3];		// penetrance given copy of diease causative allele = i
	r[0] = pr / ( (1-pd)*(1-pd) + 2*pd*(1-pd)*rr1 + pd*pd*rr2 );
	r[1] = r[0]*rr1;
	r[2] = r[0]*rr2;
	double PYdXb[2][3] = { {0,0,0},{0,0,0} };	// prob of Y=d Xb=i
	for (int i=0;i<3;++i)
		for (int j=0;j<3;++j)
			PYdXb[1][i] += r[j]*PrXaXb[j][i];
	for (int i=0;i<3;++i)
		for (int j=0;j<3;++j)
			PYdXb[0][i] += (1-r[j])*PrXaXb[j][i];
	double P00 = PYdXb[0][0]/(1-pr);
	double P10 = PYdXb[0][1]/(1-pr);
	double P20 = PYdXb[0][2]/(1-pr);
	double P01 = PYdXb[1][0]/pr;
	double P11 = PYdXb[1][1]/pr;
	double P21 = PYdXb[1][2]/pr;
	genetic_power_analysis_return_type result;
	double mean = log(n0*P00)+log(n1*P11+n1*P21)-log(n1*P01)-log(n0*P10+n0*P20);
	double variance = 1/(n1*P01) + 1/(n0*P10+n0*P20) + 1/(n0*P00) + 1/(n1*P11+n1*P21);
	double sd = sqrt(variance);
	result.power_dom = 1 - cdf_norm_p(mean, sd, qnorm(mean, sd, 1-alpha/2));
	return result;
}

// ------------------------------------------- Cochran_Armitage_test_for_trend -------------------------------------------

template <typename T1,typename T2,typename T3>
void Cochran_Armitage_test_for_trend(T1 N[][2], T2 t[], int k, T3& x2, T3& p)
// http://en.wikipedia.org/wiki/Cochran-Armitage_test_for_trend
{
	x2=0; p=1;
	int R1=0,R2=0;
	std::vector<int> C(k);
	for (int i=0;i<k;++i) { R1+=N[i][0]; R2+=N[i][1]; C[i]=N[i][0]+N[i][1]; }
	if (!R1 || !R2) return;
	int N_=R1+R2;
	double T=0,Vpart1=0,Vpart2=0;
	for (int i=0;i<k;++i)
	{	T+= t[i]*((int)N[i][0]*R2-(int)N[i][1]*R1);
		Vpart1+=t[i]*t[i]*C[i]*(N_-C[i]);
	}
	for (int i=0;i<k-1;++i)
		for (int j=i+1;j<k;++j)
			Vpart2+=t[i]*t[j]*C[i]*C[j];
	double Var=R1*R2*(Vpart1-2*Vpart2)/N_;
	if (!Var) return;
	double z=T/sqrt(Var);
	x2=z*z;
	p=cdf_chisq_1df_q(x2);
}
template void Cochran_Armitage_test_for_trend<unsigned,double,double>(unsigned N[][2], double t[], int k, double& x2, double& p);
template void Cochran_Armitage_test_for_trend<int,double,double>(int N[][2], double t[], int k, double& x2, double& p);
template void Cochran_Armitage_test_for_trend<double,double,double>(double N[][2], double t[], int k, double& x2, double& p);

template <typename T1,typename T2,typename T3>
void Cochran_Armitage_test_for_trend(T1 N[], T2 t[], int k, T3& x2, T3& p)
// http://en.wikipedia.org/wiki/Cochran-Armitage_test_for_trend
{
	x2=0; p=1;
	int R1=0,R2=0;
	std::vector<int> C(k);
	for (int i=0;i<k;++i) { R1+=N[i]; R2+=N[k+i]; C[i]=N[i]+N[k+i]; }
	if (!R1 || !R2) return;
	int N_=R1+R2;
	double T=0,Vpart1=0,Vpart2=0;
	for (int i=0;i<k;++i)
	{	T+= t[i]*((int)N[i]*R2-(int)N[k+i]*R1);
		Vpart1+=t[i]*t[i]*C[i]*(N_-C[i]);
	}
	for (int i=0;i<k-1;++i)
		for (int j=i+1;j<k;++j)
			Vpart2+=t[i]*t[j]*C[i]*C[j];
	double Var=R1*R2*(Vpart1-2*Vpart2)/N_;
	if (!Var) return;
	double z=T/sqrt(Var);
	x2=z*z;
	p=cdf_chisq_1df_q(x2);
}
template void Cochran_Armitage_test_for_trend<unsigned,double,double>(unsigned N[], double t[], int k, double& x2, double& p);
template void Cochran_Armitage_test_for_trend<int,double,double>(int N[], double t[], int k, double& x2, double& p);
template void Cochran_Armitage_test_for_trend<double,double,double>(double N[], double t[], int k, double& x2, double& p);

/*============== test data =============
STATA: paste data, then "ptrend n1 n2 a" ( http://www.stata.com/support/faqs/stat/trend.html )
a	n1	n2
1	19	1
2	31	5
3	67	21
 output: Chi2(1) for trend = 4.546,  pr>chi2 = 0.0330
int cc[3][2]={19,1,31,5,67,21};
double tr[3]={0,0.5,1}; // 0,0.5,1 or 0,1,2 or 1,2,3 are all the same
double x2,p;
Cochran_Armitage_test_for_trend(cc,tr,3,x2,p); // expect: x2=4.54646 p=0.0329869
  ======================================== */

// ------------------------------------------- chi-square test -------------------------------------------

// likelihodd ratio test for a contingency table, from sas
template <typename T1,typename T2, typename T3>
void chi_square_test_lr(T1* obs, T2* expct,int num, int df, T3& x2, T3& p)
{
	double chi2=0;
	for (int i=0;i<num;++i)
	{
		if (obs[i]==0) continue; // do nothing, tested with SAS FREQ 
		if (expct[i]) chi2 += log((double)obs[i]/(double)expct[i]) * (double)obs[i];
	}
	x2=chi2*2;
	p=cdf_chisq_q(x2,df);
}
template void chi_square_test_lr<unsigned, double, double>(unsigned* obs, double* expct,int num, int df, double& x2, double& p);
template void chi_square_test_lr<int, double, double>(int* obs, double* expct,int num, int df, double& x2, double& p);
template void chi_square_test_lr<double, double, double>(double* obs, double* expct,int num, int df, double& x2, double& p);

template <typename T1,typename T2, typename T3>
void chi_square_test_pearson(T1* obs, T2* expct,int num, int df, T3& x2, T3& p)
{
	double chi2=0;
	for (int i=0;i<num;++i)
		if (expct[i]) chi2+=((double)expct[i]-(double)obs[i])*((double)expct[i]-(double)obs[i])/(double)expct[i];
	x2=chi2;
	p=cdf_chisq_q(x2,df);
}
template void chi_square_test_pearson<unsigned, double, double>(unsigned* obs, double* expct,int num, int df, double& x2, double& p);
template void chi_square_test_pearson<int, double, double>(int* obs, double* expct,int num, int df, double& x2, double& p);
template void chi_square_test_pearson<double, double, double>(double* obs, double* expct,int num, int df, double& x2, double& p);

// http://v8doc.sas.com/sashtml/stat/chap28/sect19.htm , make more sense than wikipedia
template <typename T1,typename T2, typename T3>
void chi_square_test_yate(T1* obs, T2* expct,int num, int df, T3& x2, T3& p)
{
	double chi2=0;
	for (int i=0;i<num;++i)
		if (expct[i])
		{
			double diff = fmax(0, std::fabs((double)expct[i]-(double)obs[i])-0.5 );
			chi2 += diff*diff/(double)expct[i];
		}
	x2=chi2;
	p=cdf_chisq_q(x2,df);
}
template void chi_square_test_yate<unsigned, double, double>(unsigned* obs, double* expct,int num, int df, double& x2, double& p);
template void chi_square_test_yate<int, double, double>(int* obs, double* expct,int num, int df, double& x2, double& p);
template void chi_square_test_yate<double, double, double>(double* obs, double* expct,int num, int df, double& x2, double& p);

// ------------------------------------------- for contingency table -------------------------------------------
// all statistic function use T* so that both C array and C++ vector can be used for input. This is also the
// reason why I don't put the data inside the class, because otherwise I need to specify how to store the data.
// However, because of this user needs to make sure prepare_from() is called before using any stat functions!
// I don't put ContingencyTableStatType inside struct or class for shorter names, but it pollutes the namespace.
// All functions have tested for VAL & ASE, but pvl is note tested because SAS doesn't show it. Be careful.

/* removed because they are in libfbj_math.hpp
enum ContingencyTableStatType
{
	ChiSq_Ps, ChiSq_Yt, ChiSq_LR, UncertCR, UncertRC, UncertSm, LambdaCR, LambdaRC, LambdaSm,	// nominal
	KendllTb, StuartTc, GK_gamma, SomerDcr, SomerDrc, PearsonR, SpearmnR, JTtrendS,	MHChiSqT,	// ordinal
	BowkersQ, KappaSmp,																			// RxR
	CochranQ,																					// other
	NumStats																					// # of stats
};

template <typename T>
class ContingencyTable
{
public:
	std::vector<double> val, ase, pvl, dof; // value, Asymptotic Standard Error, p-value, degree of freedom
	void InitAllStats() { val.assign(NumStats,-2); ase.assign(NumStats,0); pvl.assign(NumStats,-1); dof.assign(NumStats,0); }
	void InitStat(ContingencyTableStatType type) { val[type]=-2; ase[type]=0; pvl[type]=-1; dof[type]=0; }
	
	unsigned nr, nc, nd;			// #row, #col, #data
	double n;						// sum of all cells
	std::vector<double> sumr, sumc; // sum of row, col
	ContingencyTable():nr(0),nc(0),nd(0),n(0) {}
	void prepare_from(T* data, unsigned r, unsigned c)
	{
		nr=r; nc=c; nd=r*c;
		sumr.assign(nr,0);
		sumc.assign(nc,0);
		n=0;
		for (unsigned ij=0,i=0; i<nr; ++i)
			for (unsigned  j=0; j<nc; ++j, ++ij)
			{
				n      +=data[ij];
				sumr[i]+=data[ij];
				sumc[j]+=data[ij];
			}
		InitAllStats();
	}
}; */

// has segmentation fault: 11
template <typename T>
void RxC_both_nominal(T* N, ContingencyTable<T>& ctab)
{
	ctab.InitStat(ChiSq_Ps);
	ctab.InitStat(ChiSq_Yt);
	ctab.InitStat(ChiSq_LR);
	ctab.InitStat(UncertCR);
	ctab.InitStat(UncertRC);
	ctab.InitStat(UncertSm);
	ctab.InitStat(LambdaCR);
	ctab.InitStat(LambdaRC);
	ctab.InitStat(LambdaSm);
	unsigned& nr=ctab.nr;
	unsigned& nc=ctab.nc;
	unsigned& nd=ctab.nd;
	double& n=ctab.n;
	
	std::vector<double> E(nr * nc, 0);
	unsigned er=nr; // effective number of rows
	unsigned ec=nc; // effective number of cols
	for (unsigned r=0; r<nr; ++r) if (ctab.sumr[r]==0) --er;
	for (unsigned c=0; c<nc; ++c) if (ctab.sumc[c]==0) --ec;
	if (er<2 || ec<2) return;
	for (unsigned ij=0,r=0; r<nr; ++r)
		for (unsigned  c=0; c<nc; ++c, ++ij)
			E[ij] = ctab.sumr[r] * ctab.sumc[c] / n;
	ctab.dof[ChiSq_Ps]=(er-1)*(ec-1);
	ctab.dof[ChiSq_Yt]=(er-1)*(ec-1);
	ctab.dof[ChiSq_LR]=(er-1)*(ec-1);
	chi_square_test_pearson(&N[0], &E[0], nd, ctab.dof[ChiSq_Ps], ctab.val[ChiSq_Ps], ctab.pvl[ChiSq_Ps]);
	chi_square_test_yate   (&N[0], &E[0], nd, ctab.dof[ChiSq_Yt], ctab.val[ChiSq_Yt], ctab.pvl[ChiSq_Yt]);
	chi_square_test_lr	   (&N[0], &E[0], nd, ctab.dof[ChiSq_LR], ctab.val[ChiSq_LR], ctab.pvl[ChiSq_LR]);
	
	// Uncertainty coefficient. treatment of zeros same as SAS, see examples at http://www4.stat.ncsu.edu/~dzhang2/st744/table3.9.lst.txt
	{
		double Hx=0, Hy=0, Hxy=0;
		for (unsigned r=0; r<nr; ++r)	Hx -= (ctab.sumr[r]/n)*log(ctab.sumr[r]/n);
		for (unsigned c=0; c<nc; ++c)	Hy -= (ctab.sumc[c]/n)*log(ctab.sumc[c]/n);
		for (unsigned ij=0,r=0; r<nr; ++r)
			for (unsigned  c=0; c<nc; ++c, ++ij)
				if (N[ij]) Hxy -= (N[ij]/n)*log(N[ij]/n);
		double v = Hx + Hy - Hxy;
		double Hx_pl_Hy = Hx+Hy;
		ctab.val[UncertCR] = v / Hy;
		ctab.val[UncertRC] = v / Hx;
		ctab.val[UncertSm] = ( 2*v )/( Hx+Hy );
		for (unsigned ij=0,r=0; r<nr; ++r)
			for (unsigned  c=0; c<nc; ++c, ++ij)
				if (N[ij])
				{
					double a;
					a = Hy*log(N[ij]/ctab.sumr[r]) + (Hx-Hxy)*log(ctab.sumc[c]/n);			ctab.ase[UncertCR] += N[ij]*a*a;
					a = Hx*log(N[ij]/ctab.sumc[c]) + (Hy-Hxy)*log(ctab.sumr[r]/n);			ctab.ase[UncertRC] += N[ij]*a*a;
					a = Hxy*log(ctab.sumr[r]*ctab.sumc[c]/(n*n)) - Hx_pl_Hy*log(N[ij]/n);	ctab.ase[UncertSm] += N[ij]*a*a;
				}
		ctab.ase[UncertCR] = sqrt( ctab.ase[UncertCR] / (n*n*Hy*Hy*Hy*Hy) );
		ctab.ase[UncertRC] = sqrt( ctab.ase[UncertRC] / (n*n*Hx*Hx*Hx*Hx) );
		ctab.ase[UncertSm] = sqrt( 4 * ctab.ase[UncertSm] / ( n * n * Hx_pl_Hy * Hx_pl_Hy * Hx_pl_Hy * Hx_pl_Hy ) );
	}
	
	// Lambda (prob improvment in predicting C|R or R|C or Symmetric)
	{
		std::vector<T> ri(nr,0);
		std::vector<T> cj(nc,0);
		std::vector<unsigned> li(nr,0);
		std::vector<unsigned> kj(nc,0);
		for (unsigned ij=0,i=0; i<nr; ++i)
			for (unsigned  j=0; j<nc; ++j, ++ij)
			{
				if (N[ij]>ri[i]) { ri[i]=N[ij]; li[i]=j; }
				if (N[ij]>cj[j]) { cj[j]=N[ij]; kj[j]=i; }
			}
		double r=0, c=0; // r = max(ctab.sumc[j]);
		unsigned l=-1, k=-1; // l is the unique value of j such that r = ctab.sumc[j], same as k.
		for (unsigned  j=0; j<nc; ++j) if (ctab.sumc[j]>r) { r=ctab.sumc[j]; l=j; }
		for (unsigned  i=0; i<nr; ++i) if (ctab.sumr[i]>c) { c=ctab.sumr[i]; k=i; }
		double e=0, f=0, g=0, h=0;
		for (unsigned  i=0; i<nr; ++i) { e+=ri[i]; if (li[i]==l) f+=ri[i]; }
		for (unsigned  j=0; j<nc; ++j) { g+=cj[j]; if (kj[j]==k) h+=cj[j]; }
		double LambdaCR_var = ((n-e)/((n-r)*(n-r)*(n-r))) * (e+r-2*f);
		double LambdaRC_var = ((n-g)/((n-c)*(n-c)*(n-c))) * (g+c-2*h);
		ctab.val[LambdaCR] = (e-r) / (n-r);
		ctab.val[LambdaRC] = (g-c) / (n-c);
		ctab.ase[LambdaCR] = sqrt(LambdaCR_var);
		ctab.ase[LambdaRC] = sqrt(LambdaRC_var);
		double w = 2*n-r-c;
		double v = 2*n-e-g;
		double x = f+h+ri[k]+cj[l];
		double y = 8*n-w-v-2*x;
		ctab.val[LambdaSm] = (w-v)/w;
		double z = 0;
		for (unsigned ij=0,i=0; i<nr; ++i)
			for (unsigned  j=0; j<nc; ++j, ++ij)
				if (j==li[i] && i==kj[j]) z += N[ij];
		double LambdaSm_var = (w*v*y - 2*w*w*(n-z) - 2*v*v*(n-N[nc*k+l])) / (w*w*w*w);
		ctab.ase[LambdaSm] = sqrt(LambdaSm_var);
	}
}

// http://v8doc.sas.com/sashtml/stat/chap28/sect18.htm
template <typename T>
void RxC_both_ordinal(T* N, ContingencyTable<T>& ctab)
{
	ctab.InitStat(KendllTb);
	ctab.InitStat(StuartTc);
	ctab.InitStat(GK_gamma);
	ctab.InitStat(SomerDcr);
	ctab.InitStat(SomerDrc);
	ctab.InitStat(PearsonR);
	ctab.InitStat(SpearmnR);
	ctab.InitStat(JTtrendS);
	ctab.InitStat(MHChiSqT);
	unsigned& nr=ctab.nr;
	unsigned& nc=ctab.nc;
	unsigned& nd=ctab.nd;
	double& n=ctab.n;

	// prepare
	std::vector<double> A(nd,0);
	std::vector<double> D(nd,0);
	double  P=0, Q=0, Ty=0, Tx=0;
	for (unsigned ij=0,i=0; i<nr; ++i)
	{
		for (unsigned  j=0; j<nc; ++j, ++ij)
		{
			for (unsigned kl=0,k=0; k<nr; ++k)
			{
				for (unsigned  l=0; l<nc; ++l, ++kl)
				{
					if (j!=l && i==k) { Tx += N[ij] * N[kl]; continue; }
					if (j==l && i!=k) { Ty += N[ij] * N[kl]; continue; }
					if (i==k || j==l) continue;
					if (i>k)	{ if (j>l) A[ij] += N[kl]; else D[ij] += N[kl]; }
					else		{ if (j<l) A[ij] += N[kl]; else D[ij] += N[kl]; }
				}
			}
		}
	}
	if (n<2) return;
	for (unsigned ij=0,i=0; i<nr; ++i)
		for (unsigned  j=0; j<nc; ++j, ++ij)
		{
			P += N[ij] * A[ij];
			Q += N[ij] * D[ij];
		}
	double P_pl_Q = P+Q;
	double P_mi_Q = P-Q;
	double ndd = 0;
	for (unsigned ij=0,i=0; i<nr; ++i)
		for (unsigned  j=0; j<nc; ++j, ++ij)
			ndd += N[ij]*(A[ij]-D[ij])*(A[ij]-D[ij]);
	double ndd_pq2n = ndd - P_mi_Q*P_mi_Q/n;
	double wr=n*n; for (unsigned i=0; i<nr; ++i) wr-=ctab.sumr[i]*ctab.sumr[i];
	double wc=n*n; for (unsigned j=0; j<nc; ++j) wc-=ctab.sumc[j]*ctab.sumc[j];
	
	// Goodman and Kruskal's gamma
	if (P_pl_Q)
	{
		ctab.val[GK_gamma] = P_mi_Q / P_pl_Q;
		double GK_gamma_var = 0;
		for (unsigned ij=0,i=0; i<nr; ++i)
			for (unsigned  j=0; j<nc; ++j, ++ij)
			{
				double a = Q*A[ij]-P*D[ij];
				GK_gamma_var += N[ij]*a*a;
			}
		GK_gamma_var = 16*GK_gamma_var/(P_pl_Q*P_pl_Q*P_pl_Q*P_pl_Q);
		double gamma_vh0 = (4 / (P_pl_Q*P_pl_Q)) * ndd_pq2n; // variance under H0 (GK_gamma=0)
		ctab.ase[GK_gamma] = sqrt(GK_gamma_var);
		ctab.pvl[GK_gamma] = cdf_norms_2sided_pv(GK_gamma/sqrt(gamma_vh0));
	}
	
	// Somers' D C|R
	if (wr)
	{
		ctab.val[SomerDcr] = P_mi_Q/wr;
		double SomerDcr_var=0;
		for (unsigned ij=0,i=0; i<nr; ++i)
			for (unsigned  j=0; j<nc; ++j, ++ij)
			{
				double a = wr*(A[ij]-D[ij]) - P_mi_Q*(n-ctab.sumr[i]);
				SomerDcr_var += N[ij]*a*a;
			}
		SomerDcr_var = 4*SomerDcr_var/(wr*wr*wr*wr);
		double SomerDcr_vh0 = (4/(wr*wr)) * ndd_pq2n;
		ctab.ase[SomerDcr] = sqrt(SomerDcr_var);
		ctab.pvl[SomerDcr] = cdf_norms_2sided_pv(ctab.val[SomerDcr]/sqrt(SomerDcr_vh0));
	}
	
	// Somers' D R|C
	if (wc)
	{
		ctab.val[SomerDrc] = P_mi_Q/wc;
		double SomerDrc_var=0;
		for (unsigned ij=0,i=0; i<nr; ++i)
			for (unsigned  j=0; j<nc; ++j, ++ij)
			{
				double a = wc*(A[ij]-D[ij]) - P_mi_Q*(n-ctab.sumc[j]);
				SomerDrc_var += N[ij]*a*a;
			}
		SomerDrc_var = 4*SomerDrc_var/(wc*wc*wc*wc);
		double SomerDrc_vh0= (4/(wc*wc)) * ndd_pq2n;
		ctab.ase[SomerDrc] = sqrt(SomerDrc_var);
		ctab.pvl[SomerDrc] = cdf_norms_2sided_pv(ctab.val[SomerDrc]/sqrt(SomerDrc_vh0));
	}

	// Kendall's Tau-b
	{
		double w = sqrt (wr*wc);
		ctab.val[KendllTb] = P_mi_Q / w;
		double KendllTb_var = 0;
		for (unsigned ij=0,i=0; i<nr; ++i)
			for (unsigned  j=0; j<nc; ++j, ++ij)
			{
				double a = 2*w*(A[ij]-D[ij]) + ctab.val[KendllTb]*(ctab.sumr[i]*wc+ctab.sumc[j]*wr);
				KendllTb_var += N[ij]*a*a;
			}
		KendllTb_var = (KendllTb_var - n*n*n*ctab.val[KendllTb]*ctab.val[KendllTb]*(wr+wc)*(wr+wc)) / (w*w*w*w);
		double KendllTb_vh0= (4/(wr*wc)) * ndd_pq2n;
		ctab.ase[KendllTb] = sqrt(KendllTb_var);
		ctab.pvl[KendllTb] = cdf_norms_2sided_pv(ctab.val[KendllTb]/sqrt(KendllTb_vh0));
	}
	
/*
	// http://en.wikipedia.org/wiki/Kendall_tau_rank_correlation_coefficient tested with STATA
	double tiedX_notY = Tx/2;
	double tiedY_notX = Ty/2;
	double concordant = P/2;
	double discordant = Q/2;
	SomerDcr = (concordant-discordant)/(concordant+discordant+tiedY_notX);
	SomerDrc = (concordant-discordant)/(concordant+discordant+tiedX_notY);
	double n0 = n*(n-1)/2;
	tau_a = (concordant-discordant)/n0; // same as STATA
	double Za = 1.5*(concordant-discordant)/sqrt(n*(n-1)*(n*2+5)/2);
	tau_a_p=cdf_norms_2sided_pv(Za); // not tested, but for 2x2 pv similar to chi-square test
	double n1=0, n2=0;
	for (unsigned r=0; r<nr; r++) n1+=ctab.sumr[r]*(ctab.sumr[r]-1)/2;
	for (unsigned c=0; c<nc; c++) n2+=ctab.sumc[c]*(ctab.sumc[c]-1)/2;
	KendllTb = (concordant-discordant)/sqrt((n0-n1)*(n0-n2)); // same as STATA
	double vt=0, vu=0, v1a=0, v1b=0, v2a=0, v2b=0;
	for (unsigned r=0; r<nr; r++) vt+=ctab.sumr[r]*(ctab.sumr[r]-1)*(ctab.sumr[r]*2+5);
	for (unsigned c=0; c<nc; c++) vu+=ctab.sumc[c]*(ctab.sumc[c]-1)*(ctab.sumc[c]*2+5);
	for (unsigned r=0; r<nr; r++) v1a+=ctab.sumr[r]*(ctab.sumr[r]-1);
	for (unsigned c=0; c<nc; c++) v1b+=ctab.sumc[c]*(ctab.sumc[c]-1);
	double v1=v1a*v1b/(2*n*(n-1));
	for (unsigned r=0; r<nr; r++) v2a+=ctab.sumr[r]*(ctab.sumr[r]-1)*(ctab.sumr[r]-2);
	for (unsigned c=0; c<nc; c++) v2b+=ctab.sumc[c]*(ctab.sumc[c]-1)*(ctab.sumc[c]-2);
	double v2=v2a*v2b/(9*n*(n-1)*(n-2));
	double v0=n*(n-1)*(n*2+5);
	double Zb=(concordant-discordant);
	if (Zb>0) Zb=fmax(0,Zb-1); else Zb=fmin(0,Zb+1); // continuity correction, I guessed & tested with STATA ktau
	double denominator = (v0-vt-vu)/18+v1+v2;
	if (denominator>0) KendllTb_p=cdf_norms_2sided_pv(Zb/sqrt(denominator));	// same as STATA ktau command
*/
	
	// Stuart's Tau-c
	{
		double m = fmin(nr,nc);
		ctab.val[StuartTc] = (m*P_mi_Q)/(n*n*(m - 1));
		double StuartTc_var = (4*m*m/((m-1)*(m-1)*n*n*n*n)) * ndd_pq2n; // var=vh0
		ctab.ase[StuartTc] = sqrt(StuartTc_var);
		ctab.pvl[StuartTc] = cdf_norms_2sided_pv(ctab.val[StuartTc]/sqrt(StuartTc_var));
	}
	
	// Spearman Rank Correlation Coefficient, same as SAS, STATA and
	// http://department.obg.cuhk.edu.hk/ResearchSupport/Spearman.asp
	{
		std::vector<double> R0(nr,0); // SCORE=TABLE
		std::vector<double> C0(nc,0);
		std::vector<double> R1(nr,0); // SCORE=RANK
		std::vector<double> C1(nc,0);
		std::vector<double> R(nr,0);
		std::vector<double> C(nc,0);
		std::vector<double> vij(nd,0);
		std::vector<double> wij(nd,0);
		std::vector<double> zij(nd,0);
		double Rb=0, Cb=0;
		for (unsigned i=0; i<nr; ++i)
		{
			R0[i] = i+1;
			for (unsigned k=0; k<i; ++k) R1[i] += ctab.sumr[k];
			R1[i] += (ctab.sumr[i])/2; // same as SAS and as the paper bellow but wrong.
			// see "Sampling Behavior of Test for Corr in 2-Way Contingency Tables. J of Am Stat Assoc. 72(358):309-315" for the wrong equation
			// and SAS manual for the correct equation: (ctab.sumr[i]+1)/2
			R[i] = R1[i] - n/2;
			Rb += R0[i]*ctab.sumr[i];
		}
		for (unsigned j=0; j<nc; ++j)
		{
			C0[j] = j+1;
			for (unsigned l=0; l<j; ++l) C1[j] += ctab.sumc[l];
			C1[j] += (ctab.sumc[j])/2; // see above
			C[j] = C1[j] - n/2;
			Cb += C0[j]*ctab.sumc[j];
		}
		Rb/=n; Cb/=n;
		double v=0;
		for (unsigned ij=0,i=0; i<nr; ++i)
			for (unsigned  j=0; j<nc; ++j, ++ij)
				v += N[ij]*R[i]*C[j];
		double F=n*n*n, G=F;
		for (unsigned i=0; i<nr; ++i) F -= ctab.sumr[i]*ctab.sumr[i]*ctab.sumr[i];
		for (unsigned j=0; j<nc; ++j) G -= ctab.sumc[j]*ctab.sumc[j]*ctab.sumc[j];
		double w=sqrt(F*G)/12;
		ctab.val[SpearmnR] = v/w;
		
		double zb=0;
		double vb=0;
		for (unsigned ij=0,i=0; i<nr; ++i)
			for (unsigned  j=0; j<nc; ++j, ++ij)
			{
				double a1=0, a2=0, a3=0, a4=0;
				for (unsigned  l=0; l<nc; ++l) a1+= N[nc*i+l]*C[l]; a1/=2;
				for (unsigned  k=0; k<nr; ++k) a2+= N[nc*k+j]*R[k]; a2/=2;
				for (unsigned  k=i+1; k<nr; ++k) for (unsigned  l=0; l<nc; ++l) a3+= N[nc*k+l]*C[l];
				for (unsigned  l=j+1; l<nc; ++l) for (unsigned  k=0; k<nr; ++k) a4+= N[nc*k+l]*R[k];
				vij[ij] = n*(R[i]*C[j]+a1+a2+a3+a4);
				wij[ij] = (n/(-96*w)) * (F*ctab.sumc[j]*ctab.sumc[j]+G*ctab.sumr[i]*ctab.sumr[i]);
				zij[ij] = w*vij[ij] - v*wij[ij];
				zb += N[ij]*zij[ij];
				vb += N[ij]*vij[ij];
			}
		zb /= n;
		vb /= n;
		double SpearmnR_Vh0=0, SpearmnR_Var=0;
		for (unsigned ij=0,i=0; i<nr; ++i)
			for (unsigned  j=0; j<nc; ++j, ++ij)
			{
				SpearmnR_Vh0 += N[ij]*(vij[ij]-vb)*(vij[ij]-vb);
				SpearmnR_Var += N[ij]*(zij[ij]-zb)*(zij[ij]-zb);
			}
		SpearmnR_Vh0 /= n*n*w*w;
		SpearmnR_Var /= n*n*w*w*w*w;
		ctab.ase[SpearmnR] = sqrt(SpearmnR_Var);
		ctab.pvl[SpearmnR] = cdf_norms_2sided_pv(ctab.val[SpearmnR]/sqrt(SpearmnR_Vh0)); // SAS, not tested
//		double SpearmnR_t = ctab.val[SpearmnR] * sqrt( (n-2) / (1 - ctab.val[SpearmnR] * ctab.val[SpearmnR]) );
//		ctab.pvl[SpearmnR] = cdf_t_2sided_pv(n-2, SpearmnR_t); // STATA
		
		// Pearson corr
		std::vector<double> bij(nd,0);
		double ssr=0, ssc=0, ssrc=0;
		for (unsigned ij=0,i=0; i<nr; ++i)
			for (unsigned  j=0; j<nc; ++j, ++ij)
			{
				ssr += N[ij]*(R0[i]-Rb)*(R0[i]-Rb);
				ssc += N[ij]*(C0[j]-Cb)*(C0[j]-Cb);
				ssrc+= N[ij]*(R0[i]-Rb)*(C0[j]-Cb);
			}
		for (unsigned ij=0,i=0; i<nr; ++i)
			for (unsigned  j=0; j<nc; ++j, ++ij)
				bij[ij] = (R0[i]-Rb)*(R0[i]-Rb)*ssc + (C0[j]-Cb)*(C0[j]-Cb)*ssr;
		double ssrssc = ssr*ssc;
		double PearsonR_w = sqrt(ssrssc);
		ctab.val[PearsonR] = ssrc / PearsonR_w;
		double PearsonR_vh0 = 0;
		for (unsigned ij=0,i=0; i<nr; ++i)
			for (unsigned  j=0; j<nc; ++j, ++ij)
				PearsonR_vh0 += N[ij]*(R0[i]-Rb)*(R0[i]-Rb)*(C0[j]-Cb)*(C0[j]-Cb);
		PearsonR_vh0 = (PearsonR_vh0 - ssrc*ssrc/n) / ssrssc;
		double PearsonR_var = 0;
		for (unsigned ij=0,i=0; i<nr; ++i)
			for (unsigned  j=0; j<nc; ++j, ++ij)
			{
				double a = PearsonR_w*(R0[i]-Rb)*(C0[j]-Cb) - bij[ij]*ssrc/(2*PearsonR_w);
				PearsonR_var += N[ij]*a*a;
			}
		PearsonR_var = PearsonR_var / (ssrssc*ssrssc);
		ctab.ase[PearsonR] = sqrt (PearsonR_var);
		ctab.pvl[PearsonR] = cdf_norms_2sided_pv(ctab.val[PearsonR]/sqrt(PearsonR_vh0));
		
		// Mantel-Haenszel chi-square test
		ctab.val[MHChiSqT]=(n-1)*ctab.val[PearsonR]*ctab.val[PearsonR];
		ctab.pvl[MHChiSqT]=cdf_chisq_1df_q(ctab.val[MHChiSqT]);
		ctab.dof[MHChiSqT]=1;
	}
	
	// Jonckheere's trend test (Jonckheere–Terpstra test)
	// http://en.wikipedia.org/wiki/Jonckheere's_trend_test
	{
		double P=0, Q=0;
		for (unsigned ij=0,i=0; i<nr; ++i)
			for (unsigned  j=0; j<nc; ++j, ++ij)
				if (N[ij])
					for (unsigned k=0;k<nr;++k)
						for (unsigned l=j+1;l<nc;++l)
						{	if		(k<i) Q+= N[ij]*N[nc*k+l];
							else if (k>i) P+= N[ij]*N[nc*k+l];  }
		ctab.val[JTtrendS] = P-Q;
		double t2=0, t3=0, u2=0, u3=0;
		for (unsigned j=0; j<nc; ++j) { t2+=ctab.sumc[j]*ctab.sumc[j]; t3+=ctab.sumc[j]*ctab.sumc[j]*ctab.sumc[j]; }
		for (unsigned i=0; i<nr; ++i) { u2+=ctab.sumr[i]*ctab.sumr[i]; u3+=ctab.sumr[i]*ctab.sumr[i]*ctab.sumr[i]; }
		double JTtrendS_var = (2*(n*n*n-t3-u3)+3*(n*n-t2-u2)+5*n)/18 + (t3-3*t2+2*n)*(u3-3*u2+2*n)/(9*n*(n-1)*(n-2)) + (t2-n)*(u2-n)/(2*n*(n-1));
		ctab.pvl[JTtrendS] = cdf_norms_2sided_pv(ctab.val[JTtrendS]/sqrt(JTtrendS_var));
	}
}

template <typename T>
void RxR_agreement(T* N, ContingencyTable<T>& ctab)
{
	ctab.InitStat(BowkersQ);
	ctab.InitStat(KappaSmp);
	unsigned& nr=ctab.nr;
	unsigned& nc=ctab.nc;
	double& n=ctab.n;

	// Bowker's test of Symmetry, tested with manual calculation.
	{
		ctab.val[BowkersQ]=0;
		for (unsigned i=0; i<nr; ++i)
			for (unsigned j=i+1; j<nr; ++j)
				if (N[nc*i+j]+N[nc*j+i]) ctab.val[BowkersQ] += ((double)N[nc*i+j]-N[nc*j+i]) * ((double)N[nc*i+j]-N[nc*j+i]) / (N[nc*i+j]+N[nc*j+i]);
		ctab.pvl[BowkersQ] = cdf_chisq_q(ctab.val[BowkersQ],nr*(nr-1)/2);
		ctab.dof[BowkersQ] = nr*(nr-1)/2;
	}
	
	// simple Kappa coefficient
	{
		double Po=0, Pe=0, A=0, B=0;
		for (unsigned i=0; i<nr; ++i) Po += N[nc*i+i]/n;
		for (unsigned i=0; i<nr; ++i) Pe += ctab.sumr[i]*ctab.sumc[i]/(n*n);
		ctab.val[KappaSmp] = (Po-Pe)/(1-Pe);
		double OneMiK = 1-ctab.val[KappaSmp];
		for (unsigned i=0; i<nr; ++i) A  += N[nc*i+i] * (n-(ctab.sumr[i]+ctab.sumc[i])*OneMiK) * (n-(ctab.sumr[i]+ctab.sumc[i])*OneMiK);
		A /= (n*n*n);
		for (unsigned i=0; i<nr; ++i)
			for (unsigned j=0; j<nc; ++j)
				if (i!=j) B+= N[nc*i+j] * (ctab.sumc[i]+ctab.sumr[j]) * (ctab.sumc[i]+ctab.sumr[j]);
		B = B * OneMiK * OneMiK / (n*n*n);
		double C = ctab.val[KappaSmp] - Pe * OneMiK;
		C*=C;
		double KappaSmp_var = (A+B-C) / ((1-Pe)*(1-Pe)*n);
		ctab.ase[KappaSmp] = sqrt(KappaSmp_var);
		double KappaSmp_vh0 = 0;
		for (unsigned i=0; i<nr; ++i)
			KappaSmp_vh0 += ctab.sumr[i] * ctab.sumc[i] * (ctab.sumr[i]+ctab.sumc[i]);
		KappaSmp_vh0 = (Pe + Pe*Pe - KappaSmp_vh0/(n*n*n)) / ((1-Pe)*(1-Pe)*n);
		ctab.pvl[KappaSmp] = cdf_norms_2sided_pv(ctab.val[KappaSmp]/sqrt(KappaSmp_vh0));
	}
}

// http://en.wikipedia.org/wiki/Cochran's_Q_test
// tested with http://support.sas.com/documentation/cdl/en/procstat/63104/HTML/default/viewer.htm#procstat_freq_sect034.htm
template <typename T>
void RxC_BinDV_CochranQ(T* N, ContingencyTable<T>& ctab)
{
	ctab.InitStat(CochranQ);
	unsigned& nr=ctab.nr;
	unsigned& nc=ctab.nc;
	double& n=ctab.n;

	{
		double a=0, b=0;
		for (unsigned j=0; j<nc; ++j) a+= (ctab.sumc[j]-n/nc)*(ctab.sumc[j]-n/nc);
		for (unsigned i=0; i<nr; ++i) b+=  ctab.sumr[i]*(nc-ctab.sumr[i]);
		ctab.val[CochranQ]=nc*(nc-1)*a/b;
		ctab.pvl[CochranQ]=cdf_chisq_q(ctab.val[CochranQ],nc-1);
	}
}

// below instantiate the above 4 functions. 
// http://stackoverflow.com/questions/495021/why-can-templates-only-be-implemented-in-the-header-file

template void RxC_both_nominal<unsigned>(unsigned* N, ContingencyTable<unsigned>& ctab);
template void RxC_both_ordinal<unsigned>(unsigned* N, ContingencyTable<unsigned>& ctab);
template void RxR_agreement<unsigned>(unsigned* N, ContingencyTable<unsigned>& ctab);
template void RxC_BinDV_CochranQ<unsigned>(unsigned* N, ContingencyTable<unsigned>& ctab);

template void RxC_both_nominal<int>(int* N, ContingencyTable<int>& ctab);
template void RxC_both_ordinal<int>(int* N, ContingencyTable<int>& ctab);
template void RxR_agreement<int>(int* N, ContingencyTable<int>& ctab);
template void RxC_BinDV_CochranQ<int>(int* N, ContingencyTable<int>& ctab);

template void RxC_both_nominal<double>(double* N, ContingencyTable<double>& ctab);
template void RxC_both_ordinal<double>(double* N, ContingencyTable<double>& ctab);
template void RxR_agreement<double>(double* N, ContingencyTable<double>& ctab);
template void RxC_BinDV_CochranQ<double>(double* N, ContingencyTable<double>& ctab);

// ------------------------------------------- Wilcoxon_signed_rank_test -------------------------------------------

template <typename T1,typename T2>
void Wilcoxon_signed_rank_test(const std::multiset<T1>& z, T2& sumW, T2& p)
{
	
	// critical values table can be taken from Biostatistical Analysis, 4/E Jerrold H. Zar 1999 appendixB12
	// or calculated using algrithm in http://www.psych.cornell.edu/darlington/wilcoxon/wilcox4.htm
	static std::map<int,std::map<int,double> > pv; // pv[N][W] = 1 tailed p value
	static std::mutex mt;
	if (pv.empty())
	{
		mt.lock();
		if (pv.empty())
		{
			// initialize data for N=1
			int num=1;								// N
			int len=2;								// length of vector
			std::vector<double> nOutcome(2,1);			// vector for "number of outcome"
			nOutcome.reserve(5051);					// reserve space for N=100
			for (int k=2;k<=100;++k)				// calcualte N up to 100
			{
				std::vector<double> vector_2=nOutcome;	// vector for shifting down N
				++num; len+=num;
				nOutcome.insert(nOutcome.end(),  num,0);
				vector_2.insert(vector_2.begin(),num,0);
				double total=pow(2., num);			// total number of outcomes
				double le=0;						// number of outcomes that S<=l
				for (int s=0;s<len;++s)
				{
					nOutcome[s]+=vector_2[s];
					le+=nOutcome[s];
					pv[num][s]=le/total;			// 1 tailed p value
				}
			}
		}
		mt.unlock();
	}
	
	std::map<T1,double> r; // r[std::fabs(z)] = rank
	typename std::multiset<T1>::const_iterator it1;
	for (it1=z.begin(); it1!=z.end(); ++it1) if (*it1!=0) ++r[std::fabs(*it1)];
	int N=0; // total number of non-zeros less than *it2
	double adjtV=0;
	typename std::map<T1,double>::iterator it2;
	for (it2=r.begin(); it2!=r.end(); ++it2)
	{
		int n = it2->second;
		if (n>1)
		{
			it2->second = N + (1+n)/2. ;
			adjtV += (n*n*n-n) / 48.;
		}
		else // n==1
			it2->second = N + 1 ;
		N += n;
	}
	// now N become the total number of non-zero observations
	if (N<2) { sumW=0; p=1; return; }
	
	// calculate W=min(W+,W-)
	double Wp = 0;					// W+
	for (it1=z.begin(); it1!=z.end(); ++it1) if (*it1 > 0) Wp+=r[*it1];
	double Wn = N*(N+1)/2 - Wp;		// W-
	double W = Wp < Wn ? Wp : Wn ;	// W = min(W+,W-)
	sumW= Wp - Wn;
	
	// calcualte p value
	//	static const boost::math::normal dist; // standard normal
	if (N<=100) // see simulation, even N=100 normal asymt is not as good as exact test but very close
		p=pv[N][ceil(W)]*2;
	else
	{
		double Zscore;
		
		// Z test method 1
		double meanW = N*(N+1)/4. ;	// mean of W
		double stdvW = sqrt( N*(N+1)*(N+N+1)/24. -adjtV ); // SD(W). -adjtV is adjustment for non-zero ties. Hollander-Wolfe (1973), Lehman(1975)
		Zscore= (W-meanW)/stdvW; // some argue to correct for continuity: -.5 http://www.fon.hum.uva.nl/Service/Statistics/Signed_Rank_Test.html */
		
		/*/ Z test method 2 http://faculty.vassar.edu/lowry/ch12a.html adjust for non-zero ties? t1err lower than method 1
		double cont = sumW > 0 ? -.5 : .5 ; // correction for continuity
		double stdvW = sqrt( N*(N+1)*(N+N+1)/6. ); // SD(W).
		Zscore = (sumW + cont) / stdvW; //*/
		
		/*/ Z test method 3 http://www.psych.cornell.edu/darlington/normscor.htm, t1err much lower than others plus weird
		double M = 0;
		for (it1=z.begin(); it1!=z.end(); ++it1)
		{
			if (*it1 > 0)		M+=gsl_cdf_ugaussian_Pinv( r[ *it1]/(z.size()+1) );
			else if (*it1 < 0 )	M+=gsl_cdf_ugaussian_Pinv( r[-*it1]/(z.size()+1) );
			else continue;
		}
		Zscore = M / sqrt(z.size()); //*/
		
		p = cdf_norms_2sided_pv(Zscore);
		// p = gsl_cdf_ugaussian_Q(std::fabs(Zscore))*2; // p = ( 1 - boost::math::cdf(dist,std::fabs(Zscore)) ) * 2;		
	}
}

template void Wilcoxon_signed_rank_test<unsigned,double>(const std::multiset<unsigned>& z, double& sumW, double& p);
template void Wilcoxon_signed_rank_test<int,double>(const std::multiset<int>& z, double& sumW, double& p);
template void Wilcoxon_signed_rank_test<double,double>(const std::multiset<double>& z, double& sumW, double& p);

// ------------------------------------------- general functions for containers -------------------------------------------


// ------------------------------------------- Values class -------------------------------------------

const std::vector<std::string> STAT::StatName = {
	"MIN","MAX","MEAN","SUM","N","TSS","TSS0","RMS","SPL_VAR","POP_VAR","SPL_SD","POP_SD","PRD","M3","SKEW","ADJ_SKEW","SEM","POP_VMR","CV","RSD","UCV","URSD","NUM_NAN","MEAN_ABS","MED_ABS","MEAN_POS","MED_POS",
	"IQR","QCD","MAD","BWMV","MEDIAN","P0.01","P0.05","P0.1","P0.5","P1","P2","P2.5","P5","P9","P10","P25","P75","P90","P91","P95","P97.5","P98","P99","P99.5","P99.9","P99.95","P99.99","SD_LO","SD_UP","SORTED",
	"ABP_MC","ABP_LF","ABP_UF","ABP_LO","ABP_UO","Adil_LF","Adil_UF","Kimber_LF","Kimber_UF",
	"NUM_STAT"
};

STAT::StatType to_StatType(const std::string& name)
{
	std::string u = to_upper_copy(name);
	if (str_startsw(u,"PCT")) u = "P" + u.substr(3);
	for (int i=0; i<STAT::NUM_STAT; ++i)
		if (u==STAT::StatName[i]) return static_cast<STAT::StatType>(i);
	exit_error(name+" is not a StatID. Available StatIDs are "+str_of_container(STAT::StatName,','));
	return STAT::NUM_STAT; // should not happen
}

std::string	to_StatName(const STAT::StatType type)
{
	return STAT::StatName[static_cast<int>(type)];
}

template <typename T>
int Values<T>::num_ge(const T& threshold)
{
	if (std::isnan(threshold)) exit_error("Values::num_ge cutoff is NaN.");
	if (num_ge_cnt.empty()) {
		for (each_element(data,it)) if (!std::isnan(*it)) ++num_ge_cnt[*it];
		int N=0;
		for (each_element(num_ge_cnt,it_vc)) { it_vc->second += N; N=it_vc->second; }
	}
	if (exist_element(num_ge_cnt,threshold)) return num_ge_cnt[threshold];
	else
	{
		auto it=num_ge_cnt.upper_bound(threshold);
		if (it == num_ge_cnt.begin()) return 0;
		--it;
		return it->second;
	}
}

template <typename T>
int Values<T>::num_le(const T& threshold)
{
	if (std::isnan(threshold)) exit_error("Values::num_le cutoff is NaN.");
	if (num_le_cnt.empty()) {
		for (each_element(data,it)) if (!std::isnan(*it)) ++num_le_cnt[*it];
		int N=0;
		for (each_element(num_le_cnt,it_vc)) { it_vc->second += N; N=it_vc->second; }
	}
	if (exist_element(num_le_cnt,threshold)) return num_le_cnt[threshold];
	else
	{
		auto it=num_le_cnt.upper_bound(threshold);
		if (it == num_le_cnt.begin()) return 0;
		--it;
		return it->second;
	}
}

template <typename T>
void Values<T>::do_min() {
	if (get(STAT::SORTED)==1) {
		if (!data.empty())	result[STAT::MIN]=data.front();
		else				result[STAT::MIN]=std::numeric_limits<double>::signaling_NaN();
	}
	else {
		double res=DBL_MAX;
		for (auto &i:data) if (!std::isnan(i) && i<res) res=i;
		if (res!=DBL_MAX)	result[STAT::MIN]=res;
		else				result[STAT::MIN]=std::numeric_limits<double>::signaling_NaN();
	}
}
template <typename T>
void Values<T>::do_max() {
	if (get(STAT::SORTED)==1) {
		auto it = rbegin();
		while (it!=rend() && std::isnan(*it)) ++it;
		if (it!=rend())	result[STAT::MAX] = *it;
		else			result[STAT::MAX] = std::numeric_limits<double>::signaling_NaN();
	}
	else {
		double res=-DBL_MAX;
		for (auto &i:data) if (!std::isnan(i) && i>res) res=i;
		if (res!=-DBL_MAX)	result[STAT::MAX]=res;
		else				result[STAT::MAX]=std::numeric_limits<double>::signaling_NaN();
	}
}
template <typename T>
void Values<T>::do_rms() { // Root mean square
	if (!exist_element(result,STAT::N)) do_n();
	if (!exist_element(result,STAT::TSS0)) get(STAT::TSS0);
	if (result[STAT::N]>0)	result[STAT::RMS]=sqrt(result[STAT::TSS0]/result[STAT::N]);
	else					result[STAT::RMS]=std::numeric_limits<double>::signaling_NaN();
}
template <typename T>
void Values<T>::do_pop_vmr() { // variance-to-mean ratio
	if (!exist_element(result,STAT::POP_VAR)) do_pop_var();
	if (result[STAT::MEAN]!=0)	result[STAT::POP_VMR]=result[STAT::POP_VAR]/result[STAT::MEAN];
	else						result[STAT::POP_VMR]=std::numeric_limits<double>::signaling_NaN();
}
template <typename T>
void Values<T>::do_cv() { // [sample] coefficient of variation
	if (!exist_element(result,STAT::SPL_SD)) get(STAT::SPL_SD);
	if (result[STAT::MEAN]!=0)	result[STAT::CV]=result[STAT::SPL_SD]/result[STAT::MEAN];
	else						result[STAT::CV]=std::numeric_limits<double>::signaling_NaN();
}
template <typename T>
void Values<T>::do_mean_abs() { // mean of absolute value
	double offset = 0;	for (each_element(data,it)) if (!std::isnan(*it)) { offset=std::fabs(*it); break; }
	double m=0; int n=0;for (each_element(data,it)) if (!std::isnan(*it)) { m+=std::fabs(*it)-offset; ++n; }
	result[STAT::N]=n;
	result[STAT::NUM_NAN]=data.size()-n;
	if (!n) { result[STAT::MEAN_ABS]=std::numeric_limits<double>::signaling_NaN(); }
	else	{ result[STAT::MEAN_ABS]=m/n+offset; }
}
template <typename T>
void Values<T>::do_med_abs() {	// median of absolute value
	std::deque<T> abs_val;	// absolute values
	int n=0;
	for (each_element(data,it)) if (!std::isnan(*it)) { abs_val.push_back(std::fabs(*it)); ++n; }
	result[STAT::N]=n;
	result[STAT::NUM_NAN]=data.size()-n;
	if (!n) { result[STAT::MED_ABS]=std::numeric_limits<double>::signaling_NaN(); }
	else	{ std::sort(abs_val.begin(),abs_val.end()); result[STAT::MED_ABS]=median(abs_val); }
}
template <typename T>
void Values<T>::do_mean_pos() { // mean of positive value
	double offset = 0;	for (each_element(data,it)) if (!std::isnan(*it)&&*it>0) { offset=*it; break; }
	double m=0; int n=0;for (each_element(data,it)) if (!std::isnan(*it)&&*it>0) { m+=*it-offset; ++n; }
	if (!n) { result[STAT::MEAN_POS]=std::numeric_limits<double>::signaling_NaN(); }
	else	{ result[STAT::MEAN_POS]=m/n+offset; }
}
template <typename T>
void Values<T>::do_med_pos() {	// median of positive value
	std::deque<T> pos_val;	// absolute values
	int n=0;
	for (each_element(data,it)) if (!std::isnan(*it)&&*it>0) { pos_val.push_back(*it); ++n; }
	if (!n) { result[STAT::MED_POS]=std::numeric_limits<double>::signaling_NaN(); }
	else	{ std::sort(pos_val.begin(),pos_val.end()); result[STAT::MED_POS]=median(pos_val); }
}
template <typename T>
void Values<T>::do_n() {
	int ian=0, nan=0;
	for (each_element(data,it))
		if (std::isnan(*it)) ++nan; else ++ian;
	result[STAT::N]=ian;
	result[STAT::NUM_NAN]=nan;
}
template <typename T>
void Values<T>::do_spl_var() {
	if (!exist_element(result,STAT::TSS)) get(STAT::TSS);
	if (result[STAT::N]>1)	result[STAT::SPL_VAR] = result[STAT::TSS] / (result[STAT::N]-1);
	else					result[STAT::SPL_VAR] = std::numeric_limits<double>::signaling_NaN();
}
template <typename T>
void Values<T>::do_pop_var() {
	if (!exist_element(result,STAT::TSS)) get(STAT::TSS);
	if (result[STAT::N]>0)	result[STAT::POP_VAR] = result[STAT::TSS] / result[STAT::N];
	else					result[STAT::POP_VAR] = std::numeric_limits<double>::signaling_NaN();
}
template <typename T>
void Values<T>::do_sd1sided() {
	double med = get(STAT::MEDIAN);
	double N = get(STAT::N);
	if (!std::isnan(med) && N>=2)
	{
		double lo_sum=0, up_sum=0;
		for (auto &Xi:data)
			if (!std::isnan(Xi))
			{
				if (Xi<med) lo_sum += 2*(Xi-med)*(Xi-med);
				if (Xi>med) up_sum += 2*(Xi-med)*(Xi-med);
			}
		result[STAT::SD_LO] = sqrt(lo_sum/(N-1));
		result[STAT::SD_UP] = sqrt(up_sum/(N-1));
	}
	else
	{
		result[STAT::SD_LO] = std::numeric_limits<double>::signaling_NaN();
		result[STAT::SD_UP] = std::numeric_limits<double>::signaling_NaN();
	}
}
template <typename T>
void Values<T>::do_prd() { // product
	int size=0;
	T p=1;
	for (each_element(data,it)) if (!std::isnan(*it)) { p *= *it; ++size; }
	result[STAT::N]=size;
	if (size) result[STAT::PRD] = p;
	else result[STAT::PRD] = std::numeric_limits<double>::signaling_NaN();
}
template <typename T>
void Values<T>::sort() {
	if (get(STAT::SORTED)==1) return;
	size_t size_1 = data.size(); // bigger
	// data.erase(remove_if(data.begin(),data.end(),std::isnan<T>),data.end()); // compiled in old gcc but not the new one because it's not a template anymore
	// http://stackoverflow.com/questions/17574242/how-to-use-isnan-as-a-predicate-function-to-stdfind-if-c11
	data.erase(remove_if(data.begin(),data.end(),[](T d){return std::isnan(d);}),data.end());
	size_t size_2 = data.size(); // smaller
	result[STAT::NUM_NAN] = size_1-size_2;
	result[STAT::N] = size_2;
	if (!data.empty()) std::sort(data.begin(),data.end());
	if (size_1!=size_2) data.insert(data.end(),size_1-size_2,std::numeric_limits<T>::signaling_NaN());
	result[STAT::SORTED] = 1;
}
template <typename T>
void Values<T>::remove_nan() {
	if (get(STAT::SORTED)==1)
	{
		while (!data.empty() && std::isnan(data.back())) data.pop_back();
	}
	else
	{
		size_t size_1 = data.size(); // bigger
		// data.erase(remove_if(data.begin(),data.end(),std::isnan<T>),data.end()); // compiled in old gcc but not the new one because it's not a template anymore
		// http://stackoverflow.com/questions/17574242/how-to-use-isnan-as-a-predicate-function-to-stdfind-if-c11
		data.erase(remove_if(data.begin(),data.end(),[](T d){return std::isnan(d);}),data.end());
		size_t size_2 = data.size(); // smaller
		result[STAT::NUM_NAN] = size_1-size_2;
		result[STAT::N] = size_2;
		if (!data.empty()) std::sort(data.begin(),data.end());
		result[STAT::SORTED] = 1;
	}
}
template <typename T>
void Values<T>::rank() {
	if (ranks.empty())
	{
		remove_nan();
		ranks.assign(data.size(),0);
		int tie_bgn=0, tie_end=0;
		for (size_t i=0; i<data.size(); ++i)
		{
			ranks[i]=i+1;
			if (i)
			{
				if (data[i]==data[i-1])
				{
					++tie_end;
					if (i+1==data.size())
					{
						double mean = (tie_bgn+tie_end+2)/2.0;
						for (int j=tie_bgn; j<=tie_end; ++j) ranks[j]=mean;
					}
				}
				else
				{
					if (tie_bgn!=tie_end)
					{
						double mean = (tie_bgn+tie_end+2)/2.0;
						for (int j=tie_bgn; j<=tie_end; ++j) ranks[j]=mean;
					}
					tie_bgn=i;
					tie_end=i;
				}
			}
		}
	}
}
template <typename T>
void Values<T>::do_kimber() {	// SIQR BP by Kimber 1990. Ref: 2011 Contributions to Outlier Detection Methods (dissertation).
	T	Q1    = get(STAT::PCT25);	// Q1
	T	Q3    = get(STAT::PCT75);	// Q3
	T	med   = get(STAT::MEDIAN);	// median
	T	SIQRl = med-Q1;
	T	SIQRu = Q3-med;
	result[STAT::Kimber_LF] = Q1 - 3.0 * SIQRl;
	result[STAT::Kimber_UF] = Q3 + 3.0 * SIQRu; //*/
}
template <typename T>
void Values<T>::do_abp() {	// adjusted boxplot
	remove_nan();						// will call sort()
	int		jend  = get(STAT::N);		// N
	if (jend<4) exit_error("Sample size (<4) is too small for AdjustedBoxPlot.");
	
	// make IQR MEDIAN MEDIAN_LB MEDIAN_UB MEDIAN_K
	T	IQR	  = get(STAT::IQR);		// IQR
	T	Q1    = get(STAT::PCT25);	// Q1
	T	Q3    = get(STAT::PCT75);	// Q3
	T	med   = get(STAT::MEDIAN);	// median
	T	SIQRl = med-Q1;
	T	SIQRu = Q3-med;
	int	medlb = std::lower_bound(data.begin(),data.end(),med) - data.begin();
	int	medub = std::upper_bound(data.begin(),data.end(),med) - data.begin();
	int	medk  = medub-medlb;		// number of data[i]=med
	
	// make mcij, values to calculate MC
	std::deque<float> mcij;
	for (int i=0;i<medub;++i)
		for (int j=medlb;j<jend;++j)
		{
			if (i==j) continue;
			if (data[i]==med && data[j]==med)
			{
				if		(i+j+1< medk)	mcij.push_back(-1);
				else if (i+j+1==medk)	mcij.push_back(0);
				else					mcij.push_back(1);
			}
			else
			{
				mcij.push_back( ((data[j]-med)-(med-data[i]))/((float)data[j]-data[i]) );
			}
		}
	
	// calculate MC from mcij
	std::sort(mcij.begin(),mcij.end());
	float MC = median(mcij);		// MedCouple
	T lf,uf;
	
	// SIQR BP by Kimber 1990. Ref: 2011 Contributions to Outlier Detection Methods (dissertation).
	result[STAT::Kimber_LF] = Q1 - 3.0 * SIQRl;
	result[STAT::Kimber_UF] = Q3 + 3.0 * SIQRu; //*/
	
	/*/ earlier ABP by G. Brys 2005. A Robustification of Independent Component Analysis. -- by the same people who created the official ABP later
	if (MC>=0) { lf = Q1 - 1.5 * exp(-3.5*MC) * IQR ; uf = Q3 + 1.5 * exp(4.0*MC) * IQR; }
	else	   { lf = Q1 - 1.5 * exp(-4.0*MC) * IQR ; uf = Q3 + 1.5 * exp(3.5*MC) * IQR; } //*/
	
	// official ABP by M. Hubert 2007. An Adjusted Boxplot for Skewed Distributions. (tested by R adjboxStats {robustbase})
	if (MC>=0) { lf = Q1 - 1.5 * exp(-4.0*MC) * IQR ; uf = Q3 + 1.5 * exp(3.0*MC) * IQR; }
	else	   { lf = Q1 - 1.5 * exp(-3.0*MC) * IQR ; uf = Q3 + 1.5 * exp(4.0*MC) * IQR; } //*/
	
	/*/ modified ABP by Yinaze Herve Dovoedo 2011. Contributions to Outlier Detection Methods (dissertation). -- worse
	lf = med - 4.0 * exp(-2.0*MC) * SIQRl;
	uf = med + 4.0 * exp( 2.0*MC) * SIQRu; //*/
	
	// modified ABP by I. H. Adil 2015. A Modified Approach for Detection of Outliers. -- not stable
	double sk = get(STAT::SKEW);
	if (sk> 3.5) sk= 3.5;
	if (sk<-3.5) sk=-3.5;
	result[STAT::Adil_LF] = Q1 - 1.5 * exp(-sk*std::abs(MC)) * IQR;
	result[STAT::Adil_UF] = Q3 + 1.5 * exp( sk*std::abs(MC)) * IQR; //*/
	
	result[STAT::ABP_MC] = MC;
	result[STAT::ABP_LF] = lf;
	result[STAT::ABP_UF] = uf;
	auto itl = std::lower_bound(data.begin(), data.end(), lf); result[STAT::ABP_LO] = itl - data.begin();
	auto itu = std::upper_bound(data.begin(), data.end(), uf); result[STAT::ABP_UO] = data.end() - itu;
}
template <typename T>
void Values<T>::do_m3() {
	if (!exist_element(result,STAT::MEAN)) get(STAT::MEAN);
	if (result[STAT::N]>1)	result[STAT::M3] = tsn(data,result[STAT::MEAN],3)/(result[STAT::N]);
	else					result[STAT::M3] = std::numeric_limits<double>::signaling_NaN();
}
template <typename T>
void Values<T>::do_skew() { // http://itl.nist.gov/div898/handbook/eda/section3/eda35b.htm , results match STATA
	if (!exist_element(result,STAT::M3))		do_m3();
	if (!exist_element(result,STAT::POP_SD))	get(STAT::POP_SD);
	if (result[STAT::N]>2)	result[STAT::SKEW] = (result[STAT::M3]/pow(result[STAT::POP_SD],3));
	else					result[STAT::SKEW] = std::numeric_limits<double>::signaling_NaN();
}
template <typename T>
void Values<T>::do_adj_skew() { // http://itl.nist.gov/div898/handbook/eda/section3/eda35b.htm
	if (!exist_element(result,STAT::SKEW))		do_skew();
	double N=result[STAT::N];
	if (!std::isnan(result[STAT::SKEW]))	result[STAT::ADJ_SKEW] = result[STAT::SKEW] * (sqrt(N*(N-1))/(N-1));
	else									result[STAT::ADJ_SKEW] = std::numeric_limits<double>::signaling_NaN();
}
template <typename T>
void Values<T>::do_sem() {	// [sample] standard error of the mean
	if (!exist_element(result,STAT::SPL_SD))	get(STAT::SPL_SD);
	if (result[STAT::N]>0)	result[STAT::SEM] = result[STAT::SPL_SD]/sqrt(result[STAT::N]);
	else					result[STAT::SEM] = std::numeric_limits<double>::signaling_NaN();
}
template <typename T>
void Values<T>::do_bwmv() {	// BiWeight MidVariance
	double M = get(STAT::MAD);
	double Q = get(STAT::MEDIAN);
	double n = 0; // nominator
	double d = 0; // denominator
	for (auto &Xi:data)
	{
		if (!std::isnan(Xi))
		{
			double Ui  = (Xi-Q)/(9*M);
			double Ui2 = Ui*Ui;
			double I   = (std::abs(Ui)<1 ? 1 : 0);
			n += pow(Xi-Q,2) * pow(1-Ui2,4) * I;
			d += (1-Ui2) * (1-5*Ui2) * I;
		}
	}
	if (d)	result[STAT::BWMV] = get(STAT::N) * n / (d * d);
	else	result[STAT::BWMV] = std::numeric_limits<double>::signaling_NaN();
}
template <typename T>
double Values<T>::get(STAT::StatType type)
{
//	if (data.empty() && type!=STAT::N && type!=STAT::SUM) return std::numeric_limits<double>::signaling_NaN();
	if (exist_element(result,type)) return result[type];
	switch (type) {
		case STAT::MIN:		do_min();		break;
		case STAT::MAX:		do_max();		break;
		case STAT::MEAN:	safe_mean_and_N(data,result[STAT::MEAN],result[STAT::N]);	break;
		case STAT::SUM:		sum_and_N(data,result[STAT::SUM],result[STAT::N]);			break;
		case STAT::N:		do_n();			break;
		case STAT::TSS:		result[STAT::TSS] =tss(data,get(STAT::MEAN));				break; // Total sum of squares around the mean
		case STAT::TSS0:	result[STAT::TSS0]=tss(data,0);								break; // Total sum of squares around 0
		case STAT::RMS:		do_rms();		break;
		case STAT::POP_VMR:	do_pop_vmr();	break;
		case STAT::CV:		do_cv();		break;
		case STAT::MEAN_ABS:do_mean_abs();	break;
		case STAT::MED_ABS: do_med_abs();	break;
		case STAT::MEAN_POS:do_mean_pos();	break;
		case STAT::MED_POS: do_med_pos();	break;
		case STAT::UCV:		result[STAT::UCV] =get(STAT::CV)*(1+0.25/get(STAT::N));		break; // [sample] unbiased coefficient of variation
		case STAT::RSD:		result[STAT::RSD] =std::abs(get(STAT::CV));					break; // [sample] relative standard deviation
		case STAT::URSD:	result[STAT::URSD]=std::abs(get(STAT::UCV));				break; // [sample] unbiased relative standard deviation
		case STAT::SPL_VAR:	do_spl_var();	break;
		case STAT::POP_VAR: do_pop_var();	break;
		case STAT::SPL_SD:	result[STAT::SPL_SD] = sqrt(get(STAT::SPL_VAR));			break;
		case STAT::POP_SD:	result[STAT::POP_SD] = sqrt(get(STAT::POP_VAR));			break;
		case STAT::PRD:		do_prd();		break;
		case STAT::M3:		do_m3();		break;
		case STAT::ADJ_SKEW:do_adj_skew();	break;
		case STAT::SKEW:	do_skew();		break;
		case STAT::SEM:		do_sem();		break;
		case STAT::BWMV:	do_bwmv();		break;
		case STAT::IQR:	  result[STAT::IQR]=get(STAT::PCT75)-get(STAT::PCT25);break; // interquartile range
		case STAT::QCD:   result[STAT::QCD]=(get(STAT::PCT75)-get(STAT::PCT25))/(get(STAT::PCT75)+get(STAT::PCT25));break; // quartile coefficient of dispersion
		case STAT::MAD:	   sort(); result[STAT::MAD]    = MAD(data);				break; // Median absolute deviation
		case STAT::MEDIAN: sort(); result[STAT::MEDIAN] = median(data);				break;
		case STAT::PCT001: sort(); result[STAT::PCT001] = quantile(data, 0.0001);	break;
		case STAT::PCT005: sort(); result[STAT::PCT005] = quantile(data, 0.0005);	break;
		case STAT::PCT01:  sort(); result[STAT::PCT01]  = quantile(data, 0.001);	break;
		case STAT::PCT05:  sort(); result[STAT::PCT05]  = quantile(data, 0.005);	break;
		case STAT::PCT1:   sort(); result[STAT::PCT1]   = quantile(data, 0.01);		break;
		case STAT::PCT2:   sort(); result[STAT::PCT2]   = quantile(data, 0.02);		break;
		case STAT::PCT2p5: sort(); result[STAT::PCT2p5] = quantile(data, 0.025);	break;
		case STAT::PCT5:   sort(); result[STAT::PCT5]   = quantile(data, 0.05);		break;
		case STAT::PCT9:   sort(); result[STAT::PCT9]   = quantile(data, 0.09);		break;
		case STAT::PCT10:  sort(); result[STAT::PCT10]  = quantile(data, 0.10);		break;
		case STAT::PCT25:  sort(); result[STAT::PCT25]  = quantile(data, 0.25);		break;
		case STAT::PCT75:  sort(); result[STAT::PCT75]  = quantile(data, 0.75);		break;
		case STAT::PCT90:  sort(); result[STAT::PCT90]  = quantile(data, 0.90);		break;
		case STAT::PCT91:  sort(); result[STAT::PCT91]  = quantile(data, 0.91);		break;
		case STAT::PCT95:  sort(); result[STAT::PCT95]  = quantile(data, 0.95);		break;
		case STAT::PCT97p5:sort(); result[STAT::PCT97p5]= quantile(data, 0.975);	break;
		case STAT::PCT98:  sort(); result[STAT::PCT98]  = quantile(data, 0.98);		break;
		case STAT::PCT99:  sort(); result[STAT::PCT99]  = quantile(data, 0.99);		break;
		case STAT::PCT995: sort(); result[STAT::PCT995] = quantile(data, 0.995);	break;
		case STAT::PCT999: sort(); result[STAT::PCT999] = quantile(data, 0.999);	break;
		case STAT::PCT9995:sort(); result[STAT::PCT9995]= quantile(data, 0.9995);	break;
		case STAT::PCT9999:sort(); result[STAT::PCT9999]= quantile(data, 0.9999);	break;
		case STAT::SORTED: return 0; // do not add result[SORTED], it's for internal use only.
		case STAT::ABP_MC: do_abp(); break;
		case STAT::ABP_LF: do_abp(); break;
		case STAT::ABP_UF: do_abp(); break;
		case STAT::ABP_LO: do_abp(); break;
		case STAT::ABP_UO: do_abp(); break;
		case STAT::Adil_LF: do_abp(); break;
		case STAT::Adil_UF: do_abp(); break;
		case STAT::Kimber_LF: do_kimber(); break;
		case STAT::Kimber_UF: do_kimber(); break;
		case STAT::SD_LO:	do_sd1sided(); break;
		case STAT::SD_UP:	do_sd1sided(); break;
		default:exit_error("Unknown result type "+itos(type)+" for Values.");
	}
	return result[type];
}

template class Values<char>;
template class Values<int8_t>;
template class Values<int16_t>;
template class Values<int32_t>;
template class Values<int64_t>;
template class Values<uint8_t>;
template class Values<uint16_t>;
template class Values<uint32_t>;
template class Values<uint64_t>;
template class Values<float>;
template class Values<double>;
template class Values<long double>;
#ifdef __APPLE__
template class Values<long>;
#endif

// --------------------- p_values or rankings ---------------------

// make ranks (1-based) for the right-values of a map, allow ties
template <typename T1, typename T2>
void rank_up_mapped_values(const std::map<T1,T2>& vMap, std::map<T1,double>& rMap) // v=value r=rank
{
	std::map<T2,double> vTOr;
	for (each_element(vMap,it_vMap)) ++vTOr[it_vMap->second];
	int N=0;
	for (each_element(vTOr,it_vTOr))
	{
		int n = it_vTOr->second;
		it_vTOr->second = N + (1+n)/2. ;
		N += n;
	}
	for (each_element(vMap,it_vMap)) rMap[it_vMap->first]=vTOr[it_vMap->second];
}

// make ranks (1-based) for the right-values of a map, allow ties
template <typename T1, typename T2, typename T3>
void num_ge_mapped_values(const std::map<T1,T2>& vMap, std::map<T1,T3>& rMap) // v=value r=rank
{
	std::map< T2,double,std::greater<T2> > vTOr;
	for (each_element(vMap,it_vMap)) ++vTOr[it_vMap->second];
	int N=0;
	for (each_element(vTOr,it_vTOr))
	{
		it_vTOr->second += N;
		N = it_vTOr->second;
	}
	for (each_element(vMap,it_vMap)) rMap[it_vMap->first]=vTOr[it_vMap->second];
}

// make ranks (1-based) for the right-values of a map, allow ties
template <typename Container>
void rank_up(const Container& vMap, std::vector<double>& rMap) // v=value r=rank
{
	typedef typename Container::value_type T;
	rMap.assign(vMap.size(),0);
	std::map<T,double> vTOr;
	for (auto &v:vMap) ++vTOr[v];
	int N=0;
	for (each_element(vTOr,it_vTOr))
	{
		int n = it_vTOr->second;
		it_vTOr->second = N + (1+n)/2. ;
		N += n;
	}
	for (size_t i=0;i<vMap.size();++i) rMap[i]=vTOr[vMap[i]];
}
template void rank_up(const std::vector<double>& vMap, std::vector<double>& rMap);
template void rank_up(const std::vector<int>& vMap, std::vector<double>& rMap);
template void rank_up(const std::deque<double>& vMap, std::vector<double>& rMap);
template void rank_up(const std::deque<int>& vMap, std::vector<double>& rMap);

// make ranks (1-based) for the right-values of a map, allow ties
template <typename Container>
void rank_dn(const Container& vMap, std::vector<double>& rMap) // v=value r=rank
{
	typedef typename Container::value_type T;
	rMap.assign(vMap.size(),0);
	std::map< T,double,std::greater<T> > vTOr;
	for (auto &v:vMap) ++vTOr[v];
	int N=0;
	for (each_element(vTOr,it_vTOr))
	{
		int n = it_vTOr->second;
		it_vTOr->second = N + (1+n)/2. ;
		N += n;
	}
	for (size_t i=0;i<vMap.size();++i) rMap[i]=vTOr[vMap[i]];
}
template void rank_dn(const std::vector<double>& vMap, std::vector<double>& rMap);
template void rank_dn(const std::vector<int>& vMap, std::vector<double>& rMap);
template void rank_dn(const std::deque<double>& vMap, std::vector<double>& rMap);
template void rank_dn(const std::deque<int>& vMap, std::vector<double>& rMap);

template <typename T1>
void FDR_from_p(const std::map<T1,double>& vMap, std::map<T1,double>& rMap) // v=value r=rank
{
	std::multimap<double,double> pv;
	for (each_element(vMap,it_vMap)) pv.insert(std::pair<double,double>(it_vMap->second,it_vMap->second));
	std::vector<double> Pi_m__i;
	int i=1,m=vMap.size();
	for (each_element(pv,it_pv),++i) Pi_m__i.push_back(it_pv->second * m / i);
	double minp=DBL_MAX;
	std::multimap<double,double>::reverse_iterator rit=pv.rbegin();
	for (i=m-1;i>=0;--i,++rit)
	{
		if (Pi_m__i[i]<minp) minp=Pi_m__i[i];
		rit->second=minp;
	}
	for (each_element(vMap,it_vMap)) rMap[it_vMap->first]=pv.find(it_vMap->second)->second;
}

/* Cross-phenotype meta-analysis. PLoS Genet 2011 7(8) e1002254.
 "This (method) give us more statistical power than relying on strategies combining association statistics,
 which would consume multiple degrees of freedom." -- I think the latter refers to Fisher's method.
 ## Bellow is the R code downloaded from http://www.cotsapaslab.info/index.php/software/cpma/
 cpma.score.exponential.nolog <- function(pvals,log=T,zero.val=NA) {
 ## this function tests deviation from the expected exponential behaviour of -log(p) for a set of associations to a SNP.
 ## modelled on suggestions from Chris Wallace/David Clayton, it's an implementation of a method proposed by Ben Voight.
 ## Chris Cotsapas 2011, based on previous code written 2009.
 ## This version avoids multiplication which eventually converges to zero for large p value series.
 ## Instead it transforms to log space and adds.
 
 ## internal function to compute likelihood of exponential distribution at a rate lambda
 int.exp.fn <- function(x,lambda=1,epsilon=0.001) {
 return( exp(-lambda * (x-epsilon)) - exp(-lambda * (x+epsilon)) )
 }
 
 pvals[pvals==0] <- zero.val
 borked.ix <- ( is.na(pvals) | is.infinite(pvals) )
 pvals <- pvals[!borked.ix]
 if(log) { pvals <- -log(pvals) }
 
 ## return the log likelihood ratio of -log(p) being exponentially
 ## distributed (cf biased exponential). Effectively, test if
 ## exponential decay parameter == 1
 
 p.obs <- sum(log(int.exp.fn(pvals,1/mean(pvals))))
 p.exp <- sum(log(int.exp.fn(pvals,1)))
 stat <- -2 * (p.obs-p.exp)
 return(std::fabs(stat)) ## this should be chi-sq,df=1
 }
 */

// likelihood of exponential distribution at a rate lambda
double lik_exp_dist(double x, double lambda, double epsilon)
{
	return exp(-lambda * (x-epsilon)) - exp(-lambda * (x+epsilon)) ;
}

// from a set of p-values, return a p-value
template <typename T1>
double CPMA_from_p(const std::map<T1,double>& pval)
{
	double p_obs=0, p_exp=0, m=0; // m = mean then 1/mean
	std::map<T1,double> logMap;
	for (each_element(pval,it_pval)) { logMap[it_pval->first]=-log(it_pval->second); m+=logMap[it_pval->first]; }
	m = pval.size()/m;
	for (each_element(logMap,it_pval)) {
		p_obs += log(lik_exp_dist(it_pval->second,m));
		p_exp += log(lik_exp_dist(it_pval->second,1));
	}
	return cdf_chisq_1df_q(std::fabs(-2 * (p_obs - p_exp)));
}

// http://en.wikipedia.org/wiki/Fisher's_method
// but also see http://www.burns-stat.com/pages/Working/perfmeasrandport.pdf page 14 for "why you should not use Fisher's method"
// use Stouffer's method instead.
// tested with http://www.ncbi.nlm.nih.gov/CBBresearch/Yu/downloads/CoinedPValues.html
template <typename T1>
double Fisher_from_p(const std::map<T1,double>& pval)
{
	double m=0;
	for (each_element(pval,it_pval)) m += log(it_pval->second);
	m *= -2;
	return cdf_chisq_q(std::fabs(m),2*pval.size());
}

// R code: pnorm(sum(qnorm(x)) / sqrt(length(x))) , where x is the vector of individual p-values.
// this gives the same results: 1 - pnorm(sum(qnorm(1-x)) / sqrt(length(x)))
// Stouffer, S. A., Suchman, E. A , DeVinney, L.C., Star, S.A., Williams, R.M. Jr (1949).
// Adjustment During Army Life. Princeton, NJ, Princeton University Press. (this is original reference)
template <typename T1>
double Stouffer_from_one_sided_p(const std::map<T1,double>& pval)
{
	double m=0;
	for (each_element(pval,it_pval)) m += qnorms(1-it_pval->second);
	m /= sqrt(pval.size());
	return cdf_norms_q(m);
}

// see 2011 Optimally weighted Z-test is a powerful method for combining probabilities in meta analysis.
template <typename T1>
double weighted_Stouffer_from_one_sided_p(const std::map<T1,double>& pval, std::map<T1,double>& wght)
{
	double m=0, denom=0;
	for (each_element(pval,it_pval))
	{
		m += wght[it_pval->first] * qnorms(1-it_pval->second);
		denom += wght[it_pval->first] * wght[it_pval->first];
	}
	m /= sqrt(denom);
	return cdf_norms_q(m);
}

// see 2011 Optimally weighted Z-test is a powerful method for combining probabilities in meta analysis.
// this function has not been tested by comparing with known programs
// different direction of association average out, make sure this is what you want, otherwise use the one_sided_p()
template <typename T1>
double Stouffer_from_two_sided_p(const std::map<T1,double>& stat, const std::map<T1,double>& pval)
{
	double m=0;
	for (each_element(pval,it_pval))
	{
		if (stat.find(it_pval->first)->second >= 0)
			m += qnorms(1-it_pval->second/2);
		else
			m += qnorms(  it_pval->second/2);
	}
	m /= sqrt(pval.size());
	// bellow is the same as cdf_norms_2sided_pv(m)
	double p = cdf_norms_q(m);
	if (p>.5) p = 1-p;
	return 2*p;
}

// this function has not been tested by comparing with known programs
template <typename T1>
double Stouffer_from_z(const std::map<T1,double>& zMap)
{
	double m=0;
	for (each_element(zMap,it_zMap)) m += it_zMap->second;
	m /= sqrt(zMap.size());
	return cdf_norms_2sided_pv(m);
}

void qq_plot(std::multiset<double> pv, const std::string& prefix, const int range, const std::string& type, const std::string& program, const bool save_data)
{
	if (range>11) exit_error("qq_plot() range must be <=11.");
	const int num_points=32;
	double power_log[num_points];
	power_log[0] =0.5;
	power_log[1] =0.1;
	power_log[2] =0.05;
	power_log[3] =0.02;
	power_log[4] =0.01;
	power_log[5] =0.005;
	power_log[6] =0.002;
	power_log[7] =0.001;
	power_log[8] =0.0005;
	power_log[9] =0.0002;
	power_log[10]=0.0001;
	power_log[11]=0.00005;
	power_log[12]=0.00002;
	power_log[13]=0.00001;
	power_log[14]=0.000005;
	power_log[15]=0.000002;
	power_log[16]=0.000001;
	power_log[17]=0.0000005;
	power_log[18]=0.0000002;
	power_log[19]=0.0000001;
	power_log[20]=0.00000005;
	power_log[21]=0.00000002;
	power_log[22]=0.00000001;
	power_log[23]=0.000000005;
	power_log[24]=0.000000002;
	power_log[25]=0.000000001;
	power_log[26]=0.0000000005;
	power_log[27]=0.0000000002;
	power_log[28]=0.0000000001;
	power_log[29]=0.00000000005;
	power_log[30]=0.00000000002;
	power_log[31]=0.00000000001;
	
	openOutFile_or_exit(data, prefix+"_data.txt");
	for (int j=0;j<num_points;++j)
	{
		if (power_log[j]*pv.size()<5) break;
		data << -log10( power_log[j] ) << '\t';
		data << -log10( quantile_Multiset(pv,power_log[j]) ) << '\n';
	}
	closefile(data);
	
	std::string out;
	if		(type=="pdf")			out="pdf";
	else if (type=="postscript")	out="eps";
	else exit_error("qq_plot() type must be pdf or postscript.");

	openOutFile_or_exit(cmdf, prefix+"_cmd.txt");
	cmdf<<"# usage: gnuplot "<<prefix<<"_cmd.txt"<<std::endl;
	cmdf<<"set terminal "<<type<<" size 7,7 font \"Times-Roman, 30\" linewidth 8"<<std::endl; // prv 8.5,8.5 for pdf with bmargin tmargin
	cmdf<<"set output '"<<prefix<<"."<<out<<"'"<<std::endl;
//	cmdf<<"set bmargin 2"<<std::endl;
//	cmdf<<"set tmargin 2"<<std::endl;
	cmdf<<"set multiplot"<<std::endl;
	cmdf<<"plot [0:"<<range<<"] [0:"<<range<<"] x with lines linecolor rgb 'grey' tit ''"<<std::endl;
	cmdf<<"plot [0:"<<range<<"] [0:"<<range<<"] '"<<prefix<<"_data.txt' with lines linecolor rgb 'blue' tit ''"<<std::endl;
	closefile(cmdf);
	
	system((program+" "+prefix+"_cmd.txt").c_str());
	if (!save_data)
	{
		system(("rm -f "+prefix+"_cmd.txt").c_str());
		system(("rm -f "+prefix+"_data.txt").c_str());
	}
}

// PO = p/(1-p) = BF * pi/(1-pi) = BF'  ==> p = BF'/(1+BF')
// with an uninformative prior (pi=0.5), pi/(1-pi)=1, so BF'=BF
double posterior_given_log10BF(double log10_BF)
{
	if (std::isnan(log10_BF)) return std::numeric_limits<double>::signaling_NaN();
	double BF = pow(10,log10_BF);
	return BF/(1+BF);
}

double posterior_given_log10BF(double log10_BF, double pi)
{
	if (std::isnan(log10_BF)) return std::numeric_limits<double>::signaling_NaN();
	if (std::isnan(pi)) return std::numeric_limits<double>::signaling_NaN();
	if (pi<0 || pi>1) exit_error("prior probability cannot be <0 or >1");
	if (pi==0) return 0;
	if (pi==1) return 1;
	double BFp = pow(10,log10_BF) * pi;
	return BFp/(BFp+1-pi);
}

// --------------------- outlier detection ---------------------

template <typename T>
double cal_outlier_byDQ(const outlier_info& ol_inf, const T& val) // return p-val
{
	if (val <= ol_inf.l_thr) return ol_inf.dq_pl;
	if (val >= ol_inf.r_thr) return ol_inf.dq_pr;
	return 1;
}

template <typename T>
void num_outlier_byDQ(const std::multiset<T>& values, outlier_info& ol_inf, double alpha)
//Dixon's Q test
//http://www.statistics4u.info/fundstat_eng/cc_outlier_tests_dixon.html#
{
	const static double crtval[101][8]={
		{ 0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000 },
		{ 0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000 },
		{ 0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000 },
		{ 0.999,0.998,0.994,0.988,0.976,0.941,0.886,0.782 },
		{ 0.964,0.949,0.921,0.889,0.847,0.766,0.679,0.561 },
		{ 0.895,0.869,0.824,0.782,0.729,0.643,0.559,0.452 },
		{ 0.822,0.792,0.744,0.698,0.646,0.563,0.484,0.387 },
		{ 0.763,0.731,0.681,0.636,0.587,0.507,0.433,0.344 },
		{ 0.799,0.769,0.724,0.682,0.633,0.554,0.480,0.386 },
		{ 0.750,0.720,0.675,0.634,0.586,0.512,0.441,0.352 },
		{ 0.713,0.683,0.637,0.597,0.551,0.477,0.409,0.325 },
		{ 0.770,0.746,0.708,0.674,0.636,0.575,0.518,0.445 },
		{ 0.739,0.714,0.676,0.643,0.605,0.546,0.489,0.420 },
		{ 0.713,0.687,0.649,0.617,0.580,0.522,0.467,0.399 },
		{ 0.732,0.708,0.672,0.640,0.603,0.546,0.491,0.422 },
		{ 0.708,0.685,0.648,0.617,0.582,0.524,0.470,0.403 },
		{ 0.691,0.667,0.630,0.598,0.562,0.505,0.453,0.386 },
		{ 0.671,0.647,0.611,0.580,0.545,0.489,0.437,0.373 },
		{ 0.652,0.628,0.594,0.564,0.529,0.475,0.424,0.361 },
		{ 0.640,0.617,0.581,0.551,0.517,0.462,0.412,0.349 },
		{ 0.627,0.604,0.568,0.538,0.503,0.450,0.401,0.339 },
		{ 0.616,0.593,0.558,0.528,0.494,0.441,0.393,0.332 },
		{ 0.606,0.582,0.548,0.518,0.485,0.432,0.384,0.324 },
		{ 0.595,0.572,0.537,0.509,0.475,0.424,0.376,0.317 },
		{ 0.585,0.561,0.527,0.499,0.466,0.415,0.367,0.309 },
		{ 0.574,0.550,0.517,0.489,0.457,0.406,0.359,0.302 },
		{ 0.567,0.543,0.510,0.482,0.451,0.400,0.354,0.297 },
		{ 0.560,0.537,0.504,0.476,0.444,0.394,0.348,0.292 },
		{ 0.553,0.530,0.497,0.469,0.438,0.388,0.343,0.288 },
		{ 0.546,0.524,0.491,0.463,0.431,0.382,0.337,0.283 },
		{ 0.539,0.517,0.484,0.456,0.425,0.376,0.332,0.278 },
		{ 0.533,0.512,0.479,0.451,0.420,0.372,0.328,0.274 },
		{ 0.528,0.506,0.474,0.446,0.415,0.367,0.324,0.271 },
		{ 0.522,0.501,0.469,0.441,0.410,0.363,0.319,0.267 },
		{ 0.517,0.495,0.464,0.436,0.405,0.358,0.315,0.264 },
		{ 0.511,0.490,0.459,0.431,0.400,0.354,0.311,0.260 },
		{ 0.507,0.486,0.455,0.427,0.396,0.351,0.308,0.257 },
		{ 0.503,0.482,0.451,0.423,0.393,0.347,0.305,0.254 },
		{ 0.498,0.477,0.446,0.420,0.389,0.344,0.301,0.252 },
		{ 0.494,0.473,0.442,0.416,0.386,0.340,0.298,0.249 },
		{ 0.490,0.469,0.438,0.412,0.382,0.337,0.295,0.246 },
		{ 0.487,0.466,0.435,0.409,0.379,0.334,0.293,0.244 },
		{ 0.484,0.463,0.432,0.406,0.376,0.331,0.290,0.241 },
		{ 0.481,0.460,0.429,0.403,0.374,0.329,0.288,0.239 },
		{ 0.478,0.457,0.426,0.400,0.371,0.326,0.285,0.236 },
		{ 0.475,0.454,0.423,0.397,0.368,0.323,0.283,0.234 },
		{ 0.472,0.451,0.420,0.394,0.365,0.321,0.281,0.232 },
		{ 0.469,0.448,0.418,0.392,0.363,0.319,0.279,0.231 },
		{ 0.466,0.445,0.415,0.389,0.360,0.316,0.276,0.229 },
		{ 0.463,0.442,0.413,0.387,0.358,0.314,0.274,0.228 },
		{ 0.460,0.439,0.410,0.384,0.355,0.312,0.272,0.226 },
		{ 0.458,0.437,0.408,0.382,0.353,0.310,0.270,0.225 },
		{ 0.455,0.435,0.406,0.380,0.351,0.308,0.269,0.223 },
		{ 0.453,0.432,0.403,0.378,0.349,0.307,0.267,0.222 },
		{ 0.451,0.430,0.401,0.376,0.347,0.305,0.266,0.220 },
		{ 0.449,0.428,0.399,0.374,0.346,0.303,0.264,0.219 },
		{ 0.446,0.426,0.397,0.371,0.344,0.301,0.262,0.217 },
		{ 0.444,0.424,0.395,0.369,0.342,0.299,0.261,0.216 },
		{ 0.442,0.421,0.392,0.367,0.340,0.298,0.259,0.214 },
		{ 0.439,0.419,0.390,0.365,0.338,0.296,0.258,0.213 },
		{ 0.437,0.417,0.388,0.363,0.336,0.294,0.256,0.211 },
		{ 0.436,0.416,0.387,0.362,0.335,0.293,0.255,0.210 },
		{ 0.434,0.414,0.385,0.360,0.333,0.291,0.254,0.209 },
		{ 0.433,0.413,0.384,0.359,0.332,0.290,0.252,0.208 },
		{ 0.431,0.411,0.382,0.357,0.330,0.288,0.251,0.207 },
		{ 0.430,0.410,0.381,0.356,0.329,0.287,0.250,0.206 },
		{ 0.428,0.409,0.380,0.355,0.327,0.286,0.249,0.205 },
		{ 0.427,0.407,0.378,0.353,0.326,0.284,0.248,0.204 },
		{ 0.425,0.406,0.377,0.352,0.324,0.283,0.246,0.203 },
		{ 0.424,0.404,0.375,0.350,0.323,0.281,0.245,0.202 },
		{ 0.422,0.403,0.374,0.349,0.321,0.280,0.244,0.201 },
		{ 0.421,0.402,0.373,0.348,0.320,0.279,0.243,0.200 },
		{ 0.419,0.400,0.371,0.347,0.319,0.278,0.242,0.199 },
		{ 0.418,0.399,0.370,0.345,0.318,0.277,0.241,0.198 },
		{ 0.416,0.397,0.368,0.344,0.317,0.276,0.240,0.197 },
		{ 0.415,0.396,0.367,0.343,0.316,0.275,0.239,0.197 },
		{ 0.414,0.395,0.366,0.342,0.314,0.274,0.238,0.196 },
		{ 0.412,0.393,0.364,0.341,0.313,0.273,0.237,0.195 },
		{ 0.411,0.392,0.363,0.339,0.312,0.272,0.236,0.194 },
		{ 0.409,0.390,0.361,0.338,0.311,0.271,0.235,0.193 },
		{ 0.408,0.389,0.360,0.337,0.310,0.270,0.234,0.192 },
		{ 0.407,0.388,0.359,0.336,0.309,0.269,0.233,0.191 },
		{ 0.406,0.387,0.358,0.335,0.308,0.268,0.232,0.191 },
		{ 0.405,0.385,0.357,0.334,0.307,0.267,0.232,0.190 },
		{ 0.404,0.384,0.356,0.333,0.306,0.266,0.231,0.189 },
		{ 0.403,0.383,0.355,0.332,0.305,0.266,0.230,0.189 },
		{ 0.401,0.382,0.354,0.330,0.304,0.265,0.229,0.188 },
		{ 0.400,0.381,0.353,0.329,0.303,0.264,0.228,0.187 },
		{ 0.399,0.379,0.352,0.328,0.302,0.263,0.228,0.186 },
		{ 0.398,0.378,0.351,0.327,0.301,0.262,0.227,0.186 },
		{ 0.397,0.377,0.350,0.326,0.300,0.261,0.226,0.185 },
		{ 0.396,0.376,0.349,0.325,0.299,0.260,0.225,0.184 },
		{ 0.395,0.375,0.348,0.324,0.298,0.259,0.225,0.184 },
		{ 0.394,0.374,0.347,0.323,0.298,0.259,0.224,0.183 },
		{ 0.393,0.373,0.346,0.322,0.297,0.258,0.223,0.183 },
		{ 0.392,0.373,0.346,0.322,0.296,0.257,0.223,0.182 },
		{ 0.391,0.372,0.345,0.321,0.295,0.256,0.222,0.181 },
		{ 0.390,0.371,0.344,0.320,0.294,0.255,0.221,0.181 },
		{ 0.389,0.370,0.343,0.319,0.294,0.255,0.220,0.180 },
		{ 0.388,0.369,0.342,0.318,0.293,0.254,0.220,0.180 },
		{ 0.387,0.368,0.341,0.317,0.292,0.253,0.219,0.179 }
	};
	const static double cvPval[9]=
	{ 0.001,0.002,0.005,0.01, 0.020,0.050,0.100,0.200,1.000 };
	ol_inf.init();
	int N=values.size();
	if (N>100 || N<3) return;
	typename std::multiset<T>::iterator fi=values.begin(); // forward iterator
	double minv = *fi; ++fi;
	double min2 = *fi; ++fi;
	double min3 = *fi;
	typename std::multiset<T>::reverse_iterator ri=values.rbegin(); // reverse iterator
	double maxv = *ri; ++ri;
	double max2 = *ri; ++ri;
	double max3 = *ri;
	ol_inf.success=true;
	double Ql,Qr;
	if (N<8)
	{
		Ql = (min2 - minv) / (maxv - minv);
		Qr = (maxv - max2) / (maxv - minv);
	}
	else if (N<11)
	{
		Ql = (min2 - minv) / (max2 - minv);
		Qr = (maxv - max2) / (maxv - min2);
	}
	else if (N<14)
	{
		Ql = (min3 - minv) / (max2 - minv);
		Qr = (maxv - max3) / (maxv - min2);
	}
	else
	{
		Ql = (min3 - minv) / (max3 - minv);
		Qr = (maxv - max3) / (maxv - min3);
	}
	int i;
	for (i=0;i<8;++i) if (Ql > crtval[N][i]) break;
	if (cvPval[i] <= alpha) { ++ol_inf.num_l_ol; ol_inf.l_thr=minv; ol_inf.dq_pl=cvPval[i]; }
	for (i=0;i<8;++i) if (Qr > crtval[N][i]) break;
	if (cvPval[i] <= alpha) { ++ol_inf.num_r_ol; ol_inf.r_thr=maxv; ol_inf.dq_pr=cvPval[i]; }
}

template <typename T>
double cal_outlier_byBP(const outlier_info& ol_inf, const T& val)
{
	if (val < ol_inf.q25pc) return (val - ol_inf.q25pc)/ol_inf.IQR;
	if (val > ol_inf.q75pc) return (val - ol_inf.q75pc)/ol_inf.IQR;
	return 0;
}

template <typename T>
void num_outlier_byBP(const std::multiset<T>& values, outlier_info& ol_inf, double k)
// box plot
{
	ol_inf.init();
	int N=values.size();
	if (N<4) return;
	ol_inf.q25pc=quantile_Multiset(values,0.25);
	ol_inf.q75pc=quantile_Multiset(values,0.75);
	ol_inf.IQR = ol_inf.q75pc - ol_inf.q25pc;
	ol_inf.l_thr=ol_inf.q25pc - k * ol_inf.IQR;
	ol_inf.r_thr=ol_inf.q75pc + k * ol_inf.IQR;
	ol_inf.success=true;
	for (typename std::multiset<T>::iterator it=values.begin(); it!=values.end(); ++it)
	{
		double fr = cal_outlier_byBP(ol_inf,*it);
		if (std::fabs(fr) > k) ++ol_inf.num_l_ol;
		else break;
	}
	for (typename std::multiset<T>::reverse_iterator it=values.rbegin(); it!=values.rend(); ++it)
	{
		double fr = cal_outlier_byBP(ol_inf,*it);
		if (std::fabs(fr) > k) ++ol_inf.num_r_ol;
		else break;
	}
}

template <typename T>
double cal_outlier_byIH(const outlier_info& ol_inf, const T& val) // return Mi
{
	return 0.6745 * ( val - ol_inf.median ) / ol_inf.mad;
}

template <typename T>
void num_outlier_byIH(const std::multiset<T>& values, outlier_info& ol_inf, double alpha) // median absolute deviation
// http://itl.nist.gov/div898/handbook/eda/section3/eda35h.htm
//+-------Iglewicz and Hoaglin's modified Z-sore to detect outliers-----+
//| modified Z-score Mi = 0.6745(Yi-Y~)/MAD                             |
//| where, median absolute deviation (MAD) = median(|Yi-Y~|)            |
//| Y~ = median of the data                                             |
//| data points whose Mi > 3.5 should be labeled as potential outliers, |
//| which has an alpha of 0.000465258 in a two-sided test.              |
//+---------------------------------------------------------------------+
{
	ol_inf.init();
	int N=values.size();
	if (N<1) return;
	ol_inf.ih_sigval=qnorms(1-alpha/2); // same as static const boost::math::normal dist; ol_inf.ih_sigval=boost::math::quantile(dist,1-alpha/2);
	ol_inf.mad=MAD_Multiset(values);
	ol_inf.median=median_Multiset(values);
	ol_inf.l_thr=ol_inf.median-ol_inf.ih_sigval*ol_inf.mad/0.6745;
	ol_inf.r_thr=ol_inf.median+ol_inf.ih_sigval*ol_inf.mad/0.6745;
	ol_inf.success=true;
	for (typename std::multiset<T>::iterator it=values.begin(); it!=values.end(); ++it)
	{
		double Mi = cal_outlier_byIH(ol_inf,*it);
		if (std::fabs(Mi) > ol_inf.ih_sigval) ++ol_inf.num_l_ol;
		else break;
	}
	for (typename std::multiset<T>::reverse_iterator it=values.rbegin(); it!=values.rend(); ++it)
	{
		double Mi = cal_outlier_byIH(ol_inf,*it);
		if (std::fabs(Mi) > ol_inf.ih_sigval) ++ol_inf.num_r_ol;
		else break;
	}
}

template <typename T>
double cal_outlier_byGr(const outlier_info& ol_inf, const T& val) // return G
{
	return (val - ol_inf.mean) / ol_inf.stddev;
}

template <typename T>
void num_outlier_byGr(const std::multiset<T>& values, outlier_info& ol_inf, double alpha) // took this number from above
{
	ol_inf.init();
	int N=values.size();
	if (N<3) return;
	ol_inf.stddev=sqrt(sample_variance(values));
	ol_inf.mean=sum(values)/values.size();
	ol_inf.gr_t2=qstudents_t(N-2,1-alpha/N); // same as: boost::math::students_t dist(N-2); ol_inf.gr_t2=quantile(dist,1-alpha/N);
	ol_inf.gr_t2*=ol_inf.gr_t2;
	ol_inf.gr_sigval=(N-1)/sqrt(N)*sqrt(ol_inf.gr_t2/(N-2+ol_inf.gr_t2));
	ol_inf.l_thr = ol_inf.mean - ol_inf.gr_sigval * ol_inf.stddev;
	ol_inf.r_thr = ol_inf.mean + ol_inf.gr_sigval * ol_inf.stddev;
	ol_inf.success=true;
	if ( std::fabs(cal_outlier_byGr(ol_inf,*(values. begin()))) > ol_inf.gr_sigval)	++ol_inf.num_l_ol;
	if ( std::fabs(cal_outlier_byGr(ol_inf,*(values.rbegin()))) > ol_inf.gr_sigval)	++ol_inf.num_r_ol;
}

template double cal_outlier_byDQ(const outlier_info& ol_inf, const double& val);
template double cal_outlier_byBP(const outlier_info& ol_inf, const double& val);
template double cal_outlier_byIH(const outlier_info& ol_inf, const double& val);
template double cal_outlier_byGr(const outlier_info& ol_inf, const double& val);
template void   num_outlier_byDQ(const std::multiset<double>& values, outlier_info& ol_inf, double alpha);
template void   num_outlier_byBP(const std::multiset<double>& values, outlier_info& ol_inf, double k);
template void   num_outlier_byIH(const std::multiset<double>& values, outlier_info& ol_inf, double alpha);
template void   num_outlier_byGr(const std::multiset<double>& values, outlier_info& ol_inf, double alpha);

Interpolate::Interpolate(std::initializer_list< std::pair<double,double> > args)
{
	for (auto iter = args.begin(); iter != args.end(); ++iter)
	{
		const double& x=iter->first;
		const double& y=iter->second;
		if (std::isnan(x) || std::isnan(y) || std::isinf(x) || std::isinf(y)) exit_error("Interpolate::setup cannot use inf/nan.");
		reference_data.insert(*iter);
	}
}

void Interpolate::setup(std::istream& file) {
	reference_data.clear();
	for (each_line(file))
	{
		double x,y;
		file >> x >> y;
		if (std::isnan(x) || std::isnan(y) || std::isinf(x) || std::isinf(y)) exit_error("Interpolate::setup cannot use inf/nan.");
		reference_data[x]=y;
	}
}

void Interpolate::setup(const std::string& filename) {
	reference_data.clear();
	openInpFile_or_exit(file,filename);
	setup(file);
	closefile(file);
}

double Interpolate::solve(double x, int mode) {
	if (reference_data.empty()) return std::numeric_limits<double>::signaling_NaN();
	if (std::isnan(x))			return std::numeric_limits<double>::signaling_NaN();
	std::map<double,double>::iterator it1 = reference_data.lower_bound(x);
	std::map<double,double>::iterator it2 = reference_data.upper_bound(x);
	if (it1==reference_data.end())
	{
		if (mode==0) return std::numeric_limits<double>::signaling_NaN();
		if (mode==1) return reference_data.rbegin()->second;
		--it1;
		--it2;
	}
	if (x < reference_data.begin()->first)
	{
		if (mode==0) return std::numeric_limits<double>::signaling_NaN();
		if (mode==1) return reference_data. begin()->second;
		++it1;
		++it2;
	}
	if (it1==it2) --it1;
	const double& x1 = it1->first;
	const double& y1 = it1->second;
	const double& x2 = it2->first;
	const double& y2 = it2->second;
	return (y2-y1)*(x-x1)/(x2-x1) + y1;
}
