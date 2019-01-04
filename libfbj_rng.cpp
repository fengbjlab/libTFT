#include <fstream>
#include <sstream>
#include <iostream>
#include <set>
#include <cstdio>
#include <algorithm>
#include <numeric>
#include <boost/date_time/local_time/local_time.hpp>
#include <boost/cstdint.hpp>
#include <unistd.h>	// need for getpid() function
#include <time.h>	// need for time() function
#include "libfbj_rng.hpp"
#include "libfbj_base.hpp"

#ifdef _CHECK_LICENSE
	#include "libfbj_program.hpp" // for license checking
	int rng_by_boost::lc = program.nudge;
#else
	int rng_by_boost::lc = 0;
#endif

rng_by_boost rng;

using namespace std;

// -------------------- class rng_by_boost --------------------

namespace {
	// http://stackoverflow.com/questions/322938/recommended-way-to-initialize-srand
	// http://www.concentric.net/~Ttwang/tech/inthash.htm
	inline uint32_t mix(uint32_t a, uint32_t b, uint32_t c)
	{
		a=a-b;  a=a-c;  a=a^(c >> 13);
		b=b-c;  b=b-a;  b=b^(a << 8);
		c=c-a;  c=c-b;  c=c^(b >> 13);
		a=a-b;  a=a-c;  a=a^(c >> 12);
		b=b-c;  b=b-a;  b=b^(a << 16);
		c=c-a;  c=c-b;  c=c^(b >> 5);
		a=a-b;  a=a-c;  a=a^(c >> 3);
		b=b-c;  b=b-a;  b=b^(a << 10);
		c=c-a;  c=c-b;  c=c^(b >> 15);
		return c;
	}
	
	inline uint32_t microsec()
	{
		boost::posix_time::ptime time_set_seed=boost::posix_time::microsec_clock::local_time();
		std::stringstream s;
		std::string temp;
		char c;
		int n;
		uint32_t microsec;
		s<<time_set_seed;
		s>>temp>>n>>c>>n>>c>>n>>c>>microsec;
		return microsec;
	}
	
	// http://stackoverflow.com/questions/2572366/how-to-use-dev-random-or-urandom-in-c
	inline uint32_t _dev_random()
	{
		ifstream file;
		file.open("/dev/random", ios::binary);
		uint32_t myRandomInteger;
		file.read((char *)&myRandomInteger, sizeof(uint32_t));
		file.close();
		return myRandomInteger;
	}
}

rng_by_boost::rng_by_boost(const std::string& f):int_one_to_max(1,INT_MAX),dice_int(*this,int_one_to_max),zeroone(rng)
{ 
	time_begin=boost::posix_time::microsec_clock::local_time();
	set_seed(f);
}

rng_by_boost::rng_by_boost():int_one_to_max(1,INT_MAX),dice_int(*this,int_one_to_max),zeroone(rng)
{
	time_begin=boost::posix_time::microsec_clock::local_time();
	set_seed();
}

rng_by_boost::rng_by_boost(unsigned s):int_one_to_max(1,INT_MAX),dice_int(*this,int_one_to_max),zeroone(rng)
{ 
	time_begin=boost::posix_time::microsec_clock::local_time();
	set_seed(s);
}

rng_by_boost::rng_by_boost(const rng_by_boost& orig):int_one_to_max(1,INT_MAX),dice_int(*this,int_one_to_max),zeroone(rng)
{
	exit_error("rng_by_boost copy constructor is forbidden because: even set with the same seed, output is diff if the other rng_by_boost is already used.");
	seed_filename=orig.seed_filename;
	time_begin=orig.time_begin;
	time_end=orig.time_end;
	set_seed(orig.readseed); 
}

rng_by_boost& rng_by_boost::operator=(const rng_by_boost& orig) 
{ 
	exit_error("rng_by_boost copy operator is forbidden because: even set with the same seed, output is diff if the other rng_by_boost is already used.");
	if (&orig!=this) {
		seed_filename=orig.seed_filename;
		time_begin=orig.time_begin;
		time_end=orig.time_end;
		set_seed(orig.readseed); }
	return *this; 
}

void rng_by_boost::set_seed(const std::string& f)
{
	seed_filename=f;
	std::fstream seedfile;
	seedfile.open(seed_filename.c_str(),std::ios::in);
	if (!seedfile.is_open()) { exit_error("Cannot open seed file "+seed_filename); }
	seedfile>>readseed;
	seedfile.close();
	seed(readseed);
	// srand (readseed); removed since it's out of the scope of this class, and it's not thread-safe
}

// seed is hash value of microsec, sec since 01/01/1970, process identifier, [and /dev/random if it's a Mac]
void rng_by_boost::set_seed()
{
// unlike Mac, in linux /dev/random could be very slow, especially for the second time it's called
#if defined(__APPLE__) || defined(__MACH__)
	uint32_t tmp_seed	= mix(368800899, getpid(), _dev_random());
	readseed			= mix(microsec(), time(NULL), tmp_seed);
#else
	readseed			= mix(microsec(), time(NULL), getpid());
#endif
	seed(readseed);
	// srand (readseed); removed since it's out of the scope of this class, and it's not thread-safe
}

void rng_by_boost::set_seed(unsigned s)
{
	readseed=s;
	seed(s);
	// srand (readseed); removed since it's out of the scope of this class, and it's not thread-safe
}

rng_by_boost::~rng_by_boost() 
{
	if (!seed_filename.empty()) //write seed for next time 
	{	
		unsigned writeseed=dice_int();
		std::fstream seedfile;
		seedfile.open(seed_filename.c_str(),std::ios::out);
		if (!seedfile.is_open()) { exit_error("Cannot write seed to "+seed_filename); }
		seedfile<<writeseed<<" <= This is the seed for next usage.\nLast seed is : "<<readseed<<endl
		<<"RNG instantiation: "<<time_begin<<endl;
		time_end=boost::posix_time::microsec_clock::local_time();
		boost::posix_time::time_duration d = time_end-time_begin;
		seedfile<<"Seed recreated at: "<<time_end<<endl
		<<"Duration         :             "<<d
		<<"\n===========================================================================\n";
		seedfile.close();
	}
}

unsigned long rng_by_boost::uniform_int_ge0leMax()				{ return dice_int();	}	// return [0,INT_MAX]
unsigned long rng_by_boost::uniform_int_ge0ltN(unsigned long n)	{ return dice_int()%n;	}	// return [0,n)
double rng_by_boost::uniform_gt0lt1()	{ for (double r=zeroone();r;r=zeroone()) return r; exit_error("uniform_gt0lt1 error"); return -1; }
double rng_by_boost::uniform_ge0lt1()	{ return zeroone(); }
int rng_by_boost::flip(const double& p) { return uniform_ge0lt1()<p; }

void rng_by_boost::memran(unsigned char* buffer, size_t length) {
	unsigned char* ptr=buffer;
	for (size_t i=0;i<length;++i,++ptr)
		*ptr = uniform_int_ge0leMax();
}

// -------------------- class rng_normal_distr --------------------

rng_normal_distr::~rng_normal_distr() { if (dist) { delete dist; delete vgen; } }
rng_normal_distr::rng_normal_distr():dist(NULL),vgen(NULL) {} // default = no allocation
rng_normal_distr& rng_normal_distr::operator=(const rng_normal_distr& othr) {
	if (this==&othr) return *this;
	if (dist) { delete dist; delete vgen; }
	dist = new boost::normal_distribution<> (*(othr.dist));
	vgen = new boost::variate_generator<RNGType&, boost::normal_distribution<> > (*(othr.vgen));
	return *this;
}
rng_normal_distr::rng_normal_distr(const rng_normal_distr& othr) {
	dist = new boost::normal_distribution<> (*(othr.dist));
	vgen = new boost::variate_generator<RNGType&, boost::normal_distribution<> > (*(othr.vgen));
}
rng_normal_distr::rng_normal_distr(double m,double s) {
	dist=new boost::normal_distribution<>(m,s);
	vgen=new boost::variate_generator<RNGType&, boost::normal_distribution<> >(rng,*dist);
}
void rng_normal_distr::setup_par(double m,double s){
	if (dist) { delete dist; delete vgen; }
	dist=new boost::normal_distribution<>(m,s);
	vgen=new boost::variate_generator<RNGType&, boost::normal_distribution<> >(rng,*dist);
}
double rng_normal_distr::gen_num() {
	if (dist) return (*vgen)();
	exit_error("Normal distribution not defined; can't generate random number ~ N(mean,sigma).");
	return 0;
}

// -------------------- class uniform_int_ge0ltN --------------------

uniform_int_ge0ltN::~uniform_int_ge0ltN() { if (dist) { delete dist; delete vgen; } }
uniform_int_ge0ltN::uniform_int_ge0ltN():dist(NULL),vgen(NULL) {}
uniform_int_ge0ltN& uniform_int_ge0ltN::operator=(const uniform_int_ge0ltN& othr) {
	if (this==&othr) return *this;
	if (dist) { delete dist; delete vgen; }
	dist = new boost::uniform_int<int> (*(othr.dist));
	vgen = new boost::variate_generator< RNGType&, boost::uniform_int<int> > (*(othr.vgen));
	return *this;
}
uniform_int_ge0ltN::uniform_int_ge0ltN(const uniform_int_ge0ltN& othr) {
	dist = new boost::uniform_int<int> (*(othr.dist));
	vgen = new boost::variate_generator< RNGType&, boost::uniform_int<int> > (*(othr.vgen));
}
uniform_int_ge0ltN::uniform_int_ge0ltN( const int N ):dist(NULL),vgen(NULL) {
	dist = new boost::uniform_int<int> (0, N-1);
	vgen = new boost::variate_generator< RNGType&, boost::uniform_int<int> > (rng,*dist);
}
void uniform_int_ge0ltN::setup_par(int N) {
	if (dist) { delete dist; delete vgen; }
	dist = new boost::uniform_int<int> (0, N-1);
	vgen = new boost::variate_generator< RNGType&, boost::uniform_int<int> > (rng,*dist);
}
int uniform_int_ge0ltN::gen_num() {
	if (dist) return (*vgen)(); 
	exit_error("N not defined; can't generate random integer in [0,N).");
	return 0;
}
std::ptrdiff_t uniform_int_ge0ltN::operator()( std::ptrdiff_t arg ) {
	if (dist) return static_cast<std::ptrdiff_t>((*vgen)());
	exit_error("N not defined; can't generate random integer in [0,N).");
	return 0;
}

// -------------------- class wheeling_select --------------------

wheeling_select::~wheeling_select(){ if (dist) { delete dist; delete vgen; } }
wheeling_select::wheeling_select():dist(NULL),vgen(NULL){}
wheeling_select& wheeling_select::operator=(const wheeling_select& othr) {
	if (this==&othr) return *this;
	if (dist) { delete dist; delete vgen; }
	dist = new boost::uniform_real<> (*(othr.dist));
	vgen = new boost::variate_generator<RNGType&, boost::uniform_real<> > (*(othr.vgen));
	cumulative=othr.cumulative;
	return *this;
}
wheeling_select::wheeling_select(const wheeling_select& othr) {
	dist = new boost::uniform_real<> (*(othr.dist));
	vgen = new boost::variate_generator<RNGType&, boost::uniform_real<> > (*(othr.vgen));
	cumulative=othr.cumulative;
}
wheeling_select::wheeling_select(const double *f,const int& n)
{	std::partial_sum(f, f + n, std::back_inserter(cumulative));
	dist = new boost::uniform_real<> (0,cumulative.back());
	vgen = new boost::variate_generator<RNGType&, boost::uniform_real<> > (rng,*dist);
}
wheeling_select::wheeling_select(const vector<double>& v)
{	std::partial_sum(v.begin(), v.end(), std::back_inserter(cumulative)); 
	dist = new boost::uniform_real<> (0,cumulative.back());
	vgen = new boost::variate_generator<RNGType&, boost::uniform_real<> > (rng,*dist);
}
void wheeling_select::setup_par(const double *f,const int& n)
{	cumulative.clear();
	std::partial_sum(f, f + n, std::back_inserter(cumulative));
	if (dist) { delete dist; delete vgen; }
	dist = new boost::uniform_real<> (0,cumulative.back());
	vgen = new boost::variate_generator<RNGType&, boost::uniform_real<> > (rng,*dist);
}
void wheeling_select::setup_par(const vector<double>& v)
{	cumulative.clear();
	std::partial_sum(v.begin(), v.end(), std::back_inserter(cumulative)); 
	if (dist) { delete dist; delete vgen; }
	dist = new boost::uniform_real<> (0,cumulative.back());
	vgen = new boost::variate_generator<RNGType&, boost::uniform_real<> > (rng,*dist);
}
int wheeling_select::gen_num() // return [0,n-1]
{	if (dist) return (std::lower_bound(cumulative.begin(), cumulative.end(), (*vgen)()) - cumulative.begin());
	exit_error("Probabilties not defined; can't generate random-integers-with-probabilities.");
	return 0;
}

// -------------------- individual functions --------------------

std::string random_string(int length)
{
	static const char alphanum[] =
	"0123456789"
	"ABCDEFGHIJKLMNOPQRSTUVWXYZ"
	"abcdefghijklmnopqrstuvwxyz";
	
	string result;
	for (int i = 0; i < length; ++i) result.push_back( alphanum[MyRandom(62)] );
	return result;
}
