#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <set>
#include "libfbj_base.hpp"
#include "libfbj_fn.hpp"
#include "libfbj_math.hpp"
#include "libfbj_program.hpp" // for license checking

/*
Each time modify rng_list, remember to: 1) cal mxneg; 2) _clear_vec 3) _setup_size !!!
INPUT: positive numbers are substract by 1, negative numbers and zero remain the same
TEST: size_valid (automatically called by tabular_file) do the following
  > check whether all negative numbers can be solved
  > update cur_size
  > check whether all REQUIRED fields exist
    violation of any one of above will return 'false', otherwise true
To get field_numbers
 1) contents_to_xxx
  > check _solvable, if not, do nothing and return false
  > If one field not exist, no insertion or insert default value. 
  > if AllRequired then throw an exception.
 2) vector<int>::const_iterator
  > if has negative field number but not yet set_size, exit program
  > iterate the field number vector (blind to the content vector, therefore, AllRequired doesn't matter here.)
  > !! It's not safe to use this method if you are not sure all fields exist !!
    It always return a value i, it would go out of boundary if inf.contents()[i]
*/

using namespace std;

struct field_numbers::field_numbers_data {
	static int				license_chk;// licensing of the class
	bool					NewEachTime;// modified only in constructor,
	bool					AllRequired;// modified only in constructor,
	bool					sizeNset;	// modified only in setup_size, current size is not set
	int						min_size;	// modified only in setup_size, minimum size required if AllRequired
	int						cur_size;	// modified only in setup_size, current size 
	int						mxneg;		// modified only in parse & push_back
	vector<pair<int,int> >	rng_list;	// modified only in parse & push_back
	vector<int>				rng_sect;	// modified only in parse & push_back
	vector<int>				fn_vec;		// created field number vector. allow dup and retain seq
	vector<int>				fn_sec;		// section
	void _parse(istream& ss);			// INPUT data, format like: 0--3,2,-5-; adj=adjustment if # >=0, normally -1
	bool _setup_size(int s);			// SETUP size, s should be .size(). return whether there's unsolved negative #
	void _clear_vec();					// clear fn_vec
	field_numbers_data():NewEachTime(false),AllRequired(true),sizeNset(true),min_size(0),cur_size(-1),mxneg(-1){}
};

int field_numbers::field_numbers_data::license_chk = program.nudge;

void field_numbers::field_numbers_data::_parse(istream& ss) // format like: 0--3,2,-5-; adj=adjustment if # >=0, normally -1
{
	//clear(); // this line removed so that parse/push_back can be called aggregately
	if (ss.peek()=='N') { ss.get(); ss.get(); NewEachTime=true;  } // potential input is N:<integer_ranges>
	if (ss.peek()=='A') { ss.get(); ss.get(); AllRequired=true;  } // or N,A:<integer_ranges>
	if (ss.peek()=='n') { ss.get(); ss.get(); NewEachTime=false; }
	if (ss.peek()=='a') { ss.get(); ss.get(); AllRequired=false; }
	for (char c=ss.peek(); (c>='0' && c<='9') || c=='-'; c=ss.peek())
	{
		pair<int,int> r;
		ss>>r.first;
		switch (ss.peek())
		{	case ',':	r.second=r.first; rng_list.push_back(r); ss.get(); continue;
			case '-':	ss.get();
				c=ss.peek();
				if ((c>='0' && c<='9') || c=='-') {	ss>>r.second; c=ss.peek(); }
				else								r.second=-1;
				rng_list.push_back(r);
				if (c==',') { ss.get(); continue; }
				break;
			default:	r.second=r.first; rng_list.push_back(r);
		}
		break;
	}
	rng_sect.push_back(rng_list.size());
	for (each_element(rng_list,it))
	{
		int l=it->first;	if (-l>mxneg) mxneg=-l;
		int r=it->second;	if (-r>mxneg) mxneg=-r;
	}
	int bkup=cur_size;
	_clear_vec();
	_setup_size(bkup);
}

bool field_numbers::field_numbers_data::_setup_size(int s) // s should be .size(). return whether there's unsolved negative # 
{	if (mxneg>s) return false; // must be the first: Don't _clear_vec if it won't work.
	_clear_vec();
	int j=0;
	for (unsigned sec=0;sec<rng_sect.size();++sec)
	{
		for (;j<rng_sect[sec];++j)
		{
			int l=rng_list[j].first;	if (l<=0) { ++l; l+=s; } // l=++l+s;
			int r=rng_list[j].second;	if (r<=0) { ++r; r+=s; } // r=++r+s;
			if (l<=r)	{ for (int i=l;i<=r;++i) fn_vec.push_back(i-1); if (r>min_size) min_size=r; }
			else		{ for (int i=l;i>=r;--i) fn_vec.push_back(i-1); if (l>min_size) min_size=l; }
		}
		fn_sec.push_back(fn_vec.size());
	}
	cur_size=s;
	sizeNset=false;
	return true;
}

void field_numbers::field_numbers_data::_clear_vec() {
	min_size=0;
	cur_size=-1;
	sizeNset=true;
	fn_vec.clear();
	fn_sec.clear();
}

void field_numbers::clear()    { d->_clear_vec(); d->mxneg=-1; d->rng_list.clear(); d->rng_sect.clear(); }
field_numbers::~field_numbers(){ delete d; }
field_numbers::field_numbers(bool n,bool a):d(NULL){ d=new field_numbers_data; d->NewEachTime=n; d->AllRequired=a; }
field_numbers::field_numbers():d(NULL){ d=new field_numbers_data; }
field_numbers::field_numbers(const field_numbers& othr):d(NULL) { d=new field_numbers_data; *d=*(othr.d);}
field_numbers& field_numbers::operator=(const field_numbers& orig) {
	if(&orig != this)	*(this->d) = *(orig.d);
	return *this;
}

void field_numbers::set_NewEachTime(bool status) { d->NewEachTime=status; }
void field_numbers::set_AllRequired(bool status) { d->AllRequired=status; }

// These functions input data to rng_list
void field_numbers::parse_or_exit(std::istream& ss, bool one_input)
{
	size_t last=d->rng_list.size();
	d->_parse(ss);
	size_t now =d->rng_list.size();
	if (last!=now) // if (last==now) exit_error("No field numbers read."); // originally if (rng_list.empty()) but wrong
		if (one_input)
			if (now-last!=1 || d->rng_list.back().first!=d->rng_list.back().second)
				exit_error("Cannot input more than one field number");
}
string field_numbers::parse_or_exit(const string& str, bool one_input) { // return remaining string
	std::stringstream ss(str);
	parse_or_exit(ss,one_input);
	string s;
	getline(ss,s,(char)EOF);
	return s;
}
void field_numbers::push_back(int i) {
	d->rng_list.push_back(pair<int,int>(i,i));
	d->rng_sect.push_back(d->rng_list.size());
	if (-i>d->mxneg) d->mxneg=-i;
	int bkup=d->cur_size;
	d->_clear_vec();
	d->_setup_size(bkup);
}
void field_numbers::push_back(int i,int j) {
	d->rng_list.push_back(pair<int,int>(i,j));
	d->rng_sect.push_back(d->rng_list.size());
	if (-i>d->mxneg) d->mxneg=-i;
	if (-j>d->mxneg) d->mxneg=-j;
	int bkup=d->cur_size;
	d->_clear_vec();
	d->_setup_size(bkup);
}
string field_numbers::print() {
	string result;
	for (vector< pair<int,int> >::const_iterator i=d->rng_list.begin(); i!=d->rng_list.end(); ++i)
	{
		if (i!=d->rng_list.begin()) result.push_back(',');
		result += itos(i->first);
		if (i->first!=i->second) result.push_back('-'); else continue;
		if (i->second!=-1)  result += itos(i->second);
	}
	return result;
}

// Test whether size is good
bool field_numbers::size_solvable(int s)
{
	if (d->mxneg>=0 && (d->NewEachTime || d->sizeNset)) return d->_setup_size(s);
	return true; // not to return s>=mxneg, because at this line negative number has been solved, no need to solve again.
}

bool field_numbers::size_valid(int s) 
{
	if (size_solvable(s)==false) return false;
	if (d->AllRequired) return s>=d->min_size;
	return true;
}

int field_numbers::min_required() { // if mxneg>=0, it should be used after size_valid(), and NewEachTime should be false, only this way the min_size was calcualted.
	if (d->mxneg>=0 && d->NewEachTime) exit_error("Field_number is set to NewEachTime, and there're negative numbers; no way to known the minimum size requirement.");
	if (d->mxneg>=0 && d->sizeNset)   exit_error("Field_number has negative numbers and size is not solved yet; no way to know the minimum size requirement.");
	if (d->AllRequired)	return d->min_size;
	else				return 0; // not to return mxneg, because at this line negative number has been solved, no need to solve again.
}

template<typename Container>
void field_numbers::contents_to_doubles(const vector<string>& srce, bool clr, Container& dest, bool add_def, const double def_val) {
	if (clr) dest.clear();
	if (!size_solvable(srce.size())) throw bad_query_fieldnumbers_unsolved(); // previously return -1;
	int size=srce.size();
	for (each_element_const(d->fn_vec,it))
	{
		if (size>*it)
		{
			try { add_to_container(dest,boost::lexical_cast<double>(srce[*it])); }
			catch (boost::bad_lexical_cast &) { if (add_def) add_to_container(dest,def_val); }
		}
		else
		{
			if (d->AllRequired) throw bad_query_fieldnumbers_absent();
			if (add_def) add_to_container(dest,def_val);
		}
	}
}
template void field_numbers::contents_to_doubles(const vector<string>& srce, bool clr,   vector<double>& dest, bool add_def, const double def_val);
template void field_numbers::contents_to_doubles(const vector<string>& srce, bool clr,    deque<double>& dest, bool add_def, const double def_val);
template void field_numbers::contents_to_doubles(const vector<string>& srce, bool clr,      set<double>& dest, bool add_def, const double def_val);
template void field_numbers::contents_to_doubles(const vector<string>& srce, bool clr, multiset<double>& dest, bool add_def, const double def_val);
template void field_numbers::contents_to_doubles(const vector<string>& srce, bool clr,   Values<double>& dest, bool add_def, const double def_val);

template<typename Container>
void field_numbers::contents_to_strings(const vector<string>& srce, bool clr, Container& dest, bool add_def, const string& def_val) {
	if (clr) dest.clear();
	if (!size_solvable(srce.size())) throw bad_query_fieldnumbers_unsolved(); // previously return -1;
	int size=srce.size();
	for (each_element_const(d->fn_vec,it))
	{
		if (size>*it)
		{
			add_to_container(dest, srce[*it]);
		}
		else
		{
			if (d->AllRequired) throw bad_query_fieldnumbers_absent();
			if (add_def) add_to_container(dest,def_val);
		}
	}
}
template void field_numbers::contents_to_strings(const vector<string>& srce, bool clr,   vector<string>& dest, bool add_def, const string& def_val);
template void field_numbers::contents_to_strings(const vector<string>& srce, bool clr,    deque<string>& dest, bool add_def, const string& def_val);
template void field_numbers::contents_to_strings(const vector<string>& srce, bool clr,      set<string>& dest, bool add_def, const string& def_val);
template void field_numbers::contents_to_strings(const vector<string>& srce, bool clr, multiset<string>& dest, bool add_def, const string& def_val);
template void field_numbers::contents_to_strings(const vector<string>& srce, bool clr,   Values<string>& dest, bool add_def, const string& def_val);

template<typename T>
void field_numbers::contents_to_a_string(const vector<string>& srce, string& dest, const T& del, bool write_endl, bool quoted, const string& def_val) {
	dest.clear();
	if (!size_solvable(srce.size())) throw bad_query_fieldnumbers_unsolved(); // previously return -1;
	if (!d->fn_vec.empty())
	{
		int size=srce.size();
		vector<int>::const_iterator it (d->fn_vec.begin());
		if (quoted)
		{
			if (size>*it) {	dest+='"'; dest+=srce[*it]; dest+='"'; }
			else		  {	dest+='"'; dest+=def_val;   dest+='"'; if (d->AllRequired) throw bad_query_fieldnumbers_absent(); }
			for (++it; it!=d->fn_vec.end(); ++it) {
				dest += del;
				if (size>*it) {	dest+='"'; dest+=srce[*it]; dest+='"'; }
				else		  {	dest+='"'; dest+=def_val;   dest+='"'; if (d->AllRequired) throw bad_query_fieldnumbers_absent(); }
			}
		}
		else
		{
			if (size>*it) {	dest+=srce[*it]; }
			else		  {	dest+=def_val;   if (d->AllRequired) throw bad_query_fieldnumbers_absent(); }
			for (++it; it!=d->fn_vec.end(); ++it) {
				dest += del;
				if (size>*it) {	dest+=srce[*it]; }
				else		  {	dest+=def_val;   if (d->AllRequired) throw bad_query_fieldnumbers_absent(); }
			}
		}
	}
	if (write_endl) dest.push_back('\n');
}
template void field_numbers::contents_to_a_string(const vector<string>& srce, string& dest, const string& del, bool write_endl, bool quoted, const string& def_val);
template void field_numbers::contents_to_a_string(const vector<string>& srce, string& dest, const char&   del, bool write_endl, bool quoted, const string& def_val);

template<typename T>
void field_numbers::contents_to_ostream(const vector<string>& srce, ostream& dest, const T& del, bool write_endl, bool quoted, const string& def_val) {
	if (!size_solvable(srce.size())) throw bad_query_fieldnumbers_unsolved(); // previously return -1;
	if (!d->fn_vec.empty())
	{
		int size=srce.size();
		vector<int>::const_iterator it (d->fn_vec.begin());
		if (quoted)
		{
			if (size>*it) {	dest << '"' << srce[*it] << '"'; }
			else		  {	dest << '"' << def_val   << '"'; if (d->AllRequired) throw bad_query_fieldnumbers_absent(); }
			for (++it; it!=d->fn_vec.end(); ++it) {
				dest << del;
				if (size>*it) {	dest << '"' << srce[*it] << '"'; }
				else		  {	dest << '"' << def_val   << '"'; if (d->AllRequired) throw bad_query_fieldnumbers_absent(); }
			}
		}
		else
		{
			if (size>*it) {	dest << srce[*it]; }
			else		  {	dest << def_val;   if (d->AllRequired) throw bad_query_fieldnumbers_absent(); }
			for (++it; it!=d->fn_vec.end(); ++it) {
				dest << del;
				if (size>*it) {	dest << srce[*it]; }
				else		  {	dest << def_val;   if (d->AllRequired) throw bad_query_fieldnumbers_absent(); }
			}
		}
	}
	if (write_endl) dest<<endl; // previously '\n' but it doesn't flush and severely delays the display to stdout
}
template void field_numbers::contents_to_ostream(const vector<string>& srce, ostream& dest, const string& del, bool write_endl, bool quoted, const string& def_val);
template void field_numbers::contents_to_ostream(const vector<string>& srce, ostream& dest, const char&   del, bool write_endl, bool quoted, const string& def_val);

////////////////////////
// rng related functions
bool field_numbers::no_input() {
	return d->rng_list.empty();
}
int field_numbers::num_inputs() { 
	return d->rng_sect.size(); // previously rng_list, changed on 05-06-2011
}
bool field_numbers::is_double() {
	if (d->rng_list.size()==2)
	{
		if (d->rng_list[0].first != d->rng_list[0].second) return false;
		if (d->rng_list[1].first != d->rng_list[1].second) return false;
		return true;
	}
	else return false;
}
bool field_numbers::is_multi() {
	if (d->rng_list.empty()) return false;
	if (d->rng_list.size()>1) return true;
	return d->rng_list.front().first != d->rng_list.front().second;
}
bool field_numbers::is_multi(int sec) {
	if (d->rng_list.empty()) return false;
	int prv = sec ? d->rng_sect[sec-1] : 0;
	int ths = d->rng_sect[sec];
	if (ths-prv>1) return true;
	return d->rng_list[prv].first != d->rng_list[prv].second;
}

////////////////////////
// vec related functions
int field_numbers::max_field_num() {
	if (d->mxneg>=0 && d->sizeNset) exit_error("Encountered an unsolved negative field number. This could be due to\n 1) insufficient number of fields to solve negative field numbers;\n 2) bug in program ( failed to call size_valid before using []).");
	int maxnum=-INT_MAX;
	for (each_element_const(d->fn_vec,it)) if (*it>maxnum) maxnum=*it;
	return maxnum;	}

bool field_numbers::contain(const int& n) const {
	if (d->mxneg>=0 && d->sizeNset) exit_error("Encountered an unsolved negative field number. This could be due to\n 1) insufficient number of fields to solve negative field numbers;\n 2) bug in program ( failed to call size_valid before using []).");
	for (each_element_const(d->fn_vec,it)) if (*it==n) return true;
	return false; }

vector<int>::const_iterator field_numbers::find(const int& n) const {
	if (d->mxneg>=0 && d->sizeNset) exit_error("Encountered an unsolved negative field number. This could be due to\n 1) insufficient number of fields to solve negative field numbers;\n 2) bug in program ( failed to call size_valid before using []).");
	for (each_element_const(d->fn_vec,it)) if (*it==n) return it;
	return end(); }

int field_numbers::operator[]( int pos ) { // return by value rather than by reference, so that it can't be changed.
	if (d->mxneg>=0 && d->sizeNset) exit_error("Encountered an unsolved negative field number. This could be due to\n 1) insufficient number of fields to solve negative field numbers;\n 2) bug in program ( failed to call size_valid before using []).");
	if (pos>=(int)d->fn_vec.size()) exit_error("Query number exceed input length of field_number.");
	if (pos<0)				exit_error("field_number_vec query number cannot be negative.");
	return d->fn_vec[pos];
}

vector<int>::const_iterator field_numbers::begin() const {
	if (d->mxneg>=0 && d->sizeNset) exit_error("Encountered an unsolved negative field number. This could be due to\n 1) insufficient number of fields to solve negative field numbers;\n 2) bug in program ( failed to call size_valid before using [])."); 
	return d->fn_vec.begin(); }

vector<int>::const_iterator field_numbers::begin(int sec) const {
	if (d->mxneg>=0 && d->sizeNset) exit_error("Encountered an unsolved negative field number. This could be due to\n 1) insufficient number of fields to solve negative field numbers;\n 2) bug in program ( failed to call size_valid before using [])."); 
	if (sec)	return d->fn_vec.begin()+d->fn_sec[sec-1];
	else		return d->fn_vec.begin(); }

vector<int>::const_iterator field_numbers::end(int sec) const {
	if (d->mxneg>=0 && d->sizeNset) exit_error("Encountered an unsolved negative field number. This could be due to\n 1) insufficient number of fields to solve negative field numbers;\n 2) bug in program ( failed to call size_valid before using [])."); 
	return d->fn_vec.begin()+d->fn_sec[sec]; }

vector<int>::const_reverse_iterator field_numbers::rbegin() const {
	if (d->mxneg>=0 && d->sizeNset) exit_error("Encountered an unsolved negative field number. This could be due to\n 1) insufficient number of fields to solve negative field numbers;\n 2) bug in program ( failed to call size_valid before using [])."); 
	return d->fn_vec.rbegin(); }

int& field_numbers::front() {
	if (d->mxneg>=0 && d->sizeNset) exit_error("Encountered an unsolved negative field number. This could be due to\n 1) insufficient number of fields to solve negative field numbers;\n 2) bug in program ( failed to call size_valid before using [])."); 
	return *(d->fn_vec.begin()); }

int& field_numbers::front(int sec) {
	if (d->mxneg>=0 && d->sizeNset) exit_error("Encountered an unsolved negative field number. This could be due to\n 1) insufficient number of fields to solve negative field numbers;\n 2) bug in program ( failed to call size_valid before using [])."); 
	if (sec)	return *(d->fn_vec.begin()+d->fn_sec[sec-1]);
	else		return *(d->fn_vec.begin()); }

vector<int>::size_type field_numbers::size() const {
	if (d->mxneg>=0 && d->sizeNset) exit_error("Encountered an unsolved negative field number. This could be due to\n 1) insufficient number of fields to solve negative field numbers;\n 2) bug in program ( failed to call size_valid before using [])."); 
	return d->fn_vec.size(); }

bool field_numbers::empty() const {
	if (d->mxneg>=0 && d->sizeNset) exit_error("Encountered an unsolved negative field number. This could be due to\n 1) insufficient number of fields to solve negative field numbers;\n 2) bug in program ( failed to call size_valid before using [])."); 
	return d->fn_vec.empty(); }

// still vec related functions, but no validity checking only to increase speed
vector<int>::const_iterator field_numbers::end() const { return d->fn_vec.end(); }
vector<int>::const_reverse_iterator field_numbers::rend() const { return d->fn_vec.rend(); }

// vec related operator, not used as of 2010 Dec 14
field_numbers operator&(const field_numbers& left, const field_numbers& right)  // previously *
{
	field_numbers result(left.d->NewEachTime|right.d->NewEachTime, left.d->AllRequired|right.d->AllRequired);
	if ( left.d->mxneg>=0 &&  left.d->sizeNset) exit_error("Encountered an unsolved negative field number. This could be due to\n 1) insufficient number of fields to solve negative field numbers;\n 2) bug in program ( failed to call size_valid before using []).");
	if (right.d->mxneg>=0 && right.d->sizeNset) exit_error("Encountered an unsolved negative field number. This could be due to\n 1) insufficient number of fields to solve negative field numbers;\n 2) bug in program ( failed to call size_valid before using []).");
	for (each_element_const(left.d->fn_vec,itl))
		for (each_element_const(right.d->fn_vec,itr))
			if (*itl==*itr) { result.push_back(*itl+1); break; }
	return result;
}

field_numbers operator|(const field_numbers& left, const field_numbers& right)
{
	field_numbers result(left.d->NewEachTime|right.d->NewEachTime, left.d->AllRequired|right.d->AllRequired);
	if ( left.d->mxneg>=0 &&  left.d->sizeNset) exit_error("Encountered an unsolved negative field number. This could be due to\n 1) insufficient number of fields to solve negative field numbers;\n 2) bug in program ( failed to call size_valid before using []).");
	if (right.d->mxneg>=0 && right.d->sizeNset) exit_error("Encountered an unsolved negative field number. This could be due to\n 1) insufficient number of fields to solve negative field numbers;\n 2) bug in program ( failed to call size_valid before using []).");
	set<int> allfld;
	for (each_element_const(left.d->fn_vec,itl))
		if (!exist_element(allfld,*itl)) { allfld.insert(*itl); result.push_back(*itl+1); }
	for (each_element_const(right.d->fn_vec,itr))
		if (!exist_element(allfld,*itr)) { allfld.insert(*itr); result.push_back(*itr+1); }
	return result;
}	

field_numbers operator+(const field_numbers& left, const field_numbers& right)
{
	field_numbers result(left.d->NewEachTime|right.d->NewEachTime, left.d->AllRequired|right.d->AllRequired);
	if ( left.d->mxneg>=0 &&  left.d->sizeNset) exit_error("Encountered an unsolved negative field number. This could be due to\n 1) insufficient number of fields to solve negative field numbers;\n 2) bug in program ( failed to call size_valid before using []).");
	if (right.d->mxneg>=0 && right.d->sizeNset) exit_error("Encountered an unsolved negative field number. This could be due to\n 1) insufficient number of fields to solve negative field numbers;\n 2) bug in program ( failed to call size_valid before using []).");
	for (each_element_const(left.d->fn_vec,itl))	result.push_back(*itl+1);
	for (each_element_const(right.d->fn_vec,itr))	result.push_back(*itr+1);
	return result;
}

field_numbers operator-(const field_numbers& left, const field_numbers& right)
{
	field_numbers result(left.d->NewEachTime|right.d->NewEachTime, left.d->AllRequired|right.d->AllRequired);
	if ( left.d->mxneg>=0 &&  left.d->sizeNset) exit_error("Encountered an unsolved negative field number. This could be due to\n 1) insufficient number of fields to solve negative field numbers;\n 2) bug in program ( failed to call size_valid before using []).");
	if (right.d->mxneg>=0 && right.d->sizeNset) exit_error("Encountered an unsolved negative field number. This could be due to\n 1) insufficient number of fields to solve negative field numbers;\n 2) bug in program ( failed to call size_valid before using []).");
	set<int> allfld;
	for (each_element_const(right.d->fn_vec,itr))	allfld.insert(*itr);
	for (each_element_const(left.d->fn_vec,itl))	if (!exist_element(allfld,*itl)) result.push_back(*itl+1);
	return result;
}

int ReadFld(const std::string& input, field_numbers& f, bool one_input, int ErrCode)
{
	string remains = f.parse_or_exit(input,one_input);
	if (remains.empty()){ return 1; }
	else if	(ErrCode>0)	{ elog.add(ErrCode); return 0; }
	else if (ErrCode<0)	{ exit_error("Fail to parse field_number, remaining string = "+remains); }
	else				{ return 0; }
	return -1; // never happen
}

void ReadArg(const std::vector<std::string>& args, size_t& argi, field_numbers& values, bool one_input, int ErrCode)
{
	std::string input;
	std::size_t found = args[argi].find('=');
	if (found!=std::string::npos)	{ NeedArg(args,0,argi); input = args[argi].substr(found+1); }
	else							{ NeedArg(args,1,argi); input = args[++argi]; }
	ReadFld(input,values,one_input,ErrCode);
}
