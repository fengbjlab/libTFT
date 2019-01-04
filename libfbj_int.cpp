#include <muParser.h>
#include <vector>
#include <deque>
#include "libfbj_base.hpp"
#include "libfbj_file.hpp"
#include "libfbj_int.hpp"
#include "libfbj_program.hpp" // for license checking

using namespace std;

struct IntRanges::IntRangesData {
	static int	license_chk;// licensing of the class
	
	// mode 1 -- integer range vector
	vector< pair<int,int> > alt_vec_;
	bool no_inv_;	// ignore inverse range, ie, n1-n2 where n2<n1
	
	// mode 2 -- math expression with/without iterators
	vector<int> max_val_,min_val_,cur_val_,incr_;
	vector<double> use_val_;
	mu::Parser p;
	string exp_s;	// expression string
	bool valid;		// curr_num_ is valid; for math expr, valid=false if start-of-a-new-round-except-the-1st-round, even no iterator exists
	char bgn_m;		// begin character for math expression
	char end_m;		// end character for math expression
	int  curr_num_;
	
	// mode 3 -- a file, content is lines of integer range vectors
	char bgn_f;
	char end_f;
	
	// flags or internal data
	bool use_alt_;	// flag for mode, use alt_vec_
	bool use_pre_;	// flag for bool, use pre_vec_
	deque<bool> pre_vec_;
	// pre_vec_ stores the pre-calculated flags for each byte, which increases the speed A LOT by
	// avoiding checking alt_vec_ for each bytes. It will be turned on if all the following is true :
	// 1) has_inf_==false
	// 2) use_alt_==true  (redundant, because if has_inf_==false math.expr will be converted to _alt_vec)
	// 3) idta->alt_vec_.size()>10
	
	// result
	bool has_inf_;	// has +infinity
	int max_,min_;	// max / min of all numbers; =INT_MAX / INT_MIN if has_inf_ and is math expr
	
	// functions
	IntRangesData():no_inv_(false),valid(false),bgn_m('{'),end_m('}'),curr_num_(INT_MIN),bgn_f(':'),end_f(':'),use_alt_(true),use_pre_(false),has_inf_(false),max_(INT_MIN),min_(INT_MAX){}
	IntRangesData(bool nr):no_inv_(nr),valid(false),bgn_m('{'),end_m('}'),curr_num_(INT_MIN),bgn_f(':'),end_f(':'),use_alt_(true),use_pre_(false),has_inf_(false),max_(INT_MIN),min_(INT_MAX){}
	void next() {
		int _old_num_ = curr_num_;
		bool next_is_valid=false;	// next try is valid
		for (int i=cur_val_.size()-1;i>=0;i--) {
			if (cur_val_[i]==max_val_[i]) { cur_val_[i]=min_val_[i]; use_val_[i]=min_val_[i]; continue; }
			cur_val_[i]+=incr_[i]; use_val_[i]=cur_val_[i]; next_is_valid=true; break;
		}
		valid=next_is_valid;
		curr_num_=p.Eval();
		if (curr_num_ < _old_num_ && has_inf_) exit_error("Int-Range requires that math_expr w/ infinities must be monotonously increasing.");
	}
	void rewind() {
		for (int i=cur_val_.size()-1;i>=0;i--) {
			cur_val_[i]=min_val_[i];
			use_val_[i]=min_val_[i];
		}
		valid=true;
		curr_num_=p.Eval();
	}
	void math_expr_to_alt_vec_() {
		alt_vec_.clear();
		for (rewind(); valid; next()) 
			alt_vec_.push_back(pair<int,int>(curr_num_,curr_num_));
		curr_num_=alt_vec_.front().first;
	}
	void clear() {
		alt_vec_.clear();
		max_val_.clear();
		min_val_.clear();
		cur_val_.clear();
		incr_.clear();
		use_val_.clear();
		exp_s.clear();
		valid=false;
		use_alt_=true;
		curr_num_=INT_MIN;
		has_inf_=false;
		min_ = INT_MAX;
		max_ = INT_MIN;
	}
};

// initialize
int IntRanges::IntRangesData::license_chk = program.nudge;
IntRanges::~IntRanges(){ delete idta; }
IntRanges::IntRanges():idta(NULL){ idta=new IntRangesData; }
IntRanges::IntRanges(bool nr):idta(NULL){ idta=new IntRangesData(nr); }
IntRanges::IntRanges(const IntRanges& othr):idta(NULL) { idta=new IntRangesData; *idta = *(othr.idta); }
IntRanges& IntRanges::operator=(const IntRanges& orig)
{
	if (&orig!=this)
	{
		*(this->idta)=*(orig.idta);
		if (!this->idta->exp_s.empty())
		{
			this->idta->p.ClearVar();
			for (unsigned i=0;i<this->idta->use_val_.size();++i) this->idta->p.DefineVar(s(char(108+i)),&(this->idta->use_val_[i]));
			try { this->idta->p.SetExpr(this->idta->exp_s); }
			catch (mu::Parser::exception_type &e) { exit_error("Error parsing math expression "+idta->exp_s+" : "+e.GetMsg()); }
		}
	}
	return *this;
}

void IntRanges::operator=(const string& in)
{
	string remains = in;
	int something_read = parse(remains);
	if (!something_read || !remains.empty()) exit_error("Error parsing "+in+" as an IntRanges object.");
}

void IntRanges::set_math_signal(char c) {
	switch (c) {
		case '{': idta->bgn_m=c; idta->end_m='}'; break;
		case '(': idta->bgn_m=c; idta->end_m=')'; break;
		case '[': idta->bgn_m=c; idta->end_m=']'; break;
		case '<': idta->bgn_m=c; idta->end_m='>'; break;
		default:  idta->bgn_m=   idta->end_m=c;   break;
	}
}

void IntRanges::set_file_signal(char c) {
	switch (c) {
		case '{': idta->bgn_f=c; idta->end_f='}'; break;
		case '(': idta->bgn_f=c; idta->end_f=')'; break;
		case '[': idta->bgn_f=c; idta->end_f=']'; break;
		case '<': idta->bgn_f=c; idta->end_f='>'; break;
		default:  idta->bgn_f=   idta->end_f=c;   break;
	}
}

void IntRanges::clear() { idta->clear(); }

// setup data
void IntRanges::push_back(int n) {
	if (n==INT_MAX) exit_error("IntRanges::push_back(a) does not allow INT_MAX.");
	idta->alt_vec_.push_back(pair<int,int>(n,n));
	if (n > idta->max_) idta->max_=n;
	if (n < idta->min_) idta->min_=n;
	idta->use_alt_=true;
	idta->curr_num_=idta->alt_vec_.front().first;
	if (idta->use_pre_) vec_deq_set(idta->pre_vec_,n,true);
}
void IntRanges::push_back_n_to_max(int n) {
	if (n==INT_MAX) exit_error("IntRanges::push_back_n_to_max(a) does not allow INT_MAX.");
	idta->alt_vec_.push_back(pair<int,int>(n,INT_MAX));
	idta->has_inf_=true;
	if (n < idta->min_) idta->min_=n;
	idta->max_=INT_MAX;
	idta->use_alt_=true;
	idta->curr_num_=idta->alt_vec_.front().first;
}
int IntRanges::parse(std::string& str) { // return # read, es become the rest
	std::istringstream ss(str);
	int to_return=parse(ss);
	str.clear();
	getline(ss,str,(char)EOF); // while (!ss.eof()) str.push_back(ss.get()); is wrong
	return to_return;
}
int IntRanges::parse(std::istream& es) // return # read
{
	clear();
	skip_whitespaces(es);
	if (es.peek()!=idta->bgn_m && es.peek()!=idta->bgn_f) {
		int n = parse_num_pairs(es,idta->alt_vec_,idta->no_inv_);
		if (n) {
			idta->has_inf_	= num_pairs_includes_inf(idta->alt_vec_);
			idta->max_		= num_pairs_max(idta->alt_vec_);
			idta->min_		= num_pairs_min(idta->alt_vec_);
			if (idta->alt_vec_.size()>10 && !idta->has_inf_) idta->use_pre_=true;
			if (idta->use_pre_) num_pairs_to_vector(idta->alt_vec_,idta->pre_vec_);
		}
		return n;
	}
	else if (es.peek()==idta->bgn_f) {
		es.get(); // get idta->bgn_f
		int n=0;
		string filename;
		osi_getline_allowing_brackets_or_quotes(es,filename,idta->end_f);
		for (Rows_in_File(inf,filename,1)) {
			if (inf.contents()[0].empty()) exit_error("read integer range from file "+filename+" error: empty string.");
			n += parse_num_pairs(inf.contents()[0],idta->alt_vec_,idta->no_inv_);
			if (!inf.contents()[0].empty()) exit_error("read integer range from file "+filename+" error: remaining string is "+inf.contents()[0]);
		}
		if (n) {
			idta->has_inf_	= num_pairs_includes_inf(idta->alt_vec_);
			idta->max_		= num_pairs_max(idta->alt_vec_);
			idta->min_		= num_pairs_min(idta->alt_vec_);
			if (idta->alt_vec_.size()>10 && !idta->has_inf_) idta->use_pre_=true;
			if (idta->use_pre_) num_pairs_to_vector(idta->alt_vec_,idta->pre_vec_);
		}
		return n;
	}
	char c;
	es>>c;
	skip_whitespaces(es);
	for (int nv=0;es.peek()=='[';++nv)
	{
		int fr,to;
		es>>c>>fr>>c;	// [ fr -
		if (es.peek()==']') to=INT_MAX;
		else es>>to;	// to
		es>>c;			// ]
		skip_whitespaces(es);
		idta->min_val_.push_back(fr);
		idta->max_val_.push_back(to);
		if (to>=fr)	idta->incr_.push_back(1);
		else		idta->incr_.push_back(-1);
		idta->cur_val_.push_back(fr);
		idta->use_val_.push_back(double(fr));
		if (fr==INT_MAX) idta->has_inf_=true;
		if (to==INT_MAX) idta->has_inf_=true;
		// p.DefineVar(s(char(97+nv)),&use_val_[nv]);
		// The above code is wrong! because address of all use_val_[i] changes after another push_back() operation
		// So only the last address will be correct. To DefineVar, have to define them all after use_val_ is complete.
	}
	for (unsigned i=0;i<idta->use_val_.size();++i) idta->p.DefineVar(s(char(108+i)),&(idta->use_val_[i]));
	osi_getline_allowing_brackets_or_quotes(es,idta->exp_s,idta->end_m); // previously while ((c=es.get())!=idta->end_m) ex.push_back(c);
	try {	idta->p.SetExpr(idta->exp_s); }
	catch (mu::Parser::exception_type &e)	{ exit_error("Error parsing math expression "+idta->exp_s+" : "+e.GetMsg()); }
	idta->rewind();
	if (idta->has_inf_) {
		for (int i=0;i<100;++i) idta->next();
		idta->rewind();
		idta->min_ = idta->curr_num_; // prv: INT_MIN;
		idta->max_ = INT_MAX;
		idta->use_alt_=false;
	}
	else {
		idta->math_expr_to_alt_vec_();
		idta->max_= num_pairs_max(idta->alt_vec_);
		idta->min_= num_pairs_min(idta->alt_vec_);
		idta->use_alt_=true;
		if (idta->alt_vec_.size()>10) idta->use_pre_=true; // idta->has_inf_ is already false
		if (idta->use_pre_) num_pairs_to_vector(idta->alt_vec_,idta->pre_vec_);
	}
	return 1;
}
void IntRanges::print(std::ostream& os) const
{
	if (idta->use_alt_) print_num_pairs(idta->alt_vec_,os);
	else os<<"{"<<idta->exp_s<<"}";
}

// access data
bool IntRanges::empty()			{ return idta->use_alt_ && idta->alt_vec_.empty(); }
bool IntRanges::has_infinity()	{ return idta->has_inf_; }
int  IntRanges::max()			{ return idta->max_; }
int  IntRanges::min()			{ return idta->min_; }
bool IntRanges::include(int n) {
	if (idta->use_alt_)
	{
		if (idta->use_pre_)	return vec_deq_get(idta->pre_vec_,n);
		else				return num_pairs_includes(idta->alt_vec_,n);
	}
	if (idta->curr_num_>n) idta->rewind();
	while (n>idta->curr_num_) { idta->next(); if (!idta->valid) break;}
	return idta->curr_num_==n;
}
bool IntRanges::is_n_to_max(int n) {
	if (!idta->use_alt_) return false;
	if (idta->alt_vec_.size()!=1) return false;
	if (idta->alt_vec_[0]!=pair<int,int>(n,INT_MAX)) return false;
	return true;
}

