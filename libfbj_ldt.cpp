#include <STLplus/containers/ntree.hpp>
#include <sstream>
#include <set>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include "libfbj_base.hpp"
#include "libfbj_file.hpp"
#include "libfbj_program.hpp"
#include "libfbj_fn.hpp"
#include "libfbj_int.hpp"
#include "libfbj_expr.hpp"
#include "libfbj_ldt.hpp"

using std::vector;
using std::string;
using std::istream;
using std::ostream;
using std::pair;
using std::set;
using stlplus::ntree;

string cond_row_selection_help_text=
"ProgramHandle, please enable regex."
"FOR ROW SELECTION\n"
" --if MC           select lines by multiple conditions MC\n"
" --renew MC        if satisfied, zero out history of --if and group number +1\n"
" --until MC        if satisfied, skip the line and all subsequent lines\n"
;

namespace {
	class row_index {
	public:
		bool			required;	// this srowbi is required (has data input)
		field_numbers	index1;
		field_numbers	index2;
		string			idxf_name;
		set<string>		indices;
		row_index():required(false),index1(false,true),index2(false,true){}
		void input(field_numbers& i1, string& s, bool cis)
		{
			index1=i1;
			index1.set_NewEachTime(false);
			index1.set_AllRequired(true);
			idxf_name=index2.parse_or_exit(s).substr(1);
			tfile_format fmt(default_tfile_format);
			fmt.set_field_nums(index2, "line(s) in "+idxf_name+" lack some index fields.", tfile_format::Continue);
			for (Rows_in_File(inf,idxf_name,&fmt))
			{
				string idx;
				index2.contents_to_a_string(inf.contents(),idx,DLMTR);
				if (cis) boost::to_lower(idx);
				indices.insert(idx);
			}
			required=true;
		}
		pair<bool,bool> match(const vector<string>& v)
		{
			if (required)
			{
				if (!index1.size_valid(v.size())) return pair<bool,bool>(false,false);
				string idx;
				index1.contents_to_a_string(v,idx,DLMTR);
				return pair<bool,bool>(true,exist_element(indices,idx));
			}
			return pair<bool,bool>(true,true);
		}
	};
	
	enum row_selection_method { on_ea=0, af_1c, oa_1c, bf_1c, ob_1c, af_ea, bf_ea };
	
	class row_selection {
	public:
		row_selection():past_sat(false),rn(-1),fr(0),md(0) { satrow_rn.push_back_n_to_max(1); satrow_fr.push_back_n_to_max(1); }
		void renew() {	past_sat=false; rn=-1; }
		void analyze(bool& expr_sat)
		{
			bool rstr_sat = expr_sat? satrow_fr.include(++fr) : false; // satisfaction restricted by occurence numbers
			if  (rstr_sat) { past_sat=true; rn=0; }	else rn += rn!=-1; // passed any restricted_satisfaction
			switch (md)
			{
				case on_ea:	expr_sat =  rstr_sat;	break;
				case af_1c:	expr_sat =  past_sat & satrow_rn.include(rn);	break;
				case oa_1c: expr_sat =  past_sat;	break;
				case bf_1c:	expr_sat = !past_sat;	break;
				case ob_1c: expr_sat = !past_sat || rstr_sat; break;
				case af_ea:	expr_sat =  past_sat & satrow_fr.include(fr) & satrow_rn.include(rn); break;
				case bf_ea:	expr_sat = !expr_sat & satrow_fr.include(fr+1); break;
				default   : exit_error("single/multiple condition: wrong P (position)");	break;
			}
		}
		void read(std::istream& ss)
		{
			// read L (line numbers)
			skip_whitespaces(ss);
			bool satrow_rn_read=false;
			satrow_rn.clear();
			if (satrow_rn.parse(ss)) satrow_rn_read=true;
			else satrow_rn.push_back_n_to_max(1);
			if (satrow_rn.min()<0) exit_error("single/multiple condition: L (line numbers) cannot be <0.");
			
			// read P (position)
			skip_whitespaces(ss);
			char c1,c2;
			switch (ss.peek())
			{
				case 'o': ss>>c1>>c2; if (c2=='a') md = oa_1c; else if (c2=='b') md = ob_1c; else exit_error("wrong P (position)"); break;
				case 'a': ss>>c1>>c2; if (c2=='f') md = af_1c; else if (c2=='e') md = af_ea; else exit_error("wrong P (position)"); break;
				case 'b': ss>>c1>>c2; if (c2=='f') md = bf_1c; else if (c2=='e') md = bf_ea; else exit_error("wrong P (position)"); break;
				default:
					if (satrow_rn_read)
					{
						satrow_fr=satrow_rn;
						satrow_rn.clear();
						satrow_rn.push_back(0);
						return;
					}
					break;
			}
			if (md==oa_1c || md==ob_1c || md==bf_1c || md==bf_ea)
				if (satrow_rn_read) exit_error("single/multiple condition: L (line numbers) is not allowed for oa/ob/bf/be.");
			
			// read N (Nth satisfaction, satrow_fr)
			skip_whitespaces(ss);
			satrow_fr.clear();
			if (satrow_fr.parse(ss)==0)
			{
				if (md==oa_1c || md==ob_1c || md==af_1c || md==bf_1c)
					satrow_fr.push_back(1);
				else
					satrow_fr.push_back_n_to_max(1);
			}
			else
			{
				if (satrow_fr.min()<1) exit_error("single/multiple condition: N (occurrence numbers) cannot be <1.");
				if (md==oa_1c || md==ob_1c || md==af_1c || md==bf_1c)
					if (satrow_fr.min()!=satrow_fr.max()) exit_error("single/multiple condition: when P=oa/ob/af/bf, N should be singular.");
			}
		}
		void pr(IntRanges& ir, ostream& os) // print to os if ir is not default [1,max]
		{
			if (!ir.is_n_to_max(1)) ir.print(os);
		}
		void write(std::ostream& os)
		{
			switch (md)
			{
				case on_ea:	                            pr(satrow_fr,os); break;
				case af_1c: pr(satrow_rn,os); os<<"af"; pr(satrow_fr,os); break;
				case oa_1c:                   os<<"oa"; pr(satrow_fr,os); break;
				case bf_1c:                   os<<"bf"; pr(satrow_fr,os); break;
				case ob_1c:                   os<<"ob"; pr(satrow_fr,os); break;
				case af_ea: pr(satrow_rn,os); os<<"ae"; pr(satrow_fr,os); break;
				case bf_ea:                   os<<"be"; pr(satrow_fr,os); break;
				default:exit_error("single/multiple condition: wrong P (position).");
			}
		}
	private:
		IntRanges satrow_rn;
		IntRanges satrow_fr;
		bool past_sat;
		int rn;
		int fr;
		int md;
	};
	
}

// ------------- single_condition -------------

struct single_condition::SingleConditionData {
	static int  license_chk;// license checking
	char		comparison;
	bool		tan;		// test as numeric fields
	bool		lao;		// logical AND operation to combine results from individual fields
	bool		ng1;		// logical negation 1
	bool		ng2;		// logical negation 2
	bool		cis;		// case insensitive
	bool		tnf;		// test number of fields
	bool		tfl;		// test field length
	int			error1,error2,error3,error4;
	expressions::math_expression	rMath;	// reference math expression
	expressions::strg_expression	rStrg;	// reference strg expression
	field_numbers flds;	// original inputed field ranges
	string		input_ref;
	set<string> refdb;
	row_index	srowbi;
	boost::regex filter;
	boost::smatch what;	// no use, but needed
	row_selection rws;
	SingleConditionData():comparison('n'),tan(false),lao(true),ng1(false),ng2(false),cis(false),tnf(false),tfl(false),flds(false,true){
		error1=elog.get_token("single conditional statement(s) cannot be tested due to non-numeric fields:");
		error2=elog.get_token("single conditional statement(s) cannot be tested due to the lack of fields.");
		error3=elog.get_token("single conditional statement(s) cannot be tested due to math expression evaluation failures.");
		error4=elog.get_token("single conditional statement(s) cannot be tested due to string expression evaluation failures.");
	}
};

int single_condition::SingleConditionData::license_chk = program.nudge;

char single_condition::type_of_test() {
	return d->comparison;
}

string single_condition::ref_string() {
	return d->input_ref;
}

int single_condition::first_field() {
	return d->flds.front();
}

single_condition::~single_condition() {
	delete d;
}

single_condition::single_condition():d(NULL) {
	d=new SingleConditionData;
}

single_condition::single_condition(const single_condition& othr):d(NULL) {
	d=new SingleConditionData;
	*d = *(othr.d);
}

single_condition& single_condition::operator=(const single_condition& orig) {
	if(&orig != this)	*(this->d) = *(orig.d);
	return *this;
}

void single_condition::initialize() {
	d->comparison='n';
	d->tan=false;
	d->lao=true;
	d->ng1=false;
	d->ng2=false;
	d->cis=false;
	d->tnf=false;
	d->tfl=false;
	d->input_ref.clear();
	d->refdb.clear();
	d->flds.clear();
}

int  single_condition::max_field_num() {
	return d->flds.max_field_num();
}

bool single_condition::size_valid(int size)	{
	return d->flds.size_valid(size);
}

void single_condition::suppress_elog() {
	d->error1 = d->error2 = d->error3 = d->error4 = -1;
}

pair<bool,bool> single_condition::evaluate(const vector<string>& sv) {
	pair<bool,bool> r=eval_without_inverse(sv);
	if (r.first && d->ng2) r.second=!r.second;
	d->rws.analyze(r.second);
	return r;
}

pair<bool,bool> single_condition::eval_without_inverse(const vector<string>& sv) {
	if (d->tnf)
	{
		if (!d->rMath.solve()) { elog.add(d->error3); return pair<bool,bool>(false,false); }
		double refnum=d->rMath.val();
		bool nf_result;
		switch (d->comparison)
		{
			case 'e': nf_result = sv.size() == refnum; break;
			case 'g': nf_result = sv.size() >= refnum; break;
			case 'G': nf_result = sv.size() >  refnum; break;
			case 'l': nf_result = sv.size() <= refnum; break;
			case 'L': nf_result = sv.size() <  refnum; break;
			case 'n': nf_result = sv.size() != refnum; break;
			default: exit_error("single condition: wrong comparison for NF (number of fields).");
		}
		if (d->ng1) nf_result = !nf_result;
		return pair<bool,bool>(true,nf_result);
	}
	if (! d->flds.size_valid(sv.size()))
	{	elog.add(d->error2);
		return pair<bool,bool>(false,false);
	}
	// Now this program does not allow mising even one field, although logically it's still possible to evaluate the condition,
	// say [|s1-eq] and there's 1 field empty but some others missing, the result could be true, but here is <false,false>.
	// This is so because field_numbers::[] do not allow it. In the future, if it really requires that some fields are allowed
	// to be missing, use flds(true,true).
	bool r=d->lao;
	if (d->tan)
	{
		if (d->comparison=='S') // is missing, no refnum
		{
			for (each_element(d->flds,it))
			{
				bool ers;
				try { boost::lexical_cast<double>(sv[*it]); ers=false; }
				catch (boost::bad_lexical_cast &)		{	ers=true; }
				if (d->ng1) ers=!ers;
				if (d->lao) r&=ers; else r|=ers;
			}
		}
		else // other tests, with refnum
		{
			if (!d->rMath.solve()) { elog.add(d->error3); return pair<bool,bool>(false,false); }
			double refnum=d->rMath.val();
			for (each_element(d->flds,it))
			{
				double n;
				if (d->tfl)
				{
					n=sv[*it].length();
				}
				else
				{
					try { n=boost::lexical_cast<double>(sv[*it]); }
					catch (boost::bad_lexical_cast &) {	elog.add(d->error1,sv[*it]); return pair<bool,bool>(false,false); }
				}
				bool ers;
				switch (d->comparison)
				{
					case 'e': ers = (n==refnum); break;
					case 'g': ers = (n>=refnum); break;
					case 'G': ers = (n> refnum); break;
					case 'l': ers = (n<=refnum); break;
					case 'L': ers = (n< refnum); break;
					case 'n': ers = (n!=refnum); break;
					default: exit_error("single condition: wrong comparison for a numeric field.");
				}
				if (d->ng1) ers=!ers;
				if (d->lao) r&=ers; else r|=ers;
			}
		}
	}
	else
	{
		vector<string> const *vp=&sv;	// vp = vector pointer; sv = source vector
		vector<string> lv;				// lv = lower-case vector
		if (!d->rStrg.solve()) { elog.add(d->error4); return pair<bool,bool>(false,false); }
		const string& refstr = *(d->rStrg.str());
		if (d->cis)
		{
			for (each_element(sv,it)) lv.push_back(boost::to_lower_copy(*it));
			vp=&lv;
		}
		if (d->comparison=='k')
		{
			pair<bool,bool> k_result=d->srowbi.match(*vp);
			if (k_result.first)
				if (d->ng1) k_result.second=!k_result.second;
			return k_result;
		}
		for (each_element(d->flds,it))
		{
			bool ers;
			switch (d->comparison)
			{
				case 'e': ers = ((*vp)[*it]==refstr); break;
				case 'g': ers = ((*vp)[*it]>=refstr); break;
				case 'G': ers = ((*vp)[*it]> refstr); break;
				case 'l': ers = ((*vp)[*it]<=refstr); break;
				case 'L': ers = ((*vp)[*it]< refstr); break;
				case 'n': ers = ((*vp)[*it]!=refstr); break;
				case 'm': ers = boost::regex_match( (*vp)[*it], d->what, d->filter ); break;
				case 's': ers = (str_startsw((*vp)[*it],refstr));		break;
				case 'd': ers = (str_endsw((*vp)[*it],refstr));			break;
				case 'h': ers = (str_has((*vp)[*it],refstr));			break;
				case 'i': ers = ( exist_element(d->refdb,(*vp)[*it]));	break;
				case 'o': ers = (!exist_element(d->refdb,(*vp)[*it]));	break;
				case 'S': ers = ((*vp)[*it].empty());					break;
				default: exit_error("single condition: wrong comparison for a non-numeric field.");
			}
			if (d->ng1) ers=!ers;
			if (d->lao) r&=ers; else r|=ers;
		}
	}
	return pair<bool,bool>(true,r);
}

bool single_condition::read (istream & in)
{
	skip_whitespaces(in);
	if (in.peek()!='[') return false;
	initialize();
	char c1,c2;
	
	// read [
	in>>c1;
	
	// read row_selection
	d->rws.read(in);
		
	// read O (operation = / !) (& is default, so must be omitted)
	in>>c1;
	if (c1=='!'||c1=='^') {	d->ng2=true;  in>>c1; }
	if (c1=='/'||c1=='|') {	d->lao=false; in>>c1; }
	
	// read T (type = l n s NF) and N (field numbers)
	if (c1=='N') in>>c1;
	switch (c1)
	{
		case 'l':d->tan=true;	d->tfl=true;	skip_whitespaces(in); d->flds.parse_or_exit(in);	break;
		case 'n':d->tan=true;					skip_whitespaces(in); d->flds.parse_or_exit(in);	break;
		case 's':d->tan=false;					skip_whitespaces(in); d->flds.parse_or_exit(in);	break;
		case 'F':d->tan=true;	d->tnf=true;														break;
		default: exit_error("single condition: wrong type of data "+s(c1)+". Must be l/n/s/NF.");
	}
	
	// read C
	in>>c1;
	if (c1=='!'||c1=='^') {	d->ng1=true; in>>c1; }
	if (c1=='~') {	d->cis=true; in>>c1; }
	c2=in.peek();
	switch (c1) {
		case '!':	if (c2=='=') {	d->comparison='n'; in>>c2; break; } exit_error("single condition: wrong C "+s(c1)+s(c2));	// not equal
		case '^':	if (c2=='=') {	d->comparison='n'; in>>c2; break; } exit_error("single condition: wrong C "+s(c1)+s(c2));	// not equal
		case '<':	if (c2=='=') {	d->comparison='l'; in>>c2; break; }											// less/equal
									d->comparison='L';		   break;											// less than
		case '=':	if (c2=='=') {	d->comparison='e'; in>>c2; break; }											// equal
									d->comparison='e';		   break;											// equal
		case '>':	if (c2=='=') {	d->comparison='g'; in>>c2; break; }											// greater/equal
									d->comparison='G';		   break;	 										// greater than
		case 'e':	if (c2=='q') {	d->comparison='e'; in>>c2; break; }											// equal
					if (c2=='w') {	d->comparison='d'; in>>c2; break; } exit_error("single condition: wrong C "+s(c1)+s(c2));	// end with
		case 'g':	if (c2=='t') {	d->comparison='G'; in>>c2; break; }											// greater than
					if (c2=='e') {	d->comparison='g'; in>>c2; break; } exit_error("single condition: wrong C "+s(c1)+s(c2));	// greater equal
		case 'h':	if (c2=='s') {	d->comparison='h'; in>>c2; break; } exit_error("single condition: wrong C "+s(c1)+s(c2));	// has
		case 'i':	if (c2=='n') {	d->comparison='i'; in>>c2; break; }															// in
					if (c2=='m') {	d->comparison='S'; in>>c2; break; } exit_error("single condition: wrong C "+s(c1)+s(c2));	// is missing
		case 'k':	if (c2=='m') {	d->comparison='k'; in>>c2; break; } exit_error("single condition: wrong C "+s(c1)+s(c2));	// key match
		case 'l':	if (c2=='t') {	d->comparison='L'; in>>c2; break; }															// less than
					if (c2=='e') {	d->comparison='l'; in>>c2; break; } exit_error("single condition: wrong C "+s(c1)+s(c2));	// less equal
		case 'n':	if (c2=='e') {	d->comparison='n'; in>>c2; break; }															// not equal
					if (c2=='i') {	d->comparison='o'; in>>c2; break; } exit_error("single condition: wrong C "+s(c1)+s(c2));	// not in
		case 'r':	if (c2=='m') {	d->comparison='m'; in>>c2; break; } exit_error("single condition: wrong C "+s(c1)+s(c2));	// regex match
		case 's':	if (c2=='w') {	d->comparison='s'; in>>c2; break; } exit_error("single condition: wrong C "+s(c1)+s(c2));	// start with
		default: exit_error("Wrong test: "+s(c1)+s(c2));
	}
	if (d->tan)
	{
		switch (d->comparison) {
			case 'd': exit_error("single condition: comparison 'ew' is not for numeric fields.");
			case 'h': exit_error("single condition: comparison 'hs' is not for numeric fields.");
			case 'i': exit_error("single condition: comparison 'in' is not for numeric fields.");
			case 'k': exit_error("single condition: comparison 'km' is not for numeric fields.");
			case 'o': exit_error("single condition: comparison 'ni' is not for numeric fields.");
			case 'm': exit_error("single condition: comparison 'rm' is not for numeric fields.");
			case 's': exit_error("single condition: comparison 'sw' is not for numeric fields.");
			default: break;
		}
		if (d->cis)	exit_error("single condition: case-insensitive (~) doesn't mean anything for numeric tests.");
	}
	
	// read R (reference)
	for (int i=0;!in.eof();) { // i is the extra [
		c1=in.get();
		if (c1=='[') { ++i; }
		if (c1==']') { if (!i) break; --i; }
		d->input_ref.push_back(c1); }
	if (d->comparison=='S' && !d->input_ref.empty()) exit_error("single condition: is_missing does not allow right operand.");
	
	// prepare
	// if ( d->cis) boost::to_lower(d->input_ref); // don't add this, because sometimes input_ref is a filename (i/o/k) or regex (m)
	// instead, I require users to take care of input_ref themselves -- it has to be in lower case
	
	if ( d->tan && d->comparison!='S')	d->rMath.set_expression(d->input_ref);
	if (!d->tan)						d->rStrg.set_expression(d->input_ref);
	if (d->comparison=='i' || d->comparison=='o')
	{
		if (d->cis)
			for (Rows_in_File(inf,d->input_ref,1)) d->refdb.insert(boost::to_lower_copy(inf[0]));
		else
			for (Rows_in_File(inf,d->input_ref,1)) d->refdb.insert(inf[0]);
	}
	if (d->comparison=='m')	d->filter.assign(d->input_ref, program.regex_syntax);
	if (d->comparison=='k')	d->srowbi.input(d->flds,d->input_ref,d->cis);
	
	// all done, return success
	return true;
}

void single_condition::write_eval(ostream& out) {
	out<<'[';
	d->rws.write(out);
	if (!d->lao) out<<'|';
	if ( d->ng2) out<<'!';
	if ( d->tnf) out<<"NF";
	else
	{
		if		(d->tfl)	out<<'l';
		else if (d->tan)	out<<'n';
		else				out<<'s';
		print_container(d->flds,out,',');
	}
	if (d->ng1) out<<'!';
	if (d->cis) out<<'~';
	switch (d->comparison) {
		case 'd':out<<"ew";break;
		case 'e':out<<"eq";break;
		case 'g':out<<"ge";break;
		case 'h':out<<"hs";break;
		case 'i':out<<"in";break;
		case 'k':out<<"km";break;
		case 'l':out<<"le";break;
		case 'm':out<<"rm";break;
		case 'n':out<<"ne";break;
		case 'o':out<<"ni";break;
		case 's':out<<"sw";break;
		case 'w':out<<"sw";break;
		case 'G':out<<"gt";break;
		case 'L':out<<"lt";break;
		case 'S':out<<"im";break;
		default :exit_error("single condition: wrong C (comparison).");
	}
	out<<d->input_ref;
	out<<']';
}

// ------------- multiple_conditions internal class -------------

namespace {
	class operator_btwn_conditions {
	public:
		bool and_operation;
		bool assigned;		// already asigned &,|
		operator_btwn_conditions():and_operation(true),assigned(false){}
//		void write_oper(ostream& out) {
//			if (and_operation)  out<<"&";
//			else				out<<"|";
//		}
	};
	
	class logic_tree_node : public single_condition, public operator_btwn_conditions {
	public:
		bool eval_type;
		pair<bool,bool> result; // valid,satisfied
		logic_tree_node():eval_type(false),result(true,true) { }
		logic_tree_node(bool t):eval_type(t),result(true,true) { }
//		friend ostream& operator<<(ostream& out, logic_tree_node& ob) {
//			if (ob.eval_type)	ob.write_eval(out);
//			else				ob.write_oper(out);
//			return out;
//		} // rm this func coz it creates a C++11 warning: unused function 'operator<<' [-Wunused-function]
	};
	
	class logic_tree  {
	private:
		bool					df;		// defined
		int				error1,error2;
		pair<bool,bool>			rs;		// result
		vector<logic_tree_node>	nodes;	//
		ntree<int>				str;	// tree structure
		row_selection			rws;
		void _read_nodes(istream& in, ntree<int>::iterator parent)
		{
			while (!in.eof())
			{
				skip_whitespaces(in);
				char c=in.peek();
				switch (c) {
					case '(': {
						in.get();
						ntree<int>::iterator sub = str.append(parent,nodes.size());
						nodes.push_back(logic_tree_node(false));
						_read_nodes(in,sub);
						if ( (c=in.get()) != ')' ) exit_error("mingle conditions: Missing )");
						continue;
					}
					case '+':
					case '&':
						in.get();
						if (!nodes[*parent].and_operation && nodes[*parent].assigned) exit_error("mingle conditions: ..&..|.. not allowed, must be grouped by ().");
						nodes[*parent].and_operation=true;
						nodes[*parent].assigned=true;
						nodes[*parent].result.second=true;
						continue;
					case '/':
					case '|':
						in.get();
						if (nodes[*parent].and_operation && nodes[*parent].assigned) exit_error("mingle conditions: ..|..&.. not allowed, must be grouped by ().");
						nodes[*parent].and_operation=false;
						nodes[*parent].assigned=true;
						nodes[*parent].result.second=false;
						continue;
					case ')':
						return;
					case EOF: return;
					default:
						break;
				}
				str.append(parent,nodes.size());
				nodes.push_back(logic_tree_node(true));
				if (!nodes.back().read(in)) return;
			}
		}
	public:
		logic_tree():df(false),rs(true,true){
			error1=elog.get_token("multiple conditions were not fully tested and result=false");
			error2=elog.get_token("multiple conditions were not fully tested but result=true");
		}
		bool is_valid()	{ return rs.first;  }
		bool is_true()	{ return rs.second; }
		bool is_defined() { return df; }
		pair<bool,bool>result() { return rs; }
		void read(const string& instr)
		{
			std::stringstream ss;
			if (instr.find('[')==string::npos && instr.find(']')==string::npos) ss.str("["+instr+"]");
			else ss.str(instr);
			
			// clear the tree object, set first node = [&]
			ntree<int>::iterator top = str.insert(0);
			nodes.clear();
			nodes.push_back(logic_tree_node(false));
			
			// read row_selection
			rws.read(ss);
			
			// read condition tree
			_read_nodes(ss,top);
			if (str.children(top)==0)										exit_error("Multiple conditions: nothing has been read.");
			if (ss.peek()!=(char)EOF) { string s; getline(ss,s,(char)EOF);	exit_error("Multiple conditions: remaining string "+s); }
			// previously if (!ss.eof()), works fine with g++ any version and clang++ w/ libstdc++, but libc++ has a bug
			// It always return false unless a get() is evoked. This is not right but as of 2013-07-10 clang++ v4.2, this
			// problem still exist. I changed it so that it still works.
			df=true;
		}
		
		void DepthFirstTraverseUp_eval(ntree<int>& t, vector<string>& v)
		{
			for (ntree<int>::postfix_iterator loc=str.postfix_begin(); loc!=str.postfix_end(); ++loc)
			{
				int ths=*loc;
				if (nodes[ths].eval_type) nodes[ths].result=nodes[ths].evaluate(v);
				if (str.depth(loc.simplify())-1) // ths is not root
				{
					int par = * (str.parent(loc.simplify()));
					if (nodes[par].and_operation)	nodes[par].result.second &= nodes[ths].result.second;
					else							nodes[par].result.second |= nodes[ths].result.second;
					nodes[par].result.first &= nodes[ths].result.first;
					// .first is always &=, therefore, root.result.first tells whether all conditions are
					// evaluatable (ie, all fields are there). When root.result.first==false, root.result.second
					// can still be true, because not all conditions are required because of the '|' operation.
				}
			}
		}
		
		void eval(vector<string>& v)
		{
			if (str.empty()) { return; }
			if (nodes.size()<=1) { return; }
			for (each_element(nodes,it))
				if (it->and_operation)	it->result=pair<bool,bool>(true,true);
				else					it->result=pair<bool,bool>(true,false);
			DepthFirstTraverseUp_eval(str,v);
			rs=nodes[*(str.root())].result;
			if (!rs.first) { if (rs.second) elog.add(error2); else elog.add(error1); }
			rws.analyze(rs.second);
		}
		
		void renew() // prv: return whether satisfied line has passed: if (rn==-1) return 0; else return 1;
		{
			rws.renew();
		}
		
		int max_field_num()
		{
			int maxnum=-INT_MAX;
			for (each_element(nodes,loc))
				if (loc->eval_type)
					if (loc->max_field_num()>maxnum) maxnum=loc->max_field_num();
			return maxnum;
		}
		
		bool size_valid(int size)
		{
			bool sv=true;
			for (each_element(nodes,loc))
				if (loc->eval_type) sv&=loc->size_valid(size);
			return sv;
		}
	};
}

// ------------- multiple_conditions -------------

struct multiple_conditions::MultipleConditionData {
	int  grp;		// number of times it's renew-ed
	bool all_f;		// all test should return false if all_f=true, used by --until
	logic_tree t;	// test
	logic_tree u;	// until
	logic_tree r;	// renew
	MultipleConditionData():grp(0),all_f(false){}
};

multiple_conditions::~multiple_conditions() {
	delete d;
}

multiple_conditions::multiple_conditions():d(NULL) {
	d=new MultipleConditionData;
}

multiple_conditions::multiple_conditions(const multiple_conditions& othr):d(NULL) {
	d=new MultipleConditionData;
	*d = *(othr.d);
}

multiple_conditions& multiple_conditions::operator=(const multiple_conditions& orig){
	if(&orig != this)	*(this->d) = *(orig.d);
	return *this;
}

void multiple_conditions::read_arguments(vector<string>& srce_opt)
{
	vector<string> dest_opt;
	int args=srce_opt.size();
	for (int argi=0;argi<args;++argi)
	{
		if (srce_opt[argi]=="--if")
		{
			if (args-argi<2) exit_error("Insufficient arguments for the --if option.");
			d->t.read(srce_opt[++argi]);
		}
		else if (srce_opt[argi]=="--renew")
		{
			if (args-argi<2) exit_error("Insufficient arguments for the --renew option.");
			d->r.read(srce_opt[++argi]);
		}
		else if (srce_opt[argi]=="--until")
		{
			if (args-argi<2) exit_error("Insufficient arguments for the --until option.");
			d->u.read(srce_opt[++argi]);
		}
		else
			dest_opt.push_back(srce_opt[argi]);
	}
	srce_opt=dest_opt;
}

void multiple_conditions::read(const string& s) {
	d->t.read(s);
}

pair<bool,bool> multiple_conditions::evaluate(vector<string>& v) {
	if (d->all_f)	return pair<bool,bool>(true,false);
	if (d->u.is_defined()) { d->u.eval(v); if (d->u.is_true()) { d->all_f=true; return pair<bool,bool>(d->u.is_valid(),false); }}
	if (d->r.is_defined()) { d->r.eval(v); if (d->r.is_true()) { ++d->grp; d->t.renew(); } }
	d->t.eval(v);
	return d->t.result();
}

bool multiple_conditions::is_defined() {
	return d->t.is_defined();
}

bool multiple_conditions::satisfied() {
	pair<bool,bool> rs=evaluate(program.main_data().current_row()); return rs.second;
}

bool multiple_conditions::is_valid()			{ return d->t.is_valid(); }
//bool multiple_conditions::is_true()			{ return d->t.is_true();  }
//pair<bool,bool>multiple_conditions::result()	{ return d->r.result();   }
int multiple_conditions::group_number()			{ return d->grp; }

int multiple_conditions::max_field_num() {
	int maxnum=-INT_MAX, n;
	n=d->t.max_field_num(); if (n>maxnum) maxnum=n;
	n=d->r.max_field_num(); if (n>maxnum) maxnum=n;
	n=d->u.max_field_num(); if (n>maxnum) maxnum=n;
	return maxnum;
}

bool multiple_conditions::size_valid(int size) {
	bool sv=true;
	sv&=d->t.size_valid(size);
	sv&=d->r.size_valid(size);
	sv&=d->u.size_valid(size);
	return sv;
}

int ReadCnd(const std::string& input, multiple_conditions& mc, int ErrCode)
{
	mc.read(input);
	return 1;
}

void ReadArg(const std::vector<std::string>& args, size_t& argi, multiple_conditions& mc, int ErrCode)
{
	std::string input;
	std::size_t found = args[argi].find('=');
	if (found!=std::string::npos)	{ NeedArg(args,0,argi); input = args[argi].substr(found+1); }
	else							{ NeedArg(args,1,argi); input = args[++argi]; }
	ReadCnd(input,mc,ErrCode);
}
