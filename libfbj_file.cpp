#include "libfbj_base.hpp"
#include "libfbj_file.hpp"
#include "libfbj_program.hpp"
#include "libfbj_progress.hpp"
#include <algorithm>
#include <deque>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>

using namespace std;

namespace {

	#define _N_DF_OPT 10

}

// =============== tfile_format ===============

struct tfile_format::tfile_format_data {

	vector<string> _DF_OPT_TEXT={
		"DF_SUCCESSIVE_DELIMITERS_AS_ONE",\
		"DF_QUOTED",\
		"DF_SKIP_COMMENT",\
		"DF_KEEP_COMMENT",\
		"DF_SKIP_BLANK",\
		"DF_READ_BLANK",\
		"DF_TRIM_LWS",\
		"DF_TRIM_TEF",\
		"DF_IRG",\
		"DF_CSV"
	};
	
	vector<string> _DF_OPT_ARGV={
		"s",\
		"quoted",\
		"skip-comment",\
		"keep-comment",\
		"skip-blank",\
		"read-blank",\
		"trim-lws",\
		"trim-tef",\
		"irg",\
		"csv"
	};
	
	vector<bool> default_DF_OPT={
		false,\
		false,\
		true,\
		false,\
		false,\
		false,\
		false,\
		false,\
		false,\
		false,\
	};
	
	vector<bool> default_FB_OPT={
		false,\
		false,\
		false,\
		false,\
		false,\
		false,\
		false,\
		false,\
		false,\
		false,\
	};
	bool default_FB_OPT_t=false;
	bool default_FB_OPT_d=false;

	static int license_chk;			// license checking
	struct fnstruct {
		field_numbers*	fnp;		// field_numbers pointer
		int				tok;		// token for elog
		int				bak;		// backup of token
		action_code		act;		// action if unsatisfied
		fnstruct(field_numbers* f,int t,action_code a):fnp(f),tok(t),bak(t),act(a){}
		// no need for the following:Once explicitly declare any constructor, the compiler stops providing the implicit default constructor.
		// fnstruct(){exit_error("not allowed");}
	};
	vector<fnstruct> fnv;
	field_numbers dfn;				// default field_number object
	RowsBuffer*	hrw;				// Def.=NULL means don't use RowsBuffer
	deque< vector<string> >* vrw;	// Def.=NULL means don't use vector_of_row
	RowNum_t		titlelines_;	// number of title lines
	int   			n_inp_dels_;	// number of delimiters
	unsigned long	progressBar;	// show bar every 1 line, also show % of max number
	unsigned long	progressPer;	// show # read out of total=progressPer lines
	unsigned long	progressPct;	// show % read out of total=progressPct lines
	string input_dels_;				// delimiters string
	string output_del_;				// delimiter string for output
	string msg_nf_dif_;				// message_when_nf_differ
	string cmt_symbol_;				// comment line start with it. Can't be empty; otherwise all lines are comment!
	vector<bool> DF_OPT;			// delimited_file option
	vector<bool> FB_OPT;			// forbid option[i]
	bool FB_OPT_t;					// forbid -t
	bool FB_OPT_d;					// forbid -d
	bool is_default;				// this object is the default_tfile_format
	tfile_format_data():dfn(false,true),hrw(NULL),vrw(NULL) {
		_set_to_original();
		_read_env_var();
		is_default = false;
	}
	void _set_to_original() {
		fnv.clear();
		dfn.clear();
		hrw=NULL;
		vrw=NULL;
		cmt_symbol_="#";
		titlelines_=0;
		progressBar=0;
		progressPer=0;
		progressPct=0;
		_set_inp_dels_to_default();
		msg_nf_dif_="Unequal number of fields between lines in ";
		DF_OPT = default_DF_OPT;
		FB_OPT = default_FB_OPT;
		FB_OPT_t=default_FB_OPT_t;
		FB_OPT_d=default_FB_OPT_d;
	}
	void _set_inp_dels_to_default() { n_inp_dels_=2; input_dels_="\t "; output_del_="\t"; if (is_default) DLMTR = output_del_; }
	void _set_inp_dels(const string&);	// setup input & output delimiters, argument could be an empty string (ie, no delimiters)
	void _set_comb_opt();				// setup combinational options
	void _read_env_var();				// setup options from environmental variable
};

string DLMTR="\t";
int tfile_format::tfile_format_data::license_chk = program.nudge;
tfile_format default_tfile_format(true);

void tfile_format::tfile_format_data::_set_inp_dels(const string& ds) {
	input_dels_=ds;
	n_inp_dels_=ds.length();
	if (n_inp_dels_)	output_del_=s(input_dels_[0]);
	else				output_del_.clear();
	if (is_default) DLMTR = output_del_;
}

void tfile_format::tfile_format_data::_set_comb_opt()
{
	if (DF_OPT[FORMAT_IRREGULAR])
	{
		DF_OPT[SUCCESSIVE_DELIMITERS_AS_ONE]=true;
		DF_OPT[TRIM_LEADING_WHITESPACES]=true;
		DF_OPT[TRIM_LAST_EMPTY_FIELDS]=true;
		DF_OPT[SKIP_BLANKS]=true;
		msg_nf_dif_.clear();
	}
	if (DF_OPT[FORMAT_CSV])
	{
		DF_OPT[SUCCESSIVE_DELIMITERS_AS_ONE]=false;
		DF_OPT[TRIM_LEADING_WHITESPACES]=false;
		DF_OPT[TRIM_LAST_EMPTY_FIELDS]=false;
		DF_OPT[QUOTED]=true;
		_set_inp_dels(",");
	}
	if (DF_OPT[SKIP_BLANKS] && DF_OPT[READ_BLANKS]) exit_error("You can't set both SKIP_BLANKS and READ_BLANKS to yes.");
	if (DF_OPT[KEEP_NOTES]) DF_OPT[SKIP_NOTES]=false;
}

void tfile_format::tfile_format_data::_read_env_var()
{
	// read _DF_OPT
	for (int i=0;i<_N_DF_OPT;++i)
		if (bool_env(_DF_OPT_TEXT[i]).exist) 
			DF_OPT[i]=bool_env(_DF_OPT_TEXT[i]).is_true;
	EnvironmentVariable s = str_env("DF_COMMENT_SYMBOL");	if (s.exist) { if (s.value.empty()) exit_error("DF_COMMENT_SYMBOL can't be empty"); else cmt_symbol_=s.value; }
	EnvironmentVariable t = str_env("DF_HEADER_ROWS");		if (t.exist) { if (t.value.empty()) exit_error("DF_HEADER_ROWS can't be empty"); else titlelines_=boost::lexical_cast<int>(t.value); }
	EnvironmentVariable d = str_env("DF_DELIMITERS");		if (d.exist) { _set_inp_dels(replace_escape_sequence_copy(d.value)); }
	if (bool_env("DF_WARN_NO_NFV").is_true) msg_nf_dif_.clear();
	_set_comb_opt();
}

string tfile_format::help_text() {
	string result="FOR FILE READING\n";
	if (!fdta->FB_OPT_d)							result+="  -d STR / -dSTR   Delimiters (1st char. for output too, '' allowed) {'"+reverse_escape_sequence_copy(fdta->input_dels_)+"'}\n";
	if (!fdta->FB_OPT_t)							result+="  -t INT / -tINT   Number of title lines {"+itos(fdta->titlelines_)+"}\n";
	if (!fdta->FB_OPT[SUCCESSIVE_DELIMITERS_AS_ONE])result+="  -s [B]           Treat successive delimiters as one {"+str_YesOrNo(fdta->DF_OPT[SUCCESSIVE_DELIMITERS_AS_ONE])+"}\n";
	if (!fdta->FB_OPT[QUOTED])						result+=" --quoted=B        Allow quoted fields, quotation marks remain as content {"+str_YesOrNo(fdta->DF_OPT[QUOTED])+"}\n";
	if (!fdta->FB_OPT[SKIP_NOTES])					result+=" --comment=S       Lines starting with S are comments {"+comment_sw()+"}\n";
	if (!fdta->FB_OPT[SKIP_NOTES])					result+=" --skip-comment=B  Skip comment lines (omit and move on) {"+str_YesOrNo(fdta->DF_OPT[SKIP_NOTES])+"}\n";
	if (!fdta->FB_OPT[KEEP_NOTES])					result+=" --keep-comment=B  Keep comment lines (show and move on) {"+str_YesOrNo(fdta->DF_OPT[KEEP_NOTES])+"}\n";
	if (!fdta->FB_OPT[SKIP_BLANKS])					result+=" --skip-blank=B    Skip blank lines {"+str_YesOrNo(fdta->DF_OPT[SKIP_BLANKS])+"}\n";
	if (!fdta->FB_OPT[READ_BLANKS])					result+=" --read-blank=B    Read blank lines if not skipped, {"+str_YesOrNo(fdta->DF_OPT[READ_BLANKS])+"}\n";
	if (!fdta->FB_OPT[TRIM_LEADING_WHITESPACES])	result+=" --trim-lws=B      Skip leading whitespaces {"+str_YesOrNo(fdta->DF_OPT[TRIM_LEADING_WHITESPACES])+"}\n";
	if (!fdta->FB_OPT[TRIM_LAST_EMPTY_FIELDS])		result+=" --trim-tef=B      Skip trailing empty fields {"+str_YesOrNo(fdta->DF_OPT[TRIM_LAST_EMPTY_FIELDS])+"}\n";
	if (!fdta->FB_OPT[FORMAT_IRREGULAR])			result+=" --irg=B           Input is irregular: s=Y trim-lws=Y trim-tef=Y skip-blank=Y\n";
	if (!fdta->FB_OPT[FORMAT_CSV])					result+=" --csv=B           Input is CSV: -d , quoted=Yes s=No trim-lws=No trim-tef=No\n";
													result+=" --progress-bar N  Show progress as % of N lines, display a bar with scale\n";
													result+=" --progress-per N  Show progress as number of lines, display at every N lines\n";
													result+=" --progress-pct N  Show progress as % of N lines\n";
													result+="  -Wno-nfv         Do not show warning: number of fields varies between lines\n";
	int has_opt=0; for (int i=0;i<_N_DF_OPT;++i) if (!fdta->FB_OPT[i]) ++has_opt;
//	if (has_opt) {
//		result+="*  =B can be omitted if B is Yes.\n";
//		result+="   Format: --xx=B / xx=B / --xx (if B=yes) / --xx B. Examples:-s=no / s=no / -s\n";
//		result+="   These defaults are factory defaults adjusted by environment variable values.\n";
//		result+="   Env. var. is named DF_XXX (all upper case, - changed to _) E.g. DF_TRIM_LWS.\n";
//		result+="   Exceptions: -s is DF_SUCCESSIVE_DELIMITERS_AS_ONE, -d is DF_DELIMITERS.\n";
//	}
	return result;
}

tfile_format::~tfile_format() {
	delete fdta;
}

tfile_format::tfile_format():fdta(NULL) {
	fdta=new tfile_format_data;
}

tfile_format::tfile_format(bool is_default):fdta(NULL) {
	fdta=new tfile_format_data;
	fdta->is_default=is_default;
}

tfile_format::tfile_format(const tfile_format& othr):fdta(NULL) {
	fdta=new tfile_format_data;
	*fdta = *(othr.fdta);
}

tfile_format& tfile_format::operator=(const tfile_format& orig) {
	if(&orig != this)	*(this->fdta) = *(orig.fdta);
	return *this; 
}

void tfile_format::reset()
{
	fdta->_set_to_original();
}

void tfile_format::set_option(_DF_OPT_NAME opt_name, bool value)
{
	fdta->DF_OPT[opt_name]=value;
	fdta->_set_comb_opt();
}

void tfile_format::forbid_nf_rpt()
{
	fdta->msg_nf_dif_.clear();
}

void tfile_format::enable_nf_rpt()
{
	fdta->msg_nf_dif_="Unequal number of fields between lines in ";
}

void tfile_format::forbid_option(_DF_OPT_NAME opt_name)
{
	fdta->FB_OPT[opt_name]=true;
}

void tfile_format::forbid_option(const string& opt_name)
{
	if		(opt_name=="-d") fdta->FB_OPT_d=true;
	else if (opt_name=="-t") fdta->FB_OPT_t=true;
	else if (opt_name=="delimiting_related") {
		fdta->FB_OPT_d=true;
		fdta->FB_OPT[SUCCESSIVE_DELIMITERS_AS_ONE]=true;
		fdta->FB_OPT[QUOTED]=true;
		fdta->FB_OPT[TRIM_LEADING_WHITESPACES]=true;
		fdta->FB_OPT[TRIM_LAST_EMPTY_FIELDS]=true;
		fdta->FB_OPT[FORMAT_IRREGULAR]=true;
		fdta->FB_OPT[FORMAT_CSV]=true;
	}
	else exit_error("Unknown argument "+opt_name+" for forbid_option().");
}

int tfile_format::ReadOpt(vector<string>& srce_opt, size_t& argi) {
	int num_read=0;
	for (size_t arg_size=srce_opt.size(); argi<arg_size; ++argi)
	{
		if		(srce_opt[argi]=="-d")
		{
			if (fdta->FB_OPT_d) exit_error("Program option -d is forbidden.");
			if (arg_size-argi<2) exit_error("Insufficient arguments for the -d option.");
			set_delimiters(srce_opt[++argi]);
			++num_read;
		}
		else if (str_startsw(srce_opt[argi],"-d"))
		{
			if (fdta->FB_OPT_d) exit_error("Program option -d is forbidden.");
			set_delimiters(srce_opt[argi].substr(2));
			++num_read;
		}
		else if (str_startsw(srce_opt[argi],"-d="))
		{
			if (fdta->FB_OPT_d) exit_error("Program option -d is forbidden.");
			set_delimiters(srce_opt[argi].substr(3));
			++num_read;
		}
		else if (srce_opt[argi]=="-t")
		{
			if (fdta->FB_OPT_t) exit_error("Program option -t is forbidden.");
			if (arg_size-argi<2) exit_error("Insufficient arguments for the -t option.");
			try { fdta->titlelines_=boost::lexical_cast<int>(srce_opt[++argi]); }
			catch (boost::bad_lexical_cast &) {	exit_error("Error reading NUMBER for the -t option."); }
			++num_read;
		}
		else if (str_startsw(srce_opt[argi],"-t="))
		{
			if (fdta->FB_OPT_t) exit_error("Program option -t is forbidden.");
			try { fdta->titlelines_=boost::lexical_cast<int>(srce_opt[argi].substr(3)); }
			catch (boost::bad_lexical_cast &) {	exit_error("Error reading NUMBER for the -t option."); }
			++num_read;
		}
		else if (str_startsw(srce_opt[argi],"--comment=")) // prv --notes-startw comment-symbol=
		{
			if (fdta->FB_OPT[SKIP_NOTES]) exit_error("Program option --comment is forbidden.");
			fdta->cmt_symbol_=srce_opt[argi].substr(10);
			if (fdta->cmt_symbol_.empty()) exit_error("--comment should not be empty.");
			++num_read;
		}
		else if (str_startsw(srce_opt[argi],"--progress-bar"))
		{
			string parameter;
			if (srce_opt[argi].size()>14 && srce_opt[argi][14]=='=')			parameter=srce_opt[argi].substr(15);
			else { if (arg_size-argi<2) exit_error("Need Arg for "+srce_opt[argi]); parameter=srce_opt[++argi]; }
			parameter.erase (std::remove(parameter.begin(), parameter.end(), ','), parameter.end());
			int m=1;
			size_t l = parameter.size()-1;
			if		(parameter[l]=='h' || parameter[l]=='H') { parameter.pop_back(); m=100; }
			else if	(parameter[l]=='k' || parameter[l]=='K') { parameter.pop_back(); m=1000; }
			else if (parameter[l]=='m' || parameter[l]=='M') { parameter.pop_back(); m=1000000; }
			try { fdta->progressBar = boost::lexical_cast<unsigned long>(parameter) * m; }
			catch (boost::bad_lexical_cast &) {	exit_error("Error reading NUMBER for the --progress-bar option."); }
			if (fdta->progressPer)	 exit_error("You can't use both --progress-bar & --progress-per.");
			if (fdta->progressPct)	 exit_error("You can't use both --progress-bar & --progress-pct.");
			++num_read;
		}
		else if (str_startsw(srce_opt[argi],"--progress-per"))
		{
			string parameter;
			if (srce_opt[argi].size()>14 && srce_opt[argi][14]=='=')			parameter=srce_opt[argi].substr(15);
			else { if (arg_size-argi<2) exit_error("Need Arg for "+srce_opt[argi]); parameter=srce_opt[++argi]; }
			parameter.erase (std::remove(parameter.begin(), parameter.end(), ','), parameter.end());
			int m=1;
			size_t l = parameter.size()-1;
			if		(parameter[l]=='h' || parameter[l]=='H') { parameter.pop_back(); m=100; }
			else if	(parameter[l]=='k' || parameter[l]=='K') { parameter.pop_back(); m=1000; }
			else if (parameter[l]=='m' || parameter[l]=='M') { parameter.pop_back(); m=1000000; }
			try { fdta->progressPer = boost::lexical_cast<unsigned long>(parameter) * m; }
			catch (boost::bad_lexical_cast &) {	exit_error("Error reading NUMBER for the --progress-per option."); }
			if (fdta->progressBar)	 exit_error("You can't use both --progress-per & --progress-bar.");
			if (fdta->progressPct)	 exit_error("You can't use both --progress-per & --progress-pct.");
			++num_read;
		}
		else if (str_startsw(srce_opt[argi],"--progress-pct"))
		{
			string parameter;
			if (srce_opt[argi].size()>14 && srce_opt[argi][14]=='=')			parameter=srce_opt[argi].substr(15);
			else { if (arg_size-argi<2) exit_error("Need Arg for "+srce_opt[argi]); parameter=srce_opt[++argi]; }
			parameter.erase (std::remove(parameter.begin(), parameter.end(), ','), parameter.end());
			int m=1;
			size_t l = parameter.size()-1;
			if		(parameter[l]=='h' || parameter[l]=='H') { parameter.pop_back(); m=100; }
			else if	(parameter[l]=='k' || parameter[l]=='K') { parameter.pop_back(); m=1000; }
			else if (parameter[l]=='m' || parameter[l]=='M') { parameter.pop_back(); m=1000000; }
			try { fdta->progressPct = boost::lexical_cast<unsigned long>(parameter) * m; }
			catch (boost::bad_lexical_cast &) {	exit_error("Error reading NUMBER for the --progress-pct option."); }
			if (fdta->progressBar)	 exit_error("You can't use both --progress-pct & --progress-bar.");
			if (fdta->progressPer)	 exit_error("You can't use both --progress-pct & --progress-per.");
			++num_read;
		}
		else if (srce_opt[argi]=="-Wno-nfv")
		{
			forbid_nf_rpt();
			++num_read;
		}
		else
		{
			bool is_df_opt=false;
			for (int i=0;i<_N_DF_OPT;++i)
			{
				string complete_argv;
				if (fdta->_DF_OPT_ARGV[i].size()==1)	complete_argv = "-"  + fdta->_DF_OPT_ARGV[i];
				else									complete_argv = "--" + fdta->_DF_OPT_ARGV[i];
				
				if (srce_opt[argi]==complete_argv)
				{
					if (fdta->FB_OPT[i]) exit_error("Program option "+complete_argv+" is forbidden.");
					if (arg_size-argi>=2)
					{
						try { fdta->DF_OPT[i]=IsYes(srce_opt[argi+1]); ++argi; }
						catch (input_exception &) { fdta->DF_OPT[i]=true; }
					}
					else fdta->DF_OPT[i]=true;
					fdta->_set_comb_opt();
					is_df_opt=true;
					break;
				}
				if (str_startsw(srce_opt[argi],complete_argv) && srce_opt[argi][complete_argv.size()]=='=')
					srce_opt[argi]=srce_opt[argi].substr(complete_argv.size()-fdta->_DF_OPT_ARGV[i].size());
				if (str_startsw(srce_opt[argi],fdta->_DF_OPT_ARGV[i]) && srce_opt[argi][fdta->_DF_OPT_ARGV[i].size()]=='=')
				{
					if (fdta->FB_OPT[i]) exit_error("Program option "+complete_argv+" is forbidden.");
					string answer=srce_opt[argi].substr(fdta->_DF_OPT_ARGV[i].size()+1);
					try { fdta->DF_OPT[i]=IsYes(answer); }
					catch (input_exception &) { exit_error("Program option "+complete_argv+" argument wrong."); }
					fdta->_set_comb_opt();
					is_df_opt=true;
					break;
				}
			}
			if (is_df_opt)	++num_read;
			else break;
		}
	}
	--argi; // Always point to the option left to the unknown option. Good for "for (argi=1; argi<arg_size; ++argi)".
	return num_read;
}

void tfile_format::read_arguments(vector<string>& srce_opt, size_t start, bool expected, bool stop_at_unknown) {
	vector<string> dest_opt;
	int tot_read = expected; // if (expected) read @ start or ended; else, find the 1st to read.
	for (size_t ended=0, argi=0, arg_size=srce_opt.size(); argi<arg_size; ++argi)
	{
		if (argi<start) { dest_opt.push_back(srce_opt[argi]); continue; }
		if (ended)		{ dest_opt.push_back(srce_opt[argi]); continue; }
		int add_read = ReadOpt(srce_opt,argi);
		if (add_read)
		{
			tot_read += add_read;
			if (tot_read && stop_at_unknown) ended=1;
		}
		else
			dest_opt.push_back(srce_opt[++argi]);
	}
	srce_opt=dest_opt;
}

int ReadArg(std::vector<std::string>& args, size_t& argi, tfile_format& f, int ErrCode){
	size_t option_id = argi;
	int num_read=f.ReadOpt(args,++argi);
	if (num_read) return num_read;
	else if	(ErrCode>0)	{	elog.add(ErrCode); return 0; }
	else if (ErrCode<0)	{	exit_error("cannot read tfile_format options in "+args[option_id]); return 0; }
	else				{	/*elog.add(0) does nothing*/ return 0; }
}

// --------------- setup ---------------

void tfile_format::set_storage_to(RowsBuffer& h)				{ fdta->hrw=&h; fdta->vrw=NULL; }
void tfile_format::set_storage_to(deque< vector<string> >& v)	{ fdta->vrw=&v; fdta->hrw=NULL; }
void tfile_format::set_titlelines(RowNum_t l)					{ fdta->titlelines_=l; }
void tfile_format::set_delimiters_to_default()					{ fdta->_set_inp_dels_to_default(); }
void tfile_format::set_delimiters(const string& ds)				{ fdta->_set_inp_dels(ds); }
void tfile_format::clear_field_nums()							{ fdta->fnv.clear(); }

void tfile_format::set_field_nums(field_numbers& f, const string& msg, action_code action_i) {
	int token=0; // token=0 does nothing
	if (!msg.empty()) token=elog.get_token(msg);
	for (each_element(fdta->fnv,it)) {
		if (it->fnp==&f) {
			it->tok=token;
			it->bak=token;
			it->act=action_i;
			return;
		}
	}
	fdta->fnv.push_back(tfile_format_data::fnstruct(&f,token,action_i));
}

void tfile_format::set_field_nums(int n, const string& msg, action_code action_i) {
	fdta->dfn.clear();
	fdta->dfn.push_back(n);
	int token=0;
	if (!msg.empty()) token=elog.get_token(msg);
	for (each_element(fdta->fnv,it)) {
		if (it->fnp==&(fdta->dfn)) {
			it->tok=token;
			it->bak=token;
			it->act=action_i;
			return;
		}
	}
	fdta->fnv.push_back(tfile_format_data::fnstruct(&(fdta->dfn),token,action_i));
}

void tfile_format::forbid_sv_rpt() { for (each_element(fdta->fnv,it)) it->tok=0; }
void tfile_format::enable_sv_rpt() { for (each_element(fdta->fnv,it)) it->tok=it->bak; }

// --------------- access ---------------
string& tfile_format::output_del()	{ return fdta->output_del_; }
string& tfile_format::comment_sw()	{ return fdta->cmt_symbol_; }
bool tfile_format::opt(_DF_OPT_NAME option_name) { return fdta->DF_OPT[option_name]; }
RowNum_t tfile_format::titlelines() { return fdta->titlelines_; }
int tfile_format::min_NumFields() {
	int r=0;
	for (each_element_const(fdta->fnv,it))
	{
		int m=it->fnp->min_required();
		if (m>r) r=m;
	}
	return r;
}

// --------------- work ---------------

bool tfile_format::is_content(char c) {
	if (is_EndOfLine(c)) return false;
	for (int i=0;i<fdta->n_inp_dels_;++i)	if (c==fdta->input_dels_[i]) return false;
	return true;
}

bool tfile_format::is_delimiter(char c) {
	for (int i=0;i<fdta->n_inp_dels_;++i)	if (c==fdta->input_dels_[i]) return true;
	return false;
}

int tfile_format::act_by_size(int s, bool& validity)  { // return action type
	validity=true;
	int overall=None;
	for (each_element_const(fdta->fnv,it))
	{
		if (!it->fnp->size_valid(s)) 
		{
			elog.add(it->tok);
			if (it->act > overall) overall=it->act;
			validity=false;
		}
	}
	return overall;
}

// =============== begin tabular_file_data ===============

struct tabular_file::tabular_file_data {
	bool						_SizeValid;		// whether the line has a valid size based on tfile_format* format
	RowNum_t					row_num;		// row number, 1st row is 0.
	int							file_num;		// file number, 1st file is 0.
	int							skip_begin;		// skip the first few lines
	int							skip_every;		// after the 1st line read, skip this number of lines before each next
	int							line_len;		// max length of a line
	vector<string>*				cont;			// content pointer, init=NULL but not after read_row / read_line, which has at least 1 field
	tfile_format				my_format;		// internally stored format, used only when the default format is used
	tfile_format*				format;			// the format to apply, I use a pointer so that user can change format during reading
	boost::iostreams::filtering_istream file;	// the file to be read
	boost::progress_display_fbj pbar;			// progress bar
	progress_time				ppct;			// progress pct
	vector<string>				my_cont;		// internally stored contents
	deque<string>				filenames;		// filenames[0]=currently opened file; [1-]=to be opened; empty()=nothing opened.
	map<int,int>				nf_count;		// num of fields count, for error output only (diff NF between rows)
	tabular_file_data():_SizeValid(false),row_num(-1),file_num(-1),skip_begin(0),skip_every(0),line_len(1024),cont(NULL),my_format(default_tfile_format),format(&my_format) {}
	void skip_rows(int n);						// skip n rows. Used only at the line begining, because it deals with comments and blank lines.
	bool read_line();							// read one line to str and then _skip_one_line
	bool read_row();							// read one row.
	bool eof();									// the currently opened file is at eof()
	bool is_blank();							// the last line is a blank line (contents().size=1 && contents()[0].empty())
	void open_file(const string& filename);		// only modify file;   require skip_begin & skip_every be setup already!
	bool open_next();							// open the next file; require skip_begin & skip_every be setup already!
};

void tabular_file::tabular_file_data::open_file(const string& filename)
{
	file.reset();
	if (is_stdin(filename))
		file.push(std::cin);
	else {
		string name_cvt;
		if (!find_file_zippedOrNot(filename,name_cvt)) exit_cannotFind(filename);
		if		(str_endsw(name_cvt,".gz"))		file.push(boost::iostreams::gzip_decompressor());
		else if (str_endsw(name_cvt,".bz2"))	file.push(boost::iostreams::bzip2_decompressor());
		else if (str_endsw(name_cvt,".zip"))	exit_error(".zip file not supported.");
		else if (str_endsw(name_cvt,".Z"))		exit_error(".Z file not supported.");
		else if (str_endsw(name_cvt,".7z"))		exit_error(".7z file not supported.");
		file.push(boost::iostreams::file_source(name_cvt));
		if (!(file.component<boost::iostreams::file_source>(file.size()-1))->is_open()) exit_cannotOpen(name_cvt);
	}
	if (skip_begin) skip_rows(skip_begin);
}

bool tabular_file::tabular_file_data::open_next() {
	if (filenames.empty()) exit_error("open_next err");	// should be called only when one file is already open!
	while (filenames.size()>1) {						// open only when there's a next
		filenames.pop_front(); file.reset();			// close the current file
		open_file(filenames.front());					// open next file using old parameters, exit if not open
		if (eof()) continue;							// file empty, try next again
		++file_num;
		row_num=-1;
		return true;									// next file is open and contains data
	}
	return false;										// there's no next
}

bool tabular_file::tabular_file_data::eof() {
	return file.peek()==EOF;
}

bool tabular_file::tabular_file_data::is_blank() {
	return cont->size()==1 && cont->front().empty();
	// didn't test cont->empty() because it won't happen; otherwise "program.test df.grep" won't pass.
}

// Use at the line begining because it deals with comments and blank lines.
void tabular_file::tabular_file_data::skip_rows(int n)
{
	string str;
	str.reserve(line_len);
	for (int i=0; i<n && !eof(); )
	{
		safeGetline(file, str);
		int size = str.size();
		if (size>line_len) line_len=size;
		if (format->opt(TRIM_LEADING_WHITESPACES)) boost::trim_left(str);
		if (format->opt(SKIP_NOTES) && str_startsw(str, format->fdta->cmt_symbol_)) { continue; }
		if (format->opt(KEEP_NOTES) && str_startsw(str, format->fdta->cmt_symbol_)) { program.outf << str << endl; continue; }
		if (format->opt(SKIP_BLANKS) && str.empty()) continue;
		++i;
	}
}

bool tabular_file::tabular_file_data::read_line() {
	cont->assign(1,string());
	string& str(cont->front());
	str.reserve(line_len);
	while (!eof()) // no need to skip_current_row because \n is extracted
	{
		safeGetline(file, str); // \n extracted.
		int size = str.size();
		if (size>line_len) line_len=size;
		if (format->opt(TRIM_LEADING_WHITESPACES)) boost::trim_left(str);
		if (format->opt(SKIP_NOTES) && str_startsw(str, format->fdta->cmt_symbol_)) { continue; }
		if (format->opt(KEEP_NOTES) && str_startsw(str, format->fdta->cmt_symbol_)) { program.outf << str << endl; continue; }
		if (format->opt(SKIP_BLANKS) && str.empty()) continue;
		break;
	}
	if (format->opt(SKIP_BLANKS) && str.empty()) { if (open_next()) return read_line(); else return false; } // stop at eof (prv no this line but works)
	if (is_blank() && !format->opt(READ_BLANKS)) { if (open_next()) return read_line(); else return false; } // stop at blank
	++row_num; // add "if (!str_startsw(cont->front(), format->fdta->cmt_symbol_))" on 2015-07-06 but removed on 2015-09-10, otherwise run_num could be -1
	if (format->fdta->progressBar) ++pbar;
	if (format->fdta->progressPct) ++ppct;
	if (format->fdta->progressPer && row_num%format->fdta->progressPer==0) cerr<<row_num<<'\r';
	return true;
}

bool tabular_file::tabular_file_data::read_row() {
	std::istream::sentry se(file, true);
	std::streambuf* sb = file.rdbuf();
	for (cont->clear(); !eof(); cont->clear())
	{
		int col=0;
		cont->assign(1,string());
		string* str = &cont->back();
		bool AfterDel=false; // if c is delimiter, it's not the first
		for(;;) {
			char c = sb->sbumpc();
			switch (c) {
				case '\n':
					goto end_reading_a_line;
				case '\r':
					if (sb->sgetc() == '\n') sb->sbumpc();
					goto end_reading_a_line;
				case EOF:
					file.setstate(std::ios::eofbit);
					goto end_reading_a_line;
				case ' ':
				case '\t':
					if (format->opt(TRIM_LEADING_WHITESPACES) && !col && str->empty()) continue;
				default:
					if (format->is_delimiter(c))
					{
						if (format->opt(SUCCESSIVE_DELIMITERS_AS_ONE) && AfterDel) continue;
						cont->resize(cont->size()+1);
						str = &cont->back();
						++col;
						AfterDel=true;
					}
					else
					{
						if (format->opt(QUOTED) && c=='\"' && AfterDel)
						{
							str->push_back(c);
							for(;;)
							{
								c = sb->sbumpc();
								if (c==EOF) exit_error("Closing quotation mark not found before end-of-file.");
								str->push_back(c);
								if (c=='\"')
								{
									if (sb->sgetc() == '\"') sb->sbumpc();
									else break;
								}
							}
						}
						else
						{
							str->push_back(c);
						}
						AfterDel=false;
					}
			}
		}
end_reading_a_line:
		// handle comments if it is a comment line
		if (format->opt(SKIP_NOTES))
			if (str_startsw(cont->front(), format->fdta->cmt_symbol_)) { continue; }
		if (format->opt(KEEP_NOTES))
			if (str_startsw(cont->front(), format->fdta->cmt_symbol_)) { print_container(*cont, program.outf, format->fdta->output_del_, true); continue; }
		
		// remove trailing empty fields
		if (format->opt(TRIM_LAST_EMPTY_FIELDS))
			while (cont->size()>1 && cont->back().empty()) cont->pop_back();
		
		// skip a blank line or return a blank line, whether to take it depends on the caller
		if (format->opt(SKIP_BLANKS) && is_blank()) continue;
		else break;
	}
	if (cont->empty()) { cont->assign(1,string()); if (open_next()) return read_row(); else return false; } // stop at eof
	if (is_blank() && !format->opt(READ_BLANKS)) { if (open_next()) return read_row(); else return false; } // stop at blank
	++nf_count[cont->size()];
	++row_num; // add "if (!str_startsw(cont->front(), format->fdta->cmt_symbol_))" on 2015-07-06 but removed on 2015-09-10, otherwise run_num could be -1
	if (format->fdta->progressBar) ++pbar;
	if (format->fdta->progressPct) ++ppct;
	if (format->fdta->progressPer && row_num%format->fdta->progressPer==0) cerr<<row_num<<'\r';
	return true;
}

// =============== end tabular_file_data ===============

// =============== begin tabular_file ===============
// Make sure: eof() should be together with open_next(). Therefore, I made a function eof(); always use it instead of d->eof().

void tabular_file::close()
{
	if (is_open())
	{
		if (tdta->format->fdta->progressPer) cerr<<"# "<<tdta->row_num+1<<" lines read.                                                            \n";
		if (tdta->format->fdta->progressBar) cerr<<"# "<<tdta->row_num+1<<" lines read.                                                            \n";
		if (tdta->format->fdta->progressPct) cerr<<"# "<<tdta->row_num+1<<" lines read.                                                            \n";
		if (tdta->nf_count.size()>1 && !tdta->format->fdta->msg_nf_dif_.empty())
		{
			stringstream ss;
			for (each_element(tdta->nf_count,itn)) ss<<' '<<itn->first<<'('<<itn->second<<')';
			elog.add(tdta->format->fdta->msg_nf_dif_+tdta->filenames.front()+" --"+ss.str(),false);
		}
		tdta->file.reset();
		tdta->filenames.clear();
		tdta->nf_count.clear();
		tdta->row_num=-1;
		tdta->cont=NULL;
	}
}

tabular_file::~tabular_file()
{
	close();
	delete tdta;
}

void tabular_file::open(const string& filename,int n, int skip_begin_n, int skip_every_m)
{
	close();
	if (n) {
		if (tdta->format != &(tdta->my_format)) exit_error("Tabular_file point to external format while you want to chagne a internal parameter.");
		tdta->my_format.set_field_nums(n,"line(s) in '"+filename+"' do not have sufficient fields.",tfile_format::Continue); }
	tdta->skip_begin=skip_begin_n;
	tdta->skip_every=skip_every_m;
	tdta->open_file(filename);
	tdta->filenames.assign(1,filename);
	if (tdta->format->fdta->progressBar) tdta->pbar.restart(tdta->format->fdta->progressBar);
	if (tdta->format->fdta->progressPct) tdta->ppct.restart(tdta->format->fdta->progressPct);
	++tdta->file_num;
	tdta->row_num=-1;
}

tabular_file::tabular_file(const string& filename,int n, int skip_begin_n, int skip_every_m):tdta(NULL)
{
	tdta=new tabular_file_data;
	open(filename,n,skip_begin_n,skip_every_m);
}

void tabular_file::open(const string& filename,tfile_format* f, int skip_begin_n, int skip_every_m)
{
	close();
	if (f) tdta->format=f;
	tdta->skip_begin=skip_begin_n;
	tdta->skip_every=skip_every_m;
	tdta->open_file(filename);
	tdta->filenames.assign(1,filename);
	if (tdta->format->fdta->progressBar) tdta->pbar.restart(tdta->format->fdta->progressBar);
	if (tdta->format->fdta->progressPct) tdta->ppct.restart(tdta->format->fdta->progressPct);
	++tdta->file_num;
	tdta->row_num=-1;
}

tabular_file::tabular_file(const string& filename,tfile_format* f, int skip_begin_n, int skip_every_m):tdta(NULL)
{
	tdta=new tabular_file_data;
	open(filename,f,skip_begin_n,skip_every_m);
}

template<typename StrContainr>
tabular_file::tabular_file(const StrContainr& namelist,int n, int skip_begin_n, int skip_every_m):tdta(NULL)
{
	tdta=new tabular_file_data;
	open(namelist,n,skip_begin_n,skip_every_m);
}
template tabular_file::tabular_file(const vector<string>& namelist, int n, int b, int e);
template tabular_file::tabular_file(const  deque<string>& namelist, int n, int b, int e);
template tabular_file::tabular_file(const    set<string>& namelist, int n, int b, int e);

template<typename StrContainr>
tabular_file::tabular_file(const StrContainr& namelist,tfile_format* f, int skip_begin_n, int skip_every_m):tdta(NULL)
{
	tdta=new tabular_file_data;
	open(namelist,f,skip_begin_n,skip_every_m);
}
template tabular_file::tabular_file(const vector<string>& namelist, tfile_format* f, int b, int e);
template tabular_file::tabular_file(const  deque<string>& namelist, tfile_format* f, int b, int e);
template tabular_file::tabular_file(const    set<string>& namelist, tfile_format* f, int b, int e);

tabular_file::tabular_file():tdta(NULL)
{
	tdta=new tabular_file_data;
}

tabular_file::tabular_file(const tabular_file& othr):tdta(NULL) {
//	tdta=new tabular_file_data;
//	*tdta = *(othr.tdta);
	exit_error("Copy constructor is prohibited for tabular_file.");
}

tabular_file& tabular_file::operator=(const tabular_file& orig) {
//	if(&orig != this)
//	{
//		*(this->tdta)	= *(orig.tdta);
//		*(this->format) = *(orig.format);
//	}
//	this won't compile coz iostream cannot be copied, tabular_file::operator= doesn't make sense anyway
	exit_error("Assignment operation is prohibited for tabular_file.");
	return *this;
}

string&			tabular_file::operator[](int p) { if (p>=0) return (*(tdta->cont))[p]; else return (*(tdta->cont))[tdta->cont->size()+p]; }
vector<string>& tabular_file::contents()		{ return *(tdta->cont); }
RowNum_t		tabular_file::RowNumber()		{ return tdta->row_num; }
int				tabular_file::FileNumber()		{ return tdta->file_num; }
int				tabular_file::NumFields()		{ return tdta->cont->size(); }
bool			tabular_file::SizeValid()		{ return tdta->_SizeValid; }
bool			tabular_file::eof()				{ return tdta->eof() && !tdta->open_next(); }
bool			tabular_file::is_open()			{ return tdta->filenames.size(); }
bool			tabular_file::is_blank()		{ return tdta->is_blank(); }
bool			tabular_file::is_header()		{ return RowNumber() < tdta->format->titlelines(); }
bool			tabular_file::is_comment()		{ return str_startsw(tdta->cont->front(), tdta->format->fdta->cmt_symbol_); }
void			tabular_file::clear_nf()		{ tdta->nf_count.clear(); }
void			tabular_file::seekg(long o)		{ tdta->file.seekg(o); }
string			tabular_file::FileName()		{ if (tdta->filenames.empty()) exit_error("no file opened"); return tdta->filenames[0]; }

void tabular_file::add_file(const string& filename) {
	if (tdta->filenames.empty())
	{
		tdta->open_file(filename);
		tdta->filenames.assign(1,filename);
		if (tdta->format->fdta->progressBar) tdta->pbar.restart(tdta->format->fdta->progressBar);
		if (tdta->format->fdta->progressPct) tdta->ppct.restart(tdta->format->fdta->progressPct);
	}
	else
		tdta->filenames.push_back(filename);
}

bool tabular_file::next() {
	if (tdta->skip_every) tdta->skip_rows(tdta->skip_every);
	if (eof()) return false;
	return true;
}

bool tabular_file::read_r() {
	if (eof()) return false;
	RowsBuffer*					hrw = tdta->format->fdta->hrw;
	deque< vector<string> >*	vrw = tdta->format->fdta->vrw;
	if		(hrw) tdta->cont=&(hrw->push_back());
	else if (vrw) { vrw->resize(vrw->size()+1); tdta->cont=&(vrw->back()); }
	else tdta->cont=&(tdta->my_cont);
	if (!tdta->read_row()) { if (hrw) hrw->pop_back(); if (vrw) vrw->pop_back(); return false; }
	for (;;)
	{
		int action=tdta->format->act_by_size(NumFields(),tdta->_SizeValid);
		switch (action) {
			case tfile_format::None:	return true;
			case tfile_format::Expand:{	int mr=tdta->format->min_NumFields();
										vec_deq_set(*(tdta->cont),mr-1,s(""));}
										return true;
			case tfile_format::Wr_Cont: write_r(program.outf,true);
			case tfile_format::Continue:if (!(next()&&tdta->read_row())) { if (hrw) hrw->pop_back(); if (vrw) vrw->pop_back(); return false; }
										continue;
			case tfile_format::Break:	return false;
			case tfile_format::Exit:	exit_error("Insufficient number of columns in the file "+FileName());
			default:					exit_error("Undefined action for tfile_format class");
		}
	}
	exit_error("pleae check tabular_file::read_r().");
	return true; // should not happen
}

bool tabular_file::read_l() {
	if (eof()) return false;
	RowsBuffer*					hrw = tdta->format->fdta->hrw;
	deque< vector<string> >*	vrw = tdta->format->fdta->vrw;
	if		(hrw)	tdta->cont=&(hrw->push_back());
	else if (vrw) { vrw->resize(vrw->size()+1); tdta->cont=&(vrw->back()); }
	else			tdta->cont=&(tdta->my_cont);
	if (!tdta->read_line()) { if (hrw) hrw->pop_back(); if (vrw) vrw->pop_back(); return false; }
	return true;
}

void tabular_file::write_r(ostream& ostr,bool write_endl) {
	print_container(*tdta->cont, ostr, tdta->format->fdta->output_del_, write_endl);
}

void tabular_file::row2str(string& ostr,bool write_endl) {
	ostr.clear();
	string& del = tdta->format->fdta->output_del_;
	vector<string>& val = *tdta->cont;
	if (!val.empty())
	{
		vector<string>::iterator it(val.begin());
		ostr+=*it;
		for (++it; it!=val.end(); ++it) { ostr+=del; ostr+=*it; }
	}
	if (write_endl) ostr.push_back('\n');
}

void tabular_file::write_r(ostream & ostr, vector<int>& f, bool write_endl) {
	if (!f.empty()) {
		vector<int>::iterator it(f.begin());
		ostr << (*tdta->cont)[*it] ;
		for (++it; it!=f.end(); ++it) ostr << tdta->format->fdta->output_del_ << (*tdta->cont)[*it] ;
	}
	if (write_endl) ostr<<endl; // previously '\n' but it doesn't flush and severely delays the display in stdout
}

void tabular_file::write_r(ostream & ostr, field_numbers& f, bool write_endl) {
	try 
	{
		f.contents_to_ostream(*(tdta->cont),ostr,tdta->format->fdta->output_del_,write_endl);
	}
	catch (bad_query_fieldnumbers_unsolved &e)	{ elog.add("lines not printed due to unsolved negative field numbers."); }
	catch (bad_query_fieldnumbers_absent &e)	{ elog.add("lines not printed due to absent must-have fields."); }
}

