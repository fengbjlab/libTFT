#include <map>
#include <list>
#include <cstdint>				// for        [u]int8/16/32/64_t
#include <boost/date_time/local_time/local_time.hpp>
#include <boost/regex.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/interprocess/sync/file_lock.hpp> // for boost::interprocess::file_lock
#include <boost/thread/thread.hpp>

#include "libfbj_program.hpp"
#include "libfbj_base.hpp"
#include "libfbj_security.hpp"
#include "_all_txt.hpp"	// need _package_dir _package_vendor _package_version _program_name _program_manual _program_history _program_finale

using namespace std;

bool bad_curl(string GotStr)
{
	if (GotStr.empty()) return true; // if the website is down, it may return nothing without any error message
	boost::to_lower(GotStr);
	if (str_has(GotStr,"400") && str_has(GotStr,"bad request")) return true;
	if (str_has(GotStr,"401") && str_has(GotStr,"authorization required")) return true;
	if (str_has(GotStr,"403") && str_has(GotStr,"permission denied")) return true;
	if (str_has(GotStr,"404") && str_has(GotStr,"not found")) return true;
	return false;
}

// use the curl command, which has the -L option, good for redirected URLs
// use two URLs, in case one URL is down.
// return 0:not_checked 1:up_to_date -1:not_up_to_date
int chk_net_vers(const string& url1, const string& url2, const string& this_version, string& read_version) // previously bool is_expired_by_net()
{
	if (url1.empty() || this_version.empty()) return 0;
	if (linux_command_which_dir("curl").empty()) exit_error("cannot find the program curl.");
	try { read_version = exec("curl -sL "+url1+" 2>/dev/null",true); } // "2>/dev/null" hide the progress, but also any error messages
	catch (const std::exception& error) { exit_error("failed to get the latest version number"); } // impossible because of the true above
	if (bad_curl(read_version))
	{
		if (!url2.empty())
		{
			try { read_version = exec("curl -sL "+url2+" 2>/dev/null",true); } // "2>/dev/null" hide the progress, but also any error messages
			catch (const std::exception& error) { exit_error("failed to get the latest version number"); } // impossible because of the true above
		}
		if (bad_curl(read_version)) exit_error("Cannot obtain the current version.");
	}
	if (read_version != this_version) return -1;
	return 1;
}

// --------------- objects ---------------

// Previously contruct program at last, now before everything because elog & lns require program.
// construct elog before lns because I want to destruct lns before elog.
// after lns is destructed, '\n' is printed at the end, so that elog can continue to print

const boost::gregorian::date	compile_date	 (boost::gregorian::from_simple_string(macro_date_to_boost(__DATE__)));
const std::string				compile_date_str (boost::gregorian::to_simple_string(compile_date));

string _package_license =
"You may not decompile, reverse-engineer, disassemble, or otherwise attempt to derive the source code for the Software. \
You may not modify the Software or create any derivative work of the Software or its accompanying documentation. \
Derivative works include but are not limited to translations. You may not alter any files or libraries in any portion of the Software. \
\n\n\
THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, \
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, \
TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE \
BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT \
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.";

ProgramHandle program;
ErrorLogger elog;
logger lns;

// http://stackoverflow.com/questions/81870/print-variable-type-in-c

template <typename T>
class type_name {
public:
    static const char *name;
};

#define DECLARE_TYPE_NAME(x) template<> const char * type_name<x>::name = #x;
#if __cplusplus >= 201103L
	#define GET_TYPE_NAME(x) (type_name<__typeof__(x)>::name)
#elif defined __GNUG__
	#define GET_TYPE_NAME(x) (type_name<typeof(x)>::name)
#else
	#error This library requires either a GNU C++ compiler or C++0x/C++11 features.
#endif

// end from http://stackoverflow.com/questions/81870/print-variable-type-in-c
// see some DECLARE_TYPE_NAME(x) in libfbj_program.cpp, or put your own in main.cpp before main()

DECLARE_TYPE_NAME(string);
DECLARE_TYPE_NAME(char);
//DECLARE_TYPE_NAME(short);
//DECLARE_TYPE_NAME(int);
//DECLARE_TYPE_NAME(long);
//DECLARE_TYPE_NAME(long long);
//DECLARE_TYPE_NAME(unsigned);
DECLARE_TYPE_NAME(float);
DECLARE_TYPE_NAME(double);
DECLARE_TYPE_NAME(long double);
DECLARE_TYPE_NAME(int8_t);
DECLARE_TYPE_NAME(int16_t);
DECLARE_TYPE_NAME(int32_t);
DECLARE_TYPE_NAME(int64_t);
DECLARE_TYPE_NAME(uint8_t);
DECLARE_TYPE_NAME(uint16_t);
DECLARE_TYPE_NAME(uint32_t);
DECLARE_TYPE_NAME(uint64_t);
#ifdef __APPLE__
DECLARE_TYPE_NAME(size_t);
#endif

template <typename T>
int ReadStr(const std::string& input, T& output, int error_code)
{
	try { output=boost::lexical_cast<T>(input); return 1; }
	catch (boost::bad_lexical_cast &)
	{
		if		(error_code>0)	{	elog.add(error_code); return 0; }
		else if (error_code<0)	{	exit_error("Cannot read \""+input+"\" as "+GET_TYPE_NAME(output)); }
		else					{	/*elog.add(0)=Do_nothing*/ return 0; }
	}
	return -1; // never happens
}
template int ReadStr<string>		(const std::string& input, string&		output, int error_code);
template int ReadStr<char>			(const std::string& input, char&		output, int error_code);
//template int ReadStr<short>		(const std::string& input, short&		output, int error_code);
//template int ReadStr<int>			(const std::string& input, int&			output, int error_code);
//template int ReadStr<long>		(const std::string& input, long&		output, int error_code);
//template int ReadStr<long long>	(const std::string& input, long long&	output, int error_code);
//template int ReadStr<unsigned>	(const std::string& input, unsigned&	output, int error_code);
template int ReadStr<float>			(const std::string& input, float&		output, int error_code);
template int ReadStr<double>		(const std::string& input, double&		output, int error_code);
template int ReadStr<long double>	(const std::string& input, long double& output, int error_code);
template int ReadStr<int8_t>		(const std::string& input, int8_t&		output, int error_code);
template int ReadStr<int16_t>		(const std::string& input, int16_t&		output, int error_code);
template int ReadStr<int32_t>		(const std::string& input, int32_t&		output, int error_code);
template int ReadStr<int64_t>		(const std::string& input, int64_t&		output, int error_code);
template int ReadStr<uint8_t>		(const std::string& input, uint8_t&		output, int error_code);
template int ReadStr<uint16_t>		(const std::string& input, uint16_t&	output, int error_code);
template int ReadStr<uint32_t>		(const std::string& input, uint32_t&	output, int error_code);
template int ReadStr<uint64_t>		(const std::string& input, uint64_t&	output, int error_code);
#ifdef __APPLE__
template int ReadStr<size_t>		(const std::string& input, size_t&		output, int error_code);
#endif

int ReadStr(const std::string& input, bool& output, int error_code)
{
	try { output=IsYes(input); return 1; }
	catch (input_exception &)
	{
		if		(error_code>0)	{	elog.add(error_code); return 0; }
		else if (error_code<0)	{	exit_error("Cannot read \""+input+"\" as bool"); }
		else					{	/*elog.add(0)=Do_nothing*/ return 0; }
	}
	return -1; // never happens
}

int ReadStr(const std::string& input, double& output, int error_code)
{
	// try { output=boost::lexical_cast<double>(input); if (!std::isnan(output)) return 1; else throw "read nan";}
	// if ( input=="NaN" ) { output = std::numeric_limits<double>::signaling_NaN(); return 1; }
	try { output=boost::lexical_cast<double>(input); return 1; }
	catch (...) // boost::bad_lexical_cast &
	{
		if		(error_code>0)	{	elog.add(error_code); return 0; }
		else if (error_code<0)	{	exit_error("Cannot read \""+input+"\" as "+GET_TYPE_NAME(output)); }
		else					{	/*elog.add(0)=Do_nothing*/ return 0; }
	}
	return -1; // never happens
}

int ReadStr(const std::string& input, IntRanges& output, int error_code)
{
	string remains = input;
	int something_read = output.parse(remains);
	if (something_read && remains.empty()) return 1;
	if		(error_code>0)	{	elog.add(error_code); return 0; }
	else if (error_code<0)	{	exit_error("Cannot read \""+input+"\" as IntRanges"); return 0; }
	else					{	/*elog.add(0)=Do_nothing*/ return 0; }
}

void NeedArg(const std::vector<std::string>& args, int need,int argi)
{
	if (((int)args.size()-argi)<(need+1))
	{
		exit_error("Insufficient arguments for the "+args[argi]+" option.");
	}
}

void ReadArg(const std::vector<std::string>& args, size_t& argi, bool& value, int ErrCode)
{
	if (argi>=args.size()) exit_error("argi is out of boundary");
	std::string input;
	std::size_t found = args[argi].find('=');
	if		(found!=std::string::npos)	{ input = args[argi].substr(found+1); }
	else if ((args.size()-argi)==1)		{ input = "Yes"; }
	else if (IsBool(args[argi+1]))		{ input = args[++argi]; }
	else								{ input = "Yes"; }
	ReadStr(input,value,ErrCode);
}

template <typename Container>
int files_matched(const std::string& match_str,Container& v,const std::string& dir)
{
	v.clear();
	boost::filesystem::path p(boost::filesystem::current_path());
	if (!dir.empty()) p = boost::filesystem::path(dir);
	if (!boost::filesystem::is_directory(p)) exit_error(dir+" is not a directory.");
	const boost::regex my_filter(match_str, program.regex_syntax);
	boost::filesystem::directory_iterator i(p);
	for (boost::filesystem::directory_iterator end_itr; i!=end_itr; ++i)
	{
		boost::smatch what;
		if( !boost::filesystem::is_regular_file(i->status()) ) continue;	// Skip if not a file
		if( !boost::regex_match( i->path().filename().string(), what, my_filter ) ) continue;
		add_to_container(v,i->path().string()); // previously i->path().filename().string()
	}
	return v.size();
}
template int files_matched(const std::string& match_str,multiset<string>& v,const std::string& dir);
template int files_matched(const std::string& match_str,  vector<string>& v,const std::string& dir);
template int files_matched(const std::string& match_str,   deque<string>& v,const std::string& dir);
template int files_matched(const std::string& match_str,     set<string>& v,const std::string& dir);

template <typename Container>
int files_matched_zipOrNot(const std::string& match_str,Container& v,const std::string& dir)
{
	v.clear();
	boost::filesystem::path p(boost::filesystem::current_path());
	if (!dir.empty()) p = boost::filesystem::path(dir);
	if (!boost::filesystem::is_directory(p)) exit_error(dir+" is not a directory.");
	const boost::regex my_filter(match_str+"(\\.gz|\\.bz2|)", program.regex_syntax);
	boost::filesystem::directory_iterator i(p);
	for (boost::filesystem::directory_iterator end_itr; i!=end_itr; ++i)
	{
		boost::smatch what;
		if( !boost::filesystem::is_regular_file(i->status()) ) continue;	// Skip if not a file
		if( !boost::regex_match( i->path().filename().string(), what, my_filter ) ) continue;
		std::string filename = i->path().string();  // previously i->path().filename().string()
		if (str_endsw(filename,".gz")) filename.resize(filename.size()-3);
		else if (str_endsw(filename,".bz2")) filename.resize(filename.size()-4);
		add_to_container(v,filename);
	}
	return v.size();
}
template int files_matched_zipOrNot(const std::string& match_str,multiset<string>& v,const std::string& dir);
template int files_matched_zipOrNot(const std::string& match_str,  vector<string>& v,const std::string& dir);
template int files_matched_zipOrNot(const std::string& match_str,   deque<string>& v,const std::string& dir);
template int files_matched_zipOrNot(const std::string& match_str,     set<string>& v,const std::string& dir);

// --------------- ProgramData struct -----------------

#include <mutex>	// g++47 needs it for <thread>! clang++4.0/g++44 doesn't.
#include <thread>	// std::this_thread::sleep_for
#include <chrono>	// std::chrono::seconds

struct ProgramHandle::ProgramData
{
	static int instantiated;// 0=not instantiated 1=instantiated
	
	// run time information
	string exe_filename;	// name of the program at running (user may changed the filename)
	string program_path;	// path to the executable program
	string package_path;	// path to the executable programs as detected by xxx.activate
	string initial_path;	// working dir at the time of entry to main()
	string configr_path;	// path to the configure files, should end with /, previously=_license_path
	string commands_log;	// log of all commands, one option per line
	string lic_fullpath;	// full path to license file
	string pdf_fullpath;	// full path to .pd
	string WebVersion_1;	// up-to-date version number URL
	string WebVersion_2;	// up-to-date version number URL
	string latestVerNum;	// up-to-date version number
	bool   latestVerWeb;	// up-to-date version checking by web

	// user's information
	string home_path;		// home dir = HOME()

	// license holder's information, read from ~/<package>/l
	string first_name;		
	string last_name;		
	string organization;	
	string email_address;	
	string license_version;	
	string fst_use_dt_str;
	string lst_use_dt_str;
	string expiration_str;
	string license_home;
	int	   days_usage;		// tiems usage of the program
	
	// program data
	vector<string>		conf_files;
	vector<string>		arguments;
	set<string>			forbid_opt;
	RowsBuffer			main_data;
	deque<string>		help_texts;
	map<string,string>	help_replc;
	bool				help_requested; // whether arguments has -h
	bool				chk_help;		// if there's -h, print help and exit when reading arg
	bool				arg_read;
	bool				echo_cmd;
	string				outf_name;
	string				outf_pref;
	
	// for progress showing
	std::thread		progress;
	static int		spin_per;
	static int		tickStop;
	static void tick()
	{
		string symbol = "-\\|/"; // "<^>v"
		for (int i=0;!tickStop;++i,i%=4)
		{
			std::this_thread::sleep_for (std::chrono::milliseconds(spin_per));
			cerr << symbol[i] << '\b';
		}
	}
	
	ProgramData():latestVerWeb(true),help_requested(false),chk_help(false),arg_read(false),echo_cmd(false) //,outf_pref(_program_name)
	{
		if (instantiated++) 
		{
			std::cerr<<"ProgramData cannot be instantiated more than once.\n";
			exit(1);
		}
		putenv(strdup("LC_ALL=C"));	// Fix the problem of not running in other machines. strdup() is just to get rid off a warning.
		user_details();				// get HOME_path, exe_filename, program_path, etc.; and get configr_path

#if defined _LIBRARY_DISTR
		forbid_opt.insert("--tft-history");
		forbid_opt.insert("--tft-EULA");
#endif
		forbid_opt.insert("--regex-syntax");
		forbid_opt.insert("--nt");
		forbid_opt.insert("--prefix");
		forbid_opt.insert("--no-web");
		forbid_opt.insert("--spin-per");
	}
	void user_details();
	void replace_text(string& text);
	string _help_text();
	string _copy_right();
};

int ProgramHandle::ProgramData::instantiated = 0;
int ProgramHandle::ProgramData::spin_per = 0;
int ProgramHandle::ProgramData::tickStop = 0;

void ProgramHandle::ProgramData::user_details()
{
	home_path = HOME();
	program_path = self_path();
	exe_filename = substr_after_rfind(program_path,"/");
	program_path = trim_after_rfind(program_path,"/");
	configr_path = home_path+_package_dir;
	package_path = linux_command_which_dir(_package_shortname+string(".activate"));
	if		(FileExists(program_path+_package_shortname+".pd"))	pdf_fullpath = program_path+_package_shortname+".pd";
	else if (FileExists(package_path+_package_shortname+".pd"))	pdf_fullpath = package_path+_package_shortname+".pd";
	else														pdf_fullpath = configr_path+_package_shortname+".pd";
}

bool remove_dir_if_empty(boost::filesystem::path dir)
{
	string dir_string = dir.string();
	if (dir_string.back()=='/') dir_string.pop_back();
	boost::filesystem::path to_remove(dir_string+"/.DS_Store");
	if (boost::filesystem::exists(to_remove)) boost::filesystem::remove(to_remove);
	if (boost::filesystem::is_empty(dir)) { boost::filesystem::remove(dir); return true; }
	else {return false;}
}

string ProgramHandle::ProgramData::_copy_right() {
	string output;
	string date(__DATE__);
	output = "LibTFT version:" + s(_package_version) + " build:" + compile_date_str + '\n';
	output+="(C) 2010-"+date.substr(7)+" "+_package_vendor+".\n";
	return output;
}

string ProgramHandle::ProgramData::_help_text() {
	string result="BASIC OPTIONS\n";
	if (!exist_element(forbid_opt, "-o"))				result += "  -o FILE          output to FILE (.gz/.bz2 allowed) {stdout}\n"; // \"\"=cout=stdout
	if (!exist_element(forbid_opt, "-q"))				result += "  -q / --quiet     run program quietly\n";
	if (!exist_element(forbid_opt, "-h"))				result += "  -h / --help      show this help text\n";
	if (!exist_element(forbid_opt, "--version"))		result += " --version         show version number and check the latest version from the web\n";
	if (!exist_element(forbid_opt, "--nt"))				result += " --nt S            set the number of threads. S may be I, auto, auto/I, auto-I. {1}\n";
	if (!exist_element(forbid_opt, "--prefix"))			result += " --prefix S        set output file prefix. {"+string(_program_name)+"}\n";
	//if (!exist_element(forbid_opt, "--no-web"))		result += " --no-web          do not check for updates\n";
	if (!exist_element(forbid_opt, "--regex-syntax"))	result += " --regex-syntax S  set regular expression syntax (basic/extended/perl/sed/grep/egrep/awk) {basic}\n";
	if (!exist_element(forbid_opt, "--spin-per"))		result += " --spin-per D      show a progress wheel spinning per D second\n";
	if (!exist_element(forbid_opt, "--tft-history"))	result += " --tft-history     show LibTFT's update history\n";
	if (!exist_element(forbid_opt, "--tft-EULA"))		result += " --tft-EULA        show LibTFT's End User License Agreement\n";
	if (!exist_element(forbid_opt, "--tft-credits"))	result += " --tft-credits     show LibTFT's Acknowledgement statement\n";
	return result;
}

void ProgramHandle::ProgramData::replace_text(string& text)
{
	boost::algorithm::replace_all(text,"__PROGRAM_NAME__",string(_program_name));
	boost::algorithm::replace_all(text,"__PACKAGE_VERSION__",string(_package_version));
	boost::algorithm::replace_all(text,"__PACKAGE_DIR__",string(_package_dir));
	boost::algorithm::replace_all(text,"__PACKAGE_VENDOR__",string(_package_vendor));
	for (auto &i:help_replc) boost::algorithm::replace_all(text,i.first,i.second);
}

// --------------- ProgramHandle class -----------------

ProgramHandle::ProgramHandle():quiet(false),nt(1),nudge(1),d(NULL) {
	d=new ProgramData;
	regex_syntax=boost::regex::basic;
	string	 regex_var=str_env("LIBTFT_REGEX_SYNTAX").value;
	if		(regex_var.empty())		;
	else if	(regex_var=="extended")	regex_syntax=boost::regex::extended;
	else if (regex_var=="perl")		regex_syntax=boost::regex::perl;
	else if (regex_var=="basic")	regex_syntax=boost::regex::basic;
	else if (regex_var=="sed")		regex_syntax=boost::regex::sed;
	else if (regex_var=="grep")		regex_syntax=boost::regex::grep;
	else if (regex_var=="egrep")	regex_syntax=boost::regex::egrep;
	else if (regex_var=="awk")		regex_syntax=boost::regex::awk;
	else							exit_error("LIBTFT_REGEX_SYNTAX should be basic / extended / perl / sed / grep / egrep / awk. Currently it is "+regex_var);
	// if (bool_env("LIBTFT_ECHO_COMMAND").exist) d->echo_cmd=bool_env("LIBTFT_ECHO_COMMAND").is_true;
	if (!openOutFile(outf,"")) { std::cerr<<"Cannot open stdout.\n"; exit(1); }
	manual = string(_program_manual);
}

ProgramHandle::ProgramHandle(const ProgramHandle& othr):d(NULL) {
	std::cerr<<"copy constructor is prohibited for ProgramHandle !"<<endl;
	std::cout<<"copy constructor is prohibited for ProgramHandle !"<<endl;
	std::clog<<"copy constructor is prohibited for ProgramHandle !"<<endl;
	exit(1);
//	d=new ProgramData;
//	*d = *(othr.d);
}

ProgramHandle& ProgramHandle::operator=(const ProgramHandle& orig) {
	std::cerr<<"assignment operator is prohibited for ProgramHandle !"<<endl;
	std::cout<<"assignment operator is prohibited for ProgramHandle !"<<endl;
	std::clog<<"assignment operator is prohibited for ProgramHandle !"<<endl;
	exit(1);
//	if(&orig != this)	*(this->d) = *(orig.d);
//	return *this;
}

ProgramHandle::~ProgramHandle() {
	if (d->spin_per) { d->tickStop=1; d->progress.join(); }
	delete d; }

string	ProgramHandle::name()			{ return _program_name;   }
string	ProgramHandle::prefix()			{ return d->outf_pref;    }
string	ProgramHandle::commands()		{ return d->commands_log; }
string	ProgramHandle::out_name()		{ return d->outf_name;    }
string	ProgramHandle::exe_name()		{ return d->exe_filename; }
string	ProgramHandle::exe_path()		{ return d->program_path; }
string	ProgramHandle::init_path()		{ return d->initial_path; }
string	ProgramHandle::conf_path()		{ return d->configr_path; }
RowsBuffer& ProgramHandle::main_data()	{ return d->main_data;    }
void	ProgramHandle::set_check_url(const string& url1, const string& url2, const string& s)
										{	d->WebVersion_1=url1;
											d->WebVersion_2=url2;
											d->latestVerNum=s;
											/*d->forbid_opt.erase("--no-web");*/ }
void	ProgramHandle::set_no_web(const bool noweb)		{ d->latestVerWeb=!noweb;  }
bool	ProgramHandle::no_web()		{ return !d->latestVerWeb;  }
void	ProgramHandle::set_prefix(const string& s)		{ d->outf_pref = s;        }
void	ProgramHandle::forbid_option(const string& s)	{ d->forbid_opt.insert(s); }
void	ProgramHandle::enable_option(const string& s)	{ d->forbid_opt.erase(s);  }
void	ProgramHandle::help_text_var(const std::string& Fr, const std::string& To) { d->help_replc[Fr]=To; }
void	ProgramHandle::check_help_at_arg(){ d->chk_help=true; }

vector<string>& ProgramHandle::arg() {
	if (!d->arg_read) { std::cerr<<"Please call program.read_arguments() before .arg()."<<endl; exit(1); }
	return d->arguments;
};

void ProgramHandle::read_arguments(int argc, char * const argv[], bool rpl_esc, bool MustHaveArguments) {
	if (d->arg_read) exit_error("program.read_arguments() can be called for only once.");
	if (lns.lock_set())
	{
		lns<<showl<<"Program started:";
		for (int argi=0;argi<argc;++argi) { lns<<' '<<argv[argi]; }
		lns<<flush_logger;
	}
	d->arg_read=true;
	if (d->echo_cmd)
	{
		cerr << '\n' << to_term("__THICK_LINE__");
		cerr << "### ";
		for (int argi=0;argi<argc;++argi)
		{
			vector<string> tokens;
			boost::split(tokens, argv[argi], boost::is_any_of(" \t\n\r\0"));
			if (argi) cerr<<' ';
			if (tokens.size()>1)	cerr<<"'"<<argv[argi]<<"'";
			else					cerr<<argv[argi];
		}
		cerr << '\n';
	}
	if (d->exe_filename.empty()) d->exe_filename=substr_after_rfind(argv[0],"/");
	if (d->program_path.empty()) 
	{
		if (str_startsw(argv[0],"/")) d->program_path=trim_after_rfind(argv[0],"/");
		else d->program_path = linux_command_which_dir(d->exe_filename);
	}
	namespace fs = boost::filesystem;
	fs::path full_path( fs::initial_path<fs::path>() ); // fs::initial_path() = fs::current_path() at the time of entry to main()
	d->initial_path=full_path.string();
	for (int argi=0;argi<argc;++argi) d->arguments.push_back(argv[argi]);
	process_arguments(rpl_esc);
	if (argc==1 && MustHaveArguments) { print_help_text(std::cout); exit(0); }
}

void ProgramHandle::show_version(std::ostream& os) {
	os<<to_term("__DOUBLE_LINE__")<<"PROGRAM : "<<_program_name<<' '<<trademark<<endl;
	if (d->latestVerWeb)
	{
		string read_version;
		int chk_net_vers_result = chk_net_vers(d->WebVersion_1,d->WebVersion_2,d->latestVerNum,read_version);
		if (!quiet)
		{
			if (chk_net_vers_result==-1)
				exit_error("This version ("+d->latestVerNum+") is retired.\n" +(str_startsw(read_version," ") ? read_version.substr(1) : ("New version is "+read_version)) );
			else if (chk_net_vers_result==1)
				os<<"# "<<name()<<" version is up-to-date."<<endl;
		}
	}
}

void ProgramHandle::process_arguments(bool rpl_esc) {
	vector<string> bkup = d->arguments;
	d->arguments.clear();
	stringstream all_commands;
	string command_line;
	string nt_input;
	string rg_input;
	double spin_per=0;
	int skip=0;
	int argc=bkup.size();
	for (int argi=0;argi<argc;++argi)
	{
		// write to all_commands, which can be written to a command log file
		string oristr=bkup[argi];
		bool is_opt = str_startsw(oristr,"-") && oristr!="-" ;
		if ( oristr.find('\'')!=string::npos )
		{
			boost::algorithm::replace_all(oristr, "\"", "\\\""); // can't do it recursively, otherwise dead lock
			boost::algorithm::replace_all(oristr, "$", "\\$");	 // same
			boost::algorithm::replace_all(oristr, "`", "\\`");	 // same
			boost::algorithm::replace_all(oristr, "\\", "\\\\"); // same
			oristr = "\"" + oristr + "\"";
		}
		else if ((oristr.find(' ') !=string::npos)|| \
				 (oristr.find('`') !=string::npos)|| \
				 (oristr.find('$') !=string::npos)|| \
				 (oristr.find('\"')!=string::npos)|| \
				 (oristr.find('~') !=string::npos)|| \
				 (oristr.find('&') !=string::npos)|| \
				 (oristr.find('\\')!=string::npos)|| \
				 (oristr.find('|') !=string::npos)|| \
				 (oristr.find(';') !=string::npos)|| \
				 (oristr.find('<') !=string::npos)|| \
				 (oristr.find('>') !=string::npos)|| \
				 (oristr.find('#') !=string::npos)|| \
				 (oristr.find('(') !=string::npos)|| \
				 (oristr.find(')') !=string::npos)|| \
				 (oristr.find('?') !=string::npos)|| \
				 (oristr.find('*') !=string::npos)|| \
				 (oristr.empty())) oristr = "'" + oristr + "'";
		if (is_opt)
		{
			if (!command_line.empty()) all_commands << command_line << " \\" << '\n';
			command_line = oristr + " ";
		}
		else 
		{
			command_line = command_line + oristr + " ";
		}
		
		// prepare arg for return
		if (skip) { --skip; continue; } // for options that has arguments, needs to skip this but should write to command line.
		string rplstr;
		if (rpl_esc)	rplstr=replace_escape_sequence_copy(bkup[argi]);
		else			rplstr=bkup[argi];
		if (rplstr=="--help") rplstr="-h";
		if (rplstr=="--quiet") rplstr="-q";
		if (exist_element(d->forbid_opt, rplstr) || exist_element(d->forbid_opt, substr_before_find(rplstr,"="))) d->arguments.push_back(rplstr);
		else if	(rplstr=="-q")				{	quiet=true; }
		else if	(rplstr=="--no-web")		{	d->latestVerWeb=false; }
		else if	(rplstr=="-h")				{	d->help_requested=true; if (d->chk_help) {print_help_text(std::cout); exit(0);} }
		else if	(rplstr=="--version")		{	show_version(std::cout); 							exit(0); }
		else if (rplstr=="--tft-history")	{	print_messages(std::cout,_program_history);			exit(0); }
		else if (rplstr=="--tft-EULA")		{	print_messages(std::cout,_package_license);			exit(0); }
		else if (rplstr=="--tft-credits")	{	print_messages(std::cout,_package_acknowledgement);	exit(0); }
		else if (rplstr=="-o")				{	if (argc-argi<2) exit_error("lack arg for -o");				skip=1; re_open_output(bkup[argi+1]); }
		else if (rplstr=="--prefix")		{	if (argc-argi<2) exit_error("lack arg for --prefix");		skip=1; ReadStr(bkup[argi+1],d->outf_pref); }
		else if (rplstr=="--spin-per")		{	if (argc-argi<2) exit_error("lack arg for --spin-per");		skip=1; ReadStr(bkup[argi+1],spin_per); }
		else if (rplstr=="--nt") 			{	if (argc-argi<2) exit_error("lack arg for --nt");			skip=1; nt_input=bkup[argi+1]; }
		else if (rplstr=="--regex-syntax")	{	if (argc-argi<2) exit_error("lack arg for --regex-syntax");	skip=1; rg_input=bkup[argi+1]; }
		else if (str_startsw(rplstr,"-q="))				{	ReadStr(bkup[argi].substr(3),quiet);		}
		else if (str_startsw(rplstr,"--quiet="))		{	ReadStr(bkup[argi].substr(8),quiet);		}
		else if (str_startsw(rplstr,"-o="))				{	re_open_output(bkup[argi].substr(3));		}
		else if (str_startsw(rplstr,"--prefix="))		{	ReadStr(bkup[argi].substr(9),d->outf_pref);	}
		else if (str_startsw(rplstr,"--spin-per="))		{	ReadStr(bkup[argi].substr(11),spin_per);	}
		else if (str_startsw(rplstr,"--nt=")) 			{	nt_input=bkup[argi].substr(5);				}
		else if (str_startsw(rplstr,"--regex-syntax=")) {	rg_input=bkup[argi].substr(15);				}
		else d->arguments.push_back(rplstr);
	}
	all_commands << command_line << '\n';
	d->commands_log = all_commands.str();
	if (spin_per) { d->spin_per=spin_per*1000; d->progress=std::thread(ProgramHandle::ProgramData::tick); }
	if (!nt_input.empty())
	{
		if (nt_input=="auto")					{ nt=std::thread::hardware_concurrency()-1; if (nt<1) nt=1; }
		else if	(str_startsw(nt_input,"auto/"))	{ ReadStr(nt_input.substr(5),nt); nt=std::thread::hardware_concurrency()/nt; if (nt<1) exit_error("not enough cores"); }
		else if (str_startsw(nt_input,"auto-"))	{ ReadStr(nt_input.substr(5),nt); nt=std::thread::hardware_concurrency()-nt; if (nt<1) exit_error("not enough cores"); }
		else if	(str_startsw(nt_input,"/"))		{ ReadStr(nt_input.substr(1),nt); nt=std::thread::hardware_concurrency()/nt; if (nt<1) exit_error("not enough cores"); }
		else if (str_startsw(nt_input,"-"))		{ ReadStr(nt_input.substr(1),nt); nt=std::thread::hardware_concurrency()-nt; if (nt<1) exit_error("not enough cores"); }
		else									{ ReadStr(nt_input,nt); if (nt<1) exit_error("Argument for --nt should be greater than or equal to 1"); }
	}
	if (!rg_input.empty())
	{
		if		(rg_input=="extended")	regex_syntax=boost::regex::extended;
		else if (rg_input=="perl")		regex_syntax=boost::regex::perl;
		else if (rg_input=="basic")		regex_syntax=boost::regex::basic;
		else if (rg_input=="sed")		regex_syntax=boost::regex::sed;
		else if (rg_input=="grep")		regex_syntax=boost::regex::grep;
		else if (rg_input=="egrep")		regex_syntax=boost::regex::egrep;
		else if (rg_input=="awk")		regex_syntax=boost::regex::awk;
		else							exit_error("--regex-syntax argument should be basic / extended / perl / sed / grep / egrep / awk. Your input is "+rg_input);
	}
}

void ProgramHandle::re_open_output(const string& s)
{
	outf.reset();
	d->outf_name=s;
	if (!openOutFile(outf,s)) exit_cannotOpen(s);
}

void ProgramHandle::check_help_request() {
	if (d->help_requested) { print_help_text(std::cout); exit(0); }
}

void ProgramHandle::print_help_text(std::ostream& os) {
	if (!manual.empty()) push_back_help(manual,true);
	if (!d->arg_read) exit_error("Please call program.read_arguments() before program.print_help_text().");
	os<<to_term("__DOUBLE_LINE__")<<"PROGRAM : "<<_program_name<<' '<<trademark<<endl;
	for (each_element(d->help_texts,it))
	{
		os<<to_term("__SINGLE_LINE__");
		string& msg = *it;
		d->replace_text(msg);
		os<<to_term(msg);
	}
	string msg;
	msg=d->_help_text();  if (!msg.empty()) { d->replace_text(msg); os<<to_term("__SINGLE_LINE__")<<to_term(msg); }
	msg=_program_finale;  if (!msg.empty()) { d->replace_text(msg); os<<to_term("__SINGLE_LINE__")<<to_term(msg); }
	msg=d->_copy_right(); if (!msg.empty()) { d->replace_text(msg); os<<to_term("__SINGLE_LINE__")<<to_term(msg); }
	os<<to_term("__DOUBLE_LINE__");
}

void ProgramHandle::print_messages(std::ostream& os, std::string msg)
{
	if (!d->arg_read) exit_error("Please call program.read_arguments() before program.print_messages().");
	d->replace_text(msg);
	os<<to_term("__DOUBLE_LINE__");
	os<<"PROGRAM : "<<_program_name<<' '<<trademark<<endl;
	os<<to_term(msg);
	os<<to_term("__DOUBLE_LINE__");
}

void ProgramHandle::push_back_help(string text, bool up_front) {
	if (str_has(text,"ProgramHandle, please enable regex."))
	{
		enable_option("--regex-syntax");
		boost::algorithm::erase_all(text,"ProgramHandle, please enable regex.");
	}
	if (up_front)	d->help_texts.push_front(text);
	else			d->help_texts.push_back (text);
}
