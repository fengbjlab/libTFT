#include <string>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <boost/date_time/local_time/local_time.hpp>
#include <boost/interprocess/sync/file_lock.hpp>
#include <chrono>
#include <thread>
#include "libfbj_base.hpp"
#include "libfbj_log.hpp"
#include "libfbj_program.hpp" // for license checking

using namespace std;

flush_logger_	flush_logger;
fatal_			fatal;
showline_		showl;
showerr_		showe;
showwarning_	showw;
cshowline_		cshowl;
writeline_		writel;
writeerr_		writee;
writewarning_	writew;

//------------------------------ logger data ----------------------------

struct logger::LoggerData
{
	static int		license_chk;						// licensing of the class
	enum			counts_lvl {ths_l,all_l};			// this level, all levels
	enum			output_cat { sh_n_wr, wr_only };	// show and write, write only
	int				counts[100][2][_NUM_LOGGER_CAT];	// [level] [counts_lvl] [counts_cat], include both output_cat
	int				cnt_sh[100][2][_NUM_LOGGER_CAT];	// [level] [counts_lvl] [counts_cat], output_cat: sh_n_wr
	int				level;								// level of output, starts from 0
	int				on_scr;								// whether to put on screen
	std::string		prefix;								// prefix in each line output
	std::fstream	logf;								// log file
	boost::posix_time::ptime time_begin;				// begin time of program running
	boost::posix_time::ptime time_end;					// end   time of program running

	// interprocess mutex
	std::string						mutex_fn;
	boost::interprocess::file_lock	flock;

	LoggerData() { initiate(); }
	void		initiate();
	int			sub_wr();
	int			sub_sh();
	void		increase_count(int c_cat, int o_cat);
	void		set_prefix();
	void		lock();
	void		unlock();
};

int logger::LoggerData::license_chk = program.nudge;

void logger::LoggerData::initiate()
{
	level=0;
	on_scr=1;
	prefix="";
	memset(counts,'\0', sizeof(counts));
	memset(cnt_sh,'\0', sizeof(cnt_sh));
}

int logger::LoggerData::sub_wr()
{
	return counts[level][ths_l][err]+counts[level][ths_l][wrn]+counts[level][ths_l][line];
}

int logger::LoggerData::sub_sh()
{
	return cnt_sh[level][ths_l][err]+cnt_sh[level][ths_l][wrn]+cnt_sh[level][ths_l][line];
}

void logger::LoggerData::increase_count(int c_cat, int o_cat)
{
	if		(o_cat==sh_n_wr)
	{
		int i;
		for (i=0;i<level;i++)	cnt_sh[i][all_l][c_cat]++;
		cnt_sh[i][all_l][c_cat]++;
		cnt_sh[i][ths_l][c_cat]++;
	}
	else if (o_cat==wr_only)
	{
		int i;
		for (i=0;i<level;i++)	counts[i][all_l][c_cat]++;
		counts[i][all_l][c_cat]++;
		counts[i][ths_l][c_cat]++;
	}
	else exit_error("wrong output_category");
}

void logger::LoggerData::set_prefix()
{
	prefix.clear();
	for (int i=0;i<level;i++)	prefix+="   ";
}

void logger::LoggerData::lock()
{
	if (!mutex_fn.empty())
	{
		flock.lock();
	}
}

void logger::LoggerData::unlock()
{
	if (!mutex_fn.empty())
	{
		flock.unlock();
	}
}

//------------------------------ logger ----------------------------

logger::logger():d(NULL)
{
	d=new LoggerData;
	string cerr_mutex_file_name = str_env("STDERR_MUTEX").value;
	if (!cerr_mutex_file_name.empty())
	{
		if (!FileExists(cerr_mutex_file_name))
		{
			std::ofstream mutex_file(cerr_mutex_file_name);
			mutex_file.close();
		}
		set_lock(cerr_mutex_file_name);
	}
}
logger::logger(const std::string& fn):d(NULL)	{ d=new LoggerData; open(fn,std::ios::out); }
logger::logger(const logger& othr):d(NULL)		{ exit_error("The logger object cannot be copied."); }
logger& logger::operator=(const logger& orig)	{ exit_error("The logger object cannot be copied."); return *this; }
logger::~logger() { close(0); delete d; }
void logger::rewind_timer() { d->time_begin=boost::posix_time::microsec_clock::local_time(); }
void logger::set_lock(const std::string& fn)
{
	d->mutex_fn=fn;
	if (!fn.empty()) d->flock=boost::interprocess::file_lock(fn.c_str());
}
bool logger::lock_set()
{
	return !d->mutex_fn.empty();
}

void logger::open(std::string logfilename, const std::ios_base::openmode& mode)
{
	flush_log();
	if (d->logf.is_open()) close(0);
	d->initiate();
	filename_change_home_path(logfilename);
	d->logf.open(logfilename.c_str(),mode);
	if (!d->logf.is_open()) exit_error("cannot open file "+logfilename);
	d->time_begin=boost::posix_time::microsec_clock::local_time();
	d->logf	<<"## Start  running at: " << d->time_begin << "  " << endl;
}

void logger::close(int is_fatal)
{
	// flush. can't rely on endsub() because if level=0, it won't be called.
	flush_log();
	
	// back to level 0
	while (d->level) endsub();
	
	// write summary
	d->lock();
	if (!program.quiet)
	{
		if (d->cnt_sh[0][LoggerData::all_l][err]) cerr<<getpid()<<": "<<"There are "<<d->cnt_sh[0][LoggerData::all_l][err]<<" error(s)!!!\n";
		if (d->cnt_sh[0][LoggerData::all_l][wrn]) cerr<<getpid()<<": "<<"There are "<<d->cnt_sh[0][LoggerData::all_l][wrn]<<" warning(s) ! \n";
	}
	if (is_fatal)
	{
		cerr<<getpid()<<": "<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
		cerr<<getpid()<<": "<<"!!!!! Program exit prematurely due to a fatal error !!!!!\n";
		cerr<<getpid()<<": "<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	}
	d->unlock();
	
	// logf
	if (d->logf.is_open())
	{
		d->logf<<"## There are "<<d->counts[0][LoggerData::all_l][err]<<" error(s)!!!\n";
		d->logf<<"## There are "<<d->counts[0][LoggerData::all_l][wrn]<<" warning(s) ! \n";		
		d->time_end=boost::posix_time::microsec_clock::local_time();
		boost::posix_time::time_duration td = d->time_end - d->time_begin;
		d->logf	<<"## Start  running at: " << d->time_begin << '\n'
				<<"## Finish running at: " << d->time_end << '\n'
				<<"## Duration         : " << td << '\n'
				<<"## ===========================================================================\n\n";
		d->logf.close();
	}
	
	// clear data
	d->initiate();
}

void logger::flush_log()
{
	if (str().size())
	{
		if (d->on_scr && !program.quiet) { d->lock(); cerr<<str()<<endl; d->unlock(); }
		if (d->logf.is_open()) { d->logf<<str()<<endl; d->logf.flush(); }
		str("");
	}
}

void logger::test()
{
	d->lock();
	if (str().size() && d->on_scr && !program.quiet) { cerr<<str()<<endl; str(""); }
	std::this_thread::sleep_for(std::chrono::seconds(10));
	d->unlock();
}

void logger::fatal()
{
	close(1);
	exit(1);
}

void logger::sub()
{
	flush_log();
	++d->level;
	d->set_prefix();
	return;
}

// Previously all_l is ths_l in this function. But don't know which strategy is more reasonable.
// ths_l: errors in lower level are not show, which doesn't make sense
// all_l: errors in lower level are shown even they are not offspring, but cousins
// but from the two, all_l make more sense. It can be interpreted as erros up to this level.
void logger::endsub()
{
	flush_log();
	if (d->counts[d->level][LoggerData::all_l][err] || d->counts[d->level][LoggerData::all_l][wrn])
	{
		if (d->logf.is_open())
			d->logf<<d->prefix<<"("<<d->counts[d->level][LoggerData::all_l][err]<<" errors, "<<d->counts[d->level][LoggerData::all_l][wrn]<<" warnings)"<<endl;
		if (!program.quiet)
		{
			d->lock();
			cerr<<getpid()<<": "<<d->prefix<<"("<<d->counts[d->level][LoggerData::all_l][err]<<" errors, "<<d->counts[d->level][LoggerData::all_l][wrn]<<" warnings)\n";
			d->unlock();
		}
	}
	if (d->level)
	{
		memset(d->counts[d->level],'\0',sizeof(int)*_NUM_LOGGER_CAT*2);
		memset(d->cnt_sh[d->level],'\0',sizeof(int)*_NUM_LOGGER_CAT*2);
		--d->level;
		d->set_prefix();
	}
}

std::string logger::logger_input_string()
{
	flush_log();
	string instr=input_line();
	if (d->logf.is_open()) d->logf<<" : "<<instr;
	return instr;
}

double logger::logger_input_double()
{
	double result;
	flush_log();
	input_number(result);
	if (d->logf.is_open()) d->logf<<" : "<<result;
	return result;
}

char logger::logger_input_char()
{
	flush_log();
	char c=input_char();
	if (d->logf.is_open()) d->logf<<" : "<<c;
	return c;
}

void logger::write(int t)
{
	flush_log();
	if (!d->logf.is_open()) return;
	switch (t) {
		case err:	d->logf<<d->prefix<<"!!! Error:   ";	break;
		case wrn:	d->logf<<d->prefix<<" !  Warning: ";	break;
		case line:	d->logf<<d->prefix;			break;
		default:	exit_error("wrong parameter.");		return;
	};
	d->increase_count(t,LoggerData::wr_only);
	d->on_scr=0;
	return;
}

// changed cerr.flush(); system("tput rev"); cerr<<"Err"; cerr.flush(); system("tput sgr0"); to cerr<<"Err";
// If command is redirected to a file (cmd [opt] > file) , tput writes some funny characters to file.
void logger::show(int t)
{
	write(t);
	switch (t) {
		case err:	*this<<getpid()<<": "<<d->prefix<<"!!! Error:   ";	break;
		case wrn:	*this<<getpid()<<": "<<d->prefix<<" !  Warning: ";	break;
		case line:	*this<<getpid()<<": "<<d->prefix;			break;
		default:	exit_error("wrong parameter.");							return;
	};
	d->increase_count(t,LoggerData::sh_n_wr);
	d->on_scr=1;
	return;
}

void logger::continue_write(int t)
{
	flush_log();
	if (!d->logf.is_open()) return;
	switch (t) {
		case err:						d->logf<<d->prefix<<"!!! Error:   ";	d->increase_count(t,LoggerData::wr_only); break;
		case wrn:						d->logf<<d->prefix<<" !  Warning: ";	d->increase_count(t,LoggerData::wr_only); break;
		case line:	if (d->sub_wr()) {	d->logf<<d->prefix;			d->increase_count(t,LoggerData::wr_only);}break;
		default:	exit_error("wrong parameter.");			return;
	};
	d->on_scr=0;
	return;
}

void logger::continue_show(int t)
{
	continue_write(t);
	switch (t) {
		case err:						*this<<getpid()<<": "<<d->prefix<<"!!! Error:   ";	d->increase_count(t,LoggerData::sh_n_wr);	break;
		case wrn:						*this<<getpid()<<": "<<d->prefix<<" !  Warning: ";	d->increase_count(t,LoggerData::sh_n_wr);	break;
		case line:	if (d->sub_sh()) {	*this<<getpid()<<": "<<d->prefix;			d->increase_count(t,LoggerData::sh_n_wr);}	break;
		default:	exit_error("wrong parameter.");			return;
	};
	d->on_scr=1;
	return;
}

