#include <iomanip>
#include <map>
#include <set>
#include <vector>
#include <chrono>
#include <thread>
#include <mutex>
#include <boost/interprocess/sync/file_lock.hpp>
#include "libfbj_elog.hpp"
#include "libfbj_program.hpp" // for license checking

using std::map;
using std::set;
using std::string;
using std::ostream;
using std::vector;

namespace  {
	struct ErrorType {
		bool to_show_count;		// whether to show number of errors
		bool to_show_message;	// whether to show the message
		unsigned count;			// number of errors occured
		set<string> assoc_str;	// associated strings for this error
		ErrorType():to_show_count(true),to_show_message(true),count(0){}
	};
	void write_elog(ostream& out, ErrorType& e, const string& message)
	{
		if (e.count && e.to_show_message)
		{
			if (e.to_show_count) out<<getpid()<<": "<<std::right<<std::setw(10)<<e.count<<" "<<message;
			else				 out<<getpid()<<":            "<<message;
			for (set<string>::iterator i=e.assoc_str.begin(); i!=e.assoc_str.end(); ++i)
			{
				out<<" ";
				if (i->empty()) out<<"[EmptyString]";
				else out<<*i;
			}
			out<<std::endl;
		}
	}
	/*void write_elog(ostream& out, ErrorType& e, const string& message)
	{
		if (e.count && e.to_show_message)
		{
			if (e.to_show_count) out<<"#"<<std::right<<std::setw(10)<<e.count<<" "<<message;
			else				 out<<"#           "<<message;
			for (set<string>::iterator i=e.assoc_str.begin(); i!=e.assoc_str.end(); ++i)
			{
				out<<" ";
				if (i->empty()) out<<"[EmptyString]";
				else out<<*i;
			}
			out<<std::endl;
		}
	}*/
}

struct ErrorLogger::ErrorLoggerData {
	static int			license_chk;// licensing of the class
	map<string,int>		tokens;		// location in vector<ErrorType>
	vector<ErrorType>	logs;		// log of all errors
	vector<string>		strs;		// same sorting as logs
	std::mutex			data_mutex;	// for logs
	ostream*			output;
	
	// interprocess mutex
	std::string						mutex_fn;
	boost::interprocess::file_lock	flock;
	
	ErrorLoggerData():logs(1),strs(1),output(&std::cerr){}
	void				lock();
	void				unlock();
};

int ErrorLogger::ErrorLoggerData::license_chk = program.nudge;

void ErrorLogger::ErrorLoggerData::lock()
{
	if (!mutex_fn.empty())
	{
		flock.lock();
	}
}

void ErrorLogger::ErrorLoggerData::unlock()
{
	if (!mutex_fn.empty())
	{
		flock.unlock();
	}	
}

ErrorLogger::~ErrorLogger() {
	if (d->output && !program.quiet) write(*(d->output));
	delete d;
}

ErrorLogger::ErrorLogger():d(NULL) {
	d=new ErrorLoggerData;
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

ErrorLogger::ErrorLogger(const ErrorLogger& othr):d(NULL) {
	std::cerr<<"copy constructor is prohibited for ErrorLogger !"<<std::endl;
	std::cout<<"copy constructor is prohibited for ErrorLogger !"<<std::endl;
	std::clog<<"copy constructor is prohibited for ErrorLogger !"<<std::endl;
	exit(1);
//	d=new ErrorLoggerData;
//	*d = *(othr.d);
}

ErrorLogger& ErrorLogger::operator=(const ErrorLogger& orig) {
	std::cerr<<"copy constructor is prohibited for ErrorLogger !"<<std::endl;
	std::cout<<"copy constructor is prohibited for ErrorLogger !"<<std::endl;
	std::clog<<"copy constructor is prohibited for ErrorLogger !"<<std::endl;
	exit(1);
//	if(&orig != this)	*(this->d) = *(orig.d);
//	return *this;
}

void ErrorLogger::set_lock(const std::string& fn) {
	d->mutex_fn=fn;
	if (!fn.empty()) d->flock=boost::interprocess::file_lock(fn.c_str());
}

bool ErrorLogger::lock_set()
{
	return !d->mutex_fn.empty();
}

void ErrorLogger::set_output(ostream& out) {
	d->output=&out;
} // you can do this: elog.set_output(cerr);

void ErrorLogger::set_output(ostream* out) {
	d->output=out; 
} // you can do this: elog.set_output(NULL);

int  ErrorLogger::get_token(const string& s, bool sh_c, bool sh_m) {
	if ((d->tokens.find(s)==d->tokens.end())) {
		d->tokens[s]=d->logs.size();
		ErrorType e;
		e.to_show_count=sh_c;
		e.to_show_message=sh_m;
		d->logs.push_back(e);
		d->strs.push_back(s);
	}
	return d->tokens[s];
}

// I previously did not check for t<0, trusting the program to throw "segmentation fault".
// But I found that it won't. It still runs and does nothing if t<0. So I add the check point.

void ErrorLogger::rewind(int t) {
	if (t<0) exit_error("Negative coordinate for elog.");
	d->logs[t].count=0;
}

void ErrorLogger::add(int t) {
	if (t<0) exit_error("Negative coordinate for elog.");
	if (!t) return;
	d->data_mutex.lock();
	++d->logs[t].count;
	d->data_mutex.unlock();
}

void ErrorLogger::add(int t, const string& a) {
	if (t<0) exit_error("Negative coordinate for elog.");
	if (!t) return;
	d->data_mutex.lock();
	++d->logs[t].count;
	d->logs[t].assoc_str.insert(a);
	d->data_mutex.unlock();
}

int  ErrorLogger::add(const string& s, bool sh) {
	int t=get_token(s,sh);
	d->data_mutex.lock();
	++d->logs[t].count;
	d->data_mutex.unlock();
	return t;
}

int  ErrorLogger::add(const string& s, const string& a, bool sh) {
	int t=get_token(s,sh);
	d->data_mutex.lock();
	++d->logs[t].count; 
	d->logs[t].assoc_str.insert(a);
	d->data_mutex.unlock();
	return t; 
}

void ErrorLogger::write(int token)
{
	if (d->output)
	{
		d->lock();
		write_elog(*(d->output),d->logs[token],d->strs[token]);
		d->unlock();
	}
}

void ErrorLogger::write()
{
	if (d->output)
	{
		d->lock();
		for (size_t i=1;i<d->logs.size();++i)
			write_elog(*(d->output),d->logs[i],d->strs[i]);
		d->unlock();
	}
}

void ErrorLogger::write(ostream& out, int token)
{
	d->lock();
	write_elog(out,d->logs[token],d->strs[token]);
	d->unlock();
}

void ErrorLogger::write(ostream& out)
{
	d->lock();
	for (size_t i=1;i<d->logs.size();++i)
		write_elog(out,d->logs[i],d->strs[i]);
	d->unlock();
}
