#include <cstring>								// for strcpy memcpy strlen
#include <set>
#include <mutex>								// g++47 needs it! clang++4.0/g++44 doesn't.
#include <thread>
#include <boost/regex.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/algorithm/string.hpp>			// for boost::split trim
#include <boost/algorithm/string/predicate.hpp> // for boost::algorithm::contains
#include <boost/format.hpp>						// for ftos ftos_FixWidth
#include "libfbj_base.hpp"
#include "libfbj_program.hpp"					// for lns

using namespace std;

//---------------- basic string related functions -----------------

void exit_error (const std::string& s)	{ lns << showe << s << fatal; }

//performance: string[] >> string::iterator >> string.find() >> boost::find_first(s1,s2) >> boost::algorithm::contains()

bool str_has(const std::string& s, const std::string& subs)
{
	return s.find(subs)!=std::string::npos;
}

bool str_startsw(const std::string& s1, const std::string& s2)
{
	int size1(s1.size()), size2(s2.size());
	if (size1<size2) return false;
	for (int i(0); i<size2; ++i) if (s1[i]!=s2[i]) return false;
	return true;
}

bool str_endsw(const std::string& s1, const std::string& s2)
{
	int size1(s1.size()), size2(s2.size()), i1(size1-size2);
	if (i1<0) return false;
	for (int i2(0); i2<size2; ++i1, ++i2) if (s1[i1]!=s2[i2]) return false;
	return true;
}

bool is_numeric(const std::string& s)
{
	try { boost::lexical_cast<double>(s);	return true; }
	catch (boost::bad_lexical_cast &)	{	return false; }
}

bool is_integer(const std::string& s)
{
	try { boost::lexical_cast<int>(s);	return true; }
	catch (boost::bad_lexical_cast &) {	return false; }
}

bool is_white_spaces(const std::string& str) // previously is_all_whiteSpaces
{
    for (std::string::const_iterator It = str.begin(); It != str.end() ; ++It)
    	if ( *It!=' ' && *It!='\t' ) return false;
	return true;
}

bool is_a_valid_name(const std::string& s)
{
	if (s.empty()) return false;
	string remains(s);
	extract_name(remains);
	return remains.empty();
}

double as_double_or_nan(const std::string& s)
{
	double v = std::numeric_limits<double>::signaling_NaN();
	read_val(s,v);
	return v;
}

int as_int_or_max(const std::string& s)
{
	int v = std::numeric_limits<int>::max();
	read_val(s,v);
	return v;
}

int as_int_or_min(const std::string& s)
{
	int v = -std::numeric_limits<int>::max();
	read_val(s,v);
	return v;
}

// ---------- information extraction from beginning -------------------

char extract_char(std::string& s)
{	if (s.empty()) throw bad_extracting();
	char c=s[0];
	s.erase(s.begin());
	return c;
}

double extract_double(std::string& s)
{
	if (s.empty()) throw bad_extracting();
	std::stringstream ss(s);
	double n;
	if (ss>>n) {
		if (ss.eof()) s.clear();		// this line is important, otherwise s unchanged
		else getline(ss,s,(char)EOF);	// previously \0 works fine, but EOF is better
	}
	else throw bad_extracting();
	return n;
}

int extract_int(std::string& s)
{
	if (s.empty()) throw bad_extracting();
	std::stringstream ss(s);
	int n;
	if (ss>>n) {
		if (ss.eof()) s.clear();		// this line is important, otherwise s unchanged
		else getline(ss,s,(char)EOF);	// previously \0 works fine, but EOF is better
	}
	else throw bad_extracting();
	return n;
}

std::string extract_alphabets(std::string& s)
{
	size_t size=s.size(), i;
	for (i=0;i<size;++i) if ( !isalpha(s[i]) ) break;
	std::string result=s.substr(0,i);
	s=s.substr(i);
	return result;
}

std::string extract_alnum(std::string& s)
{
	size_t size=s.size(), i;
	for (i=0;i<size;++i) if ( !isalnum(s[i]) ) break;
	std::string result=s.substr(0,i);
	s=s.substr(i);
	return result;
}

std::string extract_name(std::string& s)
{
	bool _al = false;	// character _ or alphabets occured
	size_t size=s.size(), i;
	for (i=0;i<size;++i)
	{
		if (isdigit(s[i])) { if (_al) continue; else break; }
		if ( isalpha(s[i]) || s[i]=='_' ) { _al=true; continue; }
		break;
	}
	std::string result=s.substr(0,i);
	s=s.substr(i);
	return result;
}

bool get_int (const vector<string>& INFO, const string& var, int& dest)
{
	for (auto& f:INFO)
		if (str_startsw(f,var+"=")) // (f.size()>var.size() && str_startsw(f,var) && f[var.size()]=='=')
		{
			try { dest=boost::lexical_cast<int>(f.substr(var.size()+1)); return true; }
			catch(...) { return false; }
		}
	return false;
}

double get_value (const vector<string>& INFO, const string& var)
{
	for (auto& f:INFO)
		if (str_startsw(f,var+"=")) // (f.size()>var.size() && str_startsw(f,var) && f[var.size()]=='=')
		{
			try { return boost::lexical_cast<double>(f.substr(var.size()+1)); }
			catch(...) { return std::numeric_limits<double>::signaling_NaN(); }
		}	
	return std::numeric_limits<double>::signaling_NaN();
}

std::string get_string (const vector<string>& INFO, const string& var)
{
	for (auto& f:INFO)
		if (str_startsw(f,var+"=")) // (f.size()>var.size() && str_startsw(f,var) && f[var.size()]=='=')
			return f.substr(var.size()+1);
	return string();
}

bool get_int_sw (const vector<string>& INFO, const string& var, int& dest)
{
	for (auto& f:INFO)
		if (str_startsw(f,var))
		{
			std::size_t found = f.find("=");
			if (found==std::string::npos) return false;
			try { dest=boost::lexical_cast<int>(f.substr(found+1)); return true; }
			catch(...) { return false; }
		}
	return false;
}

double get_value_sw (const vector<string>& INFO, const string& var)
{
	for (auto& f:INFO)
		if (str_startsw(f,var))
		{
			std::size_t found = f.find("=");
			if (found==std::string::npos) return std::numeric_limits<double>::signaling_NaN();
			try { return boost::lexical_cast<double>(f.substr(found+1)); }
			catch(...) { return std::numeric_limits<double>::signaling_NaN(); }
		}
	return std::numeric_limits<double>::signaling_NaN();
}

std::string get_string_sw (const vector<string>& INFO, const string& var)
{
	for (auto& f:INFO)
		if (str_startsw(f,var))
		{
			std::size_t found = f.find("=");
			if (found==std::string::npos) return string();
			return f.substr(found+1);
		}
	return string();
}

void put_string (vector<string>& INFO, const std::string& var, const std::string& val)
{
	for (auto& f:INFO)
		if (str_startsw(f,var+"=")) // (f.size()>var.size() && str_startsw(f,var) && f[var.size()]=='=')
		{
			f=(var+"="+val);
			return;
		}
	INFO.push_back(var+"="+val);
}

//---------------- string conversion functions --------------------

int int_width(long long num,int base)
{
	int w=0;
	if (num<0) { num=-num; ++w; }
	while (num)	{ num/=base; ++w; }
	return w;
}

int sig_digits(const std::string& input)
{
	int n=0; // number to the left of decimal point
	int d=0; // number to the right of decimal point
	std::string ns,ds;
	
	std::string ts = boost::to_lower_copy(input);
	boost::trim_left_if(ts,boost::is_any_of("0+- "));
	size_t found = ts.find('e');
	std::string ss; // significand string
	if (found!=std::string::npos)	ss=ts.substr(0, found); // remove exponential term
	else							ss=ts;					// no exponential term
	found = ss.find('.');
	if (found!=std::string::npos)
	{
		ns = ss.substr(0,found);
		ds = ss.substr(found+1);
		n=ns.size();
		if (n==0) boost::trim_left_if(ds,boost::is_any_of("0"));
		d=ds.size();
	}
	else
	{
		n=ss.size();
	}
	return n+d;
}

// c is [0-9a-z], num is [0-35], exit_error is out of range
char num2char(int num)
{
	if (num<0)  exit_error("in converting a number to character: N < 0");
	if (num>35) exit_error("in converting a number to character: N > 35");
	if (num<10) return (48+num);
	else		return (87+num); // not 97 here because num>=10
}

// c is [0-9a-z], num is [0-35], exit_error is out of range
int char2num(char c)
{
	if (c>='0'&&c<='9') return c-48;
	if (c>='a'&&c<='z') return c-87;
	exit_error("in converting a character to number: not [0-9a-z]");
	return -1;
}

int read_int(const std::string& str, std::string::size_type& loc)
{
	bool negative=false;
	bool read=false;
	if (str[loc]=='-') { negative=true; ++loc; }
	int num=0;
	for (;loc<str.size();++loc)
	{
		if (str[loc]<'0'||str[loc]>'9') break;
		num *= 10;
		num += char2num(str[loc]);
		read = true;
	}
	if (negative) num = -num;
	if (!read) exit_error("in reading an integer from a string: nothing read");
	return num;
}

// previously use long long, but cannot be compiled by gcc in cygwin, change to long.
// if really need long long, use the second function instead.
char * itoa( long long value, char* result, int base ) 
{
	// check that the base if valid
	if (base < 2 || base > 16)
	{	*result = 0; return result; 
	}
	
	char* out = result;
	long long quotient = value;
	
	do 
	{
		*out = "0123456789abcdef"[ std::abs( quotient % base ) ];
		++out;
		quotient /= base;
	} while ( quotient );
	
	// Only apply negative sign for base 10
	if ( value < 0 && base == 10) *out++ = '-';
	
	std::reverse( result, out );
	*out = 0;
	return result;
}


char * itoa_format(long long num,char align,int width,char fill,int base,char * output)
//long long is 8 byte, max=1e19
//num: number to be converted
//width: width of output std::string, not including \0
//fill: character to be filled in blank
//base: 2-36
{
	int i;
	char temstr[50];
	char target[50];
	long long temptr;
	long long tgtptr;
	tgtptr=0;
	temptr=0;
	bool negative;
	negative=false;
	if (num<0)
	{
		negative=true;
		num*=-1;
	}
	if (num==0)
	{
		temstr[temptr++]='0';
	}
	while (num>0)
	{
		temstr[temptr++]=num2char(num%base);
		num-=num%base;
		num/=base;
	}
	if (negative) temstr[temptr++]='-';

	if (temptr>=width)
	{	for (--temptr;temptr>=0;--temptr)
		{
			target[tgtptr++]=temstr[temptr];
		}
		target[tgtptr]='\0';
		strcpy(output,target);
		return output;
	}
	for (i=0;i<width;++i) target[i]=fill;
	target[i]='\0';
	switch (align)
	{
		case 'l':	tgtptr=temptr-1; break;
		case 'm':	tgtptr=width-1-(width-temptr)/2; break;
		case 'r':	tgtptr=width-1;break;
		default:	exit_error("wrong align number");
	}	
	for (i=0;i<temptr;++i)
		target[tgtptr--]=temstr[i];
	strcpy(output,target);
	return output;
}

std::string itos (long long num)
{
	char target[64];
	itoa(num,target,10);
	return target;
}

std::string itos_format (long long num,char align,int width,char fill,int base)
{
	char target[64];
	itoa_format(num,align,width,fill,base,target);
	return target;
}

std::string itos_max (long long num,long long maxnum)
{
	char target[64];
	itoa_format(num,'r',int_width(maxnum),'0',10,target);
	return target;
}

std::string itos_26based( long long quotient )
{
	bool negative=false;
	if (quotient==0) exit_error("itos_26based(v) cannot handle v=0");
	if (quotient<0) {
		negative=true;
		quotient = -quotient;
	}
	std::string out;
	do
	{	--quotient;
		out.insert(0, 1, "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[ quotient % 26 ] );
		quotient /= 26;
	} while ( quotient );
	if (negative) out.insert(0,1,'-');
	return out;
}

std::string ftos(double f, const std::string& MissingStr)
{
	if (std::isnan(f)) return MissingStr;
	char string_result[64];
	sprintf(string_result, "%g", f);
	return std::string(string_result);
}

std::string ftos(double f, int n, const std::string& MissingStr) // value_to_be_printed, num_sig_digits
{
	if (std::isnan(f)) return MissingStr;
	if (f==0) return "0";
	std::string s = boost::lexical_cast<string>(f); // str(boost::format("%g") % f);
	if (s=="nan")
	{
		return s;
	}
	else if (str_has(s,"e"))
	{
		string sd = substr_before_find( substr_before_find(s,"e"), ".");
		int d = (sd[0]=='-') ? sd.size()-1 : sd.size();
		s = str(boost::format ("%."+itos(std::max(n-d,0))+"e") % f);
		replace_all_R(s,"e-0","e-"); if (str_endsw(s,"e-")) s.resize(s.size()-2);
		replace_all_R(s,"e+0","e+"); if (str_endsw(s,"e+")) s.resize(s.size()-2);
		return s;
	}
	else
	{
		int d = sig_digits(substr_before_find(s,".")); if (f<1 && f>-1) d=0; // # sig dig before "."
		if (d>0 && d<=n)
		{
			return str(boost::format ("%."+itos(std::max(n-d,0))+"f") % f);
		}
		else if (d==0) // unlikely to happen
		{
			int d = sig_digits(s);					 // # sig dig after "."
			int z = substr_after_find(s,".").size(); // # dig     after "."
			if (d>=n)
				return str(boost::format ("%."+itos(z-(d-n))+"f") % f);
			else
			{
				s = str(boost::format ("%."+itos(std::max(n-d,0))+"e") % f);
				replace_all_R(s,"e-0","e-"); if (str_endsw(s,"e-")) s.resize(s.size()-2);
				replace_all_R(s,"e+0","e+"); if (str_endsw(s,"e+")) s.resize(s.size()-2);
				return s;
			}
		}
		else // d>n
		{
			s = str(boost::format ("%."+itos(n-1)+"e") % f);
			replace_all_R(s,"e-0","e-"); if (str_endsw(s,"e-")) s.resize(s.size()-2);
			replace_all_R(s,"e+0","e+"); if (str_endsw(s,"e+")) s.resize(s.size()-2);
			return s;
		}
	}
	return "ftos_err";
}

std::string ftos_FixWidth(double f, int n, const std::string& MissingStr) // value, width
{
	if (std::isnan(f)) return MissingStr;
	if (f==0) { return "0." + std::string(n-2, '0'); }
	string s1r,s2r;
	// s1
	{
		std::string s1 = str(boost::format("%e") % f);
		replace_all_R(s1,"e-0","e-"); if (str_endsw(s1,"e-")) s1.resize(s1.size()-2);
		replace_all_R(s1,"e+0","e+"); if (str_endsw(s1,"e+")) s1.resize(s1.size()-2);
		string sb = substr_before_find(s1,"e");
		string rs;
		if (str_has(sb,"."))
		{
			int acc = substr_after_find(sb,".").size() + (n-s1.size());
			rs = str(boost::format ("%."+itos(std::max(acc,0))+"e") % f);
		}
		else
		{
			int acc = n-s1.size()-1;
			rs = str(boost::format ("%."+itos(std::max(acc,0))+"e") % f);
		}
		replace_all_R(rs,"e-0","e-"); if (str_endsw(rs,"e-")) s1.resize(s1.size()-2);
		replace_all_R(rs,"e+0","e+"); if (str_endsw(rs,"e+")) s1.resize(s1.size()-2);
		if ((int)rs.size() >= n)
			s1r = rs;
		else
		{
			if (rs[0]=='-') s1r = "-" + string(n-rs.size(),'0') + rs.substr(1);
			else			s1r =       string(n-rs.size(),'0') + rs;
		}
	}
	// s2
	{
		std::string s2 = str(boost::format("%f") % f);
		if (str_has(s2,"."))
		{
			int acc = n - substr_before_find(s2,".").size() - 1;
			s2r = str(boost::format ("%."+itos(std::max(acc,0))+"f") % f);
		}
		else
		{
			int acc = n-s2.size()-1;
			s2r = str(boost::format ("%."+itos(std::max(acc,0))+"f") % f);
		}
	}
	int size1=s1r.size();
	int size2=s2r.size();
	if ((size1==n && size2==n) || (size1<n && size2<n))
	{
		int n1=sig_digits(s1r);
		int n2=sig_digits(s2r);
		if (n2>=n1) return s2r; else return s1r;
	}
	if (size1==n) return s1r;
	if (size2==n) return s2r;
	if (size1<n) return s1r;
	if (size2<n) return s2r;
	if (size1<size2) return s1r;
	if (size1>size2) return s2r;
	return "ftos_FixWidth_err";
}

std::string ftos_MaxWidth(double f, int n, const std::string& MissingStr) // value, width
{
	if (std::isnan(f)) return MissingStr;
	std::string s = str(boost::format("%g") % f);
	replace_all_R(s,"e-0","e-"); if (str_endsw(s,"e-")) s.resize(s.size()-2);
	replace_all_R(s,"e+0","e+"); if (str_endsw(s,"e+")) s.resize(s.size()-2);
	if ((int)s.size()<=n) return s;
	s = ftos_FixWidth(f,n);
	if (str_has(s,"e"))
	{
		string sb = substr_before_find(s,"e");
		if (str_has(sb,"."))
		{
			while (str_endsw(sb,"0")) sb.pop_back();
			if (str_endsw(sb,".")) sb.pop_back();
		}
		string sa = trim_before_find(s,"e");
		s = sb + sa;
	}
	else if (str_has(s,"."))
	{
		while (str_endsw(s,"0")) s.pop_back();
		if (str_endsw(s,".")) s.pop_back();
	}
	return s;
}

std::string ftos_MaxDigit(double f, int n, const std::string& MissingStr) // value, width
{
	if (std::isnan(f)) return MissingStr;
	std::string s = ftos(f,n);
	if (str_has(s,"e"))
	{
		string sb = substr_before_find(s,"e");
		if (str_has(sb,"."))
		{
			while (str_endsw(sb,"0")) sb.pop_back();
			if (str_endsw(sb,".")) sb.pop_back();
		}
		string sa = trim_before_find(s,"e");
		s = sb + sa;
	}
	else if (str_has(s,"."))
	{
		while (str_endsw(s,"0")) s.pop_back();
		if (str_endsw(s,".")) s.pop_back();
	}
	return s;
}

/*
std::string ftos_FixWidth(double f, int n, const std::string& MissingStr) // value, width
{
	if (std::isnan(f)) return MissingStr;
	if (f==0) { return "0." + std::string(n-2, '0'); }
	std::string s = str(boost::format("%g") % f);
	replace_all_R(s,"e-0","e-"); if (str_endsw(s,"e-")) s.resize(s.size()-2);
	replace_all_R(s,"e+0","e+"); if (str_endsw(s,"e+")) s.resize(s.size()-2);
	if (s=="nan")
	{
		return s;
	}
	else if (str_has(s,"e"))
	{
		string sb = substr_before_find(s,"e");
		string rs;
		if (str_has(sb,"."))
		{
			int acc = substr_after_find(sb,".").size() + (n-s.size());
			rs = str(boost::format ("%."+itos(std::max(acc,0))+"e") % f);
		}
		else
		{
			int acc = n-s.size()-1;
			rs = str(boost::format ("%."+itos(std::max(acc,0))+"e") % f);
		}
		replace_all_R(rs,"e-0","e-"); if (str_endsw(rs,"e-")) s.resize(s.size()-2);
		replace_all_R(rs,"e+0","e+"); if (str_endsw(rs,"e+")) s.resize(s.size()-2);
		if ((int)rs.size() >= n) return rs;
		if (rs[0]=='-') return "-" + string(n-rs.size(),'0') + rs.substr(1);
		else			return       string(n-rs.size(),'0') + rs;
	}
	else
	{
		if (str_has(s,"."))
		{
			int acc = substr_after_find(s,".").size() + (n-s.size());
			return str(boost::format ("%."+itos(std::max(acc,0))+"f") % f);
		}
		else
		{
			int acc = n-s.size()-1;
			return str(boost::format ("%."+itos(std::max(acc,0))+"f") % f);
		}
	}
	return "ftos_fixw_err";
}

std::string ftos_MaxWidth(double f, int n, const std::string& MissingStr) // value, width
{
	if (std::isnan(f)) return MissingStr;
	std::string s = str(boost::format("%g") % f);
	replace_all_R(s,"e-0","e-"); if (str_endsw(s,"e-")) s.resize(s.size()-2);
	replace_all_R(s,"e+0","e+"); if (str_endsw(s,"e+")) s.resize(s.size()-2);
	if ((int)s.size()<=n) return s;
	s = ftos_FixWidth(f,n);
	if (str_has(s,"e"))
	{
		string sb = substr_before_find(s,"e");
		if (str_has(sb,"."))
		{
			while (str_endsw(sb,"0")) sb.pop_back();
			if (str_endsw(sb,".")) sb.pop_back();
		}
		string sa = trim_before_find(s,"e");
		s = sb + sa;
	}
	else if (str_has(s,"."))
	{
		while (str_endsw(s,"0")) s.pop_back();
		if (str_endsw(s,".")) s.pop_back();
	}
	return s;
}*/

// this create a new string that will not contain any 'from' sub-string
// eg, replace // to /, if context=/// it returns /
// this is different from boost::algorithm::replace_all, which returns //
void replace_all_R(std::string& context, const std::string& from, const std::string& to)
{
	if (boost::algorithm::contains(to,from))
		exit_error("replace substring 'fr' to 'to'. 'fr' cannot be a substring of 'to'.");
    size_t foundHere;
	size_t fromsize=from.size();
    while((foundHere = context.find(from)) != std::string::npos)
    {
		context.replace(foundHere, fromsize, to);
    }
}

std::string replace_all_R_copy(const std::string& context, const std::string& from, const std::string& to)
{
	if (boost::algorithm::contains(to,from))
		exit_error("replace substring 'fr' to 'to'. 'fr' cannot be a substring of 'to'.");
	std::string result(context);
    size_t foundHere;
	size_t fromsize=from.size();
    while((foundHere = result.find(from)) != std::string::npos)
    {
		result.replace(foundHere, fromsize, to);
    }
    return result;
}

// similar to unix::echo except that it doesn't have \c \e \E and change \0num
std::string replace_escape_sequence_copy(const std::string& s)
{
	std::string o;
	for (size_t i=0; i<s.size(); ++i)
	{
		if (s[i]=='\\')
		{
			if (++i==s.size()) exit_error("parse escape sequence: nothing after \\");
			switch (s[i]) {
				case '\\':o.push_back('\\'); continue;
				case 'a': o.push_back('\a'); continue;
				case 'b': o.push_back('\b'); continue;
				case 'f': o.push_back('\f'); continue;
				case 'n': o.push_back('\n'); continue;
				case 'r': o.push_back('\r'); continue;
				case 't': o.push_back('\t'); continue;
				case 'v': o.push_back('\v'); continue;
				case 'x': // unlike unix::echo, it requires exactly 2 hex digits
				{
					char num=0; // no need for init but to suppress warnings.
					if (++i==s.size()) exit_error("parse escape sequence: there must be 2 hex digits after \\x");
					if		(s[i]>='0'&&s[i]<='9') num=s[i]-48;
					else if (s[i]>='a'&&s[i]<='f') num=s[i]-87;
					else if (s[i]>='A'&&s[i]<='F') num=s[i]-55;
					else exit_error("parse escape sequence: there must be 2 hex digits after \\x");
					num <<= 4;
					if (++i==s.size()) exit_error("parse escape sequence: there must be 2 hex digits after \\x");
					if		(s[i]>='0'&&s[i]<='9') num+=s[i]-48;
					else if (s[i]>='a'&&s[i]<='f') num+=s[i]-87;
					else if (s[i]>='A'&&s[i]<='F') num+=s[i]-55;
					else exit_error("parse escape sequence: there must be 2 hex digits after \\x");
					o.push_back(num); continue;
				}
				default:  o.push_back('\\'); break;
			}
		}
		o.push_back(s[i]);
	}
	return o;
}

#include <boost/algorithm/string/replace.hpp>	// for boost::algorithm::replace_all
std::string reverse_escape_sequence_copy(std::string s)
{
	boost::algorithm::replace_all(s ,"\\","\\\\"); // must be the first
	boost::algorithm::replace_all(s ,"\a","\\a");
	boost::algorithm::replace_all(s ,"\b","\\b");
	boost::algorithm::replace_all(s ,"\f","\\f");
	boost::algorithm::replace_all(s ,"\n","\\n");
	boost::algorithm::replace_all(s ,"\r","\\r");
	boost::algorithm::replace_all(s ,"\t","\\t");
	boost::algorithm::replace_all(s ,"\v","\\v");
	return s;
}

std::string reverse_escape_char(const char c)
{
	if (c=='\\') return "\\\\";
	if (c==0x07) return "\\a";
	if (c==0x08) return "\\b";
	if (c==0x0c) return "\\f";
	if (c==0x0a) return "\\n";
	if (c==0x0d) return "\\r";
	if (c==0x09) return "\\t";
	if (c==0x0b) return "\\v";
	return s(c);
}

std::string to_cmd(const std::string& input)
{
	std::string output = reverse_escape_sequence_copy(input);
	boost::algorithm::replace_all(output,"\"","\\\""); // no need to use \ for a single quote within double quotes
	if (output!=input || output.empty() ||
		output.find(' ')!=std::string::npos ||
		output.find('\t')!=std::string::npos ||
		output.find('\n')!=std::string::npos ||
		output.find('\\')!=std::string::npos ||
		output.find('#')!=std::string::npos ||
		output.find('<')!=std::string::npos ||
		output.find('>')!=std::string::npos ||
		output.find('(')!=std::string::npos ||
		output.find(')')!=std::string::npos ||
		output.find('|')!=std::string::npos ||
		output.find(';')!=std::string::npos ||
		output.find('&')!=std::string::npos ||
		output[0]=='!') output="\""+output+"\"";
	return output;
}

//---------------- string trimming functions ----------------

// Bellow are a set of functions to trim a string.
// Input is not a reference, because otherwise it makes char[] not possible,
// The input will not be changed.

std::string substr_after_find(std::string str, const std::string& sch, std::string::size_type index)
{
	std::string::size_type loc(str.find(sch,index));
	if (loc!=std::string::npos)
		str=str.substr(loc+sch.size());
	return str;
}

std::string substr_before_find(std::string str, const std::string& sch, std::string::size_type index)
{
	std::string::size_type loc(str.find(sch,index));
	if (loc!=std::string::npos)
		str.resize(loc);
	return str;
}

std::string trim_before_find(std::string str, const std::string& sch, std::string::size_type index)
{
	std::string::size_type loc(str.find(sch,index));
	if (loc!=std::string::npos)
		str=str.substr(loc);
	return str;
}

std::string trim_after_find(std::string str, const std::string& sch, std::string::size_type index)
{
	std::string::size_type loc(str.find(sch,index));
	if (loc!=std::string::npos)
		str.resize(loc+sch.size());
	return str;
}

std::string substr_after_rfind(std::string str, const std::string& sch)
{
	std::string::size_type loc(str.rfind(sch));
	if (loc!=std::string::npos)
		str=str.substr(loc+sch.size());
	return str;
}

std::string substr_before_rfind(std::string str, const std::string& sch)
{
	std::string::size_type loc(str.rfind(sch));
	if (loc!=std::string::npos)
		str.resize(loc);
	return str;
}

std::string trim_before_rfind(std::string str, const std::string& sch)
{
	std::string::size_type loc(str.rfind(sch));
	if (loc!=std::string::npos)
		str=str.substr(loc);
	return str;
}

std::string trim_after_rfind(std::string str, const std::string& sch)
{
	std::string::size_type loc(str.rfind(sch));
	if (loc!=std::string::npos)
		str.resize(loc+sch.size());
	return str;
}

//---------------- string printing functions ----------------

void _print_text(std::vector<std::string>& paragraph, int& hanging, std::string& prefix, int max_length, std::string& result)
{
	result+=prefix;
	std::string newprf;
	newprf.assign(prefix.size()+hanging,' ');
	int cur_length=prefix.size();
	bool first_word = true;
	
	for (int j=0; j<(int)paragraph.size(); ++j) // j must be int because negative is possible
	{
		std::string& s=paragraph[j];
		if (str_startsw(s,"__PRTXTFMT_INDENT__"))
		{
			try
			{
				int n=boost::lexical_cast<int>(s.substr(19));
				prefix.assign(n,' ');
				if (int(prefix.size())+hanging<0) hanging = -prefix.size();
				newprf.assign(prefix.size()+hanging,' ');
			}
			catch (boost::bad_lexical_cast &) { exit_error("__PRTXTFMT_INDENT__ parameter should be an integer."); }
			continue;
		}
		if (str_startsw(s,"__PRTXTFMT_HANGING__"))
		{
			try
			{
				int n=boost::lexical_cast<int>(s.substr(20));
				hanging=n;
				if (int(prefix.size())+hanging<0) hanging = -prefix.size();
				newprf.assign(prefix.size()+hanging,' ');
			}
			catch (boost::bad_lexical_cast &) { exit_error("__PRTXTFMT_HANGING__ parameter should be an integer."); }
			continue;
		}
		if (first_word)
		{
			if ( max_length-cur_length >= (int)s.size() )
			{
				result += s;
				cur_length += s.size();
				first_word = false;
			}
			else
			{
				result += s.substr(0,max_length-cur_length);
				s=s.substr(max_length-cur_length);  --j;
				result.push_back('\n');
				result+=newprf;
				cur_length=newprf.size();
				first_word = true;
			}
		}
		else
		{
			if ( max_length-cur_length >= (int)s.size()+1 )
			{
				result.push_back(' ');
				result += s;
				cur_length += s.size()+1;
				first_word = false;
			}
			else
			{
				--j;
				result.push_back('\n');
				result+=newprf;
				cur_length=newprf.size();
				first_word = true;
			}
		}
	}
}

void _print_test_read_priority(std::vector<std::string>& paragraph, int& priority)
{
	if (paragraph.size())
	{
		if (str_startsw(paragraph[0],"__PRTXTFMT_PRIORITY__"))
		{
			try	{	priority=boost::lexical_cast<int>(paragraph[0].substr(21));	}
			catch (boost::bad_lexical_cast &) { exit_error("__PRTXTFMT_PRIORITY__ parameter should be an integer."); }
			paragraph.erase(paragraph.begin());
		}
	}
}

std::string print_text(std::string& text, const int hanging, const std::string prefix, const int max_length, const int pri_thr)
{
	int priority = 0 ; // 0 is top priority
	int			newhng = hanging;
	std::string newprf = prefix;
	std::string result;
	std::vector<std::string> paragraph;
	std::string word;
	bool begin_of_word=true;
	if (text[text.size()-1]!='\n') text.push_back('\n'); // with this, no need to "if (text[text.size()-1]!='\n') .." after for()
	for (size_t i=0; i<text.size() ; ++i)
	{
		if ( text[i]==' ' || text[i]=='\t' )
		{
			if (begin_of_word) word.push_back(text[i]);
			else { paragraph.push_back(word); word.clear(); begin_of_word=true; } // word can't be empty
		}
		else if (text[i]=='\n')
		{
			if (!begin_of_word) { paragraph.push_back(word); word.clear(); }
			_print_test_read_priority(paragraph,priority);
			if (priority <= pri_thr)
			{
				_print_text(paragraph,newhng,newprf,max_length,result);
				result.push_back('\n');
			}
			paragraph.clear();
			word.clear();
			begin_of_word=true;
		}
		else
		{
			word.push_back(text[i]);
			begin_of_word=false;
		}
	}
	return result;
}

std::string print_a_line(const std::string& pattern, int line_length) {
	if (line_length<=0) return "";
	std::string out_str;
	int i;
	for (i=pattern.length(); i<=line_length; i+=pattern.length()) out_str+=pattern;
	i-=pattern.length();
	if (i!=line_length) out_str+=pattern.substr(0,line_length-i);
	return out_str;
}

std::string print_a_line(const std::pair<std::string,int>& pattern, int line_length) {
	if (line_length<=0) return "";
	std::string out_str;
	int i;
	for (i=pattern.second; i<=line_length; i+=pattern.second) out_str+=pattern.first;
	i-=pattern.second;
	if (i!=line_length) out_str+=pattern.first.substr(0,line_length-i);
	return out_str;
}

std::string to_term(std::string msg)
{
	int lines,columns;
	get_term_dim(lines,columns);
	boost::algorithm::replace_all(msg,"__ROUND_BIG__",print_a_line(line_pattern::round_big,columns));
	boost::algorithm::replace_all(msg,"__ROUND_SMALL__",print_a_line(line_pattern::round_small,columns));
	boost::algorithm::replace_all(msg,"__SQUARED__",print_a_line(line_pattern::squared,columns));
	boost::algorithm::replace_all(msg,"__DIAMOND__",print_a_line(line_pattern::diamond,columns));
	boost::algorithm::replace_all(msg,"__ZIGZAG_BIG__",print_a_line(line_pattern::zigzag_big,columns));
	boost::algorithm::replace_all(msg,"__ZIGZAG_ZIP__",print_a_line(line_pattern::zigzag_zip,columns));
	boost::algorithm::replace_all(msg,"__ZIGZAG_SMALL__",print_a_line(line_pattern::zigzag_small,columns));
	boost::algorithm::replace_all(msg,"__WAVE__",print_a_line(line_pattern::wave,columns));
	boost::algorithm::replace_all(msg,"__DASH__",print_a_line(line_pattern::dash,columns));
	boost::algorithm::replace_all(msg,"__LONG_DASH__",print_a_line(line_pattern::long_dash,columns));
	boost::algorithm::replace_all(msg,"__DASH_DOT__",print_a_line(line_pattern::dash_dot,columns));
	boost::algorithm::replace_all(msg,"__DOT__",print_a_line(line_pattern::dot,columns));
	boost::algorithm::replace_all(msg,"__DOT_SPARSE__",print_a_line(line_pattern::dot_sparse,columns));
	boost::algorithm::replace_all(msg,"__DOT_DOUBLE__",print_a_line(line_pattern::dot_double,columns));
	boost::algorithm::replace_all(msg,"__SINGLE_LINE__",print_a_line(line_pattern::single_line,columns));
	boost::algorithm::replace_all(msg,"__DOUBLE_LINE__",print_a_line(line_pattern::double_line,columns));
	boost::algorithm::replace_all(msg,"__THICK_LINE__",print_a_line(line_pattern::thick_line,columns));
	boost::algorithm::replace_all(msg,"__ASTERISK__",print_a_line(line_pattern::asterisk,columns));
	boost::algorithm::replace_all(msg,"__CHAIN_LIKE__",print_a_line(line_pattern::chain_like,columns));
	boost::algorithm::replace_all(msg,"__DOUBLE_HELIX__",
								  print_a_line("_   _.-'```''-..",columns)+"\n"+
								  print_a_line(" '.'_.-|`|`':-. ",columns)+"\n"+
								  print_a_line(".  '.  | |  | |'",columns)+"\n"+
								  print_a_line("/'.  '.| |  | |_",columns)+"\n"+
								  print_a_line("_/ '.  `.!__!-` ",columns)+"\n"+
								  print_a_line("     `'-.___.-'`",columns));
	std::string outs = print_text(msg,0,"",columns);
	return outs;
}

//------------ system environment variables and commands --------------
/*
EnvironmentVariable bool_env(const std::string& name) {
	EnvironmentVariable result;
	if (name.empty()) return result;
	string env=exec("echo $"+name,true);
	if (env.empty()) return result;
	if (str_endsw(env,"\n")) env.pop_back();
	if (env.empty()) return result;
	result.exist=true;
	result.value=env;
	if (IsBool(env))
		result.is_true = IsYes(result.value);
	else
		exit_error("Unknown value \""+env+"\" for $"+name+". Candidates are 0/1/t/f/y/n/true/false/yes/no/on/off.");
	return result;
}

EnvironmentVariable str_env(const std::string& name) {
	EnvironmentVariable result;
	if (name.empty()) return result;
	string env=exec("echo $"+name,true);
	if (env.empty()) return result;
	if (str_endsw(env,"\n")) env.pop_back();
	if (env.empty()) return result;
	result.exist=true;
	result.value=env;
	return result;
}
*/

EnvironmentVariable bool_env(const std::string& name) {
	EnvironmentVariable result;
	if (name.empty()) return result;
	char* env=getenv(name.c_str());
	if (env==NULL) return result;
	result.exist=true;
	result.value=env;
	if (IsBool(result.value))
		result.is_true = IsYes(result.value);
	else
		exit_error("Unknown value \""+result.value+"\" for $"+name+". Candidates are 0/1/t/f/y/n/true/false/yes/no/on/off.");
	return result;
}

EnvironmentVariable str_env(const std::string& name) {
	EnvironmentVariable result;
	if (name.empty()) return result;
	char* env=getenv(name.c_str());
	if (env==NULL) return result;
	result.exist=true;
	result.value=env;
	return result;
}

void get_term_dim(int& lines, int& columns)
{
	// First try "stty size".
	// If run within a pipe, stty cannot see the terminal and give an error "stty: standard input: Invalid argument".
	string term_dim;
	try { term_dim = exec("stty size 2>/dev/null",false); }
	catch (const std::exception& error) { } // it will fail if run within a pipe, so do not do anything.
	if (!term_dim.empty())
	{
		try {
			lines	= extract_int(term_dim); extract_char(term_dim);
			columns	= extract_int(term_dim)-1;
			return;
		}
		catch (...) { }
	} //*/

	/*/ Then try "tput".
	// Strange that it can see the terminal even within a pipe in a command line, but not within a pipe inside a program.
	// plus not all machines have it. So it's value is quite limited. I keep it here anyway.
	try
	{
		string str_lines = exec("tput lines 2>/dev/null");  lines = extract_int(str_lines);
		string str_columns = exec("tput cols 2>/dev/null"); columns = extract_int(str_columns);
		return ;
	}
	catch (...) { } //*/

	/*/ Then try environment variable. But it doesn't work within a pipe inside a program.
	try {
		string str_lines = str_env("LINES").value;		lines = extract_int(str_lines);
		string str_columns = str_env("COLUMNS").value;	columns = extract_int(str_columns);
		return ;
	}
	catch (...) { } //*/

	// if everything failed
	lines = 25;
	columns = DefTermWidth;
}

std::string HOME()
{
	static std::string home_path;
	static std::mutex mt;
	if (home_path.empty())
	{
		mt.lock();
		if (home_path.empty())
		{
			home_path = str_env("HOME").value;
			if (!str_endsw(home_path,"/"))  home_path.push_back('/');
		}
		mt.unlock();
	}
	return home_path;
}

std::string linux_command_which_dir(const std::string& prg, const std::string& var, bool check_exe) // previously: char const * const prg
{
	std::string path_var = str_env(var).value;
	vector<string>	path_str;
	boost::split(path_str, path_var, boost::is_any_of(":"), boost::token_compress_on);
	for (each_element(path_str,it))
	{
		if (it->empty()) continue; // should not happen due to boost::token_compress_on
		if (!str_endsw(*it,"/")) it->push_back('/');
		if (check_exe)	{ if (access((*it+prg).c_str(), X_OK) == 0)	return *it; }
		else			{ if (FileExists(*it+prg))					return *it; }
	}
	return "";
}

//---------- file system functions -------------

bool FileExists(std::string strFilename)
{
	boost::filesystem::path p(strFilename);
	if (!boost::filesystem::exists(p)) return false;
	if (boost::filesystem::is_directory(p)) return false;
	return true;
}

bool DirExists(std::string strFilename)
{
	boost::filesystem::path p(strFilename);
	if (!boost::filesystem::exists(p)) return false;
	if (!boost::filesystem::is_directory(p)) return false;
	return true;
}

void mkdir(const std::string& DirName)
{
	boost::filesystem::path p(DirName);
	if (!boost::filesystem::exists(p))
	{
		boost::system::error_code ec;
		boost::filesystem::create_directories(p,ec);
		if (ec) exit_error("cannot make directory "+DirName);
	}
	if (!boost::filesystem::is_directory(p))
	{
		exit_error(DirName+" already exist but is a file");
	}
}

// to create a folder:
// method 1: if (mkdir(config.c_str(), 0700)!=0) exit_error("Permission denied."); // 0700 works, 0600/0660/0666 can't write files inside
// method 2: boost::filesystem::path p(config); if (!boost::filesystem::create_directory(p))   exit_error("Failed"); // not recursive
// method 3: boost::filesystem::path p(config); if (!boost::filesystem::create_directories(p)) exit_error("Failed"); //     recursive

void filename_change_home_path(std::string& filename)
{
	if		(str_startsw(filename,"~/"))	 { filename=filename.substr(2); filename.insert(0,HOME()); }
	else if (str_startsw(filename,"$HOME/")) { filename=filename.substr(6); filename.insert(0,HOME()); }
	while (str_has(filename,"$"))
	{
		string l = substr_before_find(filename,"$");
		string r = substr_after_find(filename,"$");
		string n = extract_name(r);
		string t = str_env(n).value;
		filename = l + t + r;
	}
}

bool find_file_zippedOrNot(std::string filename,std::string& trueName)
{
	trueName.clear();
	filename_change_home_path(filename);
	if (FileExists(filename))		{ trueName=filename;		return true; }
	if (FileExists(filename+".gz"))	{ trueName=filename+".gz";  return true; }
	if (FileExists(filename+".bz2")){ trueName=filename+".bz2"; return true; }
	if (str_endsw(filename,".gz"))	{ trueName=filename.substr(0,filename.size()-3); return FileExists(trueName); }
	if (str_endsw(filename,".bz2"))	{ trueName=filename.substr(0,filename.size()-4); return FileExists(trueName); }
	return false;
}

bool find_file_zippedOrNot(const std::string& filename)
{
	std::string trueName;
	return find_file_zippedOrNot(filename,trueName);
}

bool is_stdin(const std::string& filename) { return (filename.empty()||filename=="cin"||filename=="stdin"); }
std::string label_stdin() { return "stdin"; }

// the OS macro is from http://sourceforge.net/p/predef/wiki/OperatingSystems/ which is the gold standard.
// the method is from http://stackoverflow.com/questions/1023306/finding-current-executables-path-without-proc-self-exe
#if defined(__linux__)
std::string self_path() {
	char path[10240];
	int len = readlink("/proc/self/exe", path, sizeof(path)-1);			// previously "size_t len ="  but return type is ssize_t (signed)
	if (len == -1) exit_error("Failed to find the executable file.");
	path[len] = '\0';
	std::string result(path);
	return result;
}
#elif defined(__APPLE__) || defined(__MACH__)
#include <cstdlib>			// for realpath
#include <mach-o/dyld.h>	// for _NSGetExecutablePath
std::string self_path() {
	char path[10240];
	uint32_t size = sizeof(path);
	if ( _NSGetExecutablePath(path, &size) != 0 ) {
		// if size not enough, _NSGetExecutablePath rewrite size (includes the \0 at the end) and return !=0.
		char* newp = (char*)malloc(size);
		if ( newp==NULL ) exit_error("Failed to find the executable file.");
		if (_NSGetExecutablePath(newp, &size)!=0) exit_error("Failed to find the executable file.");
		char* real = realpath(newp, NULL);
		if ( real==NULL ) exit_error("Failed to find the executable file.");
		std::string result(real);
		free(real);
		free(newp);
		return result;
	}
	else {
		char* real = realpath(path, NULL);
		if ( real==NULL ) exit_error("Failed to find the executable file.");
		std::string result(real);
		free(real);
		return result;
	}
}
#elif defined(__FreeBSD__) || defined(__NetBSD__) || defined(__OpenBSD__) || defined(__bsdi__) || defined(__DragonFly__)
std::string self_path() {
	char path[10240];
	int len = readlink("/proc/curproc/file", path, sizeof(path)-1);	// previously "size_t len ="  but return type is ssize_t (signed)
	if (len == -1) exit_error("Failed to find the executable file.");
	path[len] = '\0';
	std::string result(path);
	return result;
}
#elif defined(sun) || defined(__sun)
std::string self_path() {
	char path[10240];
	int len = readlink("/proc/self/path/a.out", path, sizeof(path)-1);	// previously "size_t len ="  but return type is ssize_t (signed)
	if (len == -1) exit_error("Failed to find the executable file.");
	path[len] = '\0';
	std::string result(path);
	return result;
}
#else // MSDOS __MSDOS__ _MSDOS __DOS__ _WIN16 _WIN32 _WIN64 OS2 _OS2 __OS2__ __TOS_OS2__ __ANDROID__ __CYGWIN__
#error this operating system is not supported
#endif

// http://stackoverflow.com/questions/478898/how-to-execute-a-command-and-get-output-of-command-within-c
// http://stackoverflow.com/questions/7807755/reading-popen-results-in-c
// http://stackoverflow.com/questions/15058876/how-to-capture-the-exit-code-and-stderr-of-the-command-that-is-run-in-c
// Need to consider whether cmd should contain 2>/dev/null. Recommend that libraries have it but programs don't.
#include <cstdio>
struct exec_exception : public std::exception { virtual char const* what() const throw() { return "exec() error"; } };
std::string exec(const std::string& cmd, bool nothrow) {
	FILE* pipe = popen(cmd.c_str(), "r");
	if (!pipe) return "Pipe error for "+cmd;
	char buffer[128];
	std::string result = "";
	while(!feof(pipe)) {
		if(fgets(buffer, 128, pipe) != NULL)
			result += buffer;
	}
	int exit_code = pclose(pipe); // 0=success 1=fail. Some say -1=fail too. So I use WEXITSTATUS(exit_code)!=0 as failure.
	if (WEXITSTATUS(exit_code)!=0 && !nothrow) throw exec_exception();
	return result;
}

// ----- boost::iostreams file open/close -----

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>

bool openInpFile(boost::iostreams::filtering_istream& f, std::string n)
{
	filename_change_home_path(n);
	if (is_stdin(n))
		f.push(std::cin);
	else
	{
		std::string newname;
		if (!find_file_zippedOrNot(n,newname)) return false;
		if		(str_endsw(newname,".gz"))	f.push(boost::iostreams::gzip_decompressor());
		else if (str_endsw(newname,".bz2"))	f.push(boost::iostreams::bzip2_decompressor());
		else if (str_endsw(newname,".zip"))	exit_error(".zip file not supported.");
		else if (str_endsw(newname,".rar"))	exit_error(".rar file not supported.");
		else if (str_endsw(newname,".Z"))	exit_error(".Z file not supported.");
		else if (str_endsw(newname,".7z"))	exit_error(".7z file not supported.");
		f.push(boost::iostreams::file_source(newname));
		if (!(f.component<boost::iostreams::file_source>(f.size()-1)->is_open())) return false;
	}
	return true;
}

// I have tried std::ios::app, but the result is weird: the file can be gunzip/bunzip2 in command line and is fine;
// but when I use boost::iostreams to read it again, the appended content is not accessible for .gz, while number of fields
// is wrong for .bz2. In this trial I used boost 1.46, zlib 1.2.5, bzip2 1.0.6, gcc 4.2
bool openOutFile(boost::iostreams::filtering_ostream& f, std::string n)
{
	filename_change_home_path(n);
	if (str_has(n,"/")) mkdir(substr_before_rfind(n,"/"));
	if		(str_endsw(n,".gz"))	f.push(boost::iostreams::gzip_compressor());
	else if (str_endsw(n,".bz2"))	f.push(boost::iostreams::bzip2_compressor());
	else if (str_endsw(n,".zip"))	exit_error(".zip file not supported.");
	else if (str_endsw(n,".rar"))	exit_error(".rar file not supported.");
	else if (str_endsw(n,".Z"))		exit_error(".Z file not supported.");
	else if (str_endsw(n,".7z"))	exit_error(".7z file not supported.");
	if (n.empty() || n=="cout" || n=="stdout") f.push(std::cout);
	else if (n=="cerr")			f.push(std::cerr);
	else if (n=="clog")			f.push(std::clog);
	else					{	f.push(boost::iostreams::file_sink(n, std::ios::binary));
		if (!(f.component<boost::iostreams::file_sink>(f.size()-1)->is_open())) return false; }
	return true;
}

// ---- other read/write utilities ----

// return whether the string is found.
// If found, file loc=end of string; otherwise, file.eof()=true;
int read_until(std::istream& file,const std::string& s)
{
	size_t i,l=s.length();
	while (!file.eof())
	{
		for (i=0; i<l ; ++i)
		{
			char c=file.get();
			if (file.eof() || c!=s[i]) break;
			//  file.eof() is required, otherwise read_until(inpf,s(EOF)) will always return true!
		}
		if (i==l) return 1;
	}
	return 0;
}

// write text to file named filename, lock the file before writing so that other processes cannot open it
/* To lock the file for reading/writing at the system level, so that no other process can access the file.
 Neither fstream nor FILE lock the file even in write/append mode.
 flockfile() is not suitable because it's at the application level (other threads cannot access the file locked by one thread).
 boost::interprocess::file_lock is only advisory lock (applications can choose to ignore it, as it does by MacOSX TextEdit).
 default to try 1 day, sleep 1 sec after each failure.
 the following method is not thread safe, valgrind --tool=drd gives some errors */
#include <boost/thread/thread.hpp> // for boost::this_thread::sleep
#include <boost/interprocess/sync/file_lock.hpp> // for boost::interprocess::file_lock
#include <boost/date_time/posix_time/posix_time_types.hpp>

// to compile this function, -lboost_thread -lrt on Linux and -lboost_thread on Mac. Recommend to use the next function.
void write_file_try_lock(const std::string& text, const std::string& filename, const std::ios_base::openmode mode, int num_tries, int sleep_sec, int incr)
{
	if (!boost::filesystem::exists(boost::filesystem::path(filename))) { FILE* file=fopen(filename.c_str(),"w"); fclose(file); }
	boost::interprocess::file_lock flock(filename.c_str());
	for ( int i=0, to_sleep=sleep_sec; i<num_tries; ++i, boost::this_thread::sleep(boost::posix_time::seconds(to_sleep+=incr)) )
	{
		if (flock.try_lock())
		{
			std::fstream file;
			file.open(filename.c_str(), mode);
			if (!file.is_open()) exit_error("Cannot open "+filename+" even after locking.");
			file << text;
			file.close();
			flock.unlock();
			return;
		}
	}
	exit_error("Cannot open "+filename+" after many trials.");
}

void write_file_lock(const std::string& text, const std::string& filename, const std::ios_base::openmode mode)
{
	if (!boost::filesystem::exists(boost::filesystem::path(filename))) { FILE* file=fopen(filename.c_str(),"w"); fclose(file); }
	boost::interprocess::file_lock flock(filename.c_str());
	flock.lock();
	std::fstream file;
	file.open(filename.c_str(), mode);
	if (!file.is_open()) exit_error("Cannot open "+filename+" even after locking.");
	file << text;
	file.close();
	flock.unlock();
}

// --------- date and time functions ----------------

// NTP client uaing Boost::Asio.
// 2010-05-11
// Jansson Consulting
// http://www.p-jansson.com
// Dedicated to the Public Domain.
// link with the boost_date_time and boost_system libraries.

std::string macro_date_to_boost(const std::string& in) { // "Nov 20 2011" => 2011-11-20
//	std::string in(__DATE__);
	std::map<std::string,std::string> month_abr2num;
	month_abr2num["jan"]="01";
	month_abr2num["feb"]="02";
	month_abr2num["mar"]="03";
	month_abr2num["apr"]="04";
	month_abr2num["may"]="05";
	month_abr2num["jun"]="06";
	month_abr2num["jul"]="07";
	month_abr2num["aug"]="08";
	month_abr2num["sep"]="09";
	month_abr2num["oct"]="10";
	month_abr2num["nov"]="11";
	month_abr2num["dec"]="12";
	return in.substr(7)+"-"+month_abr2num[boost::to_lower_copy(in.substr(0,3))]+"-"+boost::replace_first_copy(in.substr(4, 2)," ","0");
}

#include <boost/asio.hpp>

class NtpClient
{
public:
	static boost::posix_time::ptime GetTime();
	static boost::posix_time::ptime GetTime( const char* ntpServer );
};

boost::posix_time::ptime NtpClient::GetTime()
{
	return GetTime("pool.ntp.org");
}

boost::posix_time::ptime NtpClient::GetTime( const char* ntpServer )
{
	using boost::asio::ip::udp;
	boost::asio::io_service io_service;
	
	udp::resolver resolver(io_service);
	udp::resolver::query query(udp::v4(), ntpServer, "ntp");
	udp::endpoint receiver_endpoint = *resolver.resolve(query);
	
	udp::endpoint sender_endpoint;
	
	boost::uint8_t data[48] = {
		0x1B,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
		0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
	};
	
	udp::socket socket(io_service);
	socket.open(udp::v4());
	
	socket.send_to(
				   boost::asio::buffer(data),
				   receiver_endpoint);
	socket.receive_from(
						boost::asio::buffer(data),
						sender_endpoint);
	
	typedef boost::uint32_t u32;
	const u32 iPart(
					static_cast<u32>(data[40]) << 24
					| static_cast<u32>(data[41]) << 16
					| static_cast<u32>(data[42]) << 8
					| static_cast<u32>(data[43])
					);
	const u32 fPart(
					static_cast<u32>(data[44]) << 24
					| static_cast<u32>(data[45]) << 16
					| static_cast<u32>(data[46]) << 8
					| static_cast<u32>(data[47])
					);
	
	using namespace boost::posix_time;
	const ptime pt( boost::gregorian::date(1900,1,1),
				   milliseconds( iPart*1.0E3 + fPart*1.0E3/0x100000000ULL ) );
	return pt;
}

bool date_incorrect()
{
	static bool result = true;
	static bool untest = true;
	static std::mutex mt;
	if (untest)
	{
		mt.lock();
		if (untest)
		{
			boost::gregorian::date network_date (NtpClient::GetTime().date());
			boost::gregorian::date_duration dd = network_date - today_s_date();
			result = (dd.days()>1 || dd.days()<0);
			untest = false;
		}
		mt.unlock();
	}
	return result;
	// cannot test (network_date!=today_s_date) because network_date is zone 0, today_s_date is local,
	// they may be different. But network_date should always be 1 greater or equal to today_s_date.
}

boost::posix_time::ptime			duration_start(){ return boost::posix_time::microsec_clock::local_time(); }
boost::posix_time::time_duration	duration_end(const boost::posix_time::ptime& time_begin)
{
	boost::posix_time::ptime time_end  = boost::posix_time::microsec_clock::local_time();
	boost::posix_time::time_duration d = time_end - time_begin;
	return d;
}

// =================== file processing ===================

bool is_EndOfLine(char c) {
	return (c=='\n' || c=='\r' || c==EOF);
}

bool is_blank_row(std::istream & file) {
	return is_EndOfLine(file.peek());
}

bool skip_EndOfLine(std::istream & file) { // EOL not extracted yet, but sould be right at it
	char c=file.peek();
	if (c=='\n') {	file.get(); if (file.peek()=='\r') file.get(); return true; }
	if (c=='\r') {	file.get(); if (file.peek()=='\n') file.get(); return true; }
	return false;
}

bool skip_one_line(std::istream & file) {
	for (char c=file.get(); c!=EOF; c=file.get())
	{
		if (c=='\n') {	if (file.peek()=='\r') file.get(); return true; }
		if (c=='\r') {	if (file.peek()=='\n') file.get(); return true; }
	}
	return false; // happens at the EOF
}

bool skip_spaces(std::istream & file, int n) {
	for (int i=0;i<n;++i)
	{
		if (file.peek()!=' ') return false;
		file.get();
	}
	return true;
}

bool skip_whitespaces(std::istream & file) {
	char c; while ((c=file.peek())==' ' || c=='\t') file.get();
	return true;
}

bool skip_blank_lines(std::istream & file) {
	char c; while ((c=file.peek())=='\n' || c=='\r') file.get();
	return true;
}

bool skip_comments(std::istream & file) {
	while (file.peek()=='#') skip_one_line(file);
	return true;
}

bool next_row(std::istream & file) {
	for (char c=file.peek();;c=file.peek())
	{
		if (c=='#') { skip_one_line(file); continue; }
		// previously skip_whitespaces, not ideal if I want to allow for empty fields at the beginning!
		// if (skip_whitespaces(file) && is_blank_row(file)) return false;
		if (is_blank_row(file)) return false;
		return true;
	}
	return true; // should not happen
}

// EndOfLine is not extracted, return whether field exist. field starts from 0.
bool get_field(std::istream & file, int field, std::string& os, const std::string& dlt, bool skip_eol)
{
	os.clear();
	int i=0;
	while (!is_EndOfLine(file.peek()))
	{
		char c=file.get();
		if (dlt.find(c)!=std::string::npos) { ++i; continue; }
		if (i==field) os.push_back(c);
	}
	if (skip_eol) skip_EndOfLine(file);
	return i>=field;
}

// EndOfLine is not extracted, return whether read something. field starts from 0.
bool get_fields(std::istream & file, std::vector<std::string>& ov, const std::string& dlt, bool skip_eol)
{
	ov.assign(1, std::string());
	std::vector<std::string>::reverse_iterator it = ov.rbegin();
	while (!is_EndOfLine(file.peek()))
	{
		char c=file.get();
		if (dlt.find(c)!=std::string::npos) { ov.push_back(std::string()); it=ov.rbegin(); continue; }
		it->push_back(c);
	}
	if (skip_eol) skip_EndOfLine(file);
	return (ov.size()>1 || !ov[0].empty());
}


// =================== OS Independent file reading functions ====================

/////////////////////////////////////////////////////////////////////////////////
// <iostream> and <string> provides the following functions:
//
// *** when a delimiter is found, it is extracted and discarded ***
// #include <iostream>
// istream& istream::getline( char* buffer, streamsize num );
// istream& istream::getline( char* buffer, streamsize num, char delim );
// #include <string>
// istream& std::getline( istream& is, string& s );
// istream& std::getline( istream& is, string& s, Char delimiter );
//
// *** the delimiter is not extracted ***
// int      istream::get();
// istream& istream::get( char& ch );
// istream& istream::get( char* buffer, streamsize num );
// istream& istream::get( char* buffer, streamsize num, char delim );
// istream& istream::get( streambuf& buffer );
// istream& istream::get( streambuf& buffer, char delim );
//
// These are subject to non-portability due to different EOL characters.
// The following functions are safe for various EOL, with the same naming scheme:
// getline() extract EOL but not get(); copy_line() does but copy_row() does not.
/////////////////////////////////////////////////////////////////////////////////

// http://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf
std::istream& safeGetline(std::istream& is, std::string& t)
{
	t.clear();
	
	// The characters in the stream are read one-by-one using a std::streambuf.
	// That is faster than reading them one-by-one using the std::istream.
	// Code that uses streambuf this way must be guarded by a sentry object.
	// The sentry object performs various tasks,
	// such as thread synchronization and updating the stream state.
	
	std::istream::sentry se(is, true);
	std::streambuf* sb = is.rdbuf();
	
	for(;;) {
		int c = sb->sbumpc();
		switch (c) {
			case '\n':
				return is;
			case '\r':
				if(sb->sgetc() == '\n')
					sb->sbumpc();
				return is;
			case EOF:
				// Also handle the case when the last line has no line ending
				if(t.empty())
					is.setstate(std::ios::eofbit);
				return is;
			default:
				t += (char)c;
		}
	}
}

std::istream& safeGetline_add(std::istream& is, std::string& t)
{
	std::istream::sentry se(is, true);
	std::streambuf* sb = is.rdbuf();
	
	for(;;) {
		int c = sb->sbumpc();
		switch (c) {
			case '\n':
				return is;
			case '\r':
				if(sb->sgetc() == '\n')
					sb->sbumpc();
				return is;
			case EOF:
				// Also handle the case when the last line has no line ending
				if(t.empty())
					is.setstate(std::ios::eofbit);
				return is;
			default:
				t += (char)c;
		}
	}
}

// extract one char, return char by value
int osi_getchar(std::istream & file)
{
	int c=file.get();
	if (c=='\n') {	if (file.peek()=='\r') file.get(); }
	if (c=='\r') {	if (file.peek()=='\n') file.get(); }
	return c;
}

// extracts EOL, modify string
void osi_getline(std::istream & file,std::string & os)
{
	os.clear();
	//	while (!is_EndOfLine(file.peek())) os.push_back(file.get());
	//	skip_EndOfLine(file);
	for (char c=file.get(); c!=EOF; c=file.get())
	{
		if (c=='\n') {	if (file.peek()=='\r') file.get(); break; }
		if (c=='\r') {	if (file.peek()=='\n') file.get(); break; }
		os.push_back(c);
	}
}

// EndOfLine is not extracted, modify string
void osi_get(std::istream & file,std::string & os)
{
	os.clear();
	while (!is_EndOfLine(file.peek())) os.push_back(file.get());
}

// extracts EOL, return string by value
std::string osi_getline(std::istream & file)
{
	std::string os;
	osi_getline(file,os);
	return os;
}

// EndOfLine is not extracted, return string by value
std::string osi_get(std::istream & file)
{
	std::string os;
	osi_get(file,os);
	return os;
}

// extracts EOL, modify string, using es as EOL character
void osi_getline(std::istream & file,std::string& os,char es)
{	char c;
	os.clear();
	while ((c=file.peek())!=es && c!=EOF)
	{
		file.get();
		os.push_back(c);
	}
	if (c==es) file.get();
}

// EndOfLine is not extracted, modify string, using es as EOL character
void osi_get(std::istream & file,std::string& os,char es)
{	char c;
	os.clear();
	while ((c=file.peek())!=es && c!=EOF)
	{
		file.get();
		os.push_back(c);
	}
}

// copy one line from one file to another, EndOfLine is not extracted
void osi_copy_row(std::istream & ifle, std::ostream & ofle)
{
	char c;
	while (!is_EndOfLine(c=ifle.peek()))
	{
		ifle.get();
		ofle<<c;
	}
}

// copy one line from one file to another, extracts EOL
void osi_copy_line(std::istream & ifle, std::ostream & ofle)
{
	char c;
	while (!is_EndOfLine(c=ifle.get())) ofle<<c;
	if (c=='\n') {	if (ifle.peek()=='\r') ifle.get(); }
	if (c=='\r') {	if (ifle.peek()=='\n') ifle.get(); }
	ofle<<endl;
}

// extracts EOL, modify string, using endc as EOL character
void osi_getline_allowing_brackets_or_quotes(std::istream & file, std::string& outs, char endc)
{
	outs.clear();
	int double_quoted=0;
	int single_quoted=0;
	int round_bracket=0;
	int squre_bracket=0;
	int curly_bracket=0;
	int angle_bracket=0;
	for (char c=file.get(); c!=EOF; outs.push_back(c), c=file.get()) {
		if (c=='\'') single_quoted ^= 1;
		if (c=='\"') double_quoted ^= 1;
		if (single_quoted || double_quoted) continue;
		switch (c) {
			case '(': ++round_bracket; continue;
			case '[': ++squre_bracket; continue;
			case '{': ++curly_bracket; continue;
			case '<': ++angle_bracket; continue;
			case ')': if (round_bracket) { --round_bracket; continue; } break;
			case ']': if (squre_bracket) { --squre_bracket; continue; } break;
			case '}': if (curly_bracket) { --curly_bracket; continue; } break;
			case '>': if (angle_bracket) { --angle_bracket; continue; } break;
			default: break;
		}
		if (round_bracket || squre_bracket || curly_bracket || angle_bracket) continue;
		if (c==endc) break;
	}
}

// read one line of file and tokenize it in UNIX-style, EOL is extracted, robust to diff EOL / no EOL
// 1) "\c" => \c, except when c = \`$" which results in c only.
// 2)  \c  =>  c, even when c = t n SPACE
// 3) '\c' => \c
// 4) if last char is \ (no preseeding \) then go to next line (Unix do this also to `")
void tokenize(std::istream& file, std::vector<std::string>& result)
{
	int single_q=0;		// whether is single quoted
	int double_q=0;		// whether is double quoted
	int sth_read=0;		// something's been read
	std::string word;
	for (char c=osi_getchar(file); !file.eof(); c=osi_getchar(file))
	{
		if		(c=='\n' || c=='\r')
		{
			if (single_q || double_q) exit_error("tokenizer, line ends before closing a pair of quotation marks");
			else {
				if (sth_read) { result.push_back(word); word.clear(); } // word could be empty!
				sth_read=0;
			}
			break;
		}
		else if (c==' ' || c=='\t') {
			if (single_q || double_q) word.push_back(c);
			else {
				if (sth_read) { result.push_back(word); word.clear(); } // word could be empty!
				sth_read=0;
			}
		}
		else if (c=='\'') {
			if (double_q) word.push_back(c);
			else if (single_q) single_q=0;
			else single_q=1;
			++sth_read;
		}
		else if (c=='\"') {
			if (single_q) word.push_back(c);
			else if (double_q) double_q=0;
			else double_q=1;
			++sth_read;
		}
		else if (c=='\\') {
			if		(single_q)
				word.push_back(c);
			else if (double_q) {
				char c2=osi_getchar(file);
				if (file.eof()) exit_error("tokenizer, there must be a character after \\");
				if (c2!='\\' && c2!='`' && c2!='$' && c2!='\"') word.push_back('\\');
				word.push_back(c2);
			}
			else {
				char c2=osi_getchar(file);
				if (file.eof()) exit_error("tokenizer, there must be a character after \\");
				if (c2!='\n' && c2!='\r') word.push_back(c2);
			}
			++sth_read;
		}
		else {
			word.push_back(c);
			++sth_read;
		}
	}
	if (sth_read) { result.push_back(word); word.clear(); }
	if (single_q || double_q) exit_error("tokenizer, line ends before closing a pair of quotation marks");
}

// -------------------- serie --------------------

/* old format of serie: no options. Deprecated.
 int fill_number(std::string& str,int val,int width,char fill,char starting_char)
 {
 std::string search_str(s(starting_char));
 for (int n=0;;++n)
 {
 std::size_t loc=str.find(search_str);
 if ( loc == std::string::npos ) return n;
 str.replace(loc,1,itos_format(val,'r',width,fill,10));
 }
 exit_error("It won't happen.");
 return -1;
 }*/

// a serie's format: @#[/option1/option2/..]//  (# is int & >0)

// loc should point to xxx/.., return xxx, loc point to /
std::string _serie_read_opt(const std::string& str, size_t& loc, char starting_char)
{
	size_t bkp=loc;
	std::string opt;
	for (;loc<str.size() && str[loc]!='/';++loc)
		if (str[loc]==starting_char)		{ loc=bkp; return ""; }
		else opt.push_back(str[loc]);
	if (loc==str.size() || str[loc]!='/')	{ loc=bkp; return ""; }
	return opt;
}

int serie_fill_string(std::string& str,int seq, int val, std::map<int,std::vector<std::string> >& rplmap,char starting_char)
{
	std::string search_str(s(starting_char) + itos(seq));
	for (size_t loc=0,n=0;;++loc)
	{
		size_t bkp=str.find(search_str,loc);
		if ( bkp == std::string::npos ) return n;
		loc = bkp + search_str.size();
		if (loc==str.size() || str[loc]!='/')	{ loc=bkp; continue; }
		int width=0;
		char fill=' ';
		char align='l'; // left / middle / right
		char item='a';
		int item_num=-1;
		for (;;)
		{
			if (++loc==str.size())				{ loc=bkp; break; }
			if (str[loc]=='/')					{ break; }
			std::string opt=_serie_read_opt(str,loc,starting_char);
			if (opt.empty())					{ loc=bkp; break; }
			// make use of the opt bellow
			switch (opt[0]) {
				case 'w':
					try { width=boost::lexical_cast<unsigned>(opt.substr(1)); }
					catch (boost::bad_lexical_cast &) {	exit_error("parsing serie option wN/ : N is not an unsigned"); }
					break;
				case 'p':
					if (opt.size()>1)	fill=opt[1];
					else exit_error("parsing serie option pC/ : C not specified");
					break;
				case 'a':
					if (opt.size()>1)	align=tolower(opt[1]);
					else exit_error("parsing serie option aC/ : C not specified");
					if (align!='l' && align!='m' && align!='r') exit_error("parsing serie option aA/ : A must be l/m/r");
					break;
				case 'i':
					if (opt.size()>1)	item=tolower(opt[1]);
					else exit_error("parsing serie option iC/ : C not specified");
					if 		(item=='a') item_num=-1;
					else if (item>='1'&&item<='9') item_num=boost::lexical_cast<int>(s(item))-1;
					else exit_error("parsing serie option iC/ : C must be a or 1-9.");
					break;
				default:
					exit_error("parsing serie option: unknown option "+s(opt[0]));
					break;
			}
		}
		if (loc==bkp) continue;
		std::string repl;
		if (item=='a')
		{
			repl=rplmap[seq][val];
		}
		else
		{
			std::string whole=rplmap[seq][val];
			vector<string>	items;
			boost::split(items, whole, boost::is_any_of("\t"));
			if (item_num>=(int)items.size()) exit_error("not sufficient number of items for serie_fill_string");
			repl=items[item_num];
		}
		if ((int)repl.size()<width)
		{
			switch (align) {
				case 'l':repl.insert(repl.end(),  width-repl.size(),fill);	break;
				case 'm':repl=print_in_middle(repl,s(fill),width);			//break;
				case 'r':repl.insert(repl.begin(),width-repl.size(),fill);	break;
				default: exit_error("parsing serie option: unknown alignment "+s(align));
			}
		}
		str.replace(bkp,loc-bkp+1,repl);
		loc = bkp + repl.size() -1;
		++n;
	}
	exit_error("It won't happen.");
	return -1;
}

int serie_fill_number(std::string& str,int seq, double val,int width,char fill,char starting_char)
{
	std::string search_str(s(starting_char) + itos(seq));
	for (size_t loc=0,n=0;;++loc)
	{
		size_t bkp=str.find(search_str,loc);
		if ( bkp == std::string::npos ) return n;
		loc = bkp + search_str.size();
		if (loc==str.size() || str[loc]!='/')	{ loc=bkp; continue; }
		char align='r'; // left / middle / right
		for (;;)
		{
			if (++loc==str.size())				{ loc=bkp; break; }
			if (str[loc]=='/')					{ break; }
			std::string opt=_serie_read_opt(str,loc,starting_char);
			if (opt.empty())					{ loc=bkp; break; }
			// make use of the opt bellow
			switch (opt[0]) {
				case 'w':
					try { width=boost::lexical_cast<unsigned>(opt.substr(1)); }
					catch (boost::bad_lexical_cast &) {	exit_error("parsing serie option wN/ : N is not an unsigned"); }
					break;
				case 'p':
					if (opt.size()>1)	fill=opt[1];
					else exit_error("parsing serie option pC/ : C not specified");
					break;
				case 'a':
					if (opt.size()>1)	align=tolower(opt[1]);
					else exit_error("parsing serie option pC/ : C not specified");
					if (align!='l' && align!='m' && align!='r') exit_error("parsing serie option aA/ : A must be l/m/r");
					break;
				default:
					exit_error("parsing serie option: unknown option "+s(opt[0]));
					break;
			}
		}
		if (loc==bkp) continue;
		std::string repl;
		if (width<0 || fill=='d')	repl=ftos(val);
		else						repl=itos_format((int)val,align,width,fill,10);
		str.replace(bkp,loc-bkp+1,repl);
		loc = bkp + repl.size() -1;
		++n;
	}
	exit_error("It won't happen.");
	return -1;
}

// loc must be valid, test whether loc is at @ which indicates a serie
// return the serie's ID, which must >0. If 0, not a serie.
// after return, the loc is at the last /, or unchanged if not a serie
int _serie_number(const std::string& str, size_t& loc, char starting_char)
{
	// not starting with @, not a serie
	if (str[loc]!=starting_char) return 0;
	size_t bkp = loc++;
	// nothing after @, not a serie
	if (loc==str.size())			{ loc=bkp; return 0; }
	// if no number next to @, not a serie
	if (str[loc]<'0'||str[loc]>'9')	{ loc=bkp; return 0; }
	// has number, read it
	int num = read_int(str,loc);
	// next to num must be /
	if (loc==str.size() || str[loc]!='/') { loc=bkp; return 0; }
	// next to /, read options xxx/, or end sign /
	for (;;)
	{
		if (++loc==str.size())		{ loc=bkp; return 0; }
		if (str[loc]=='/')			{ break; }
		if (_serie_read_opt(str,loc,starting_char).empty())
		{ loc=bkp; return 0; }
	}
	return num;
}

int number_of_series(const std::string& str,char starting_char)
{
	std::set<int> numset;
	for (size_t loc=0; loc<str.size(); ++loc)
	{
		int n = _serie_number(str,loc,starting_char);
		if (n) numset.insert(n);
	}
	return numset.size();
}

int max_serie_ID(const std::string& str,char starting_char)
{
	int maxn=0;
	for (size_t loc=0; loc<str.size(); ++loc)
	{
		int n = _serie_number(str,loc,starting_char);
		if (n>maxn) maxn=n;
	}
	return maxn;
}

// ------------- functions for data input -----------------

bool IsNo(const std::string& answer)
{
	std::string s=to_lower_copy(answer);
	if (s=="0")			return true;
	if (s=="f")			return true;
	if (s=="n")			return true;
	if (s=="no")		return true;
	if (s=="false")		return true;
	return false;
}

bool IsYes(const std::string& answer)
{
	std::string s=to_lower_copy(answer);
	if (s=="1")			return true;
	if (s=="t")			return true;
	if (s=="y")			return true;
	if (s=="yes")		return true;
	if (s=="true")		return true;
	if (s=="0")			return false;
	if (s=="f")			return false;
	if (s=="n")			return false;
	if (s=="no")		return false;
	if (s=="false")		return false;
	throw input_exception();
}

bool IsBool(const std::string& answer)
{
	try { IsYes(answer); return true; }
	catch (input_exception &) { return false; }
}

std::string str_YesOrNo(const bool Yes)
{
	if (Yes)	return "Yes";
	else		return "No";
}

void clear_cin()
{
	if (std::cin.fail()) { std::cin.clear(); std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); }
	else if (std::cin.rdbuf()->sungetc()!=EOF) std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

// std::numeric_limits<std::streamsize>::max() is previously INT _ MAX
std::string input_line(const std::string& s)
{
	clear_cin();
	if (s.size()) std::cout<<s;
	std::string instr;
	std::getline(std::cin,instr);
	return instr;
}

// in Unix systems, use ^d to input EOF, but make sure the last line is ended with EOL,
// otherwise the text remains on the screen and not yet go to cin.
std::string input_paragraph(const std::string& s)
{
	clear_cin();
	if (s.size()) std::cout<<s;
	std::string instr;
	std::getline(std::cin,instr,(char)EOF);
	return instr;
}

bool input_yes(const std::string& s)
{
	std::cout << s + " (y/n): ";
	for (;;)
	{
		try
		{
			std::string instr = input_line();
			return IsYes(instr);
		}
		catch (input_exception &)
		{
			std::cout<<"Please input 1/0/yes/no/y/n/true/false/t/f/on/off (case insensitive): ";
			continue;
		}
	}
}

// here s cannot be default to ="" because of the previous function.
// otherwise this_func("text") will be interpreted as this_func(true)
bool input_yes(const bool dflt, const std::string& s)
{
	std::string dflt_s = dflt ? "y" : "n";
	std::cout << s + " (y/n) ["+dflt_s+"]: ";
	for (;;)
	{
		try
		{
			std::string instr = input_line();
			if (instr.empty()) return dflt;
			else return IsYes(instr);
		}
		catch (input_exception &)
		{
			std::cout<<"Please input 1/0/yes/no/y/n/true/false/t/f/on/off (case insensitive) or press ENTER for "<<dflt_s<<": ";
			continue;
		}
	}
}

char input_char(const std::string& s)
{
	std::cout << s << ": ";
	for (;;)
	{
		std::string instr=replace_escape_sequence_copy(input_line());
		if (instr.size()==1) return instr[0];
		else
		{
			std::cout<<"Please input exactly 1 character: ";
			continue;
		}
	}
}

char input_char(const char dflt, const std::string& s)
{
	std::cout << s << " ["<<dflt<<"]: ";
	for (;;)
	{
		std::string instr=replace_escape_sequence_copy(input_line());
		if		(instr.size()==1)	return instr[0];
		else if (instr.empty())		return dflt;
		else
		{
			std::cout<<"Please input at most 1 character (default="<<dflt<<"): ";
			continue;
		}
	}
}

// Two ways to get data from cin:   1) cin>>data  2) getline(cin,phrase); But they have different properties:
// the former doesn't extract the final EOL but the latter does; and getline can take spaces or tabs.
// So if these two methods are used together (eg. cin>>data then getline()) it will cause a problem.
// To solve this problem, each time before getline(), getchar() or get(), clear the buffer with the following code:
//    if (std::cin.fail()) { std::cin.clear(); std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); }
//    else if (std::cin.rdbuf()->sungetc()!=EOF) std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
// Explanation: If sungetc() is not EOF, it means last input is not completely extracted, which could be due to
// a) cin>>xxx read all inputs and leave an EOL in buffer, or b) cin>>xxx read only parts of the inputs.
// if cin>>xxx did not read anything, it's failbit is set and can be detected by fail(). All these circumstances
// need to be treated with std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
// This method should be robust to both Linux/Mac.
// There's a problem: if user inputs nothing but press enter for cin>>xxx, it never ends. So never use cin>>xxx.

/* debug
 int main (int argc, char * const argv[])
 {
 int i;
 cin>>i;
 string s=input_line("please input a line:");
 cout << "your input integer : "<< i << '\n';
 cout << "your input line is : \"" << s << "\"" << '\n';
 return 0;
 }//*/


struct progress_timer::progress_timer_data {
	boost::posix_time::ptime time_begin;
	boost::posix_time::ptime time_last;
	boost::posix_time::ptime time_now;
	boost::posix_time::time_duration td_now_from_begin;
	boost::posix_time::time_duration td_now_from_last;
	boost::posix_time::time_duration td_end_from_now;
	unsigned long	_total;
	unsigned long	_count;
	double			_factor;
	string			_prefix;
	progress_timer_data():_total(0),_count(0),_factor(0),_prefix("# Total runtime ") {}
};

progress_timer::~progress_timer(){ delete d; }
progress_timer::progress_timer() { d=new progress_timer_data; }
progress_timer::progress_timer(const progress_timer& othr) { d=new progress_timer_data; *d=*(othr.d);}
progress_timer& progress_timer::operator=(const progress_timer& orig) { if(&orig!=this) *(this->d)=*(orig.d); return *this; }

void progress_timer::prefix(const string& s)
{
	d->_prefix=s;
}

void progress_timer::start(unsigned long N) {
	if (N==0) exit_error("progress_time total cannot be 0");
	d->_total=N;
	d->_factor=100.0/d->_total;
	d->time_begin = d->time_last = boost::posix_time::microsec_clock::local_time();
}

void progress_timer::finish() {
	cerr<<"                                                                                \r";
	cerr<<d->_prefix<<d->td_now_from_begin<<endl;
}

unsigned long progress_timer::operator+=( unsigned long increment )
{
	d->_count += increment;
	d->time_now=boost::posix_time::microsec_clock::local_time();
	d->td_now_from_last = d->time_now - d->time_last;
	if (d->td_now_from_last.total_seconds() && d->_total!=d->_count)
	{
		double pct = d->_count * d->_factor;
		double mul = (d->_total-d->_count) / (double)d->_count;
		d->td_now_from_begin = d->time_now - d->time_begin;
		long seconds_now_from_begin=d->td_now_from_begin.total_seconds();
		d->td_end_from_now = boost::posix_time::seconds(seconds_now_from_begin * mul);
		cerr<<"                                                                                \r";
		cerr<<' '<<pct<<"% done; remaining "<<d->td_end_from_now<<'\r';
		d->time_last= d->time_now;
	}
	return d->_count;
}
