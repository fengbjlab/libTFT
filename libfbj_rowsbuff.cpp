#include <deque>
#include "libfbj_rowsbuff.hpp"
#include "libfbj_base.hpp" // for exit_error() only

using namespace std;

struct RowsBuffer::RowsBufferData {
	bool			retro_only; // allow access future rows
	RowNum_t		stable_sze;	// stable size, eg, if require row -5 then stable_sze=6. 0 is unlimited.
	RowNum_t		current_rn;	// row number of the current row, 1st=0, nothing_read_yet=-1
	RowNum_t		rows_added; // number of rows added
	deque<RowData>	data;
	RowsBufferData():retro_only(true),stable_sze(1),current_rn(-1),rows_added(0) {}
};

RowsBuffer::~RowsBuffer()
{
	d->data.clear();
	delete d;
}

RowsBuffer::RowsBuffer():d(NULL)
{
	d=new RowsBufferData;
}

RowsBuffer::RowsBuffer(const RowsBuffer& othr):d(NULL) {
	d=new RowsBufferData;
	*d = *(othr.d);
}

RowsBuffer& RowsBuffer::operator=(const RowsBuffer& orig) {
	if(&orig != this)	*(this->d) = *(orig.d);
	return *this;
}

RowNum_t& RowsBuffer::rowID()			{ return d->current_rn;   }
RowNum_t  RowsBuffer::rows_added() const{ return d->rows_added;   }
void RowsBuffer::keep_all_input()		{ d->stable_sze = 0;      }
void RowsBuffer::allow_future_rows()	{ d->retro_only = false;  }
size_t	RowsBuffer::size() const		{ return d->data.size();  }
bool	RowsBuffer::empty() const		{ return d->data.empty(); }
bool	RowsBuffer::has_room() const	{ return d->stable_sze==0 || (RowNum_t)d->data.size()<d->stable_sze; }
RowsBuffer::iterator				RowsBuffer::begin()   const { return d->data.begin();  }
RowsBuffer::iterator				RowsBuffer::end()	  const { return d->data.end();    }
RowsBuffer::const_iterator			RowsBuffer::cbegin()  const { return d->data.cbegin(); }
RowsBuffer::const_iterator			RowsBuffer::cend()    const { return d->data.cend();   }
RowsBuffer::reverse_iterator		RowsBuffer::rbegin()  const { return d->data.rbegin(); }
RowsBuffer::reverse_iterator		RowsBuffer::rend()    const { return d->data.rend();   }
RowsBuffer::const_reverse_iterator	RowsBuffer::crbegin() const { return d->data.crbegin();}
RowsBuffer::const_reverse_iterator	RowsBuffer::crend()   const { return d->data.crend();  }

void RowsBuffer::set_oldest(RowNum_t row) // row should be <=0, 0 means the current row
{
	if (row>0) { exit_error("Previous row number should be <= 0"); }
	if (1-row > d->stable_sze) d->stable_sze = 1-row;
}

void RowsBuffer::clear()
{
	d->data.clear();
	d->current_rn=-1;
	d->rows_added=0;
}

void RowsBuffer::clear_ExceptTheLastRow()
{
	if (d->data.size()>1)
	{
		while (d->data.size()>1) d->data.pop_front();
		d->current_rn=0;
		d->rows_added=1;
	}
}

RowsBuffer::RowData& RowsBuffer::push_back()
{
	if (d->data.size()==1 && d->stable_sze==1)
		d->data.front().clear();
	else
	{
		d->data.push_back(RowData());
		if (d->stable_sze && (RowNum_t)d->data.size()>d->stable_sze) { d->data.pop_front(); }
	}
	d->current_rn = d->data.size()-1;
	++d->rows_added;
	return d->data[d->current_rn];
}

void RowsBuffer::push_back(RowsBuffer::RowData& input)
{
	if (d->data.size()==1 && d->stable_sze==1)
		d->data.front()=input;
	else
	{
		d->data.push_back(input);
		if (d->stable_sze && (RowNum_t)d->data.size()>d->stable_sze) { d->data.pop_front(); }
	}
	++d->rows_added;
	d->current_rn = d->data.size()-1;
}

void RowsBuffer::pop_back()
{
	d->data.pop_back();
	d->current_rn = (RowNum_t)(d->data.size())-1;
	--d->rows_added;
}

// Below are functions to access data.

string& RowsBuffer::operator[]( size_t col )
{
	if (empty()) exit_error("RowsBuffer has no data.");
	RowNum_t target = d->current_rn;
	if (target<0 || target>=(RowNum_t)d->data.size()) throw RowsBufferException();
	if (d->data[target].size()<=col)				  throw RowsBufferException();
	return d->data[target][col];
}

RowsBuffer::RowData& RowsBuffer::current_row()
{
	if (empty()) exit_error("RowsBuffer has no data.");
	RowNum_t target = d->current_rn;
	if (target<0 || target>=(RowNum_t)d->data.size()) throw RowsBufferException();
	return d->data[target];
}

string& RowsBuffer::operator() (RowNum_t row, size_t col)
{
	if (row>0 && d->retro_only) exit_error("Accessing future rows is not allowed.");
	if (empty()) exit_error("RowsBuffer has no data.");
	RowNum_t target = d->current_rn + row;
	if (target<0 || target>=(RowNum_t)d->data.size())	throw RowsBufferException();
	if (d->data[target].size()<=col)					throw RowsBufferException();
	return d->data[target][col];
}

RowsBuffer::RowData& RowsBuffer::target_row(RowNum_t row)
{
	if (row>0 && d->retro_only) exit_error("Accessing future rows is not allowed.");
	if (empty()) exit_error("RowsBuffer has no data.");
	RowNum_t target = d->current_rn + row;
	if (target<0 || target>=(RowNum_t)d->data.size()) throw RowsBufferException();
	return d->data[target];
}
