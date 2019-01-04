#include "libfbj_genepi.hpp"
#include "libfbj_math.hpp"
#include "libfbj_file.hpp"
#include "libfbj_base.hpp"
#include "libfbj_program.hpp"
#include <cmath>
#include <boost/assign/list_of.hpp>					// for map_list_of
#include <boost/lexical_cast.hpp>					// lexical_cast
#include <boost/algorithm/string.hpp>				// for is_any_of
#include <boost/range/algorithm/replace_if.hpp>		// for replace_if
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <algorithm>
#include <mutex> // g++47 needs it! clang++4.0/g++44 doesn't.
#include <thread>

namespace genepi {

	using namespace std;

	// data for genome
	string MainPath = "./"; // previously "$BING/work/data/h g 1 9/"; now has to call set_path(), controled by path_set.
	string NpNmFile = "_NCBI_NP_NM";
	string SynSymTr = "_NCBI_Synonyms_Symbol";
	string RsqSymTr = "_NCBI_refSeq_Symbol";
	string txID_gID = "_NCBI_refSeq_GeneID";
	string read_thr = "_NCBI_readthroughGeneID";
	string pseudoid = "_NCBI_pseudoGeneID";
	string gID_MFOT = "_NCBI_MFOT";
	string uGID_can = "_appris_refGene";
	string SEQ_PATH = "seq/";
	int chrM_num=-1;
	int chrX_num=-1;
	int chrY_num=-1;
	int chrXYnum=-1;
	int chrA_min=-1;
	int chrA_max=-1;
	int PAR1_x_b=-1;
	int PAR1_x_e=-1;
	int PAR1_y_b=-1;
	int PAR1_y_e=-1;
	int PAR2_x_b=-1;
	int PAR2_x_e=-1;
	int PAR2_y_b=-1;
	int PAR2_y_e=-1;
	int MHC_chr=0;
	int MHC_beg=0;
	int MHC_end=0;
	bool path_set=false;
	map<int,string>		chr_num2name_;
	map<string,int>		chr_name2num_;
	map<string,string>	chr_2StdName_;
	string GenomeName;
	
	// data for genes
	string	GDB_name = "refGeneLite";
	set<string>	incl_g;	// gene symbol filter
	set<string>	incl_t;	// transcript ID filter
	bool	canonic = false; // prv true, which is bad
	bool	filt_tx = false; // filter by _NCBI_MFOT _appris_refGene _NCBI_readthroughGeneID _NCBI_pseudoGeneID

	void set_path(string& s, string genome)
	{
		if (s.empty()) exit_error("genepi::set_path(s) argument s cannot be empty.");
		if (!str_endsw(s,"/")) s.push_back('/');
		if (genome.empty()) genome=boost::to_lower_copy(s);
		else				genome=boost::to_lower_copy(genome);
		if (genome=="./")
		{
			boost::filesystem::path full_path(boost::filesystem::current_path());
			genome=boost::to_lower_copy(full_path.string());
		}
		int found_genome=0;
		if (str_has(genome,"hg19")||str_has(genome,"grch37"))
		{
			++found_genome;
			chrM_num=23;
			chrX_num=24;
			chrY_num=25;
			chrXYnum=26;
			chrA_min=1;
			chrA_max=22;
			PAR1_x_b=60001;
			PAR1_x_e=2699520;
			PAR1_y_b=10001;
			PAR1_y_e=2649520;
			PAR2_x_b=154931044;
			PAR2_x_e=155260560;
			PAR2_y_b=59034050;
			PAR2_y_e=59363566;
			MHC_chr=6;
			MHC_beg=25702009;
			MHC_end=33378772;
			GenomeName="GRCh37";
		}
		if (str_has(genome,"hg38")||str_has(genome,"grch38"))
		{
			++found_genome;
			chrM_num=23;
			chrX_num=24;
			chrY_num=25;
			chrXYnum=26;
			chrA_min=1;
			chrA_max=22;
			PAR1_x_b=10001;
			PAR1_x_e=2781479;
			PAR1_y_b=10001;
			PAR1_y_e=2781479;
			PAR2_x_b=155701383;
			PAR2_x_e=156030895;
			PAR2_y_b=56887903;
			PAR2_y_e=57217415;
			MHC_chr=6;
			MHC_beg=25701784;
			MHC_end=33410995;
			GenomeName="GRCh38";
		}
		if (found_genome==0) exit_error("cannot infer genome from input ("+s+","+genome+") or the genome is not supported.");
		if (found_genome>1) exit_error("inferred multiple known genomes from input ("+s+","+genome+")");
		MainPath = s;
		NpNmFile = MainPath+NpNmFile;
		SynSymTr = MainPath+SynSymTr;
		RsqSymTr = MainPath+RsqSymTr;
		txID_gID = MainPath+txID_gID;
		read_thr = MainPath+read_thr;
		pseudoid = MainPath+pseudoid;
		gID_MFOT = MainPath+gID_MFOT;
		uGID_can = MainPath+uGID_can;
		SEQ_PATH = MainPath+SEQ_PATH;
		if (!str_endsw(SEQ_PATH,"/")) SEQ_PATH.push_back('/');
		for (Rows_in_File(in,SEQ_PATH+"chr_std_name",2)) chr_2StdName_[in[0]]=in[1];
		for (Rows_in_File(in,SEQ_PATH+"chr_number",2))
		{
			int chrNo=0;
			if (read_val_ge(in[1],chrNo,1))
			{
				chr_name2num_[in[0]]=chrNo;
				chr_num2name_[chrNo]=in[0];
			}
			else exit_error("failed to read chromosome number from "+in[1]+" in "+SEQ_PATH+"chr_number");
		}
		if (chr_num2name_.rbegin()->first!=(int)chr_num2name_.size()) exit_error("chomosome numbers are not continguous in "+SEQ_PATH+"chr_number");
		path_set = true;
	}
	
	string get_genome()
	{
		if (!path_set) exit_error("genepi path not set");
		return GenomeName;
	}
	
	void read_arguments(vector<string>& srce_opt, bool to_set_path)
	{
		// read arguments
		vector<int> to_rm(srce_opt.size(),0);
		for (size_t argi=1; argi<srce_opt.size(); ++argi)
		{
			if		(str_startsw(srce_opt[argi],"--genome"))		{	to_rm[argi]=1; ReadArg(srce_opt,argi,MainPath);	to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--pr2tx"))			{	to_rm[argi]=1; ReadArg(srce_opt,argi,NpNmFile);	to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--syn2sym"))		{	to_rm[argi]=1; ReadArg(srce_opt,argi,SynSymTr);	to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--tx2sym"))		{	to_rm[argi]=1; ReadArg(srce_opt,argi,RsqSymTr);	to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--tx2gid"))		{	to_rm[argi]=1; ReadArg(srce_opt,argi,txID_gID);	to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--read-through"))	{	to_rm[argi]=1; ReadArg(srce_opt,argi,read_thr);	to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--pseudo"))		{	to_rm[argi]=1; ReadArg(srce_opt,argi,pseudoid);	to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--main-tx-1"))		{	to_rm[argi]=1; ReadArg(srce_opt,argi,gID_MFOT);	to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--main-tx-2"))		{	to_rm[argi]=1; ReadArg(srce_opt,argi,uGID_can);	to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--gdb"))			{	to_rm[argi]=1; ReadArg(srce_opt,argi,GDB_name);	to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--seq"))			{	to_rm[argi]=1; ReadArg(srce_opt,argi,SEQ_PATH);	to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--genes"))			{	to_rm[argi]=1; ReadSet(srce_opt,argi,incl_g);	to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--txs"))			{	to_rm[argi]=1; ReadSet(srce_opt,argi,incl_t);	to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--canonical"))		{	to_rm[argi]=1; ReadArg(srce_opt,argi,canonic);	to_rm[argi]=1; }
			else if (str_startsw(srce_opt[argi],"--filt-tx"))		{	to_rm[argi]=1; ReadArg(srce_opt,argi,filt_tx);	to_rm[argi]=1; }
		}
		vector<string> dest_opt;
		for (size_t argi=0; argi<srce_opt.size(); ++argi) { if (!to_rm[argi]) dest_opt.push_back(srce_opt[argi]); }
		srce_opt=dest_opt;
		
		// set data. set_path should be the first.
		if (to_set_path) set_path(MainPath,"");
	}

	void add_help_text_var()
	{
		program.help_text_var("_Default_genome",MainPath);
		program.help_text_var("_Default_pr2tx",NpNmFile);
		program.help_text_var("_Default_syn2sym",SynSymTr);
		program.help_text_var("_Default_tx2sym",RsqSymTr);
		program.help_text_var("_Default_tx2gid",txID_gID);
		program.help_text_var("_Default_readthr",read_thr);
		program.help_text_var("_Default_pseudo",pseudoid);
		program.help_text_var("_Default_main_tx_1",gID_MFOT);
		program.help_text_var("_Default_main_tx_2",uGID_can);
		program.help_text_var("_Default_gdb",GDB_name);
		program.help_text_var("_Default_seq",SEQ_PATH);
		program.help_text_var("_Default_genes",str_of_container(incl_g,',',false));
		program.help_text_var("_Default_txs",str_of_container(incl_t,',',false));
		program.help_text_var("_Default_canonical",str_YesOrNo(canonic));
		program.help_text_var("_Default_filt_tx",str_YesOrNo(filt_tx));
	}
	
	string help_text() {
		string result="FOR GENOME\n";
		result += " --genome DIR      genome database location  {"+MainPath+"}\n";
		result += " --gdb S           gene database (refGeneLite, refgeneFull, ensGeneLite, ensGeneFull) {"+GDB_name+"}\n";
		result += " --seq P           sequence path {"+SEQ_PATH+"}\n";
		result += " --genes Ss        include these genes only {"+str_of_container(incl_g,',',false)+"}\n";
		result += " --txs Ss          include these transcripts only {"+str_of_container(incl_t,',',false)+"}\n";
//		result += " --canonical B     include canonical transcripts only {"+str_YesOrNo(canonic)+"}\n";
//		result += " --filt-tx B       filter transcripts internally {"+str_YesOrNo(filt_tx)+"}\n";
		return result;
	}
	
	// I should add if (!path_set) exit_error("genepi path not set"); but I take a chance that chr_num is cal by read_chr_num().
	bool is_chrM		(int chr_num) { return chr_num==chrM_num; }
	bool is_chrX		(int chr_num) { return chr_num==chrX_num; }
	bool is_chrY		(int chr_num) { return chr_num==chrY_num; }
	bool is_chrXY		(int chr_num) { return chr_num==chrXYnum; }
	bool is_autosomal	(int chr_num) { return chr_num>=chrA_min && chr_num<=chrA_max; }
	bool is_autoOrPAR	(int chr_num, int bp) { return expected_ploidy(chr_num,bp,2)==2; }
	bool is_MHC(int chr_num, int bp)
	{
		if (!MHC_chr || !MHC_beg || !MHC_end) return false;
		if (chr_num!=MHC_chr) return false;
		if (bp<MHC_beg) return false;
		if (bp>MHC_end) return false;
		return true;
	}
	int expected_ploidy(int chr_num, int bp, int sex) // input sex=1(m)/2(f); output ploidy=2/1/0(NA)
	{
		if (chr_num<=0)				return 0;	// bad <=0
		if (is_autosomal(chr_num))	return 2;	// auto (previsouly <=22)
		if (is_chrXY(chr_num))		return 2;	// chrXY (previsouly ==25)
		if (is_chrM(chr_num))		return 1;	// chrM (previsouly ==26)
		if (is_chrX(chr_num))					// chrX (previsouly ==23)
		{
			if (sex!=1 && sex!=2)	return 0;
			if (sex==2)				return 2;
			if (bp<  PAR1_x_b)		return 1;
			if (bp<= PAR1_x_e)		return 2;	// PAR1
			if (bp<  PAR2_x_b)		return 1;
			if (bp<= PAR2_x_e)		return 2;	// PAR2
			return 1;
		}
		if (is_chrY(chr_num))					// chrY (previsouly ==24)
		{
			if (sex!=1 && sex!=2) return 0;
			if (bp< PAR1_y_b) { if (sex==2) return 0; else return 1; }
			if (bp<=PAR1_y_e) { if (sex==2) return 2; else return 2; } // PAR1
			if (bp< PAR2_y_b) { if (sex==2) return 0; else return 1; }
			if (bp<=PAR2_y_e) { if (sex==2) return 2; else return 2; } // PAR2
							  { if (sex==2) return 0; else return 1; }
		}
		return 0;
	}

	bool read_chr_str(const std::string& input, std::string& cs)
	{
		if (!path_set) exit_error("genepi path not set");
		cs=to_upper_copy(input);
		if (str_startsw(cs, "CHROMOSOME_")) cs=cs.substr(11);
		if (str_startsw(cs, "CHROMOSOME"))	cs=cs.substr(10);
		if (str_startsw(cs, "CHROM_"))		cs=cs.substr(6);
		if (str_startsw(cs, "CHROM"))		cs=cs.substr(5);
		if (str_startsw(cs, "CHR_"))		cs=cs.substr(4);
		if (str_startsw(cs, "CHR"))			cs=cs.substr(3);
		if (str_has(cs,"P")) cs=substr_before_find(cs,"P");
		if (str_has(cs,"Q")) cs=substr_before_find(cs,"Q");
		if (str_has(cs,".")) cs=substr_before_find(cs,".");
		if (exist_element(chr_2StdName_,cs)) cs=chr_2StdName_[cs];
		if (exist_element(chr_name2num_,cs)) return true;
		else					{ cs=input; return false; }
	}
	
	// if success, return 1-26; otherwise, 0.
	int read_chr_num(const std::string& input)
	{
		string cs;
		if (read_chr_str(input,cs)) return chr_name2num_[cs];
		return 0;
	}

	bool read_chr_num(const std::string& input, int& chrNum)
	{
		chrNum = read_chr_num(input);
		return chrNum;
	}
	
	string convert_chr_num(int chr_num)
	{
		if (!path_set) exit_error("genepi path not set");
		if (exist_element(chr_num2name_,chr_num)) return chr_num2name_[chr_num];
		exit_error("unknown chromosome number "+itos(chr_num));
		return "impossible";
	}
	
	bool merge_region(int start, int end, vector< pair<int,int> >& regions)
	{
		bool merged=false;
		for (each_element(regions,it))
		{
			int& x = it->first;
			int& y = it->second;
			if (end>=x-1 && start<=y+1)
			{
				x = std::min(start,x);
				y = std::max(end,y);
				merged=true;
				break;
			}
		}
		if (!merged) regions.push_back(pair<int,int>(start,end));
		return merged;
	}
	
	void collapse_regions(vector< pair<int,int> >& regions)
	{
		for (int merged=1; merged; )
		{
			merged=0;
			vector< pair<int,int> > tmp;
			for (each_element(regions,itr))
			{
				int& start = itr->first;
				int& end = itr->second;
				merged += merge_region(start,end,tmp);
			}
			regions=tmp;
		}
	}
	
	bool ChrRegions::empty() {
		return (regions.empty() && !whole);
	}

	bool ChrRegions::not_read() {
		return (!did_setup);
	}

	bool ChrRegions::add(const std::vector<string>& in, bool zero_based, int padding, bool omit_unk) {
		string chr;
		int start=-1,end=-1;
		
		if (in.size()==1)
		{
			if (!str_has(in[0],":")) exit_error("Unkown file format for chromosome regions.");
			chr = substr_before_find(in[0],":");
			string s = substr_after_find(in[0],":");
			if (str_has(s,"-"))
			{
				if (!read_val_ge(substr_before_find(s,"-"),start,1)) exit_error("failed to read position "+substr_before_find(s,"-"));
				if (!read_val_ge(substr_after_find (s,"-"),end,  1)) exit_error("failed to read position "+substr_after_find(s,"-"));
			}
			else
			{
				if (!read_val_ge(s,start,1)) exit_error("failed to read position "+s);
				end = start;
			}
		}
		else if (in.size()==2)
		{
			chr = in[0];
			if (!read_val_ge(in[1],start,1)) exit_error("failed to read position "+in[1]);
			end = start;
		}
		else if (in.size()>=3)
		{
			chr = in[0];
			if (zero_based)
			{
				if (!read_val_ge(in[1],start,0)) exit_error("failed to read position "+in[1]);
				++start;
			}
			else
			{
				if (!read_val_ge(in[1],start,1)) exit_error("failed to read position "+in[1]);
			}
			if (!read_val_ge(in[2],end,1)) exit_error("failed to read position "+in[2]);
		}
		else
		{
			return false;
		}
		
		int chrNum = read_chr_num(chr); if (!chrNum) { if (omit_unk) return false; else exit_error("Can't read chr "+chr); }
		start -= padding; if (start<1) start=1;
		end   += padding;
		merge_region(start,end,unfinished[chrNum]);
		return true;
	}
	
	void ChrRegions::add(const int chrNum, const int start, const int end)
	{
		merge_region(start,end,unfinished[chrNum]);
	}
	
	void ChrRegions::finishing()
	{
		for (each_element(unfinished,itc))
		{
			const int& chr=itc->first;
			collapse_regions(itc->second);
			for (each_element(itc->second,itr)) regions[chr].insert(*itr);
		}
		whole=false;
		did_setup=true;
		unfinished.clear();
	}
	
	template <typename T>
	void ChrRegions::setup(const T& filename, bool zero_based, int padding, bool omit_unk) {
		regions.clear();
		unfinished.clear();
		tfile_format fmt;
		fmt.forbid_nf_rpt();
		fmt.set_delimiters("\t");
		int num_tracks=0;
		for (Rows_in_File(in,filename,&fmt))
		{
			if (str_startsw(in[0],"browser position")) continue;
			if (str_startsw(in[0],"track name=")) { ++num_tracks; continue; }
			add(in.contents(),zero_based,padding,omit_unk);
		}
		finishing();
		if (num_tracks>1) exit_error("There were multiple tracks in a BED file.");
	}
	template void ChrRegions::setup(const vector<string>& name, bool zero_based, int padding, bool omit_unk);
	template void ChrRegions::setup(const  deque<string>& name, bool zero_based, int padding, bool omit_unk);
	template void ChrRegions::setup(const    set<string>& name, bool zero_based, int padding, bool omit_unk);
	template void ChrRegions::setup(const         string& name, bool zero_based, int padding, bool omit_unk);
	
	void ChrRegions::write(std::ostream& file, bool zero_based) {
		for (each_element(regions,itc))
		{
			const int& chr=itc->first;
			for (each_element(itc->second,itr))
			{
				const int& start = itr->first;
				const int& end = itr->second;
				file << convert_chr_num(chr) << '\t' << (zero_based ? start-1 : start) << '\t' << end << '\n';
			}
		}
	}
	
	void ChrRegions::write(int desert_bp, const string& prefix)
	{
		int round=0;
		for (each_element(regions,itc))
		{
			int prev_bp=-1;
			std::fstream outf;
			if (!openfile_successfully(outf, prefix+itos(++round)+".bed", ios::out)) exit_error("Cannot open "+prefix+itos(round)+".bed");
			const int& chr=itc->first;
			for (each_element(itc->second,itr))
			{
				const int& start = itr->first;
				const int& end = itr->second;
				if (prev_bp==-1) prev_bp=end;
				if (start-prev_bp>desert_bp)
				{
					outf.close();
					if (!openfile_successfully(outf, prefix+itos(++round)+".bed", ios::out)) exit_error("Cannot open "+prefix+itos(round)+".bed");
				}
				prev_bp=end;
				outf << convert_chr_num(chr) << '\t' << start-1 << '\t' << end << '\n';
			}
			outf.close();
		}
	}
	
	bool ChrRegions::contain(const int chrNum, const int bp)
	{
		// if (chrNum<=0 || chrNum>MaxChrNum()) exit_error("chromosome number "+itos(chrNum)+" is invalid."); // no check to save time and the need for SEQ_PATH
		if (whole) return true;
		if (!exist_element(regions,chrNum)) return false;
		map<int,int>::iterator itr = regions[chrNum].upper_bound(bp);
		if (itr!=regions[chrNum].begin()) --itr;
		const int& start = itr->first;
		const int& end = itr->second;
		return (bp>=start && bp<=end);
	}
	
	bool ChrRegions::overlap(const int chrNum, const int bp1, const int bp2)
	{
		// if (chrNum<=0 || chrNum>MaxChrNum()) exit_error("chromosome number "+itos(chrNum)+" is invalid."); // no check to save time and the need for SEQ_PATH
		if (whole) return true;
		if (!exist_element(regions,chrNum)) return false;
		map<int,int>::iterator itr = regions[chrNum].upper_bound(bp2);
		if (itr!=regions[chrNum].begin()) --itr;
		const int& start = itr->first;
		const int& end = itr->second;
		return (bp2>=start && bp1<=end);
	}

	int ChrRegions::total_bp()
	{
		if (whole) return -1;
		int total=0;
		for (auto &chr:regions)
			for (auto &seg:chr.second)
				total += (seg.second-seg.first+1);
		return total;
	}
	
	// ------------------ genes ------------------

	map<string,gene_info>	gene_db;
	map< string, tx_set_t > gene_bySymbol;
	map< string, tx_set_t > gene_byTranscript;

	void gene_cal_hgvs(gene_info& g)
	{
		if (g.cdsStart == g.cdsEnd)
		{
			// return; // add on 2014-06-12 but changed to the following on 2015-09-01
			g.xSrelative.assign(g.exonCount, 0);
			g.xErelative.assign(g.exonCount, 0);
			if (g.strand=='+')
			{
				int previous=0;
				for (int i=0; i<g.exonCount; ++i)
				{
					g.xSrelative[i]=previous+1;
					previous += (g.exonEnds[i]-g.exonStarts[i]);
					g.xErelative[i]=previous;
				}
			}
			else
			{
				int previous=0;
				for (int i=g.exonCount-1; i>=0; --i)
				{
					g.xErelative[i]=previous+1;
					previous += (g.exonEnds[i]-g.exonStarts[i]);
					g.xSrelative[i]=previous;
				}
			}
		}
		else
		{
			int i;
			g.xSrelative.assign(g.exonCount, 0);
			g.xErelative.assign(g.exonCount, 0);
			if (g.strand=='+')
			{
				for (i=0; i<g.exonCount; ++i)
					if (g.cdsStart>=g.exonStarts[i] && g.cdsStart<g.exonEnds[i]) break;
				if (i==g.exonCount) exit_error(g.name+"'s cdsStart is not located inside any exon.");
				//																			//        -->>>>
				g.xSrelative[i] = g.exonStarts[i] - g.cdsStart;								// start  012345 : 0-2=-2
				g.xErelative[i] = g.exonEnds[i]   - g.cdsStart;								// end    123456 : 6-2=4
				if (g.xSrelative[i]==0) g.xSrelative[i]=1;									//        >>>>>>
				for (int j=i+1; j<g.exonCount; ++j)
				{																			//               >>>>>>
					g.xSrelative[j] = g.xErelative[j-1] + 1;								// 4+1=5
					g.xErelative[j] = g.xErelative[j-1] + (g.exonEnds[j] - g.exonStarts[j]);// 4+(6-0)=10
				}
				for (int j=i-1; j>=0; --j)
				{
					if (g.xSrelative[j+1]==1)												// ------ >>>>>>
					{																		//s2   e8
						g.xErelative[j] = -1;												//     -1
						g.xSrelative[j] = g.exonStarts[j] - g.exonEnds[j];					//-6 (2-8)
					}
					else																	// ------ -->>>>
					{
						g.xErelative[j] = g.xSrelative[j+1] - 1;								// -2-1=-3
						g.xSrelative[j] = g.xSrelative[j+1] - (g.exonEnds[j] - g.exonStarts[j]);// -2-(6-0)=-8
					}
				}
			}
			else
			{
				for (i=0; i<g.exonCount; ++i)
					if (g.cdsEnd>g.exonStarts[i] && g.cdsEnd<=g.exonEnds[i]) break;
				if (i==g.exonCount) exit_error(g.name+"'s cdsEnd is not located inside any exon.");
				//																			//       <<<<--
				g.xSrelative[i] = g.cdsEnd - g.exonStarts[i];								// start 012345 : 4-0=4
				g.xErelative[i] = g.cdsEnd - g.exonEnds[i];									// end   123456 : 4-6=-2
				if (g.xErelative[i]==0) g.xErelative[i]=1;									//       <<<<<<
				for (int j=i+1; j<g.exonCount; ++j)
				{
					if (g.xErelative[j-1]==1)												//       <<<<<< ------
					{																		//             s2   e8
						g.xSrelative[j] = -1;												//             -1
						g.xErelative[j] = g.exonStarts[j] - g.exonEnds[j];					//                  -6 (2-8)
					}
					else																	//       <<<<-- ------
					{
						g.xSrelative[j] = g.xErelative[j-1] - 1;										// -3 (-2-1)
						g.xErelative[j] = g.xErelative[j-1] - (g.exonEnds[j] - g.exonStarts[j]);		//      -8 (-2-(6-0))
					}
				}
				for (int j=i-1; j>=0; --j)
				{																			// <<<<<<
					g.xErelative[j] = g.xSrelative[j+1] + 1;								// 4+1=5
					g.xSrelative[j] = g.xSrelative[j+1] + (g.exonEnds[j] - g.exonStarts[j]);// 4+(6-0)=10
				}
			}
		}
	}
	
	int genomic_to_cds_location(const genepi::gene_info& g, int chr, int bp)
	{
		if (g.chrNumPlink==chr && g.txStart<bp && bp<=g.txEnd && g.cdsStart!=g.cdsEnd)
		{
			int i;
			for (i=0; i<g.exonCount; ++i)
				if (bp>g.exonStarts[i] && bp<=g.exonEnds[i]) break;
			if (i==g.exonCount) return 0; // intronic
			
			if (g.strand=='+')
			{
				if		(g.xSrelative[i]>0 && g.xErelative[i]>0) return g.xSrelative[i] + ( bp - g.exonStarts[i] -1 );
				else if (g.xSrelative[i]<0 && g.xErelative[i]<0) return g.xSrelative[i] + ( bp - g.exonStarts[i] -1 );
				else if (bp > g.cdsStart)						 return bp - g.cdsStart;
				else											 return bp - g.cdsStart -1;
			}
			else
			{
				if		(g.xSrelative[i]>0 && g.xErelative[i]>0) return g.xSrelative[i] - ( bp - g.exonStarts[i] -1 );
				else if (g.xSrelative[i]<0 && g.xErelative[i]<0) return g.xSrelative[i] - ( bp - g.exonStarts[i] -1 );
				else if (bp <= g.cdsEnd)						 return g.cdsEnd - bp +1;
				else											 return g.cdsEnd - bp;
			}
		}
		else
			return -std::numeric_limits<int>::max(); // not protein-coding OR not inside the transcript
	}
	
	int genomic_to_rna_location(const genepi::gene_info& g, int chr, int bp)
	{
		if (g.chrNumPlink==chr && g.txStart<bp && bp<=g.txEnd)
		{
			int i;
			if (g.strand=='+')
			{
				int previous=0;
				for (i=0; i<g.exonCount; ++i)
				{
					if (bp>g.exonStarts[i] && bp<=g.exonEnds[i]) return bp-g.exonStarts[i]+previous;
					previous += (g.exonEnds[i]-g.exonStarts[i]);
				}
			}
			else
			{
				int previous=1;
				for (i=g.exonCount-1; i>=0; --i)
				{
					if (bp>g.exonStarts[i] && bp<=g.exonEnds[i]) return g.exonEnds[i]-bp+previous;
					previous += (g.exonEnds[i]-g.exonStarts[i]);
				}
			}
			return 0; // intronic
		}
		else
			return -std::numeric_limits<int>::max(); // not inside the transcript
	}
	
	bool overlap_StoE(const genepi::gene_info& g, int chr, int b1, int b2, int up, int dn)
	{
		if (g.chrNumPlink!=chr) return false;
		int start;
		int end;
		if (g.strand=='-')	{ start = g.txStart - dn; end = g.txEnd + up; }
		else				{ start = g.txStart - up; end = g.txEnd + dn; }
		return (b1<=(start+1) && b2>=end);
	}

	bool overlap_gene(const genepi::gene_info& g, int chr, int b1, int b2, int up, int dn)
	{
		if (g.chrNumPlink!=chr) return false;
		int start;
		int end;
		if (g.strand=='-')	{ start = g.txStart - dn; end = g.txEnd + up; }
		else				{ start = g.txStart - up; end = g.txEnd + dn; }
		return (start<b2 && b1<=end);
	}

	bool overlap_exon(const genepi::gene_info& g, int chr, int b1, int b2, int up, int dn, int in)
	{
		if (g.chrNumPlink!=chr) return false;
		int start;
		int end;
		if (g.strand=='+')	{ start = g.txStart - up; end = g.txEnd + dn; }
		else				{ start = g.txStart - dn; end = g.txEnd + up; }
		if (start<b2 && b1<=end)
		{
			if (g.strand=='+')
				for (int i=0; i<g.exonCount; ++i)
				{
					if (i==0 && i==g.exonCount-1)	{ int e1=g.exonStarts[i]-up+1, e2=g.exonEnds[i]+dn; if (e1<=b2 && b1<=e2) return true; }
					else if	(i==0)					{ int e1=g.exonStarts[i]-up+1, e2=g.exonEnds[i]+in; if (e1<=b2 && b1<=e2) return true; }
					else if (i==g.exonCount-1)		{ int e1=g.exonStarts[i]-in+1, e2=g.exonEnds[i]+dn; if (e1<=b2 && b1<=e2) return true; }
					else							{ int e1=g.exonStarts[i]-in+1, e2=g.exonEnds[i]+in; if (e1<=b2 && b1<=e2) return true; }
				}
			else
				for (int i=g.exonCount-1; i>=0; --i)
				{
					if (i==0 && i==g.exonCount-1)	{ int e1=g.exonStarts[i]-dn+1, e2=g.exonEnds[i]+up; if (e1<=b2 && b1<=e2) return true; }
					else if	(i==0)					{ int e1=g.exonStarts[i]-dn+1, e2=g.exonEnds[i]+in; if (e1<=b2 && b1<=e2) return true; }
					else if (i==g.exonCount-1)		{ int e1=g.exonStarts[i]-in+1, e2=g.exonEnds[i]+up; if (e1<=b2 && b1<=e2) return true; }
					else							{ int e1=g.exonStarts[i]-in+1, e2=g.exonEnds[i]+in; if (e1<=b2 && b1<=e2) return true; }
				}
			return false;
		}
		else return false;
	}
	// to test this function:
	// cout<<genepi::overlap_exon(genepi::gene_byTranscript["NM_000059"].begin()->second,13,32889581,32889582,up_reg,dn_reg,in_reg)<<endl; // 1
	// cout<<genepi::overlap_exon(genepi::gene_byTranscript["NM_000059"].begin()->second,13,32889580,32889581,up_reg,dn_reg,in_reg)<<endl; // 0
	// cout<<genepi::overlap_exon(genepi::gene_byTranscript["NM_000059"].begin()->second,13,32889816,32889817,up_reg,dn_reg,in_reg)<<endl; // 1
	// cout<<genepi::overlap_exon(genepi::gene_byTranscript["NM_000059"].begin()->second,13,32889817,32889818,up_reg,dn_reg,in_reg)<<endl; // 0
	// cout<<genepi::overlap_exon(genepi::gene_byTranscript["NM_007294"].begin()->second,17,41196311,41196312,up_reg,dn_reg,in_reg)<<endl; // 1
	// cout<<genepi::overlap_exon(genepi::gene_byTranscript["NM_007294"].begin()->second,17,41196310,41196311,up_reg,dn_reg,in_reg)<<endl; // 0
	// cout<<genepi::overlap_exon(genepi::gene_byTranscript["NM_007294"].begin()->second,17,41197831,41197832,up_reg,dn_reg,in_reg)<<endl; // 1
	// cout<<genepi::overlap_exon(genepi::gene_byTranscript["NM_007294"].begin()->second,17,41197832,41197833,up_reg,dn_reg,in_reg)<<endl; // 0
	
	bool overlap_cds_(const genepi::gene_info& g, int chr, int b1, int b2)
	{
		if (g.chrNumPlink!=chr) return false;
		if (g.cdsStart==g.cdsEnd) return false;
		if (g.cdsStart<b2 && b1<=g.cdsEnd)
		{
			for (int i=0; i<g.exonCount; ++i)
			{
				if (g.cdsStart>=g.exonEnds[i])  continue;// no CDS in this exon
				if (g.cdsEnd<g.exonStarts[i]+1) continue;// no CDS in this exon
				int e1 = std::max(g.exonStarts[i],g.cdsStart)+1;
				int e2 = std::min(g.exonEnds[i],g.cdsEnd);
				if (e1<=b2 && b1<=e2) return true;
			}
			return false;
		}
		else return false;
	}

	// This position in within the last exon of gene g. Return 1(yes) 0(no)
	bool in_last_exon(const genepi::gene_info& g, int chr, int bp)
	{
		if (g.chrNumPlink==chr && g.txStart<bp && bp<=g.txEnd)
		{
			int i;
			for (i=0; i<g.exonCount; ++i)
				if (bp>g.exonStarts[i] && bp<=g.exonEnds[i]) break;
			if (i==g.exonCount) return false; // intronic
			
			if (g.strand=='+')	return (++i==g.exonCount);
			else				return (++i==1);
		}
		else
			return false; // not in the gene
	}

	// This position in within the first exon of gene g. Return 1(yes) 0(no)
	bool in_first_exon(const genepi::gene_info& g, int chr, int bp)
	{
		if (g.chrNumPlink==chr && g.txStart<bp && bp<=g.txEnd)
		{
			int i;
			for (i=0; i<g.exonCount; ++i)
				if (bp>g.exonStarts[i] && bp<=g.exonEnds[i]) break;
			if (i==g.exonCount) return false; // intronic
			
			if (g.strand=='+')	return (++i==1);
			else				return (++i==g.exonCount);
		}
		else
			return false; // not in the gene
	}

	// return 1 (success) 2 (failed)
	bool nearest_copy(int chr, int bp, const string& TxID, gene_info* g)
	{
		g = NULL;
		int nearest = std::numeric_limits<int>::max();
		for (auto &tx:gene_byTranscript[TxID])
		{
			if (tx.second.chrNumPlink!=chr) continue;
			int Distance = std::abs(bp-tx.second.txStart-1);
			if (overlap_gene(tx.second,chr,bp,bp)) Distance=0;
			if (Distance < nearest) { g = &tx.second; nearest=Distance; }
			// do nothing if ==, so that it choose the one with the longest cdsLen / txLen
		}
		return g!=NULL;
	}
	
	// xS, xE are distances to exonStart exonEnd, independent of strand. Positive if inside an exon, negative if upstream to Tx.
	int cds_to_genomic_location(const gene_info& g, int loc, string& chr, int& bp, int& xS, int& xE)
	{
		chr="Error"; bp=-1;
		if (loc==0) return 0; // input wrong: there is no nucleotide 0
		if (g.cdsStart==g.cdsEnd) return 1;
		if (g.strand=='+')
		{
			int i=0; // just to avoid warning: ‘i’ may be used uninitialized
			if (loc<g.xSrelative[0])
			{
				i=0;
				if (g.xSrelative[i]==1) loc+=1; // unlikely
			}
			else
			{
				for (i=0; i<g.exonCount; ++i) if (loc >= g.xSrelative[i] && loc <= g.xErelative[i]) break;
				if (i==g.exonCount) return 0; // input wrong
			}
			xS = loc - g.xSrelative[i];
			xE = g.xErelative[i] - loc;
			chr=g.chr;
			if		(g.xSrelative[i]>0 && loc>0) bp = g.exonStarts[i] + (loc - g.xSrelative[i]) + 1;
			else if (g.xSrelative[i]<0 && loc<0) bp = g.exonStarts[i] + (loc - g.xSrelative[i]) + 1;
			else								 bp = g.exonStarts[i] + (loc - g.xSrelative[i]);
		}
		else
		{
			int i=g.exonCount-1; // just to avoid warning: ‘i’ may be used uninitialized
			if (loc<g.xErelative[g.exonCount-1])
			{
				i=g.exonCount-1;
				if (g.xErelative[g.exonCount-1]==1) loc+=1; // unlikely
			}
			else
			{
				for (i=0; i<g.exonCount; ++i) if (loc <= g.xSrelative[i] && loc >= g.xErelative[i]) break;
				if (i==g.exonCount) return 0; // input wrong
			}
			xS = g.xSrelative[i] - loc;
			xE = loc - g.xErelative[i];
			chr=g.chr;
			if		(g.xSrelative[i]>0 && loc>0) bp = g.exonStarts[i] - (loc - g.xSrelative[i]) + 1;
			else if (g.xSrelative[i]<0 && loc<0) bp = g.exonStarts[i] - (loc - g.xSrelative[i]) + 1;
			else								 bp = g.exonStarts[i] - (loc - g.xSrelative[i]);
		}
		return 1;
	}

	int rna_to_genomic_location(const gene_info& g, int loc, string& chr, int& bp, int& xS, int& xE)
	{
		chr="Error"; bp=-1;
		if (loc==0) return 0; // input wrong: there is no nucleotide 0
		int i;
		int lb_size=0;
		int ub_size=0;
		if (g.strand=='+')
		{
			if (loc<0)
			{
				i=0;
				loc+=1;
				lb_size = 0;
				ub_size = g.exonEnds[i]-g.exonStarts[i];
			}
			else
			{
				for (i=0; i<g.exonCount; ++i)
				{
					lb_size = ub_size;
					ub_size += g.exonEnds[i]-g.exonStarts[i];
					if (loc <= ub_size) break;
				}
				if (i==g.exonCount) return 0; // input wrong
			}
			chr= g.chr;
			xS = loc - lb_size -1;
			xE = ub_size - loc;
			bp = g.exonStarts[i] + (loc - lb_size);
		}
		else
		{
			if (loc<0)
			{
				i=g.exonCount-1;
				loc+=1;
				lb_size = 0;
				ub_size = g.exonEnds[i]-g.exonStarts[i];
			}
			else
			{
				for (i=g.exonCount-1; i>=0; --i)
				{
					lb_size = ub_size;
					ub_size += g.exonEnds[i]-g.exonStarts[i];
					if (loc <= ub_size) break;
				}
				if (i<0) return 0; // input wrong
			}
			chr= g.chr;
			xS = ub_size - loc;
			xE = loc - lb_size -1;
			bp = g.exonStarts[i] + (ub_size - loc) + 1;
		}
		return 1;
	}
	
	// ---------------- begin of about chromosomes ----------------
	
	typedef unsigned char BYTE;
	map<int,int>						_chrlen_bp;	// chrlen_bp[chrNum], chrNum is in PLINK convention: 1-22,X,Y,XY,M
	map<int,int>						_chrlen_bp_cumulative;
	int									_MaxChrNum=-1;
	int									_GenomeLen=-1;
	std::vector< std::vector<BYTE> >	_chr_ref;
	std::vector<int>					_ref_read;
	ChrRegions							_chr_Ns;
	std::mutex							_chr_mt;	// mutex

	void _prepare_chrlen()
	{
		if (!path_set) exit_error("genepi path not set");
		_chr_mt.lock();
		if (_chrlen_bp.empty())
		{
			filename_change_home_path(SEQ_PATH);
			if (SEQ_PATH.empty()) exit_error("path/to/sequence/ not defined.");
			if (!DirExists(SEQ_PATH)) exit_error(SEQ_PATH+" not exist.");
			if (!str_endsw(SEQ_PATH,"/")) SEQ_PATH+='/';
			_MaxChrNum=-1;
			for (Rows_in_File(in,SEQ_PATH+"chr_len",2))
			{
				int c; if (!read_chr_num(in[0],c)) exit_error("unknown chromosome in "+SEQ_PATH+"chr_len");
				int l; if (!read_val_ge(in[1],l,1)) exit_error("cannot read "+in[1]+" in "+SEQ_PATH+"chr_len as basepair");
				_chrlen_bp[c]=l;
				if (c>_MaxChrNum) _MaxChrNum=c;
			}
			_GenomeLen=0;
			for (auto &c:_chrlen_bp)
			{
				_chrlen_bp_cumulative[c.first]=_GenomeLen;
				_GenomeLen += c.second;
			}
			_chr_ref.assign(_MaxChrNum,std::vector<BYTE>());
			_ref_read.assign(_MaxChrNum,0);
			_chr_Ns.setup(SEQ_PATH+"Ns",false,0,false);
		}
		_chr_mt.unlock();
	}
	
	int chrlen_bp           (int chrNum) { if (_chrlen_bp.empty()) _prepare_chrlen(); return _chrlen_bp[chrNum]; }
	int chrlen_bp_cumulative(int chrNum) { if (_chrlen_bp.empty()) _prepare_chrlen(); return _chrlen_bp_cumulative[chrNum]; }
	int total_genome_len()				 { if (_chrlen_bp.empty()) _prepare_chrlen(); return _GenomeLen; }
	int MaxChrNum()						 { if (_chrlen_bp.empty()) _prepare_chrlen(); return _MaxChrNum; }
	
	// Reads DNA sequences stored in files created by the program read_fasta.
	// http://stackoverflow.com/questions/15138353/reading-the-binary-file-into-the-vector-of-unsigned-chars
	string DNA_seq(int chrNum, int bp, int len)
	{
		if (_chrlen_bp.empty()) _prepare_chrlen();
		
		// check errors
		if (chrNum<=0 || chrNum>MaxChrNum()) exit_error("chromosome number "+itos(chrNum)+" is invalid.");
		if (bp<1) exit_error("basepair position must be greater than zero.");
		if (len<1) exit_error("DNA_seq query length must be at least 1 bp.");
		if (bp>chrlen_bp(chrNum)) exit_error("query position ("+itos(bp)+") exceeds the length ("+itos(chrlen_bp(chrNum))+") of chromosome "+itos(chrNum));
		if (bp+len-1>chrlen_bp(chrNum)) exit_error("query beyond end of chromosome");
		
		// prepare for this chr
		int chr_ID = chrNum-1;
		if (!_ref_read[chr_ID])
		{
			_chr_mt.lock();
			if (!_ref_read[chr_ID])
			{
				string chrName = convert_chr_num(chrNum);
				ifstream file(SEQ_PATH+chrName, ios::in|ios::binary);
				file.unsetf(std::ios::skipws);
				std::streampos fileSize;
				file.seekg(0, std::ios::end);
				fileSize = file.tellg();
				file.seekg(0, std::ios::beg);
				_chr_ref[chr_ID].reserve(fileSize);
				_chr_ref[chr_ID].insert(_chr_ref[chr_ID].begin(),
										std::istream_iterator<BYTE>(file),
										std::istream_iterator<BYTE>());
				_ref_read[chr_ID]=1;
			}
			_chr_mt.unlock();
		}
		
		int loc = (bp-1)%4;		// allele code is in which bit  within the byte = 0/1/2/3
		int byte =(bp-1)/4;		// allele code is in which byte within the data = 0,1,2.. + start_coordinate
		string result;
		for (int i=0;i<len;++i)
		{
			if (_chr_Ns.contain(chrNum,bp))
				result+='N';
			else
			{
				char a = _chr_ref[chr_ID][byte];
				switch (loc) {
					case 0: a>>=6; break;
					case 1: a>>=4; break;
					case 2: a>>=2; break;
					case 3:        break;
					default: exit_error("wrong loc"); break;
				}
				a &= 0x03;
				switch (a) {
					case 0x00: result+='A'; break;
					case 0x01: result+='C'; break;
					case 0x02: result+='G'; break;
					case 0x03: result+='T'; break;
					default: exit_error("wrong a"); break;
				}
			}
			if (++loc>3) { loc=0; ++byte; }
			++bp;
		}
		return result;
	}
	
	string cytoband(const int chrNum, const int bp1, const int bp2)
	{
		if (_chrlen_bp.empty()) _prepare_chrlen();
		
		static std::map< int, std::map< int, std::pair<int,string> > >	regions;
		static std::mutex mt;
		if (regions.empty())
		{
			mt.lock();
			if (regions.empty())
			{
				if (!find_file_zippedOrNot(SEQ_PATH+"cytoBand.txt.gz")) exit_error("cytoBand.txt.gz not found.");
				for (Rows_in_File(in,SEQ_PATH+"cytoBand.txt.gz",4))
				{
					if (str_startsw(in[0],"chr")) in[0]=in[0].substr(3);
					int in_chrNum = read_chr_num(in[0]); if (!in_chrNum) { exit_error("Can't read chr "+in[0]); }
					int in_bp1; if (!read_val_ge(in[1],in_bp1,0)) exit_error("cannot read bp1 "+in[1]); ++in_bp1;
					int in_bp2; if (!read_val_ge(in[2],in_bp2,1)) exit_error("cannot read bp2 "+in[2]);
					regions[in_chrNum][in_bp1]=make_pair(in_bp2,in[0]+in[3]);
				}
				if (regions.empty()) exit_error("nothing was read from "+SEQ_PATH+"cytoBand.txt.gz");
			}
			mt.unlock();
		}

		
		if (!exist_element(regions,chrNum)) return "NA";
		string str1;
		{
			map<int,std::pair<int,string>>::iterator itr = regions[chrNum].upper_bound(bp1);
			if (itr!=regions[chrNum].begin()) --itr;
			const int& start = itr->first;
			const int& end = itr->second.first;
			if (bp1>=start && bp1<=end) str1=itr->second.second;
		}
		string str2;
		{
			map<int,std::pair<int,string>>::iterator itr = regions[chrNum].upper_bound(bp2);
			if (itr!=regions[chrNum].begin()) --itr;
			const int& start = itr->first;
			const int& end = itr->second.first;
			if (bp2>=start && bp2<=end) str2=itr->second.second;
		}
		if (str1.empty() && str2.empty()) return "NA";
		if (str1==str2)		return str1;
		if (str1.empty())	return str2;
		if (str2.empty())	return str1;
		return str1+"-"+str2;
	}
	
	// ---------------- end of about chromosomes ----------------
	
	string RNA_seq(const genepi::gene_info& g, bool convert_to_U)
	{
		string result;
		for (size_t i=0;i<g.exonStarts.size();++i)
			result += DNA_seq(g.chrNumPlink,g.exonStarts[i]+1,g.exonEnds[i]-g.exonStarts[i]);
		if (g.strand=='-') genepi::dna_rc(result,true);
		if (convert_to_U)
		{
			std::replace( result.begin(), result.end(), 't', 'u' );
			std::replace( result.begin(), result.end(), 'T', 'U' );
		}
		return result;
	}

	string UT5_seq(const genepi::gene_info& g, bool convert_to_U)
	{
		if (g.cdsStart == g.cdsEnd) return string();
		return RNA_seq(g,convert_to_U).substr(0,g.u5Len);
	}
	
	string UT3_seq(const genepi::gene_info& g, bool convert_to_U)
	{
		if (g.cdsStart == g.cdsEnd) return string();
		return RNA_seq(g,convert_to_U).substr(g.u5Len+g.cdsLen);
	}

	string CDS_seq(const genepi::gene_info& g)
	{
		if (g.cdsStart == g.cdsEnd) return string();
		string cds;
		for (int i=0; i<g.exonCount; i++)
		{
			int s = g.exonStarts[i];
			int e = g.exonEnds[i];
			if (g.cdsStart > s)	s = g.cdsStart;
			if (g.cdsEnd < e)	e = g.cdsEnd;
			if (s>=e) continue; // not s>e  , assuming all start coordinates are 0-based (UCSC).
			cds += DNA_seq(g.chrNumPlink,s+1,e-s);
		}
		if (g.strand=='-') genepi::dna_rc(cds,true);
		return cds;
	}
	
	string AA_seq(const genepi::gene_info& g)
	{
		if (g.cdsStart == g.cdsEnd) return string();
		string pr = translate_cds(CDS_seq(g),true);
		if (pr.empty()) { lns << showe << "Problem in translating CDS for "+g.name+" "+g.name2+": empty CDS" << flush_logger; return string(); }
		if (pr[pr.size()-1]=='*') pr.pop_back(); // It's possible if cdsEnd is "incmpl"
		if (pr.empty()) { lns << showe << "Problem in translating CDS for "+g.name+" "+g.name2+": only *" << flush_logger; return string(); }
		if (str_startsw(g.name,"NM_")||str_startsw(g.name2,"NM_")) pr[0]='M'; // Always starts with M even codes for different AA
		return pr;
	}
	// on 2014-07-01 I found one problem: Problem in translating CDS for NM_001010890 PRAMEF9: only *

	string Syn2Sym(const string& input, bool NoChangeIfNotExist)
	{
		if (!path_set) exit_error("genepi path not set");
		static TranslateVariableSameType<string> tr_Syn2Sym;
		static std::mutex mt;
		if (tr_Syn2Sym.reference_data.empty())
		{
			mt.lock();
			if (tr_Syn2Sym.reference_data.empty())
			{
				if (SynSymTr.empty()) exit_error("NCBI_Synonyms_to_Symbol_translation_file not defined.");
				if (!find_file_zippedOrNot(SynSymTr)) exit_error("NCBI_Synonyms_to_Symbol_translation_file not exist.");
				tfile_format format;
				format.set_delimiters("\t");
				for (Rows_in_File(in,SynSymTr,&format))
				{
					if (exist_element(tr_Syn2Sym.reference_data,in[0]) && tr_Syn2Sym.reference_data[in[0]]!=in[1])
						exit_error("ambiguous record for "+in[0]+" in "+SynSymTr);
					tr_Syn2Sym.reference_data[in[0]]=in[1];
				}
				if (tr_Syn2Sym.reference_data.empty()) exit_error("nothing was read from "+SynSymTr);
			}
			mt.unlock();
		}
		if		(NoChangeIfNotExist)		return tr_Syn2Sym.tr_if_exist(input);
		else if (tr_Syn2Sym.exist(input))	return tr_Syn2Sym.tr(input);
		else								return "";
	}
	
	string Rsq2Sym(const string& input, bool NoChangeIfNotExist)
	{
		if (!path_set) exit_error("genepi path not set");
		static TranslateVariableSameType<string> tr_Rsq2Sym;
		static std::mutex mt;
		if (tr_Rsq2Sym.reference_data.empty())
		{
			mt.lock();
			if (tr_Rsq2Sym.reference_data.empty())
			{
				if (RsqSymTr.empty()) exit_error("_NCBI_refSeq_to_Symbol_translation_file not defined.");
				if (!find_file_zippedOrNot(RsqSymTr)) exit_error("_NCBI_refSeq_to_Symbol_translation_file not exist.");
				tfile_format format;
				format.set_delimiters("\t");
				for (Rows_in_File(in,RsqSymTr,&format))
				{
					if (exist_element(tr_Rsq2Sym.reference_data,in[0]) && tr_Rsq2Sym.reference_data[in[0]]!=in[1])
						exit_error("ambiguous record for "+in[0]+" in "+RsqSymTr);
					tr_Rsq2Sym.reference_data[in[0]]=in[1];
				}
				if (tr_Rsq2Sym.reference_data.empty()) exit_error("nothing was read from "+RsqSymTr);
			}
			mt.unlock();
		}
		if		(NoChangeIfNotExist)		return tr_Rsq2Sym.tr_if_exist(input);
		else if (tr_Rsq2Sym.exist(input))	return tr_Rsq2Sym.tr(input);
		else								return "";
	}
	
	string NP_2_NM(const string& input, bool NoChangeIfNotExist)
	{
		if (!path_set) exit_error("genepi path not set");
		static TranslateVariableSameType<string> tr_NP_2_NM;
		static std::mutex mt;
		if (tr_NP_2_NM.reference_data.empty())
		{
			mt.lock();
			if (tr_NP_2_NM.reference_data.empty())
			{
				if (NpNmFile.empty()) exit_error("NCBI_NP_to_NM_translation_file not defined.");
				if (!find_file_zippedOrNot(NpNmFile)) exit_error("NCBI_NP_to_NM_translation_file not exist.");
				tr_NP_2_NM.setup(NpNmFile);
				if (tr_NP_2_NM.reference_data.empty()) exit_error("nothing was read from "+NpNmFile);
			}
			mt.unlock();
		}
		if		(NoChangeIfNotExist)		return tr_NP_2_NM.tr_if_exist(input);
		else if (tr_NP_2_NM.exist(input))	return tr_NP_2_NM.tr(input);
		else								return "";
	}

	// input: code remain g; output: chr loc info
	int read_hgvs_position(const char& code, string& remain, gene_info& g, string& chr, int& loc, string& info)
	{
		int loc_Tx=0;	// location within transcript
		int in_3prime=0;// offset in 3' UTR (c.*) or gene flanking region (r.*).
		int in_intron=0;// offset in intron. These two will not be both non-zero.
		
		if (!remain.empty())
		{
			if (str_startsw(remain,"IVS")) // outdated HGVS, but program it anyway.
			{
				remain.erase(0,3);
				size_t intron_number;
				if (!extract_int(remain,intron_number)) { info+=":ErrIvsNumber(NaN)"; return 0; }
				if (intron_number >= g.exonStarts.size()) { info+=":ErrIvsNumber"; return 0; }
				if (remain.empty()) { info+=":ErrIvsOffset"; return 0; }
				if (remain[0]!='+'&&remain[0]!='-') { info+=":ErrIvsOffset"; return 0; }
				if (!extract_int(remain,in_intron)) { info+=":ErrIvsOffset(NaN)"; return 0; }
				if (in_intron==0) { info+=":ErrIvsOffset"; return 0; }
				if (g.strand=='-')
				{
					if (in_intron>0)	loc_Tx = g.xSrelative[g.exonStarts.size()-intron_number];
					else				loc_Tx = g.xErelative[g.exonStarts.size()-intron_number-1];
				}
				else
				{
					if (in_intron>0)	loc_Tx = g.xErelative[intron_number-1];
					else				loc_Tx = g.xSrelative[intron_number];
				}
			}
			else if (remain[0]=='*')
			{
				if		(code=='c')				{ if (g.strand=='-') loc = g.cdsStart+1; else loc = g.cdsEnd; }
				else if (code=='r'||code=='n')	{ if (g.strand=='-') loc = g.txStart+1;  else loc = g.txEnd; }
				remain.erase(0,1);
				if (!extract_int(remain,in_3prime)) { info+=":Err*Offset(NaN)"; return 0; }
				if (in_3prime==0) { info+=":Err*Offset(zero)"; return 0; }
				if (in_3prime<0) { info+=":Err*Offset(negative)"; return 0; } // unlikely but legitimate. remember to ++in_3prime before adding a negative in_3prime.
				if (g.strand=='-') in_3prime*=-1;
				loc += in_3prime;
				chr = g.chr;
				return 1;
			}
			else
			{
				if (!extract_int(remain,loc_Tx)) { info+=":ErrPos(NaN)"; return 0; }
				if (!remain.empty() && (remain[0]=='+'||remain[0]=='-'))
					if (!extract_int(remain,in_intron)) { info+=":ErrPos(NaN)"; return 0; }
			}
		}
		else
		{
			info+=":ErrPos";
			return 0;
		}
		
		if		(code=='g'||code=='m')
		{
			chr=g.chr;
			loc=loc_Tx;
		}
		else if (code=='r'||code=='n')
		{
			int xS, xE;
			if (!rna_to_genomic_location(g,loc_Tx,chr,loc,xS,xE))
			{
				info+=":ErrRNA";
				return 0;
			}
			if (in_intron>0 && g.strand!='-' && xE!=0) { info+=":ErrExonBoundary"; return 0; }
			if (in_intron>0 && g.strand=='-' && xS!=0) { info+=":ErrExonBoundary"; return 0; }
			if (in_intron<0 && g.strand!='-' && xS!=0) { info+=":ErrExonBoundary"; return 0; }
			if (in_intron<0 && g.strand=='-' && xE!=0) { info+=":ErrExonBoundary"; return 0; }
			if (g.strand=='-') in_intron*=-1; loc += in_intron;
		}
		else if (code=='c')
		{
			int xS, xE;
			if (!cds_to_genomic_location(g,loc_Tx,chr,loc,xS,xE))
			{
				info+=":ErrCDS";
				return 0;
			}
			if (in_intron>0 && g.strand!='-' && xE!=0) { info+=":ErrExonBoundary"; return 0; }
			if (in_intron>0 && g.strand=='-' && xS!=0) { info+=":ErrExonBoundary"; return 0; }
			if (in_intron<0 && g.strand!='-' && xS!=0) { info+=":ErrExonBoundary"; return 0; }
			if (in_intron<0 && g.strand=='-' && xE!=0) { info+=":ErrExonBoundary"; return 0; }
			if (g.strand=='-') in_intron*=-1; loc += in_intron;
		}
		else
		{
			info+=":ErrLvl("+s(code)+")"; // n. p.
			return 0;
		}
		return 1;
	}
	
	bool is_RefSeq(const string& tx)
	{
		return (str_startsw(tx,"NP_") || str_startsw(tx,"XP_") || str_startsw(tx,"NM_") || str_startsw(tx,"NR_") || str_startsw(tx,"XM_") || str_startsw(tx,"XR_")) ;
	}
	
	// remove .version_number from transcript.version_number
	string rmTxVer_copy(const string& tx)
	{
		size_t found = tx.find('.');
		if (found!=string::npos) return tx.substr(0,found);
		else return tx;
	}
	void rmTxVer(string& tx)
	{
		tx = rmTxVer_copy(tx);
	}
	
	void RNA2DNA(string& seq)
	{
		std::replace( seq.begin(), seq.end(), 'u', 't');
		std::replace( seq.begin(), seq.end(), 'U', 'T');
	}
	
	string RNA2DNA_copy(const string& seq)
	{
		string cp = seq;
		RNA2DNA(cp);
		return cp;
	}
	
	// only the first instance of gene_byTranscript is chosen, but RefSeq is confusing in this regards and there's no solutions.
	// Require: nucleotide must be in upper case; c./p./etc must be in lower case; No space; no unknown position;
	int read_hgvs(const string& input, string& gene, char& code, string& chr, int& stt, int& end, string& fr, string& to, string& info)
	{
		// gene_byTranscript must be already in memory. Requires RefSeq but no checking.
		if (gene_byTranscript.empty())						exit_error("read_hgvs() cannot work without read_refGene().");
		if (!is_RefSeq(gene_byTranscript.begin()->first))	exit_error("read_hgvs() cannot work without read_refGene().");
		
		// initial value
		gene.clear(); code='-'; chr="Error"; fr="."; to="."; stt=-1; end=-1;
		string variant; // things after :
		std::tuple<int,int,int,int> priority (0,0,0,0);

		// read gene & variant
		vector<string> transcript;
		boost::split(transcript,input,boost::is_any_of(","));
		for (auto &hgvs:transcript)
		{
			// prepare
			if (hgvs.empty()) continue;
			vector<string> ifs;
			boost::split(ifs,hgvs,boost::is_any_of(":"));
			if (ifs.size()<=1) continue;
			
			// read gene. Priority: RefSeq > Other Tx > Tx from Symbol > the first field. This list is best if gene_byTranscript is built by read_refGene().
			string	rd_g;	// read gene
			int		src =0; // source of this reading
			if (rd_g.empty()) { for (auto &i:ifs) if (is_RefSeq(i))											{ rd_g=NP_2_NM(rmTxVer_copy(i),true);		  src=3; break; } }
			if (rd_g.empty()) { for (auto &i:ifs) if (exist_element(gene_byTranscript,rmTxVer_copy(i)))		{ rd_g=rmTxVer_copy(i);						  src=2; break; } }
			if (rd_g.empty()) { for (auto &i:ifs) if (exist_element(gene_bySymbol,i))						{ rd_g=gene_bySymbol[i].begin()->second.name; src=1; break; } }
			if (rd_g.empty()) { for (auto &i:ifs) { i=Syn2Sym(i,true); if (exist_element(gene_bySymbol,i))	{ rd_g=gene_bySymbol[i].begin()->second.name; src=1; break; } } }
			if (rd_g.empty()) { rd_g=rmTxVer_copy(ifs[0]); }
			
			// read variant. Priority: c.xxx === n.xxx >> p.xxx == r.xxx >> g.xxx == m.xxx
			string rd_v;
			for (auto &i:ifs)
			{
				if (str_startsw(i,"c.")) {	rd_v=i; break; }
				if (str_startsw(i,"n.")) {	rd_v=i; break; }
				if (str_startsw(i,"p.")) 	rd_v=i;
				if (str_startsw(i,"r."))	rd_v=i;
				if (str_startsw(i,"g.")) {	if (rd_v.empty()) rd_v=i; }
				if (str_startsw(i,"m.")) {	if (rd_v.empty()) rd_v=i; }
			}
			
			// store results if successful. Use RefSeq whenever possible.
			if (!rd_g.empty() && !rd_v.empty())
			{
				int cdsLen = 0;	if (exist_element(gene_byTranscript,rd_g)) cdsLen =		gene_byTranscript[rd_g].begin()->second.cdsLen;
				int txLen = 0;	if (exist_element(gene_byTranscript,rd_g)) txLen  =		gene_byTranscript[rd_g].begin()->second.txLen;
				int idLen = 0;	if (exist_element(gene_byTranscript,rd_g)) idLen  = 100-gene_byTranscript[rd_g].begin()->second.name.size();
				std::tuple<int,int,int,int> this_priority (src,cdsLen,txLen,idLen);
				if (this_priority>priority || gene.empty()) { gene = rd_g; variant = rd_v; priority = this_priority; }
			}
		}
		info = gene + ":" + variant;

		// read code and remain
		if (variant.size()<3)	{ info+=":ErrHGVS"; return 0; }
		if (variant[1]!='.')	{ info+=":ErrHGVS"; return 0; }
		code = variant[0];
		string remain = variant.substr(2);
		if (str_has(remain,"?")) { info+=":ErrQMrk"; return 0; } // does not support question marks: BRCA1:c.-232-?_134+?dup
		if (remain[0]=='[') { info+=":Err2+ch"; return 0; } // does not support multiple changes: c.[76A>C; 83G>C] / c.[76A>C];[83G>C] or masaicism c.[83G=/83G>C] or Chimerism c.[=//83G>C]
		if (remain[0]=='t') { info+=":ErrTrLoc"; return 0; } // does not support translocations: t(X;4)(p21.2;q35)(c.857+101_857+102)

		if (!exist_element(gene_byTranscript,gene)) { info+=":ErrTx"; return 0; } // doesn't work for m.xxx, which doesn't need a gene
		if (gene_byTranscript[gene].size()>1)		{ info+=":DupTx"; return 0; }
		gene_info& g = gene_byTranscript[gene].begin()->second;
		info = g.name2 + ":" + info;
		
		// read Start End Ref Alt
		switch (code) {
			case 'r':
			case 'n':
			case 'm':
			case 'g':
			case 'c': {
				// NM_000051:c.G8158C / NM_000051:c.G-5C / NM_000051:c.G4909+1A / NM_000051:c.G8851-1T / NM_000051:c.C*5T
				if (remain[0]=='A'||remain[0]=='C'||remain[0]=='G'||remain[0]=='T'|| (remain[0]=='U' && (code=='r'||code=='n')) )
				{
					string rd;
					if (!remain.empty()) rd = RNA2DNA_copy(s(extract_char(remain))); else { info+=":ErrImpossible"; return 0; }
					if (!read_hgvs_position(code,remain,g,chr,stt,info)) return 0; else end = stt;
					if (!remain.empty()) to = RNA2DNA_copy(s(extract_char(remain))); else { info+=":ErrAlt"; return 0; }
					if ((code=='c' || code=='r' || code=='n') && g.strand=='-') { if (!dna_rc(rd,false)) info+=":ErrNT"; if (!dna_rc(to,false)) info+=":ErrNT"; }
					fr = boost::to_upper_copy(DNA_seq(read_chr_num(chr),stt,1));
					if (fr!=rd) { info+=":ErrRef"; return 0; }
					if (to[0]=='R') { if (fr[0]=='A') to[0]='G'; else if (fr[0]=='G') to[0]='A'; }
					if (to[0]=='Y') { if (fr[0]=='C') to[0]='T'; else if (fr[0]=='T') to[0]='C'; }
					if (to[0]=='S') { if (fr[0]=='G') to[0]='C'; else if (fr[0]=='C') to[0]='G'; }
					if (to[0]=='W') { if (fr[0]=='A') to[0]='T'; else if (fr[0]=='T') to[0]='A'; }
					if (to[0]=='K') { if (fr[0]=='G') to[0]='T'; else if (fr[0]=='T') to[0]='G'; }
					if (to[0]=='M') { if (fr[0]=='A') to[0]='C'; else if (fr[0]=='C') to[0]='A'; }
					if (str_has(info,":Err")) return 0;
				}
				else if (str_has(remain,"delins")) // c.112_117delinsTG / BRCA1:c.5564_5572delins8 (not allowed by HGVS standards but I programed it anyway)
				{
					string left = substr_before_find(remain,"delins");
					vector<string> ifs;
					boost::split(ifs,left,boost::is_any_of("_"));
					if (ifs.size()>0) { if (!read_hgvs_position(code,ifs[0],g,chr,stt,info)) return 0; } else { info+=":ErrNo1stPos"; return 0; }
					if (ifs.size()>1) { if (!read_hgvs_position(code,ifs[1],g,chr,end,info)) return 0; } else end = stt;
					if (ifs.size()>2 || str_has(left,"?") || str_has(left,"(")) { info+=":ErrUknPos"; return 0; }
					to = RNA2DNA_copy(substr_after_find(remain,"delins"));
					if (is_integer(to)) to = string(boost::lexical_cast<int>(to),'N');
					if (to.empty()) { info+=":ErrIns"; return 0; }
					if ((code=='c' || code=='r' || code=='n') && g.strand=='-') { std::swap(stt,end); if (!dna_rc(to,false)) info+=":ErrNT"; }
					fr = boost::to_upper_copy(DNA_seq(read_chr_num(chr),stt,(end-stt+1)));
					if (str_has(info,":Err")) return 0;
				}
				else if (str_has(remain,"del") && str_has(remain,"ins")) // c.112_117delAGGTCAinsTG / c.112_117delAGGTCAins2 (not allowed by HGVS standards but I programed it anyway)
				{
					string left = substr_before_find(remain,"del");
					string right = substr_after_find(remain,"del");
					string rd;
					vector<string> ifs;
					boost::split(ifs,left,boost::is_any_of("_"));
					if (ifs.size()>0) { if (!read_hgvs_position(code,ifs[0],g,chr,stt,info)) return 0; } else { info+=":ErrNo1stPos"; return 0; }
					if (ifs.size()>1) { if (!read_hgvs_position(code,ifs[1],g,chr,end,info)) return 0; } else { end = stt; }
					if (ifs.size()>2 || str_has(left,"?") || str_has(left,"(")) { info+=":ErrUknPos"; return 0; }
					rd = RNA2DNA_copy(substr_before_find(right,"ins"));
					to = RNA2DNA_copy(substr_after_find(right,"ins"));
					if (is_integer(rd))
					{
						// check whether 7644-7636+1 == 9 in ATM:c.7636_7644del9
						if (ifs.size()!=1)
						{
							int num_del1 = boost::lexical_cast<int>(rd);
							int num_del2 = abs(end - stt) + 1;
							if (num_del2!=num_del1) info+=":ErrNumNtDel";
						}
						else
						{
							int num_del1 = boost::lexical_cast<int>(rd);
							if ((code=='c' || code=='r' || code=='n') && g.strand=='-')	end=stt-(num_del1-1);
							else														end=stt+(num_del1-1);
						}
						// rd here is duplicated information, just clear it. So ATM:c.7636_7644del9 become ATM:c.7636_7644del.
						rd.clear();
					}
					if (rd.size()>1 && ifs.size()==1)
					{
						if ((code=='c' || code=='r' || code=='n') && g.strand=='-')	end=stt-(rd.size()-1);
						else														end=stt+(rd.size()-1);
					}
					if (is_integer(to)) to = string(boost::lexical_cast<int>(to),'N');
					if (to.empty()) { info+=":ErrIns"; return 0; }
					if ((code=='c' || code=='r' || code=='n') && g.strand=='-') { std::swap(stt,end); if (!dna_rc(rd,false)) info+=":ErrNT"; if (!dna_rc(to,false)) info+=":ErrNT"; }
					fr = boost::to_upper_copy(DNA_seq(read_chr_num(chr),stt,(end-stt+1)));
					if (!rd.empty() && fr!=rd) { info+=":ErrRef"; return 0; }
					if (str_has(info,":Err")) return 0;
				}
				else if (str_has(remain,"inv")) // c.203_506inv
				{
					string left = substr_before_find(remain,"inv");
					vector<string> ifs;
					boost::split(ifs,left,boost::is_any_of("_"));
					if (ifs.size()>0) { if (!read_hgvs_position(code,ifs[0],g,chr,stt,info)) return 0; } else { info+=":ErrNo1stPos"; return 0; }
					if (ifs.size()>1) { if (!read_hgvs_position(code,ifs[1],g,chr,end,info)) return 0; } else { info+=":ErrNo2ndPos"; return 0; }
					if (ifs.size()>2 || str_has(left,"?") || str_has(left,"(")) { info+=":ErrUknPos"; return 0; }
					if ((code=='c' || code=='r' || code=='n') && g.strand=='-') { std::swap(stt,end); }
					fr = boost::to_upper_copy(DNA_seq(read_chr_num(chr),stt,(end-stt+1)));
					to = fr; std::reverse(to.begin(),to.end());
					if (str_has(info,":Err")) return 0;
				}
				else if (str_has(remain,"del")) // BRCA1:c.185_186delCT / BRCA2:c.517-23delTA / ATM:c.7636_7644del9 / ATM:c.7636_7644del / ATM:c.7636del9.
				{
					string left = substr_before_find(remain,"del");
					vector<string> ifs;
					boost::split(ifs,left,boost::is_any_of("_"));
					if (ifs.size()>0) { if (!read_hgvs_position(code,ifs[0],g,chr,stt,info)) return 0; } else { info+=":ErrNo1stPos"; return 0; }
					if (ifs.size()>1) { if (!read_hgvs_position(code,ifs[1],g,chr,end,info)) return 0; } else { end = stt; }
					if (ifs.size()>2 || str_has(left,"?") || str_has(left,"(")) { info+=":ErrUknPos"; return 0; }
					string rd = RNA2DNA_copy(substr_after_find(remain,"del"));
					if (is_integer(rd))
					{
						// check whether 7644-7636+1 == 9 in ATM:c.7636_7644del9
						if (ifs.size()!=1)
						{
							int num_del1 = boost::lexical_cast<int>(rd);
							int num_del2 = abs(end - stt) + 1;
							if (num_del2!=num_del1) info+=":ErrNumNtDel";
						}
						else
						{
							int num_del1 = boost::lexical_cast<int>(rd);
							if ((code=='c' || code=='r' || code=='n') && g.strand=='-')	end=stt-(num_del1-1);
							else														end=stt+(num_del1-1);
						}
						// rd here is duplicated information, just clear it. So ATM:c.7636_7644del9 become ATM:c.7636_7644del.
						rd.clear();
					}
					if (rd.size()>1 && ifs.size()==1)
					{
						if ((code=='c' || code=='r' || code=='n') && g.strand=='-')	end=stt-(rd.size()-1);
						else														end=stt+(rd.size()-1);
					}
					if ((code=='c' || code=='r' || code=='n') && g.strand=='-') { std::swap(stt,end); if (!dna_rc(rd,false)) info+=":ErrNT"; }
					fr = boost::to_upper_copy(DNA_seq(read_chr_num(chr),stt,(end-stt+1)));
					if (!rd.empty() && fr!=rd) info+=":ErrRef";
					to = boost::to_upper_copy(DNA_seq(read_chr_num(chr),--stt,1));
					fr = to+fr;
					if (str_has(info,":Err")) return 0;
				}
				else if (str_has(remain,"dup")) // NM_024675:c.1947dupA
				{
					string left = substr_before_find(remain,"dup");
					vector<string> ifs;
					boost::split(ifs,left,boost::is_any_of("_"));
					if (ifs.size()>0) { if (!read_hgvs_position(code,ifs[0],g,chr,stt,info)) return 0; } else { info+=":ErrNo1stPos"; return 0; }
					if (ifs.size()>1) { if (!read_hgvs_position(code,ifs[1],g,chr,end,info)) return 0; } else end = stt;
					if (ifs.size()>2 || str_has(left,"?") || str_has(left,"(")) { info+=":ErrUknPos"; return 0; }
					string rd = RNA2DNA_copy(substr_after_find(remain,"dup"));
					if ((code=='c' || code=='r' || code=='n') && g.strand=='-') { std::swap(stt,end); if (!dna_rc(rd,false)) info+=":ErrNT"; }
					fr = boost::to_upper_copy(DNA_seq(read_chr_num(chr),stt,(end-stt+1)));
					to = fr + fr;
					if (!rd.empty() && fr!=rd) { info+=":ErrRef"; return 0; }
					if (str_has(info,":Err")) return 0;
				}
				else if (str_has(remain,"["))
				{
					string left = substr_before_find(remain,"[");
					vector<string> ifs;
					boost::split(ifs,left,boost::is_any_of("_"));
					if (ifs.size()>0) { if (!read_hgvs_position(code,ifs[0],g,chr,stt,info)) return 0; } else { info+=":ErrNo1stPos"; return 0; }
					if (ifs.size()>1) { if (!read_hgvs_position(code,ifs[1],g,chr,end,info)) return 0; } else end = stt;
					if (ifs.size()>2 || str_has(left,"?") || str_has(left,"(")) { info+=":ErrUknPos"; return 0; }
					if (ifs.size()>1)	remain = ifs[1] + trim_before_find(remain,"[");
					else				remain = ifs[0] + trim_before_find(remain,"[");
					if (remain[0]=='[') // NM_000051:c.*8_*9[4]
					{
						if ((code=='c' || code=='r' || code=='n') && g.strand=='-') { std::swap(stt,end); }
						fr = boost::to_upper_copy(DNA_seq(read_chr_num(chr),stt,(end-stt+1)));
						remain.erase(0,1);
						int dup;
						try { dup = extract_int(remain) / fr.size(); } catch (bad_extracting &e) { info+=":ErrDupNum(NaN)"; return 0; }
						to=fr; for (int i=0;i<dup;++i) to+=fr;
						if (remain!="]") { info+=":Err2+Rpt"; return 0; } // does not support multiple repeated elements g.456TG[4]TA[9]TG[3] (or g.456_465[4]466_489[9]490_499[3])
					}
					else // NM_000051:c.*8TA[4]
					{
						string rd = RNA2DNA_copy(substr_before_find(remain,"["));
						end = stt + rd.size() -1;
						if ((code=='c' || code=='r' || code=='n') && g.strand=='-') { std::swap(stt,end); if (!dna_rc(rd,false)) info+=":ErrNT"; }
						fr = boost::to_upper_copy(DNA_seq(read_chr_num(chr),stt,(end-stt+1)));
						if (fr!=rd) info+=":ErrRef";
						remain=substr_after_find(remain,"[");;
						int dup ;
						try { dup = extract_int(remain) / fr.size(); } catch (bad_extracting &e) { info+=":ErrDupNum(NaN)"; return 0; }
						to=fr; for (int i=0;i<dup;++i) to+=fr;
						if (remain!="]") info+=":Err2+Rpt"; // does not support multiple repeated elements g.456TG[4]TA[9]TG[3] (or g.456_465[4]466_489[9]490_499[3])
						if (str_has(info,":Err")) return 0;
					}
					if (str_has(info,":Err")) return 0;
				}
				else if (str_has(remain,"ins")) // NM_000051:c.7176_7177insT / BRCA1:c.80_81ins6 ("not allowed" by HGVS but I program it anyway) / ATM:c.5713insT (not standard but I allow it)
				{
					string left = substr_before_find(remain,"ins");
					vector<string> ifs;
					boost::split(ifs,left,boost::is_any_of("_"));
					if (ifs.size()>0) { if (!read_hgvs_position(code,ifs[0],g,chr,stt,info)) return 0; } else { info+=":ErrNo1stPos"; return 0; }
					if (ifs.size()>1) { if (!read_hgvs_position(code,ifs[1],g,chr,end,info)) return 0; } else { end=0; } // else { info+=":ErrNo2ndPos"; return 0; }
					if (ifs.size()>2 || str_has(left,"?") || str_has(left,"(")) { info+=":ErrUknPos"; return 0; }
					to = RNA2DNA_copy(substr_after_find(remain,"ins"));
					if (is_integer(to)) to = string(boost::lexical_cast<int>(to),'N');
					if (to.empty()) { info+=":ErrIns"; return 0; }
					if (end==0) { if ((code=='c' || code=='r' || code=='n') && g.strand=='-') end=stt-1; else end=stt+1; }
					if ((code=='c' || code=='r' || code=='n') && g.strand=='-') { std::swap(stt,end); if (!dna_rc(to,false)) info+=":ErrNT"; }
					if (end!=stt+1) { info+=":ErrEnd!=Start+1"; return 0; }
					fr = boost::to_upper_copy(DNA_seq(read_chr_num(chr),stt,1));
					to = fr+to;
					end = stt;
					if (str_has(info,":Err")) return 0;
				}
				else if (str_has(remain,"con")) // g.123_678conNG_012232.1:g.9456_10011
				{
					info+=":ConvNotProgramed"; return 0; // this version doesn't support conversions
				}
				else // NM_000051:c.8158G>C / NM_000051:c.-5G>C / NM_000051:c.4909+1G>A / NM_000051:c.8851-1G>T / NM_000051:c.*5C>T / NM_000051:c.IVS1-5T>C / NM_000051:c.IVS2+1G>C
				{
					string ch, rd;
					if (!read_hgvs_position(code,remain,g,chr,stt,info)) return 0; else end = stt;
					if (!remain.empty()) rd = RNA2DNA_copy(s(extract_char(remain))); else { info+=":ErrIncompleteHGVS"; return 0; }
					if (!remain.empty()) ch = s(extract_char(remain)); else { info+=":Err'>'"; return 0; }
					if (!remain.empty()) to = RNA2DNA_copy(s(extract_char(remain))); else { info+=":ErrAlt"; return 0; }
					if ((code=='c' || code=='r' || code=='n') && g.strand=='-') { if (!dna_rc(rd,false)) info+=":ErrNT"; if (!dna_rc(to,false)) info+=":ErrNT"; }
					fr = boost::to_upper_copy(DNA_seq(read_chr_num(chr),stt,1));
					if (fr!=rd) { info+=":ErrRef"; return 0; }
					if (to[0]=='R') { if (fr[0]=='A') to[0]='G'; else if (fr[0]=='G') to[0]='A'; }
					if (to[0]=='Y') { if (fr[0]=='C') to[0]='T'; else if (fr[0]=='T') to[0]='C'; }
					if (to[0]=='S') { if (fr[0]=='G') to[0]='C'; else if (fr[0]=='C') to[0]='G'; }
					if (to[0]=='W') { if (fr[0]=='A') to[0]='T'; else if (fr[0]=='T') to[0]='A'; }
					if (to[0]=='K') { if (fr[0]=='G') to[0]='T'; else if (fr[0]=='T') to[0]='G'; }
					if (to[0]=='M') { if (fr[0]=='A') to[0]='C'; else if (fr[0]=='C') to[0]='A'; }
					if (str_has(info,":Err")) return 0;
				}
				break; }
			case 'p':
				chr=g.chr;
				if (remain=="?" || remain=="(=)" || remain=="=" || remain=="0" || remain=="0?")
				{
					fr=remain;
					to=remain;
					break;
				}
				if (remain.size()<3) { info+=":ErrIncompleteHGVS"; return 0; }
				if (remain[0]=='(' && remain[remain.size()-1]==')') remain=remain.substr(1,remain.size()-2); // p.(Arg3381Ser) => p.Arg3381Ser
				if (remain.size()<3) { info+=":ErrIncompleteHGVS"; return 0; }
				if (remain[1]>='0' && remain[1]<='9') // p.W26* / p.*110Glnext*17 / p.R83=
				{
					try {
						fr = s(extract_char(remain));
						stt = end = extract_int(remain);
						to = remain;
					} catch (bad_extracting &e) { info+=":ErrPos"; return 0; }
					if (fr=="*" && str_has(to,"ext")) to=substr_before_find(to,"ext");
				}
				else // p.Trp26Ter / p.Ter110GlnextTer17 / p.Arg83=
				{
					if (remain.find('+')!=string::npos || remain.find('-')!=string::npos) { info+=":p.xxxWith+-NotProgramed"; return 0; } // need coding
					try {
						fr = extract_alphabets(remain);
						stt = end = extract_int(remain);
						to = remain;
					} catch (bad_extracting &e) { info+=":ErrPos"; return 0; }
					if (fr=="Ter" && str_has(to,"ext")) to=substr_before_find(to,"ext");
				}
				if (exist_element(AA_3to1,fr)) fr = AA_3to1[fr]; else { info+=":ErrAAfrom"; return 0; }
				if (exist_element(AA_3to1,to)) to = AA_3to1[to]; else { if (to=="=") to=fr; else { info+=":ErrAAto"; return 0; } }
				break;
			default:
				info+=":ErrLvl("+s(code)+")"; // n.
				return 0;
		}
		return 1;
	}

	void _gene_cds_sizes(gene_info& g, vector<int>& out)
	{
		out.clear();
		for (int i=0; i<g.exonCount; i++)
		{
			int s = g.exonStarts[i];
			int e = g.exonEnds[i];
			if (g.cdsStart > s)	s = g.cdsStart;
			if (g.cdsEnd < e)	e = g.cdsEnd;
			if (s>=e) continue; // not s>e
			out.push_back(e-s); // not e-s+1, assuming all start coordinates are 0-based (UCSC).
		}
	}
	
	int _gene_CDS_len(gene_info& g)
	{
		vector<int> cds_sizes;
		_gene_cds_sizes(g,cds_sizes);
		int len=0;
		for (each_element(cds_sizes, it)) len += *it;
		return len;
	}

	int _gene_UTR5_len(gene_info& g)
	{
		if (g.cdsStart==g.cdsEnd) return g.txLen;
		int len=0;
		if (g.strand=='+')
		{
			for (int i=0; i<g.exonCount; ++i)
			{
				if (g.cdsStart < g.exonEnds[i])
				{
					len += g.cdsStart - g.exonStarts[i];
					break;
				}
				else
				{
					len += g.exonEnds[i] - g.exonStarts[i];
				}
			}
		}
		else
		{
			for (int i=g.exonCount-1; i>=0; --i)
			{
				if (g.cdsEnd > g.exonStarts[i])
				{
					len += g.exonEnds[i] - g.cdsEnd;
					break;
				}
				else
				{
					len += g.exonEnds[i] - g.exonStarts[i];
				}
			}
		}
		return len;
	}
	
	void _offset_cds(gene_info& g)
	{
		if (g.cdsStart==g.cdsEnd || g.cdsOffset==0) return ;
		if (g.strand=='+')
		{
			for (int i=0; i<g.exonCount; ++i)
			{
				if (g.cdsStart < g.exonEnds[i])
				{
					if (g.cdsStart+g.cdsOffset<g.exonEnds[i]) g.cdsStart += g.cdsOffset;
					else if (i+1<g.exonCount) g.cdsStart = g.exonStarts[i+1]+(g.cdsOffset-(g.exonEnds[i]-g.cdsStart));
					else exit_error("cannot offset cds for "+g.name);
					break;
				}
			}
		}
		else
		{
			for (int i=g.exonCount-1; i>=0; --i)
			{
				if (g.cdsEnd > g.exonStarts[i])
				{
					if (g.cdsEnd-g.cdsOffset>g.exonStarts[i]) g.cdsEnd -= g.cdsOffset;
					else if (i-1>=0) g.cdsEnd = g.exonEnds[i-1]-(g.cdsOffset-(g.cdsEnd-g.exonStarts[i]));
					else exit_error("cannot offset cds for "+g.name);
					break;
				}
			}
		}
	}
	
	void _restore_cds(gene_info& g)
	{
		if (g.cdsStart==g.cdsEnd || g.cdsOffset==0) return ;
		if (g.strand=='+')
		{
			for (int i=0; i<g.exonCount; ++i)
			{
				if (g.cdsStart < g.exonEnds[i])
				{
					if (g.cdsStart-g.cdsOffset>=g.exonStarts[i]) g.cdsStart -= g.cdsOffset;
					else if (i-1>=0) g.cdsStart = g.exonEnds[i-1]-1-(g.cdsOffset-(g.cdsStart-g.exonStarts[i]+1));
					else exit_error("cannot restore cds for "+g.name);
					break;
				}
			}
		}
		else
		{
			for (int i=g.exonCount-1; i>=0; --i)
			{
				if (g.cdsEnd > g.exonStarts[i])
				{
					if (g.cdsEnd+g.cdsOffset<=g.exonEnds[i]) g.cdsEnd += g.cdsOffset;
					else if (i+1<g.exonCount) g.cdsEnd = g.exonStarts[i+1]+1+(g.cdsOffset-(g.exonEnds[i]-g.cdsEnd+1));
					else exit_error("cannot restore cds for "+g.name);
					break;
				}
			}
		}
	}

	// read file wget http://hgdownload.cse.ucsc.edu/goldenPath/h g 1 9/database/refGene.txt.gz
	// all start positions are 0-based while end positions are 1-based
	// remove records whose chromosome is chr#_xxx, which are haplotypes not assigned to ref genome, see
	// http://vega.sanger.ac.uk/info/data/MHC_Homo_sapiens.html
	// http://genomeref.blogspot.com/
	// http://genome.ucsc.edu/cgi-bin/hgTracks?hgsid=191600823&chromInfoPage=
	// cdsStartStat/cdsEndStat is described in
	// https://lists.soe.ucsc.edu/pipermail/genome/2005-December/009184.html
	// On Apr 17 2012, there're 42217 data lines in the file, where the break down of cdsStartStat is:
	// cmpl   	35189	83.3456
	// incmpl 	56		0.139754
	// unk    	6972	16.5147
	// cdsEndStat is similar (35186/59/6972, respectively).
	// So I use this to choose transcript: cmpl >  unk > incmpl, and cdsStartStat > cdsEndStat.
	// 2nd criteria is longest cds
	// 3rd criteria is longest cDNA
	// as per http://www.biostars.org/post/show/19029/what-is-rationale-behind-ucsc-canonical-transcripts/
	// input format of refGene.txt
	// ---------------------------------------------------------------------------------------------------------------------
	// field		example		SQL type							info	description
	// ---------------------------------------------------------------------------------------------------------------------
	// bin			112			smallint(5) unsigned				range	Indexing field to speed chromosome range queries
	// name			NM_005300	varchar(255)						values	Name of gene (usually transcript_id from GTF)
	// chrom		chrX		varchar(255)						values	Reference sequence chromosome or scaffold
	// strand		+			char(1)								values	+ or - for strand
	// txStart		41548225	int(10) unsigned					range	Transcription start position
	// txEnd		41556530	int(10) unsigned					range	Transcription end position
	// cdsStart		41554886	int(10) unsigned					range	Coding region start
	// cdsEnd		41556032	int(10) unsigned					range	Coding region end
	// exonCount	3			int(10) unsigned					range	Number of exons
	// exonStarts	41548225,41548998,41554880,	longblob					Exon start positions
	// exonEnds		41548337,41549088,41556530,	longblob					Exon end positions
	// score		0			int(11)								range	score
	// name2		GPR34		varchar(255)						values	Alternate name (e.g. gene_id from GTF)
	// cdsStartStat	cmpl		enum('none','unk','incmpl','cmpl')	values	enum('none','unk','incmpl','cmpl')
	// cdsEndStat	cmpl		enum('none','unk','incmpl','cmpl')	values	enum('none','unk','incmpl','cmpl')
	// exonFrames	-1,-1,0,	longblob									Exon frame {0,1,2}, or -1 if no frame for exon
	// ---------------------------------------------------------------------------------------------------------------------
	
	// input format of ensGene.txt
	// ---------------------------------------------------------------------------------------------------------------------
	// field		example				SQL type							info	description
	// ---------------------------------------------------------------------------------------------------------------------
	// bin			585					smallint(5) unsigned				range	Indexing field to speed chromosome range queries
	// name			ENST00000456328		varchar(255)						values	Name of gene (usually transcript_id from GTF)
	// chrom		chr1				varchar(255)						values	Reference sequence chromosome or scaffold
	// strand		+					char(1)								values	+ or - for strand
	// txStart		11868				int(10) unsigned					range	Transcription start position
	// txEnd		14409				int(10) unsigned					range	Transcription end position
	// cdsStart		14409				int(10) unsigned					range	Coding region start
	// cdsEnd		14409				int(10) unsigned					range	Coding region end
	// exonCount	3					int(10) unsigned					range	Number of exons
	// exonStarts	11868,12612,13220,	longblob									Exon start positions
	// exonEnds		12227,12721,14409,	longblob									Exon end positions
	// score		0					int(11)								range	score
	// name2		ENSG00000223972		varchar(255)						values	Alternate name (e.g. gene_id from GTF)
	// cdsStartStat	none				enum('none','unk','incmpl','cmpl')	values	enum('none','unk','incmpl','cmpl')
	// cdsEndStat	none				enum('none','unk','incmpl','cmpl')	values	enum('none','unk','incmpl','cmpl')
	// exonFrames	-1,-1,-1,			longblob									Exon frame {0,1,2}, or -1 if no frame for exon
	// ---------------------------------------------------------------------------------------------------------------------
	
	// filename can also be ensGene.txt because it has the same format
	void read_refGene(const string& filename, bool read_tx, bool filter)
	{
		gene_bySymbol.clear();
		gene_byTranscript.clear();
		gene_db.clear();
		
		map< string, char >			geneIDtoStrand;	// geneID -> strand
		map< string, string >		refseq2geneID;	// RefSeq -> geneID
		map< string, string >		NatReadthrough;	// geneID -> geneSymbol
		map< string, string >		pseudogene_gID;	// geneID -> geneSymbol
		map< string,set<string> >	MostFrqOccurTx;	// geneID -> RefSeq
		map< string,set<string> >	SProtCanonical;	// geneID -> RefSeq

		if (filter)
		{
			if (!path_set) exit_error("genepi path not set");
			
			set<string> to_be_read; // RefSeq ID
			for (Rows_in_File(row,filename,16))
			{
				if (!read_chr_num(row[2])) continue;
				to_be_read.insert(row[1]);
			}
			
			//string strand_file="$BING/work/data/h g 1 9/_CCDS_GeneID_Strand";
			//lns<<showl<<"Reading "<<strand_file<<flush_logger;
			//for (Rows_in_File(in,strand_file,2))
			//	geneIDtoStrand[in[0]]=in[1][0];
			
			if (!txID_gID.empty())
			{
				lns<<showl<<"Reading "<<txID_gID<<flush_logger;
				for (Rows_in_File(in,txID_gID,2))
					refseq2geneID[in[0]]=in[1];
			}
			
			if (!read_thr.empty())
			{
				lns<<showl<<"Reading "<<read_thr<<flush_logger;
				for (Rows_in_File(in,read_thr,2))
					NatReadthrough[in[0]]=in[1];
			}
			
			if (!pseudoid.empty())
			{
				lns<<showl<<"Reading "<<pseudoid<<flush_logger;
				for (Rows_in_File(in,pseudoid,2))
					pseudogene_gID[in[0]]=in[1];
			}
			
			if (!gID_MFOT.empty())
			{
				lns<<showl<<"Reading "<<gID_MFOT<<flush_logger;
				for (Rows_in_File(in,gID_MFOT,2))
					if (exist_element(to_be_read,in[1]))
						MostFrqOccurTx[in[0]].insert(in[1]);
			}
			
			if (!uGID_can.empty())
			{
				lns<<showl<<"Reading "<<uGID_can<<flush_logger;
				for (Rows_in_File(in,uGID_can,2))
				{
					string& txID = in[1];
					string& gID = in[0];
					if (!exist_element(refseq2geneID,txID)) continue;
					if (refseq2geneID[txID]!=gID) continue; // above added on 2018-10-09. Appris has 9718 NM_014693, but NM_014693 is a readthrough in the new NCBI with GeneID 110599583.
					if (exist_element(to_be_read,in[1]) && !exist_element(MostFrqOccurTx,in[0])) SProtCanonical[in[0]].insert(in[1]);
				}
			}
		}

		if (filter) lns<<showl<<"Reading "<<filename<<" ... "<<flush_logger;
		for (Rows_in_File(row,filename,16))
		{
			if (!read_chr_num(row[2])) continue;
			if (filter)
			{
				// Remove definitely unwanted Tx (not the most frequently observed or is a readthrough).
				// "on the wrong strand" is not used because they may be gene duplications. To bring it back, just un-comment the first line below.
				if (exist_element(refseq2geneID,row[1]))
				{
					string & GeneID = refseq2geneID[row[1]];
				//	if (exist_element(geneIDtoStrand,GeneID) && geneIDtoStrand[GeneID]!=row[3][0]) continue;
					if (exist_element(MostFrqOccurTx,GeneID) && !exist_element(MostFrqOccurTx[GeneID],row[1])) continue;
					if (exist_element(SProtCanonical,GeneID) && !exist_element(SProtCanonical[GeneID],row[1])) continue;
					if (exist_element(NatReadthrough,GeneID)) continue;
					if (exist_element(pseudogene_gID,GeneID)) continue;
				}
				else
				{
					// "continue" will remove the RefSeq that are not in NCBI, which are discontinued by NCBI but remained in UCSC. It's OK to remove them.
					// Otherwise, they will be kept in the output. But it's likely that I cannot translate the RefSeq ID to an official Symbol or GeneID.
					continue;
				}
			}
			
			// read data
			gene_info newg;
			if (!read_chr_str(row[2],newg.chr)) exit_error("cannot read "+row[2]+" as a chromosome");
			newg.chrNumPlink = read_chr_num(newg.chr); if (!newg.chrNumPlink) exit_error("cannot read "+newg.chr+" as a chromosome");
			newg.strand = row[3][0]; if (newg.strand!='+'&&newg.strand!='-') exit_error("strand should be '+' or '-'.");
			newg.txStart = boost::lexical_cast<int>(row[4]);
			newg.txEnd = boost::lexical_cast<int>(row[5]);
			newg.cdsStart = boost::lexical_cast<int>(row[6]);
			newg.cdsEnd = boost::lexical_cast<int>(row[7]);
			newg.exonCount = boost::lexical_cast<int>(row[8]);
			newg.txLen = 0;
			newg.cdsOffset = 0;
			newg.col16 = row[15];
			try
			{
				for (int i=0; i<newg.exonCount; i++)
				{
					newg.exonStarts.push_back(extract_int(row[9]));
					newg.exonEnds.push_back(extract_int(row[10]));
					extract_char(row[9]);
					extract_char(row[10]);
					newg.txLen += newg.exonEnds.back() - newg.exonStarts.back();
				}
			} catch (bad_extracting &e) { exit_error("failed to read exon start/end positions for "+row[1]); }
			try
			{
				vector<string> INFO;
				boost::split(INFO,row[15],boost::is_any_of(";"));
				if (get_int(INFO,"cdsOffset",newg.cdsOffset))
				{
					_offset_cds(newg);
				}
				else
				{
					vector<int> exonFrames;
					for (auto &i:INFO)
						if (!str_has(i,"=") && str_endsw(i,","))
						{
							while (!i.empty())
							{
								int f=extract_int(i); extract_char(i);
								exonFrames.push_back(f);
							}
							break;
						}
					if ((int)exonFrames.size()==newg.exonCount)
					{
						if (newg.strand=='+')
						{
							for (int i=0; i<newg.exonCount; i++)
							{
								int& f = exonFrames[i];
								if 		(f==-1) continue;
								else if (f==0) newg.cdsOffset=0;
								else if (f==1) newg.cdsOffset=2;
								else if (f==2) newg.cdsOffset=1;
								else exit_error("unknown exonFrames "+itos(f));
								break;
							}
						}
						else
						{
							for (int i=0; i<newg.exonCount; i++)
							{
								int& f = exonFrames[newg.exonCount-1-i];
								if 		(f==-1) continue;
								else if (f==0) newg.cdsOffset=0;
								else if (f==1) newg.cdsOffset=2;
								else if (f==2) newg.cdsOffset=1;
								else exit_error("unknown exonFrames "+itos(f));
								break;
							}
						}
						_offset_cds(newg);
					}
				}
			} catch (bad_extracting &e) { exit_error("failed to read exon frame for "+row[1]); }
			newg.cdsLen = _gene_CDS_len(newg);
			newg.u5Len = _gene_UTR5_len(newg);
			newg.u3Len = newg.txLen-newg.u5Len-newg.cdsLen;
			newg.txRgn = newg.txEnd - newg.txStart;
			newg.cdsRgn = newg.cdsEnd - newg.cdsStart;
			newg.name = row[1];
			if (filter)
			{
				newg.name2=Rsq2Sym(row[1],false);
				if (newg.name2.empty()) newg.name2=row[12];
			}
			else
			{
				newg.name2=row[12];
				if (row[12]!=row[0] && !is_integer(row[0]) && row[0]!=".") newg.name2+="("+row[0]+")";
			}
			if (newg.strand=='-')	newg.name3 = newg.name + "_CHR" + newg.chr + "_" + itos(newg.txEnd);
			else					newg.name3 = newg.name + "_CHR" + newg.chr + "_" + itos(newg.txStart+1);
			if (newg.strand=='-')	newg.name4 = newg.name2+ "_CHR" + newg.chr + "_" + itos(newg.txEnd);
			else					newg.name4 = newg.name2+ "_CHR" + newg.chr + "_" + itos(newg.txStart+1);
			newg.symbol=row[0];
			if (!incl_g.empty() && !exist_element(incl_g,newg.name2)) continue;
			if (!incl_t.empty() && !exist_element(incl_t,newg.name )) continue;

			newg.cdsStartStat=row[13];
			newg.cdsEndStat=row[14];
			int cdsStartStat=0, cdsEndStat=0; // "none" is 0 automatically
			if		(row[13]=="incmpl")	cdsStartStat=1;
			else if (row[13]=="unk")	cdsStartStat=2;
			else if (row[13]=="cmpl")	cdsStartStat=3;
			if		(row[14]=="incmpl")	cdsEndStat=1;
			else if (row[14]=="unk")	cdsEndStat=2;
			else if (row[14]=="cmpl")	cdsEndStat=3;
			if (newg.strand=='+')	newg.cdsStat = cdsStartStat*10 + cdsEndStat;
			else					newg.cdsStat = cdsEndStat*10 + cdsStartStat;
			newg.lenStat = std::pair<int,int>(newg.cdsLen,newg.txLen);
			gene_cal_hgvs(newg);
			// select transcript by cdsStat & lenStat
			if (exist_element(gene_bySymbol,newg.name2))
			{
				gene_info& oldg = gene_bySymbol[newg.name2].begin()->second;
				if (newg.cdsStat > oldg.cdsStat)
				{
					gene_bySymbol[newg.name2].clear();
					gene_bySymbol[newg.name2].insert(pair<lenStat_t,gene_info>(newg.lenStat,newg));
				}
				else if (newg.cdsStat == oldg.cdsStat)
				{
					if (read_tx)
						gene_bySymbol[newg.name2].insert(pair<lenStat_t,gene_info>(newg.lenStat,newg));
					else if (newg.lenStat == oldg.lenStat)
						gene_bySymbol[newg.name2].insert(pair<lenStat_t,gene_info>(newg.lenStat,newg));
					else if (newg.lenStat  > oldg.lenStat)
					{
						gene_bySymbol[newg.name2].clear();
						gene_bySymbol[newg.name2].insert(pair<lenStat_t,gene_info>(newg.lenStat,newg));
					}
					else continue;
				}
				else continue;
			}
			else
				gene_bySymbol[newg.name2].insert(pair<lenStat_t,gene_info>(newg.lenStat,newg));
		}
		for (auto &gene:gene_bySymbol) for (auto &tx:gene.second) gene_byTranscript[tx.second.name].insert(pair<lenStat_t,gene_info>(tx.second.lenStat,tx.second));
		for (auto &gene:gene_bySymbol) for (auto &tx:gene.second) gene_db[tx.second.name]=tx.second;
		if (filter) lns<<showl<<gene_bySymbol.size()<<" genes read."<<flush_logger;
	}
	
	// input format of hg19 knownGene.txt (knownGene discontinued for GRCh38):
	// ---------------------------------------------------------------------------------------------------------------------------------
	// field		example				SQL type			info	description
	// ---------------------------------------------------------------------------------------------------------------------------------
	// name			uc001aaa.3			varchar(255)		values	Name of gene
	// chrom		chr1				varchar(255)		values	Reference sequence chromosome or scaffold
	// strand		+					char(1)				values	+ or - for strand
	// txStart		11873				int(10) unsigned	range	Transcription start position (or end position for minus strand item)
	// txEnd		14409				int(10) unsigned	range	Transcription end position (or start position for minus strand item)
	// cdsStart		11873				int(10) unsigned	range	Coding region start (or end position if for minus strand item)
	// cdsEnd		11873				int(10) unsigned	range	Coding region end (or start position if for minus strand item)
	// exonCount	3					int(10) unsigned	range	Number of exons
	// exonStarts	11873,12612,13220,	longblob	 				Exon start positions (or end positions for minus strand item)
	// exonEnds		12227,12721,14409,	longblob	 				Exon end positions (or start positions for minus strand item)
	// proteinID						varchar(40)			values	UniProt display ID, UniProt accession, or RefSeq protein ID
	// alignID		uc001aaa.3			varchar(255)		values	Unique identifier (GENCODE transcript ID for GENCODE Basic)
	// ---------------------------------------------------------------------------------------------------------------------------------
	void read_knownGene(string prefix, bool read_tx)
	{
		gene_bySymbol.clear();
		gene_byTranscript.clear();
		gene_db.clear();
		
		if (!str_endsw(prefix,"/")) prefix+='/';
		tfile_format fmt(default_tfile_format);
		fmt.set_field_nums(5,"lines missing required field",tfile_format::Continue);
		fmt.set_delimiters("\t");
		lns<<showl<<"Reading "<<prefix<<"kgXref.txt"<<flush_logger;
		map<string,string> knownToGeneSymbol;
		for (Rows_in_File(in,prefix+"kgXref.txt",&fmt))
		{	knownToGeneSymbol[in[0]]=in[4]; }
		
		lns<<showl<<"Reading "<<prefix<<"knownCanonical.txt"<<flush_logger;
		set<string> knownCanonical, knownGeneSymbol;
		for (Rows_in_File(in,prefix+"knownCanonical.txt",6))
		{	knownCanonical.insert(in[4]);
			knownGeneSymbol.insert(knownToGeneSymbol[in[4]]); }
		
		lns<<showl<<"Reading "<<prefix<<"knownGene.txt"<<" ... "<<flush_logger;
		for (Rows_in_File(row,prefix+"knownGene.txt",12))
		{
			// decide whether to read
			string geneSymbol = knownToGeneSymbol[row[0]];
			if (!exist_element(knownCanonical,row[0]))
				if (geneSymbol.empty() || exist_element(knownGeneSymbol,geneSymbol))
					continue;
			// knownCanon	knownSymbol	SymbolEmpty	toRead	Notes
			// y			y			n			yes		knownCanon must has knownSymbol
			// n			y			n			no		knownSymbol doesn't contain [blank]
			// n			n			y			no
			// n			n			n			yes
			// The second "if" is important, otherwise some gene cannot be found, eg CYP3A5, DDX11L1, etc.
			// This is because some geneSymbol doesn't associate with a knownCanonical. As of 2013-03-18,
			// I got 27993 knownGene_Canonical with that "if", but 26601 without it.
			
			// read data
			gene_info newg;
			if (!read_chr_str(row[1],newg.chr)) exit_error("cannot read "+row[2]+" as a chromosome");
			newg.chrNumPlink = read_chr_num(newg.chr); if (!newg.chrNumPlink) exit_error("cannot read "+newg.chr+" as a chromosome");
			newg.strand = row[2][0]; if (newg.strand!='+'&&newg.strand!='-') exit_error("strand should be '+' or '-'.");
			newg.txStart = boost::lexical_cast<int>(row[3]);
			newg.txEnd = boost::lexical_cast<int>(row[4]);
			newg.cdsStart = boost::lexical_cast<int>(row[5]);
			newg.cdsEnd = boost::lexical_cast<int>(row[6]);
			newg.exonCount = boost::lexical_cast<int>(row[7]);
			newg.txLen = 0;
			newg.cdsOffset = 0;
			newg.col16 = ".";
			try
			{
				for (int i=0; i<newg.exonCount; i++)
				{
					newg.exonStarts.push_back(extract_int(row[8]));
					newg.exonEnds.push_back(extract_int(row[9]));
					extract_char(row[8]);
					extract_char(row[9]);
					newg.txLen += newg.exonEnds.back() - newg.exonStarts.back();
				}
			} catch (bad_extracting &e) { exit_error("failed to read exon start/end positions for "+row[0]); }
			newg.proteinID = row[10];
			newg.alignID = row[11];
			newg.cdsLen = _gene_CDS_len(newg);
			newg.u5Len = _gene_UTR5_len(newg);
			newg.u3Len = newg.txLen-newg.u5Len-newg.cdsLen;
			newg.name  = row[0]; // not the same as RefSeq! but ID must be stored
			newg.name2 = geneSymbol;
			if (!incl_g.empty() && !exist_element(incl_g,newg.name2)) continue;
			if (!incl_t.empty() && !exist_element(incl_t,newg.name )) continue;

			if (newg.strand=='-')	newg.name3 = newg.name + "_CHR" + newg.chr + "_" + itos(newg.txEnd);
			else					newg.name3 = newg.name + "_CHR" + newg.chr + "_" + itos(newg.txStart+1);
			if (newg.strand=='-')	newg.name4 = newg.name2+ "_CHR" + newg.chr + "_" + itos(newg.txEnd);
			else					newg.name4 = newg.name2+ "_CHR" + newg.chr + "_" + itos(newg.txStart+1);
			newg.cdsStartStat = "NA";
			newg.cdsEndStat = "NA";
			newg.cdsStat = 0;
			newg.lenStat = std::pair<int,int>(newg.cdsLen,newg.txLen);
			gene_cal_hgvs(newg);
			
			// Check duplications. knownCanonical does have dup tx for the same gene.
			// Besides, those geneSymbol without knownCanonical may have duplicates too.
			// see https://lists.soe.ucsc.edu/pipermail/genome/2010-April/021963.html
			if (exist_element(gene_bySymbol,newg.name2))
			{
				gene_info& oldg = gene_bySymbol[newg.name2].begin()->second;
				if (newg.cdsStat > oldg.cdsStat)
				{
					gene_bySymbol[newg.name2].clear();
					gene_bySymbol[newg.name2].insert(pair<lenStat_t,gene_info>(newg.lenStat,newg));
				}
				else if (newg.cdsStat == oldg.cdsStat)
				{
					if (read_tx)
						gene_bySymbol[newg.name2].insert(pair<lenStat_t,gene_info>(newg.lenStat,newg));
					else if (newg.lenStat == oldg.lenStat)
						gene_bySymbol[newg.name2].insert(pair<lenStat_t,gene_info>(newg.lenStat,newg));
					else if (newg.lenStat  > oldg.lenStat)
					{
						gene_bySymbol[newg.name2].clear();
						gene_bySymbol[newg.name2].insert(pair<lenStat_t,gene_info>(newg.lenStat,newg));
					}
					else continue;
				}
				else continue;
			}
			else
				gene_bySymbol[newg.name2].insert(pair<lenStat_t,gene_info>(newg.lenStat,newg));
		}
		
		for (auto &gene:gene_bySymbol) for (auto &tx:gene.second) gene_byTranscript[tx.second.name].insert(pair<lenStat_t,gene_info>(tx.second.lenStat,tx.second));
		for (auto &gene:gene_bySymbol) for (auto &tx:gene.second) gene_db[tx.second.name]=tx.second;
		lns<<showl<<gene_bySymbol.size()<<" genes read."<<flush_logger;
	}
	
	void read_genes()
	{
		read_UCSCgenes(MainPath, GDB_name, !canonic, filt_tx);
	}
	
	void read_UCSCgenes(string dir, const std::string& type, bool read_tx, bool filter)
	{
		if (!str_endsw(dir,"/")) dir+='/';
		if		(str_has(type,"refGene"))	read_refGene  (dir+type+".txt",read_tx,filter);
		else if (str_has(type,"ensGene"))	read_refGene  (dir+type+".txt",read_tx,filter);
		else if (str_has(type,"knownGene")) read_knownGene(dir,read_tx);
		else exit_error("Unkown UCSC gene DB type. Must be refGene/ensGene/knownGene.");
	}
	
	void write_knownGene(const string& filename, bool sort_tx_len)
	{
		if (sort_tx_len)
		{
			for (auto &gene:gene_bySymbol)
			{
				tx_set_t new_set;
				for (auto &tx:gene.second)
				{
					lenStat_t new_len = std::pair<int,int>(tx.first.second,tx.first.first);
					new_set.insert(pair<lenStat_t,gene_info>(new_len,tx.second));
				}
				gene.second = new_set;
			}
		}

		lns<<showl<<"Writing "<<filename<<flush_logger;
		openOutFile_or_exit(out,filename);
		for (auto &gene:gene_bySymbol)
			for (auto &tx:gene.second)
			{
				gene_info& g = tx.second;
				out<<g.name<<"\t"<<g.chr<<'\t'<<g.strand<<'\t'<<g.txStart<<'\t'<<g.txEnd<<'\t'<<g.cdsStart<<'\t'<<g.cdsEnd<<'\t';
				out<<g.exonCount<<'\t';
				print_container(g.exonStarts,out,',');
				out<<",\t";
				print_container(g.exonEnds,out,',');
				out<<",\t";
				out<<g.proteinID<<'\t'<<g.alignID<<endl;
			}
		closefile(out);
	}
	
	void write_universal(const string& filename, bool sort_tx_len)
	{
		if (sort_tx_len)
		{
			for (auto &gene:gene_bySymbol)
			{
				tx_set_t new_set;
				for (auto &tx:gene.second)
				{
					lenStat_t new_len = std::pair<int,int>(tx.first.second,tx.first.first);
					new_set.insert(pair<lenStat_t,gene_info>(new_len,tx.second));
				}
				gene.second = new_set;
			}
		}
		
		lns<<showl<<"Writing "<<filename<<flush_logger;
		openOutFile_or_exit(out,filename);
		for (auto &gene:gene_bySymbol)
			for (auto &tx:gene.second)
			{
				gene_info& g = tx.second;
				_restore_cds(g);
				if (g.symbol.empty())	out<<".\t";
				else					out<<g.symbol<<"\t";
				out<<g.name;
				out<<"\t"<<g.chr<<'\t'<<g.strand<<'\t'<<g.txStart<<'\t'<<g.txEnd<<'\t';
				out<<g.cdsStart<<'\t'<<g.cdsEnd<<'\t';
				out<<g.exonCount<<'\t';
				print_container(g.exonStarts,out,',');
				out<<",\t";
				print_container(g.exonEnds,out,',');
				out<<",\t.\t";
				out<<g.name2;
				out<<'\t'<<g.cdsStartStat<<'\t'<<g.cdsEndStat<<'\t'<<(g.col16.empty()?".":g.col16)<<endl;
			}
		closefile(out);
	}
	
	// ------------------ DNA sequence ------------------

	// The following is not standard AAs, so they are removed.
	// AA in some species coded for by codons that are usually interpreted as stop codons: U O
	// chemical or crystallographic analysis cannot conclusively determine the residue: B Z J X
	// Previously X mean unkown amino acid "Unk". Since many people use it to represent Ter,
	// on 2015-11-19 I removed ("Xaa",'X') ("Unk",'X') and add ("X",'*') ('X',"Ter") instead of ("X",'X') ('X',"Unk").
	std::map< std::string, char > AA_3to1 = boost::assign::map_list_of\
	("Ala",'A')
	("Arg",'R')
	("Asn",'N')
	("Asp",'D')
	("Cys",'C')
	("Glu",'E')
	("Gln",'Q')
	("Gly",'G')
	("His",'H')
	("Ile",'I')
	("Leu",'L')
	("Lys",'K')
	("Met",'M')
	("Phe",'F')
	("Pro",'P')
	("Ser",'S')
	("Thr",'T')
	("Trp",'W')
	("Tyr",'Y')
	("Val",'V')
	("Ter",'*')
//	("Sec",'U')
//	("Pyl",'O')
//	("Asx",'B')
//	("Glx",'Z')
//	("Xle",'J')
	("A",'A')
	("R",'R')
	("N",'N')
	("D",'D')
	("C",'C')
	("E",'E')
	("Q",'Q')
	("G",'G')
	("H",'H')
	("I",'I')
	("L",'L')
	("K",'K')
	("M",'M')
	("F",'F')
	("P",'P')
	("S",'S')
	("T",'T')
	("W",'W')
	("Y",'Y')
	("V",'V')
//	("U",'U')
//	("O",'O')
//	("B",'B')
//	("Z",'Z')
//	("J",'J')
	("X",'*')
	("*",'*');
	
	std::map< char, std::string > AA_1to3 = boost::assign::map_list_of\
	('A',"Ala")
	('R',"Arg")
	('N',"Asn")
	('D',"Asp")
	('C',"Cys")
	('E',"Glu")
	('Q',"Gln")
	('G',"Gly")
	('H',"His")
	('I',"Ile")
	('L',"Leu")
	('K',"Lys")
	('M',"Met")
	('F',"Phe")
	('P',"Pro")
	('S',"Ser")
	('T',"Thr")
	('W',"Trp")
	('Y',"Tyr")
	('V',"Val")
	('*',"Ter")
	('X',"Ter");

	string tr_AA_1to3(const char& in)
	{
		if (!exist_element(AA_1to3,in)) exit_error("Unknown amino acid "+s(in));
		return AA_1to3[in];
	}
	
	string tr_AA_1to3(const string& in)
	{
		string out;
		for (auto &c:in) {
			if (!exist_element(AA_1to3,c)) exit_error("Unknown amino acid "+s(c));
			out += AA_1to3[c]; }
		return out;
	}
	
	// Below use * to represent stop codon instead of X, as recommended by hgvs.org; X is used to represent unknown AA (wikipedia).
	std::map< std::string, char > genet_code = boost::assign::map_list_of\
	("TTT",'F') ("TTC",'F') ("TTA",'L') ("TTG",'L') ("CTT",'L') ("CTC",'L') ("CTA",'L') ("CTG",'L') ("ATT",'I') ("ATC",'I') ("ATA",'I') ("ATG",'M')\
	("GTT",'V') ("GTC",'V') ("GTA",'V') ("GTG",'V') ("TCT",'S') ("TCC",'S') ("TCA",'S') ("TCG",'S') ("CCT",'P') ("CCC",'P') ("CCA",'P') ("CCG",'P')\
	("ACT",'T') ("ACC",'T') ("ACA",'T') ("ACG",'T') ("GCT",'A') ("GCC",'A') ("GCA",'A') ("GCG",'A') ("TAT",'Y') ("TAC",'Y') ("TAA",'*') ("TAG",'*')\
	("CAT",'H') ("CAC",'H') ("CAA",'Q') ("CAG",'Q') ("AAT",'N') ("AAC",'N') ("AAA",'K') ("AAG",'K') ("GAT",'D') ("GAC",'D') ("GAA",'E') ("GAG",'E')\
	("TGT",'C') ("TGC",'C') ("TGA",'*') ("TGG",'W') ("CGT",'R') ("CGC",'R') ("CGA",'R') ("CGG",'R') ("AGT",'S') ("AGC",'S') ("AGA",'R') ("AGG",'R')\
	("GGT",'G') ("GGC",'G') ("GGA",'G') ("GGG",'G')\
	// now RNA
	("UUU",'F') ("UUC",'F') ("UUA",'L') ("UUG",'L') ("CUU",'L') ("CUC",'L') ("CUA",'L') ("CUG",'L') ("AUU",'I') ("AUC",'I') ("AUA",'I') ("AUG",'M')\
	("GUU",'V') ("GUC",'V') ("GUA",'V') ("GUG",'V') ("UCU",'S') ("UCC",'S') ("UCA",'S') ("UCG",'S') ("CCU",'P') ("CCC",'P') ("CCA",'P') ("CCG",'P')\
	("ACU",'T') ("ACC",'T') ("ACA",'T') ("ACG",'T') ("GCU",'A') ("GCC",'A') ("GCA",'A') ("GCG",'A') ("UAU",'Y') ("UAC",'Y') ("UAA",'*') ("UAG",'*')\
	("CAU",'H') ("CAC",'H') ("CAA",'Q') ("CAG",'Q') ("AAU",'N') ("AAC",'N') ("AAA",'K') ("AAG",'K') ("GAU",'D') ("GAC",'D') ("GAA",'E') ("GAG",'E')\
	("UGU",'C') ("UGC",'C') ("UGA",'*') ("UGG",'W') ("CGU",'R') ("CGC",'R') ("CGA",'R') ("CGG",'R') ("AGU",'S') ("AGC",'S') ("AGA",'R') ("AGG",'R')\
	("GGU",'G') ("GGC",'G') ("GGA",'G') ("GGG",'G');

	std::string translate_cds(const std::string& DNA, bool stop_at_stop)
	{
		std::string result;
		for (size_t loc=0, len=DNA.size(); len-loc >= 3; loc+=3)
		{
			string codon=to_upper_copy(DNA.substr(loc,3));
			if (exist_element(genet_code,codon))
			{
				char aa = genet_code[codon];
				result.push_back(aa);
				if (aa=='*' && stop_at_stop) break;
			}
			else
			{
				lns << showe << "Unable to translate CDS to protein sequence due to an unknown codon " << codon << fatal;
			}
		}
		return result;
	}
	
	// n N - w s W S are omitted.
	std::map<char,char> IUPAC_complement = boost::assign::map_list_of\
		('a','t')
		('c','g') 
		('g','c') 
		('t','a') 
		('A','T')
		('C','G') 
		('G','C') 
		('T','A') 
		('m','k')
		('k','m') 
		('r','y') 
		('y','r') 
		('b','v') 
		('d','h') 
		('h','d') 
		('v','b') 
		('M','K')
		('K','M') 
		('R','Y') 
		('Y','R') 
		('B','V') 
		('D','H') 
		('H','D') 
		('V','B')
		('-','-')
		('.','.')
		('n','n')
		('N','N');
	
	
	// http://en.wikipedia.org/wiki/Nucleic_acid_notation exclude U because it's DNA code
	bool dna_complement(std::string& seq, bool ExitIfErr)
	{
		bool no_error=true;
		for (size_t i=0;i<seq.size();i++)
		{
			std::map<char,char>::iterator it = IUPAC_complement.find(seq[i]);
			if (it!=IUPAC_complement.end()) seq[i]=it->second;
			else if (ExitIfErr) exit_error("wrong DNA nucleotide "+s(seq[i]));
			else no_error=false;
		}
		return no_error;
	}
	
	std::string dna_complement_copy(const std::string& seq, bool ExitIfErr)
	{
		string res=seq;
		dna_complement(res,ExitIfErr);
		return res;
	}
	
	bool dna_rc(std::string& seq, bool ExitIfErr)
	{
		std::reverse(seq.begin(),seq.end());
		return dna_complement(seq,ExitIfErr);
	}
	
	std::string dna_rc_copy(const std::string& seq, bool ExitIfErr)
	{
		string res=seq;
		dna_rc(res,ExitIfErr);
		return res;
	}
	
	// ------------------ genotypes ------------------

	void gtp_par::read(const string& s) // FORMAT examples, GT:AD:DP:GQ:PL from RGC, GT:ADS:DS:GP from IMPUTE2
	{
		clear();
		vector<string> fields;
		boost::split(fields,s,boost::is_any_of(":"));
		for (size_t i=0;i<fields.size();++i)
		{
			if (fields[i]=="DP") DP_col=i+1;
			if (fields[i]=="GQ") GQ_col=i+1;
			if (fields[i]=="GP") GP_col=i+1;
		}
	}

	string genotype::to_linkage(const char c) {
		if (c==mss)	return "0 0";
		if (c==HoR)	return "1 1";
		if (c==HoA)	return "2 2";
		if (c==P01)	return "2 1";
		if (c==P10)	return "1 2";
		if (c==Het)	return "1 2";
		if (c==HpM)	return "0 0";
		if (c==HpR)	return "1 1";
		if (c==HpA)	return "2 2";
		return "0 0"; // c==CNA
	}

	string genotype::to_vcf(const char c) {
		if (c==mss)	return "./.";
		if (c==HoR)	return "0|0";
		if (c==HoA)	return "1|1";
		if (c==P01)	return "1|0";
		if (c==P10)	return "0|1";
		if (c==Het)	return "0/1";
		if (c==HpM)	return ".|.";
		if (c==HpR)	return "0|0";
		if (c==HpA)	return "1|1";
		return "."; // c==CNA
	}

	bool   genotype::valid(const char c) { return c!=CNA; }
	bool   genotype::is_missing(const char c) { return c==mss || c==HpM; }
	bool   genotype::has_missing(const string& in) { return in.find(mss)!=std::string::npos || in.find(HpM)!=std::string::npos; }
	
	double genotype::prob_1(const char c) { if (c==Het) return 0.5; else if (c==HoA||c==P01||c==HpA) return 1; else if (c==HoR||c==P10||c==HpR) return 0; else return std::numeric_limits<double>::signaling_NaN(); }
	double genotype::prob1_(const char c) { if (c==Het) return 0.5; else if (c==HoA||c==P10||c==HpA) return 1; else if (c==HoR||c==P01||c==HpR) return 0; else return std::numeric_limits<double>::signaling_NaN(); }
	
	string	genotype::help_text() {
		string result="FOR GENOTYPE READING\n";
		result += " --filt-DP I       Exclude genotypes if DP<I (0=no_filter) {"+itos(DP_cut)+"}\n";
		result += " --filt-GQ I       Exclude genotypes if GQ<I (0=no_filter) {"+itos(GQ_cut)+"}\n";
		result += " --filt-GP D       Exclude genotypes if GP<D (0=no_filter) {"+ftos(GP_cut)+"}\n";
		result += " --filt-domGP B    Exclude genotypes by GP with a dominant model {"+str_YesOrNo(GP_dom)+"}\n";
		result += " --filt-recGP B    Exclude genotypes by GP with a recessive model {"+str_YesOrNo(GP_rec)+"}\n";
		result += " --no-miss         Treat missing genotype as homozygous ref. allele\n";
		return result;
	}

	void genotype::read_arguments(vector<string>& srce_opt, size_t start, bool expected, bool stop_at_unknown) {
		vector<string> dest_opt;
		int tot_read = expected; // if (expected) read @ start or ended; else, find the 1st to read.
		for (size_t ended=0, argi=0, args=srce_opt.size(); argi<args; ++argi)
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

	int genotype::ReadOpt(vector<string>& srce_opt, size_t& argi) {
		int num_read=0;
		for (size_t args=srce_opt.size(); argi<args; ++argi)
		{
			if (srce_opt[argi]=="--no-miss")
			{
				NoMiss=true;
				++num_read;
			}
			else if (str_startsw(srce_opt[argi],"--no-miss="))
			{
				string input = srce_opt[argi].substr(10);
				try { NoMiss=IsYes(input); }
				catch (input_exception &) { exit_error("Cannot read \""+input+"\" as a boolean value."); }
				++num_read;
			}
			else if (str_startsw(srce_opt[argi],"--filt-DP"))
			{
				std::size_t found = srce_opt[argi].find('=');
				string input;
				if (found!=std::string::npos) input = srce_opt[argi].substr(found+1);
				else if (args-argi<2) exit_error("Insufficient argument for the --filt-DP option.");
				else input = srce_opt[++argi];
				try { DP_cut = boost::lexical_cast<int>(input); }
				catch (...) { exit_error("Cannot read "+input+" as an integer."); }
				++num_read;
			}
			else if (str_startsw(srce_opt[argi],"--filt-GQ"))
			{
				std::size_t found = srce_opt[argi].find('=');
				string input;
				if (found!=std::string::npos) input = srce_opt[argi].substr(found+1);
				else if (args-argi<2) exit_error("Insufficient argument for the --filt-GQ option.");
				else input = srce_opt[++argi];
				try { GQ_cut = boost::lexical_cast<int>(input); }
				catch (...) { exit_error("Cannot read "+input+" as an integer."); }
				++num_read;
			}
			else if (str_startsw(srce_opt[argi],"--filt-GP"))
			{
				std::size_t found = srce_opt[argi].find('=');
				string input;
				if (found!=std::string::npos) input = srce_opt[argi].substr(found+1);
				else if (args-argi<2) exit_error("Insufficient argument for the --filt-GP option.");
				else input = srce_opt[++argi];
				try { GP_cut = boost::lexical_cast<double>(input); }
				catch (...) { exit_error("Cannot read "+input+" as an integer."); }
				++num_read;
			}
			else if (str_startsw(srce_opt[argi],"--filt-domGP="))
			{
				string input = srce_opt[argi].substr(13);
				try { GP_dom=IsYes(input); }
				catch (input_exception &) { exit_error("Cannot read \""+input+"\" as a boolean value."); }
				++num_read;
			}
			else if (str_startsw(srce_opt[argi],"--filt-recGP="))
			{
				string input = srce_opt[argi].substr(13);
				try { GP_rec=IsYes(input); }
				catch (input_exception &) { exit_error("Cannot read \""+input+"\" as a boolean value."); }
				++num_read;
			}
			else
			{
				break;
			}
		}
		--argi; // Always point to the option left to the unknown option. Good for "for (argi=1; argi<args; ++argi)".
		return num_read;
	}

	char   genotype::fix_ploidy(const char c, const int exp_ploidy, int& ploidy_err)
	{
		int obs_ploidy=ploidy(c);
		if (exp_ploidy==0) return CNA;
		if (exp_ploidy==1)
		{
			if (obs_ploidy==1) return c;
			if (obs_ploidy==2)
			{
				if (c==mss)	{				return HpM; }
				if (c==HoR)	{				return HpR; }
				if (c==HoA)	{				return HpA; }
				if (c==P01)	{ ++ploidy_err; return HpM; }
				if (c==P10)	{ ++ploidy_err;	return HpM; }
				if (c==Het)	{ ++ploidy_err;	return HpM; }
			}
			exit_error(" cannot convert N/A to haploidy.");
		}
		if (exp_ploidy==2)
		{
			if (obs_ploidy==2) return c;
			if (obs_ploidy==1)
			{
				// previously exit_error(" cannot convert haploidy to diploidy.");
				// Most likely explanation: IMPUTE2 imputes variants on chrX PAR but treats it as non-PAR, which is a mistake.
				// Therefore, the genotype call is likely wrong, hence cannot be converted to a diploidy.
				// Therefore, I will make them all missing. And I wil not count it as ploidy_err for this sample as it's due to software error.
				return mss;
			}
			exit_error(" cannot convert N/A to diploidy.");
		}
		exit_error(" do not recognize expected_ploidy "+itos(exp_ploidy));
		return CNA;
	}
	bool genotype::vcf_gtp(const string& s)
	{
		if (s.find(':')!=string::npos) return true;
		if (s.size()==3 && (s[1]=='/'||s[1]=='|')) return true;
		return false;
	}
	char   genotype::to_char(const string& s, const bool to_filter, int& GQ, int& DP, double& GP, const gtp_par& par)
	{
		GQ=-1; DP=-1; GP=-1;
		string geno = substr_before_find(s,":");
		if (geno=="./." || geno==".|.") return NoMiss ? HoR : mss;
		
		// quality control
		bool is_vcf = false;
		if (vcf_gtp(s)) // if VCF
		{
			is_vcf = true;
			vector<string> fields;
			boost::split(fields,s,boost::is_any_of(":"));
			bool bad = false;
			if (par.DP_col && DP_cut) // if need to do QC by DP
			{
				if (fields.size()>=par.DP_col && read_val_ge(fields[par.DP_col-1],DP,0))
				{
					if (DP<DP_cut) bad=true;
					DP_tot += DP;
					++DP_cnt;
				}
			}
			if (par.GQ_col && GQ_cut) // if need to do QC by GQ
			{
				if (fields.size()>=par.GQ_col && read_val_ge(fields[par.GQ_col-1],GQ,0))
				{
					if (GQ<GQ_cut) bad=true;
					GQ_tot += GQ;
					++GQ_cnt;
				}
			}
			if (par.GP_col && GP_cut) // if need to do QC by GP
			{
				if (fields.size()>=par.GP_col)
				{
					vector<string> valstr;
					boost::split(valstr,fields[par.GP_col-1],boost::is_any_of(","));
					if (valstr.size()==3)
					{
						vector<double> values;
						for (auto &x:valstr)
						{
							double v;
							if (read_val_noNaN(x,v)) values.push_back(v);
						}
						if (values.size()==3)
						{
							if 		(GP_dom) GP=std::max( values[0],values[1]+values[2]);
							else if (GP_rec) GP=std::max( values[0]+values[1],values[2]);
							else			 GP=std::max({values[0],values[1],values[2]});
							if (GP<GP_cut) bad=true;
							GP_tot += GP;
							++GP_cnt;
						}
					}
				}
			}
			if (to_filter && bad) boost::replace_if(geno, ::isalnum, '.'); // geno.clear();
		}
		
		// read genotype
		if (geno.empty() || geno=="NA" || geno=="N/A" || geno=="#N/A")
		{
			if (is_vcf) exit_error("Cannot read genotype from "+s+" because I don't know the ploidy.");
			else		return NoMiss ? HoR : mss;
		}
		if (geno==".")
		{
			if (is_vcf) return NoMiss ? HpR : HpM;
			else		return NoMiss ? HoR : mss;
		}
		if (geno.size()==3) // diploidy
		{
			char a1 = geno[0];
			char a2 = geno[2];
			if (a1=='.' && NoMiss) a1='0';
			if (a2=='.' && NoMiss) a2='0';
			if (a1=='.' || a2=='.') return mss;
			if (a1=='1' && a2=='1') return HoA;
			if (a1!='1' && a2!='1') return HoR;
			if (geno[1]=='/')		return Het;
			if (a1=='1')			return P01;
			if (a2=='1')			return P10;
		}
		else if (geno.size()==1)
		{
			if (is_vcf) // haploidy
			{
				char a1 = geno[0];
				if		(a1=='0')	return HpR;
				else if	(a1=='1')	return HpA;
				else if (a1=='.')	return NoMiss ? HpR : HpM;
			}
			else // diploidy alt count
			{
				char a1 = geno[0];
				if		(a1=='0')	return HoR;
				else if	(a1=='1')	return Het;
				else if (a1=='2')	return HoA;
				else if (a1=='.')	return NoMiss ? HoR : mss;
			}
		}
		exit_error("Error reading genotype "+s);
		return mss; // won't happen
	}
	
	void genotype::to_void(string& s)
	{
		if (vcf_gtp(s)) // VCF
		{
			string::iterator it;
			for (it = s.begin(); it < s.end(); ++it)
			{
				if (*it==':') break;
				if (isdigit(*it)) *it='.';
			}
		}
		else if (s.size()==3) // diploidy
		{
			s="./.";
		}
		else if (s.size()==1) // haploidy or #copy of ALT
		{
			s=".";
		}
	}
	
	void genotype::to_rewr(char g, string& s)
	{
		if (!usable(g))
		{
			if (vcf_gtp(s)) // VCF
			{
				string::iterator it;
				for (it = s.begin(); it < s.end(); ++it)
				{
					if (*it==':') break;
					if (isdigit(*it)) *it='.';
				}
			}
			else if (s.size()==3) // diploidy
			{
				s="./.";
			}
			else if (s.size()==1) // haploidy or #copy of ALT
			{
				s=".";
			}
		}
		else if (NoMiss) // if g is not missing, s is missing, and user don't want any missing
		{
			if (vcf_gtp(s) && s[0]=='.') s=to_vcf(g);
			// else could be haploidy / #copy_of_alt. Don't know what to do.
		}
	}

	void genotype::to_dn(string& s)
	{
		if (vcf_gtp(s)) // VCF, robust to any ploidy
		{
			for (size_t i=0;i<s.size();++i)
			{
				if		(s[i]==':') break;
				else if (s[i]=='|') s[i]='/';
				else if (i)			s[i]='0';
				else				s[i]='1';
			}
		}
		else if (s.size()==3) // diploidy
		{
			s="1/0";
		}
		else if (s.size()==1) // haploidy or #copy of ALT
		{
			s="1";
		}
	}

	void genotype::to_ref(string& s)
	{
		if (vcf_gtp(s)) // VCF
		{
			string::iterator it;
			for (it = s.begin(); it < s.end(); ++it)
			{
				if (*it==':') break;
				if (*it!='|' && *it!='/') *it='0';
			}
		}
		else if (s.size()==3) // diploidy
		{
			s="0/0";
		}
		else if (s.size()==1) // haploidy or #copy of ALT
		{
			s="0";
		}
	}
	
	void genotype::to_swap(string& s) // VCF only; not robust to multi-ALT variants after "vSPLIT -a=0".
	{
		if (vcf_gtp(s)) // VCF
		{
			string::iterator it;
			for (it = s.begin(); it < s.end(); ++it)
			{
				if		(*it==':') break;
				else if	(*it=='1') *it='0';
				else if	(*it=='0') *it='1';
			}
		}
	}
	
	int		genotype::mean_GQ() { if (GQ_cnt) return (double)GQ_tot/GQ_cnt; else return -1; }
	int 	genotype::mean_DP() { if (DP_cnt) return (double)DP_tot/DP_cnt; else return -1; }
	double	genotype::mean_GP() { if (GP_cnt) return (double)GP_tot/GP_cnt; else return -1; }
	void genotype::rewind()
	{
		DP_tot=0;
		GQ_tot=0;
		GP_tot=0;
		DP_cnt=0;
		GQ_cnt=0;
		GP_cnt=0;
	}
	
	int genotype::num_alt(const char c)
	{
		if (c==HoR)		return 0;
		if (c==HoA)		return 2;
		if (c==P01)		return 1;
		if (c==P10)		return 1;
		if (c==Het)		return 1;
		if (c==HpR)		return 0;
		if (c==HpA)		return 1;
		if (c==mss)		exit_error("did not test usability of genotype");
		if (c==HpM)		exit_error("did not test usability of genotype");
		if (c==CNA)		exit_error("did not test usability of genotype");
		exit_error("genepi::genotype::num_alt() doesn't recognize "+s(c));
		return -999999; // won't happen
	}
	int genotype::num_alt_recessive(const char c)
	{
		if (c==HoR)		return 0;
		if (c==HoA)		return 2;
		if (c==P01)		return 1;
		if (c==P10)		return 1;
		if (c==Het)		return 1;
		if (c==HpR)		return 0;
		if (c==HpA)		return 2; // the only difference from above
		if (c==mss)		exit_error("did not test usability of genotype");
		if (c==HpM)		exit_error("did not test usability of genotype");
		if (c==CNA)		exit_error("did not test usability of genotype");
		exit_error("genepi::genotype::num_alt_recessive() doesn't recognize "+s(c));
		return -999999; // won't happen
	}
	int genotype::ploidy(const char c)
	{
		if (c==mss)		return 2;
		if (c==HoR)		return 2;
		if (c==HoA)		return 2;
		if (c==P01)		return 2;
		if (c==P10)		return 2;
		if (c==Het)		return 2;
		if (c==HpM)		return 1;
		if (c==HpR)		return 1;
		if (c==HpA)		return 1;
		if (c==CNA)		return 0;
		exit_error("genepi::genotype::ploidy() doesn't recognize "+s(c));
		return -999999; // won't happen
	}
	bool genotype::denovo(const char g_ch, const char g_pa, const char g_ma, bool allow_hom)
	{
		if (valid(g_pa)) { if (is_missing(g_pa)) return false; if (num_alt(g_pa)!=0) return false; }
		if (valid(g_ma)) { if (is_missing(g_ma)) return false; if (num_alt(g_ma)!=0) return false; }
		if (valid(g_ch)) { if (is_missing(g_ch)) return false;
			// now pa,ma is non-missing and is 0
			int g=num_alt(g_ch);
			if		(g==1)				return true;
			else if (g==2 && allow_hom) return true;
			else						return false; }
		return false;
	}
	bool genotype::MenErr_auto(const char g_ch, const char g_pa, const char g_ma) // including denovo!
	{
		if (g_ch==mss) return false;
		if (g_ch==HoR)								return  (g_pa==HoA || g_ma==HoA);
		if (g_ch==HoA)								return  (g_pa==HoR || g_ma==HoR);
		if (g_ch==Het || g_ch==P01 || g_ch==P10)	return ((g_pa==HoR && g_ma==HoR)||(g_pa==HoA && g_ma==HoA));
		exit_error("MenErr_auto() encounter a problem: the child's ploidy is wrong");
		return true;
	}
	bool genotype::MenErr_chrX(const char g_ch, const char g_pa, const char g_ma) // including denovo!
	{
		if (g_ch==mss || g_ch==HpM) return false;
		if (ploidy(g_ch)==2) // PAR or female
		{
			if (ploidy(g_pa)==2) // PAR
			{
				if (g_ch==HoR)								return  (g_pa==HoA || g_ma==HoA);
				if (g_ch==HoA)								return  (g_pa==HoR || g_ma==HoR);
				if (g_ch==Het || g_ch==P01 || g_ch==P10)	return ((g_pa==HoR && g_ma==HoR)||(g_pa==HoA && g_ma==HoA));
				return false; // will not happen
			}
			else // ploidy=1, no checking for CNA
			{
				if (g_ch==HoR)								return  (g_pa==HpA || g_ma==HoA);
				if (g_ch==HoA)								return  (g_pa==HpR || g_ma==HoR);
				if (g_ch==Het || g_ch==P01 || g_ch==P10)	return ((g_pa==HpR && g_ma==HoR)||(g_pa==HpA && g_ma==HoA));
				return false; // will not happen
			}
		}
		else // male, ploidy of g_pa must be 1 too, this X is from mom. no checking for CNA.
		{
			if (g_ch==HpR)	return (g_ma==HoA);
			if (g_ch==HpA)	return (g_ma==HoR);
			return false; // will not happen
		}
		exit_error("MenErr_auto() encounter a problem: the child's ploidy is wrong"); // will not happen
		return true;
	}
	bool genotype::MenErr_chrY(const char g_ch, const char g_pa, const char g_ma) // including denovo!
	{
		if (g_ch==CNA || g_ch==HpM || g_ch==mss) return false;
		if (ploidy(g_ch)==2) // PAR
		{
			if (g_ch==HoR)								return  (g_pa==HoA || g_ma==HoA);
			if (g_ch==HoA)								return  (g_pa==HoR || g_ma==HoR);
			if (g_ch==Het || g_ch==P01 || g_ch==P10)	return ((g_pa==HoR && g_ma==HoR)||(g_pa==HoA && g_ma==HoA));
			return false; // will not happen
		}
		else // male, ploidy of g_pa must be 1 too, this Y from dad. no need to check for CNA.
		{
			if (g_ch==HpR)	return (g_pa==HpA);
			if (g_ch==HpA)	return (g_pa==HpR);
			return false; // will not happen
		}
	}
	bool genotype::MenErr_chrM(const char g_ch, const char g_pa, const char g_ma) // including denovo!
	{
		if (ploidy(g_ch)!=1) exit_error("MenErr_auto() encounter a problem: the child's ploidy is wrong");
		if (g_ch==HpM) return false;
		if (g_ch==HpR) return (g_ma==HpA);
		if (g_ch==HpA) return (g_ma==HpR);
		return false;
	}
	void genotype::swap_allele(char& c)
	{
		if		(c==HoR) c=HoA;
		else if (c==HoA) c=HoR;
		else if (c==P01) c=P10;
		else if (c==P10) c=P01;
		else if (c==HpR) c=HpA;
		else if (c==HpA) c=HpR;
	}
	void genotype::dephase(char& c)
	{
		if		(c==P01) c=Het;
		else if (c==P10) c=Het;
	}
	
	// Does not maintain missing! Because otherwise some individuals may be excluded from analysis.
	// But keep NA. Because these people are supposed to be excluded anyway.
	void genotype::set_all_ref(char& c)
	{
		int p = ploidy(c);
		if		(p==2) c=HoR;
		else if (p==1) c=HpR;
	}

	void genotype::set_missing(char& c)
	{
		int p = ploidy(c);
		if		(p==2) c=mss;
		else if (p==1) c=HpM;
	}

	char genotype::read_vcf(const string& s, const int chr_num, const int bp, const int sex, const bool to_filter, int& GQ, int& DP, double& GP, const gtp_par& par)
	{
		int ploidy_err=0;
		return fix_ploidy(to_char(s,to_filter,GQ,DP,GP,par),expected_ploidy(chr_num,bp,sex),ploidy_err);
	}
	
	char genotype::read_vcf(const string& s, const bool to_filter, int& GQ, int& DP, double& GP, const gtp_par& par)
	{
		return to_char(s,to_filter,GQ,DP,GP,par);
	}

	char genotype::read(const string& s, const int chr_num, const int bp, const int sex, const gtp_par& par)
	{
		int ploidy_err=0;
		int GQ,DP;
		double GP;
		return fix_ploidy(to_char(s,true,GQ,DP,GP,par),expected_ploidy(chr_num,bp,sex),ploidy_err);
	}
	
	char genotype::read(const string& s, const gtp_par& par)
	{
		int GQ,DP;
		double GP;
		return to_char(s,true,GQ,DP,GP,par);
	}

	void genotype::count(const char g, int& i, int& m, double& n, double& x, const double wt)
	{
		if (valid(g))
		{
			++i;
			if (is_missing(g)) ++m;
			else { n+=ploidy(g)*wt; x+=num_alt(g)*wt; }
		}
	}
	
	// Default values for DP (8) and GQ (20) was taken from Carson 2014 (Effective filtering strategies to improve data quality .. BMC Bioinfo. 15:125. doi: 10.1186/1471-2105-15-125)
	// Previously DP_cut=20 (the default value by GATK) and GQ_cut=30 (the default value by GATK & Somatic Variant Caller).
	// GATK comments, "by GQ = 26 isn't very good", at http://gatkforums.broadinstitute.org/discussion/1268/what-is-a-vcf-and-how-should-i-interpret-it
	// In RGC, 30/12 separate the concordance/discordance genotypes.
	// In my experience with the TCGA-KG data, 30/12 results in a bad QQ plot, while 40/8 is fine. I also used 50/20, the QQ plot looks the same.
	// The Dutch pipeline use 30/7, which suggests that GQ is quite important.
	int     genotype::GQ_tot=0;
	int     genotype::DP_tot=0;
	double	genotype::GP_tot=0;
	int     genotype::GQ_cnt=0;
	int     genotype::DP_cnt=0;
	int     genotype::GP_cnt=0;
	bool	genotype::NoMiss=false;
	int		genotype::DP_cut=10;
	int		genotype::GQ_cut=40;
	double	genotype::GP_cut=0.8;	// This threshold is obtained from https://www.med.unc.edu/pgc/files/resources/visualization-tools/imputation_ricopili
	bool	genotype::GP_dom=false; // QC in dominant model,  can only be used for analysis under a dominant model
	bool	genotype::GP_rec=false; // QC in recessive model, can only be used for analysis under a recessive model

	// The following 3 functions are over-simplified because it returns one MAF instead of probability of 3 genotypes.
	// However, they have the benefit of ease of use and doesn't inflate type-I-error for association tests require HWE.
	// calculate MAF in cases given population MAF (af) and a penetrance model (penetr)
	double up_af(const double& af, const vector<double>& penetr)
	{
		if (penetr.size()!=3) exit_error("penetrance should be 3 floating point numbers.");
		double Paa = (1-af)*(1-af);
		double PaA = 2*af*(1-af);
		double PAA = af*af;
		double den_A = penetr[0]*Paa + penetr[1]*PaA + penetr[2]*PAA;	// denominator for cases
		double PaA_A = penetr[1]*PaA / den_A;
		double PAA_A = penetr[2]*PAA / den_A;
		double piA = (PaA_A + 2*PAA_A)/2;
		return piA;
	}
	
	// calculate MAF in controls given population MAF (af) and a penetrance model (penetr)
	double dn_af(const double& af, const vector<double>& penetr)
	{
		if (penetr.size()!=3) exit_error("penetrance should be 3 floating point numbers.");
		double Paa = (1-af)*(1-af);
		double PaA = 2*af*(1-af);
		double PAA = af*af;
		double den_U = (1-penetr[0])*Paa + (1-penetr[1])*PaA + (1-penetr[2])*PAA;
		double PaA_U = (1-penetr[1])*PaA / den_U;
		double PAA_U = (1-penetr[2])*PAA / den_U;
		double piU = (PaA_U + 2*PAA_U)/2;
		return piU;
	}

	void up_dn_af(const double& af, const vector<double>& penetr, double& piA, double& piU)
	{
		if (penetr.size()!=3) exit_error("penetrance should be 3 floating point numbers.");
		double Paa = (1-af)*(1-af);
		double PaA = 2*af*(1-af);
		double PAA = af*af;
		double den_A = penetr[0]*Paa + penetr[1]*PaA + penetr[2]*PAA;	// denominator for cases
		double PaA_A = penetr[1]*PaA / den_A;
		double PAA_A = penetr[2]*PAA / den_A;
		piA = (PaA_A + 2*PAA_A)/2;
		double den_U = (1-penetr[0])*Paa + (1-penetr[1])*PaA + (1-penetr[2])*PAA;
		double PaA_U = (1-penetr[1])*PaA / den_U;
		double PAA_U = (1-penetr[2])*PAA / den_U;
		piU = (PaA_U + 2*PAA_U)/2;
	}

	// The following two functions are more accurate in that it takes 3 probabilities and output another 3, not depending on HWE.
	// calculate MAF in cases given population MAF (af) and a penetrance model (penetr)
	std::tuple<double,double,double> up_af(const std::tuple<double,double,double>& af, const vector<double>& penetr)
	{
		if (penetr.size()!=3) exit_error("penetrance should be 3 floating point numbers.");
		double Paa = std::get<0>(af);
		double PaA = std::get<1>(af);
		double PAA = std::get<2>(af);
		double den_A = penetr[0]*Paa + penetr[1]*PaA + penetr[2]*PAA;	// denominator for cases
		double Paa_A = penetr[0]*Paa / den_A;
		double PaA_A = penetr[1]*PaA / den_A;
		double PAA_A = penetr[2]*PAA / den_A;
		return std::make_tuple(Paa_A,PaA_A,PAA_A);
	}
	
	// calculate MAF in controls given population MAF (af) and a penetrance model (penetr)
	std::tuple<double,double,double> dn_af(const std::tuple<double,double,double>& af, const vector<double>& penetr)
	{
		if (penetr.size()!=3) exit_error("penetrance should be 3 floating point numbers.");
		double Paa = std::get<0>(af);
		double PaA = std::get<1>(af);
		double PAA = std::get<2>(af);
		double den_U = (1-penetr[0])*Paa + (1-penetr[1])*PaA + (1-penetr[2])*PAA;
		double Paa_U = (1-penetr[0])*Paa / den_U;
		double PaA_U = (1-penetr[1])*PaA / den_U;
		double PAA_U = (1-penetr[2])*PAA / den_U;
		return std::make_tuple(Paa_U,PaA_U,PAA_U);
	}

	void up_dn_af(const std::vector<double>& af, const vector<double>& rr, const double& prevalence, std::vector<double>& cs, std::vector<double>& ct, char model)
	{
		// store temporary result to speed up
		vector< vector<double> > tmp_A(af.size(),vector<double>(af.size(),0));
		vector< vector<double> > tmp_B(af.size(),vector<double>(af.size(),0));
		
		// calculate r0 (baseline risk of a haplotype w/o any causal variants)
		double P_x_RR = 0;
		if		(model=='m')
		{
			for (size_t i=0;i<af.size();++i)
				for (size_t j=0;j<af.size();++j)
				{
					tmp_A[i][j] = af[i]*af[j];
					tmp_B[i][j] = tmp_A[i][j] * rr[i]*rr[j];
					P_x_RR += tmp_B[i][j];
				}
		}
		else if (model=='d')
		{
			for (size_t i=0;i<af.size();++i)
				for (size_t j=0;j<af.size();++j)
				{
					tmp_A[i][j] = af[i]*af[j];
					tmp_B[i][j] = tmp_A[i][j] * std::max(rr[i],rr[j]);
					P_x_RR += tmp_B[i][j];
				}
		}
		else if (model=='r')
		{
			for (size_t i=0;i<af.size();++i)
				for (size_t j=0;j<af.size();++j)
				{
					tmp_A[i][j] = af[i]*af[j];
					tmp_B[i][j] = tmp_A[i][j] * std::min(rr[i],rr[j]);
					P_x_RR += tmp_B[i][j];
				}
		}
		else if (model=='a')
		{
			for (size_t i=0;i<af.size();++i)
				for (size_t j=0;j<af.size();++j)
				{
					tmp_A[i][j] = af[i]*af[j];
					tmp_B[i][j] = tmp_A[i][j] * (rr[i]+rr[j]);
					P_x_RR += tmp_B[i][j];
				}
		}
		else
			exit_error("Model for up_dn_af() should be m/d/r/a.");
		double r0 = prevalence / P_x_RR;
		
		// calculate denominator. Inverse to speed up; 0.5 to make sum(cs)=sum(ct)=1.
		double inv_den_A = 0.5/prevalence;
		double inv_den_U = 0.5/(1-prevalence);
	
		// calcualte results
		cs.assign(af.size(),0);
		ct.assign(af.size(),0);
		for (size_t i=0;i<af.size();++i)
			for (size_t j=0;j<af.size();++j)
			{
				tmp_B[i][j] *= r0;
				double p =  tmp_B[i][j] * inv_den_A;
				cs[i] += p;
				cs[j] += p;
				double q = (tmp_A[i][j] - tmp_B[i][j]) * inv_den_U;
				ct[i] += q;
				ct[j] += q;
			}
	}

	void up_dn_af(const std::vector<double>& af, const vector<double>& rr, const double& prevalence, std::vector<double>& cs, std::vector<double>& ct)
	{
		// store temporary result to speed up
		cs.assign(af.size(),0);
		ct.assign(af.size(),0);
		
		// calculate r0 (baseline risk of a haplotype w/o any causal variants)
		double P_x_RR = 0;
		for (size_t i=0;i<af.size();++i)
		{
			cs[i] = af[i] * rr[i];
			P_x_RR += cs[i];
		}
		double r0 = prevalence / P_x_RR;
		
		// calculate denominator. Inverse to speed up
		double inv_den_A = 1/prevalence;
		double inv_den_U = 1/(1-prevalence);
		
		// calcualte results
		for (size_t i=0;i<af.size();++i)
		{
			cs[i] *= r0;
			ct[i] = (af[i] - cs[i]) * inv_den_U;
			cs[i] *= inv_den_A;
		}
	}

	void up_dn_af(const std::vector<double>& af, const vector<double>& pn, std::vector<double>& cs, std::vector<double>& ct)
	{
		cs.assign(af.size(),0);
		ct.assign(af.size(),0);

		// calculate prevalence
		double prevalence = 0;
		for (size_t i=0;i<af.size();++i)
		{
			cs[i] = af[i] * pn[i];
			prevalence += cs[i];
		}
		
		// calculate denominator for affected and unaffected. Inverse to speed up
		double inv_den_A = 1/prevalence;
		double inv_den_U = 1/(1-prevalence);
		
		// calcualte results
		for (size_t i=0;i<af.size();++i)
		{
			ct[i] = (af[i] - cs[i]) * inv_den_U;
			cs[i] *=  inv_den_A;
		}
	}

	// ------------------ linkage analysis ------------------

	HLOD_result cal_HLOD(const std::multiset<double>& lodmset)
	{
		double minAlpha=0, maxAlpha=1+std::numeric_limits<double>::epsilon()*20, LOD=0;
		for (each_element_const(lodmset,it)) LOD+=*it;
		double maxhlod=0;
		double maxa=0;
		for (;;)
		{
			double incremnt = (maxAlpha-minAlpha)/20;
			for (double a=minAlpha; a<=maxAlpha; a+=incremnt)
			{
				double hlod=0;
				for (each_element_const(lodmset,it)) hlod+=log10(a*pow(10,*it)+1-a);
				if (hlod>maxhlod) { maxhlod=hlod; maxa=a; }
			}
			if (incremnt<=0.001) break;
			minAlpha = std::max(minAlpha, maxa - incremnt);
			maxAlpha = std::min(maxAlpha, maxa + incremnt);
		}
		return HLOD_result(LOD,maxa,maxhlod);
	}
	
	HLOD_result cal_HLOD(const std::vector<double>& lodvec)
	{
		double minAlpha=0, maxAlpha=1+std::numeric_limits<double>::epsilon()*20, LOD=0;
		unsigned size=lodvec.size();
		for (unsigned i=0;i<size;i++) LOD+=lodvec[i];
		double maxhlod=0;
		double maxa=0;
		for (;;)
		{
			double incremnt = (maxAlpha-minAlpha)/20;
			for (double a=minAlpha; a<=maxAlpha; a+=incremnt)
			{
				double hlod=0;
				for (unsigned i=0;i<size;i++) hlod += log10( a*pow(10,lodvec[i]) + 1-a );
				if (hlod>maxhlod) { maxhlod=hlod; maxa=a; }
			}
			if (incremnt<=0.001) break;
			minAlpha = std::max(minAlpha, maxa - incremnt);
			maxAlpha = std::min(maxAlpha, maxa + incremnt);
		}
		return HLOD_result(LOD,maxa,maxhlod);
	}
	
	// calculate HLOD conditional on another locus, prbvec is the probability of linkage for the other locus/loci
	// note: L+/L- = LR = 10^LOD ; a1=alpha_when_link_to_other_loci ; a2=alpha_when_NOT_link_to_other_loci
	HLOD_result cal_HLOD_conditional(const std::vector<double>& lodvec, const std::vector<double>& prbvec)
	{
		double LOD=0;
		double mina1=0, maxa1=1+std::numeric_limits<double>::epsilon()*20;
		double mina2=0, maxa2=1+std::numeric_limits<double>::epsilon()*20;
		unsigned size=lodvec.size();
		for (unsigned i=0;i<size;i++) LOD+=lodvec[i];
		double besthlod=0,besta1=0,besta2=0,bestma=0;
		for (;;)
		{
			double incr1 = (maxa1-mina1)/20;
			double incr2 = (maxa2-mina2)/20;
			for (double a1=mina1; a1<=maxa1; a1+=incr1)
			{
				for (double a2=mina2; a2<=maxa2; a2+=incr2)
				{
					double hlod=0,a=0;
					for (unsigned i=0;i<size;i++)
					{
						double LR=pow(10,lodvec[i]);	// likelihood ratio for the studied locus, = L+/L-
						double PPL=prbvec[i];			// posterior probability of linkage to the other loci
						double ai=a1*PPL + a2*(1-PPL);	// alpha_i,  probability of linkage to the studied locus
						hlod += log10( ai*LR + 1-ai );	// same as traditional HLOD
						a += ai;
					}
					if (hlod>besthlod) { besthlod=hlod; besta1=a1; besta2=a2; bestma=a/size; } // ma is mean(alpha_i)
				}
			}
			if (incr1<=0.001) break;
			mina1 = std::max(mina1, besta1 - incr1);
			maxa1 = std::min(maxa1, besta1 + incr1);
			mina2 = std::max(mina2, besta2 - incr2);
			maxa2 = std::min(maxa2, besta2 + incr2);
		}
		return HLOD_result(LOD,bestma,besthlod);
	}
	
	// calculate HLOD with individual alpha for each fam, prbvec is the probability of linkage for the target locus
	HLOD_result cal_HLOD_iAlpha(const std::vector<double>& lodvec, const std::vector<double>& prbvec)
	{
		double minAlpha=0, maxAlpha=1+std::numeric_limits<double>::epsilon()*20, LOD=0;
		unsigned size=lodvec.size();
		for (unsigned i=0;i<size;i++) LOD+=lodvec[i];
		double maxhlod=0;
		double maxa=0;
		for (;;)
		{
			double incremnt = (maxAlpha-minAlpha)/20;
			for (double a=minAlpha; a<=maxAlpha; a+=incremnt)
			{
				double hlod=0;
				for (unsigned i=0;i<size;i++)
				{
					double ai=a*prbvec[i];
					hlod+=log10(ai*pow(10,lodvec[i])+1-ai);
				}
				if (hlod>maxhlod) { maxhlod=hlod; maxa=a; }
			}
			if (incremnt<=0.001) break;
			minAlpha = std::max(minAlpha, maxa - incremnt);
			maxAlpha = std::min(maxAlpha, maxa + incremnt);
		}
		return HLOD_result(LOD,maxa,maxhlod);
	}
	
	// calculate summation HLOD (overall HLOD for several loci as a whole, not useful for testing a specific locus)
	HLOD_result cal_HLOD_summation(const std::vector<double>& lodvec, const std::vector< std::vector<double> >& lodoth, const std::vector<double>& aj)
	{
		double minAlpha=0, maxAlpha=1+std::numeric_limits<double>::epsilon()*20, LOD=0;
		unsigned size=lodvec.size();
		for (unsigned i=0;i<size;i++) LOD+=lodvec[i];
		double maxhlod=0;
		double maxa=0;
		for (;;)
		{
			double incremnt = (maxAlpha-minAlpha)/20;
			for (double a=minAlpha; a<=maxAlpha; a+=incremnt)
			{
				double hlod=0;
				for (unsigned i=0;i<size;i++)
				{
					double l = 1 + a * pow(10,lodvec[i]) - a;
					for (unsigned j=0;j<aj.size();j++) l += aj[j] * ( pow(10,lodoth[j][i]) -1 );
					hlod+=log10(l);
				}
				if (hlod>maxhlod) { maxhlod=hlod; maxa=a; }
			}
			if (incremnt<=0.01) break;
			minAlpha = std::max(minAlpha, maxa - incremnt);
			maxAlpha = std::min(maxAlpha, maxa + incremnt);
		}
		return HLOD_result(LOD,maxa,maxhlod);
	}
	
	HLOD_result cal_HLOD_summation(const std::vector<double>& lodvec, const std::vector< std::vector<double> >& lodoth, std::vector<double>& result_aj, double incr_aj)
	{
		double minAlpha=0, maxAlpha=1+std::numeric_limits<double>::epsilon()*20;
		std::vector<double> aj(lodoth.size(),0.0);
		HLOD_result best_r;
		for (;;)
		{
			HLOD_result r = cal_HLOD_summation(lodvec,lodoth,aj);
			if (r.HLOD > best_r.HLOD) { best_r=r; result_aj=aj; }
			bool valid=false;	// next try is valid
			for (unsigned i=0; i<aj.size(); i++) {
				if (aj[i]+incr_aj > maxAlpha) { aj[i]=minAlpha; continue; }
				aj[i]+=incr_aj; valid=true; break;
			}
			if (!valid) break;
		}
		return best_r;
	} //*/
	
	/*/ calculate summation HLOD (overall HLOD for several loci as a whole, not useful for testing a specific locus).  assume sum(a)=1 !
	HLOD_result cal_HLOD_summation(const std::vector<double>& lodvec, const std::vector< std::vector<double> >& lodoth, const std::vector<double>& aj)
	{
		double minAlpha=0, maxAlpha=1-std::accumulate( aj.begin(), aj.end(), 0.0 ), LOD=0;
		unsigned size=lodvec.size();
		for (unsigned i=0;i<size;i++) LOD+=lodvec[i];
		double maxhlod=0;
		double maxa=0;
		for (;;)
		{
			double incremnt = (maxAlpha-minAlpha)/20;
			for (double a=minAlpha; a<=maxAlpha; a+=incremnt)
			{
				double hlod=0;
				for (unsigned i=0;i<size;i++)
				{
					double l = 1 + a * pow(10,lodvec[i]) - a;
					for (unsigned j=0;j<aj.size();j++) l += aj[j] * ( pow(10,lodoth[j][i]) -1 );
					hlod+=log10(l);
				}
				if (hlod>maxhlod) { maxhlod=hlod; maxa=a; }
			}
			if (incremnt<=0.01) break;
			minAlpha = std::max(minAlpha, maxa - incremnt);
			maxAlpha = std::min(maxAlpha, maxa + incremnt);
		}
		return HLOD_result(LOD,maxa,maxhlod);
	}
	
	HLOD_result cal_HLOD_summation(const std::vector<double>& lodvec, const std::vector< std::vector<double> >& lodoth, std::vector<double>& result_aj, double incr_aj)
	{
		std::vector<double> aj(lodoth.size(),0.0);
		HLOD_result max_r;
		for (;;)
		{
			HLOD_result r = cal_HLOD_summation(lodvec,lodoth,aj);
			if (r.HLOD > max_r.HLOD) { max_r=r; result_aj=aj; }
			double sum=std::accumulate( aj.begin(), aj.end(), 0.0 );
			if (sum + incr_aj > 1)
			{
				bool not_added=true;
				for (unsigned j=0;j<aj.size()-1;j++)
				{
					if (aj[j]!=0)
					{
						aj[j]=0;
						aj[j+1]+=incr_aj;
						not_added=false;
					}
				}
				if (not_added) break;
			}
			else
				aj[0]+=incr_aj;
		}
		return max_r;
	} //*/
	
	double LOD_to_P(const double LOD)
	// Province MA. (2001) The Signiﬁcance of Not Finding a Gen. AJHG 69:660-3
	// careful: this assume df=1, i.e. only theta varied to maximize likelihood. If >1 parameters varied (eg, alpha), use the other function.
	{
		int sign_of_LOD = LOD>=0?1:-1;
		//		return gsl_cdf_ugaussian_Q( sign_of_LOD * sqrt(4.605170186*fabs(LOD)) );	// works fine
		return cdf_norms_q( sign_of_LOD * sqrt(4.605170186*fabs(LOD)) );
	}

	double P_to_LOD(const double P) // reverse of the above
	{
		double LOD = 0;
		if (P<0.5)
		{
			double z = qnorms(P);
			double chi = z*z;
			LOD = chi * 0.217147240951626;
		}
		else if (P>0.5)
		{
			double z = qnorms(P);
			double chi = z*z;
			LOD = chi * -0.217147240951626;
		}
		return LOD;
	}
	
	double LOD_to_P(const double LOD, const unsigned df)
	//http://www.sph.umich.edu/csg/abecasis/LAMP/tour/association.html
	// Convert LOD score to p-value: LOD*2ln(10) = chi_square(df=2)
	{
		if (LOD<0) exit_error("This LOD_to_P function can only handle positive LOD. For negative ones please use another function.");
		//		return gsl_cdf_chisq_P(LOD*4.605170186,df);		// works fine
		return cdf_chisq_p(LOD*4.605170186,df);
	}
	
	NPL_result::NPL_result(double d, double c, double z):delta(d),chisq(c),Zmean(z)
	{
		lod_chisq  = delta>=0 ? chisq*0.2171472409516				:	-chisq*0.2171472409516;
//		pval_chisq = delta>=0 ? (1-gsl_cdf_chisq_P(chisq,1))/2		:	1-(1-gsl_cdf_chisq_P(chisq,1))/2;		// works fine
//		pval_Zmean = delta>=0 ? gsl_cdf_ugaussian_Q(fabs(Zmean))	:	1-gsl_cdf_ugaussian_Q(fabs(Zmean));	}	// works fine
		pval_chisq = delta>=0 ? (1-cdf_chisq_1df_p(chisq))/2		:	1-(1-cdf_chisq_1df_p(chisq))/2;
		pval_Zmean = delta>=0 ? cdf_norms_q(fabs(Zmean))			:	1-cdf_norms_q(fabs(Zmean));
	}
	
	NPL_result cal_NPL(const std::vector<double>& Z_scores, double minDelta, double maxDelta) // return <delta,chisq>
	{
		if (minDelta == maxDelta) return NPL_result();
		unsigned size=Z_scores.size();
		if (size==0) return NPL_result();
		maxDelta+=std::numeric_limits<double>::epsilon()*20;
		minDelta-=std::numeric_limits<double>::epsilon()*20;
		double sumZ = 0.0;
		for (unsigned i=0;i<size;i++) sumZ+=Z_scores[i];
		double Zmean=sumZ/sqrt(size);
		double maxf=-DBL_MAX;	// max sum
		double maxb=0;			// delta value to obtain maxf
		for (;;)
		{
			double incremnt = (maxDelta-minDelta)/20;
			for (double delta=minDelta; delta<=maxDelta; delta+=incremnt)
			{
				double sum = 0.0;
				for (unsigned i=0;i<size;i++)
					sum += 2.0 * log(1.0 + delta * Z_scores[i]);
				if (sum>maxf) { maxf=sum; maxb=delta; }
			}
			if (incremnt<=0.00001) break;
			minDelta = std::max(minDelta, maxb - incremnt);
			maxDelta = std::min(maxDelta, maxb + incremnt);
		}
		return NPL_result(maxb,maxf,Zmean);
	}
	
	NPL_result cal_NPL(const std::multiset<double>& Z_scores, double minDelta, double maxDelta) // return <delta,chisq>
	{
		if (minDelta == maxDelta) return NPL_result();
		unsigned size=Z_scores.size();
		if (size==0) return NPL_result();
		maxDelta+=std::numeric_limits<double>::epsilon()*20;
		minDelta-=std::numeric_limits<double>::epsilon()*20;
		double sumZ = 0.0;
		for (each_element_const(Z_scores,it)) sumZ+=*it;
		double Zmean=sumZ/sqrt(size);
		double maxf=-DBL_MAX;	// max sum
		double maxb=0;			// delta value to obtain maxf
		for (;;)
		{
			double incremnt = (maxDelta-minDelta)/20;
			for (double delta=minDelta; delta<=maxDelta; delta+=incremnt)
			{
				double sum = 0.0;
				for (each_element_const(Z_scores,it))
					sum += 2.0 * log(1.0 + delta * *it);
				if (sum>maxf) { maxf=sum; maxb=delta; }
			}
			if (incremnt<=0.00001) break;
			minDelta = std::max(minDelta, maxb - incremnt);
			maxDelta = std::min(maxDelta, maxb + incremnt);
		}
		return NPL_result(maxb,maxf,Zmean);
	}
	
	double kosambi2rf(double k)
	{
		double morgan,r;
		if (k<0) exit_error("Centimorgan cannot be negative"); // return(-10000);
		morgan=k/100;
		r=0.5*(exp(4*morgan)-1)/(exp(4*morgan)+1);
		return r;
	}
	
	double haldane2rf(double h)
	{
		double morgan,r;
		if (h<0) exit_error("Centimorgan cannot be negative"); // return(-10000);
		morgan=(h)/100;
		r=0.5*(1-exp(-2*morgan));
		return r;
	}

	// ------------------ linkage disequilibrium ------------------

	LDmodel::LDmodel():is_set(false) {}
	
	void LDmodel::set(double p, double q, double Dp) // p for A2, q for B2
	{
		is_set = true;
		p1 = 1-p;
		p2 = p;
		q1 = 1-q;
		q2 = q;
		set(Dp);
	}

	void LDmodel::set(double Dp)
	{
		if (!is_set) exit_error("LDmodel not set.");
		Dprime = Dp;
		if (Dprime<0)	Dmax = std::min(p1*q1,p2*q2);
		else			Dmax = std::min(p1*q2,p2*q1);
		D = Dprime * Dmax;
		X11 = p1*q1 + D;
		X12 = p1*q2 - D;
		X21 = p2*q1 - D;
		X22 = p2*q2 + D;
	}

	double LDmodel::cond_p(int B)
	{
		if (!is_set) exit_error("LDmodel not set.");
		if (B)	return X22;
		else	return X21;
	}

	
	template <typename T1,typename T2>
	void cal_LD(T1 x11, T1 x12, T1 x21, T1 x22, T2& Dprime, T2& Rsquare)
	// http://en.wikipedia.org/wiki/Linkage_disequilibrium
	{
		double p1,p2,q1,q2,D;
		p1 = x11 + x12;
		p2 = x21 + x22;
		q1 = x11 + x21;
		q2 = x12 + x22;
		if (p1==0 || p2==0 || q1==0 || q2==0) { Dprime=0; Rsquare=0; return; } // when a column or row is 0, D cannot be calculated
		D=x11-p1*q1;
		if (D>=0)	Dprime = D /  std::min(p1*q2,p2*q1); // D/Dmax
		else		Dprime = D / -std::min(p1*q1,p2*q2); // D/Dmin
		Rsquare=D*D/p1/p2/q1/q2;
	}
	template void cal_LD<unsigned,double>(unsigned x11,unsigned x12, unsigned x21, unsigned x22, double& Dprime, double& Rsquare);
	template void cal_LD<int,double>(int x11,int x12, int x21, int x22, double& Dprime, double& Rsquare);
	template void cal_LD<double,double>(double x11,double x12, double x21, double x22, double& Dprime, double& Rsquare);

	template <typename T1,typename T2>
	void cal_LD_from_hap_counts(T1 x11,T1 x12, T1 x21, T1 x22, T2& Dprime, T2& Rsquare)
	{
		double t=x11+x12+x21+x22;
		if (t)	cal_LD(x11/t,x12/t,x21/t,x22/t,Dprime,Rsquare);
		else {	Dprime=0; Rsquare=0; }
	}	
	template void cal_LD_from_hap_counts<unsigned,double>(unsigned x11,unsigned x12, unsigned x21, unsigned x22, double& Dprime, double& Rsquare);
	template void cal_LD_from_hap_counts<int,double>(int x11,int x12, int x21, int x22, double& Dprime, double& Rsquare);
	template void cal_LD_from_hap_counts<double,double>(double x11,double x12, double x21, double x22, double& Dprime, double& Rsquare);
	
	template <typename T>
	void cubex (const T inputarray[3][3], map<string,double>& result)
	// http://www.oege.org/software/cubex/ modified so that result["Dprime"] & result["Rsquared"] is the most probable solution.
	// When no solution exists, result["Dprime"] & result["Rsquared"] doesn't exist. Access to them will return 0, which are good default values.
	// Therefore, no need to check their existence. This happens when at least one of the SNPs is non-polymorphic.
	// The most probable solution = alpha/beta/gamma that has the lowest xxxx-chisq (min_chisq). But sig chisq means genotype out of HWE!
	// Not all alpha/beta/gamma solution exist: it exist only when xxxxposs=1
	// I also changed result["xxxx-e####"] to a double variable for faster computation.
	// result["Error"]=1 seems never will happen, from my simulations.
	{
		result.clear();
		const T& n1111 = inputarray[0][0];
		const T& n1112 = inputarray[0][1];
		const T& n1122 = inputarray[0][2];
		const T& n1211 = inputarray[1][0];
		const T& n1212 = inputarray[1][1];
		const T& n1222 = inputarray[1][2];
		const T& n2211 = inputarray[2][0];
		const T& n2212 = inputarray[2][1];
		const T& n2222 = inputarray[2][2];
		
		// calcualte df for x2_good_fit_HWE, based on http://www.oege.org/software/cubex/
		int df=5;
		if (n1111+n1112+n1122==0 || n2211+n2212+n2222==0 || n1111+n1211+n2211==0 || n1122+n1222+n2222==0) df=6;
		for (unsigned r=0;r<3;r++)
			for (unsigned c=0;c<3;c++)
				if (inputarray[r][c]==0) df--;
		//if (df<=0) df=1; // This is what I guess
		result["df"]=df;
		
		double n = (n1111 + n1112 + n1122 + n1211 + n1212 + n1222 + n2211 + n2212 + n2222); // prv int
		result["n"] = n;
		double p = ((n1111+n1112+n1122)*2.0+(n1211+n1212+n1222))/(2.0 * n);
		double q = ((n1111+n1211+n2211)*2.0+(n1112+n1212+n2212))/(2.0 * n);
		double n11 = (2*n1111 + n1112 + n1211); // prv int
		//	int n12 = (2*n1122 + n1112 + n1222);
		//	int n21 = (2*n2211 + n2212 + n1211);
		//	int n22 = (2*n2222 + n2212 + n1222);
		double a0 = -n11*p*q;
		double a1 = -n11*(1.0 - 2.0*p - 2.0*q) - n1212*(1.0 - p - q) + 2.0*n*p*q;
		double a2 = 2.0*n*(1.0 - 2.0*p - 2.0*q) - 2.0*n11 - n1212;
		double a3 = 4.0 * n;
		double minhap = n11 / (2.0 * n);
		double maxhap = (n11 + n1212) / (2.0 * n);
		result["minhap"] = minhap;
		result["maxhap"] = maxhap;
		if ( p < 1.0 && p > 0.0 )
			result["hwchisnp1"] = \
			pow(n1111 + n1112 + n1122 - p*p*n , 2) / (p*p*n) + \
			pow(n1211 + n1212 + n1222 - 2*p*(1-p)*n , 2) / (2*p*(1-p)*n) + \
			pow(n2211 + n2212 + n2222 - (1-p)*(1-p)*n , 2) / ((1-p)*(1-p)*n);
		else
			result["hwchisnp1"] = 0.0;
		if ( q < 1.0 && q > 0.0 )
			result["hwchisnp2"] = \
			pow(n1111 + n1211 + n2211 - q*q*n , 2) / (q*q*n) + \
			pow(n1112 + n1212 + n2212 - 2*q*(1-q)*n , 2) / (2 * q * (1-q)*n) + \
			pow(n1122 + n1222 + n2222 - (1-q)*(1-q)*n , 2) / ((1-q)*(1-q)*n);
		else
			result["hwchisnp2"] = 0.0;
		
		double a = a3;
		double b = a2;
		double c = a1;
		double dee = a0;
		double xN = -b/(3.0*a);
		double d2 = (pow(b,2)-3.0*a*c)/(9*pow(a,2));
		double yN = a * pow(xN,3) + b * pow(xN,2) + c * xN + dee;
		double yN2 = pow(yN,2);
		double h2 = 4 * pow(a,2) * pow(d2,3);
		double f11,f12,f21,f22;
		double chisq1111,chisq1112,chisq1122,chisq1211,chisq1212,chisq1222,chisq2211,chisq2212,chisq2222;
		double alpha_e1111,alpha_e1112,alpha_e1122,alpha_e1211,alpha_e1212,alpha_e1222,alpha_e2211,alpha_e2212,alpha_e2222;
		double beta_e1111,beta_e1112,beta_e1122,beta_e1211,beta_e1212,beta_e1222,beta_e2211,beta_e2212,beta_e2222;
		double gamma_e1111,gamma_e1112,gamma_e1122,gamma_e1211,gamma_e1212,gamma_e1222,gamma_e2211,gamma_e2212,gamma_e2222;
		double D, Dmax, Dprime, rsquared;
		double min_chisq=DBL_MAX;
		result["realnoproblem"] = 0;
		
		if ( fabs(yN2-h2) <= 0.00001 )
			result["realnoproblem"] = 1;
		result["Error"] = 0;
		if ( yN2 > h2 ) // option 1
		{
			double number1 = 0.0;
			double number2 = 0.0;
			if ( (1.0/(2.0*a)*(-yN + pow((yN2 - h2),0.5))) < 0 )
				number1 = -pow(-(1.0/(2.0*a)*(-yN + pow((yN2 - h2),0.5))),1.0/3.0);
			else
				number1 = pow((1.0/(2.0*a)*(-yN + pow((yN2 - h2),0.5))),1.0/3.0);
			if ( (1.0/(2.0*a)*(-yN - pow((yN2 - h2),0.5))) < 0 )
				number2 = -pow(-(1.0/(2.0*a)*(-yN - pow((yN2 - h2),0.5))),1.0/3.0);
			else
				number2 = pow((1.0/(2.0*a)*(-yN - pow((yN2 - h2),0.5))),1.0/3.0);
			double alpha = xN + number1 + number2;
			result["alpha"] = alpha;
			result["beta"] = -9999999; // "Not a real root";
			result["gamma"] = -9999999; // "Not a real root";
			result["p"] = p;
			result["q"] = q;
			result["pq"] = p * q;
			if ( result["alpha"] >= minhap - 0.00001 && result["alpha"] <= maxhap + 0.00001 )
			{
				result["alphaposs"] = 1;
				f11 = result["alpha"];
				f12 = p - f11;
				f21 = q - f11;
				f22 = 1 - (f11 + f12 + f21);
				D = (f11 * f22) - (f12 * f21);
				if ( D >= 0.0 )	Dmax = min(p * (1.0-q), q * (1.0-p));
				else			Dmax = min(p*q,(1-p)*(1-q));
				Dprime = D / Dmax;
				rsquared = (D*D) / (p * (1-p) * q * (1-q));
				result["alphaDprime"] = Dprime;
				result["alpharsquared"] = rsquared;
				result["Dprime"] = Dprime;
				result["Rsquared"] = rsquared;
				result["alphaf11"] = f11;
				result["alphaf12"] = f12;
				result["alphaf21"] = f21;
				result["alphaf22"] = f22;
				if ( min(min(f11,f12),min(f21,f22)) < -0.0000001 )
					result["alphaposs"] = 0;
				else
				{
					alpha_e1111 = n * f11*f11;
					alpha_e1112 = 2 * n * f11 * f12;
					alpha_e1122 = n * f12*f12;
					alpha_e1211 = 2 * n * f11 * f21;
					alpha_e1212 = 2 * n * f12 * f21 + 2 * n * f11 * f22 ;
					alpha_e1222 = 2 * n * f12 * f22;
					alpha_e2211 = n * f21*f21;
					alpha_e2212 = 2 * n * f21 * f22;
					alpha_e2222 = n * f22*f22;
					if ( alpha_e1111 > 0.0 ) chisq1111 = pow(n1111 - alpha_e1111,2)/alpha_e1111;	else chisq1111 = 0.0;
					if ( alpha_e1112 > 0.0 ) chisq1112 = pow(n1112 - alpha_e1112,2)/alpha_e1112;	else chisq1112 = 0.0;
					if ( alpha_e1122 > 0.0 ) chisq1122 = pow(n1122 - alpha_e1122,2)/alpha_e1122;	else chisq1122 = 0.0;
					if ( alpha_e1211 > 0.0 ) chisq1211 = pow(n1211 - alpha_e1211,2)/alpha_e1211;	else chisq1211 = 0.0;
					if ( alpha_e1212 > 0.0 ) chisq1212 = pow(n1212 - alpha_e1212,2)/alpha_e1212;	else chisq1212 = 0.0;
					if ( alpha_e1222 > 0.0 ) chisq1222 = pow(n1222 - alpha_e1222,2)/alpha_e1222;	else chisq1222 = 0.0;
					if ( alpha_e2211 > 0.0 ) chisq2211 = pow(n2211 - alpha_e2211,2)/alpha_e2211;	else chisq2211 = 0.0;
					if ( alpha_e2212 > 0.0 ) chisq2212 = pow(n2212 - alpha_e2212,2)/alpha_e2212;	else chisq2212 = 0.0;
					if ( alpha_e2222 > 0.0 ) chisq2222 = pow(n2222 - alpha_e2222,2)/alpha_e2222;	else chisq2222 = 0.0;
					result["alpha-chisq"] = (chisq1111+chisq1112+chisq1122+chisq1211+chisq1212+chisq1222+chisq2211+chisq2212+chisq2222);
					min_chisq=result["alpha-chisq"];
					result["min_chisq"]=min_chisq;
				}
			}
			else
				result["alphaposs"] = 0;
			result["betaposs"] = 0;
			result["gammaposs"] = 0;
		}
		else if ( yN2 == h2 ) // option 2
		{
			double delta = pow((yN/2.0*a),(1.0/3.0));
			result["alpha"] = xN + delta;
			result["beta"] = xN + delta;
			result["gamma"] = xN - 2.0*delta;
			result["p"] = p;
			result["q"] = q;
			result["pq"] = p * q;
			if ( result["alpha"] >= minhap - 0.00001 && result["alpha"] <= maxhap + 0.00001 )
			{
				result["alphaposs"] = 1;
				f11 = result["alpha"];
				f12 = p - f11;
				f21 = q - f11;
				f22 = 1 - (f11 + f12 + f21);
				D = (f11 * f22) - (f12 * f21);
				if ( D >= 0.0 )	Dmax = min(p * (1-q), q * (1-p));
				else			Dmax = min(p * q    , (1-p)*(1-q));
				Dprime = D / Dmax;
				rsquared = (D * D) / (p * (1-p) * q * (1-q));
				result["alphaDprime"] = Dprime;
				result["alpharsquared"] = rsquared;
				result["alphaf11"] = f11;
				result["alphaf12"] = f12;
				result["alphaf21"] = f21;
				result["alphaf22"] = f22;
				if ( min(min(f11,f12),min(f21,f22)) < -0.0000001 )
					result["alphaposs"] = 0;
				else
				{
					alpha_e1111 = n * f11*f11;
					alpha_e1112 = 2 * n * f11 * f12;
					alpha_e1122 = n * f12*f12;
					alpha_e1211 = 2 * n * f11 * f21;
					alpha_e1212 = 2 * n * f12 * f21 + 2 * n * f11 * f22 ;
					alpha_e1222 = 2 * n * f12 * f22;
					alpha_e2211 = n * f21*f21;
					alpha_e2212 = 2 * n * f21 * f22;
					alpha_e2222 = n * f22*f22;
					if ( alpha_e1111 > 0.0 ) chisq1111 = pow(n1111 - alpha_e1111,2)/alpha_e1111;	else chisq1111 = 0.0;
					if ( alpha_e1112 > 0.0 ) chisq1112 = pow(n1112 - alpha_e1112,2)/alpha_e1112;	else chisq1112 = 0.0;
					if ( alpha_e1122 > 0.0 ) chisq1122 = pow(n1122 - alpha_e1122,2)/alpha_e1122;	else chisq1122 = 0.0;
					if ( alpha_e1211 > 0.0 ) chisq1211 = pow(n1211 - alpha_e1211,2)/alpha_e1211;	else chisq1211 = 0.0;
					if ( alpha_e1212 > 0.0 ) chisq1212 = pow(n1212 - alpha_e1212,2)/alpha_e1212;	else chisq1212 = 0.0;
					if ( alpha_e1222 > 0.0 ) chisq1222 = pow(n1222 - alpha_e1222,2)/alpha_e1222;	else chisq1222 = 0.0;
					if ( alpha_e2211 > 0.0 ) chisq2211 = pow(n2211 - alpha_e2211,2)/alpha_e2211;	else chisq2211 = 0.0;
					if ( alpha_e2212 > 0.0 ) chisq2212 = pow(n2212 - alpha_e2212,2)/alpha_e2212;	else chisq2212 = 0.0;
					if ( alpha_e2222 > 0.0 ) chisq2222 = pow(n2222 - alpha_e2222,2)/alpha_e2222;	else chisq2222 = 0.0;
					result["alpha-chisq"] = (chisq1111+chisq1112+chisq1122+chisq1211+chisq1212+chisq1222+chisq2211+chisq2212+chisq2222);
					if (result["alpha-chisq"] < min_chisq)
					{
						min_chisq=result["alpha-chisq"];
						result["min_chisq"]=min_chisq;
						result["Dprime"] = Dprime;
						result["Rsquared"] = rsquared;
					}
				}
			}
			else
				result["alphaposs"] = 0;
			if ( result["beta"] >= minhap - 0.00001 && result["beta"] <= maxhap + 0.00001 )
			{
				result["betaposs"] = 1;
				f11 = result["beta"];
				f12 = p - f11;
				f21 = q - f11;
				f22 = 1 - (f11 + f12 + f21);
				D = (f11 * f22) - (f12 * f21);
				if ( D >= 0.0 )	Dmax = min(p * (1.0-q), q * (1.0-p));
				else			Dmax = min(p*q,(1-p)*(1-q));
				Dprime = D / Dmax;
				rsquared = (D * D) / (p * (1-p) * q * (1-q));
				result["betaDprime"] = Dprime;
				result["betarsquared"] = rsquared;
				result["betaf11"] = f11;
				result["betaf12"] = f12;
				result["betaf21"] = f21;
				result["betaf22"] = f22;
				if ( min(min(f11,f12),min(f21,f22)) < -0.0000001 )
					result["betaposs"] = 0;
				else
				{
					beta_e1111 = n * f11*f11;
					beta_e1112 = 2 * n * f11 * f12;
					beta_e1122 = n * f12*f12;
					beta_e1211 = 2 * n * f11 * f21;
					beta_e1212 = 2 * n * f12 * f21 + 2 * n * f11 * f22 ;
					beta_e1222 = 2 * n * f12 * f22;
					beta_e2211 = n * f21*f21;
					beta_e2212 = 2 * n * f21 * f22;
					beta_e2222 = n * f22*f22;
					if ( beta_e1111 > 0.0 )	chisq1111 = pow(n1111 - beta_e1111,2)/beta_e1111;	else chisq1111 = 0.0;
					if ( beta_e1112 > 0.0 )	chisq1112 = pow(n1112 - beta_e1112,2)/beta_e1112;	else chisq1112 = 0.0;
					if ( beta_e1122 > 0.0 )	chisq1122 = pow(n1122 - beta_e1122,2)/beta_e1122;	else chisq1122 = 0.0;
					if ( beta_e1211 > 0.0 )	chisq1211 = pow(n1211 - beta_e1211,2)/beta_e1211;	else chisq1211 = 0.0;
					if ( beta_e1212 > 0.0 )	chisq1212 = pow(n1212 - beta_e1212,2)/beta_e1212;	else chisq1212 = 0.0;
					if ( beta_e1222 > 0.0 )	chisq1222 = pow(n1222 - beta_e1222,2)/beta_e1222;	else chisq1222 = 0.0;
					if ( beta_e2211 > 0.0 )	chisq2211 = pow(n2211 - beta_e2211,2)/beta_e2211;	else chisq2211 = 0.0;
					if ( beta_e2212 > 0.0 )	chisq2212 = pow(n2212 - beta_e2212,2)/beta_e2212;	else chisq2212 = 0.0;
					if ( beta_e2222 > 0.0 )	chisq2222 = pow(n2222 - beta_e2222,2)/beta_e2222;	else chisq2222 = 0.0;
					result["beta-chisq"] = (chisq1111+chisq1112+chisq1122+chisq1211+chisq1212+chisq1222+chisq2211+chisq2212+chisq2222);
					if (result["beta-chisq"] < min_chisq)
					{
						min_chisq=result["beta-chisq"];
						result["min_chisq"]=min_chisq;
						result["Dprime"] = Dprime;
						result["Rsquared"] = rsquared;
					}
				}
			}
			else
				result["betaposs"] = 0;
			if ( result["gamma"] >= minhap - 0.00001 && result["gamma"] <= maxhap + 0.00001 )
			{
				result["gammaposs"] = 1;
				f11 = result["gamma"];
				f12 = p - f11;
				f21 = q - f11;
				f22 = 1 - (f11 + f12 + f21);
				D = (f11 * f22) - (f12 * f21);
				if ( D >= 0.0 )	Dmax = min(p * (1.0-q), q * (1.0-p));
				else			Dmax = min(p*q,(1-p)*(1-q));
				Dprime = D / Dmax;
				rsquared = (D * D) / (p * (1-p) * q * (1-q));
				result["gammaDprime"] = Dprime;
				result["gammarsquared"] = rsquared;
				result["gammaf11"] = f11;
				result["gammaf12"] = f12;
				result["gammaf21"] = f21;
				result["gammaf22"] = f22;
				if ( min(min(f11,f12),min(f21,f22)) < -0.0000001 )
					result["gammaposs"] = 0;
				else
				{
					gamma_e1111 = n * f11*f11;
					gamma_e1112 = 2 * n * f11 * f12;
					gamma_e1122 = n * f12*f12;
					gamma_e1211 = 2 * n * f11 * f21;
					gamma_e1212 = 2 * n * f12 * f21 + 2 * n * f11 * f22 ;
					gamma_e1222 = 2 * n * f12 * f22;
					gamma_e2211 = n * f21*f21;
					gamma_e2212 = 2 * n * f21 * f22;
					gamma_e2222 = n * f22*f22;
					if ( gamma_e1111 > 0.0 ) chisq1111 = pow(n1111 - gamma_e1111,2)/gamma_e1111;	else chisq1111 = 0.0;
					if ( gamma_e1112 > 0.0 ) chisq1112 = pow(n1112 - gamma_e1112,2)/gamma_e1112;	else chisq1112 = 0.0;
					if ( gamma_e1122 > 0.0 ) chisq1122 = pow(n1122 - gamma_e1122,2)/gamma_e1122;	else chisq1122 = 0.0;
					if ( gamma_e1211 > 0.0 ) chisq1211 = pow(n1211 - gamma_e1211,2)/gamma_e1211;	else chisq1211 = 0.0;
					if ( gamma_e1212 > 0.0 ) chisq1212 = pow(n1212 - gamma_e1212,2)/gamma_e1212;	else chisq1212 = 0.0;
					if ( gamma_e1222 > 0.0 ) chisq1222 = pow(n1222 - gamma_e1222,2)/gamma_e1222;	else chisq1222 = 0.0;
					if ( gamma_e2211 > 0.0 ) chisq2211 = pow(n2211 - gamma_e2211,2)/gamma_e2211;	else chisq2211 = 0.0;
					if ( gamma_e2212 > 0.0 ) chisq2212 = pow(n2212 - gamma_e2212,2)/gamma_e2212;	else chisq2212 = 0.0;
					if ( gamma_e2222 > 0.0 ) chisq2222 = pow(n2222 - gamma_e2222,2)/gamma_e2222;	else chisq2222 = 0.0;
					result["gamma-chisq"] = (chisq1111+chisq1112+chisq1122+chisq1211+chisq1212+chisq1222+chisq2211+chisq2212+chisq2222);
					if (result["gamma-chisq"] < min_chisq)
					{
						min_chisq=result["gamma-chisq"];
						result["min_chisq"]=min_chisq;
						result["Dprime"] = Dprime;
						result["Rsquared"] = rsquared;
					}
				}
			}
			else
				result["gammaposs"] = 0;
		}
		else if ( yN2 < h2 ) // option 3
		{
			double h = pow(h2, 0.5);
			double theta = ((acos(-yN/h))/3.0);
			// delta = pow((yN/2.0*a),(1.0/3.0)) # is it correct to reuse this?
			double delta = pow(d2,0.5);
			result["alpha"] = xN + 2.0 * delta * cos(theta);
			result["beta"] = xN + 2.0 * delta * cos(2.0 * M_PI/3.0 + theta);
			result["gamma"] = xN + 2.0 * delta * cos(4.0 * M_PI/3.0 + theta);
			result["p"] = p;
			result["q"] = q;
			result["pq"] = p * q;
			if ( result["alpha"] >= minhap - 0.00001 && result["alpha"] <= maxhap + 0.00001 )
			{
				result["alphaposs"] = 1;
				f11 = result["alpha"];
				f12 = p - f11;
				f21 = q - f11;
				f22 = 1 - (f11 + f12 + f21);
				D = (f11 * f22) - (f12 * f21);
				if ( D >= 0.0 )	Dmax = min(p * (1.0-q), q * (1.0-p));
				else			Dmax = min(p*q,(1-p)*(1-q));
				Dprime = D / Dmax;
				rsquared = (D * D) / (p * (1-p) * q * (1-q));
				result["alphaDprime"] = Dprime;
				result["alpharsquared"] = rsquared;
				result["alphaf11"] = f11;
				result["alphaf12"] = f12;
				result["alphaf21"] = f21;
				result["alphaf22"] = f22;
				if ( min(min(f11,f12),min(f21,f22)) < -0.0000001 )
					result["alphaposs"] = 0;
				else
				{
					alpha_e1111 = n * f11*f11;
					alpha_e1112 = 2 * n * f11 * f12;
					alpha_e1122 = n * f12*f12;
					alpha_e1211 = 2 * n * f11 * f21;
					alpha_e1212 = 2 * n * f12 * f21 + 2 * n * f11 * f22 ;
					alpha_e1222 = 2 * n * f12 * f22;
					alpha_e2211 = n * f21*f21;
					alpha_e2212 = 2 * n * f21 * f22;
					alpha_e2222 = n * f22*f22;
					if ( alpha_e1111 > 0.0 ) chisq1111 = pow(n1111 - alpha_e1111,2)/alpha_e1111;	else chisq1111 = 0.0;
					if ( alpha_e1112 > 0.0 ) chisq1112 = pow(n1112 - alpha_e1112,2)/alpha_e1112;	else chisq1112 = 0.0;
					if ( alpha_e1122 > 0.0 ) chisq1122 = pow(n1122 - alpha_e1122,2)/alpha_e1122;	else chisq1122 = 0.0;
					if ( alpha_e1211 > 0.0 ) chisq1211 = pow(n1211 - alpha_e1211,2)/alpha_e1211;	else chisq1211 = 0.0;
					if ( alpha_e1212 > 0.0 ) chisq1212 = pow(n1212 - alpha_e1212,2)/alpha_e1212;	else chisq1212 = 0.0;
					if ( alpha_e1222 > 0.0 ) chisq1222 = pow(n1222 - alpha_e1222,2)/alpha_e1222;	else chisq1222 = 0.0;
					if ( alpha_e2211 > 0.0 ) chisq2211 = pow(n2211 - alpha_e2211,2)/alpha_e2211;	else chisq2211 = 0.0;
					if ( alpha_e2212 > 0.0 ) chisq2212 = pow(n2212 - alpha_e2212,2)/alpha_e2212;	else chisq2212 = 0.0;
					if ( alpha_e2222 > 0.0 ) chisq2222 = pow(n2222 - alpha_e2222,2)/alpha_e2222;	else chisq2222 = 0.0;
					result["alpha-chisq"] = (chisq1111+chisq1112+chisq1122+chisq1211+chisq1212+chisq1222+chisq2211+chisq2212+chisq2222);
					if (result["alpha-chisq"] < min_chisq)
					{
						min_chisq=result["alpha-chisq"];
						result["min_chisq"]=min_chisq;
						result["Dprime"] = Dprime;
						result["Rsquared"] = rsquared;
					}
				}
			}
			else result["alphaposs"] = 0;
			if ( result["beta"] >= minhap - 0.00001 && result["beta"] <= maxhap + 0.00001 )
			{
				result["betaposs"] = 1;
				f11 = result["beta"];
				f12 = p - f11;
				f21 = q - f11;
				f22 = 1 - (f11 + f12 + f21);
				D = (f11 * f22) - (f12 * f21);
				if ( D >= 0.0 )	Dmax = min(p * (1.0-q), q * (1.0-p));
				else			Dmax = min(p*q,(1-p)*(1-q));
				Dprime = D / Dmax;
				rsquared = (D * D) / (p * (1-p) * q * (1-q));
				result["betaDprime"] = Dprime;
				result["betarsquared"] = rsquared;
				result["betaf11"] = f11;
				result["betaf12"] = f12;
				result["betaf21"] = f21;
				result["betaf22"] = f22;
				if ( min(min(f11,f12),min(f21,f22)) < -0.0000001 )
					result["betaposs"] = 0;
				else
				{
					beta_e1111 = n * f11*f11;
					beta_e1112 = 2 * n * f11 * f12;
					beta_e1122 = n * f12*f12;
					beta_e1211 = 2 * n * f11 * f21;
					beta_e1212 = 2 * n * f12 * f21 + 2 * n * f11 * f22;
					beta_e1222 = 2 * n * f12 * f22;
					beta_e2211 = n * f21*f21;
					beta_e2212 = 2 * n * f21 * f22;
					beta_e2222 = n * f22*f22;
					if ( beta_e1111 > 0.0 ) chisq1111 = pow(n1111 - beta_e1111,2)/beta_e1111;	else chisq1111 = 0.0;
					if ( beta_e1112 > 0.0 ) chisq1112 = pow(n1112 - beta_e1112,2)/beta_e1112;	else chisq1112 = 0.0;
					if ( beta_e1122 > 0.0 ) chisq1122 = pow(n1122 - beta_e1122,2)/beta_e1122;	else chisq1122 = 0.0;
					if ( beta_e1211 > 0.0 ) chisq1211 = pow(n1211 - beta_e1211,2)/beta_e1211;	else chisq1211 = 0.0;
					if ( beta_e1212 > 0.0 ) chisq1212 = pow(n1212 - beta_e1212,2)/beta_e1212;	else chisq1212 = 0.0;
					if ( beta_e1222 > 0.0 ) chisq1222 = pow(n1222 - beta_e1222,2)/beta_e1222;	else chisq1222 = 0.0;
					if ( beta_e2211 > 0.0 ) chisq2211 = pow(n2211 - beta_e2211,2)/beta_e2211;	else chisq2211 = 0.0;
					if ( beta_e2212 > 0.0 ) chisq2212 = pow(n2212 - beta_e2212,2)/beta_e2212;	else chisq2212 = 0.0;
					if ( beta_e2222 > 0.0 ) chisq2222 = pow(n2222 - beta_e2222,2)/beta_e2222;	else chisq2222 = 0.0;
					result["beta-chisq"] = (chisq1111+chisq1112+chisq1122+chisq1211+chisq1212+chisq1222+chisq2211+chisq2212+chisq2222);
					if (result["beta-chisq"] < min_chisq)
					{
						min_chisq=result["beta-chisq"];
						result["min_chisq"]=min_chisq;
						result["Dprime"] = Dprime;
						result["Rsquared"] = rsquared;
					}
				}
			}
			else result["betaposs"] = 0;
			if ( result["gamma"] >= minhap - 0.00001 && result["gamma"] <= maxhap + 0.00001 )
			{
				result["gammaposs"] = 1;
				f11 = result["gamma"];
				f12 = p - f11;
				f21 = q - f11;
				f22 = 1 - (f11 + f12 + f21);
				D = (f11 * f22) - (f12 * f21);
				if ( D >= 0.0 )	Dmax = min(p * (1.0-q), q * (1.0-p));
				else			Dmax = min(p*q,(1-p)*(1-q));
				Dprime = D / Dmax;
				rsquared = (D * D) / (p * (1-p) * q * (1-q));
				result["gammaDprime"] = Dprime;
				result["gammarsquared"] = rsquared;
				result["gammaf11"] = f11;
				result["gammaf12"] = f12;
				result["gammaf21"] = f21;
				result["gammaf22"] = f22;
				if ( min(min(f11,f12),min(f21,f22)) < -0.0000001 )
					result["gammaposs"] = 0;
				else
				{
					gamma_e1111 = n * f11*f11;
					gamma_e1112 = 2 * n * f11 * f12;
					gamma_e1122 = n * f12*f12;
					gamma_e1211 = 2 * n * f11 * f21;
					gamma_e1212 = 2 * n * f12 * f21 + 2 * n * f11 * f22 ;
					gamma_e1222 = 2 * n * f12 * f22;
					gamma_e2211 = n * f21*f21;
					gamma_e2212 = 2 * n * f21 * f22;
					gamma_e2222 = n * f22*f22;
					if ( gamma_e1111 > 0.0 ) chisq1111 = pow(n1111 - gamma_e1111,2)/gamma_e1111;	else chisq1111 = 0.0;
					if ( gamma_e1112 > 0.0 ) chisq1112 = pow(n1112 - gamma_e1112,2)/gamma_e1112;	else chisq1112 = 0.0;
					if ( gamma_e1122 > 0.0 ) chisq1122 = pow(n1122 - gamma_e1122,2)/gamma_e1122;	else chisq1122 = 0.0;
					if ( gamma_e1211 > 0.0 ) chisq1211 = pow(n1211 - gamma_e1211,2)/gamma_e1211;	else chisq1211 = 0.0;
					if ( gamma_e1212 > 0.0 ) chisq1212 = pow(n1212 - gamma_e1212,2)/gamma_e1212;	else chisq1212 = 0.0;
					if ( gamma_e1222 > 0.0 ) chisq1222 = pow(n1222 - gamma_e1222,2)/gamma_e1222;	else chisq1222 = 0.0;
					if ( gamma_e2211 > 0.0 ) chisq2211 = pow(n2211 - gamma_e2211,2)/gamma_e2211;	else chisq2211 = 0.0;
					if ( gamma_e2212 > 0.0 ) chisq2212 = pow(n2212 - gamma_e2212,2)/gamma_e2212;	else chisq2212 = 0.0;
					if ( gamma_e2222 > 0.0 ) chisq2222 = pow(n2222 - gamma_e2222,2)/gamma_e2222;	else chisq2222 = 0.0;
					result["gamma-chisq"] = (chisq1111+chisq1112+chisq1122+chisq1211+chisq1212+chisq1222+chisq2211+chisq2212+chisq2222);
					if (result["gamma-chisq"] < min_chisq)
					{
						min_chisq=result["gamma-chisq"];
						result["min_chisq"]=min_chisq;
						result["Dprime"] = Dprime;
						result["Rsquared"] = rsquared;
					}
				}
			}
			else result["gammaposs"] = 0;
		}
		else result["Error"] = 1;
	}
	template void cubex<unsigned> (const unsigned inputarray[3][3], map<string,double>& result);
	template void cubex<int> (const int inputarray[3][3], map<string,double>& result);
	template void cubex<double> (const double inputarray[3][3], map<string,double>& result);
	
	template <typename T>
	void cubex (const T inputarray[9], map<string,double>& result)
	{
		T gtp[3][3];
		memcpy(gtp,inputarray,sizeof(T)*9);
		cubex (gtp, result);
	}
	template void cubex<unsigned> (const unsigned inputarray[9], map<string,double>& result);
	template void cubex<int> (const int inputarray[9], map<string,double>& result);
	template void cubex<double> (const double inputarray[9], map<string,double>& result);

	// ------------------ end ------------------

}
