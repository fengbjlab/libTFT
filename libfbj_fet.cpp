#include <iostream>
#include "libfbj_fet.hpp"

void FisherExactTest::minimize(int fix_r,int fix_c,int next)
{
	int r,c,i,j,l,move,last;//row,column,,,logic_bit,move,last cell
	last=FET_r*FET_c-FET_c-1;
	for (l=FET_c*fix_r+fix_c+next;l<last;++l)
	{	r=l/FET_c;
		c=l%FET_c;
		for (i=r+1;i<FET_r;++i)
		{	if (FET_a[r][c]==0) break;
			for (j=c+1;j<FET_c;++j)
			{	if (FET_a[i][j]==0) continue;
				move=FET_a[r][c]<FET_a[i][j]?FET_a[r][c]:FET_a[i][j];
				FET_a[r][c]-=move;
				FET_a[i][j]-=move;
				FET_a[r][j]+=move;
				FET_a[i][c]+=move;
				if (FET_a[r][c]==0) break;
			}
		}
	}
}

int FisherExactTest::increment()
{
	int l,r,c,i,j;
	for (l=FET_r*FET_c-FET_c-2;l>-1;--l)
	{	r=l/FET_c;
		c=l%FET_c;
		for (i=r+1;i<FET_r;++i)
		{	for (j=c+1;j<FET_c;++j)
			{	if (FET_a[r][j]==0 || FET_a[i][c]==0) continue;
				++FET_a[r][c];
				++FET_a[i][j];
				--FET_a[r][j];
				--FET_a[i][c];
				minimize(r,c,1);
				return 1;
			}
		}
	}
	return 0;
}

FACTORIAL FisherExactTest::pvalue(FACTORIAL cf)
{
	for (int i=0;i<FET_r;++i)
		for (int j=0;j<FET_c;++j)
			cf.remove(FET_a[i][j]);
	return cf;
}

void FisherExactTest::test()
{
	//calculate marginal total
	for (int i=0;i<FET_r;++i)
	{
		FET_rm[i]=0;
		for (int j=0;j<FET_c;++j)
			FET_rm[i]+=FET_a[i][j];
	}
	for (int j=0;j<FET_c;++j)
	{
		FET_cm[j]=0;
		for (int i=0;i<FET_r;++i)
			FET_cm[j]+=FET_a[i][j];
	}
	//calculate grand total
	FET_gt=0;
	for (int i=0;i<FET_r;++i)
		for (int j=0;j<FET_c;++j)
			FET_gt+=FET_a[i][j];
	
	//calculate common factor
	FACTORIAL cf;
	for (int i=0;i<FET_r;++i) cf.add(FET_rm[i]);
	for (int j=0;j<FET_c;++j) cf.add(FET_cm[j]);
	cf.remove(FET_gt);

	//set the threshold of p value
	FACTORIAL FET_th=pvalue(cf);
	FET_th.cal_log10();
	long double FET_th_log10prob = FET_th.result();
	long double FET_th_prob = pow(10,FET_th_log10prob);
	//std::cerr<<' '<<FET_th_log10prob<<' '<<FET_th_prob<<std::endl;
	
	//begin traverse all possible contingent tables
	FET_pv_both_tail=0;
	FET_pv_this_tail=0;
	FET_pv_othr_tail=0;
	int FET_b[FET_MAX_R][FET_MAX_C]; // backup array
	memcpy(FET_b, FET_a, sizeof(int)*FET_MAX_R*FET_MAX_C);
	bool same=false, past_same=false;
	minimize(0,0,0);
	do
	{
		if (!same)	same = !memcmp(FET_b, FET_a, sizeof(int)*FET_MAX_R*FET_MAX_C);
		else		past_same = true;
		//std::cerr<<FET_a[0][0]<<' '<<FET_a[0][1]<<' '<<FET_a[1][0]<<' '<<FET_a[1][1]<<' '<<same<<' '<<past_same;
		FACTORIAL p = pvalue(cf);
		p.cal_log10(); // added on 2016-04-27, solved the problem of big numbers
		long double log10prob = p.result();
		long double prob = pow(10,log10prob);
		//std::cerr<<' '<<log10prob<<' '<<prob;
		if (mid_p)
		{
			long double to_add = prob;
			if (same && !past_same) to_add/=2;
			if (past_same)		FET_pv_othr_tail += prob;
			else if (!same)		FET_pv_this_tail += prob;
			else			{	FET_pv_this_tail += to_add; FET_pv_othr_tail += to_add; }
			if 		(prob  < FET_th_prob) FET_pv_both_tail += prob;
			else if (prob == FET_th_prob) FET_pv_both_tail += prob/2;
		}
		else
		{
			if (past_same)		FET_pv_othr_tail += prob;
			else if (!same)		FET_pv_this_tail += prob;
			else			{	FET_pv_this_tail += prob; FET_pv_othr_tail += prob; }
			if 	(prob <= FET_th_prob)	FET_pv_both_tail += prob;
		}
		//std::cerr<<' '<<FET_pv_this_tail<<' '<<FET_pv_othr_tail<<' '<<FET_pv_both_tail<<std::endl;
	} while (increment());
	// Although the results may be correct when printed, the internal representation may be wrong due to incrementation of small errors.
	// For example, calculated from a table 1734 3 182 0, the cerr outputs are
	//   1734 3 182 0 1 0 -0.129897 0.741487 0.741487 0.741487 0.741487
	//   1735 2 181 1 1 1 -0.632004 0.233344 0.741487 0.974831 0.974831
	//   1736 1 180 2 1 1 -1.61387 0.0243291 0.741487 0.99916 0.99916
	//   1737 0 179 3 1 1 -3.07552 0.000840382 0.741487 1 1
	// Note that 0.741487+0.233344+0.0243291+0.000840382=1.000000482
	// This will create an error if the p-value (1.000000482) is used as an input for other functions.
	// Therefore, I added the following 3 lines to make sure this error does not happen.
	if (FET_pv_this_tail>1) FET_pv_this_tail=1;
	if (FET_pv_othr_tail>1) FET_pv_othr_tail=1;
	if (FET_pv_both_tail>1) FET_pv_both_tail=1;
	modified=false;
	memcpy(FET_a, FET_b, sizeof(int)*FET_MAX_R*FET_MAX_C);
}

