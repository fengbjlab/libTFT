#include "libfbj_base.hpp"
#include "libfbj_math.hpp"
#include "libfbj_regress.hpp"

int multifit_logistic_validate (const Eigen::MatrixXd & X,
								const Eigen::VectorXd & y,
								const int p_start) // if X[][0] is constant term, skip it by p_start=1, otherwise 0
{
	const int ni = X.rows();		// number of individuals (observations)
	const int np = X.cols();		// number of parameters
	
	// check validity of vector y
	int responders=0;
	for (int i=0;i<ni;i++)
	{
		double response = y(i);
		if (response)
		{
			if (response!=1) exit_error("Logistic regression: response should be either 0 or 1");
			responders++;
		}
	}
	if (responders==0 || responders==ni) return 0; // exit_error("Logistic regression: outcome does not vary");
	
	// check validity of matrix X
	for (int j=p_start;j<np;j++)
	{
		double min0= INFINITY, min1= INFINITY;
		double max0=-INFINITY, max1=-INFINITY;
		for (int i=0;i<ni;i++)
		{
			double value = X(i,j);
			if (y(i)!=0) {	if (value>max1) max1=value; if (value<min1) min1=value; }
			else		 {	if (value>max0) max0=value; if (value<min0) min0=value; }
		}
		if (max0==min0 && max1==min1 && max0==max1) return 0; // exit_error("Logistic regression: one independent variable does not vary");
		if (max1 < min0)	return 0; // exit_error("Logistic regression: 1 variable predicts data perfectly");
		if (min1 > max0)	return 0; // exit_error("Logistic regression: 1 variable predicts data perfectly");
		if (max0==min0 || max1==min1) return 0; // exit_error("Logistic regression: 1 variable predicts an outcome perfectly.");
		// The last 3 situations would be excluded by STATA. Although I can still calculate, but the p-value seems not right.
		// At these situations, I recommend linear regression because they can be handled correctly.
	}
	return 1;
}

// Fit: logit(p) = X c, where X is an MxN matrix of M observations for N variables, p is probability of response.
void multifit_logistic (const Eigen::MatrixXd & X,			// predictor matrix
						const Eigen::VectorXd & y,			// response vector
						Eigen::VectorXd & c,				// coefficients
						Eigen::MatrixXd & cov,				// covariance
						double & chisq,						// likelihood ratio test
						multifit_logistic_workspace& work)	// workspace
{
	if		(X.rows() != y.size())						exit_error("number of observations in y does not match rows of matrix X");
	else if (X.cols() != c.size())						exit_error("number of parameters c does not match columns of matrix X");
	else if (cov.rows() != cov.cols())					exit_error("covariance matrix is not square");
	else if (c.size()!=cov.rows())						exit_error("number of parameters does not match size of covariance matrix");
	else if (X.rows()!=work.ni || X.cols()!=work.np)	exit_error("size of workspace does not match size of observation matrix");
	else
    {
		const int ni = X.rows();		// number of individuals (observations)
		const int np = X.cols();		// number of parameters
		chisq=0;						// no effect if there's no data validation bellow
		Eigen::VectorXd & p = work.p;	// prob of response
		Eigen::VectorXd & V = work.V;	// logit(p)
		Eigen::VectorXd & t = work.t;	// temp vector
		Eigen::MatrixXd & X_= work.X_;	// X~
		double n1 = y.sum();			// number of response, assume y=0/1, no checking!!!
		double llk0 = n1 * log(n1/ni) + (ni-n1) * log(1-n1/ni);	// log likelihood of NULL model
		double last_llk=-INFINITY;								// last log likelihood value
		c.setZero(np);											// set initial value of all coef to 0
		
		for (int iter=0; iter < 20; iter++)						// for a limited number of iterations
		{
			// calculate p and V
			t = X * c;
			for (int i=0; i<ni; i++)
			{
				p(i) = 1 / (1+exp(-t(i)));
				V(i) = p(i) * (1-p(i));
			}
			
			// b +=  inv( t(X) * X_ ) * t(X) * ( Y - p )
			for (int i=0;i<ni;i++)
				for (int j=0;j<np;j++)
					X_(i,j) = V(i) * X(i,j);
			c = c + (X.transpose() * X_).inverse() * X.transpose() * ( y - p );

			// check for convergence by log-likelihood
			double llk=0;
			for (int i=0;i<ni;i++)
				if (y(i))	llk+=log(p(i));
				else		llk+=log(1-p(i));
			if (llk-last_llk < 0.0001) break;
			last_llk=llk;
		}
		
		cov = (X_.transpose() * X).inverse();
		chisq = 2 * (last_llk - llk0);
    }
}

int multifit_linear_validate(const Eigen::MatrixXd & X,
							 const Eigen::VectorXd & y,
							 const int p_start) // if X[][0] is constant term, skip it by p_start=1, otherwise 0
{
	const int ni = X.rows();		// number of individuals (observations)
	const int np = X.cols();		// number of parameters
	
	// check validity of matrix X
	for (int j=p_start;j<np;j++)
	{
		double minV= INFINITY;
		double maxV=-INFINITY;
		for (int i=0;i<ni;i++)
		{
			double value = X(i,j);
			if (value>maxV) maxV=value;
			if (value<minV) minV=value;
		}
		if (maxV==minV) return 0; // exit_error("Linear regression: one independent variable does not vary");
	}
	return 1;
}

void multifit_linear(const Eigen::MatrixXd		& X,
					 const Eigen::VectorXd		& y,
					 Eigen::VectorXd			& c,
					 Eigen::MatrixXd			& cov,
					 double						& chisq,
					 multifit_linear_workspace	& work )
{
	if		(X.rows() != y.size())	exit_error ("number of observations in y does not match rows of matrix X");
	else if (X.cols() != c.size())	exit_error ("number of parameters c does not match columns of matrix X");
	else if (cov.rows()!=cov.cols())exit_error ("covariance matrix is not square");
	else if (c.size()!=cov.rows())	exit_error ("number of parameters does not match size of covariance matrix");
	else if (X.rows()!=work.ni)		exit_error ("size of workspace does not match size of observation matrix");
	else if (X.cols()!=work.np)		exit_error ("size of workspace does not match size of observation matrix");
	else if (X.cols()>=X.rows())	exit_error ("number of observations must > number of parameters for linear regression.");
	else
	{
		work.t = (X.transpose() * X).inverse() * X.transpose();
		work.H = X * work.t;
		work.e = (work.I - work.H) * y;
		c = work.t * y;
		chisq = work.e.squaredNorm();
		double sigma_square = work.e.transpose() * work.e;
		sigma_square /= (work.ni-work.np);
		cov = sigma_square * (X.transpose() * X).inverse();
	}
}

// Do regression with a constance term. Return p-value for the first beta.
// Here one-sided means P(T>=t), so if Beta_1 is <0, then the p-value is large.

double pv_1st_linear(const Eigen::MatrixXd	& X, const Eigen::VectorXd	& y, bool one_sided)
{
	const int ni = X.rows();		// number of individuals (observations)
	const int np = X.cols();		// number of parameters
	Eigen::VectorXd c(np);			// beta
	Eigen::MatrixXd cov(np,np);		// covariance matrix
	if (multifit_linear_validate(X,y,1))
	{
		multifit_linear_workspace work(ni,np);
		double chisq;
		multifit_linear(X,y,c,cov,chisq,work);
		double standard_err = sqrt(cov(1,1));
		double dl = ni-(np-1)-1;						// degrees of freedom for t statistics
		double t = c(1)/standard_err;
		if (one_sided)	return cdf_t_q(t,dl);
		else			return cdf_t_2sided_pv(t,dl);
	}
	return std::numeric_limits<double>::signaling_NaN();
}

double pv_1st_logistic(const Eigen::MatrixXd	& X, const Eigen::VectorXd	& y, bool one_sided)
{
	const int ni = X.rows();		// number of individuals (observations)
	const int np = X.cols();		// number of parameters
	Eigen::VectorXd c(np);			// beta
	Eigen::MatrixXd cov(np,np);		// covariance matrix
	if (multifit_logistic_validate(X,y,1))
	{
		multifit_logistic_workspace work(ni,np);
		double chisq;
		multifit_logistic(X,y,c,cov,chisq,work);
		double standard_err=sqrt(cov(1,1));
		double Z = c(1)/standard_err;
		if (one_sided)	return cdf_norms_q(Z);
		else			return cdf_norms_2sided_pv(Z);
	}
	return std::numeric_limits<double>::signaling_NaN();
}

double or_1st_logistic(const Eigen::MatrixXd	& X, const Eigen::VectorXd	& y, bool one_sided)
{
	const int ni = X.rows();		// number of individuals (observations)
	const int np = X.cols();		// number of parameters
	Eigen::VectorXd c(np);			// beta
	Eigen::MatrixXd cov(np,np);		// covariance matrix
	if (multifit_logistic_validate(X,y,1))
	{
		double chisq;
		multifit_logistic_workspace work(ni,np);
		multifit_logistic(X,y,c,cov,chisq,work);
		return exp(c(1));
	}
	return std::numeric_limits<double>::signaling_NaN();
}

