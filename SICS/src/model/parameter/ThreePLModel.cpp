/*
 * ThreePLModel.cpp
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#include <model/parameter/ThreePLModel.h>
#include <string>

ThreePLModel::ThreePLModel()
{
	parameterSet = NULL;
	probabilityMatrix = NULL;
	nodes = NULL;
} 

void ThreePLModel::setEstimationNodes(QuadratureNodes* n) { this->nodes = n; }

void ThreePLModel::successProbability(DimensionModel *dimensionModel, QuadratureNodes * quadNodes)
{
	int q = 0;
	double a_d, d_d, c_d, theta_d; // d stands from "double"

	if ( dimensionModel != NULL )
		q = quadNodes->size();

	if(typeid(*dimensionModel)==typeid(UnidimensionalModel))
	{
		if(probabilityMatrix == NULL)
			//Creates the matrix if it is not already created
			probabilityMatrix = new Matrix<double>(q,items);

		for (int k = 0; k < q; k++)
		{
			for ( int i = 0; i < items; i++ )
			{
				// 3PL Success Probability Function
				theta_d = (*quadNodes->getTheta())(0,k);
				a_d = parameterSet[0][0][i];
				d_d = parameterSet[1][0][i];
				c_d = parameterSet[2][0][i];

				(*probabilityMatrix)(k,i) = successProbability ( theta_d, a_d, d_d, c_d );
			}
		}
	}
}

double *** ThreePLModel::getParameterSet() { return (this->parameterSet); }

void ThreePLModel::setParameterSet(double ***) { this->parameterSet = parameterSet; }

double ThreePLModel::successProbability(double theta, double a, double d, double c)
{
	long double exponential = (Constant::NORM_CONST)*(a*theta+d);

	if ( exponential > Constant::MAX_EXP )
		exponential = Constant::MAX_EXP;

	else if ( exponential < -(Constant::MAX_EXP*1.0) )
		exponential = -Constant::MAX_EXP;

	exponential = exp(-exponential) ;
	double ec = exp(c);

	return ( (ec/(1+ec)) + (1 - (ec/(1+ec))) * (1/(1+exponential)) );
}

adouble ThreePLModel::successProbabilityAD(double theta, adouble a, adouble d, adouble c)
{
	adouble exponential = (Constant::NORM_CONST)*(a*theta+d);

	if ( exponential > Constant::MAX_EXP )
		exponential = Constant::MAX_EXP;

	else if ( exponential < -(Constant::MAX_EXP*1.0) )
		exponential = -Constant::MAX_EXP;

	exponential = CppAD::exp(-exponential) ;
	adouble ec = CppAD::exp(c);

	return ( (ec/(1+ec)) + (1 - (ec/(1+ec))) * (1/(1+exponential)) );
}

void ThreePLModel::getParameters(double * parameters)
{
	int i = 0;

	for (int j = 0; j < items; j++)
		parameters[i++] = parameterSet[0][0][j];
	for (int j = 0; j < items; j++)
		parameters[i++] = parameterSet[1][0][j];
	for (int j = 0; j < items; j++)
		parameters[i++] = parameterSet[2][0][j];
}

void ThreePLModel::setParameters(double * parameters)
{
	int i = 0;

	for (int j = 0; j < items; j++)
		this->parameterSet[0][0][j] = parameters[i++];
	for (int j = 0; j < items; j++)
		this->parameterSet[1][0][j] = parameters[i++];
	for (int j = 0; j < items; j++)
		this->parameterSet[2][0][j] = parameters[i++];
}

double ThreePLModel::successProbability(double theta, double * zita)
{
	long double exponential = (Constant::NORM_CONST)*(zita[0]*theta+zita[1]);

	if ( exponential > Constant::MAX_EXP )
		exponential = Constant::MAX_EXP;

	else if ( exponential < -(Constant::MAX_EXP*1.0) )
		exponential = -Constant::MAX_EXP;

	exponential = exp(-exponential) ;
	double ec = exp(zita[2]);

	return ( (ec/(1+ec)) + (1 - (ec/(1+ec))) * (1/(1+exponential)) );
}

double ThreePLModel::getProbability(int node, int item) { return ((*probabilityMatrix)(node, item)); }


double ThreePLModel::itemLogLik2 (double* args, double* pars, int nargs, int npars)
{
	int nP = 0;
	unsigned int q, items;
	int index = 0;
	double *theta, *r, *f;
	double a, b, c;
	double sum=0;
	long double tp , tq;
	
	index = pars[npars-1];
	
	// Obtain q
	q = pars[nP ++]; // q is obtained and npars is augmented
	
	// Obtain I
	items = pars[nP ++];
	theta = new double[q];
	r = new double[q];
	f = new double[q];
	
	// Obtain theta
	for (unsigned int k=0; k<q; k++)
		theta[k] = pars[nP ++];
	
	// Obtain f
	for (unsigned int k=0; k<q; k++)
		f[k] = pars[nP ++];

	// Obtain r that becomes a vector
	for (unsigned int k=0; k<q; k++)
	{
		nP += index;
		r[k] = pars[nP];
		nP += (items-index); 
	}
	
	// Obtain a
	a = args[0];
	// Obtain b
	b = args[1];
	// Obtain c
	c = args[2];

	if(abs(a)>5)
		a = 0.851;

	double dd = 0;
	dd = -b/a;

	if(abs(dd)>5)
		b = 0;

	if(abs(c)>5)
		c = 0.3;

	for (unsigned int k = 0; k < q; ++k)
	{
		tp = (ThreePLModel::successProbability ( theta[k], a,b,c));
		
		if (tp < 1e-08)
			tp = 1e-08;
		
		tq = 1-tp;
		
		if (tq < 1e-08)
			tq = 1e-08;
		
		sum += (r[k]*log(tp))+(f[k]-r[k])*log(tq);
	}

	args[0] = a;
	args[1] = b;
	args[2] = c;

	delete[] theta;
	delete[] f;
	delete[] r;

	return (-sum);
}

adouble ThreePLModel::itemLogLikAD(std::vector<adouble> args, double* pars, int nargs, int npars)
{
	double *theta, *r, *f;
	unsigned int nP, q, items, index;
	adouble a, b, c, tp, tq, sum;

	sum = nP = index = 0;

	a = args[0];
	b = args[1];
	c = args[2];

	q = pars[nP++];
	items = pars[nP++];
	index = pars[npars - 1];

	theta = new double[q];
	r = new double[q];
	f = new double[q];

	// Obtain theta
	for (unsigned int k=0; k<q; k++)
		theta[k] = pars[nP++];

	// Obtain f
	for (unsigned int k=0; k<q; k++)
		f[k] = pars[nP++];

	// Obtain r that becomes a vector
	for (unsigned int k=0; k<q; k++)
	{
		nP += index;
		r[k] = pars[nP];
		nP += (items-index);
	}

	adouble bound1 = 5;
	if(CppAD::abs(a) > bound1)
	{
		a = 0.851;
		cout << "→" << endl;
	}
	adouble dd = -b/a;
	if(CppAD::abs(dd) > bound1)
	{
		b = 0;
		cout << "↓" << endl;
	}
	if(CppAD::abs(c) > bound1)
	{
		c = 0.3;
		cout << "ø" << endl;
	}

	for (unsigned int k = 0; k < q; ++k)
	{
		tp = (ThreePLModel::successProbabilityAD(theta[k], a, b, c));

		if (tp < 1e-08)
			tp = 1e-08;
		
		tq = 1 - tp;

		if (tq < 1e-08)
			tq = 1e-08;

		sum += (r[k]*log(tp))+(f[k]-r[k])*log(tq);
	}

	// TODO: Analize
	args[0] = a;
	args[1] = b;
	args[2] = c;

	delete[] theta;
	delete[] f;
	delete[] r;

	return (-sum);
}

void ThreePLModel::itemGradientAD(double* args, double* pars, int nargs, int npars, double* gradient)
{
	std::vector<adouble> args_t(3);

	for(int i = 0; i < 3; i++)
		args_t[i] = args[i];

	CppAD::Independent(args_t);
	size_t m = 1;
	std::vector<CppAD::AD<double>> Y(m);
	Y[0] = itemLogLikAD(args_t,pars,nargs,npars);

	CppAD::ADFun<double> f(args_t, Y);

	std::vector<double> jac(m * 3);
	std::vector<double> x(3);

	for(int i = 0; i < 3; i++)
		x[i] = args[i];

	jac  = f.Jacobian(x);

	memset(gradient, 0, sizeof(double) * 3);

	for(int i = 0; i < 3; i++)
		gradient[i] = jac[i];
}

double ThreePLModel::successProbability_cPrime(double theta, double a, double b, double c)
{
	long double cPrime = exp(c)/(1+exp(c));
	return (successProbability ( theta, a, b, cPrime ));
}

void ThreePLModel::printParameterSet(ostream& out)
{
	out << "\"a\" \"b\" \"c\"" << endl;

	for (int i = 0; i < items; i++)
		out << parameterSet[0][0][i] << " "
		    << parameterSet[1][0][i] << " "
		    << parameterSet[2][0][i] << endl;
}

