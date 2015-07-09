/*
 * TwoPLModel.cpp
 *
 *  Created on: 18 Jun 2014
 *      Author: cesandovalp
 */

#include <model/parameter/TwoPLModel.h>

TwoPLModel::TwoPLModel()
{
	parameterSet = NULL;
	probabilityMatrix = NULL;
}

void TwoPLModel::untransform()
{
    for (int i = 0; i < itemModel->getDataset()->countItems(); ++i)
        parameterSet[1][0][i] /= -parameterSet[0][0][i];
}

void TwoPLModel::getParameters(double * parameters)
{
    int i = 0;

    for (int j = 0; j < items; j++)
        parameters[i++] = parameterSet[0][0][j];
    for (int j = 0; j < items; j++)
        parameters[i++] = parameterSet[1][0][j];
}

void TwoPLModel::setParameters(double * parameters)
{
    int i = 0;

    for (int j = 0; j < items; j++)
        this->parameterSet[0][0][j] = parameters[i++];
    for (int j = 0; j < items; j++)
        this->parameterSet[1][0][j] = parameters[i++];
}

inline void TwoPLModel::successProbability(DimensionModel *dimensionModel, QuadratureNodes *quadNodes)
{
	int q = 0;
	double a_d, d_d, theta_d; // d stands from "double"

	if (dimensionModel != NULL)
		q = quadNodes->size();

	if (typeid(*dimensionModel) == typeid(UnidimensionalModel))
	{
		if (probabilityMatrix == NULL)
			//Creates the matrix if it is not already created
			probabilityMatrix = new Matrix<double>(q, items);
		
		for (int k = 0; k < q; k++)
		{
			for (int i = 0; i < items; i++)
			{
				theta_d = (*quadNodes->getTheta())(0, k);
				a_d = parameterSet[0][0][i];
				d_d = parameterSet[1][0][i];
				(*probabilityMatrix)(k, i) = successProbability(theta_d, a_d,
						d_d);
			}
		}
	}
}

inline double TwoPLModel::successProbability(double theta, double a, double d)
{
	long double exponential = (Constant::D_CONST * ((a * theta) + d));

	if (exponential > Constant::MAX_EXP)
		exponential = Constant::MAX_EXP;
	else if (exponential < -(Constant::MAX_EXP))
		exponential = -Constant::MAX_EXP;

	return (1 / (1 + exp(-exponential)));
}

adouble TwoPLModel::successProbabilityAD(double theta, adouble a, adouble d)
{
	adouble exponential = (Constant::D_CONST * ((a * theta) + d));
	if (exponential > Constant::MAX_EXP) {
		exponential = Constant::MAX_EXP;
	}

	else if (exponential < -(Constant::MAX_EXP)) {
		exponential = -Constant::MAX_EXP;
	}
	return (1 / (1 + CppAD::exp(-exponential)));
}

inline double TwoPLModel::successProbability(double theta, double * zita)
{
	long double exponential = (Constant::D_CONST * ((zita[0] * theta) + zita[1]));

	if (exponential > Constant::MAX_EXP)
		exponential = Constant::MAX_EXP;
	else if (exponential < -(Constant::MAX_EXP)) 
		exponential = -Constant::MAX_EXP;

	delete [] zita;

	return (1 / (1 + exp(-exponential)));
}

double *** TwoPLModel::getParameterSet() { return (this->parameterSet); }

void TwoPLModel::setParameterSet(double *** pair) { this->parameterSet = pair; }

double TwoPLModel::getProbability(int node, int item) { return ((*probabilityMatrix)(node, item)); }

void TwoPLModel::itemGradientAD(double* args, double* pars, int nargs, int npars, double* gradient)
{
	std::vector<adouble> args_t(2);

	for(int i = 0; i < 2; i++)
		args_t[i] = args[i];

	CppAD::Independent(args_t);
	size_t m = 1;
	std::vector<CppAD::AD<double>> Y(m);
	Y[0] = itemLogLikAD(args_t,pars,nargs,npars);

	CppAD::ADFun<double> f(args_t, Y);

	std::vector<double> jac(m * 2);
	std::vector<double> x(2);

	for(int i = 0; i < 2; i++)
		x[i] = args[i];

	jac  = f.Jacobian(x);

	memset(gradient, 0, sizeof(double) * 2);

	for(int i = 0; i < 2; i++)
		gradient[i] = jac[i];
}

double TwoPLModel::itemLogLik (double* args, double* pars, int nargs, int npars)
{
	int nP = 0;
	unsigned int q, items;
	int index = 0;
	double *theta, *r, *f;
	double a, b;
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
	
	if(abs(a)>5)
		a = 0.851;

	double dd = 0;
	dd = -b/a;
	
	if(abs(dd)>5)
		b = 0;

	for (unsigned int k = 0; k < q; ++k)
	{
		tp = (TwoPLModel::successProbability ( theta[k], a,b));

		if (tp < 1e-08)
			tp = 1e-08;

		tq = 1-tp;

		if (tq < 1e-08)
			tq = 1e-08;

		sum += (r[k]*log(tp))+(f[k]-r[k])*log(tq);
	}

	// TODO: Analize
	args[0] = a;
	args[1] = b;

	delete[] theta;
	delete[] f;
	delete[] r;

	return (-sum);
}

adouble TwoPLModel::itemLogLikAD(std::vector<adouble> args, double* pars, int nargs, int npars)
{
	double *theta, *r, *f;
	unsigned int nP, q, items, index;
	adouble a, b, tp, tq, sum;

	sum = nP = index = 0;

	a = args[0];
	b = args[1];
	
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

	for (unsigned int k = 0; k < q; ++k)
	{
		tp = (TwoPLModel::successProbabilityAD(pars[k + 2], a, b));

		if (tp == 0)
			tp = 1e-08;
		
		tq = 1 - tp;

		if (tq == 0)
			tq = 1e-08;

		sum += (r[k]*log(tp))+(f[k]-r[k])*log(tq);
	}

	// TODO: Analize
	args[0] = a;
	args[1] = b;

	delete[] theta;
	delete[] f;
	delete[] r;

	return (-sum);
}

void TwoPLModel::printParameterSet(ostream& out)
{
	cout << "\"a\" \"b\" \"c\"" << endl;
	
	for (int i = 0; i < items; i++)
		cout << parameterSet[0][0][i] << " "
		     << parameterSet[1][0][i] << " "
		     << 0 << endl;
}

