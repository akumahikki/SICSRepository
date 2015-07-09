/*
 * OnePLModel.cpp
 *
 *  Created on: Nov 16, 2014
 *      Author: anmrodriguezre
 */

#include <model/parameter/OnePLModel.h>

OnePLModel::OnePLModel()
{
	parameterSet = NULL;
	probabilityMatrix = NULL;
}

inline void OnePLModel::successProbability(DimensionModel *dimensionModel, QuadratureNodes * quadNodes)
{
	int q = 0;

	if (dimensionModel != NULL)
		q = quadNodes->size();

	if (typeid(*dimensionModel) == typeid(UnidimensionalModel))
	{
		if (probabilityMatrix == NULL)
			//Creates the matrix if it is not already created
			probabilityMatrix = new Matrix<double>(q, items);

		for (int k = 0; k < q; k++)
			for (int i = 0; i < items; i++)
				// Rasch Success Probability Function
				(*probabilityMatrix)(k, i) =
					successProbability((*quadNodes->getTheta())(0, k), (parameterSet[0][0][i]));
	}
}

inline double OnePLModel::successProbability(double theta, double b)
{
	long double exponential = ((theta) - b);

	if (exponential > Constant::MAX_EXP)
		exponential = Constant::MAX_EXP;
	else if (exponential < -(Constant::MAX_EXP))
		exponential = -Constant::MAX_EXP;

	return (1 / (1.0 + exp(-exponential)));
}

adouble OnePLModel::successProbabilityAD(double theta, adouble b)
{
	adouble exponential = ((theta) - b);

	if (exponential > Constant::MAX_EXP)
		exponential = Constant::MAX_EXP;
	else if (exponential < -(Constant::MAX_EXP))
		exponential = -Constant::MAX_EXP;

	return (1 / (1.0 + CppAD::exp(-exponential)));
}

inline double OnePLModel::successProbability(double theta, double * zita)
{
	long double exponential = ((theta) - zita[0]);

	if (exponential > Constant::MAX_EXP)
		exponential = Constant::MAX_EXP;
	else if (exponential < -(Constant::MAX_EXP))
		exponential = -Constant::MAX_EXP;

	delete [] zita;

	return (1 / (1.0 + exp(-exponential)));
}

void OnePLModel::setParameterSet(double*** par) { this->parameterSet = par; }

double*** OnePLModel::getParameterSet() { return (this->parameterSet); }

void OnePLModel::getParameters(double * parameters)
{
	for (int i = 0; i < items; i++)
		parameters[i] = parameterSet[0][0][i];
}

void OnePLModel::setParameters(double * parameters)
{
	for (int i = 0; i < items; i++)
		this->parameterSet[0][0][i] = parameters[i];
}

double OnePLModel::getProbability(int node, int item) { return ((*probabilityMatrix)(node, item)); }

void OnePLModel::printParameterSet(ostream& out)
{
	out << "\"a\" \"b\" \"c\"" << "\n";

	for (int k = 0; k < items; k++)
		out << 1 << " " << (parameterSet[0][0][k]) << " " << 0 << "\n";
}

void OnePLModel::itemGradientAD(double* args, double* pars, int nargs, int npars, double* gradient)
{
	std::vector<adouble> args_t(1);

	args_t[0] = args[0];

	CppAD::Independent(args_t);
	size_t m = 1;
	std::vector<CppAD::AD<double>> Y(m);
	Y[0] = itemLogLikAD(args_t,pars,nargs,npars);

	CppAD::ADFun<double> f(args_t, Y);

	std::vector<double> jac(m);
	std::vector<double> x(1);

	x[0] = args[0];

	jac  = f.Jacobian(x);

	memset(gradient, 0, sizeof(double));

	gradient[0] = jac[0];
}

double OnePLModel::itemLogLik (double* args, double* pars, int nargs, int npars)
{
	int nP = 0;
	unsigned int q, items;
	int index = 0;
	double *theta, *r, *f;
	double a;
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

	for (unsigned int k = 0; k < q; ++k)
	{
		tp = (OnePLModel::successProbability ( theta[k], a));
		
		if (tp < 1e-08)
			tp = 1e-08;

		tq = 1-tp;

		if (tq < 1e-08)
			tq = 1e-08;

		sum += (r[k]*log(tp))+(f[k]-r[k])*log(tq);
	}

	// TODO: Analize
	args[0] = a;

	delete[] theta;
	delete[] f;
	delete[] r;

	return (-sum);
}

adouble OnePLModel::itemLogLikAD(std::vector<adouble> args, double* pars, int nargs, int npars)
{
	double *theta, *r, *f;
	unsigned int nP, q, items, index;
	adouble a, b, tp, tq, sum;

	sum = nP = index = 0;

	a = args[0];
	
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
		tp = (OnePLModel::successProbabilityAD(pars[k + 2], a));

		if (tp == 0)
			tp = 1e-08;
		
		tq = 1 - tp;

		if (tq == 0)
			tq = 1e-08;

		sum += (r[k]*log(tp))+(f[k]-r[k])*log(tq);
	}

	// TODO: Analize
	args[0] = a;

	delete[] theta;
	delete[] f;
	delete[] r;

	return (-sum);
}