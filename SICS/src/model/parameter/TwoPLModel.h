/*
 * TwoPLModel.h
 *
 *  Created on: 18 Jun 2014
 *      Author: cesandovalp
 */

#ifndef TWOPLMODEL_H_
#define TWOPLMODEL_H_

#include <model/parameter/ParameterModel.h>
#include <type/Constant.h>
#include <cmath>
#include <typeinfo>
#include <model/dimension/UnidimensionalModel.h>
#include <model/dimension/MultidimensionalModel.h>
#include <model/dimension/MultiUniDimModel.h>
#include <model/item/PolytomousModel.h>
#include <model/item/DichotomousModel.h>
#include <util/util.h>

class TwoPLModel: public ParameterModel
{

public:

	TwoPLModel();

	// Methods
	inline void successProbability(DimensionModel *, QuadratureNodes *);
	inline static double successProbability(double , double , double );
	static adouble successProbabilityAD(double, adouble, adouble);
	double successProbability(double, double*);
	static double patternProbability();
	static void itemHessian(double*, double*, int, int, double*);
	static void itemGradient (double*, double*, int, int, double*);
	static void itemGradientAD (double*, double*, int, int, double*);
	static double itemLogLik (double*, double*, int, int);
	static adouble itemLogLikAD(std::vector<adouble>, double*, int, int);

	double *** getParameterSet();
	double getProbability(int, int);
	void getParameters(double *);
	void setParameterSet(double***);
	void setParameters(double *);
	void printParameterSet(ostream&);
};

#endif /* TWOPLMODEL_H_ */
