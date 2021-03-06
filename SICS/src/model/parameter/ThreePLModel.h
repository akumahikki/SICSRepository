/*
 * ThreePLModel.h
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#ifndef THREEPLMODEL_H_
#define THREEPLMODEL_H_

#include <typeinfo>
#include <model/parameter/ParameterModel.h>
#include <model/item/DichotomousModel.h>
#include <model/item/PolytomousModel.h>
#include <model/dimension/UnidimensionalModel.h>
#include <model/dimension/MultidimensionalModel.h>
#include <model/dimension/MultiUniDimModel.h>
#include <type/PatternMatrix.h>
#include <type/Constant.h>
#include <cmath>
#include <util/util.h>
#include <type/QuadratureNodes.h>
/**
 * Model for the 3pl model, uses parameters a d c.
 * unidimensional
 * */
class ThreePLModel: public ParameterModel
{

private:

	QuadratureNodes* nodes;

public:

	static double itemLogLik2(double*, double* , int, int);
	static adouble itemLogLikAD(std::vector<adouble>, double*, int, int);
	static void itemGradient2(double*, double*, int, int, double*);
	static void itemGradientAD (double*, double*, int, int, double*);
	static double successProbability(double, double, double, double);
	static adouble successProbabilityAD(double, adouble, adouble, adouble);
	double successProbability(double, double *);
	static double successProbability_cPrime (double, double, double, double);
	virtual void transform();
	virtual void untransform();

	// Constructor
	ThreePLModel();

	// Methods
	void setEstimationNodes(QuadratureNodes * );
	void successProbability(DimensionModel *, QuadratureNodes *);
	// Getters and Setters
	double*** getParameterSet() ;
	void getParameters(double *);
	void setParameters(double *);
	void setParameterSet(double ***);
	void printParameterSet(ostream&);
	double getProbability(int, int);
};

#endif /* THREEPLMODEL_H_ */
