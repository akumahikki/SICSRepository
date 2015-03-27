/*
 * EM.cpp
 *
 *  Created on: May 29, 2014
 *      Author: mirt
 */

#include <estimation/classical/EMEstimation.h>
#include <util/util.h>

EMEstimation::EMEstimation() {
	iterations = 0;
	model = NULL;
	profiler = NULL;
	f = NULL;
	r = NULL;
	logger = NULL;
	optim = NULL;
	convergenceSignal = false;
	optim = new Optimizer();
	quadNodes = NULL;
	estimator = NULL;
}

EMEstimation::~EMEstimation() {
	if (f != NULL) {
		delete f;
	}
	if (r != NULL) {
		delete r;
	}
	if (logger != NULL) {
		delete logger;
	}
	if (optim != NULL) {
		delete optim;
	}
}

void EMEstimation::setProfiler(Trace* t) {
	profiler = t;
	//model->parameterModel->setProfiler(t);
}

/**
 * Sets the model to be estimated, currently only supports 3PL model
 */
void EMEstimation::setModel(Model* Model) {
	int q;
	int It;
	this->model = Model;
	q = quadNodes->size();
	It = model->getItemModel()->countItems();

	f = new Matrix<double>(1, q);
	r = new Matrix<double>(q, It);

	//Discriminate by models
	if (Model->Modeltype() == Constant::THREE_PL) {
		estimator = new EM3PL(model, quadNodes, f, r); //Initializes estimator
		return;
	}

	if (Model->Modeltype() == Constant::RASCH_A1) {
		estimator = new EM1PL(model, quadNodes, f, r); //Initializes estimator
		return;
	}

	if (Model->Modeltype() == Constant::TWO_PL) {
		estimator = new EM2PL(model, quadNodes, f, r); //Initializes estimator with Cristian's 2PL Model
		return;
	}

	if (Model->Modeltype() == Constant::RASCH_A_CONSTANT) {
		estimator = new EM1PLAC(model, quadNodes, f, r);
		return;
	}

}
/**
 * Sets the initial values for the estimation, use this for inputting a matrix as initial values
 */
void EMEstimation::setInitialValues(double*** parameterSet) {
	estimator->setInitialValues(parameterSet, model);
}
/**
 * Sets the initial values according to a method of calculating the values
 * Possible methods :
 * ANDRADE,
 * OSPINA,
 * RANDOM,
 *
 * The default method is OSPINA , this is the fastest method according to the SICS calculations
 */
void EMEstimation::setInitialValues(int method) {
	estimator->setInitialValues(method, model);
}
/**
 * Main loop of the EM estimation
 * orchestrates the parameters for the estimation, and holds the estimation
 * for the iterations.
 *
 * TODO : read maxiterations as a input parameter , idea : calculate the max iterations depending on the items
 * TODO : Output last estimation onto the json for recovery in the program.
 */
void EMEstimation::estimate() {
	estimator->transform();
	iterations = 0;

	double ** args_hist;
	int nargs;
	int size = 3 * model->getItemModel()->getDataset()->countItems();
	(args_hist) = new double*[3];
	(args_hist)[0] = new double[size];
	(args_hist)[1] = new double[size];
	(args_hist)[2] = new double[size];

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < size; j++)
			args_hist[i][j] = 0;

	for (;;)
	{
		//cout << iterations << endl;
		estimator->stepE();
		estimator->stepM(&args_hist, &nargs);
		estimator->stepRamsay(&args_hist, &nargs, size, iterations > 5 && (iterations) % 3 == 0);

		convergenceSignal = model->itemParametersEstimated;

		if (iterations++ > Constant::MAX_EM_ITERS || convergenceSignal)
			break;
	}
	estimator->untransform();
	model->printParameterSet(cout);
	//	cout << "Total time from estimation " << profiler->dr("estim") << endl
	//			<< "E step time : " << profiler->dr("Et") << endl
	//			<< "M step time : " << profiler->dr("Mt") << endl;
}

/**Returns the iterations that took the estimation to obtain an answer*/
int EMEstimation::getIterations() const {
	return (iterations);
}

QuadratureNodes* EMEstimation::getQuadratureNodes() const {
	return (quadNodes);
}

void EMEstimation::setQuadratureNodes(QuadratureNodes* nodes) {
	this->quadNodes = nodes;
}
