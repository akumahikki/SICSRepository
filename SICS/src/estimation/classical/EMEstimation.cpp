/*
 * EM.cpp
 *
 *  Created on: May 29, 2014
 *      Author: mirt
 */

#include "EMEstimation.h"

EMEstimation::EMEstimation() {

	model = NULL;
	f = NULL;
	r = NULL;
	logger = NULL;
	optim = NULL;
	optim = new Optimizer();

}

EMEstimation::~EMEstimation() {
	if( f!= NULL ) {
		delete f;
	}
	if( r!= NULL ) {
		delete r;
	}
	if (logger != NULL ) {
		delete logger;
	}
	if (optim != NULL ) {
			delete optim;
	}
}
/*
 * Model is set.
 */
void EMEstimation::setModel(Model* Model){
	int q;
	int I;
	this->model=Model;
	q = model->getDimensionModel()->getLatentTraitSet()->getTheta()->nC();
	I = model->getItemModel()->countItems();

	f = new Matrix<double> (1,q);
	r = new Matrix<double> (q,I);

}
/*
 * Sets or selects the initial values
 */
//Sets
void EMEstimation::setInitialValues(map<Parameter, Matrix<double>*> parameterSet){
	model->getParameterModel()->setParameterSet(parameterSet);
}
//Selects
void EMEstimation::setInitialValues(string method){
	/*TODO
	 * Possible methods
	 * ANDRADE
	 * OSPINA
	 * RANDOM
	 *
	 * The default method is OSPINA
	 */
}

void EMEstimation::stepE () {
	/*
	 * What we need
	 * q
	 * a pattern iterator
	 * item number
	 * success probability matrix
	 * thetas
	 * weights
	 * parameter set
	 */
	//Dataset by patterns

	PatternMatrix* data = dynamic_cast <PatternMatrix * > (model->getItemModel()->getDataset()) ;
	//Pattern iterator is data->iterator
	//Item number
	const double items = data->countItems();
	//Success probability matrix is obtained via pm->getProbability(int,int)
	ParameterModel* pm = model->getParameterModel();
	//Thetas
	Matrix<double>* thetas = model->getDimensionModel()->getLatentTraitSet()->getTheta();
	//Amount of nodes
	const int q = thetas->nC();
	//Weights
	Matrix<double>* weights = model->getDimensionModel()->getLatentTraitSet()->getWeight();
	//A Matrix
	Matrix<double>* A = model->getParameterModel()->getParameterSet()[a];
	//B Matrix
	Matrix<double>* B = model->getParameterModel()->getParameterSet()[b];
	//C Matrix
	Matrix<double>* C = model->getParameterModel()->getParameterSet()[c];
	//Auxiliar array for the nodes
	long double faux[q];
	long double sum = 0.0;
	int pat = 0; // N of patterns
	//Restart f and r to zero
	f->reset();
	r->reset();

	for (data->resetIterator();data->checkEnd();data->iterate()){
		//Initialize faux in 1 to later calculate the productory
		for ( int k = 0; k < q; k++ ) {
					faux[k] = 1;
				}
		//Calculate g*(k) for all the k's
		//first calculate the P for each k and store it in the array f aux
		for (int k = 0; k < q; k++) {
			//Calculate the p (iterate over the items in the productory)
			for ( unsigned int i = 0; i < items; i++ ) {
				double prob = pm->getProbability(k,i);
				if(!data->getCurrentBitSet()[i]){
					prob = 1-prob;
				}
				faux[k] = faux[k]*prob;
			}

		//At this point the productory is calculated and faux[k] is equivalent to p(u_j,theta_k)
		//Now multiply by the weight

		faux[k] = faux[k] * (*weights)(0,k);
		}

		//compute the total of the p*a' (denominator of the g*)
		sum = 0.0;
		for (int k = 0; k < q; k++) {
			sum += faux[k];
		}

		for (int k = 0; k < q; k++) {
			faux[k] = faux[k] / sum;	//This is g*_j_k

			//Multiply the f to the frequency of the pattern
			faux[k] = ((long double) data->getCurrentFrequency()) * faux[k];
			(*f)(1,k) = faux[k] = (*f)(0,k);//TODO Possible mistake here JLP

			//Now selectively add the faux to the r
			for ( unsigned int i = 0; i < items; i++ ) {
				if (data->getCurrentBitSet()[i]) {
					(*r)(k,i)= (*r)(k,i) + faux[k];
				} // if
			} // for
		} // for


	}
} //end E step

void EMEstimation::stepM(){
	/*
	 * TODO Step M for EML3M
	 */
	//Step M implementation using the BFGS Algorithm
	/*
	 * What we need
	 * fptr the pointer to loglik
	 * gprt the pointer to gradient
	 * hessptr the pointer to hessian matrix calculatrix
	 * args the a,b, and c
	 * pars, the other parameters q and stuff
	 * nargs, npars, sizes.
	 */
	//fptr
	double (*fptr)(double*,double*,int,int);
	void (*gptr)(double*,double*,int,int,double*);
	void (*hptr)(double*,double*,int,int,double*);
	fptr = &ThreePLModel::logLikelihood;
	gptr = &ThreePLModel::gradient;
	hptr = NULL;
	int I = model->getItemModel()->getDataset()->countItems();
	int q = model->getDimensionModel()->getLatentTraitSet()->getTheta()->nC();
	double args[3*I];
	double pars[2+2*q+q*I];
	int nargs = 3*I;
	int npars = 2+2*q+q*I;
	//filling args
	int nA=0;
	// Obtain a
	//A Matrix
	Matrix<double>* A = model->getParameterModel()->getParameterSet()[a];
	//B Matrix
	Matrix<double>* B = model->getParameterModel()->getParameterSet()[b];
	//C Matrix
	Matrix<double>* C = model->getParameterModel()->getParameterSet()[c];
		for (int i=0; i<I; i++) {
			args[nA] = (*A)(0,i);
			nA++;
		}

		// Obtain b
		for (int i=0; i<I; i++) {
			args[nA] = (*B)(0,i);
			nA++;
		}

		// Obtain c
		for (int i=0; i<I; i++) {
			args[nA] = (*C)(0,i);
			nA++;
		}
	//Filling pars
	int nP=0;
	// Obtain q
	pars[nP] = q;
	nP++;
	// Obtain I
	pars[nP] = I;
	nP++;
	// Obtain theta
	//Thetas
	Matrix<double>* thetas = model->getDimensionModel()->getLatentTraitSet()->getTheta();
	for (int k=0; k<q; k++) {
		pars[nP]=(*thetas)(0,q);//TODO correct indexing on this and nearby matrices
		nP++;
	}
	// Obtain f
	for (int k=0; k<q; k++) {
		pars[nP]=(*f)(0,k);
		nP++;
	}
	// Obtain r
	for (int k=0; k<q; k++) {
		for (int i=0; i<I; i++) {
			pars[nP]=(*r)(k,i);
			nP++;
		}
	}
	nargs = nA;
	npars = nP;
	optim->searchOptimal(fptr,gptr,hptr,args,pars,nargs,npars);
	 // Now pass the optimals to the Arrays.
	nA = 0;
	// Obtain a
	for (int i=0; i<I; i++) {
		(*A)(0,i) = args [nA ++];
	}
	// Obtain b
	for (int i=0; i<I; i++) {
		(*B)(0,i) = args [nA ++];
	}
	// Obtain c
	for (int i=0; i<I; i++) {
		(*C)(0,i) = args [nA ++];
	}

	//And set the parameter sets
	//model->parameterModel->
}
void EMEstimation::estimate(){
	/*
	 * TODO Estimate
	 */
}
//Check if all the conditions are met for running the model, can report an error to a logger
void EMEstimation::checkRunningConditions(){

}
	//Sets the optimization algorithm
void EMEstimation::setOptimizationAlgorithm(string algorithm){

}
	//Sets the reporter for the trace
void EMEstimation::setTrace(string filename){

}
void EMEstimation::setTrace(Trace trace){

}


