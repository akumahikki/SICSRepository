/*
 * MultiUniDimModel.h
 *
 *  Created on: 18 Jun 2014
 *      Author: jlgpisa
 */

#ifndef MULTIUNIDIMMODEL_H_
#define MULTIUNIDIMMODEL_H_

#include <model/dimension/DimensionModel.h>

class MultiUniDimModel: public DimensionModel {
public:
	// Constructor
	MultiUniDimModel();

	// Methods
	int getNumDimensions ();
	vector<double> getDimVector();

	// Getters and Setters
	const LatentTraitSet* getLatentTraitSet() const;
	void setLatentTraitSet(LatentTraitSet* latentTraitSet);

	// Destructor
	virtual ~MultiUniDimModel();
};

#endif /* MULTIUNIDIMMODEL_H_ */
