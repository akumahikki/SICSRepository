/*
 * UnidimensionalModel.h
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#ifndef UNIDIMENSIONALMODEL_H_
#define UNIDIMENSIONALMODEL_H_

#include <model/dimension/DimensionModel.h>

class UnidimensionalModel : public DimensionModel {
public:
	// Constructor
	UnidimensionalModel();

	// Methods
	int getNumDimensions ();
	vector<double> getDimVector();

	// Getters and Setters
	LatentTraitSet* getLatentTraitSet() const;
	void setLatentTraitSet(LatentTraitSet* latentTraitSet);

	// Destructor
	virtual ~UnidimensionalModel();
};

#endif /* UNIDIMENSIONALMODEL_H_ */
