/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * SubPopModel.h
 *
 *  Created on: Mar 22, 2012
 *      Author: gareth
 */

#ifndef SUBPOPMODEL_H_
#define SUBPOPMODEL_H_

#include "cuda_accel/cuda_accel.h"

struct SubPopModel : public CudaSubPopModel
{
    SubPopModel() { type = HW; theta_bar = 0; }
    SubPopModel(Type t, double f=0) { type = t; theta_bar = f; }

    double prob(bool homozygous, double p1, double p2) const;
};

#endif /* SUBPOPMODEL_H_ */
