/*
 * (c) Crown copyright 2009-2015.
 *
 * This software is distributed under the MIT License.
 * (See accompanying file LICENSE.txt or copy at
 *  http://opensource.org/licenses/MIT)
 */
/*
 * MutModel.h
 *
 *  Created on: Mar 22, 2012
 *      Author: gareth
 */

#ifndef MUTMODEL_H_
#define MUTMODEL_H_

struct MutModel
{
    enum Type { NO_MUTATIONS = 0, CRIME_IS_ANCESTOR, REF_IS_ANCESTOR };

    MutModel(Type t = NO_MUTATIONS, double rate = 0) : type(t), mutation_rate(rate) {}

    Type type;
    double mutation_rate;
};

#endif /* MUTMODEL_H_ */
