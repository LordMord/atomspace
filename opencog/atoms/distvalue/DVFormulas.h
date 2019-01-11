/*
 * opencog/truthvalue/DVFormulas.h
 *
 * Copyright (C) 2018 SingularityNet
 * All Rights Reserved
 *
 * Written by Roman Treutlein
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License v3 as
 * published by the Free Software Foundation and including the exceptions
 * at http://opencog.org/wiki/Licenses
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program; if not, write to:
 * Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#ifndef _OPENCOG_DV_FORMULAS_H
#define _OPENCOG_DV_FORMULAS_H

#include <memory>
#include <string>
#include <vector>
#include <limits>

#include <opencog/atoms/distvalue/ConditionalDV.h>
#include <opencog/atoms/distvalue/DistributionalValue.h>

/** \addtogroup grp_atomspace
 *	@{
 */

namespace opencog
{

/**
 * This class contains Formulas for working with normal and conditional
 * DistributionalValues. These Formulas are mainly used in PLN.
 */
class DVFormulas
{

public:
	//Given a NdimBin for each Dimension/Interval it contains return it's
	//starting point(min) or ending point(max)
	static DVec get_key_min(const NdimBin&);
	static DVec get_key_max(const NdimBin&);

	//Given P(X,Y) and P(X) calculate P(Y|X)
	//using the formulae P(Y|X) = P(X,Y) / P(X)
	static ConditionalDVPtr joint_to_cdv(DistributionalValuePtr,
	                                     DistributionalValuePtr,
	                                     int);

	//Given P(X,Y) calculate P(X) for an idx of 0
	//                    or P(Y) for an idx of 1
	//using the formula P(X) = sum over all y P(X,y)
	static DistributionalValuePtr sum_joint(DistributionalValuePtr, int);


	//Given P(X) and P(Y) caluclate P(X & Y)
	//
	static DistributionalValuePtr conjunction(DistributionalValuePtr,
	                                          DistributionalValuePtr);
	static DistributionalValuePtr disjunction(DistributionalValuePtr,
	                                          DistributionalValuePtr);

	static ConditionalDVPtr consequent_disjunction_elemination(ConditionalDVPtr,
	                                                           ConditionalDVPtr);
};

} // namespace opencog

/** @}*/
#endif // _OPENCOG_TRUTH_VALUE_H
