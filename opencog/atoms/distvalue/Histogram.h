/*
 * opencog/truthvalue/DistributionalValue.cc
 *
 * Copyright (C) 2019 SingularityNet
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


#ifndef _OPENCOG_HISTOGRAM_H
#define _OPENCOG_HISTOGRAM_H

#include <vector>

#include <opencog/util/Counter.h>

#include <boost/numeric/interval.hpp>

using namespace boost::numeric;
using namespace interval_lib;

namespace opencog
{

typedef std::vector<double> DVec;
//A Left-open Interval
//struct Interval
//{
//	double left;
//	double right;
//};
//typedef Interval Bin;
//A N-dimensional Bin
//typedef std::vector<Bin> NBin;
typedef interval<DVec> NBin;

//A Sequence of N-dimensional Bins
typedef std::vector<NBin> NBinSeq;

template<typename CT>
class Histogram
	: public Counter<NBin,CT>
{
public:
	static DVec center_of_bin(const NBin&);

	DVec get_mean() const;
	bool has_bin(const NBin&) const;
	NBinSeq get_bins() const;
};

}

#endif // _OPENCOG_HISTOGRAM_H
