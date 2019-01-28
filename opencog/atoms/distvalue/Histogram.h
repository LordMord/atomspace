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

#include <opencog/atoms/distvalue/Interval.h>

namespace opencog
{
using namespace interval;

template<typename CT>
class Histogram
	: public Counter<NBin,CT>
{
public:
	DVec get_mean() const;
	bool has_bin(const NBin&) const;
	NBinSeq get_bins() const;

	double get_count(const NBin&) const;
	double get_contained_count(const NBin&) const;

	double get_frequency_for(const double&) const;
	double get_frequency(const NBin&) const;

	Interval minmax_count() const;

	Histogram<CT> merge(const Histogram<CT>&) const;

	std::string to_string(const std::string& indent) const;
};

}

#endif // _OPENCOG_HISTOGRAM_H
