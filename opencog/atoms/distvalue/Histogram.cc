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

#include<sstream>
#include <opencog/util/exceptions.h>

#include <opencog/atoms/distvalue/Histogram.h>

using namespace opencog;
using namespace interval;

template class Histogram<double>;

template<typename CT>
DVec Histogram<CT>::get_mean() const
{
	DVec res;
	for (auto elem : *this)
		res = res + median(elem.first) * elem.second;
	return res / this->total_count();
}

template<typename CT>
bool Histogram<CT>::has_bin(const NBin &k) const
{
	auto it = this->find(k);
	return it != this->end();
}

template<typename CT>
NBinSeq Histogram<CT>::get_bins() const
{
	NBinSeq res(0);
	res.reserve(this->size());
	for (auto elem : *this)
		res.push_back(elem.first);
	return res;
}

template<typename CT>
double Histogram<CT>::get_count(const NBin &b) const
{
	auto it = this->find(b);
	if (it == this->end())
		throw RuntimeException(TRACE_INFO, "No Key for this value.");
	return it->second;
}

//Get the Count of a Bin that might not be in the DV explicitly
//by a weighted sum of all the Bins that are in the DV
//weighted by the overlapp of the given Bin with the Bins of the DV
template<typename CT>
double Histogram<CT>::get_contained_count(const NBin &b) const
{
	double res = 0;
	for (auto v : *this)
	{
		double weigth = conditional_probability(b,v.first);
		res += v.second * weigth;
	}
	return res;
}

template<typename CT>
double Histogram<CT>::get_frequency_for(const double& cnt) const
{
	return cnt / this->total_count();
}

template<typename CT>
double Histogram<CT>::get_frequency(const NBin &b) const
{
	return get_frequency_for(get_count(b));
}

//Get the lowest and highest count of all Interals
template<typename CT>
Interval Histogram<CT>::minmax_count() const
{
	double min = std::numeric_limits<double>::max();
	double max = 0;
	for (auto elem : *this)
	{
		if (min >= elem.second)
			min = elem.second;
		if (max <= elem.second)
			max = elem.second;
	}
	return Interval{min,max};
}

//Merge 2 Histograms into 1
template<typename CT>
Histogram<CT> Histogram<CT>::merge(const Histogram<CT> &other) const
{
	Histogram<CT> res;
	auto it1 = this->begin();
	auto it2 = other.begin();
	auto it1end = this->end();
	auto it2end = other.end();

	bool keep1 = false;
	bool keep2 = false;

	NBin bin1;
	NBin bin2;
	DVec lower1;
	DVec lower2;
	DVec upper1;
	DVec upper2;
	double count1;
	double count2;

	while (true)
	{
		if (keep1)
		{
			keep1 = false;
		}
		else
		{
			bin1 = it1->first;
			count1 = it1->second;
		}

		lower1 = lower(bin1);
		upper1 = upper(bin1);

		if (keep2)
		{
			keep2 = false;
		}
		else
		{
			bin2 = it2->first;
			count2 = it2->second;
		}

		lower2 = lower(bin2);
		upper2 = upper(bin2);

		if (upper1 < lower2)
		{
			res[bin1] = count1;
			it1++;
			continue;
		}
		if (lower1 > upper2)
		{
			res[bin2] = count2;
			it2++;
			continue;
		}
		NBin is = intersect(bin1,bin2);
		res[is] = conditional_probability(bin1,is) * count1
				+ conditional_probability(bin2,is) * count2;
		if (lower1 < lower2)
		{
			NBin lowerbin = createNBin(lower1,lower2);
			double prob = conditional_probability(bin1,lowerbin);
			if (prob != 0) res[lowerbin] = prob * count1;
		}
		else
		{
			NBin lowerbin = createNBin(lower2,lower1);
			double prob = conditional_probability(bin2,lowerbin);
			if (prob != 0) res[lowerbin] = prob * count2;
		}
		if (upper1 < upper2)
		{
			NBin upperbin = createNBin(upper1,upper2);
			count2 = conditional_probability(bin2,upperbin) * count2;
			it1++;
			if (it1 == it1end)
			{
				res[upperbin] = count2;
				for (; it2 != it2end;it2++)
					res[it2->first] = it2->second;
				return res;
			}
			else
			{
				bin2 = upperbin;
				keep2 = true;
			}
		}
		else
		{
			NBin upperbin = createNBin(upper2,upper1);
			count1 = conditional_probability(bin1,upperbin) * count1;
			it2++;
			if (it2 == it2end)
			{
				res[upperbin] = count1;
				for (it1++; it1 != it1end; it1++)
					res[it1->first] = it1->second;
				return res;
			}
			else
			{
				bin1 = upperbin;
				keep1 = true;
			}
		}
	}
}

template<typename CT>
std::string Histogram<CT>::to_string(const std::string& indent) const
{
	std::stringstream ss;
	if (this->size() == 0)
		ss << "Empty Histogram" << std::endl;
	for (auto elem : *this)
	{
		ss << interval::to_string(elem.first,indent)
		   << " Count: "
		   << elem.second
		   << " Frequency: "
		   << get_frequency_for(elem.second)
		   << std::endl;
	}
	return ss.str();
}
