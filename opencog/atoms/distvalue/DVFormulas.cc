/*
 * opencog/truthvalue/DVFormulas.cc
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

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <iomanip>
#include <vector>

#include <opencog/util/numeric.h>
#include <opencog/atoms/distvalue/DVFormulas.h>
#include <opencog/atomspace/AtomSpace.h>

using namespace opencog;

#if 0

//Given a key return the min of all it's Intervals
DVec DVFormulas::get_key_min(const NdimBin &k)
{
	std::vector<double> res;
	for (auto interval : k)
		res.push_back(interval[0]);
	return res;
}

//Given a key return the max of all it's Intervals
DVec DVFormulas::get_key_max(const NdimBin &k)
{
	std::vector<double> res;
	for (auto interval : k)
		if (interval.size() == 1)
			res.push_back(interval[0]);
		else
			res.push_back(interval[1]);

	return res;
}

#endif

//(A,B,C) + (B,C) => (B,C) -> A
//idx is the position of consequent in the joint Distribution
//(A,B,C) A is at idx 0
//Result has the same total_count as dv2
ConditionalDVPtr DVFormulas::joint_to_cdv(DistributionalValuePtr dv1,
                                          DistributionalValuePtr dv2,
                                          uint idx)
{
	size_t dv1dims = dv1->value().dimensions();
	size_t dv2dims = dv2->value().dimensions();

	size_t dv1size = dv1->value().size();
	size_t dv2size = dv2->value().size();

	CDVrep res = CDVrep(dv1size / dv2size,dv2dims);

	if (dv1dims <= 1)
		throw RuntimeException(TRACE_INFO,"Can't divide non Joint DV.");
	if (dv1dims - 1 != dv2dims)
		throw RuntimeException(TRACE_INFO,"The Divisor DV has to have exaclty 1 less dimensions then then dividend. This is not the case.");

	for (auto elem : dv1->value())
	{
		auto size = sizeof(double);
		double * hs = (double*)malloc(size * dv2dims);

		//get the consequent
		double * h = (double*)malloc(size);
		memcpy(h,elem.pos+idx,size);
		//remove it from the condition
		int j = 0;
		for (uint i = 0; i < dv2dims; i++)
		{
			if (i != idx)
			{
				memcpy(hs+j,elem.pos+i,size);
				j++;
			}
		}

		DVec hsv = res.arrayToVec(hs);

		if (dv2->get_contained_mean(hsv) != 0)
		{
			double count = dv1->get_mean_for(elem.value) /
							dv2->get_contained_mean(hsv) *
							dv2->total_count();

			CHist<double> val = CHist<double>(dv2size,1);
			val.insert(h,count);
			res.insert(hs,val);
		}

	}
	return ConditionalDV::createCDV(res);
}

#if 0

//(A,B,C) => (A,C)
//idx is the position of the Element to sum out of the joint-dv
DistributionalValuePtr DVFormulas::sum_joint(DistributionalValuePtr dv, int pos)
{
	DVCounter res;
	for (auto elem : dv->value())
	{
		NdimBin key = elem.first;
		key.erase(key.begin() + pos);
		res[key] += elem.second;
	}
	return DistributionalValue::createDV(res);
}

#define EPSILON 1e-16

//Create a Conjuction from 2 DVs
DistributionalValuePtr
DVFormulas::conjunction(DistributionalValuePtr dv1,
                        DistributionalValuePtr dv2)
{
	DVCounter res;
	double count = std::min(dv1->total_count(),dv2->total_count());

	//We start at the begining of the map with Keys of the lowest Value
	auto it1 = dv1->value().begin();
	auto it2 = dv2->value().begin();

	//Weighting factor representing how much of a given DV has be used already
	double m1 = 1;
	double m2 = 1;

	while (not is_within(m1,0.0,EPSILON) && not is_within(m2,0.0,EPSILON))
	{
		DVec v1 = get_key_min(it1->first);
		DVec v2 = get_key_min(it2->first);
		//We check which key represents a lower Truthness/Value
		//This is a fuzzy conjunction so we want to take the min of that
		if (v1 < v2)
		{
			//We get the mean for that Key
			double mean = dv1->get_mean_for(it1->second);
			//Multiply by the count to de-normalize that
			//Weighted by how much of the other DV we already used
			res[it1->first] += count * mean * m2;
			//Update m1 to refelect that we have used "mean" of this DV
			m1 -= mean;
			it1++;
		}
		else
		{
			//Same as above just flipped dv1/it1/m1 and dv2/it2/m2
			double mean = dv2->get_mean_for(it2->second);
			res[it2->first] += count * mean * m1;
			m2 -= mean;
			it2++;
		}
	}

	return DistributionalValue::createDV(res);
}

//Create a disjuction from 2 DVs
DistributionalValuePtr
DVFormulas::disjunction(DistributionalValuePtr dv1,
                        DistributionalValuePtr dv2)
{
	DVCounter res;
	double count = std::min(dv1->total_count(),dv2->total_count());

	//We start at the end of the map with Keys of the highest Value
	auto it1 = dv1->value().end();
	auto it2 = dv2->value().end();

	it1--;
	it2--;

	//Weighting factor representing how much of a given DV has be used already
	double m1 = 1;
	double m2 = 1;

	while (not is_within(m1,0.0,EPSILON) && not is_within(m2,0.0,EPSILON))
	{
		if (m1 < 0 || m2 < 0)
			throw RuntimeException(TRACE_INFO,"This should not happen.");
		DVec v1 = get_key_max(it1->first);
		DVec v2 = get_key_max(it2->first);
		//We check which key represents a higher Truthness/Value
		//This is a fuzzy disjunction so we want to take the max of that
		if (v1 > v2)
		{
			double mean = dv1->get_mean_for(it1->second);
			//Multiply by the count to de-normalize that
			//Weighted by how much of the other DV we already used
			res[it1->first] += count * mean * m2;
			//Update m1 to refelect that we have used "mean" of this DV
			m1 -= mean;
			it1--;
		}
		else
		{
			//Same as above just flipped dv1/it1/m1 and dv2/it2/m2
			double mean = dv2->get_mean_for(it2->second);
			res[it2->first] += count * mean * m1;
			m2 -= mean;
			it2--;
		}
	}

	return DistributionalValue::createDV(res);
}

//A -> (B || C) + A -> B => A -> C
ConditionalDVPtr
DVFormulas::consequent_disjunction_elemination(ConditionalDVPtr cdv1,
                                               ConditionalDVPtr cdv2)
{
	CDVrep res;
	for (auto elem : cdv2->value())
	{
		DistributionalValuePtr v1 = cdv1->get_unconditional(elem.first);
		DistributionalValuePtr v2 = DistributionalValue::createDV(elem.second);

		double count = std::min(v1->total_count(),v2->total_count());

		DVCounter partres;
		for (auto elem2 : elem.second)
		{
			double m1 = v1->get_contained_mean(elem2.first);
			double m2 = v2->get_mean(elem2.first);
			partres[elem2.first] = count * ((m1 - m2) / (1 - m2));
		}
		res[elem.first] = partres;
	}
	return ConditionalDV::createCDV(res);
}

#endif