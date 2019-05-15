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

//(A,B,C) + (B,C) => (B,C) -> A
//idx is the position of consequent in the joint Distribution
//(A,B,C) A is at idx 0
//Result has the same total_count as dv2
ConditionalDVPtr DVFormulas::joint_to_cdv(DistributionalValuePtr dv1,
                                          DistributionalValuePtr dv2,
                                          int idx)
{
	size_t dv1dims = dv1->value().dims();
	size_t dv2dims = dv2->value().dims();

	size_t dv1size = dv1->value().size();
	size_t dv2size = dv2->value().size();

	CDVrep res = CDVrep(dv1size / dv2size,dv2dims);

	if (dv1dims <= 1)
		throw RuntimeException(TRACE_INFO,"Can't divide non Joint DV.");
	if (dv1dims - 1 != dv2dims)
		throw RuntimeException(TRACE_INFO,"The Divisor DV has to have exaclty 1 less dimensions then then dividend. This is not the case.");


	DVecSeq keys = dv1->value().get_posvec();
	for (unsigned int i = 0; i < keys.size(); i++)
	{
		keys[i].erase(keys[i].begin() + idx);
	}
	std::sort(keys.begin(),keys.end());
	keys.erase(std::unique(keys.begin(),keys.end()),keys.end());
	DistributionalValuePtr dv2remap = dv2->remap(keys);

	for (auto elem : dv1->value())
	{
		DVec hs = elem.pos;

		DVec h = DVec{hs[idx]};
		hs.erase(hs.begin() + idx);

		std::cout << ::to_string(hs) << " hs,h " << ::to_string(h) << std::endl;
		//std::cout << "dv2m: " << dv2->get_contained_mean(hs) << std::endl;
		//Counter<DVec,double> means = dv2->get_contained_means(hs);
		if (dv2remap->get_mean(hs) != 0)
		{
			double count = dv1->get_mean_for(elem.value) /
							dv2remap->get_mean(hs) *
							dv2remap->total_count();

			//std::cout << "dv1m: " << dv1->get_mean_for(elem.value) << std::endl;
			//std::cout << "dv2m: " << dv2->get_contained_mean(hs) << std::endl;
			//std::cout << "dv2t: " << dv2->total_count() << std::endl;

			CTHist<double> val = CTHist<double>(dv2size,1);
			val.insert(h,count);
			res.insert(hs,val);
			//res.print();
			//std::cout << std::endl;
		}
	}
	return ConditionalDV::createCDV(res);
}

//(A,B,C) => (A,C)
//idx is the position of the Element to sum out of the joint-dv
DistributionalValuePtr DVFormulas::sum_joint(DistributionalValuePtr dv, int pos)
{
	CTHist<double> res = CTHist<double>(dv->value().size(),dv->value().dims()-1);
	for (auto elem : dv->value())
	{
		DVec key = elem.pos;
		key.erase(key.begin() + pos);
		res.insert(key,elem.value);
	}
	return DistributionalValue::createDV(res);
}


#define EPSILON 1e-16

//Create a Conjuction from 2 DVs
DistributionalValuePtr
DVFormulas::conjunction(DistributionalValuePtr dv1,
                        DistributionalValuePtr dv2)
{
	CTHist<double> res = CTHist<double>(dv1->value().size() + dv2->value().size(),
										dv1->value().dims());
	double count = std::min(dv1->total_count(),dv2->total_count());

	const std::vector<CoverTreeNode<double>> & nodes1 = dv1->value().nodes();
	const std::vector<CoverTreeNode<double>> & nodes2 = dv2->value().nodes();

	std::vector<int> idxs1(nodes1.size());
	std::size_t c1(0);
    std::generate(std::begin(idxs1), std::end(idxs1), [&]{ return c1++; });

	std::vector<int> idxs2(nodes2.size());
	std::size_t c2(0);
    std::generate(std::begin(idxs2), std::end(idxs2), [&]{ return c2++; });

    auto cmp1 = [&](int i1,int i2) -> bool
				{
				   return nodes1[i1].pos < nodes1[i2].pos;
				};
    auto cmp2 = [&](int i1,int i2) -> bool
				{
				   return nodes2[i1].pos < nodes2[i2].pos;
				};

	std::sort(idxs1.begin(),idxs1.end(),cmp1);
	std::sort(idxs2.begin(),idxs2.end(),cmp2);

	//We start at the begining of the map with Keys of the lowest Value
	auto it1 = idxs1.begin();
	auto it2 = idxs2.begin();
	const CoverTreeNode<double> * n1 = &nodes1[*it1];
	const CoverTreeNode<double> * n2 = &nodes2[*it2];

	//Weighting factor representing how much of a given DV has be used already
	double m1 = 1;
	double m2 = 1;

	while (not is_within(m1,0.0,EPSILON) && not is_within(m2,0.0,EPSILON))
	{
		//We check which key represents a lower Truthness/Value
		//This is a fuzzy conjunction so we want to take the min of that
		if (n1->pos < n2->pos)
		{
			//We get the mean for that Key
			double mean = dv1->get_mean_for(n1->value);
			//Multiply by the count to de-normalize that
			//Weighted by how much of the other DV we already used
			res.insert(n1->pos, count * mean * m2);
			//Update m1 to refelect that we have used "mean" of this DV
			m1 -= mean;
			it1++;
			n1 = &nodes1[*it1];
		}
		else
		{
			//Same as above just flipped dv1/it1/m1 and dv2/it2/m2
			double mean = dv2->get_mean_for(n2->value);
			res.insert(n2->pos, count * mean * m1);
			m2 -= mean;
			it2++;
			n2 = &nodes2[*it2];
		}
	}

	return DistributionalValue::createDV(res);
}

//Create a disjunction from 2 DVs
DistributionalValuePtr
DVFormulas::disjunction(DistributionalValuePtr dv1,
                        DistributionalValuePtr dv2)
{
	CTHist<double> res = CTHist<double>(dv1->value().size() + dv2->value().size(),
										dv1->value().dims());
	double count = std::min(dv1->total_count(),dv2->total_count());

	const std::vector<CoverTreeNode<double>> & nodes1 = dv1->value().nodes();
	const std::vector<CoverTreeNode<double>> & nodes2 = dv2->value().nodes();

	std::vector<int> idxs1(nodes1.size());
	std::size_t c1(0);
    std::generate(std::begin(idxs1), std::end(idxs1), [&]{ return c1++; });

	std::vector<int> idxs2(nodes2.size());
	std::size_t c2(0);
    std::generate(std::begin(idxs2), std::end(idxs2), [&]{ return c2++; });

    auto cmp1 = [&](int i1,int i2) -> bool
				{
				   return nodes1[i1].pos < nodes1[i2].pos;
				};
    auto cmp2 = [&](int i1,int i2) -> bool
				{
				   return nodes2[i1].pos < nodes2[i2].pos;
				};

	std::sort(idxs1.begin(),idxs1.end(),cmp1);
	std::sort(idxs2.begin(),idxs2.end(),cmp2);

	//We start at the begining of the map with Keys of the lowest Value
	auto it1 = idxs1.end() - 1;
	auto it2 = idxs2.end() - 1;
	const CoverTreeNode<double> * n1 = &nodes1[*it1];
	const CoverTreeNode<double> * n2 = &nodes2[*it2];

	//Weighting factor representing how much of a given DV has be used already
	double m1 = 1;
	double m2 = 1;

	while (not is_within(m1,0.0,EPSILON) && not is_within(m2,0.0,EPSILON))
	{
		//We check which key represents a lower Truthness/Value
		//This is a fuzzy conjunction so we want to take the min of that
		if (n1->pos > n2->pos)
		{
			//We get the mean for that Key
			double mean = dv1->get_mean_for(n1->value);
			//Multiply by the count to de-normalize that
			//Weighted by how much of the other DV we already used
			res.insert(n1->pos, count * mean * m2);
			//Update m1 to refelect that we have used "mean" of this DV
			m1 -= mean;
			it1--;
			n1 = &nodes1[*it1];
		}
		else
		{
			//Same as above just flipped dv1/it1/m1 and dv2/it2/m2
			double mean = dv2->get_mean_for(n2->value);
			res.insert(n2->pos, count * mean * m1);
			m2 -= mean;
			it2--;
			n2 = &nodes2[*it2];
		}
	}

	return DistributionalValue::createDV(res);
}

#if 0
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
