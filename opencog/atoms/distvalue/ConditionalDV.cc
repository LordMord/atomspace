/*
 * opencog/atoms/distvalue/ConditionalDV.cc
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
#include <sstream>

#include <opencog/atoms/distvalue/ConditionalDV.h>
#include <opencog/atoms/distvalue/DistributionalValue.h>

using namespace opencog;

ConditionalDV::ConditionalDV()
	: Value(CONDITIONAL_DISTRIBUTIONAL_VALUE) , _value(0,0)
{}

ConditionalDV::ConditionalDV(const CDVrep &rep)
	: Value(CONDITIONAL_DISTRIBUTIONAL_VALUE) , _value(rep)
{}

ConditionalDV::ConditionalDV(const DVecSeq &conds,
                             const std::vector<DistributionalValuePtr> &dvs)
	: Value(CONDITIONAL_DISTRIBUTIONAL_VALUE) , _value(conds.size(),conds[0].size())
{
	auto it1 = conds.begin();
	auto it2 = dvs.begin();
	auto end1	= conds.end();
	auto end2	= dvs.end();
	if (conds.empty())
		throw RuntimeException(TRACE_INFO,"Conds may not be empty.");

	if (dvs.empty())
		throw RuntimeException(TRACE_INFO,"DVs may not be empty.");

	if (dvs.size() != conds.size())
		throw RuntimeException(TRACE_INFO,"DVs and Conds must be the same lenght.");

	for (;(it1 != end1) && (it2 != end2); ++it1, ++it2)
	{
		_value.insert(*it1,(*it2)->_value);
	}
}

ConditionalDVPtr ConditionalDV::createCDV()
{
	return std::make_shared<const ConditionalDV>();
}

ConditionalDVPtr ConditionalDV::createCDV(const CDVrep &rep)
{
	return std::make_shared<const ConditionalDV>(rep);
}

ConditionalDVPtr ConditionalDV::createCDV(const DVecSeq &conds,
                                          const std::vector<DistributionalValuePtr> &dvs)
{
	return std::make_shared<const ConditionalDV>(conds,dvs);
}

//Get all the Conditions
DVecSeq ConditionalDV::get_conditions() const
{
	return _value.get_posvec();
}

//Get all the DVs without the Conditions
std::vector<DistributionalValuePtr> ConditionalDV::get_unconditionals() const
{
	std::vector<DistributionalValuePtr> res;
	for (auto node : _value)
	{
		res.push_back(DistributionalValue::createDV(node.value));
	}
	return res;
}

ConditionalDVPtr ConditionalDV::remap(const DVecSeq &k) const
{
	CDVrep rep = _value.remap(k);
	return createCDV(rep);
}

/*
 * Get a DV that is the weighted combination of all contained DVs based
 * on a Distribution of the Condition
 */
DistributionalValuePtr ConditionalDV::get_unconditional(DistributionalValuePtr condDist) const
{
	auto size = std::pow(2,std::log2(_value.size()+1) +
	                     std::log2(_value[0].value.size()+1)) - 1;
	auto dims = _value[0].value.dims();
	CTHist<double> res = CTHist<double>(size,dims);
	double sum = 0;

	DVecSeq keys = condDist->_value.get_posvec();
	CDVrep remaped = _value.remap(keys);

	for (auto v : condDist->_value)
	{
		double val = condDist->get_mean_for(v.value);
		sum += total_count() / remaped.get(v.pos).total_count();
		res += remaped.get(v.pos) * val;
	}
	res /= sum;
	return std::make_shared<const DistributionalValue>(res);
}

double ConditionalDV::total_count() const
{
	return _value.total_count();
}

double ConditionalDV::avg_count() const
{
	double res = 0;
	int count = 0;
	for (auto elem : _value)
	{
		count++;
		res += elem.value.total_count();
	}
	return res / count;
}

//Given a Distribution of the Condition calculate a Joint Probability distribution
DistributionalValuePtr ConditionalDV::get_joint_probability(DistributionalValuePtr base) const
{
	auto s1 = base->_value.size();
	auto d1 = base->_value.dims();
	auto s2 = _value.begin()->value.size();
	auto d2 = _value.begin()->value.dims();
	CTHist<double> res = CTHist<double>(s1*s2,d1+d2);
	DVecSeq ivsBASE = base->_value.get_posvec();
	//std::cout << base->to_string();
	//std::cout << "ivsBase.size: " << ivsBASE.size() << std::endl;

	ConditionalDVPtr remaped = remap(ivsBASE);

	for (auto k1 : ivsBASE) {
		DistributionalValuePtr uncond = DistributionalValue::createDV(remaped->value().get(k1));
		DVecSeq ivsTHIS = uncond->_value.get_posvec();
		//uncond->value().print();
		//std::cout << "ivsThis.size: " << ivsTHIS.size() << std::endl;
		for (auto k2 : ivsTHIS) {
			DVec k;
			k.insert(k.end(),k1.begin(),k1.end());
			k.insert(k.end(),k2.begin(),k2.end());

			//Res count based on base count
			auto v1 = base->_value.get(k1);
			auto v2 = uncond->get_mean(k2);
			std::cout << "v1: " << v1 << std::endl;
			std::cout << "v2: " << v2 << std::endl;
			//std::cout << "Inserting: " << ::to_string(k) << std::endl;
			res.insert(k,v1 * v2);
			//res.print();
			//std::cout << std::endl;
		}
	}
	return DistributionalValue::createDV(res);
}

// A->C + B->C => (A,B)->C
ConditionalDVPtr ConditionalDV::merge(ConditionalDVPtr cdv2) const
{
	if (_value.dims() != cdv2->_value.dims())
		throw RuntimeException(TRACE_INFO,"Can't merge CDVs with different number of dimensions.");

	CDVrep res = CDVrep(std::max(_value.size(),cdv2->_value.size()),_value.dims());
	for (auto elem1 : _value)
	{
		for (auto elem2 : cdv2->_value)
		{
			DVec k;
			k.insert(k.end(),elem1.pos.begin(),elem1.pos.end());
			k.insert(k.end(),elem2.pos.begin(),elem2.pos.end());
			CTHist<double> tmp = CTHist<double>::merge(elem1.value,elem2.value);
			res.insert(k,tmp);
		}
	}
	return createCDV(res);
}


std::string ConditionalDV::to_string(const std::string& indent) const
{
	std::stringstream ss;
	if (_value.size() == 0)
		ss << "Empty ConditionalDV" << std::endl;
	for (auto elem : _value)
	{
		ss << indent << "{";
		for (double coord : elem.pos)
		{
				ss << coord << ";";
		}
		ss.seekp(-1,std::ios_base::end);
		ss << "} DV: "
		   << std::endl
		   << DistributionalValue(elem.value).to_string(indent + "    ")
		   << std::endl;
	}
	return ss.str();
}

bool ConditionalDV::operator==(const Value& other) const
{
	if (CONDITIONAL_DISTRIBUTIONAL_VALUE != other.get_type()) return false;
	const ConditionalDV* cov = (const ConditionalDV*) &other;

	if (_value != cov->_value)
		return false;
	return true;
}
