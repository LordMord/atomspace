/*
 * opencog/truthvalue/DistributionalValue.cc
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
#include <tuple>
#include <iomanip>

#include <opencog/atoms/distvalue/DistributionalValue.h>
#include <opencog/atomspace/AtomSpace.h>

#include <boost/range/combine.hpp>

using namespace opencog;

count_t DistributionalValue::DEFAULT_K = 800.0;

DistributionalValue::DistributionalValue()
	: Value(DISTRIBUTIONAL_VALUE) , _value(15,1)
{
}

DistributionalValue::DistributionalValue(const CHist<double> &hist)
	: Value(DISTRIBUTIONAL_VALUE) //, _value(hist)
{
	_value = hist;
}

//Create a DV from the Parameters of a SimpleTV
//Not recommended as it results in a DV with only 1 singleton set
//which in some calculations is unusalbe as the chance of overlaps in the
//Sets/Bins is small
DistributionalValue::DistributionalValue(double mode,double conf)
	: Value(DISTRIBUTIONAL_VALUE) , _value(1,1)
{
	confidence_t cf = std::min(conf, 0.9999998);
	double count = (DEFAULT_K * cf / (1.0 - cf));
	//DV Bin with count 0 is undefined
	count = std::max(count, 0.0000002);
	_value.insert(DVec{mode},count);
}

DistributionalValuePtr DistributionalValue::createDV(const CHist<double> &hist)
{
	auto res = std::make_shared<const DistributionalValue>(hist);
	return res;
}

DistributionalValuePtr DistributionalValue::createDV(double mode,
                                                     double conf)
{
	return std::make_shared<const DistributionalValue>(mode,conf);
}

DistributionalValuePtr DistributionalValue::TRUE_TV()
{
	static DistributionalValuePtr instance;
	if (instance == nullptr)
	{
		CHist<double> hist = CHist<double>(1,1);
		hist.insert(DVec{1.0},1.0);
		instance = std::make_shared<const DistributionalValue>(hist);
	}
    return instance;
}
DistributionalValuePtr DistributionalValue::FALSE_TV()
{
	static DistributionalValuePtr instance;
	if (instance == nullptr)
	{
		CHist<double> hist = CHist<double>(1,1);
		hist.insert(DVec{0.0},1.0);
		instance = std::make_shared<const DistributionalValue>(hist);
	}
    return instance;
}
DistributionalValuePtr DistributionalValue::DEFAULT_TV()
{
	static DistributionalValuePtr instance;
	if (instance == nullptr)
	{
		//TODO: Should this have a size of 15?
		CHist<double> hist = CHist<double>(0,1);
		instance = std::make_shared<const DistributionalValue>(hist);
	}
    return instance;
}

//Add Evidence for the provided key i.e increment the count of this Bin
//Returns a new DVPtr with the update count
void DistributionalValue::add_evidence(const DVec& pos)
{
	_value.insert(pos,1.0);
}

//Combines this DV with the provided DV
//and returns this new one ass a result
DistributionalValuePtr DistributionalValue::merge(DistributionalValuePtr other) const
{
	CHist<double> hist = CHist<double>::merge(_value,other->_value);
	return createDV(hist);
}

std::vector<double> DistributionalValue::bin_modes() const
{
	std::vector<double> probs;
	for (auto elem : _value)
	{
		probs.push_back(get_mode_for(elem.value));
	}
	return probs;
}

std::vector<double> DistributionalValue::bin_means() const
{
	std::vector<double> probs;
	for (auto elem : _value)
	{
		probs.push_back(get_mean_for(elem.value));
	}
	return probs;
}

//Get the variance for all Keys
std::vector<double> DistributionalValue::bin_vars() const
{
	std::vector<double> probs;
	for (auto elem : _value)
	{
		probs.push_back(get_var_for(elem.value));
	}
	return probs;
}

double DistributionalValue::get_mode_for(double ai) const
{
	return (ai - 1) / (total_count() - _value.count_bins());
}

double DistributionalValue::get_mean_for(double ai) const
{
	if (total_count() == 0)
		return 0;
	else
		return ai / total_count();
}

//Get the variance for a Count
double DistributionalValue::get_var_for(double ai) const
{
		double a0 = total_count();
		return ai*(a0-ai) / a0*a0*(a0+1);
}

double DistributionalValue::total_count() const
{
	return _value.total_count();
}

double DistributionalValue::get_confidence() const
{
	int c = total_count();
	return to_conf(c);
}

double DistributionalValue::to_conf(int c)
{
	return c / (c + DEFAULT_K);
}

int DistributionalValue::to_count(double cf)
{
	return (cf * DEFAULT_K / (1 - cf));
}

double DistributionalValue::get_mode(const DVec &val) const
{
	return get_mode_for(_value.get(val));
}
double DistributionalValue::get_mean(const DVec &val) const
{
	return get_mean_for(_value.get(val));
}
double DistributionalValue::get_contained_mean(const DVec &val) const
{
	return get_mean_for(_value.get_avg(val));
}
double DistributionalValue::get_var(const DVec &val) const
{
	return get_var_for(_value.get(val));
}

std::string DistributionalValue::to_string(const std::string& indent) const
{
	std::stringstream ss;
	if (_value.size() == 0)
		ss << "Empty DistributionalValue" << std::endl;

	for (auto elem : _value)
	{
		ss << Node<double>::to_string(_value,elem)
		   << " Mean: "
		   << get_mean_for(elem.value)
		   << std::endl;
	}
	return ss.str();
}

bool DistributionalValue::operator==(const Value& other) const
{
	if (DISTRIBUTIONAL_VALUE != other.get_type()) return false;
	const DistributionalValue* dov = (const DistributionalValue*) &other;

	return _value == dov->_value;
}

std::string oc_to_string(const DistributionalValuePtr& dvp,
                         const std::string& indent)
{
	return dvp->to_string(indent);
}
