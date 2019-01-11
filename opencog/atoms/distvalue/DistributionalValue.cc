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

#include <opencog/util/numeric.h>

#include <opencog/atoms/distvalue/DistributionalValue.h>
#include <opencog/atomspace/AtomSpace.h>

#include <boost/range/combine.hpp>

using namespace opencog;

count_t DistributionalValue::DEFAULT_K = 800.0;

DistributionalValue::DistributionalValue()
	: Value(DISTRIBUTIONAL_VALUE)
{}

DistributionalValue::DistributionalValue(const Histogram<double> &dvctr)
	: Value(DISTRIBUTIONAL_VALUE)
{
	for (auto elem : dvctr)
	{
		if (elem.second <= 0)
			throw std::invalid_argument("DV cant have an element with count 0.");
	}
	_value = dvctr;
}

//Create a DV from the Parameters of a SimpleTV
//Not recommended as it results in a DV with only 1 singleton set
//which in some calculations is unusalbe as the chance of overlaps in the
//Sets/Bins is small
DistributionalValue::DistributionalValue(double mode,double conf)
	: Value(DISTRIBUTIONAL_VALUE)
{
	confidence_t cf = std::min(conf, 0.9999998);
	double count = (DEFAULT_K * cf / (1.0 - cf));
	//DV Bin with count 0 is undefined
	count = std::max(count, 0.0000002);
	NBin k{Bin{mode,mode}};
	_value[k] = count;
}

DistributionalValuePtr DistributionalValue::createDV(const Histogram<double> &dvctr)
{
	return std::make_shared<const DistributionalValue>(dvctr);
}

DistributionalValuePtr DistributionalValue::createDV(double mode,
                                                     double conf)
{
	return std::make_shared<const DistributionalValue>(mode,conf);
}

DistributionalValuePtr
DistributionalValue::UniformDistributionalValue(const NBin &k,int c)
{
	Histogram<double> dvctr;
	dvctr[k] = c;
	return createDV(dvctr);
}


DistributionalValuePtr
DistributionalValue::UniformDistributionalValue(const NBinSeq &ks,int c)
{
	Histogram<double> dvctr;
	for (auto k : ks)
	{
		dvctr[k] = c;
	}
	return createDV(dvctr);
}

DistributionalValuePtr DistributionalValue::TRUE_TV()
{
	static DistributionalValuePtr instance;
	if (instance == nullptr)
	{
		NBin v1{Bin{1.0}};
		NBin v2{Bin{0.0}};
		Histogram<double> dvc;
		dvc[v1] = 1;
		dvc[v2] = 0;
		instance = std::make_shared<const DistributionalValue>(dvc);
	}
    return instance;
}
DistributionalValuePtr DistributionalValue::FALSE_TV()
{
	static DistributionalValuePtr instance;
	if (instance == nullptr)
	{
		NBin v1{Bin{1.0}};
		NBin v2{Bin{0.0}};
		Histogram<double> dvc;
		dvc[v1] = 0;
		dvc[v2] = 1;
		instance = std::make_shared<const DistributionalValue>(dvc);
	}
    return instance;
}
DistributionalValuePtr DistributionalValue::DEFAULT_TV()
{
	static DistributionalValuePtr instance;
	if (instance == nullptr)
	{
		NBin v1{Bin{1.0}};
		NBin v2{Bin{0.0}};
		Histogram<double> dvc;
		dvc[v1] = 0;
		dvc[v2] = 0;
		instance = std::make_shared<const DistributionalValue>(dvc);
	}
    return instance;
}

bool DistributionalValue::is_uniform() const
{
	double val = _value.begin()->second;
	for (auto p : _value)
	{
		if (val != p.second)
		{
			return false;
		}
	}
	return true;
}

//Add Evidence for the provided key i.e increment the count of this Bin
//Returns a new DVPtr with the update count
DistributionalValuePtr DistributionalValue::add_evidence(const NBin &h) const
{
	Histogram<double> newdvc = _value;
	newdvc[h] += 1;
	return createDV(newdvc);
}

//Combines this DV with the provided DV
//and returns this new one ass a result
DistributionalValuePtr DistributionalValue::merge(DistributionalValuePtr other) const
{
	Histogram<double> newdvc = _value;
	newdvc += other->_value;
	return createDV(newdvc);
}

//Flip all the counts
DistributionalValuePtr DistributionalValue::negate() const
{
	Interval minmax = minmax_count();
	double total = minmax.left + minmax.right;
	Histogram<double> res;
	for (auto elem : _value)
	{
		res[elem.first] = total - elem.second;
	}
	return createDV(res);
}

//Get the lowest and highest count of all Interals
Interval DistributionalValue::minmax_count() const
{
	Interval minmax = Bin{std::numeric_limits<double>::max(),0};
	for (auto elem : _value)
	{
		if (minmax.left >= elem.second)
			minmax.left = elem.second;
		if (minmax.right <= elem.second)
			minmax.right = elem.second;
	}
	return minmax;
}

std::vector<double> DistributionalValue::bin_modes() const
{
	std::vector<double> probs;
	for (auto elem : _value)
	{
		probs.push_back(get_mode_for(elem.second));
	}
	return probs;
}

double DistributionalValue::get_mode_for(double ai) const
{
	return (ai - 1) / (total_count() - _value.size());
}

std::vector<double> DistributionalValue::bin_means() const
{
	std::vector<double> probs;
	for (auto elem : _value)
	{
		probs.push_back(get_mean_for(elem.second));
	}
	return probs;
}

double DistributionalValue::get_mean_for(double ai) const
{
	double count = total_count();
	if (count == 0)
		return 0;
	else
		return ai / total_count();
}

//Condenses the DV into a single strenght value
//This should be avoided where ever possible as it loses a lot of information
double DistributionalValue::get_fstord_mean() const
{
	double res = 0;
	for (auto elem : _value)
	{
		//TODO: Is mode correct or should i use mean?
		DVec mof = Histogram<double>::center_of_bin(elem.first);
		double mode = std::accumulate(mof.begin(),mof.end(),0.0)/mof.size();
		res = res + mode * get_mode_for(elem.second);
	}
	return res;
}

//Get the variance for all Keys
std::vector<double> DistributionalValue::bin_vars() const
{
	std::vector<double> probs;
	for (auto elem : _value)
	{
		probs.push_back(get_var_for(elem.second));
	}
	return probs;
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

double DistributionalValue::get_count(const NBin &h) const
{
	auto pos = _value.find(h);
	if (pos != _value.end())
		return pos->second;
	throw RuntimeException(TRACE_INFO, "No Key for this value.");
}

//Find out how much of Key1 is contained in Key2
//The Keys can be considered continuous uniform distributions
//This calculates the conditional_probabilty of Key2 given Key1
double DistributionalValue::conditional_probabilty(const NBin &ks1,const NBin &ks2)
{
	Bin k1,k2;
	//Start with the assumption that 100% of Key1 is in Key2
	double sum = 1;
	for (auto zipped : boost::combine(ks1,ks2))
	{
		boost::tie(k1,k2) = zipped;
		//Two Bins
		//Calculate the overlapp
		Bin res;
		//Find the start of the Overlapp
		if (k1.left > k2.left)
			res.left = k1.left;
		else
			res.left = k2.left;

		//Find the end of the Overlapp
		if (k1.right < k2.right)
			res.right = k1.right;
		else
			res.right = k2.right;

		//Check that the Overlapp starts before it ends
		if (res.left <= res.right)
			//Figure out how much smaler the Overlapp is then K1
			sum *= (res.right - res.left) / (k1.right - k1.left);
		else
			return 0;
	}
	return sum;
}

//Get the Count of a Key that might not be in the DV explicitly
//by a weighted sum of all the Keys that are in the DV
//weighted by the overlapp of the given Key with the Keys of the DV
double DistributionalValue::get_contained_count(const NBin &h) const
{
	double res = 0;
	for (auto v : _value)
	{
		double weigth = DistributionalValue::conditional_probabilty(h,v.first);
		res += v.second * weigth;
	}
	return res;
}

double DistributionalValue::get_mode(const NBin &val) const
{
	return get_mode_for(get_count(val));
}
double DistributionalValue::get_mean(const NBin &val) const
{
	return get_mean_for(get_count(val));
}
double DistributionalValue::get_contained_mean(const NBin &val) const
{
	return get_mean_for(get_contained_count(val));
}
double DistributionalValue::get_var(const NBin &val) const
{
	return get_var_for(get_count(val));
}

std::string DistributionalValue::to_string(const std::string& indent) const
{
	std::stringstream ss;
	if (_value.size() == 0)
		ss << "Empty DistributionalValue" << std::endl;
	for (auto elem : _value)
	{
		ss << indent << "{";
		for (auto bin : elem.first)
		{
			ss << "["
			   << bin.left
			   << ","
			   << bin.right
			   << ");";
		}
		ss.seekp(-1,std::ios_base::end);
		ss << "}"
		   << " Count: "
		   << std::setprecision(18)
		   << elem.second
		   << " Mean: "
		   << get_mean_for(elem.second)
		   << std::endl;
	}
	return ss.str();
}

bool DistributionalValue::operator==(const Value& other) const
{
	if (DISTRIBUTIONAL_VALUE != other.get_type()) return false;
	const DistributionalValue* dov = (const DistributionalValue*) &other;

	if (_value.size() != dov->_value.size()) return false;

	for (auto elem : _value) {
		double v1 = elem.second;
		double v2 = dov->_value.get(elem.first);
		if (not is_approx_eq_ulp(v1,v2,24))
			return false;
	}
	return true;
}

std::string oc_to_string(const DistributionalValuePtr& dvp,
                         const std::string& indent)
{
	return dvp->to_string(indent);
}
