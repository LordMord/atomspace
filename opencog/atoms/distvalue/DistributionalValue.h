/*
 * opencog/truthvalue/DistributionalValue.h
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

#ifndef _OPENCOG_DISTRIBUTIONAL_VALUE_H
#define _OPENCOG_DISTRIBUTIONAL_VALUE_H

#include <memory>
#include <string>
#include <vector>
#include <limits>

#include <opencog/atoms/value/Value.h>
#include <opencog/atoms/distvalue/CTHist.h>
#include <opencog/util/Counter.h>

/** \addtogroup grp_atomspace
 *	@{
 */

namespace opencog
{

typedef double strength_t;
typedef double confidence_t;
typedef double count_t;

class DistributionalValue;

typedef std::shared_ptr<const DistributionalValue> DistributionalValuePtr;

class AtomSpace;
class ConditionalDV;

/**
 * This class is used for unconditional distributions which are stored as
 * dirichlet distributions where each class is the bin of a multi dimensional
 * histogram
 */
class DistributionalValue
	: public Value
{
	CTHist<double> _value;
	friend ConditionalDV;

public:
	static count_t DEFAULT_K;

	DistributionalValue();
	DistributionalValue(const CTHist<double>&);
	DistributionalValue(double, double);

	const CTHist<double>& value() const { return _value; }

	static DistributionalValuePtr TRUE_TV();
	static DistributionalValuePtr FALSE_TV();
	static DistributionalValuePtr DEFAULT_TV();
	static DistributionalValuePtr createDV(double, double);
	static DistributionalValuePtr createDV(const CTHist<double>&);

	//Utility functions to convert between confidence and count
	//Should this be moved elsewhere???
	static double to_conf(int c);
	static int to_count(double);

	bool is_uniform() const;

	std::map<DVec,double> bin_modes() const;
	std::map<DVec,double> bin_means() const;
	std::map<DVec,double> bin_vars() const;

	double get_mode_for(double) const;
	double get_mean_for(double) const;
	double get_var_for(double) const;

	void add_evidence(const DVec&);
	DistributionalValuePtr merge(DistributionalValuePtr) const;
	DistributionalValuePtr remap(const DVecSeq&) const;
	DistributionalValuePtr mirrorLinf() const;

	double total_count() const;
	double get_confidence() const;

	double get_mean(const DVec&) const;
	double get_mode(const DVec&) const;
	double get_var(const DVec&) const;

	virtual bool operator==(const Value& rhs) const;

	using Value::to_string;
	std::string to_string(const std::string&) const;
};

static inline DistributionalValuePtr DistributionalValueCast(const ValuePtr& pa)
	{ return std::dynamic_pointer_cast<const DistributionalValue>(pa); }

static inline ValuePtr ValueCast(const DistributionalValuePtr& dv)
{
	// This should have worked!?
	// return std::const_pointer_cast<Value>(tv);

	// This, too, should have worked!?
	// return std::shared_ptr<Value>(tv, const_cast<Value*>(tv.get()));

	// This works...
	return std::shared_ptr<Value>(dv, (Value*) dv.get());
}

std::string oc_to_string(const DistributionalValue&, const std::string&);

} // namespace opencog

/** @}*/
#endif // _OPENCOG_DISTRIBUTIONAL_VALUE_H
