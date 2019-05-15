/*
 * opencog/atoms/distvalue/ConditionalDV.h
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

#ifndef _OPENCOG_CONDITIONAL_GDTV_H
#define _OPENCOG_CONDITIONAL_GDTV_H

#include <memory>
#include <string>
#include <vector>
#include <limits>

#include <opencog/util/exceptions.h>
#include <opencog/atoms/value/Value.h>
#include <opencog/atoms/distvalue/DistributionalValue.h>

/** \addtogroup grp_atomspace
 *	@{
 */

namespace opencog
{

class DistributionalValue;
typedef std::shared_ptr<const DistributionalValue> DistributionalValuePtr;

class ConditionalDV;
typedef std::shared_ptr<const ConditionalDV> ConditionalDVPtr;

typedef CTHist<CTHist<double>> CDVrep;

/*
 * This class is used for conditional distributions. It's represented by a
 * histogram where each bin does not contain a count but a DistributionalValue
 */
class ConditionalDV
	: public Value
{
	CDVrep _value;

	// Disallow assignment -- truth values are immutable!
	ConditionalDV& operator=(const ConditionalDV& rhs) {
		throw RuntimeException(TRACE_INFO, "Cannot modify truth values!");
	}

public:
	ConditionalDV();
	ConditionalDV(const CDVrep&);
	ConditionalDV(const DVecSeq &conds,
                      const std::vector<DistributionalValuePtr> &dvs);

	const CDVrep& value() const { return _value; }

	static ConditionalDVPtr createCDV();
	static ConditionalDVPtr createCDV(const CDVrep&);
	static ConditionalDVPtr createCDV(const DVecSeq &conds,
                                          const std::vector<DistributionalValuePtr> &dvs);

	ConditionalDVPtr merge(ConditionalDVPtr) const;

	DVecSeq get_conditions() const;
	std::vector<DistributionalValuePtr> get_unconditionals() const;

	CTHist<double> get_unconditionalP(const DVec&) const;
	ConditionalDVPtr remap(const DVecSeq&) const;
	DistributionalValuePtr get_unconditional(DistributionalValuePtr) const;

	DistributionalValuePtr get_joint_probability(DistributionalValuePtr) const;

	//Consequent-Disjuction-Elimination
	//ConditionalDVPtr CDE(ConditionalDVPtr) const;

	double total_count() const;
	double avg_count() const;

	virtual bool operator==(const Value& rhs) const;

	using Value::to_string;
	std::string to_string(const std::string&) const;

};

static inline ConditionalDVPtr ConditionalDVCast(const ValuePtr& pa)
	{ return std::dynamic_pointer_cast<const ConditionalDV>(pa); }

static inline ValuePtr ValueCast(const ConditionalDVPtr& cdv)
{
	// This should have worked!?
	// return std::const_pointer_cast<Value>(tv);

	// This, too, should have worked!?
	// return std::shared_ptr<Value>(tv, const_cast<Value*>(tv.get()));

	// This works...
	return std::shared_ptr<Value>(cdv, (Value*) cdv.get());
}
} // namespace opencog

/** @}*/
#endif // _OPENCOG_TRUTH_VALUE_H
