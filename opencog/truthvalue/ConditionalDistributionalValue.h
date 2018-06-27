/*
 * opencog/truthvalue/ConditionalDistributionalValue.h
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

#ifndef _OPENCOG_CONDITIONAL_DISTRIBUTIONAL_VALUE_H
#define _OPENCOG_CONDITIONAL_DISTRIBUTIONAL_VALUE_H

#include <memory>
#include <string>
#include <vector>
#include <limits>

#include <opencog/util/exceptions.h>
#include <opencog/atoms/base/ProtoAtom.h>

#include <opencog/truthvalue/DistributionalValue.h>
/** \addtogroup grp_atomspace
 *  @{
 */

namespace opencog
{

class ConditionalDistributionalValue;
typedef std::shared_ptr<const ConditionalDistributionalValue> ConditionalDistributionalValuePtr;

typedef std::map<Handle,HandleCounter> DistributionalValuerep;

class ConditionalDistributionalValue
    : public ProtoAtom
{
    DistributionalValuerep value;

    // Disallow assignment -- truth values are immutable!
    ConditionalDistributionalValue& operator=(const ConditionalDistributionalValue& rhs) {
        throw RuntimeException(TRACE_INFO, "Cannot modify truth values!");
    }

public:
    ConditionalDistributionalValue();
    ConditionalDistributionalValue(DistributionalValuerep);

    DistributionalValuePtr getUnconditional(Handle);
    DistributionalValuePtr getUnconditional(DistributionalValuePtr);

};

} // namespace opencog

/** @}*/
#endif // _OPENCOG_TRUTH_VALUE_H
