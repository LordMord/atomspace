/*
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

#include <opencog/atoms/distvalue/CHistIter.h>

#include <opencog/atoms/distvalue/CHist.h>

using namespace opencog;

template <typename c_typ, typename pointer, typename ref>
iterator_template<c_typ,pointer,ref>::iterator_template(uint i, uint d, ref & h)
		: idx(i) , dir(d) , hist(h) {}

template <typename c_typ, typename pointer, typename ref>
iterator_template<c_typ,pointer,ref>&
iterator_template<c_typ,pointer,ref>::operator++()
{
	idx = hist.next(idx,dir);
	return *this;
}

template <typename c_typ, typename pointer, typename ref>
iterator_template<c_typ,pointer,ref>
iterator_template<c_typ,pointer,ref>::operator++(int)
{
	auto result(*this);
	++(*this);
	return result;
}

template <typename c_typ, typename pointer, typename ref>
bool iterator_template<c_typ,pointer,ref>::operator==(iterator_template<c_typ,pointer,ref> other) const
{
	return (idx == other.idx && dir == other.dir && hist == other.hist);
}

template <typename c_typ, typename pointer, typename ref>
bool iterator_template<c_typ,pointer,ref>::operator!=(iterator_template<c_typ,pointer,ref> other) const
{
	return (idx != other.idx || dir != other.dir || hist != other.hist);
}

template <typename c_typ, typename pointer, typename ref>
const Node<c_typ>& iterator_template<c_typ,pointer, ref>::operator*()
{
	return hist.nodes[idx];
}

template <typename c_typ, typename pointer, typename ref>
pointer iterator_template<c_typ,pointer, ref>::operator->()
{
	return &(hist.nodes[idx]);
}
