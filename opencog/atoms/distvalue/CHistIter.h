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

#ifndef _OPENCOG_CHISTITER_H
#define _OPENCOG_CHISTITER_H

#include <iterator>

namespace opencog
{

template <typename c_typ>
class CHist;

template <typename c_typ>
struct Node;

template <typename c_typ,typename pointer,typename ref>
class iterator_template :
	public std::iterator<std::forward_iterator_tag, // iterator_category
						 Node<c_typ> ,                     // value_type
						 std::ptrdiff_t,            // difference_type
						 pointer,                   // pointer
						 Node<c_typ>&>                     // reference

{
	uint idx;
	uint dir;
	ref & hist;
public:
	iterator_template<c_typ,pointer,ref>(uint i, uint d, ref & h);

	iterator_template<c_typ,pointer,ref>& operator++();

	iterator_template<c_typ,pointer,ref> operator++(int);

	bool operator==(iterator_template<c_typ,pointer,ref> other) const;

	bool operator!=(iterator_template<c_typ,pointer,ref> other) const;

	const Node<c_typ>& operator*();

	pointer operator->();
};

template<typename c_typ>
using iterator = iterator_template<c_typ, Node<c_typ> *, CHist<c_typ>>;

template<typename c_typ>
using const_iterator = iterator_template<c_typ, const Node<c_typ> *
											  , const CHist<c_typ>>;

template class iterator_template<double,Node<double> *, CHist<double>>;

template class iterator_template<double,const Node<double> *, const CHist<double>>;

template class iterator_template<CHist<double>
								 , Node<CHist<double>> *
								 , CHist<CHist<double>>>;

template class iterator_template<CHist<double>,
								 const Node<CHist<double>> *,
								 const CHist<CHist<double>>>;

} // namespace opencog

#endif // _OPENCOG_CHISTITER_H
