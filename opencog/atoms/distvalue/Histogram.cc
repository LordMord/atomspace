/*
 * opencog/truthvalue/DistributionalValue.cc
 *
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


#include <opencog/atoms/distvalue/Histogram.h>

using namespace opencog;

//Calculates the Middle of all intervals in a Key
template<typename CT>
DVec Histogram<CT>::center_of_bin(const NBin &k)
{
	auto size = k.size();
	DVec res(size);
	for (int i = 0; i < size; i++)
		res[i] = (k[i].left + k[i].right)/2;
	return res;
}

template<typename CT>
DVec Histogram<CT>::get_mean() const
{
	DVec res;
	for (auto elem : *this)
		res = res + center_of_bin(elem->first);
	return res / this->size();
}

template<typename CT>
bool Histogram<CT>::has_bin(const NBin &k) const
{
	auto it = this->find(k);
	return it != this->end();
}

template<typename CT>
NBinSeq Histogram<CT>::get_bins() const
{
	auto size = this->size();
	NBinSeq res(size);
	for (int i = 0; i < size; i++)
		res[i] = this->operator[](i).first;
	return res;
}
