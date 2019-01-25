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


#ifndef _OPENCOG_CHIST_H
#define _OPENCOG_CHIST_H

#include <vector>

namespace opencog
{

typedef std::pair<double,double> elem;

template <size_t len>
struct Node {
	int *pos;
	double count;
};

template <size_t len>
class CHist
{
	int size;
	int count;
	int subs;
	int dimensions;
	std::vector<pos_t> limits;
	std::vector<Node<pos_t>> nodes;
	int (*cmp)(const pos_t&,const pos_t&);
	double (*dist)(const pos_t&,const pos_t&);

	int parent(int);
	int child(int, int);
	int height(int);
	int child_opposite(int);

	void move_to(int, int);
	void shift_down(int, int);
	void rotate_up(int);
	void rebalance(int);
	void make_space(int,int);

	int min_max_heights(int, int*, int*);

	int max_height_child(int);

	void merge(int, pos_t, double);

	void insertFill(pos_t,double);
	void insertMerge(pos_t,double);

	void dumpP(int idx);
	void print(int, int);
public:
	CHist(int,int,int (*)(const pos_t&, const pos_t&),
	      double (*)(const pos_t&, const pos_t&));

	void insert(pos_t,double);

	static CHist<pos_t> merge(CHist<pos_t>, CHist<pos_t>);

	void dump();

	void print();
};

}

#endif // _OPENCOG_CHIST_H
