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


#include <opencog/atoms/distvalue/CHist.h>
#include <opencog/util/exceptions.h>
#include <stdio.h>
#include <cmath>
#include <sstream>
#include <limits>

using namespace opencog;

template class CHist<double>;
template class CHist<elem>;

elem operator*(const elem& e,const double& d)
{
	return elem(e.first * d,e.second * d);
}
elem operator/(const elem& e,const double& d)
{
	return elem(e.first / d,e.second / d);
}
elem operator+(const elem& e,const double& d)
{
	return elem(e.first + d,e.second + d);
}
elem operator+(const elem& e1,const elem& e2)
{
	return elem(e1.first + e2.first,e1.second + e2.second);
}

std::string to_string(const double& d)
{
	std::stringstream ss;
	ss << d;
	return ss.str();
}
std::string to_string(const elem& e)
{
	std::stringstream ss;
	ss << "(" << e.first << "," << e.second << ")";
	return ss.str();
}

template <typename pos_t>
CHist<pos_t>::CHist(int s,int d,int (*cmpf)(const pos_t&, const pos_t&),
					double (*distf)(const pos_t&, const pos_t&))
	: size(s) , count(0) , subs(pow(2,d)) , dimensions(d)
{
	nodes.resize(s);
	limits.resize(subs);
	cmp = cmpf;
	dist = distf;
}

template <typename pos_t>
int CHist<pos_t>::parent(int idx)
{
	if (idx == 0)
		throw RuntimeException(TRACE_INFO,"No parent for root node.");
	return (idx - 1) / subs;
}

template <typename pos_t>
int CHist<pos_t>::child(int idx,int child)
{
	return idx * subs + child;
}


template <typename pos_t>
int CHist<pos_t>::height(int idx)
{
	if (idx >= size || nodes[idx].count == 0)
		return 0;
	int max = 0;
	for (int i = 1; i <= subs; i ++)
	{
		int tmp = height(child(idx,i));
		if (tmp > max)
			max = tmp;
	}
	return max+1;
}

//Move a sub-tree starting at idx to towards
//Don't use if towards is in the subtree startin at idx
//use shift down instead
//FIXME: Add check against this ^
template <typename pos_t>
void CHist<pos_t>::move_to(int idx,int towards)
{
	if (nodes[idx].count == 0 || idx >= size)
		return;

	nodes[towards] = nodes[idx];
	nodes[idx].count = 0;
	for (int i = 1; i <= subs; i ++)
	{
		move_to(child(idx,i),child(towards,i));
	}
}

//Move a subtree starting at idx dowards in direction dir
//idx: Index of element to shift_down
//dir: Direction of shift, not an index
template <typename pos_t>
void CHist<pos_t>::shift_down(int idx,int dir)
{
	if (nodes[idx].count == 0 || idx >= size)
		return;

	//First Shift down the child in the specified direction
	//As otherwise nodes could get moved into this subtree from other children
	shift_down(child(idx,dir),dir);

	int towards = child(idx,dir);
	for (int i = 1; i <= subs; i ++)
	{
		if (i == dir) //We already did this.
			continue;
		move_to(child(idx,i),child(towards,i));
	}
	nodes[towards] = nodes[idx];
	nodes[idx].count = 0;
}

template <typename pos_t>
int CHist<pos_t>::child_opposite(int child_idx)
{
	return child_idx + subs - 1 - 2 * (child_idx - 1);
}

template <typename pos_t>
void CHist<pos_t>::rotate_up(int idx)
{
	if (nodes[idx].count == 0 || idx >= size)
		return;

	int p = parent(idx);

	int child_idx = idx - p * subs;
	int child_opp = child_opposite(child_idx);

	//Move down the parrents Children except idx
	for (int i = 1; i <= subs; i ++)
	{
		if (i == child_idx)
			continue;
		shift_down(child(p,i),child_opp);
	}
	//Move down the parrent
	nodes[child(p,child_opp)] = nodes[p];

	//Move child_opp to other side of tree
	move_to(child(idx,child_opp),child(child(p,child_opp),child_idx));
	//Set it to 0 in original location
	nodes[child(idx,child_opp)].count = 0;

	//Move idx Upwards
	nodes[p] = nodes[idx];
	//Incase there is nothing to shift up into this position
	nodes[idx].count = 0;

	//Move Idx's children Upwards except child_opp
	for (int i = 1; i <= subs; i ++)
	{
		if (i == child_opp)
			continue;
		move_to(child(idx,i),child(p,i));
	}

}

template <typename pos_t>
int CHist<pos_t>::min_max_heights(int idx, int *min, int *max)
{
	*min = size;
	*max = 0;
	int max_idx = -1;
	for (int i = 1; i <= subs; i ++)
	{
		int tmp = height(child(idx,i));
		if (tmp > *max)
		{
			*max = tmp;
			max_idx = child(idx,i);
		}
		if (tmp < *min)
			*min = tmp;
	}
	return max_idx;
}

template <typename pos_t>
int CHist<pos_t>::max_height_child(int idx)
{
	int max = 0;
	int max_idx = -1;
	for (int i = 1; i <= subs; i ++)
	{
		int tmp = height(child(idx,i));
		if (tmp > max)
		{
			max = tmp;
			max_idx = child(idx,i);
		}
	}
	return max_idx;
}

template <typename pos_t>
void CHist<pos_t>::rebalance(int idx)
{
	while (true)
	{
		int min;
		int max;
		int max_idx = min_max_heights(idx,&min,&max);
		int child_idx = max_idx - idx * subs;
		int child_opp = child_opposite(child_idx);

		//std::cout << " idx: " << idx << " min: " << min
	   // 	      << " max: " << max << " maxidx: " << max_idx << std::endl;

		if (2 <= (max-min))
		{
			if (child(max_idx,child_opp) == max_height_child(max_idx))
				rotate_up(child(max_idx,child_opp));
			rotate_up(max_idx);
		}

		//End if we reached the root
		if (0 == idx)
			return;

		//Check if we need to reblance the parent next
		idx = parent(idx);
	}
}

//Used to rotate the tree to make space for a new entry
//at a specificy location
//We can only rotate in a specific direction to not mess up the ordering
template <typename pos_t>
void CHist<pos_t>::make_space(int idx, int dir)
{
	while (true)
	{
		int min;
		int max;
		int max_idx = min_max_heights(idx,&min,&max);

		if (1 <= (max-min))
		{
			rotate_up(max_idx);
			return;
		}

		//End if we reached the root
		if (0 == idx)
			return;

		//Go up 1 layer only if it is in the specified direction
		int p = parent(idx);
		if (idx != child(p,dir))
			return;
		idx = p;
	}
}

template <typename pos_t>
void CHist<pos_t>::merge(int idx, pos_t pos,double c)
{
	Node<pos_t> &n = nodes[idx];
	n.pos = (n.pos * n.count + pos * c) / (n.count + c);
	n.count += c;
}

template <typename pos_t>
void CHist<pos_t>::insert(pos_t pos,double c)
{
	if (count == 0)
		std::fill(limits.begin(),limits.end(),pos);
	else
		for (int i = 1; i <= subs; i++)
			if (i == cmp(limits[i-1],pos))
				limits[i-1] = pos;

	if (size == count)
		insertMerge(pos,c);
	else
		insertFill(pos,c);
}

template <typename pos_t>
void CHist<pos_t>::insertFill(pos_t pos,double c)
{
	int i;
	double mindist = std::numeric_limits<double>::infinity();
	int minidx = -1;
	for (i = 0; i < size; )
	{
		if (nodes[i].count == 0)
		{
			nodes[i].pos = pos;
			nodes[i].count = c;
			count++;

			if (i == 0)
				return;

			rebalance(parent(i));
			return;
		}

		if (nodes[i].pos == pos)
		{
			nodes[i].count += c;
			return;
		}
		else
		{
			double tmp = dist(nodes[i].pos,pos);
			if (tmp < mindist)
			{
				mindist = tmp;
				minidx = i;
			}
			int child_idx = cmp(nodes[i].pos,pos);
			i = child(i,child_idx);
		}
	}

	int idx = parent(i);
	if (cmp(nodes[idx].pos,pos) == cmp(nodes[parent(idx)].pos,nodes[idx].pos))
	{
		int p = parent(idx);
		int dir = idx - p * subs;
		make_space(p,dir);
		if (nodes[idx].count == 0)
		{
			count++;
			merge(idx,pos,c);
			return;
		}
	}

	merge(minidx,pos,c);
}

template <typename pos_t>
void CHist<pos_t>::insertMerge(pos_t pos,double c)
{
	int i;
	double mindist = std::numeric_limits<double>::infinity();
	int minidx = -1;
	for (i = 0; i < size; )
	{
		if (nodes[i].pos == pos)
		{
			merge(i,pos,c);
			return;
		}

		double tmp = dist(nodes[i].pos,pos);

		if (tmp < mindist)
		{
			mindist = tmp;
			minidx = i;
		}

		int child_idx = cmp(nodes[i].pos,pos);
		int child_i = child(i,child_idx);
		if (nodes[child_i].count == 0)
		{
			merge(minidx,pos,c);
			return;
		}
		else
		{
			i = child_i;
		}
	}
	merge(minidx,pos,c);
}

template <typename pos_t>
void CHist<pos_t>::dump()
{
	std::cout << subs << std::endl;
	for (int i = 0; i < subs; i++)
		std::cout << to_string(limits[i]) << "\n";
	dumpP(0);
}

template <typename pos_t>
void CHist<pos_t>::dumpP(int idx)
{
	int child_l = child(idx,1);
	if (child_l < size && nodes[child_l].count != 0)
		dumpP(child_l);

	std::cout << to_string(nodes[idx].pos)
			  << "," << nodes[idx].count << std::endl;

	int child_r = child(idx,2);
	if (child_r < size && nodes[child_r].count != 0)
		dumpP(child_r);
}

template <typename pos_t>
void CHist<pos_t>::print()
{
	std::cout << "limits: ";
	for (int i = 0; i < subs; i++)
		std::cout << "i: " << (i+1) << "," << to_string(limits[i]) << " ";
	std::cout << std::endl;
	print(0,0);
}

template <typename pos_t>
void CHist<pos_t>::print(int idx, int d)
{
	int i;

	if (size <= idx)
		return;

	for (i=0; i<d; i++)
		printf("  ");
	printf("%i: ", idx);

	//if (nodes[idx].count == 0)
	//{
	//	printf("\n");
	//	return;
	//}

	std::cout << "pos: " << to_string(nodes[idx].pos)
			  << " count: " << nodes[idx].count << std::endl;
	for (int i = 1; i <= subs; i ++)
	{
		print(child(idx,i),d+1);
	}
}


template <typename pos_t>
CHist<pos_t> CHist<pos_t>::merge(CHist<pos_t> h1, CHist<pos_t> h2)
{
	if (h1.size != h2.size)
		throw RuntimeException(TRACE_INFO,"CHists must be same size");
	if (h1.subs != h2.subs)
		throw RuntimeException(TRACE_INFO,"CHists must be same dimensions");

	CHist<pos_t> res = CHist(h1.size,h1.dimensions,h1.cmp,h1.dist);

	for (int i = 0; i < h1.size; i++)
	{
		Node<pos_t> n1 = h1.nodes[i];
		res.insert(n1.pos,n1.count);
		Node<pos_t> n2 = h2.nodes[i];
		res.insert(n2.pos,n2.count);
	}


	for (int i = 0; i < res.subs; i++)
		if ((i+1) == res.cmp(h1.limits[i],h2.limits[i]))
			res.limits[i] = h2.limits[i];
		else
			res.limits[i] = h1.limits[i];

	return res;
}

