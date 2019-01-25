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
#include <algorithm>

using namespace opencog;

double* CHist::mul(double * ds, double d)
{
	for (uint i = 0; i < dimensions; i++)
	{
		*(ds+i) = *(ds+i) * d;
	}
	return ds;
}

double* CHist::add(double * ds1, double * ds2)
{
	for (uint i = 0; i < dimensions; i++)
		*(ds1+i) = *(ds1+i) + *(ds2+i);
	return ds1;
}

double* CHist::div(double * ds, double d)
{
	for (uint i = 0; i < dimensions; i++)
		*(ds+i) = *(ds+i) / d;
	return ds;
}

CHist::CHist(uint s,uint d)
	: size(s) , count(0) , subs(pow(2,d)) , dimensions(d)
{
	if (d >= 32)
		throw RuntimeException(TRACE_INFO,"More then 31 Dimensions not supported.");
	nodes.resize(s);
	limits.resize(subs);
}

CHist::~CHist()
{
	for (Node elem : nodes)
		delete[] elem.pos;
}

uint CHist::cmp(double *p1, double *p2)
{
	uint res = 0;
	for (uint i = 0; i < dimensions; i++)
	{
		if (*(p1+i) < *(p2+i))
			res = res | (1 << i);
	}
	return (res+1);
}

double CHist::dist(double *p1, double *p2)
{
	if (dimensions == 1)
		return fabs(*p1 - *p2);
	else
	{
		double res = 0;
		for (uint i = 0; i < dimensions; i++)
			res += pow(*(p1+i) - *(p2+i),2);
		return sqrt(res);
	}
}

uint CHist::parent(uint idx)
{
	if (idx == 0)
		throw RuntimeException(TRACE_INFO,"No parent for root node.");
	return (idx - 1) / subs;
}

uint CHist::child(uint idx,uint child)
{
	return idx * subs + child;
}


uint CHist::height(uint idx)
{
	if (idx >= size || nodes[idx].count == 0)
		return 0;
	uint max = 0;
	for (uint i = 1; i <= subs; i ++)
	{
		uint tmp = height(child(idx,i));
		if (tmp > max)
			max = tmp;
	}
	return max+1;
}

//Move a sub-tree starting at idx to towards
//Don't use if towards is in the subtree startin at idx
//use shift down instead
//FIXME: Add check against this ^
void CHist::move_to(uint idx,uint towards)
{
	if (nodes[idx].count == 0 || idx >= size)
		return;

	nodes[towards] = nodes[idx];
	nodes[idx].pos = nullptr;
	nodes[idx].count = 0;
	for (uint i = 1; i <= subs; i ++)
	{
		move_to(child(idx,i),child(towards,i));
	}
}

//Move a subtree starting at idx dowards in direction dir
//idx: Index of element to shift_down
//dir: Direction of shift, not an index
void CHist::shift_down(uint idx,uint dir)
{
	if (nodes[idx].count == 0 || idx >= size)
		return;

	//First Shift down the child in the specified direction
	//As otherwise nodes could get moved uinto this subtree from other children
	shift_down(child(idx,dir),dir);

	uint towards = child(idx,dir);
	for (uint i = 1; i <= subs; i ++)
	{
		if (i == dir) //We already did this.
			continue;
		move_to(child(idx,i),child(towards,i));
	}
	nodes[towards] = nodes[idx];
	nodes[idx].pos = nullptr;
	nodes[idx].count = 0;
}

uint CHist::child_opposite(uint child_idx)
{
	return child_idx + subs - 1 - 2 * (child_idx - 1);
}

void CHist::rotate_up(uint idx)
{
	if (nodes[idx].count == 0 || idx >= size)
		return;

	uint p = parent(idx);

	uint child_idx = idx - p * subs;
	uint child_opp = child_opposite(child_idx);

	//Move down the parrents Children except idx
	for (uint i = 1; i <= subs; i ++)
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
	nodes[child(idx,child_opp)].pos = nullptr;
	nodes[child(idx,child_opp)].count = 0;

	//Move idx Upwards
	nodes[p] = nodes[idx];
	//Incase there is nothing to shift up into this position
	nodes[idx].pos = nullptr;
	nodes[idx].count = 0;

	//Move Idx's children Upwards except child_opp
	for (uint i = 1; i <= subs; i ++)
	{
		if (i == child_opp)
			continue;
		move_to(child(idx,i),child(p,i));
	}

}

uint CHist::min_max_heights(uint idx, uint *min, uint *max)
{
	*min = size;
	*max = 0;
	uint max_idx = -1;
	for (uint i = 1; i <= subs; i ++)
	{
		uint tmp = height(child(idx,i));
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

uint CHist::max_height_child(uint idx)
{
	uint max = 0;
	uint max_idx = -1;
	for (uint i = 1; i <= subs; i ++)
	{
		uint tmp = height(child(idx,i));
		if (tmp > max)
		{
			max = tmp;
			max_idx = child(idx,i);
		}
	}
	return max_idx;
}

void CHist::rebalance(uint idx)
{
	while (true)
	{
		uint min;
		uint max;
		uint max_idx = min_max_heights(idx,&min,&max);
		uint child_idx = max_idx - idx * subs;
		uint child_opp = child_opposite(child_idx);

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
void CHist::make_space(uint idx, uint dir)
{
	while (true)
	{
		uint min;
		uint max;
		uint max_idx = min_max_heights(idx,&min,&max);

		if (1 <= (max-min))
		{
			rotate_up(max_idx);
			return;
		}

		//End if we reached the root
		if (0 == idx)
			return;

		//Go up 1 layer only if it is in the specified direction
		uint p = parent(idx);
		if (idx != child(p,dir))
			return;
		idx = p;
	}
}

void CHist::mergeNode(uint idx, double * pos,double c)
{
	Node &n = nodes[idx];
	n.pos = div(add(mul(n.pos,n.count),mul(pos,c)),(n.count + c));
	n.count += c;

	if (std::none_of(limits.begin(),limits.end(),
	                 [pos](double *l){ return pos == l;}))
		delete[] pos;
}

void CHist::insert(std::vector<double> posv,double c)
{
	if (posv.size() != dimensions)
		throw RuntimeException(TRACE_INFO,"Position Vector needs to be the same lenght as the number of dimensions!");

	auto size = sizeof(double)*dimensions;
	double * pos = (double*)malloc(size);
	memcpy(pos,&posv[0],size);

	insertP(pos,c);
}

void CHist::insertP(double * pos,double c)
{
	if (count == 0)
		std::fill(limits.begin(),limits.end(),pos);
	else
		for (uint i = 1; i <= subs; i++)
			if (i == cmp(limits[i-1],pos))
				limits[i-1] = pos;

	//std::cout << "s,c: " << size << "," << count << "," << to_string(pos) << "\n";
	if (size == count)
		insertMerge(pos,c);
	else
		insertFill(pos,c);
}

void CHist::insertFill(double * pos,double c)
{
	uint i;
	double mindist = std::numeric_limits<double>::infinity();
	uint minidx = -1;
	for (i = 0; i < size; )
	{
		if (nodes[i].count == 0)
		{
			count++;
			nodes[i].pos = pos;
			nodes[i].count = c;

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
			//std::cout << "DistF: " << tmp << std::endl;
			if (tmp < mindist)
			{
				mindist = tmp;
				minidx = i;
			}
			uint child_idx = cmp(nodes[i].pos,pos);
			i = child(i,child_idx);
		}
	}

	uint idx = parent(i);
	if (cmp(nodes[idx].pos,pos) == cmp(nodes[parent(idx)].pos,nodes[idx].pos))
	{
		uint p = parent(idx);
		uint dir = idx - p * subs;
		//std::cout << "Before make_space\n";
		//print();
		make_space(p,dir);
		//std::cout << "After make_space\n";
		//print();
		if (nodes[idx].count == 0)
		{
			count++;
			nodes[idx].pos = pos;
			nodes[idx].count = c;
			return;
		}
	}

	mergeNode(minidx,pos,c);
}

void CHist::insertMerge(double * pos,double c)
{
	uint i;
	double mindist = std::numeric_limits<double>::infinity();
	uint minidx = -1;
	for (i = 0; i < size; )
	{
		if (nodes[i].pos == pos)
		{
			mergeNode(i,pos,c);
			return;
		}

		double tmp = dist(nodes[i].pos,pos);

		//std::cout << "Dist: " << tmp << std::endl;

		if (tmp < mindist)
		{
			mindist = tmp;
			minidx = i;
		}

		uint child_idx = cmp(nodes[i].pos,pos);
		uint child_i = child(i,child_idx);
		if (nodes[child_i].count == 0)
		{
			mergeNode(minidx,pos,c);
			return;
		}
		else
		{
			i = child_i;
		}
	}
	mergeNode(minidx,pos,c);
}

void CHist::dump()
{
	std::cout << subs << std::endl;
	for (uint i = 0; i < subs; i++)
		std::cout << to_string(limits[i]) << "\n";
	dumpP(0);
}

void CHist::dumpP(uint idx)
{
	uint child_l = child(idx,1);
	if (child_l < size && nodes[child_l].count != 0)
		dumpP(child_l);

	std::cout << to_string(nodes[idx].pos)
			  << "," << nodes[idx].count << std::endl;

	uint child_r = child(idx,2);
	if (child_r < size && nodes[child_r].count != 0)
		dumpP(child_r);
}

void CHist::print()
{
	std::cout << "limits: ";
	for (uint i = 0; i < subs; i++)
		std::cout << "i: " << (i+1) << "," << to_string(limits[i]) << " ";
	std::cout << std::endl;
	print(0,0);
	std::cout << std::endl;
}

void CHist::print(uint idx, uint d)
{
	uint i;

	if (size <= idx || nodes[idx].count == 0)
		return;

	for (i=0; i<d; i++)
		printf("  ");
	printf("%i: ", idx);

	std::cout << "pos: " << to_string(nodes[idx].pos)
			  << " count: " << nodes[idx].count << std::endl;
	for (uint i = 1; i <= subs; i ++)
	{
		print(child(idx,i),d+1);
	}
}


CHist CHist::merge(CHist h1, CHist h2)
{
	if (h1.size != h2.size)
		throw RuntimeException(TRACE_INFO,"CHists must be same size");
	if (h1.subs != h2.subs)
		throw RuntimeException(TRACE_INFO,"CHists must be same dimensions");

	CHist res = CHist(h1.size,h1.dimensions);

	for (uint i = 0; i < h1.size; i++)
	{
		auto size = sizeof(double)*res.dimensions;

		Node n1 = h1.nodes[i];
		double *n1pos = (double*)malloc(size);
		memcpy(n1pos,n1.pos,size);
		res.insertP(n1pos,n1.count);

		Node n2 = h2.nodes[i];
		double *n2pos = (double*)malloc(size);
		memcpy(n2pos,n2.pos,size);
		res.insertP(n2pos,n2.count);
	}


	for (uint i = 0; i < res.subs; i++)
		if ((i+1) == res.cmp(h1.limits[i],h2.limits[i]))
			res.limits[i] = h2.limits[i];
		else
			res.limits[i] = h1.limits[i];

	return res;
}

std::string CHist::to_string(double * pos)
{
	std::stringstream ss;
	if (dimensions == 1)
	{
		ss << *pos;
	}
	else
	{
		ss << "(";
		for (uint i = 0; i < dimensions; i++)
			ss << *(pos+i) << ",";
		ss.seekp(-1,std::ios_base::end);
		ss << ")";
	}
	return ss.str();
}
