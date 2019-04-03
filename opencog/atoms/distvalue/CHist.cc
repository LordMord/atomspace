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
#include <opencog/util/numeric.h>
#include <stdio.h>
#include <cmath>
#include <sstream>
#include <limits>
#include <algorithm>

using namespace opencog;

template <typename c_typ>
double* CHist<c_typ>::mul(double * ds, double d) const
{
	for (uint i = 0; i < _dimensions; i++)
	{
		*(ds+i) = *(ds+i) * d;
	}
	return ds;
}

template <typename c_typ>
double* CHist<c_typ>::add(double * ds1, double * ds2) const
{
	for (uint i = 0; i < _dimensions; i++)
		*(ds1+i) = *(ds1+i) + *(ds2+i);
	return ds1;
}

template <typename c_typ>
double* CHist<c_typ>::div(double * ds, double d) const
{
	for (uint i = 0; i < _dimensions; i++)
		*(ds+i) = *(ds+i) / d;
	return ds;
}

template <typename c_typ>
bool CHist<c_typ>::eq(double * ds1, double * ds2) const
{
	for (uint i = 0; i < _dimensions; i++)
		if (!is_approx_eq_ulp(*(ds1+i),*(ds2+i),24))
			return false;
	return true;
}

template <typename c_typ>
bool CHist<c_typ>::eq(const Node<c_typ> &n1,const Node<c_typ> &n2) const
{
	return (eq(n1.pos,n2.pos) && n1.count == n2.count);
}

template <typename c_typ>
CHist<c_typ>::CHist(uint s,uint d)
	: _size(s) , _count_elems(0) , _total_count(0)
	, _subs(pow(2,d)) , _dimensions(d)
{
	if (d >= 32)
		throw RuntimeException(TRACE_INFO,"More then 31 Dimensions not supported.");
	nodes.resize(s,Node<c_typ>());
	limits.resize(_subs);

	auto size = sizeof(double)*_dimensions;
	for (uint i = 0; i < _subs; i++)
	{
		limits[i] = (double*)malloc(size);
		if(limits[i] == nullptr)
			throw RuntimeException(TRACE_INFO,"Malloc Failed");
	}
}

template <typename c_typ>
CHist<c_typ>::~CHist()
{
	for (Node<c_typ> elem : nodes)
		delete[] elem.pos;
	for (auto elem : limits)
		delete[] elem;
}

template <typename c_typ>
uint CHist<c_typ>::cmp(double *p1, double *p2) const
{
	uint res = 0;
	for (uint i = 0; i < _dimensions; i++)
	{
		if (*(p1+i) <= *(p2+i))
			res = res | (1 << i);
	}
	//printf("Cmp p1,p2,res: %f,%f,%i \n",*p1,*p2,res);
	return (res+1);
}

template <typename c_typ>
double CHist<c_typ>::dist(double *p1, double *p2) const
{
	if (_dimensions == 1)
		return fabs(*p1 - *p2);
	else
	{
		double res = 0;
		for (uint i = 0; i < _dimensions; i++)
			res += pow(*(p1+i) - *(p2+i),2);
		return sqrt(res);
	}
}

template <typename c_typ>
uint CHist<c_typ>::parent(uint idx) const
{
	if (idx == 0)
		throw RuntimeException(TRACE_INFO,"No parent for root node.");
	return (idx - 1) / _subs;
}

template <typename c_typ>
uint CHist<c_typ>::child(uint idx,uint child) const
{
	return idx * _subs + child;
}


template <typename c_typ>
uint CHist<c_typ>::height(uint idx) const
{
	if (idx >= _size || nodes[idx].count == 0)
		return 0;
	uint max = 0;
	for (uint i = 1; i <= _subs; i ++)
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
template <typename c_typ>
void CHist<c_typ>::move_to(uint idx,uint towards)
{
	if (nodes[idx].pos == nullptr || idx >= _size)
		return;

	nodes[towards] = nodes[idx];
	nodes[idx].pos = nullptr;
	nodes[idx].count = 0;
	for (uint i = 1; i <= _subs; i ++)
	{
		move_to(child(idx,i),child(towards,i));
	}
}

//Move a subtree starting at idx dowards in direction dir
//idx: Index of element to shift_down
//dir: Direction of shift, not an index
template <typename c_typ>
void CHist<c_typ>::shift_down(uint idx,uint dir)
{
	if (nodes[idx].pos == nullptr || idx >= _size)
		return;

	//First Shift down the child in the specified direction
	//As otherwise nodes could get moved uinto this subtree from other children
	shift_down(child(idx,dir),dir);

	uint towards = child(idx,dir);
	for (uint i = 1; i <= _subs; i ++)
	{
		if (i == dir) //We already did this.
			continue;
		move_to(child(idx,i),child(towards,i));
	}
	nodes[towards] = nodes[idx];
	nodes[idx].pos = nullptr;
	nodes[idx].count = 0;
}

template <typename c_typ>
uint CHist<c_typ>::child_opposite(uint child_idx) const
{
	return child_idx + _subs - 1 - 2 * (child_idx - 1);
}

template <typename c_typ>
void CHist<c_typ>::rotate_up(uint idx)
{
	if (nodes[idx].pos == nullptr || idx >= _size)
		return;

	uint p = parent(idx);

	uint child_idx = idx - p * _subs;
	uint child_opp = child_opposite(child_idx);

	//Move down the parrents Children except idx
	for (uint i = 1; i <= _subs; i ++)
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
	for (uint i = 1; i <= _subs; i ++)
	{
		if (i == child_opp)
			continue;
		move_to(child(idx,i),child(p,i));
	}

}

template <typename c_typ>
uint CHist<c_typ>::min_max_heights(uint idx, uint *min, uint *max) const
{
	*min = _size;
	*max = 0;
	uint max_idx = -1;
	for (uint i = 1; i <= _subs; i ++)
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

template <typename c_typ>
uint CHist<c_typ>::max_height_child(uint idx) const
{
	uint max = 0;
	uint max_idx = -1;
	for (uint i = 1; i <= _subs; i ++)
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

template <typename c_typ>
uint CHist<c_typ>::get_dir(uint p,uint c) const
{
	return c - p * _subs;
}

template <typename c_typ>
void CHist<c_typ>::rebalance(uint idx)
{
	while (true)
	{
		uint min;
		uint max;
		uint max_idx = min_max_heights(idx,&min,&max);
		uint child_idx = get_dir(idx,max_idx);
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
template <typename c_typ>
void CHist<c_typ>::make_space(uint idx, uint dir)
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

template <typename c_typ>
void CHist<c_typ>::mergeNode(uint idx, double * pos, c_typ c)
{
	Node<c_typ> &n = nodes[idx];
	double nc = get_count(n);
	double ac = get_count(c);
	n.pos = div(add(mul(n.pos,nc),mul(pos,ac)),(nc + ac));
	update_count(n,nc+ac);

	delete[] pos;
}


template <typename c_typ>
double * CHist<c_typ>::vecToArray(DVec vec) const
{
	if (vec.size() != _dimensions)
		throw RuntimeException(TRACE_INFO,"Vector needs to be the same lenght as the number of dimensions!");

	auto size = sizeof(double)*_dimensions;
	double * arr = (double*)malloc(size);
	memcpy(arr,&vec[0],size);

	return arr;
}

template <typename c_typ>
void CHist<c_typ>::insert(DVec posv,c_typ c)
{
	return insert(vecToArray(posv),c);
}

template <typename c_typ>
void CHist<c_typ>::insert(double * pos,c_typ c)
{
	//Update limits if necesary
	auto size = sizeof(double)*_dimensions;
	if (_count_elems == 0)
		for (uint i = 1; i <= _subs; i++)
			memcpy(limits[i-1],pos,size);
	else
		for (uint i = 1; i <= _subs; i++)
			if (i == cmp(limits[i-1],pos))
				memcpy(limits[i-1],pos,size);

	//std::cout << "s,c: " << _size << "," << _count_elems << "," << to_string(pos) << "\n";
	if (_size == _count_elems)
		insertMerge(pos,c);
	else
		insertFill(pos,c);
	_total_count++;
}

template <typename c_typ>
void CHist<c_typ>::insertFill(double * pos,c_typ c)
{
	uint i;
	double mindist = std::numeric_limits<double>::infinity();
	uint minidx = -1;
	for (i = 0; i < _size; )
	{
		if (nodes[i].count == 0)
		{
			_count_elems++;
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
		uint dir = get_dir(p,idx);
		//std::cout << "Before make_space\n";
		//print();
		make_space(p,dir);
		//std::cout << "After make_space\n";
		//print();
		if (nodes[idx].count == 0)
		{
			_count_elems++;
			nodes[idx].pos = pos;
			nodes[idx].count = c;
			return;
		}
	}

	mergeNode(minidx,pos,c);
}

template <typename c_typ>
void CHist<c_typ>::insertMerge(double * pos,c_typ c)
{
	uint i;
	double mindist = std::numeric_limits<double>::infinity();
	uint minidx = -1;
	for (i = 0; i < _size; )
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


//Actually mirroring about the L-infinity mean
//template <typename c_typ>
//CHist<c_typ> CHist<c_typ>::negate() const
//{
//	DVec minmax = minmax_count();
//	double total = minmax[0] + minmax[1];
//	CHist res = copy();
//
//	for (uint i = 0; i < _size; i ++)
//	{
//		res.nodes[i].count = total - res.nodes[i].count;
//	}
//	return res;
//}
//
////Get the lowest and highest count of all Interals
//template <typename c_typ>
//DVec CHist<c_typ>::minmax_count() const
//{
//	DVec minmax = DVec{std::numeric_limits<double>::max(),0};
//	for (auto elem : nodes)
//	{
//		if (minmax[0] >= get_count(elem))
//			minmax[0] = get_count(elem);
//		if (minmax[1] <= get_count(elem))
//			minmax[1] = get_count(elem);
//	}
//	return minmax;
//}


template <typename c_typ>
uint CHist<c_typ>::child_loop(uint idx,uint dir) const
{
	while (true)
	{
		uint child_idx = child(idx,dir);
		if (child_idx >= _size || nodes[child_idx].pos == nullptr)
			return idx;
		idx = child_idx;
	}
}

template <typename c_typ>
uint CHist<c_typ>::next(uint idx,uint &dir) const
{
	while (dir <= _subs)
	{
		uint child_idx = child(idx,dir);
		if (child_idx < _size && nodes[child_idx].pos != nullptr)
		{
			dir = 2;
			return child_loop(child_idx,1);
		}
		dir++;
	}

	uint p = parent(idx);
	if (p == 0 && child(p,_subs) == idx)
	{
		dir = 0;
		return 0;
	}
	while (child(p,_subs) == idx)
	{
		idx = p;
		p = parent(idx);
		if (p == 0 && child(p,_subs) == idx)
		{
			dir = 0;
			return 0;
		}
	}

	dir = 1 + get_dir(p,idx);

	if (dir-1 == _subs/2)
		return p;
	else
		return next(p,dir);
}

template <typename c_typ>
c_typ CHist<c_typ>::get(DVec pos) const
{
	return get(vecToArray(pos));
}

template <typename c_typ>
c_typ CHist<c_typ>::get(double * pos) const
{
	uint idx = 0;
	while (idx < _size || nodes[idx].pos == nullptr)
	{
		if (nodes[idx].pos == pos)
			return nodes[idx].count;
		uint dir = cmp(nodes[idx].pos,pos);
		idx = child(idx,dir);
	}
	throw RuntimeException(TRACE_INFO,"This Position is not an element of this Histogram.");
}

template <typename c_typ>
c_typ CHist<c_typ>::get_avg(DVec pos) const
{
	return get_avg(vecToArray(pos));
}

template <typename c_typ>
c_typ CHist<c_typ>::get_avg(double * pos) const
{
	double mindist = std::numeric_limits<double>::infinity();
	uint minidx = -1;

	uint idx = 0;
	while (idx < _size && nodes[idx].pos != nullptr)
	{
		if (eq(nodes[idx].pos,pos))
			return nodes[idx].count;

		double tmp = dist(nodes[idx].pos,pos);
		if (tmp < mindist)
		{
			mindist = tmp;
			minidx = idx;
		}
		uint dir = cmp(nodes[idx].pos,pos);
		idx = child(idx,dir);
	}

	Node<c_typ> n1 = nodes[minidx];
	double dist1 = dist(n1.pos,pos);
	double s1 = 1/dist1;

	c_typ sum1 = n1.count * s1;
	double sum2 = s1;

	for (uint i = 1; i <= _subs; i++)
	{
		uint neighbor_idx = neighbor(minidx,i);
		if (neighbor_idx == (uint)-1)
			continue;

		Node<c_typ> n2 = nodes[neighbor_idx];
		double dist2 = dist(n2.pos,pos);
		double s2 = 1/dist2;
		sum1 += n2.count * s2;
		sum2 += s2;

		std::cout << Node<c_typ>::to_string(*this,n2) << "dist: "<< dist2 << std::endl;
	}

	return sum1 / sum2;
}

template <typename c_typ>
uint CHist<c_typ>::neighbor(uint idx,uint dir) const
{
	uint opp = child_opposite(dir);
	uint c = child(idx,dir);
	if (c < _size && nodes[c].pos != nullptr)
	{
		return child_loop(c,opp);
	}

	uint p = parent(idx);
	uint pdir = get_dir(p,idx);
	if (pdir == opp)
		return p;

	if (pdir == dir)
	{
		while (true)
		{
			if (p == 0)
				return -1;
			idx = p;
			p = parent(idx);
			uint pdir = get_dir(p,idx);

			if (pdir == dir)
				return child_loop(p,opp);
		}
	}

	uint opp_pdir = child_opposite(pdir);
	c = child(p,dir);
	return child_loop(c,opp_pdir);

}

template <typename c_typ>
void CHist<c_typ>::dump() const
{
	std::cout << _subs << std::endl;
	for (uint i = 0; i < _subs; i++)
		std::cout << to_string(limits[i]) << "\n";
	dumpP(0);
}

template <typename c_typ>
void CHist<c_typ>::dumpP(uint idx) const
{
	uint child_l = child(idx,1);
	if (child_l < _size && nodes[child_l].pos != nullptr)
		dumpP(child_l);

	std::cout << to_string(nodes[idx].pos)
			  << "," << nodes[idx].count << std::endl;

	uint child_r = child(idx,2);
	if (child_r < _size && nodes[child_r].pos != nullptr)
		dumpP(child_r);
}

template <typename c_typ>
void CHist<c_typ>::print() const
{
	std::cout << "limits: ";
	for (uint i = 0; i < _subs; i++)
		std::cout << "i: " << (i+1) << "," << to_string(limits[i]) << " ";
	std::cout << std::endl;
	print(0,0);
	std::cout << std::endl;
}

template <typename c_typ>
void CHist<c_typ>::print(uint idx, uint d) const
{
	uint i;

	if (_size <= idx || nodes[idx].pos == nullptr)
		return;

	for (i=0; i<d; i++)
		printf("  ");
	printf("%i: ", idx);

	std::cout << "pos: " << to_string(nodes[idx].pos)
			  << " count: " << nodes[idx].count << std::endl;
	for (uint i = 1; i <= _subs; i ++)
	{
		print(child(idx,i),d+1);
	}
}

template <typename c_typ>
CHist<c_typ> CHist<c_typ>::copy() const
{
	CHist<c_typ> res = CHist<c_typ>(_size,_dimensions);

	auto size = sizeof(double) * _dimensions;

	for (uint i = 0; i < _size; i++)
	{
		memcpy(res.nodes[i].pos,nodes[i].pos,size);
		res.nodes[i].count = nodes[i].count;
	}

	for (uint i = 0; i < _subs; i++)
		memcpy(res.limits[i],limits[i],size);

	return res;
}

template <typename c_typ>
CHist<c_typ> CHist<c_typ>::merge(const CHist<c_typ> &h1, const CHist<c_typ> &h2)
{
	if (h1._size != h2._size)
		throw RuntimeException(TRACE_INFO,"CHists must be same size");
	if (h1._subs != h2._subs)
		throw RuntimeException(TRACE_INFO,"CHists must be same dimensions");

	CHist res = CHist(h1._size,h1._dimensions);

	auto size = sizeof(double)*res._dimensions;

	for (uint i = 0; i < h1._size; i++)
	{
		Node<c_typ> n1 = h1.nodes[i];
		double *n1pos = (double*)malloc(size);
		memcpy(n1pos,n1.pos,size);
		res.insert(n1pos,n1.count);

		Node<c_typ> n2 = h2.nodes[i];
		double *n2pos = (double*)malloc(size);
		memcpy(n2pos,n2.pos,size);
		res.insert(n2pos,n2.count);
	}


	for (uint i = 0; i < res._subs; i++)
		if ((i+1) == res.cmp(h1.limits[i],h2.limits[i]))
		{
			memcpy(res.limits[i],h2.limits[i],size);
		}
		else
		{
			memcpy(res.limits[i],h1.limits[i],size);
		}

	return res;
}

template <typename c_typ>
bool CHist<c_typ>::operator==(const CHist<c_typ> &other) const
{
	if (_size != other._size || _subs != other._subs ||
	    _count_elems != other._count_elems || _total_count != other._total_count)
		return false;

	for (uint i = 0; i < _subs; i++)
		if (!eq(limits[i],other.limits[i]))
			return false;

	for (uint i = 0; i < _size; i++)
		if (!eq(nodes[i],other.nodes[i]))
			return false;

	return true;
}

template <typename c_typ>
CHist<c_typ>& CHist<c_typ>::operator+=(const CHist<c_typ>& val)
{
	if (_size != val.size() && _dimensions != val.dimensions())
		throw RuntimeException(TRACE_INFO,"Wrong size or dimensions!");
	for (uint i = 0; i < _size; i++)
		nodes[i].count += val.nodes[i].count;
	return *this;
}

template <typename c_typ>
CHist<c_typ>& CHist<c_typ>::operator+=(const double& val)
{
	for (auto node : nodes)
		node.count += val;
	return *this;
}

template <typename c_typ>
CHist<c_typ>& CHist<c_typ>::operator-=(const double& val)
{
	for (auto node : nodes)
		node.count -= val;
	return *this;
}

template <typename c_typ>
CHist<c_typ>& CHist<c_typ>::operator*=(const double& val)
{
	for (auto node : nodes)
		node.count *= val;
	return *this;
}

template <typename c_typ>
CHist<c_typ>& CHist<c_typ>::operator/=(const double& val)
{
	for (auto node : nodes)
		node.count /= val;
	return *this;
}

template <typename c_typ>
std::ostream& operator<<(std::ostream& os, const CHist<c_typ>& chist)
{
	os << chist.print() << std::endl;
	return os;
}

template <typename c_typ>
std::string CHist<c_typ>::to_string(double * pos) const
{
	std::stringstream ss;
	if (_dimensions == 1)
	{
		ss << *pos;
	}
	else
	{
		ss << "(";
		for (uint i = 0; i < _dimensions; i++)
			ss << *(pos+i) << ",";
		ss.seekp(-1,std::ios_base::end);
		ss << ")";
	}
	return ss.str();
}

template <typename c_typ>
std::string Node<c_typ>::to_string(const CHist<c_typ> &h,Node<c_typ> n)
{
	std::stringstream ss;
	ss << "Node pos,count: " << h.to_string(n.pos) << "," << n.count << " ";
	return ss.str();
}

double get_count(const Node<double>& val)
{
	return val.count;
}

double get_count(const Node<CHist<double>>& val)
{
	return val.count.total_count();
}

void update_count(Node<double>& val, double c)
{
	val.count = c;
}

void update_count(Node<CHist<double>>& val, double c)
{
	double tc = val.count.total_count();

	for (auto i = val.count.begin(); i != val.count.end(); i++)
	{
		double ec = i->count;
		i->count = ec / tc * c;
	}
}

double get_count(const double& val)
{
	return val;
}

double get_count(const CHist<double>& val)
{
	return val.total_count();
}

void update_count(double& val, double c)
{
	val = c;
}

void update_count(CHist<double>& val, double c)
{
	double tc = val.total_count();

	for (auto i = val.begin(); i != val.end(); i++)
	{
		double ec = i->count;
		i->count = ec / tc * c;
	}
}
