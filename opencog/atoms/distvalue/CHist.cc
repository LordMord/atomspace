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

#include <boost/algorithm/string.hpp>

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
	if (ds1 == nullptr && ds2 == nullptr)
		return true;
	if (ds1 == nullptr || ds2 == nullptr)
		return false;
	for (uint i = 0; i < _dimensions; i++)
		if (!is_approx_eq_ulp(*(ds1+i),*(ds2+i),24))
			return false;

	return true;
}

template <typename c_typ>
bool CHist<c_typ>::eq(const Node<c_typ> &n1,const Node<c_typ> &n2) const
{
	return (eq(n1.pos,n2.pos) && n1.value == n2.value);
}

template <typename c_typ>
CHist<c_typ>::CHist(uint s,uint d)
	: _size(s) , _count_elems(0) , _total_count(0)
	, _subs(pow(2,d)) , _dimensions(d)
{
	_levels = sizeToLevels(s,_subs);
	if (_size != levelsToSize(_levels,_subs))
		throw RuntimeException(TRACE_INFO,"Not a prooper size.");

	if (d >= 32)
		throw RuntimeException(TRACE_INFO,"More then 31 Dimensions not supported.");
	nodes.resize(s);
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
CHist<c_typ>::CHist(const CHist<c_typ> & other)
	: _size(other._size) , _count_elems(other._count_elems)
	, _total_count(other._total_count) , _subs(other._subs)
	, _levels(other._levels) , _dimensions(other._dimensions)
{
	auto size = sizeof(double)*_dimensions;

	nodes.resize(_size);
	for (uint i = 0; i < _size; i++)
	{
		if (other.nodes[i].pos == nullptr)
			continue;

		nodes[i].pos = (double*)malloc(size);
		if(nodes[i].pos == nullptr)
			throw RuntimeException(TRACE_INFO,"Malloc Failed");
		memcpy(nodes[i].pos,other.nodes[i].pos,size);

		nodes[i].value = other.nodes[i].value;
	}

	limits.resize(_subs);
	for (uint i = 0; i < _subs; i++)
	{
		limits[i] = (double*)malloc(size);
		if(limits[i] == nullptr)
			throw RuntimeException(TRACE_INFO,"Malloc Failed");
		memcpy(limits[i],other.limits[i],size);
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
uint CHist<c_typ>::levelsToSize(uint levels,uint subs)
{
	return (std::pow(subs,levels) - 1) / (subs - 1);
}

template <typename c_typ>
uint CHist<c_typ>::sizeToLevels(uint size,uint subs)
{
	return log(size * (subs - 1) + 1) / log(subs);
}

template <typename c_typ>
CHist<c_typ>& CHist<c_typ>::operator=(const CHist<c_typ>& other)
{
	if (this == &other)
		return *this;

	for (Node<c_typ> elem : nodes)
		delete[] elem.pos;
	for (auto elem : limits)
		delete[] elem;


	for (uint i = 0; i < _size ; i++)
		nodes[i].pos = nullptr;
	for (uint i = 0; i < _subs ; i++)
		limits[i] = nullptr;

	_size       = other._size;
	_count_elems = other._count_elems;
	_total_count = other._total_count;
	_subs       = other._subs;
	_dimensions = other._dimensions;
	_levels     = other._levels;

	auto size = sizeof(double)*_dimensions;

	nodes.resize(_size);
	for (uint i = 0; i < _size; i++)
	{
		if (other.nodes[i].pos == nullptr)
			continue;

		nodes[i].pos = (double*)malloc(size);
		if(nodes[i].pos == nullptr)
			throw RuntimeException(TRACE_INFO,"Malloc Failed");
		memcpy(nodes[i].pos,other.nodes[i].pos,size);

		nodes[i].value = other.nodes[i].value;
	}

	limits.resize(_subs);
	for (uint i = 0; i < _subs; i++)
	{
		limits[i] = (double*)malloc(size);
		if(limits[i] == nullptr)
			throw RuntimeException(TRACE_INFO,"Malloc Failed");
		memcpy(limits[i],other.limits[i],size);
	}

	return *this;
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
int CHist<c_typ>::height(uint idx) const
{
	if (idx >= _size || nodes[idx].pos == nullptr)
		return 0;
	int max = 0;
	for (uint i = 1; i <= _subs; i ++)
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
template <typename c_typ>
void CHist<c_typ>::move_to(uint idx,uint towards)
{
	if (nodes[idx].pos == nullptr || idx >= _size)
		return;

	nodes[towards] = nodes[idx];
	nodes[idx].pos = nullptr;
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
}

template <typename c_typ>
uint CHist<c_typ>::opposite_dir(uint dir) const
{
	//Old didn't match cmp
	//return child_idx + _subs - 1 - 2 * (child_idx - 1);
	return ((dir - 1) + _subs / 2) % _subs + 1;
}

template <typename c_typ>
void CHist<c_typ>::rotate_up(uint idx)
{
	if (nodes[idx].pos == nullptr || idx >= _size)
		return;

	if (idx == 0)
		throw RuntimeException(TRACE_INFO,"Can't rotate_up root node.");

	uint p = parent(idx);
	uint cdir = get_dir(p,idx);
	uint cdir_opp = opposite_dir(cdir);

	//Move down the parrents Children except idx
	for (uint i = 1; i <= _subs; i ++)
	{
		if (i == cdir)
			continue;
		if (i == cdir_opp)
			shift_down(child(p,i),cdir_opp);
		else
			move_to(child(p,i),child(child(p,cdir_opp),i));
	}
	//Move down the parrent
	nodes[child(p,cdir_opp )] = nodes[p];

	//Move child_opp to other side of tree
	move_to(child(idx,cdir_opp),child(child(p,cdir_opp),cdir));
	//Set it to 0 in original location
	nodes[child(idx,cdir_opp)].pos = nullptr;

	//Move idx Upwards
	nodes[p] = nodes[idx];
	//Incase there is nothing to shift up into this position
	nodes[idx].pos = nullptr;

	//Move Idx's children Upwards except child_opp
	for (uint i = 1; i <= _subs; i ++)
	{
		if (i == cdir_opp)
			continue;
		move_to(child(idx,i),child(p,i));
	}

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
void CHist<c_typ>::rebalance(uint idx,uint dir)
{
	while (true)
	{
		uint c1 = child(idx,dir);
		uint c2 = child(idx,opposite_dir(dir));
		int h1 = height(c1);
		int h2 = height(c2);

		if (abs(h1-h2) > 1)
		{
			uint max_idx;
			if (h1>h2)
				max_idx = c1;
			else
				max_idx = c2;
			uint child_idx = get_dir(idx,max_idx);
			uint child_opp = opposite_dir(child_idx);

			if (child(max_idx,child_opp) == max_height_child(max_idx))
				rotate_up(child(max_idx,child_opp));
			rotate_up(max_idx);
		}

		//End if we reached the root
		if (0 == idx)
			return;

		//Check if we need to reblance the parent next
		uint p = parent(idx);
		dir = get_dir(p,idx);
		idx = p;
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
		uint c1 = child(idx,dir);
		uint c2 = child(idx,opposite_dir(dir));
		int h1 = height(c1);
		int h2 = height(c2);

		if (abs(h1-h2) > 0)
		{
			if (h1>h2)
				rotate_up(c1);
			else
				rotate_up(c2);
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
void CHist<c_typ>::mergeNode(uint idx, double * pos,const c_typ & c)
{
	Node<c_typ> &n = nodes[idx];
	double nc = get_count(n);
	double ac = get_count(c);
	n.pos = div(add(mul(n.pos,nc),mul(pos,ac)),(nc + ac));
	merge_count(n.value,c);

	delete[] pos;
}


template <typename c_typ>
double * CHist<c_typ>::vecToArray(const DVec& vec) const
{
	if (vec.size() != _dimensions)
		throw RuntimeException(TRACE_INFO,"Vector needs to be the same lenght as the number of dimensions!");

	auto size = sizeof(double)*_dimensions;
	double * arr = (double*)malloc(size);
	memcpy(arr,&vec[0],size);

	return arr;
}

template <typename c_typ>
DVec CHist<c_typ>::arrayToVec(const double * arr) const
{
	return DVec(arr,arr + _dimensions);
}

template <typename c_typ>
void CHist<c_typ>::insert(DVec posv,const c_typ & c)
{
	return insert(vecToArray(posv),c);
}

template <typename c_typ>
void CHist<c_typ>::insert(double * pos,const c_typ & c)
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

	if (_size == _count_elems)
		insertMerge(pos,c);
	else
		insertFill(pos,c);
	_total_count += get_count(c);
}

template <typename c_typ>
void CHist<c_typ>::insertFill(double * pos, const c_typ & c_val)
{
	uint i;
	double mindist = std::numeric_limits<double>::infinity();
	uint minidx = -1;
	for (i = 0; i < _size; )
	{
		if (nodes[i].pos == nullptr)
		{
			_count_elems++;
			nodes[i].pos = pos;
			nodes[i].value = c_val;

			if (i == 0)
				return;
			uint p = parent(i);
			uint dir = get_dir(p,i);
			rebalance(p,dir);
			return;
		}

		if (eq(nodes[i].pos,pos))
		{
			merge_count(nodes[i].value,c_val);
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
			uint child_idx = cmp(nodes[i].pos,pos);
			i = child(i,child_idx);
		}
	}

	uint idx = parent(i);
	uint p = parent(idx);
	uint dir = get_dir(p,idx);
	make_space(p,dir);

	if (nodes[idx].pos == nullptr)
	{
		_count_elems++;
		nodes[idx].pos = pos;
		nodes[idx].value = c_val;
		return;
	}

	mergeNode(minidx,pos,c_val);
}

template <typename c_typ>
void CHist<c_typ>::insertMerge(double * pos,const c_typ & c)
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

		if (tmp < mindist)
		{
			mindist = tmp;
			minidx = i;
		}

		uint child_idx = cmp(nodes[i].pos,pos);
		uint child_i = child(i,child_idx);
		if (nodes[child_i].pos == nullptr)
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
uint CHist<c_typ>::child_loop(uint idx,uint dir,bool check) const
{
	while (true)
	{
		uint child_idx = child(idx,dir);
		if (child_idx >= _size || (check && nodes[child_idx].pos == nullptr))
			return idx;
		idx = child_idx;
	}
}

template <typename c_typ>
uint CHist<c_typ>::next(uint idx,uint &dir) const
{
	uint res = nextP(idx,dir);
	if (nodes[res].pos == nullptr)
		return next(res,dir);
	return res;
}

template <typename c_typ>
uint CHist<c_typ>::nextP(uint idx,uint &dir) const
{
	if (dir == 0)
		return idx;

	while (dir <= _subs)
	{
		uint child_idx = child(idx,dir);
		//TODO: Improvment if one child is > size then all should be
		if (child_idx < _size)
		{
			dir = 2;
			return child_loop(child_idx,1,false);
		}
		dir++;
	}

	uint p = parent(idx);
	while (child(p,_subs) == idx)
	{
		if (p == 0 && child(p,_subs) == idx)
		{
			dir = 0;
			return 0;
		}
		idx = p;
		p = parent(idx);
	}

	dir = 1 + get_dir(p,idx);

	if (dir-1 == _subs/2)
		return p;
	else
		return nextP(p,dir);
}

template <typename c_typ>
DVecSeq CHist<c_typ>::get_posvvec() const
{
	auto tmp = get_posavec();
	DVecSeq res;
	std::transform(tmp.begin(),tmp.end(),std::back_inserter(res),
	               [this](const double * arr) { return this->arrayToVec(arr); });
	return res;
}

template <typename c_typ>
std::vector<double*> CHist<c_typ>::get_posavec() const
{
	std::vector<double*> res;
	for (auto elem : nodes)
		if (elem.pos != nullptr)
			res.push_back(elem.pos);
	return res;
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
		if (eq(nodes[idx].pos,pos))
			return nodes[idx].value;
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
			return nodes[idx].value;

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
	uint mdir = cmp(pos,n1.pos);
	double s1 = 1/dist1;

	c_typ sum1 = n1.value * s1;
	double sum2 = s1;

	for (uint i = 1; i <= _subs; i++)
	{
		if (i == mdir)
			continue;

		uint neighbor_idx = neighbor(minidx,i);
		if (neighbor_idx == (uint)-1)
			continue;

		Node<c_typ> n2 = nodes[neighbor_idx];
		double dist2 = dist(n2.pos,pos);
		double s2 = 1/dist2;
		sum1 += n2.value * s2;
		sum2 += s2;
	}
	return sum1 / sum2;
}

template <>
CHist<double> CHist<CHist<double>>::get_avg(double * pos) const
{
	double mindist = std::numeric_limits<double>::infinity();
	uint minidx = -1;

	uint idx = 0;
	while (idx < _size && nodes[idx].pos != nullptr)
	{
		if (eq(nodes[idx].pos,pos))
			return nodes[idx].value;

		double tmp = dist(nodes[idx].pos,pos);
		if (tmp < mindist)
		{
			mindist = tmp;
			minidx = idx;
		}
		uint dir = cmp(nodes[idx].pos,pos);
		idx = child(idx,dir);
	}

	Node<CHist<double>> n1 = nodes[minidx];
	double dist1 = dist(n1.pos,pos);
	uint mdir = cmp(pos,n1.pos);
	double s1 = 1/dist1;

	auto dims = n1.value.dimensions();
	uint levels = _levels + n1.value.levels();
	auto size = levelsToSize(levels,pow(2,dims));
	CHist<double> sum1 = CHist<double>(size,dims);
	sum1 += n1.value * s1;
	double sum2 = s1;

	for (uint i = 1; i <= _subs; i++)
	{
		if (i == mdir)
			continue;

		uint neighbor_idx = neighbor(minidx,i);
		if (neighbor_idx == (uint)-1)
			continue;

		Node<CHist<double>> n2 = nodes[neighbor_idx];
		double dist2 = dist(n2.pos,pos);
		double s2 = 1/dist2;
		sum1 += n2.value * s2;
		sum2 += s2;
	}
	return sum1 / sum2;
}
template <typename c_typ>
uint CHist<c_typ>::neighbor(uint idx,uint dir) const
{
	uint opp = opposite_dir(dir);
	uint c = child(idx,dir);
	if (c < _size && nodes[c].pos != nullptr)
	{
		return child_loop(c,opp);
	}

	if (idx == 0)
		return -1;

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

	uint opp_pdir = opposite_dir(pdir);
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
			  << "," << nodes[idx].value << std::endl;

	uint child_r = child(idx,2);
	if (child_r < _size && nodes[child_r].pos != nullptr)
		dumpP(child_r);
}

template <typename c_typ>
std::string CHist<c_typ>::to_string() const
{
	std::stringstream ss;
	ss << "\nlimits: ";

	for (uint i = 0; i < _subs; i++)
		ss << "i: " << (i+1) << "," << to_string(limits[i]) << " ";

	ss << std::endl
	   << to_string(0,0)
	   << std::endl;

	return ss.str();
}

template <typename c_typ>
std::string CHist<c_typ>::to_string(uint idx, uint d) const
{
	std::stringstream ss;
	uint i;

	if (_size <= idx || nodes[idx].pos == nullptr)
		return "";

	for (i=0; i<d; i++)
		ss << "  ";
	ss << idx << ": ";

	std::stringstream sstmp;
	sstmp << nodes[idx].value;
	std::string tmp = sstmp.str();
	boost::replace_all(tmp,"\n","\n    ");

	ss << "pos: " << to_string(nodes[idx].pos)
	   << " count: " << tmp  << std::endl;

	for (uint i = 1; i <= _subs; i ++)
	{
		ss << to_string(child(idx,i),d+1);
	}

	return ss.str();
}

template <typename c_typ>
void CHist<c_typ>::print() const
{
	std::cout << to_string();
}

template <typename c_typ>
CHist<c_typ> CHist<c_typ>::copy() const
{
	CHist<c_typ> res = CHist<c_typ>(_size,_dimensions);

	auto size = sizeof(double) * _dimensions;

	for (uint i = 0; i < _size; i++)
	{
		memcpy(res.nodes[i].pos,nodes[i].pos,size);
		res.nodes[i].value = nodes[i].value;
	}

	for (uint i = 0; i < _subs; i++)
		memcpy(res.limits[i],limits[i],size);

	return res;
}

template <typename c_typ>
CHist<c_typ> CHist<c_typ>::merge(const CHist<c_typ> &h1, const CHist<c_typ> &h2)
{
	//if (h1._size != h2._size)
    //	throw RuntimeException(TRACE_INFO,"CHists must be same size");
	if (h1._subs != h2._subs)
		throw RuntimeException(TRACE_INFO,"CHists must be same dimensions");

	CHist res = CHist(std::max(h1._size,h2._size),h1._dimensions);

	auto size = sizeof(double)*res._dimensions;

	auto it1 = h1.nodes.begin();
	auto it2 = h2.nodes.begin();
	auto end1 = h1.nodes.end();
	auto end2 = h2.nodes.end();
	bool first = false;

	while (it1 != end1 && it2 != end2)
	{
		first = !first;
		if (first)
		{
			if (it1 == end1)
				continue;
			Node<c_typ> n1 = *it1;
			if (n1.pos != nullptr)
			{
				double *n1pos = (double*)malloc(size);
				memcpy(n1pos,n1.pos,size);
				res.insert(n1pos,n1.value);
			}
			it1++;
		}
		else
		{
			if (it2 == end2)
				continue;
			Node<c_typ> n2 = *it2;
			if (n2.pos != nullptr)
			{
				double *n2pos = (double*)malloc(size);
				memcpy(n2pos,n2.pos,size);
				res.insert(n2pos,n2.value);
			}
			it2++;
		}
	}


	if (h1._count_elems != 0 && h2._count_elems != 0)
		for (uint i = 0; i < res._subs; i++)
		{
			if (h1._count_elems == 0)
				memcpy(res.limits[i],h2.limits[i],size);
			else if (h2._count_elems == 0)
				memcpy(res.limits[i],h1.limits[i],size);
			else if ((i+1) == res.cmp(h1.limits[i],h2.limits[i]))
				memcpy(res.limits[i],h2.limits[i],size);
			else
				memcpy(res.limits[i],h1.limits[i],size);
		}

	return res;
}

template <typename c_typ>
bool CHist<c_typ>::operator==(const CHist<c_typ> &other) const
{
	//std::cout << "Operator==\n";
	if (_size != other._size || _subs != other._subs ||
	    _count_elems != other._count_elems || _total_count != other._total_count)
		return false;

	//std::cout << "Check1\n";
	//std::cout << *this
//			  << other;

	for (uint i = 0; i < _subs; i++)
		if (!eq(limits[i],other.limits[i]))
			return false;
	//std::cout << "Check2\n";

	for (uint i = 0; i < _size; i++)
		if (!eq(nodes[i],other.nodes[i]))
			return false;
	//std::cout << "Check3\n";

	return true;
}

template <typename c_typ>
CHist<c_typ>& CHist<c_typ>::operator+=(const CHist<c_typ>& val)
{
	*this = merge(*this,val);
	return *this;
}

template <typename c_typ>
CHist<c_typ>& CHist<c_typ>::operator+=(const double& val)
{
	_total_count = 0;
	for (auto& node : nodes)
	{
		node.value += val;
		_total_count += get_count(node);
	}
	return *this;
}

template <typename c_typ>
CHist<c_typ>& CHist<c_typ>::operator-=(const double& val)
{
	_total_count = 0;
	for (auto& node : nodes)
	{
		node.value -= val;
		_total_count += get_count(node);
	}
	return *this;
}

template <typename c_typ>
CHist<c_typ>& CHist<c_typ>::operator*=(const double& val)
{
	for (auto& node : nodes)
		node.value *= val;
	_total_count *= val;
	return *this;
}

template <typename c_typ>
CHist<c_typ>& CHist<c_typ>::operator/=(const double& val)
{
	for (auto& node : nodes)
		node.value /= val;
	_total_count /= val;
	return *this;
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
	if (n.pos == nullptr)
		return "Empty Node ";
	std::stringstream ss;
	ss << "Node pos,count: " << h.to_string(n.pos) << "," << n.value << " ";
	return ss.str();
}

namespace opencog
{

double get_count(const Node<double>& val)
{
	return val.value;
}

double get_count(const Node<CHist<double>>& val)
{
	return val.value.total_count();
}

void update_count(Node<double>& val, double c)
{
	val.value = c;
}

void update_count(Node<CHist<double>>& val, double c)
{
	double tc = val.value.total_count();

	for (auto i = val.value.begin(); i != val.value.end(); i++)
	{
		double ec = i->value;
		i->value = ec / tc * c;
	}
}

void merge_count(double& val, double c)
{
	val += c;
}
void merge_count(CHist<double>& val, CHist<double> c)
{
	val = CHist<double>::merge(val,c);
}
void merge_count(Node<double>& val, double c)
{
	val.value += c;
}
void merge_count(Node<CHist<double>>& n1, Node<CHist<double>> n2)
{
	n1.value = CHist<double>::merge(n1.value,n2.value);
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
		double ec = i->value;
		i->value = ec / tc * c;
	}
}

std::string to_string(const std::vector<double>& v)
{
	std::stringstream ss;
	for (auto e : v)
		ss << e << ",";
	return ss.str();
}

}


