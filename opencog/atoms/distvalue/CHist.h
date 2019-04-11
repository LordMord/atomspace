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
#include <string>
#include <ostream>

#include <boost/operators.hpp>

#include <opencog/atoms/distvalue/CHistIter.h>

namespace opencog
{

typedef unsigned int uint;
typedef std::vector<double> DVec;
typedef std::vector<DVec> DVecSeq;

template <typename c_typ>
class CHist;

template <typename c_typ>
struct Node {
	double *pos;
	c_typ value;

	static std::string to_string(const CHist<c_typ>&, Node);
};

template <typename c_typ>
class CHist : boost::arithmetic2<CHist<c_typ>,double>
			, boost::addable1<CHist<c_typ>>
{
	friend iterator<c_typ>;
	friend const_iterator<c_typ>;
	//TODO: Size is also stored in the nodes vector
	//Should we use a vector? If so get ride of _size
	uint _size;
	uint _count_elems;
	uint _total_count;
	uint _subs;
	uint _dimensions;
	std::vector<double*> limits;
	std::vector<Node<c_typ>> nodes;

	/*
	 * Return the Direction from one Pos to the Other
	 */
	uint cmp(double*,double*) const;

	/*
	 * Calculate the distance between 2 positions
	 */
	double dist(double*,double*) const;

	/*
	 * Given the index of a node
	 * find the index of it's parent
	 */
	uint parent(uint) const;

	/*
	 * Given the index of a node
	 * and a direction to the child
	 * find the index of this child
	 */
	uint child(uint, uint) const;

	/*
	 *  Find the Neighbor (closest other node)
	 *  in a given direction
	 */
	uint neighbor(uint idx,uint dir) const;

	/*
	 * Given the index of a node
	 * find the height of the sub-tree
	 * starting at the index
	 */
	uint height(uint) const;

	/*
	 * Given a Node and one of it's Children
	 * Calculate the Direction from the Parent to the Child
	 */
	uint get_dir(uint,uint) const;

	/*
	 * Given a Direction to a Child
	 * Find the Opposite Direction
	 */
	uint child_opposite(uint) const;

	/*
	 * Given 2 indices move the sub-tree starting at the first
	 * to the location of the second. Overwrites anything at
	 * the target location.
	 * Target Location should not be of the sub-tree to be moved.
	 */
	void move_to(uint, uint);

	/*
	 * Move down the subtree at the given Index
	 * in the given direction
	 */
	void shift_down(uint, uint);

	/*
	 * Given the Index of a Node. Rotate the tree such that
	 * the Node ends up in the position of it's parent.
	 */
	void rotate_up(uint);

	/*
	 * Check if the Tree needs to be rebalanced at the given Index
	 * and all it's parents. And if needed does the rebalancing.
	 */
	void rebalance(uint);

	/*
	 * If we want to add a node but there is no space in the tree anymore
	 * see if we can rotate the tree such that a space opens up in the
	 * requried position without destroying the order.
	 */
	void make_space(uint,uint);

	/*
	 * Helper function that given an Index stores the min and max
	 * height of the Nodes Childrend in argument 2 and 3
	 * and returns the Index of the Child with max height
	 */
	uint min_max_heights(uint, uint*, uint*) const;

	/*
	 * Helper function to find the Index of the Child of the given Node
	 * that has the highest height
	 */
	uint max_height_child(uint) const;

	/*
	 * Given an Index a postiona and count
	 * Merge the Position of the Node with the given
	 * weighted by the counts
	 */
	void mergeNode(uint, double*,const c_typ&);

	/*
	 * Insert a value into the Tree trying to fill it
	 * and rebalance if required
	 */
	void insertFill(double*,const c_typ&);

	/*
	 * For when the Tree is alredy full.
	 * Find the closest Node to the provided value an merge them
	 */
	void insertMerge(double*,const c_typ&);

	/*
	 * Helpers for the dump/print/insert function
	 */
	void dumpP(uint idx) const;
	std::string to_string(uint, uint) const;

	/*
	 * Helpers for calculation with Positions
	 */
	double* mul(double *, double) const;
	double* div(double *, double) const;
	double* add(double *, double *) const;
	bool eq(double *, double *) const;
	bool eq(const Node<c_typ>&, const Node<c_typ>&) const;

	/*
	 * Get Both the minimum and maximum counts
	 */
	DVec minmax_count() const;

	/*
	 * Helper for Next Function
	 */
	uint nextP(uint,uint&) const;

public:

	CHist(uint s = 0,uint d = 0);

	CHist(const CHist &);

	~CHist();

	uint size() const {return _size;}
	uint dimensions() const {return _dimensions;}
	uint total_count() const {return _total_count;}
	uint count_bins() const {return _count_elems;}

	uint child_loop(uint,uint) const;
	uint next(uint,uint&) const;

	iterator<c_typ> begin()
	{
		if (_count_elems == 0)
			return end();
		return iterator<c_typ>(child_loop(0,1),2,*this);
	}
	iterator<c_typ> end() {return iterator<c_typ>(0,0,*this);}

	const_iterator<c_typ> begin() const
	{
		if (_count_elems == 0)
			return end();
		return const_iterator<c_typ>(child_loop(0,1),2,*this);
	}

	const_iterator<c_typ> end() const
	{return const_iterator<c_typ>(0,0,*this);}

	const_iterator<c_typ> cbegin() const
	{
		if (_count_elems == 0)
			return cend();
		return const_iterator<c_typ>(child_loop(0,1),2,*this);
	}

	const_iterator<c_typ> cend() const
	{return const_iterator<c_typ>(0,0,*this);}

	double* vecToArray(const DVec&) const;
	DVec arrayToVec(const double *) const;

	/*
	 * Insert a Value into the Histogram
	 */
	void insert(DVec, const c_typ&);
	void insert(double*,const c_typ&);

	/*
	 * Get all positions.
	 */
	DVecSeq get_posvvec() const;
	std::vector<double*> get_posavec() const;

	/*
	 * Provide a Position and get the Count at that Position
	 * throws Exception if Position is not in the Tree
	 */
	c_typ get(DVec) const;
	c_typ get(double*) const;

	/*
	 * Provide a Position and get the Count at that Position
	 * If the position is not in the Tree take the average of the
	 * closest point and one on the other side
	 */
	c_typ get_avg(DVec) const;
	c_typ get_avg(double*) const;

	/*
	 * Mirror the Histogram about the L-infinity mean
	 * TODO: Should this be called negate.
	 */
	//CHist negate() const;

	/*
	 * Merge 2 Histograms into 1
	 */
	static CHist merge(const CHist&, const CHist&);

	/*
	 * Copty a Histogram
	 */
	CHist copy() const;

	/*
	 * Dump the Histogram for displaying it via Python Script
	 */
	void dump() const;

	/*
	 * Print the Histogram Tree
	 */
	void print() const;

	/*
	 * Convert the Histogram to a Tree
	 */
	std::string to_string() const;

	/*
	 * Convert a Position to a String for printing
	 */
	std::string to_string(double*) const;

	bool operator==(const CHist<c_typ>&) const;
	bool operator!=(const CHist<c_typ>& other) const {return !(*this == other);};

	CHist<c_typ>& operator+=(const CHist<c_typ>&);

	CHist<c_typ>& operator+=(const double&);
	CHist<c_typ>& operator-=(const double&);
	CHist<c_typ>& operator*=(const double&);
	CHist<c_typ>& operator/=(const double&);

	CHist<c_typ>& operator=(const CHist<c_typ>& other);

	Node<c_typ>& operator[](int idx)
	{
		return nodes[idx];
	}

	Node<c_typ> operator[](int idx) const
	{
		return nodes[idx];
	}

	friend std::ostream& operator<<(std::ostream& os, const CHist<c_typ>& chist)
	{
		os << chist.to_string() << std::endl;
		return os;
	}
};

template class CHist<double>;
template class CHist<CHist<double>>;

template class Node<double>;
template class Node<CHist<double>>;

double get_count(const double& val);
double get_count(const CHist<double>& val);
double get_count(const Node<double>& val);
double get_count(const Node<CHist<double>>& val);

void update_count(double& val, double c);
void update_count(CHist<double>& val, double c);
void update_count(Node<double>& val, double c);
void update_count(Node<CHist<double>>& val, double c);

void merge_count(double& val, double c);
void merge_count(CHist<double>& val, CHist<double>c);
void merge_count(Node<double>& val, double c);
void merge_count(Node<CHist<double>>& val, Node<CHist<double>>c);

}

#endif // _OPENCOG_CHIST_H
