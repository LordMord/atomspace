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

#include <opencog/atoms/distvalue/CHistIter.h>

namespace opencog
{

typedef unsigned int uint;
typedef std::vector<double> DVec;

class CHist;

struct Node {
	double *pos;
	double count;

	static std::string to_string(const CHist&, Node);
};

class CHist
{
	friend iterator;
	friend const_iterator;
	//TODO: Size is also stored in the nodes vector
	//Should we use a vector? If so get ride of _size
	uint _size;
	uint _count_elems;
	uint _total_count;
	uint _subs;
	uint _dimensions;
	std::vector<double*> limits;
	std::vector<Node> nodes;

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
	void mergeNode(uint, double*, double);

	/*
	 * Insert a value into the Tree trying to fill it
	 * and rebalance if required
	 */
	void insertFill(double*,double);

	/*
	 * For when the Tree is alredy full.
	 * Find the closest Node to the provided value an merge them
	 */
	void insertMerge(double*,double);

	/*
	 * Helpers for the dump/print/insert function
	 */
	void dumpP(uint idx) const;
	void print(uint, uint) const;

	double* vecToArray(DVec) const;

	/*
	 * Helpers for calculation with Positions
	 */
	double* mul(double *, double) const;
	double* div(double *, double) const;
	double* add(double *, double *) const;
	bool eq(double *, double *) const;
	bool eq(const Node&, const Node&) const;

	/*
	 * Get Both the minimum and maximum counts
	 */
	DVec minmax_count() const;

public:

	CHist(uint,uint);

	~CHist();

	uint size() const {return _size;}
	uint dimensions() const {return _dimensions;}
	uint total_count() const {return _total_count;}
	uint count_bins() const {return _count_elems;}

	uint child_loop(uint,uint) const;
	uint next(uint,uint&) const;

	iterator begin() {return iterator(child_loop(0,1),2,*this);}
	iterator end() {return iterator(0,0,*this);}

	const_iterator begin() const {return const_iterator(child_loop(0,1),2,*this);}
	const_iterator end() const {return const_iterator(0,0,*this);}

	const_iterator cbegin() const {return const_iterator(child_loop(0,1),2,*this);}
	const_iterator cend() const {return const_iterator(0,0,*this);}

	/*
	 * Insert a Value into the Histogram
	 */
	void insert(DVec, double);
	void insert(double*,double);

	/*
	 * Provide a Position and get the Count at that Position
	 * throws Exception if Position is not in the Tree
	 */
	double get(DVec) const;
	double get(double*) const;

	/*
	 * Provide a Position and get the Count at that Position
	 * If the position is not in the Tree take the average of the
	 * closest point and one on the other side
	 */
	double get_avg(DVec) const;
	double get_avg(double*) const;

	/*
	 * Mirror the Histogram about the L-infinity mean
	 * TODO: Should this be called negate.
	 */
	CHist negate() const;

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
	 * Convert a Position to a String for printing
	 */
	std::string to_string(double*) const;

	bool operator==(const CHist&) const;
	bool operator!=(const CHist& other) const {return !(*this == other);};
};

}

#endif // _OPENCOG_CHIST_H
