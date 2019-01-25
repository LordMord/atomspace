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

namespace opencog
{

typedef unsigned int uint;

struct Node {
	double *pos;
	double count;
};

class CHist
{
	uint size;
	uint count;
	uint subs;
	uint dimensions;
	std::vector<double*> limits;
	std::vector<Node> nodes;
	uint cmp(double*,double*);
	double dist(double*,double*);

	/*
	 * Given the index of a node
	 * find the index of it's parent
	 */
	uint parent(uint);

	/*
	 * Given the index of a node
	 * and a direction to the child
	 * find the index of this child
	 */
	uint child(uint, uint);

	/*
	 * Given the index of a node
	 * find the height of the sub-tree
	 * starting at the index
	 */
	uint height(uint);

	/*
	 * Given a Direction to a Child
	 * Find the Opposite Direction
	 */
	uint child_opposite(uint);

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
	uint min_max_heights(uint, uint*, uint*);

	/*
	 * Helper function to find the Index of the Child of the given Node
	 * that has the highest height
	 */
	uint max_height_child(uint);

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
	void dumpP(uint idx);
	void print(uint, uint);
	void insertP(double *, double);

	/*
	 * Convert a Position to a String for printing
	 */
	std::string to_string(double*);

	/*
	 * Helpers for calculation with Positions
	 */
	double* mul(double *, double);
	double* div(double *, double);
	double* add(double *, double *);

public:
	CHist(uint,uint);

	~CHist();

	/*
	 * Insert a Value into the Histogram
	 */
	void insert(std::vector<double>, double);

	/*
	 * Merge 2 Histograms into 1
	 */
	static CHist merge(const CHist, const CHist);

	/*
	 * Dump the Histogram for displaying it via Python Script
	 */
	void dump();

	/*
	 * Print the Histogram Tree
	 */
	void print();
};

}

#endif // _OPENCOG_CHIST_H
