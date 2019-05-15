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


#ifndef _OPENCOG_COVERTREE_H
#define _OPENCOG_COVERTREE_H

#include <ostream>
#include "Utils.h"

#include <boost/operators.hpp>

using namespace std;

namespace opencog
{

template <typename val_t>
class CoverTreeNode : boost::addable<CoverTreeNode<val_t>>
{
public:
	CoverTreeNode() {}
	CoverTreeNode(DVec x,val_t v) : pos(x) , value(v) {}
	DVec pos;
	val_t value;

	std::vector<int> children;

	bool operator!=(const CoverTreeNode<val_t>& other) const;
	bool operator==(const CoverTreeNode<val_t> & other) const;
};

template <typename val_t>
using CTNSeq = std::vector<CoverTreeNode<val_t>>;
template <typename val_t>
using CTNPSeq = std::vector<CoverTreeNode<val_t>*>;

typedef CoverTreeNode<double> CTN1;

template <typename val_t>
class CoverTree
{
protected:
	int _root_idx;
	int _root_level;
	int _elem_count;
	double _total_count;
	size_t _dims;
	//Is this a nearest Ancestor Tree
	bool _nat;
	std::vector<CoverTreeNode<val_t>> _nodes;

	/*
	 * Calculate the Distance between 2 nodes
	 */
	static double dist(const CoverTreeNode<val_t> & n1, const CoverTreeNode<val_t> & n2);

	/*
	 * Calculate the cover distance for a given level
	 */
	static double covdist(int level);

	/*
	 * Calculate the separation distance for a given level
	 */
	static double sepdist(int level);

	/*
	 * Get the Maximumdistance of a node to it's children
	 */
	double maxdist(const CoverTreeNode<val_t> &) const;

	/*
	 * Remove and return a Leaf Node thats a descendants of the given Node
	 */
	int popLeaf(CoverTreeNode<val_t> & n);

	/*
	 * This Helper Implements the actual algorithm for finding the
	 * NearestNeighbor
	 */
	const CoverTreeNode<val_t>* findNearestNeighbor_(const CoverTreeNode<val_t> & x,
			      								     const CoverTreeNode<val_t> & p,
											         const CoverTreeNode<val_t> * y) const;

	int findNearestNeighbor_(const CoverTreeNode<val_t> & x, int p, int y,
							 int level,int & ret_level, int & parent);

	/*
	 * Insert a Node x into the Node p at level:level;
	 */
	void insert(int node_idx, int & p_id, int level);

	/*
	 * Insert x Recursivly into p at level: level;
	 */
	void insert_rec(int node_idx, CoverTreeNode<val_t> & p, int level);

	/*
	 * Insert x Recursivly into p at level: level;
	 * in an ordered way , rebalance afterwards if neccesary
	 */
	//void insert_ord(CoverTreeNode<val_t> & x, CoverTreeNode<val_t> & p, int level);

	/*
	 * Rebalance according to Nearest ancestor requirments
	 */
	//void rebalance(CoverTreeNode<val_t> & p, CoverTreeNode<val_t> & x, int level);

	/*
	 * Main algorithm to recursivly rebalance starting from p;
	 */
	//void rebalance_rec(CoverTreeNode<val_t> & p, CoverTreeNode<val_t> & pq,
	//                   typename CTNSeq<val_t>::iterator & q,
	//                   CoverTreeNode<val_t> & x,
	//                   CTNSeq<val_t> & moveset, CTNSeq<val_t> & stayset,
	//                   int level);

public:


	CoverTree();
	CoverTree(int dims);
	CoverTree(CoverTreeNode<val_t> n,int dims);


	//CoverTreeNode<val_t> root() const {return _root;};
	int elem_count() const {return _elem_count;};
	size_t dims() const {return _dims;};
	double total_count() const {return _total_count;};
	const std::vector<CoverTreeNode<val_t>> & nodes() const {return _nodes;};

	/*
	 * Get a Vector of all positions
	 */
	DVecSeq get_posvec() const;

	/*
	 * Get the value at a given Position
	 */
	val_t get(DVec pos) const;

	/*
	 * This Function is the public Interface to find the NearestNeighbor
	 */
	const CoverTreeNode<val_t>* findNearestNeighbor(const DVec & pos) const;
	const CoverTreeNode<val_t>* findNearestNeighbor(const CoverTreeNode<val_t> & x) const;

	/*
	 * Insert a Node into the Tree
	 * Public Interface
	 */
	void insert(DVec,val_t);
	void insert(CoverTreeNode<val_t> x);


	std::string to_string() const;

	/*
	 * Print the Tree to std::cout
	 */
	void print() const;

	/*
	 * Returns a list of all Descendants
	 */
	void descendants(const CoverTreeNode<val_t> &, std::vector<int> & res) const;

	/*
	 * Merges 2 Trees together into a new one
	 */
	static CoverTree<val_t> merge(const CoverTree<val_t>&, const CoverTree<val_t>&);

	/*
	 * Check if the Tree is still valid.
	 */
	bool is_valid() const;
    bool is_valid_rec(int idx,int level) const;

	CoverTree<val_t>& operator+=(const CoverTree<val_t>&);

	bool operator!=(const CoverTree<val_t>& other) const;
	bool operator==(const CoverTree<val_t>& other) const;

	CoverTreeNode<val_t>& operator[](int idx)
	{
		return _nodes[idx];
	}

	CoverTreeNode<val_t> operator[](int idx) const
	{
		return _nodes[idx];
	}

	typename CTNSeq<val_t>::iterator begin() {return _nodes.begin();}
	typename CTNSeq<val_t>::iterator end() {return _nodes.end();}
	typename CTNSeq<val_t>::const_iterator begin() const {return _nodes.begin();}
	typename CTNSeq<val_t>::const_iterator end() const {return _nodes.end();}

	//CoverTree<val_t>& operator+=(const double&);
	//CoverTree<val_t>& operator-=(const double&);
	CoverTree<val_t>& operator*=(const double&);
	CoverTree<val_t>& operator/=(const double&);

	friend std::ostream& operator<<(std::ostream& os, const CoverTree<val_t>& t)
	{
		os << t.to_string() << std::endl;
		return os;
	}
};

double get_count(const double v);

double get_count(const CoverTree<double>& ct);

void update_count(double & v,const double n);

void update_count(CoverTree<double>& ct,const double n);

}

#endif // _OPENCOG_COVERTREE_H
