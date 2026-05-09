// Exclude.h
// Copyright Justin Dufresne and Terrance Howard, 2013

#ifndef EXCLUDE_H
#define EXCLUDE_H

#include "JamesBond.h"
#include <limits.h>

class Exclude {
public:
	Exclude() : ind1(-1), ind2(-1) {}
	Exclude(int ind1, int ind2) : ind1(ind1), ind2(ind2) {}
	bool operator==(const Exclude& e) const;
	bool operator!=(const Exclude& e) const;
	void print();
	int ind1;
	int ind2;
};

class Node {
public:
	Node(int index);
	void clearTree();
	int makeTree(Node** particles, Bond* bonds, int2* bondMap, int bondstart, int bondend);
	void add(Node* n);
	bool inTree;
	int index;
	int cap;
	int numBonds;
	Node** bonds;
};

// makeExcludes(Bond* bonds, int* bondMap, int num, int numBonds, String exList)
// @param    list of sorted cell bonds; corresponding bond map; number of particles; number of bonds;
//           string formated like so "EXCLUDE 1-2 1-3 1-4"; number of excludes
// @return   Array of Excludes
// This algorithm finds the central particle in every bond tree,
// then creates a list of exclusions for the particle pairs 
// defined in exList. For example, 1-2 means that there should
// be an exclusion between the central particle and every 
// particle it is directly bonded to. 1-3 means that there should
// be an exclusion between the central particle and every particle
// it is two bonds away from
Exclude* makeExcludes(Bond* bonds, int2* bondMap, int num, int numBonds,
		String exList, int& numExcludes);
void getExcludes(int root, Node* curr, Exclude* result, int depth, int& capacity,
		int& numExcludes, bool sentinel, bool* done);

#endif
