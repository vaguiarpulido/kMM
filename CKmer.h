/*
 * Kmer.h
 *
 *  Created on: 23 de sept. de 2015
 *      Author: vanessa
 */
#include<stdio.h>
#include <iostream>
#include <string>
#include <map>
using namespace std;

#ifndef CKMER_H_
#define CKMER_H_

class CKmer {
public:
	CKmer(int size);
	virtual ~CKmer();
	bool generateKmerList();
private:
	int order;
	char next(char current);
	void successor(string* current, int* pos);
	map<string, int> kmerList;
};

#endif /* CKMER_H_ */
