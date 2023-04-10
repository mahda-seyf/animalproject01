#include<iostream>

using namespace std;

class QuarterTree {

public:
	Virus_node *Root;
public:
	QuarterTree() {
		Root = NULL;
	}
	void Virus_Insert(char * word);
	Virus_node* createNode(char val);
	int Virus_Count(Virus_node*node);
	int Virus_input(Virus_node*  Root);
	void Virus_tvs(char * sequence);
};

class Virus_node {
public:
	char Data;
	Virus_node *aChild;
	Virus_node *gChild;
	Virus_node *cChild;
	Virus_node *tChild;
};

