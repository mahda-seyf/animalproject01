#include <fstream>
#include <iostream>
#include <vector>
#include <string.h>
#include "Virus.h"

using namespace std;

bool get_Vrs(char * file_name, vector<string> &r, int vsize) {
	ifstream in;
	in.open(file_name);
	if (!in.is_open()) {
		cout << "couldt\t";
		return false;
	}
	char * word = new char[(vsize + 1)]; 
	while (in.peek() != EOF)
	{
		in.getline(word, (vsize + 1), '\n'); 
		r.push_back(word);
	}
	in.clear();
	in.close();
	
	word = NULL; 
	delete word;
	
	return true;
}

void ch_virus(char input) {
	while (input != 'y' && input != 'n' && input != 'Y' && input != 'N')
	{
		cout << "Enter y or n " << endl;
		cin >> input;
		cin.clear();
		cin.ignore(100, '\n');
	}
}

void Virus_gen :: InsertNode(char * word) { 

	if (!Root) {	
		Root = createNode(word[0]);
	}
	Virus_genome *ptr; 
	ptr = Root;
	int i = 1; 

	while (i < strlen(word)) {
		if (word[i] == 'a') { 
			if (ptr->aChild) {
				ptr = ptr->aChild;
			} 
			else {
				ptr->aChild = createNode(word[i]);
				ptr = ptr->aChild;
			} 
			i++;
		}
		else if (word[i] == 'g') {
			if (ptr->gChild) {
				ptr = ptr->gChild;
			}	
			else {
				ptr->gChild = createNode(word[i]);
				ptr = ptr->gChild;
			}	
			i++;
		}
		else if (word[i] == 'c') { 
			if (ptr->cChild) {
				ptr = ptr->cChild;
			}
			else {
				ptr->cChild = createNode(word[i]);
				ptr = ptr->cChild;
			}
			i++;
		}
		else if (word[i] == 't') {
			if (ptr->tChild) {
				ptr = ptr->tChild;
			}
			else {
				ptr->tChild = createNode(word[i]);
				ptr = ptr->tChild;
			}
			i++;
		}
	}
}

int Virus_gen::CountNodes(Virus_genome* Root) {
	if (Root == NULL) { return 0; } 
	else { 
		return (1 + CountNodes(Root->aChild) + CountNodes(Root->gChild) + CountNodes(Root->cChild) + CountNodes(Root->tChild));
	}
}


Virus_genome * Virus_gen::createNode(char val) {
	Virus_genome *newnode = new Virus_genome();
	newnode->Data = val; 
	newnode->aChild = NULL; 
	newnode->gChild = NULL;
	newnode->cChild = NULL;
	newnode->tChild = NULL;
	return newnode;
}

int Virus_gen::ComputeHeight(Virus_genome*  Root) {
	int aHeight, gHeight, cHeight, tHeight;
	if (Root == NULL) { return 0; } 

	aHeight = ComputeHeight(Root->aChild); 
	gHeight = ComputeHeight(Root->gChild);
	cHeight = ComputeHeight(Root->cChild); 
	tHeight = ComputeHeight(Root->tChild); 

	int left = fmax(aHeight, gHeight);
	int right = fmax(cHeight, tHeight);
	int largest = fmax(left, right);
	

	return largest + 1;
}

void Virus_gen::TraverseTree(char * sequence) {
	if (!Root) { return; }
	int size = 0;
	size = strlen(sequence);
	Virus_genome *ptr; 
	ptr = Root;
	int i = 1; 
	while (i < (size - 1)) {
		cout << "Checking   : " << sequence[i] << endl;
		if (sequence[i] == 'a') {
			if (ptr->aChild) { ptr = ptr->aChild; i++; }
			else { cout << endl << "broken: " << sequence[i] << "iteration " << i << endl; break; }
		}
		else if (sequence[i] == 'g') {
			if (ptr->gChild) { ptr = ptr->gChild; i++; }
			else { cout << endl << "Sequence: " << sequence[i] << "iteration " << i << endl; break; }
		}
		else if (sequence[i] == 'c') { 
			if (ptr->cChild) { ptr = ptr->cChild; i++; }
			else { cout << endl << "Sequence: " << sequence[i] << "iteration " << i << endl; break; }
		}
		else if (sequence[i] == 't') { 
			if (ptr->tChild) { ptr = ptr->tChild; i++; }
			else { cout << endl << "Sequence " << sequence[i] << "iteration " << i << endl; break; }
		}
		else {
			cout << "sequence broken." << endl; break;
		}
	}
	if (i == (size - 1)) 
	{
		cout << endl << "found: " << endl << endl << sequence << endl << endl;
	}
	else { cout << endl << "Unable find: " << endl << endl << sequence << endl << endl; }
}
