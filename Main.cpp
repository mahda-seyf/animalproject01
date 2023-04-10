#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <tchar.h>
#include <stdlib.h>
#include <time.h>
#include <vector>

using namespace std;


int Main(int argc, const char * argv[])
{
	int vsize = 0;
	int len = 0;
	char input;

	cout << "sequence" << endl;
	while (len < 1 || len > 1000000000) {
		cout << "genome " << endl;
		cin >> len;
		cin.clear();
		cin.ignore(1000000000, '\n');
	}

	cout << "sequence" << endl;
	while (vsize < 1 || vsize > 10000) {
		cout << "length" << endl;
		cin >> vsize;
		cin.clear();
		cin.ignore(1000000000, '\n');
	}

	cout << "could not" << len << " characters" << endl;

	Genome_inp a, g, c, t;

	ifstream in;
	in.open("genome.txt"); 
	if (!in.is_open()) {
		cout << "could not \t";
	}
	
	char * word = new char[vsize]; 
	
	word[vsize] = '\0'; 


	for (int i = 1; i < (len - vsize -1); i++)  
	{
		for (int j = 0; j < (vsize - 1); j++) 
			word[j] = word[j + 1];
		in >> word[(vsize - 1)];
		word[vsize] = '\0';
		if (i % (vsize - 1) == 0) 
		{
			if (word[0] == 'a')a.InsertNode(word); 
			else if (word[0] == 'g')g.InsertNode(word);
			else if (word[0] == 'c')c.InsertNode(word);
			else if (word[0] == 't')t.InsertNode(word);
		}
	}
	
	in.clear();
	in.close(); 
	
	word = NULL;
	delete word; 
	 
	cout << "a" << a.Cell_len(a.Root) << "a" << a.CountNodes(a.Root) << endl;	
	cout << "g" << g.Cell_len(g.Root) << "g" << g.CountNodes(g.Root) << endl;
	cout << "c" << c.Cell_len(c.Root) << "c" << c.CountNodes(c.Root) << endl;
	cout << "t" << t.Cell_len(t.Root) << "t" << t.CountNodes(t.Root) << endl;


	vector<string> rds; 
	get_reads("virus.txt", rds, vsize); 
	char *wo = new char[(vsize + 1)]; 
	for (int k = 0; k < rds.size(); k++) 
	{
		cout << wo << endl;
		cout << wo[0];
		cout << "sequence: " << rds[k] << endl;
		string rd = rds[k]; 
		strcpy_s(wo, (vsize + 1), rd.c_str());
		if (wo[0] == 'a')a.TraverseTree(wo);
		else if (wo[0] == 'g') g.TraverseTree(wo); 
		else if (wo[0] == 'c') c.TraverseTree(wo);
		else if (wo[0] == 't')t.TraverseTree(wo);
		else break;
	}

	cout << "sequence" << endl;

	char * test_sequence = new char[vsize]; 
											
	string test = "";
	cout << "length " << vsize << "Enter y or n." << endl;
	cin >> input;
	cin.ignore(100, '\n');
	Cell_check(input); 
	
	while (true) {
		if (input == 'y' || input == 'Y') {
			cout << "Note" << vsize << "accepted." << endl;
			cout << "sequence: " << endl;
			cin >> test; 
			
			while (test.length() != vsize)
			{
				cout << "size : " << vsize << " characters" << endl;
				cin >> test; 
			}

			for (int k = 0; k < vsize; k++) {
				test[k] = tolower(test[k]);
			}

			strcpy_s(test_sequence, (vsize + 1), test.c_str());
			
			if (test_sequence[0] == 'a') { (a.TraverseTree(test_sequence)); } 
			else if (test_sequence[0] == 'g') { g.TraverseTree(test_sequence); } 
			else if (test_sequence[0] == 'c') { c.TraverseTree(test_sequence); }
			else if (test_sequence[0] == 't') { t.TraverseTree(test_sequence); }
		}
		cout << " length " << vsize << " Enter y or n." << endl;
		cin >> input; 
		cin.ignore(100, '\n');
		Cell_check(input);
		if (input == 'n' || input == 'N') { break; }
	}

	delete wo;
	wo = NULL; 
	
	system("pause");
	return true;

    unsigned popSize =  10;
    unsigned nGen    =  3;
    unsigned founderSize=10;
  
  
    Animal osim(genomeFile,mapFile);
    osim.popFounders(founderSize,haplotype);
   
    Animal osim1(genomeFile,mapFile);
    Animal osim2(genomeFile,mapFile);
    Animal osim3(genomeFile,mapFile);

    osim.popSample(popSize,nGen);
    
    Cell_Matrix outGeno;
    outGeno=osim.getGenotypes();
    outFile << outGeno;
    cout<<"genotyped."<<endl;
    cout<<endl<<endl;
    
    
    Cell_Vector outPhen=osim.getPhenotype(10, 0.5);
    cout<<outPhen<<endl;

    cout<<endl<<endl<<endl<<endl;
    cout<< "\n";
	