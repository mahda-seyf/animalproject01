#include <iostream>
#include <string>
#include <algorithm>
#include "Animal.h"

using namespace std;

struct mutantInfo
{
    unsigned whichOri;
    float    whichStartPos;
    float    whichEndPos;
    float    mutPos;
    ArrayXf  val;
};


struct chrs
{
    ArrayXf    typehp;
    ArrayXi    ori;
    ArrayXf    pos;
};

class ANimal
{    
public:
    ANimal(void)
    {
        number = ++countId;
        unsigned numchrsPair=G.get_num_chrom();
        GenomePat.resize(numchrsPair);
        GenomeMat.resize(numchrsPair);
    };
    
    vector<chrs> GenomePat;
    vector<chrs> GenomeMat;
    
    int number, sireId, damId;
    
    static vector<ANimal*> founders;
    static vector<mutantInfo*> mutants;
 
    void smpf();
    void smpf(vector<string> tokens1,vector<string> tokens2);
    void initF ();
    void rHaps();                   
    void inerHaps(vector<string> tokens,vector<chrs> &Genome);  
    void inptHaps(vector<string> tokens);
    void inaps(vector<string> tokens);
 
    void sampleNonFounder(ANimal& father, ANimal& mother);
    void sampleMyPosOri(ANimal& father, ANimal& mother);
    void sampleOnePosOri(ANimal& individual,vector<chrs> &Genome);

    void sampleMyMutation();
    void sampleOneMutation(vector<chrs> &Genome);

    void getMyHaps();
    void getOneHaps(vector<chrs> &Genome);
    ArrayXf getMyHapSeg(int i,int myOri,int start,int numcopy);
    
    void display();
    unsigned displayNumPos();
    ArrayXf myGenotype;
    void getMyGenotype();
 
    static Genome_info G;
    static unsigned countchrs;
    static unsigned countId;
    static default_random_engine randGen;
};
#endif  