#include "Animal.h"

Animal::Animal(string paramFile, string mapFile){
    
    if (genomeInfoDone==false)
    {
        ANIMAL::G.num_breeds=1;

        ParmMap parameters;
        parameters.inputParms(paramFile);

        vector<double> seed_v      =   parameters["s"];
        int seed= int(seed_v[0]);
        ANIMAL::randGen.seed(seed);

        vector<double> numChr_v      =   parameters["n"];
        vector<double> nLoci_v       =   parameters["nL"];
        vector<double> chrLength_v   =   parameters["ch"];
        vector<double> mutRate_v     =   parameters["mu"];

        int numChr= int(numChr_v[0]);
        ANIMAL::G.set_num_chrom(numChr);
        double mutRate = mutRate_v[0];
        ANIMAL::G.mutRate=mutRate;
    
        int numLoci;
        double chrLength;
    
        int start=0;

        ifstream datafile;
        datafile.open(mapFile.c_str());
        if(!datafile) {
            cerr << "Couldn't open: " << mapFile << endl;
            exit (-1);
        }

        std::string inputStr;
        vector<string> tokens;
    
        for(auto &chromosome : ANIMAL::G){
        
            numLoci=nLoci_v[start];
            chrLength=chrLength_v[start];
            chromosome.set_num_loci(numLoci);
            chromosome.chr_length = chrLength;
            start++;
        
            getline(datafile,inputStr);
        
            boost::split(tokens, inputStr, boost::is_any_of(" "));
       
            unsigned i=0;
            for(auto &locus : chromosome){
                locus.locusType="M";
                locus.alleleFreq=0.0;
                locus.map_pos=getDouble(tokens[i]);
                i++;
            }
        }
    
        if(!ANIMAL::G.mapPosDone) ANIMAL::G.mkMapPos();
        makeGenomeInfoDone();
    
        datafile.clear();
        datafile.close();
    }
}
Animal::Animal(string paramFile, string mapFile){
    
    if (genomeInfoDone==false)
    {
        ANIMAL::G.num_breeds=1;

        ParmMap parameters;
        parameters.inputParms(paramFile);

        vector<double> seed_v      =   parameters["seed"];
        int seed= int(seed_v[0]);
        ANIMAL::randGen.seed(seed);

        vector<double> numChr_v      =   parameters["n"];
        vector<double> nLoci_v       =   parameters["nL"];
        vector<double> chrLength_v   =   parameters["ch"];
        vector<double> mutRate_v     =   parameters["mu"];

        int numChr= int(numChr_v[0]);
        ANIMAL::G.set_num_chrom(numChr);
        double mutRate = mutRate_v[0];
        ANIMAL::G.mutRate=mutRate;
    
        int numLoci;
        double chrLength;
    
        int start=0;

        ifstream datafile;
        datafile.open(mapFile.c_str());
        if(!datafile) {
            cerr << "Couldn't open " << mapFile << endl;
            exit (-1);
        }

        std::string inputStr;
        vector<string> tokens;
    
        for(auto &chromosome : ANIMAL::G){
        
            numLoci=nLoci_v[start];
            chrLength=chrLength_v[start];
            chromosome.set_num_loci(numLoci);
            chromosome.chr_length = chrLength;
            start++;
        
            getline(datafile,inputStr);
        
            boost::split(tokens, inputStr, boost::is_any_of(" "));
       
            unsigned i=0;
            for(auto &locus : chromosome){
                locus.locusType="Marker";
                locus.alleleFreq=0.0;
                locus.map_pos=getDouble(tokens[i]);
                i++;
            }
        }
    
        if(!ANIMAL::G.mapPosDone) ANIMAL::G.mkMapPos();
        makeGenomeInfoDone();
    
        datafile.clear();
        datafile.close();
    }
}
void ANIMAL::posigenomecls(ANIMAL& individual, vector<chromosome> &Genome)
{
    chromosome *Cell_current;
    unsigned numChromosomePair = G.get_num_chrom();
    

    ArrayXf tempPos;
    ArrayXi tempOri;
    tempPos.resize(ORIGIN_MAX);
    tempOri.resize(ORIGIN_MAX);
    
    for(unsigned i=0;i<numChromosomePair;i++)
    {
        double chrLength = G[i].chr_length;
        unsigned binomialN = chrLength*3 + 1;
        vector<float> rPos;
        binomial_distribution<int> Binom(binomialN,chrLength/binomialN);
        int numCrossover=Binom(randGen);
        
        uniform_real_distribution<float> u(0,1);
        for (unsigned k=0; k<numCrossover; k++)
        {
            rPos.push_back(chrLength*u(randGen));
        }
        rPos.push_back(chrLength);
        sort(rPos.begin(),rPos.end());
        
        unsigned   Cell_position=0;
        unsigned   startPosMe=0;
        
        Cell_current = (u(randGen)<0.5)?&individual.GenomePat[i]:&individual.GenomeMat[i];
        
        for(unsigned j=0;j<rPos.size();j++){
            
            unsigned numCopy=(Cell_current->pos < rPos[j]).count()-Cell_position;
            tempPos.block(startPosMe,0,numCopy,1)=Cell_current->pos.block(Cell_position,0,numCopy,1);
            tempOri.block(startPosMe,0,numCopy,1)=Cell_current->ori.block(Cell_position,0,numCopy,1);
            
            Cell_current=(Cell_current==&individual.GenomePat[i])?&individual.GenomeMat[i]:&individual.GenomePat[i];
            
            Cell_position=(Cell_current->pos < rPos[j]).count();
            startPosMe+=numCopy;
            tempPos(startPosMe)=rPos[j];
            tempOri(startPosMe)=Cell_current->ori(Cell_position-1);
            startPosMe++;
        } 
        
        unsigned keep=0;
        Genome[i].pos.resize(startPosMe);
        Genome[i].ori.resize(startPosMe);
        Genome[i].pos[0]=tempPos[0];
        Genome[i].ori[0]=tempOri[0];
        
        for(unsigned m=1;m < startPosMe-1; m++)
        {
            if(tempOri[m]!=tempOri[m-1])
            {
                keep=keep+1;
                Genome[i].pos[keep]=tempPos[m];
                Genome[i].ori[keep]=tempOri[m];
            }
        }
        unsigned PosOriSize=keep+1;
        Genome[i].pos.conservativeResize(PosOriSize);
        Genome[i].ori.conservativeResize(PosOriSize);
    }
}
void ANIMAL::animalcellpos()
{
    unsigned numChromosomePair = G.get_num_chrom();
    for(unsigned i=0;i<numChromosomePair;i++){
        GenomePat[i].ori.resize(1);
        GenomePat[i].ori[0]=countChromosome;
        GenomePat[i].pos.resize(1);
        GenomePat[i].pos[0]=0;
        
        GenomeMat[i].ori.resize(1);
        GenomeMat[i].ori[0]=countChromosome+1;
        GenomeMat[i].pos.resize(1);
        GenomeMat[i].pos[0]=0;
        
    }
    countChromosome += 2;
}
void ANIMAL::genomefund()
{
    sireId=0;
    damId=0;
    animalcellpos();
    initFounderHaps();
}

void ANIMAL::genomefund(vector<string> tokens1,vector<string> tokens2)
{
    sireId=0;
    damId=0;
    animalcellpos();
    fund(tokens1,GenomePat);
    fund(tokens2,GenomeMat);
}




void ANIMAL::initFounderHaps()
{
    unsigned numChromosomePair = G.get_num_chrom();
    
    for(unsigned i=0;i<numChromosomePair;i++)
    {
        unsigned numLoci=G[i].get_num_loci();
        GenomePat[i].haplotype.resize(numLoci);
        GenomeMat[i].haplotype.resize(numLoci);
        
        for(unsigned j=0;j<numLoci;j++)
        {
            binomial_distribution<int> Binom(1,G[i][j].alleleFreq);
            GenomePat[i].haplotype[j]= Binom(randGen);
            GenomeMat[i].haplotype[j]= Binom(randGen);
        }
    }
}


void ANIMAL::fund(vector<string> tokens,vector<chromosome> &Genome)
{
    unsigned numChromosomePair = G.get_num_chrom();
    
    for(unsigned i=0;i<numChromosomePair;i++)
    {
        unsigned numLoci=G[i].get_num_loci();
        Genome[i].haplotype.resize(numLoci);
        for(unsigned j=0;j<numLoci;j++){
            Genome[i].haplotype[j] = stod(tokens[j]);
            
            G[i][j].alleleFreq    += stod(tokens[j]);
        }
     }
}


void ANIMAL::sampleMyPosOri(ANIMAL& father, ANIMAL& mother)
{
    sireId  = father.myId;
    damId   = mother.myId;
    posigenomecls(father,GenomePat);
    posigenomecls(mother,GenomeMat);
}

void ANIMAL::sampleMyMutation()
{
    if(!G.mapPosDone) G.mkMapPos();
    sampleOneMutation(GenomePat);
    sampleOneMutation(GenomeMat);
}

void ANIMAL::sampleNonFounder(ANIMAL& father, ANIMAL& mother)
{
    sireId=father.myId;
    damId=mother.myId;
    sampleMyPosOri(father,mother);
    sampleMyMutation();
};

ArrayXf ANIMAL::getMyHapSeg(int i,int myOri,int start,int numcopy){
    
    if(myOri< 2*founders.size()){
        
        unsigned whichFounder=unsigned(myOri/2);
        chromosome* GenomePatOrMatInThisFounder= (myOri%2)? &ANIMAL::founders[whichFounder]->GenomeMat[i]: &ANIMAL::founders[whichFounder]->GenomePat[i];        
        
        return(GenomePatOrMatInThisFounder->haplotype.block(start,0,numcopy,1));
        
    }else{
        
        mutantInfo* thisMutant  =ANIMAL::mutants[myOri-2*founders.size()];
        auto myMutPos    =thisMutant->mutPos;
        auto myStartPos  =thisMutant->whichStartPos;
        auto myEndPos    =thisMutant->whichEndPos;
        auto myNewOri    =thisMutant->whichOri;
        
        unsigned myStartLocus   =(G[i].MapPos < myStartPos).count();
        unsigned myNEWnumcopy   =(G[i].MapPos < myEndPos).count()-myStartLocus;
        unsigned myMutationLoci =(G[i].MapPos <= myMutPos).count()-myStartLocus;  

        if(!(thisMutant->val).size()){
            
            
            thisMutant->val.resize(myNEWnumcopy);
            thisMutant->val=getMyHapSeg(i,myNewOri,myStartLocus,myNEWnumcopy); 
            
            unsigned myMutIndex=myMutationLoci-1;  
            
            thisMutant->val(myMutIndex)= 1-thisMutant->val(myMutIndex);
            
        }
        
        unsigned myStartKeep=start-myStartLocus;
        return((thisMutant->val).block(myStartKeep,0,numcopy,1));
        
        
    }
    
}

void ANIMAL::getMyHaps()
{
    getgenome(GenomePat);
    getgenome(GenomeMat);
}


void ANIMAL::getgenome(vector<chromosome> &Genome)
{
    unsigned numChromosomePair = G.get_num_chrom();
    
    for(unsigned i=0;i<numChromosomePair;i++)
    {
        
        unsigned numLoci=G[i].get_num_loci();
        Genome[i].haplotype.resize(numLoci);
        
        unsigned numOri = Genome[i].ori.size();
        
        for(unsigned segment=0;segment<numOri;segment++)
        {
            unsigned   numcopy;
            unsigned   start=(G[i].MapPos<Genome[i].pos[segment]).count();
            int        myOri=Genome[i].ori[segment];
            
            if(segment!= numOri-1)
            {
                numcopy=(G[i].MapPos<Genome[i].pos[segment+1]).count()-start;
            }else
            {
                numcopy=(G[i].MapPos<G[i].chr_length).count()-start;
            }
            
            Genome[i].haplotype.block(start,0,numcopy,1)
            =getMyHapSeg(i,myOri,start,numcopy); 
        }
    }
}


void ANIMAL::genometype()
{
    unsigned numChromosomePair = G.get_num_chrom();
    unsigned numTotalLoci=0;
    
    for(unsigned i=0;i<numChromosomePair;i++){
        numTotalLoci+=G[i].get_num_loci();
    }
    myGenotype.resize(numTotalLoci);
    
    unsigned start=0;
    for(unsigned i=0;i<numChromosomePair;i++){
        
        unsigned numLoci=G[i].get_num_loci();
        myGenotype.block(start,0,numLoci,1)= GenomePat[i].haplotype+GenomeMat[i].haplotype;
        start+=numLoci;
    }
    
    cout<<"Gi"<< this->myId <<endl;
}


void ANIMAL::display()
{
    IOFormat CommaInitFmt(StreamPrecision, DontAlignCols, ", ", ", ", "", "", " << ", ";");
    unsigned numChromosomePair = G.get_num_chrom();
    for(unsigned i=0;i<numChromosomePair;i++){
        cout<<"positions"<<i<<" for "<<endl;
        cout<< GenomePat[i].pos.format(CommaInitFmt)<<endl;
        cout<<"origins"<<i<<" for "<<endl;
        cout<< GenomePat[i].ori.format(CommaInitFmt)<<endl;
        cout<<"positions"<<i<<" for "<<endl;
        cout<< GenomeMat[i].pos.format(CommaInitFmt)<<endl;
        cout<<"origins"<<i<<" for "<<endl;
        cout<< GenomeMat[i].ori.format(CommaInitFmt)<<endl;
    }
};

unsigned ANIMAL::positionfinde()
{
    unsigned sumSizePat=0;
    unsigned sumSizeMat=0;
    for(unsigned i=0;i<G.get_num_chrom();i++){
        sumSizePat+=GenomePat[i].pos.size();
        sumSizeMat+=GenomeMat[i].pos.size();
    }
    return((sumSizePat+sumSizeMat)/2);
};
