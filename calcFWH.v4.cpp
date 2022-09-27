#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <map>
#include <algorithm>
#include <cmath>

int sampleCount;
int markerCount;


int calcSegSiteD(std::vector<std::string>& seqVec, int& ancAlleleSeqPos){
    int segSiteD;
    segSiteD = 99999999;
    std::map<std::string,int> alleleCountMap;
    std::string alleleA;
    for(int z=0; z<seqVec.size();z++){
        if (seqVec[z] != "N" && z != ancAlleleSeqPos){
            if( alleleCountMap.find(seqVec[z]) == alleleCountMap.end()){
                alleleCountMap[seqVec[z]] =1;
            }
            else{
                alleleCountMap[seqVec[z]] +=1;
            }
        }
    }
    if( ancAlleleSeqPos >= 0 ){
        alleleA = seqVec[ancAlleleSeqPos];
        if( alleleA == "N"){alleleA = "NA";}
    }
    if( ancAlleleSeqPos < 0){ alleleA = "N";}
    if( alleleA != "NA"){
        if( alleleCountMap.size() >2 ){std::cout<<"INPUT ERROR: the number of alleles greater than two "<<"\n";std::exit(EXIT_FAILURE);}
        else{
                std::map<std::string,int>::iterator b;
                for( b= alleleCountMap.begin(); b!= alleleCountMap.end(); b++ ){
                    if( alleleA =="N" ){
                        if( b->second < segSiteD ){
                            segSiteD = b-> second;
                        }
                    }
                    else{
                        if (b-> first != alleleA){
                            segSiteD = b->second;
                        }
                    }
                }
        }
    }
    return segSiteD;

}

std::map <int,std::vector<int>> mapFToMap(const std::string& mapFile){

    std::map<int,std::vector<int>> chromPosMap;
        std::ifstream source;
        source.open(mapFile);
        std::string line;
        while (std::getline(source, line)){
            std::stringstream ss ( line );
            std::string word;
            std::vector<std::string> vectorLine;
            while ( std::getline ( ss, word, ' ' ) ){
                vectorLine.push_back(word);
            }
    
    if(vectorLine.size() != 2){std::cout<<"INPUT ERROR: the total numbers of column is not equal to 2 "<<line<<"\n";std::exit(EXIT_FAILURE);}
        if(chromPosMap.size() == 0){
        std::vector<int>posVec;
        posVec.push_back(std::stoi(vectorLine[1]));
        chromPosMap[std::stoi(vectorLine[0])] = posVec;
    }
    else{
        chromPosMap[std::stoi(vectorLine[0])].push_back(std::stoi(vectorLine[1]));
    }
    }


    return chromPosMap;
}

std::vector<std::string> phylipToVect(const std::string& phylipFile){

        std::vector<std::string>seqVec;
        std::ifstream source;
        source.open(phylipFile);
        std::string line;
    int lineCount = 0;
        while (std::getline(source, line)){
        lineCount++;
            std::stringstream ss ( line );
            std::string word;
            std::vector<std::string> outputLine;
            std::vector<std::string> vectorLine;
            while ( std::getline ( ss, word, ' ' ) ){
                    vectorLine.push_back(word);
                }
        if (lineCount == 1){
            sampleCount = std::stoi(vectorLine[0]);
            markerCount = std::stoi(vectorLine[1]);
        }
        else{
        if(vectorLine.size() != 2){std::cout<<"INPUT ERROR: the total numbers of column is not equal to 2 "<<vectorLine[0]<<"\n";std::exit(EXIT_FAILURE);}
            seqVec.push_back(vectorLine[1]);
        }

        }
    return seqVec;
}

std::map<int,int> segToSfs(std::vector<int>& segCountVec){
    
    std::map<int,int>sfsArray;
    for(int i=1;i<sampleCount;i++){
        int freq = std::count(segCountVec.begin(),segCountVec.end(),i);
        sfsArray[i] = freq;
    }
    return sfsArray;

}


std::vector<float> calcPopMeasures(std::map<int,int>& sfsArray, int& ancAlleleSeqPos){

    int popSampleCount;
    float a1 = 0.0;
    float a2 = 0.0;
    float thetaL_num = 0.0;
    int numPi = 0 ;
    int numH = 0;
    int segSites = 0;
    popSampleCount = sampleCount;
    if( ancAlleleSeqPos >= 0){
        popSampleCount -= 1;
    }
    int denomPi = popSampleCount * ( popSampleCount - 1 );
    for( int i=0; i<(sampleCount -1) ; i++ ){
        int di = 0;
	a1 += 1/ (float) (i+1);
	a2 += 1/((float) (i+1)*(float) (i+1));
        if ( sfsArray.find(i) != sfsArray.end() ){
            di = sfsArray[i];
            segSites += di;
        }
        numPi += 2 * di * i * ( popSampleCount - i);
        numH += 2*di*i*i;
	thetaL_num += i*di;
    }

    float b1 = (float)(popSampleCount+1)/(float)(3*popSampleCount-3);
    float b2 = 2*((float)(popSampleCount*popSampleCount)+(float)(popSampleCount)+3)/(9*(float)(popSampleCount)*(float)(popSampleCount-1));
    float c1 = b1 - (1/a1);
    float c2 = b2 - ((float) (popSampleCount+2)/( a1 * (float) popSampleCount))+(a2/(a1*a1));
    float e1 = c1/a1;
    float e2 = c2/(a1*a1+a2);
    float numSeg = (float) segSites;
    float thetaW = (float) numSeg/ a1;
    float pi = (float) numPi/(float) denomPi;
    float tajimaDNum = pi - thetaW;
    float tajimaDDeno = e1*numSeg + e2*numSeg*(numSeg-1);
    float tajimaD = tajimaDNum/sqrt(tajimaDDeno);

    float thetaL = thetaL_num/(float) (popSampleCount-1);
    float thetaH = (float) numH/(float) denomPi;
    float FWH = pi - thetaH;
    float Hnum = pi - thetaL;
    float Hden_LHS = ((popSampleCount-2)/6*(popSampleCount-1))*thetaW;
    float theta_square = (numSeg * (numSeg-1))/(a1*a1+a2);
    float Hden_RHSnum = ((18*popSampleCount*popSampleCount*(3*popSampleCount+2)*a2*a2)-
		    (88*popSampleCount*popSampleCount*popSampleCount+9*popSampleCount*popSampleCount-13*popSampleCount+
		     6)) * theta_square;
    int Hden_RHSden = 9*popSampleCount*(popSampleCount-1)*(popSampleCount-1);
    float Hden_RHS = Hden_RHSnum/(float) Hden_RHSden;
    float Hden = sqrt(Hden_RHS+Hden_LHS);
    float H = Hnum/Hden;  

    std::vector<float> popMeasures{numSeg, thetaW, pi, tajimaD, FWH, H};
    
    return popMeasures;


}

void printPopMeasures(std::vector<float>& popMeasures, int& startWindow){
    std::cout<<startWindow<<"\t";
    for( auto & measure: popMeasures ){
        std::cout<<measure<<"\t";
        }
    std::cout<<"\n";

}


void vectToArray(std::vector<std::string>& seqVec, std::map <int,std::vector<int>>& mapFMap, int& windowSize, int& ancAlleleSeqPos){

    std::vector<int> segCountVec;
    std::vector<std::string> tmpSeqVec;
    int startWindow;
    startWindow = windowSize;
    std::map<int,std::vector<int>>::iterator i;
    for(i= mapFMap.begin(); i!= mapFMap.end(); i++){
        for(int p=0; p< i->second.size(); p++){
            for(int s=0; s<seqVec.size();s++){
		std::string seq;
		std::string base;
		seq = seqVec[s];
		base = seq[p];
                tmpSeqVec.push_back(base);
            }
            int segSite = calcSegSiteD(tmpSeqVec, ancAlleleSeqPos);
            tmpSeqVec.clear();
            if(i->second[p] <= startWindow){
                segCountVec.push_back(segSite);
            }
            else{
                std::map<int,int> sfsArray = segToSfs(segCountVec);
                std::vector<float> popMeasures = calcPopMeasures(sfsArray, ancAlleleSeqPos);
                segCountVec.clear();
                segCountVec.push_back(segSite);
                printPopMeasures(popMeasures, windowSize);
                while(i->second[p] > startWindow){
                startWindow += windowSize;
                }
            
            }
        
        }
        if(segCountVec.size()>0){
        
                std::map<int,int> sfsArray = segToSfs(segCountVec);
                std::vector<float> popMeasures = calcPopMeasures(sfsArray, ancAlleleSeqPos);
                segCountVec.clear();
                printPopMeasures(popMeasures, startWindow);
        
        }
        }
}



int main(int argc, char **argv){
    std::string phylipFile(argv[1]);
    std::string mapFile(argv[2]);
    int windowSize = std::stoi(argv[3]);
    int ancAlleleSeqPos = std::stoi(argv[4]);
    std::map <int,std::vector<int>> mapFMap = mapFToMap(mapFile); //this is the output associative array of this function {1:(1,2,3),2:(3,4,5)}
    std::vector<std::string>seqVec = phylipToVect(phylipFile); //this is the output vector of this function (("A","T","G"),("A","C","T")), the length should be equal to the sample size
    vectToArray(seqVec, mapFMap, windowSize, ancAlleleSeqPos);
}

