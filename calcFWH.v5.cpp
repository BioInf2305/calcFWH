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


std::map<std::string,std::vector<int>> indToMap(const std::string& indFile){

    	std::map<std::string,std::vector<int>> popMap;
	//the following variable link the individual from .indi file to its respective genotype position in each row in geno file
	int lineNumber = 0;

        std::ifstream source;
        source.open(indFile);
        std::string line;
	//read the file
        while (std::getline(source, line)){
	    lineNumber +=1 ;
            std::stringstream ss ( line );
            std::string word;
            std::vector<std::string> vectorLine;
	    //store its content in another vector, easy to retrieve the desire column
            while ( std::getline ( ss, word, ' ' ) ){
                vectorLine.push_back(word);
            }
	std::string popName = vectorLine[2];

        if(popMap.find(popName) == popMap.end()){
		//if the population id is not in popmap file,create empty vector as value
		std::vector<int>lineCountVec;
		popMap[popName] = lineCountVec;
	 	}
	popMap[popName].push_back(lineNumber);
    }


    return popMap;
}

std::vector<int> snpToVec(const std::string& snpFile, int& windowSize){
	
	std::vector<int> snpVec;
	//the following variable link the marker record from .snp file to its respective row in geno file
	int lineNumber = -1;

	int startWindow = windowSize;
	std::string chrom = "NA";
	std::ifstream source;
	source.open(snpFile);
	std::string line;
	while(std::getline(source,line)){
		lineNumber +=1 ;
		std::stringstream ss (line);
		std::string word;
		std::vector<std::string>vectorLine;
		while(std::getline(ss, word, ' ')){
			vectorLine.push_back(word);
		}
		//if the chromosome name changes, add the last record as the row number
		if (vectorLine[0] != chrom){
			chrom = vectorLine[0];
			startWindow = windowSize;
			snpVec.push_back(lineNumber);
		}
		int pos = std::stoi(vectorLine[3]);
		if(pos>startWindow){
		int insideWhile = 0;
		while(pos > startWindow){
			insideWhile = 1;
			startWindow += windowSize;
			}	
		if(insideWhile == 1){
			snpVec.push_back(lineNumber);
			}

		}

			
	}
	
	return snpVec;

}

std::vector<int> countGenotypes(std::string& popString){
	
	std::vector<int>countGenoVec;

	std::string::difference_type refCount = std::count(popString.begin(), popString.end(), '2');
    	std::string::difference_type hetCount = std::count(popString.begin(), popString.end(), '1');
    	std::string::difference_type altCount = std::count(popString.begin(), popString.end(), '0');

	int refAlleleCount = refCount * 2 + hetCount;
	int altAlleleCount = altCount * 2 + hetCount;

	countGenoVec.push_back(refAlleleCount);
	countGenoVec.push_back(altAlleleCount);

	return countGenoVec;

}

std::map<std::string,int> countDerivedAllele(std::string& line, std::map<std::string,std::vector<int>& popMap, const std::string& outPopName){

	std::map<std::string,int> locusDerivedCountMap;

	if(outPopName != "NA"){
		std::vector<int>posVec;
		std::string outPopGeno="";
		posVec = popMap[outPopName];
		for(int i : line){
			outPopGeno += line[i]
		}
		std::vector<int>countGenoVec = countGenotypes(outPopGeno);
		std::int derivedAllele = (countGenoVec[1] >= countGenoVec[0]) ? 1 : 0;

	}
	else{
		std::vector<int>countGenoVec = countGenotypes(line);
		std::int derivedAllele = (countGenoVec[1] >= countGenoVec[0]) ? 1 : 0;
	}

	std::map<std::string,std::vector<int>:: iterator p;
	for(p = popMap.begin(); p!=popMap.end(); p++){
		std::string popSeq = "";
		if(p->first != outPopName){
			for(int i=0; i<p->second.size();i++){
				popSeq += line[p->second[i]];
			}
			std::vector<int>countGenoVec = countGenotypes(popSeq);
			locusDerivedCountMap[p->first] = countGenoVec[derivedAllele];
		}
		
	}

	return locusDerivedCountMap;
}

void writeOutput(const std::string& genoFile, std::map<std::string,std::vector<int>>& popMap, std::vector<int>& snpVec, const std::string& outPopName){
	
	FILE *fp[popMap.size()];
	//the following variable compare the element from snpvec with that of the line count in the geno file
	int snpVecCompare = 0;

	int lineCount = 0;
	//following map store pop as key and the vector of derived alleles as its value
	std::map<std::string,std::vector<int>> popDerivedCountMap;

	std::vec<std::string>popVec;

	std::map<std::string,std::vector<int>>::iterator p;
	for(p=popMap.begin(); p!=popMap.end();p++){
		std::vector<int>derivedAlleleCountVec;
		popDerivedCountMap[p->first] = derivedAlleleCountVec;
		popVec.push_back(p->first);
	}
	
	//open file pointer to write and append the output for each population
	
	for(int i=0;i<popVec.size();i++){
		fp[i] = (popVec[i]+".out.txt","a");
	}

	std::ifstream source;
	source.open(genoFile);
	std::string line;
	while(std::getline(source,line)){
		lineCount +=1;
		std::map<std::string,int> locusDerivedCountMap = countDerivedAllele(line, popMap, outPopName)
		


}


int main(int argc, char **argv){
    std::string genoFile(argv[1]);
    std::string snpFile(argv[2]);
    std::string indFile(argv[3]);
    int windowSize = std::stoi(argv[4]);
    std::string outPopName(argv[5]);
    //the following function stores the pop and its sample position in each row of the geno file 
    std::map<std::string,std::vector<int>> popMap = indToMap(indFile);
    std::vector<int>snpVec = snpToVec(snpFile, windowSize);
    writeOutput(genoFile,popMap,snpTovec, outPopName);
}
