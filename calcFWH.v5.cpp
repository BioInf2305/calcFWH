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

//following function output the map with pop name as a key and a vector of its indices (position in geno line) as value
//geno file 
//11112222
//so first four genotype belongs to pop1 and then the next genotypes belong to pop2
//output map--> {pop1,{0,1,2,3},pop2,{4,5,6,7}}

std::map<std::string,std::vector<int>> indToMap(const std::string& indFile){

    	std::map<std::string,std::vector<int>> popMap;
   
   //the following variable link the individual from .indi file to its respective genotype position in each row in geno file
	int lineNumber = -1;

        std::ifstream source;
        source.open(indFile);
        std::string line;

	    //read .indi file
        while (std::getline(source, line)){

                lineNumber +=1 ;
                std::stringstream ss ( line );
                std::vector<std::string> vectorLine;

            //store its content in another vector, easy to retrieve the desire column
                for ( std::string line; ss >> line; ){

                    vectorLine.push_back(line);
                }
                std::string popName = vectorLine[2];
                if(popMap.find(popName) == popMap.end()){

            //if the population id is not in popmap file,create empty vector as value
                    std::vector<int> lineCountVec;

                    lineCountVec.push_back(lineNumber);

	 	        }
	            popMap[popName].push_back(lineNumber);
		
	
    }


    return popMap;
}

//folllowing function output the map with line number as key and vector of chrom name, and start position and end position of the windows as values
//{1,{chrm1,0,50},2,{chrm1,0,50}}

std::map<int,std::vector<std::string>> snpToMap(const std::string& snpFile, int& windowSize){
	
	std::map<int,std::vector<std::string>> snpLineMap;
	
    //the following variable link the marker record from .snp file to its respective row in geno file, it is also the key in map
	int lineNumber = 0;

	int startWindow = 0;

	int pos = 0;
	std::string minWindow = "-9";
	std::string maxWindow = "-9";
	std::string chrom = "NA";
    
    //folowing vector stores the splitted content of each line
	std::vector<std::string>vectorLine;

	std::ifstream source;
	source.open(snpFile);
	std::string line;

    //start reading the file
	while(std::getline(source,line)){

		lineNumber ++ ;
		std::stringstream ss (line);
        
        //clear the content of the previous line
		vectorLine.clear();

        //split line by white space delimiter
		for (std::string line; ss >> line;){
			vectorLine.push_back(line);
		}
		
        //if the chromosome name changes, add the last record in the map
		if (vectorLine[1] != chrom && chrom != "NA" ){
			chrom = vectorLine[1];
			startWindow = windowSize;
			std::vector<std::string>tmpVec {vectorLine[1], minWindow, maxWindow};

            //substract -1 as +1 has been added to line number
			snpLineMap[ lineNumber - 1 ] = tmpVec;
		}
		
		int pos = std::stoi(vectorLine[3]);
		chrom = vectorLine[1];

        // if the position increases the startwindow, add the lineNumber to the map
		if(pos>startWindow){

			int insideWhile = 0;

            // here it is assumed that there might be gap of more than the windowSize between two consecutive SNPs, therefore keep adding windowSize
			while(pos > startWindow){

				insideWhile = 1;
				startWindow += windowSize;
				}	
			if(insideWhile == 1){

            //minWindow is the closest position from the pos with remainder zero 
				minWindow = std::to_string(pos - pos % windowSize);

            //maxWindow is just the minWindow + windowSize 
				maxWindow = std::to_string(std::stoi(minWindow)+windowSize);

				std::vector<std::string>tmpVec {vectorLine[1], minWindow, maxWindow};
				snpLineMap[lineNumber] = tmpVec;
				}

		}

			
	}

	//add the last record in map if it exceeds the maxWindow size variable stored in the previously added record
	if( pos > stoi(maxWindow) ){

		std::vector<std::string>tmpVec {vectorLine[1],minWindow, maxWindow};
		snpLineMap[lineNumber] = tmpVec;
	}
	
	return snpLineMap;

}

//following function sorts the snp line map
//{line number, {chrm, pos}}
//input --> {3,{chrm2,10},4,{chrm2,15},2,{chrm1,20},1,{chrm1,10}}
//output --> {1,{chrm1,10},2,{chrm1,20},3,{chrm2,10},4,{chrm2,15}}

std::vector<int> snpLineMapToSortedLineNumVec(std::map<int,std::vector<std::string>>& snpLineMap){

	std::vector<int>sortedLineNumVec;

	std::map<int,std::vector<std::string>>:: iterator i;
	for(i = snpLineMap.begin(); i != snpLineMap.end(); i++){
		sortedLineNumVec.push_back(i->first);
	}
	std::sort(sortedLineNumVec.begin(), sortedLineNumVec.end());
	return sortedLineNumVec;

}

//following function count the genotype 

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

std::map<std::string,int> countDerivedAllele(std::string& line, std::map<std::string,std::vector<int>>& popMap, const std::string& outPopName){

	std::map<std::string,int> locusDerivedCountMap;
	int derivedAllele;

	if(outPopName != "NA"){
		std::vector<int>posVec;
		std::string outPopGeno="";
		posVec = popMap[outPopName];
		for(int i : posVec){
			outPopGeno += line[i];
		}
		std::vector<int>countGenoVec = countGenotypes(outPopGeno);
		derivedAllele = (countGenoVec[1] >= countGenoVec[0]) ? 1 : 0;

	}
	else{

		std::vector<int>countGenoVec = countGenotypes(line);
		derivedAllele = (countGenoVec[1] >= countGenoVec[0]) ? 0 : 1;
	}
	std::map<std::string,std::vector<int>>:: iterator p;
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


std::map<int,int> alleleCntToSfs(std::vector<int>& derivedAlleleCountVec, int& sampleCountPop){

	std::map<int,int>popSfs;

	std::sort(derivedAlleleCountVec.begin(), derivedAlleleCountVec.end());

	for(int i=0; i<sampleCountPop+1;i++){
		int siteCount;
		siteCount = std::count(derivedAlleleCountVec.begin(), derivedAlleleCountVec.end(), i);
		popSfs[i] = siteCount;
	}
	
	return popSfs;

}


std::vector<float> calcPopMeasures(std::map<int,int>& sfsArray, int& samplePopCount){

    float a1 = 0.0;
    float a2 = 0.0;
    float thetaL_num = 0.0;
    int numPi = 0 ;
    int numH = 0;
    int segSites = 0;
    int popSampleCount = 2 * samplePopCount;
    int denomPi = popSampleCount * ( popSampleCount - 1 );
    for( int i=1; i<popSampleCount ; i++ ){
	a1 += 1/ (float) i;
	a2 += 1/((float) i * (float) i);
        segSites += sfsArray[i];
        numPi += 2 * sfsArray[i] * i * ( popSampleCount - i);
        numH += 2* sfsArray[i] * i *i;
	thetaL_num += i * sfsArray[i];
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


void writeOutput(std::map<std::string,std::vector<int>>& popDerivedCountMap, std::map<std::string, std::vector<int>>& popMap, std::vector<std::string>& popVec, std::map<int,std::vector<std::string>>& snpLineMap, int& lineCount){

	//create file pointer variable

	std::ofstream os;

	for(std::string pop : popVec ){
		int sampleSize = popMap[pop].size();
		std::map<int,int>popSfs = alleleCntToSfs(popDerivedCountMap[pop], sampleSize);
		std::vector<float>popMeasures = calcPopMeasures(popSfs, sampleSize);
		std::string popName = pop;
		const char *popF = popName.c_str();
		std::vector<std::string>tmpMap = snpLineMap[lineCount];
		os.open(popF, std::ofstream::out | std::ofstream::app);
		os<<tmpMap[0]<<"\t"<<tmpMap[1]<<"\t"<<tmpMap[2]<<"\t";
		for(float measure: popMeasures ){
			os<< measure <<"\t";
		}
		os<<"\n";
		os.close();
	}

}

void readGenoWriteOutput(const std::string& genoFile, std::map<std::string,std::vector<int>>& popMap, std::map<int,std::vector<std::string>>& snpLineMap, const std::string& outPopName){
	

	
	//following method retrieve the sorted vector of line number
	
	std::vector<int>sortedLineNumVec = snpLineMapToSortedLineNumVec(snpLineMap);


	//the following variable compare the element from snpvec with that of the line count in the geno file
	int snpVecCompare = 0;

	int lineCount = 0;
	
	//following map store pop as key and the vector of derived alleles as its value
	
	std::map<std::string,std::vector<int>> popDerivedCountMap;

	std::vector<std::string>popVec;

	std::map<std::string,std::vector<int>>::iterator p;
	for(p=popMap.begin(); p!=popMap.end();p++){
		std::vector<int>derivedAlleleCountVec;
		popDerivedCountMap[p->first] = derivedAlleleCountVec;
		popVec.push_back(p->first);
	}
	
	std::ifstream source;
	source.open(genoFile);
	std::string line;
	while(std::getline(source,line)){
		lineCount +=1;
		std::map<std::string,int> locusDerivedCountMap = countDerivedAllele(line, popMap, outPopName);
		if ( snpVecCompare+1 != sortedLineNumVec.size() ) {
			if( lineCount >= sortedLineNumVec[snpVecCompare+1] ){
				writeOutput(popDerivedCountMap, popMap, popVec, snpLineMap, sortedLineNumVec[snpVecCompare]);
				snpVecCompare ++;
				std::map<std::string,std::vector<int>>::iterator p;
				for(p=popMap.begin(); p!=popMap.end();p++){
					popDerivedCountMap[p->first].clear();
				}
			
			}
		}
		for(std::string i : popVec){
			popDerivedCountMap[i].push_back(locusDerivedCountMap[i]);
		}
	
	}
	if (popDerivedCountMap[popVec[0]].size() > 0){
		writeOutput(popDerivedCountMap, popMap, popVec, snpLineMap, sortedLineNumVec[snpVecCompare]);
	}
}



int main(int argc, char **argv){
    std::string genoFile(argv[1]);
    std::string snpFile(argv[2]);
    std::string indFile(argv[3]);
    int windowSize = std::stoi(argv[4]);
    std::string outPopName(argv[5]);
    //the following function stores the pop and its sample position in each row of the geno file 
    std::map<std::string,std::vector<int>> popMap = indToMap(indFile);
    std::map<int,std::vector<std::string>>snpLineMap = snpToMap(snpFile, windowSize);
    readGenoWriteOutput(genoFile,popMap,snpLineMap, outPopName);
}
