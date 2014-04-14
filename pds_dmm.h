//
//  pds_dmm.h
//  pds_dmm
//
//  Created by Patrick Schloss on 11/7/12.
//  Copyright (c) 2012 University of Michigan. All rights reserved.
//

#ifndef pds_dmm_pds_dmm_h
#define pds_dmm_pds_dmm_h

/**************************************************************************************************/

#include <stdlib.h>     /* exit, EXIT_FAILURE */

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <map>
#include <algorithm>
#include <limits>
#include <string.h>

using namespace std;

/**************************************************************************************************/

template<typename T>
string toString(const T&x){
	
    stringstream output;
    output << x;
    return output.str();
	
}

/**************************************************************************************************/

inline void gobble(istream& f){
	
	char d;
    while(isspace(d=f.get()))		{;}
	f.putback(d);
	//
}

/**************************************************************************************************/

inline string getline(ifstream& fileHandle) {
    string line = "";
    
    while (fileHandle)	{
        //get next character
        char c = fileHandle.get(); 
        
        //are you at the end of the line
        if ((c == '\n') || (c == '\r') || (c == '\f') || (c == EOF)){  break;	}	
        else {		line += c;		}
    }
    
    return line;
}

/**************************************************************************************************/

inline double max(double a, double b){
    
    if(a>b) {    return a;  }
    else    {   return b;   }
}

/**************************************************************************************************/

inline void readSharedFile(string sharedFileName, vector<vector<int> >& sharedVector, vector<string>& otuNames, vector<string>& sampleNames){
    
    ifstream sharedFile(sharedFileName.c_str());
    
    string header = getline(sharedFile);
    string colHead;
    
    stringstream line(header);
    line >> colHead;
    line >> colHead;
    line >> colHead;
    
    while(line){
        line >> header;
        otuNames.push_back(header);
        gobble(line);
    }
    
    string label;
    string sample;
    int numOTUs;
    
    
    while(sharedFile){
        
        sharedFile >> label >> sample >> numOTUs;
        
        vector<int> counts(numOTUs, 0);
        for(int i=0;i<numOTUs;i++){
            sharedFile >> counts[i];
        }
        
        sharedVector.push_back(counts);
        sampleNames.push_back(sample);
        
        gobble(sharedFile);
    }

}

/**************************************************************************************************/

inline void readDesignFile(string designFileName, vector<string> sampleNames, vector<vector<double> >& partitions){
    
    int numSamples = (int)sampleNames.size();
    int numPartitions = 0;

    map<string, int> partitionIndex;
    
    string sample, partition;
    
    ifstream designFile(designFileName.c_str());
    for(int i=0;i<numSamples;i++){
                                  
        designFile >> sample >> partition;

        if(sample != sampleNames[i]){   cout << sample << "!=" << sampleNames[i] << endl;   }
        else{
            if(partitionIndex.count(partition) == 1){
                partitions[partitionIndex[partition]][i] = 1;
            }
            else{
                numPartitions++;
                partitions.resize(numPartitions);
                partitions[numPartitions-1].resize(numSamples);
                partitionIndex[partition] = numPartitions-1;
                
                partitions[partitionIndex[partition]][i] = 1;
                
            }
        }
    }
    designFile.close();
}

/**************************************************************************************************/

inline vector<double> generateDesignFile(int numPartitions, string fileRoot){

    vector<double> piValues(numPartitions, 0);
    
    ifstream postFile((fileRoot + toString(numPartitions) + "mix.posterior").c_str());
    ofstream designFile((fileRoot + "mix.design").c_str());

    vector<string> titles(numPartitions);
    
    for(int i=0;i<numPartitions;i++){   postFile >> titles[i];  }
    
    double posterior;
    string sampleName;
    int numSamples = 0;
    
    while(postFile){
        double maxPosterior = 0.0000;
        int maxPartition = -1;
        
        postFile >> sampleName;
        
        for(int i=0;i<numPartitions;i++){

            postFile >> posterior;
            if(posterior > maxPosterior){
                maxPosterior = posterior;
                maxPartition = i;
            }
            piValues[i] += posterior;
            
        }
        
        designFile << sampleName << '\t' << titles[maxPartition] << endl;
        
        numSamples++;
        gobble(postFile);
    }
    for(int i=0;i<numPartitions;i++){
        piValues[i] /= (double)numSamples;
    }
    
    
    postFile.close();
    designFile.close();
    
    return piValues;
}

/**************************************************************************************************/

struct summaryData {

    string name;
    double refMean, difference;
    vector<double> partMean, partLCI, partUCI;
    
};

/**************************************************************************************************/

inline bool summaryFunction(summaryData i, summaryData j){ return i.difference > j.difference;   }

/**************************************************************************************************/

inline void generateSummaryFile(int numPartitions, string fileRoot){
    
    vector<summaryData> summary;
    
    vector<double> pMean(numPartitions, 0);
    vector<double> pLCI(numPartitions, 0);
    vector<double> pUCI(numPartitions, 0);
    
    string name, header;
    double mean, lci, uci;
    
   
    vector<double> piValues = generateDesignFile(numPartitions, fileRoot);
    
    ifstream referenceFile((fileRoot + "1mix.relabund").c_str());
    ifstream partitionFile((fileRoot + toString(numPartitions) + "mix.relabund").c_str());

    header = getline(referenceFile);
    header = getline(partitionFile);
    stringstream head(header);
    string dummy, label;
    head >> dummy;
    vector<string> thetaValues(numPartitions, "");
    for(int i=0;i<numPartitions;i++){
        head >> label >> dummy >> dummy;
        thetaValues[i] = label.substr(label.find_last_of('_')+1);
    }
    
    
    vector<double> partitionDiff(numPartitions, 0.0000);

    while(referenceFile){
        referenceFile >> name >> mean >> lci >> uci;
        
        summaryData tempData;
        tempData.name = name;
        tempData.refMean = mean;
        
        double difference = 0.0000;
        
        partitionFile >> name;
        for(int j=0;j<numPartitions;j++){
            partitionFile >> pMean[j] >> pLCI[j] >> pUCI[j];
            difference += abs(mean - pMean[j]);
            partitionDiff[j] += abs(mean - pMean[j]);;
        }

        tempData.partMean = pMean;
        tempData.partLCI = pLCI;
        tempData.partUCI = pUCI;
        tempData.difference = difference;
        summary.push_back(tempData);
        
        gobble(referenceFile);
        gobble(partitionFile);
    }
    referenceFile.close();
    partitionFile.close();

    
    int numOTUs = (int)summary.size();
    
    sort(summary.begin(), summary.end(), summaryFunction);
    
    
    ofstream parameterFile((fileRoot + "mix.parameters").c_str());
    parameterFile.setf(ios::fixed, ios::floatfield);
    parameterFile.setf(ios::showpoint);

    double totalDifference =  0.0000;
    parameterFile << "Part\tDif2Ref_i\ttheta_i\tpi_i\n";
    for(int i=0;i<numPartitions;i++){
        parameterFile << i+1 << '\t' << setprecision(2) << partitionDiff[i] << '\t' << thetaValues[i] << '\t' << piValues[i] << endl;
        totalDifference += partitionDiff[i];
    }
    parameterFile.close();
    
    ofstream summaryFile((fileRoot + "mix.summary").c_str());
    summaryFile.setf(ios::fixed, ios::floatfield);
    summaryFile.setf(ios::showpoint);
    
    
    summaryFile << "OTU\tP0.mean";
    for(int i=0;i<numPartitions;i++){
        summaryFile << "\tP" << i+1 << ".mean\tP" << i+1 << ".lci\tP" << i+1 << ".uci";
    }
    summaryFile << "\tDifference\tCumFraction" << endl;
    
    double cumDiff = 0.0000;
    
    for(int i=0;i<numOTUs;i++){
        summaryFile << summary[i].name << setprecision(2) << '\t' << summary[i].refMean;
        for(int j=0;j<numPartitions;j++){
            summaryFile  << '\t' << summary[i].partMean[j] << '\t' << summary[i].partLCI[j] << '\t' << summary[i].partUCI[j];
        }
        
        cumDiff += summary[i].difference/totalDifference;
        summaryFile << '\t' << summary[i].difference << '\t' << cumDiff << endl;
    }
    summaryFile.close();

}

/**************************************************************************************************/

#endif
