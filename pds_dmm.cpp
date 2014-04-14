//
//  pds_dmm.cpp
//  pds_dmm
//
//  Created by Patrick Schloss on 11/7/12.
//  Copyright (c) 2012 University of Michigan. All rights reserved.
//

/**************************************************************************************************/

#include "pds_dmm.h"
#include "qFinderDMM.h"

/**************************************************************************************************/


int main(int argc, char *argv[]){
    
    srand( (unsigned)time( NULL ) );
    
    cout.setf(ios::fixed, ios::floatfield);
    cout.setf(ios::showpoint);
    
    string sharedFileName, designFileName;
    vector<vector<int> > sharedMatrix;
    vector<string> otuNames;
    vector<string> sampleNames;

    int minNumPartitions = 5;
    int maxNumPartitions = 100;
    int optimizeGap = 3;
    
    if(argc > 1) {
        for(char **p=argv+1;p<argv+argc;p++) {
            if(strcmp(*p,"-shared")==0) {
                if(++p>=argv+argc){}
                istringstream f(*p);
                if(!(f >> sharedFileName)){}
                if(sharedFileName=="") {
                    cerr << "Error: must provide a shared file." << endl;
                }
            }
            else if(strcmp(*p,"-design")==0) {
                if(++p>=argv+argc){}
                istringstream f(*p);
                if(!(f >> designFileName)){}
            }
            else if(strcmp(*p,"-minpartitions")==0) {
                if(++p>=argv+argc){}
                istringstream f(*p);
                if(!(f >> minNumPartitions)){}
            }
            else if(strcmp(*p,"-maxpartitions")==0) {
                if(++p>=argv+argc){}
                istringstream f(*p);
                if(!(f >> maxNumPartitions)){}
            }
            else if(strcmp(*p,"-optimize")==0) {
                if(++p>=argv+argc){}
                istringstream f(*p);
                if(!(f >> optimizeGap)){}
            }
            else{   
                cout << "you entered the wrong parameter" << endl;
            }
        }
    }
    

    double minLaplace = 1e10;
    int minPartition = 0;
    
    readSharedFile(sharedFileName, sharedMatrix, otuNames, sampleNames);


    if(designFileName==""){
        string fileRoot = sharedFileName.substr(0,sharedFileName.find_last_of(".")+1);
        ofstream fitData((fileRoot+"mix.fit").c_str());
        fitData.setf(ios::fixed, ios::floatfield);
        fitData.setf(ios::showpoint);
     
        cout << "K\tNLE\t\tlogDet\tBIC\t\tAIC\t\tLaplace" << endl;
        fitData << "K\tNLE\tlogDet\tBIC\tAIC\tLaplace" << endl;

        for(int numPartitions=1;numPartitions<=maxNumPartitions;numPartitions++){
            qFinderDMM findQ(sharedMatrix, numPartitions);
            
            double laplace = findQ.getLaplace();
            cout << numPartitions << '\t';
            cout << setprecision (2) << findQ.getNLL() << '\t' << findQ.getLogDet() << '\t';
            cout << findQ.getBIC() << '\t' << findQ.getAIC() << '\t' << laplace;
            
            fitData << numPartitions << '\t';
            fitData << setprecision (2) << findQ.getNLL() << '\t' << findQ.getLogDet() << '\t';
            fitData << findQ.getBIC() << '\t' << findQ.getAIC() << '\t' << laplace << endl;

            if(laplace < minLaplace){
                minPartition = numPartitions;
                minLaplace = laplace;
                cout << "***";
            }
            cout << endl;
            
            findQ.printZMatrix(fileRoot+toString(numPartitions)+"mix.posterior", sampleNames);
            findQ.printRelAbund(fileRoot+toString(numPartitions)+"mix.relabund", otuNames);

            if(optimizeGap != -1 && (numPartitions - minPartition) >= optimizeGap && numPartitions >= minNumPartitions){ break;  }
        }
        fitData.close();

        generateSummaryFile(minPartition, fileRoot);
    }
    else{
        string fileRoot = designFileName.substr(0,designFileName.find_last_of(".")+1);
        vector<vector<double> > partitions;

        readDesignFile(designFileName, sampleNames, partitions);
        qFinderDMM findQ(sharedMatrix, partitions);
        
        double laplace = findQ.getLaplace();

        ofstream fitData((fileRoot+"fit").c_str());
        fitData.setf(ios::fixed, ios::floatfield);
        fitData.setf(ios::showpoint);
        
        cout << "K\tNLE\t\tlogDet\tBIC\t\tAIC\t\tLaplace" << endl;
        fitData << "K\tNLE\tlogDet\tBIC\tAIC\tLaplace" << endl;

        cout << partitions.size() << '\t';
        cout << setprecision (2) << findQ.getNLL() << '\t' << findQ.getLogDet() << '\t';
        cout << findQ.getBIC() << '\t' << findQ.getAIC() << '\t' << laplace << endl;
        
        fitData << partitions.size() << '\t';
        fitData << setprecision (2) << findQ.getNLL() << '\t' << findQ.getLogDet() << '\t';
        fitData << findQ.getBIC() << '\t' << findQ.getAIC() << '\t' << laplace << endl;
        fitData.close();
    }
    
}

/**************************************************************************************************/
