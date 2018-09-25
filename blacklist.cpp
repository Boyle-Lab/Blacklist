//
//  blacklist.cpp
//  blacklist generation for ENCODE
//
//  Author: Alan Boyle
//  Blacklist Copyright (c) 2018 Alan Boyle
//  This program comes with ABSOLUTELY NO WARRANTY.
//  This is free software, and you are welcome to redistribute it
//  under certain conditions.

#include <iostream>
#include <fstream>
#include <vector>
#include <stdint.h>
#include <math.h>
#include <dirent.h>
#include <string>
#include <algorithm>
#include <iterator>
#include "api/BamReader.h"
using namespace std;


class SequenceData {
	public:
		std::vector<int> binsInput;
		std::vector<int> binsSpikes;
		std::vector<int> binsMultimapping;
		std::vector<int> binsTemp;
		int totalReads;
		int zeroMulti = 0;

		SequenceData(std::string bamFile, std::string bamIndexFile) {
			this->bamFile = bamFile;
			this->bamIndexFile = bamIndexFile;
		}

		SequenceData() {
		}

		void setBamFile(std::string bamFile) {
			this->bamFile = bamFile;
		}

		std::string getBamFile() {
			return bamFile;
		}

		void setBamIndexFile(std::string bamIndexFile) {
			this->bamIndexFile = bamIndexFile;
		}

		void getInputBins(std::vector<uint8_t>& mappability, int binSize, int binOverlap, std::string refName)
		{
			this->binsInput.clear();
			this->binsSpikes.clear();
			this->binsMultimapping.clear();
			this->binsTemp.clear();

			BamTools::BamReader reader;
			BamTools::BamAlignment al;
			std::vector<int> tempCounts (mappability.size(), 0);
			std::vector<int> multiCounts (mappability.size(), 0);
			totalReads = 0;
			int testCntr = 0;
			int zmultiCntr = 0;

			reader.Open(this->bamFile);
			reader.OpenIndex(this->bamIndexFile); // Significant speedup here

 		  	// This restricts output to a specific chromosome
			if(reader.GetReferenceID(refName) == -1) {
				//cout << "No Reads to load!\n"; //DEBUG
			} else {
			  	reader.SetRegion(reader.GetReferenceID(refName), 0, reader.GetReferenceID(refName), mappability.size());

 				 // Keep vector of locations and counts of reads
 			  	// filter on mappability based on each read length
				while(reader.GetNextAlignmentCore(al)) {
			        	if(mappability[al.Position] > 0 && mappability[al.Position] <= al.Length) { // Filter mappability
			        		tempCounts[al.Position] = tempCounts[al.Position] + 1; // Count reads at each position that are uniquely mappable
							if(tempCounts[al.Position] == 1) { // For lambda calculation
								testCntr++;
							}
 		      			} else { // A Multimapping read
		            			multiCounts[al.Position] = multiCounts[al.Position] + 1; // Count reads at each position that are not uniquely mappable
						zmultiCntr++;
					}
					totalReads++;
				}
			}

			reader.Close();

			if(zmultiCntr < 10000) {
				zeroMulti++;
				//cerr << "No multimapping" << endl;
			}

			//Poisson lambda threshold calculation
			float lambda = (float)totalReads / (float)testCntr;
			float lambdaThreshold = lambda + 4.0 * sqrt(lambda);
			int thresh = (int)ceil(lambdaThreshold);
			if(totalReads == 0) {
				thresh = 0;
			}


			int readCntr;	//collapsed reads per bin
			int spikeCntr;	//Sikes above thresh per bin
			int multiCntr; //Multimapping reads per bin
			for(int i = 0; i < tempCounts.size() - binSize; i+=binOverlap) {
				readCntr = 0;
				spikeCntr = 0;
				multiCntr = 0;
				for(int j = 0; j < binSize; j++) {
					if(tempCounts[i+j] > 0) {
						readCntr += tempCounts[i+j];
					}
					if(tempCounts[i+j] > thresh) {
						spikeCntr++;
					}
					if(multiCounts[i+j] > 0) {
						multiCntr += multiCounts[i+j];
					}
				}

				this->binsInput.push_back(readCntr);
				this->binsSpikes.push_back(spikeCntr);
				this->binsMultimapping.push_back(multiCntr);
				this->binsTemp.push_back(0);
			}

		}

	private:
		std::string bamFile;
		std::string bamIndexFile;


};

// Takes as input the mappability at each base and returns
//  binned counts of mappability
void getMappabilityBins(std::vector<int>& retVal, std::vector<uint8_t>& mappability, int binSize, int binOverlap)
{
	std::vector<int> v;
	int uniqueCntr;
	int uniqueLength = 36; //This is arbitraty and defines how long a read needs
				// to be to be considered unique
				// Should be set to something actually calculated in the uint8 files

	for(int i = 0; i <= mappability.size() - binSize; i+=binOverlap) {
		uniqueCntr = 0;
		for(int j = 0; j < binSize; j++) {
			if(mappability[i+j] > 0 && mappability[i+j] <= uniqueLength) {
				uniqueCntr++;
			}
		}
		v.push_back(uniqueCntr);
	}

	retVal.swap(v);
}

int getdir(std::string dirname, std::vector<string> & files, std::string filetype)
{
	DIR *dir;
	int pos;
	struct dirent *ent;
	dir = opendir(dirname.c_str());
	if (dir != NULL) {
	  	while ((ent = readdir (dir)) != NULL) {
			pos = strlen(ent->d_name) - 4;
			if (! strcmp(&ent->d_name[pos], filetype.c_str())) {
				//printf("%s\n", ent->d_name); //DEBUG
				files.push_back(string(ent->d_name));
			}
	  	}
  		closedir (dir);
	} else {
 	 	/* could not open directory */
		cerr << "Unable to read input files!" << endl;
		exit(1);
	}
}

double quantile(std::vector<double> v, double q) {
	sort(v.begin(), v.end());
	double h = ((v.size() - 1) * q) + 1;
	if(h >= v.size()) {
		h = v.size() - 1;
	}
	return v[floor(h)] + ((h - floor(h)) * (v[floor(h) + 1] - v[floor(h)]));
}

int quantile(std::vector<int> v, double q) {
	sort(v.begin(), v.end());
	double h = (((double)v.size() - 1) * q) + 1;
	return floor((double)v[floor(h)] + ((h - floor(h)) * ((double)v[floor(h) + 1] - (double)v[floor(h)])));
}

// Function to find rank
void rankify(std::vector<double>& A) {
	int n = A.size();
	std::vector<double> R(n,0);
	std::vector<std::tuple<double, int>> T;
	int r = 1;

	// Create array of tuples storing value and index
	for(int j = 0; j < n; j++) {
		T.push_back(std::make_tuple(A[j], j));
	}

	// Sort tubples by data value
	std::sort(begin(T), end(T), [](auto const &t1, auto const &t2) {
        	return get<0>(t1) < get<0>(t2); // or use a custom compare function
	});

	int i = 0;
	int index, j;
	while(i < n) {
		j = i;

		// Get elements of same rank
		while(j < n && std::get<0>(T[j]) == std::get<0>(T[j+1])) {
			j++;
		}

		int m = j - i + 1;

		for(j = 0; j < m; j++) {
			// For each equal element use .5
			index = std::get<1>(T[i+j]);
			R[index] = r + (m-1)*0.5;
		}

		// Increment rank and index
		r+=m;
		i+=m;
	}

	A.swap(R);
}

void quantileNormalize(std::vector<std::vector<double>>& data) {
	int cellCount = data.size();
	int binCount = data[0].size();

	//First calculate rank means
	std::vector<double> rankedMean(binCount,0);
	for(int cellID = 0; cellID < cellCount; cellID++) {
		std::vector<double> x(binCount,0);
		for(int binID = 0; binID < binCount; binID++) {
			x[binID] = data[cellID][binID];
		}

		sort(x.begin(), x.end());

		for(int binID = 0; binID < binCount; binID++) {
			rankedMean[binID] += x[binID];
		}
	}
	for(int binID = 0; binID < binCount; binID++) {
		rankedMean[binID] /= (double)cellCount;
	}

	//calculate half value for ties
	std::vector<double> rankedMeanTie(binCount-1,0);
	for(int binID = 0; binID < (binCount-1); binID++) {
		rankedMeanTie[binID] = ((rankedMean[binID]+rankedMean[binID+1])/2);
	}

	//Iterate through each cell line
	for(int s = 0; s < cellCount; s++) {
		std::vector<double> bins(binCount,0);
		for(int p = 0; p < binCount; p++) {
			bins[p] = data[s][p];
		}
		rankify(bins);

		std::vector<double> binsQuantileNormalized(binCount, 0);
		for(int p = 0; p < binCount; p++) {
			if(std::fmod(bins[p],1) != 0) {
				binsQuantileNormalized[p] = rankedMeanTie[(int)floor(bins[p])-1];
			} else {
				binsQuantileNormalized[p] = rankedMean[(int)(bins[p]-1)];
			}

			data[s][p] = binsQuantileNormalized[p];
		}
	}

}



std::vector<double> getAbnormalRegions(std::vector<SequenceData> inputData, std::vector<int> binsMap, int type, bool normalize) {
	// Types:
	//  1 - Read
	//		- normalize = Reads / mapability
	//  2 - Multimapping
	//		- normalize = multimapping reads / total reads
	//  3 - Spike (not used)

	std::vector<double> normForQuantile;
	std::vector<double> result(binsMap.size());
	double normTemp;
	double quantileVal = 0;
	std::vector<std::vector<double>> data(inputData.size(), std::vector<double>(binsMap.size()));

	for(int j = 0; j < inputData.size(); j++) { //process by column

		// Generalize this a little
		if(type == 1) { //reads
			inputData[j].binsTemp = inputData[j].binsInput;
		} else if (type == 2) { //multimapping
			inputData[j].binsTemp = inputData[j].binsMultimapping;
		} else if (type == 3) { //spikes
			inputData[j].binsTemp = inputData[j].binsSpikes;
		}

		for(int i = 0; i < binsMap.size(); i++) { //track all rows in this column
			if(normalize) {
				if (type == 2) { //multimapping
					normTemp = (double)inputData[j].binsTemp[i] / (double)inputData[j].binsInput[i];
				} else {
					normTemp = (double)inputData[j].binsTemp[i] / (double)binsMap[i];
				}
			} else {
				normTemp = (double)inputData[j].binsTemp[i] / (double)inputData[j].totalReads * (double)1000000;
			}

			// Save it if we didn't divide by zero
			if(!isnan(normTemp) && !isinf(normTemp) && normTemp > 0) {
				data[j][i] = normTemp;
			} else {
				data[j][i] = 0.0;
			}
		}
	}

	quantileNormalize(data);

	//Now we collapse rows
	std::vector<double> means(inputData.size());
	for(int i = 0; i < binsMap.size(); i++) { // over each row
		for(int j = 0; j < inputData.size(); j++) { //over each column
			means[j] = data[j][i];
		}
		result[i] = quantile(means, 0.5); // 90% of cells have at least this value
	}

	return result;
}

int main(int argc, char* argv[])
{
	std::vector<uint8_t> mappability;
	std::vector<int> binsMap;
	std::vector<string> inputFileList;
	std::vector<SequenceData> inputData;
	std::string mappabilityFile;
	std::string bamFile;
	std::string bamIndexFile;
	std::string refName;
	int binCount;
	std::vector<double> readNormList;
	std::vector<double> multiList;

	//Parameters
	int binSize = 1000;
	int binOverlap = 100;

	if(argc < 2) {
		cout << "Blacklist is used to generate the ENCODE blacklists for various species." << endl;
		cout << "Usage is ./Blacklist <chr>" << endl;
		cout << "The program requires an input/ folder containing indexed bam files." << endl;
		cout << "The program requires a mappability/ folder containing Umap mappability files." << endl;
		exit(0);
	} else {
		refName = argv[1];
	}

	// Load mappability file -- we need this multiple times
	mappabilityFile = "mappability/" + refName + ".uint8.unique";
	uint8_t temp;
	std::ifstream infile;
	infile.open(mappabilityFile.c_str());
	if(!infile.good()) {
		cerr << "Unable to read mappability files!" << endl;
		exit(1);
	}
	infile.unsetf(ios_base::skipws);
	while(infile >> temp && !infile.eof()) {
		mappability.push_back(temp);
	}

	getMappabilityBins(binsMap, mappability, binSize, binOverlap);

	//Here we need to accumulate putative sites based on spikes and input levels
	getdir("input/", inputFileList, ".bam");
	for(int i = 0; i < inputFileList.size(); i++) {
		bamFile = "input/" + inputFileList[i];
		bamIndexFile = "input/" + inputFileList[i] + ".bai";
		inputData.push_back(SequenceData(bamFile, bamIndexFile));
		//cerr << "Processing " << bamFile << endl;
		inputData.back().getInputBins(mappability, binSize, binOverlap, refName);
	}

	//Fix error of one larger bins
	if(binsMap.size() > inputData[0].binsTemp.size()) {
		binsMap.pop_back();
	}

	readNormList = getAbnormalRegions(inputData, binsMap, 1, 1);
	multiList = getAbnormalRegions(inputData, binsMap, 2, 0);

	// for testing purposes, this loop will output the scores for every genomic region
	//for(int i = 0; i < binsMap.size(); i++) {
	//	cout << refName << "\t" << i*binOverlap << "\t" << i*binOverlap + binSize << "\t" << readNormList[i] << "\t" << multiList[i] << "\n";
	//}

	// Output the raw count information
	int tempFirst;
	int tempLast;
	int inRegion = 0;
	int tempHit = 0;
	int hitCode = 0;
	int hitCounter = 0;
	int miss = 0;

	// Generate the threshold levels for weak and strong hits
	double readWeakThresh = quantile(readNormList, 0.95);
	double readStrongThresh = quantile(readNormList, 0.99);
	double multiWeakThresh = quantile(multiList, 0.95);
	double multiStrongThresh = quantile(multiList, 0.99);
	double minThresh = *std::min_element( std::begin(readNormList), std::end(readNormList) ); // We can bridge regions that have 0 signal as well.
			      // Optional to use flag these as blacklist.

	//cerr << readWeakThresh << " " << readStrongThresh << " " << multiWeakThresh << " " << multiStrongThresh << " " << minThresh << endl;

	for(int i = 0; i < binsMap.size(); i++) {
		if(readNormList[i] >= readWeakThresh || multiList[i] >= multiWeakThresh || readNormList[i] <= minThresh) {
			//tracking for overlapping regions
			miss = 0;
			tempLast = i;

			// If this is a new region, record it
			if(inRegion == 0) {
				tempFirst = i;
				inRegion = 1;
			}

			// Check to see if this bin passes a threshold
			if(readNormList[i] >= readStrongThresh) {
				tempHit = 1;
				hitCode = hitCode | 1;
				hitCounter++;
			}

			if( multiList[i] > multiStrongThresh) {
				tempHit = 1;
				hitCode = hitCode | 2;
				hitCounter++;
			}
		} else {
			if(miss < (binSize/binOverlap + 5)) { // bridge over adjacent bins plus 100 * 200 = 20kb
				miss++;
			} else { //nothing in this distance
				inRegion = 0;

				// If we hit a threshold we output the whole region
				if(tempHit == 1) {
					cout << refName << "\t" << tempFirst*binOverlap << "\t" << tempLast*binOverlap + binSize << "\t" << hitCode << "\t" << hitCounter << "\n";
					tempHit = 0;
					hitCode = 0;
					hitCounter = 0;
				}
			}
		}
	}

	// If we were in a region when when hit the end of the chromosome, output region
	// This may go past chromosome end!!
	if(tempHit == 1) {
		cout << refName << "\t" << tempFirst*binOverlap << "\t" << tempLast*binOverlap + binSize << "\t" << hitCode << "\t" << hitCounter << "\n";
	}

}
