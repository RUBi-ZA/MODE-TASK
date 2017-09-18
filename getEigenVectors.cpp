// getEigenVectors.cpp
// Extracts the eigenvectors for a specfied normal mode
// Author: Caroline Ross: caroross299@gmail.com
// August 2017

//Run individually for each mode

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include <ap.h>
#include <time.h>

using namespace std;

// Get current date/time, format is YYYY-MM-DD HH:mm:ss
const std::string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);

    strftime(buf, sizeof(buf), "%Y-%m-%d %X", &tstruct);

    return buf;
}

int main(int argc, char *argv[])
{
	//Init vars
	string protomerVT, totalRes, firstMode, outdir = "output";
	int direction = 1; //direction of overlap correction - 1 is default - user will determine this from Conformational Analysis
	bool hasVt = false, hasRes = false, hasMode = false;
	
	// Welcome message
	cout<< "============================================================\n"<<endl;
	
	cout<< "\t:-) >>------->"<<argv[0]<<"<-------<< (-:\t\n"<<endl;

	cout<< "\tAuthor(s): Caroline Ross (caroross299@gmail.com)\t\t\t\t"<<endl;
	cout<< "\tResearch Unit in Bioinformatics (RUBi)\t\t"<<endl;
	cout<< "\tRhodes University, 2017\t\t\t\t"<<endl;
	cout<< "\tDistributed under GNU GPL 3.0\t\t\t\n"<<endl;
	cout<< "\thttps://github.com/michaelglenister/NMA-TASK\t\n"<<endl;

	cout<< "============================================================"<<endl;

	// Begin parameter handling
	// Add more else if statements for further parameters
	int i;
    for(i=0; i<argc; ++i)
    {
		if (strcmp(argv[i], "-h") == 0)
		{
			cout<<"usage: getEigenVectors [-h] [--vt VTMATRIX] [--mode MODEVALUE]"<<endl;
			cout<<"		[--outdir DIRECTORY] [--direction INTEGER]"<<endl;
			cout<<"arguments:"<<endl;
			cout<<" -h, --help		Show this help message and exit"<<endl;
			cout<<" --vt			VT matrix file"<<endl;
			cout<<" --mode			Mode value"<<endl;
			cout<<" --outdir		Directory to generate output to"<<endl;
			cout<<" --direction		Direction of overlap correction (Default -1)"<<endl;
			return -1;
		}
		else if(strcmp(argv[i], "--vt") == 0)
		{
			protomerVT = argv[i+1];
			hasVt = true;
		}
	    else if(strcmp(argv[i], "--mode") == 0)
		{
			firstMode = argv[i+1];
			hasMode = true;
		}
		else if(strcmp(argv[i], "--outdir") == 0)
		{
			outdir = argv[i+1];
			//hasOutdir = true;
		}

		else if(strcmp(argv[i], "--direction") == 0)
		{
			direction = atof(argv[i+1]);
			if (direction!=1&&direction!=-1)
			{	
				direction =1;
				cout<< "\n***********************************\nWARNING DEFAULT DIRECTION = 1 HAS BEEN USED\nDirection must be specified as -1 or 1\n***********************************\n"<<endl;}//if
			//hasOutdir = true;
		}
    }

	if(!hasVt)
	{
		cout<<"A VT matrix file is required, use '-h' to view help"<<endl;
		return -1;
	}
	else if(!hasMode)
	{
		cout<<"A mode value is required, use '-h' to view help"<<endl;
		return -1;
	}

	char newline = '\n';
	ifstream vfile (protomerVT.c_str());// This is the PDB file that will be coarse grained. 
	if (vfile.good() == false) {
        	cout<<newline<<"**************************************"<<newline<<"ERROR: Specified VectorMatrix file does not exist. Exiting..."<<newline<<"**************************************"<<newline;
		return -1;
        }//if
        vfile.close();

	string protomerMode = outdir + "/EVectors"+firstMode.c_str()+".txt"; 
	
	int total = atoi(totalRes.c_str());
	int mode1 = atoi(firstMode.c_str())-1;


	// End parameter handling

	// Start cronometer
	const int ONE_HOUR = 60 * 60;
	const int ONE_MINUE = 60;

	int hour;
	int min;
	int sec;
	std::cout << "Started at: " << currentDateTime() << std::endl;
	clock_t tStart = clock();

	
	
	ifstream mfile (protomerVT.c_str());

	string lineC;
	total = 0;
	while (!mfile.eof())
	{
		getline (mfile,lineC);

		if(!lineC.empty())
		{
			total++;
		}
	}//while
	
	
	mfile.close();
	mfile.clear();

	ifstream afile (protomerVT.c_str());
	double eigenVectors[total/3][3];
	int vectorIndex = 0;
	string line;
	string element;
	
	double ve;

	int countLines = 0;

	while (!afile.eof())
	{
		getline (afile,line);

		if(!line.empty())
		{
			countLines++;
		}
		if (countLines-1==mode1)
		{
			istringstream iss(line);
			for (int i = 0; i<total/3; i++)
			{
				for (int j= 0; j<3; j++)
				{
					iss>>element;
					ve = atof(element.c_str());
					eigenVectors[i][j]=ve;
				}//for
			}//for
			break;
		}//if
	}//while
	afile.close();
	afile.clear();
	//convert all vectors to unit vectors
	double magnitude;

	for (int i = 0; i<total/3; i++)
	{
		magnitude=0;
		for (int j= 0; j<3; j++)
		{
			magnitude = magnitude+(eigenVectors[i][j]*eigenVectors[i][j]);
		}//for j1

		
		magnitude = sqrt(magnitude);

		for (int j= 0; j<3; j++)
		{
			eigenVectors[i][j]=direction*eigenVectors[i][j]/magnitude;
		}//for j2
	
	}//for i

	// write ModeVectors to file. 
	ofstream outputFileW;
	outputFileW.open(protomerMode.c_str());

	for (int i = 0; i<total/3; i++)
	{
		for (int j= 0; j<3; j++)
		{
			outputFileW<<eigenVectors[i][j]<<" ";
		}//for j1

		outputFileW<<endl;
	}//for i 

	outputFileW.close();

	// End cronometer
	std::cout << "Completed at: " << currentDateTime() << std::endl;
	int time_target=(clock() - tStart)/CLOCKS_PER_SEC;

	hour=time_target/ONE_HOUR;
	time_target-=hour*ONE_HOUR;
	min=time_target/ONE_MINUE;
	time_target-=min*ONE_HOUR;
	sec=time_target;
	if (min<10 && sec<10)
	{
		printf("- Total time: %d:0%d:0%d\n",hour,min,sec);
	}
	else if (min<10)
	{
		printf("- Total time: %d:0%d:%d\n",hour,min,sec);
	}
	else if (sec<10)
	{
		printf("- Total time: %d:%d:0%d\n",hour,min,sec);
	}
	else
	{
		printf("- Total time: %d:%d:%d\n",hour,min,sec);
	}

	return 0;
}//main
