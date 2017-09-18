
// ANM.cpp
// Performs Normal Mode Analysis using the ANM tool
// Author: Caroline Ross: caroross299@gmail.com
// August 2017


//Performs Normal Mode Analysis using the ANM tool
//The script uses the alglib library and the following must be included for the script to run. These are contained in src folder in a cpp folder of the alglib package. When complying this code use g++ -I path/cpp/src as an option
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <string>
#include <sstream>
#include <limits>
#include <specialfunctions.h>
#include <specialfunctions.cpp>
#include <iostream>
#include <linalg.h>
#include <linalg.cpp>
#include <alglibinternal.h>
#include <alglibinternal.cpp>
#include <alglibmisc.h>
#include <alglibmisc.cpp>
#include <ap.h>
#include <ap.cpp>
#include <vector>

#include <stdio.h>
#include <time.h>


using namespace std;
using namespace alglib;


int countAtoms()//counts Carbon Atoms (Beta carbons for all residues but Alpha carbons for Glycine)
{
	vector< vector<int> > row;
}//countAtoms


vector< vector<double> > getCoOrds(string pdbInput,string atype)// gets the x,y,z cordinates for all the carbon atoms in PDB file. Returns 3 X numAtoms Matrix
{	
	vector< vector<double> > C;
	vector<double> atomC;
	
        ifstream mfile (pdbInput.c_str());// This is the PDB file that will be coarse grained. 
	string li;

// file is open and not empty
	
	while (! mfile.eof() )
	{       
		getline (mfile,li);
		istringstream iss(li);
		string atom, temp, type,res,x,y,z;
		iss>>atom;
		iss>>temp;
		iss>>type;
		iss>>res;
		iss>>temp;
		iss>>temp;
		iss>>x;
		iss>>y;
		iss>>z;
		if (atom=="ATOM")
		{	
			if((res=="GLY"&&type=="CA")||type==atype)// This selects only CA atoms for CB atoms use if((res=="GLY"&&type=="CA")||type=="CB") such that CA are selected in the case of Glycine
			{
				atomC.push_back(atof(x.c_str()));
				atomC.push_back(atof(y.c_str()));
				atomC.push_back(atof(z.c_str()));
				C.push_back(atomC);
				atomC.clear();				
			}//if GLY
			
		
		}//if ATOM
	}//while
    mfile.close();
	return C;
}//cords

vector< vector<double> > getForceConstants(vector<double> atom1,vector<double> atom2, double cutoff)// returns a 9x9 vector for the interaction between two nodes
{
	vector< vector<double> > dv2k;
	vector<double> dv2kROW;

	//calculate distance squared
	double dist2 =0.0;
	for (int i =0; i<3; i++)
	{
		dist2 = dist2+((atom2[i]-atom1[i])*(atom2[i]-atom1[i]));
	}//for to instantiate distance
	
	double dv2ka;
	
	//if outside cutoff or atom1=atom2
	if (dist2>cutoff)//parameter
	{
		for(int i=0;i<3;i++)
		{		
			dv2kROW.push_back(0.0);
			dv2kROW.push_back(0.0);
			dv2kROW.push_back(0.0);
			dv2k.push_back(dv2kROW);

			//reset dv2kROW
			dv2kROW.clear();
		}//for to populate row by row: times 3 rows
							
	}//if
	else 
	{
		for (int i=0; i<3; i++)
		{
			for (int j=0; j<3; j++)
			{
				dv2ka = -((atom2[i]-atom1[i])*(atom2[j]-atom1[j]))/dist2;
				dv2kROW.push_back(dv2ka);  

			}//for
			dv2k.push_back(dv2kROW);
			dv2kROW.clear();
		}//for
	}//else
	return dv2k;

}// getForceConstants

vector< vector<double> > getHessian(vector< vector<double> > C, double cutoff)// this calls the getForceConstants to set up the Hessian Matrix, we manipulate these force constants and the populate the hessian
{
	vector< vector<double> > Hessian; //Full 3Nx3N Hessian
	
	vector<double> rowX;
	vector<double> rowY;
	vector<double> rowZ;
	
	vector< vector<double> > interaction; //9x9 vector for each atom-atom interaction

	vector< vector<double> > diagonal; //holds 9x9 vector of values on a diagonal
	vector<double> diagonalROW; //holds 9x9 vector of values on a diagonal

	vector<double> atom1;
	vector<double> atom2;

	for (int i=0; i < C.size(); i++)
	{
		//reset the diagonal
		for (int d = 0; d<3; d++)
		{
			diagonalROW.push_back(0.0);
			diagonalROW.push_back(0.0);
			diagonalROW.push_back(0.0);
			diagonal.push_back(diagonalROW);
			diagonalROW.clear();
		}// set each element = 0
		
		//get x,y,z of atom1
		atom1 = C[i];

		//Nested loop for atom-atom interactions
		for (int j=0; j < C.size(); j++)
		{
			//get x,y,z of interacting atom			
			atom2 = C[j];
			if(i!=j)
			{
				interaction = getForceConstants(atom1,atom2,cutoff);
				for(int ir = 0; ir<3; ir++)
				{
					rowX.push_back(interaction[0][ir]);
					rowY.push_back(interaction[1][ir]);
					rowZ.push_back(interaction[2][ir]);
				}//itterate RowX,Y,Z

				//itterate summation of digonals
				for (int ir = 0; ir<3; ir++)
				{
					for (int ic = 0; ic<3; ic++)
					{
						diagonal[ir][ic] = diagonal[ir][ic]-interaction[ir][ic];
								
					}// increase each element
				}// increase each element
				interaction.clear();
			}//if i!=j
			else
			{
				for(int ir = 0; ir<3; ir++)
				{
					rowX.push_back(0.0);
					rowY.push_back(0.0);
					rowZ.push_back(0.0);
				}//itterate RowX,Y,Z
			}//set diagonal
		
		}//for atoms
		
		//Update diagonal
		for(int ir = 0; ir<3; ir++)
		{
			rowX[i*3+ir] = diagonal[0][ir];
			rowY[i*3+ir] = diagonal[1][ir];
			rowZ[i*3+ir] = diagonal[2][ir];
		}//itterate RowX,Y,Z
		//Reset diagonal
		diagonal.clear();
		
		//Update Hessian by row
		Hessian.push_back(rowX);
		rowX.clear();
		Hessian.push_back(rowY);
		rowY.clear();
		Hessian.push_back(rowZ);
		rowZ.clear();
		
	}//for atoms
	
	return Hessian;
}// getHessian

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
	double cutoff = 15;// this must be user define
	string pdbInput, outdir = "output";
	string atype="X";
	bool hasPdb = false;

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
			cout<<"Help"<<endl;
			return -1;
		}
		else if(strcmp(argv[i], "--pdb") == 0)
		{
			pdbInput = argv[i+1];
			hasPdb = true;
		}
		else if(strcmp(argv[i], "--cutoff") == 0)
		{
			cutoff = atof(argv[i+1]);
			//hasCutoff = true;
		}
		else if(strcmp(argv[i], "--outdir") == 0)
		{
			outdir = argv[i+1];
			//hasOutdir = true;
		}
             
	
                else if(strcmp(argv[i], "--atomType") == 0)
		{
			atype = argv[i+1];
			//hasOutdir = true;
		}
	
	}
	
	if(!hasPdb)
	{
		cout<<"A PDB file is required, use '-h' to view help"<<endl;
		return -1;
	}
		if (atype=="X") {
        	cout<<endl<<"**************************************"<<endl<<"ERROR: ATOM TYPE MUST BE SPECIFIED"<<endl<<"Input Options:"<<endl<<"CA: to select alpha carbon atoms"<<endl<<"CB: to select beta carbon atoms"<<endl<<"Exiting..."<<endl<<"**************************************"<<endl;
		return -1;
        }//if
        ifstream mfile (pdbInput.c_str());// This is the PDB file that will be coarse grained. 
	if (mfile.good() == false) {
        	cout<<endl<<"**************************************"<<endl<<"ERROR: Specified PDB file does not exist. Exiting..."<<endl<<"**************************************"<<endl;
		return -1;
        }//if
        mfile.close();

	cutoff = cutoff * cutoff;

	string eigenvalueMatrixFile = outdir + "/W_values.txt";
	string eigenvalueVTFile = outdir + "/VT_values.txt";
	string eigenvalueUFile = outdir + "/U_values.txt";
	
	// End parameter handling

	// Start cronometer
	const int ONE_HOUR = 60 * 60;
	const int ONE_MINUE = 60;

	int hour;
	int min;
	int sec;
	std::cout << "Started at: " << currentDateTime() << std::endl;
	clock_t tStart = clock();

        if (atype!="CA" && atype!="CB")
        {
        cout<<endl<<"**************************************"<<endl<<"Unrecognised atom type"<<endl<<"Input Options:"<<endl<<"CA: to select alpha carbon atoms"<<endl<<"CB: to select beta carbon atoms"<<endl<<"**************************************"<<endl;
        return -1;
        }//if
       
	string str2=".pdb";
        size_t found = pdbInput.find(str2);
        if (found==std::string::npos)
        {
         cout<<endl<<"**************************************"<<endl<<"ERROR: INPUT FILE MUST BE A PDB FILE!"<<endl<<"**************************************"<<endl;
         return -1;
        }//if not a pdbfile


	vector< vector<double> > C = getCoOrds(pdbInput,atype);
	vector< vector<double> > Hessian = getHessian(C, cutoff);

	int size = Hessian.size();
	alglib::real_2d_array Hes;
	Hes.setlength(size,size);
	
	for(int i =size-1; i>=0; i--)
	{
		for(int j =size-1; j>=0; j--)
		{
			Hes[i][j]= Hessian[i][j];
			Hessian[i].erase(Hessian[i].begin()+j);
		}			
	}

	Hessian.clear();
	cout<<"Starting Decomposition"<<endl;
	alglib::real_1d_array w;
	alglib::real_2d_array u;
	alglib::real_2d_array vt;
	alglib::rmatrixsvd(Hes,size,size,2, 2,0,w,u,vt);

	ofstream outputFileW;
 
	int r = vt.rows();
	int c = vt.cols();
 
	outputFileW.open(eigenvalueMatrixFile.c_str());// this is the eigenvalue matrix.

	for (int i=0; i<r; i++)
	{
		double e = w(i);
		outputFileW<<i+1<<" "<<e<<endl;
	}
	outputFileW.close();

	ofstream outputFileVT;
	outputFileVT.open(eigenvalueVTFile.c_str());

	for (int i=0; i<r; i++)
	{
        for (int j=0; j<c; j++)
        {
            double e = vt(i,j); //this is the eigenvector matrix - eigenvectors are rows. 
            outputFileVT<<e<<" ";
        }//for U

        outputFileVT<<endl;
	}// for vt
	outputFileVT.close();

	ofstream outputFileU;
	outputFileU.open(eigenvalueUFile.c_str());

	for (int i=0; i<r; i++)
	{
        for (int j=0; j<c; j++)
        {
            double e = u(i,j); //this is the eigenvector matrix - eigenvectors are columns, U and VT are square matrics. 
            outputFileU<<e<<" ";
        }//for U

        outputFileU<<endl;
	}// for vt
	outputFileU.close();


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
