
/*
THE DATA should be equilibrated apriori or define the Skiprows variable while reading the file. 

Generates the autocorrelation distribution (Correlation.out) and a "Results.out" which contains the average, error and autocorrelation time.  

*/

#include "Autocorrelation_Analysis2.h"

int main(void)
{

InputData=readIn2dData("energy_T=0.08.dat",0);
cout<<"The input data has "<<n<<" rows and "<<m<<" columns\n";

cout<<"Performing the autocorrelation analysis...\n";
autocorrelation();

cout<<"Printing the output\n";

PrintResults();




}


void PrintResults()
{
	//Creating the output file with the average, error and autocorrelation time.
	FILE *ptr_file;
	ptr_file=fopen("Results.out","w");
	int i;
	for (i=0;i<m;i++)
		{
		fprintf(ptr_file, " %lf %lf %lf \n" ,Results[0][i],Results[1][i],Results[2][i]);
		}

	fclose(ptr_file);
}

Vec average(){
int i,j;
Vec Average(m);

for (j=0;j<m;j++)
{
	double Result=0;
	for(i=0; i<n; i++)
		{
		Result=Result+InputData[i][j];
		}
	Average[j]=Result/n;
	
}

return(Average);

}

void autocorrelation(){
//Based on my python function.
	FILE *ptr_file;
	ptr_file=fopen("Correlation.out","w");
	int i,j, k;
	double sum;
	double val;
	int sign;
        double tol=0.0001;
        Vec Average, Error(m), CorrelationTime(m),Correlation(n);
        
	Average=average();
	Results.push_back(Average);
	
	for (j=0;j<m;j++)
	{
	
		k=0;
		val=10.0;
		sign=1;
		CorrelationTime[j]=0.5;  //To add the 0.5 of the autocorrelation time.
		
		while (val>tol && sign==1)
		{
			sum=0;
			for(i=0; i<n-k; i++)
			{
				sum+=(InputData[i][j]-Average[j])*(InputData[i+k][j]-Average[j]);

			}
			Correlation[k]=sum/(n-k);
			val=Correlation[k];
			if (k>0 && (Correlation[k]-Correlation[k-1])>0)
			{
				sign=0; //The tolerance should be satisfied only in the positive side.
			}	
			//fprintf(ptr_file, " %lf\n" ,correlation[k]);
			CorrelationTime[j]=CorrelationTime[j]+Correlation[k]/Correlation[0];
			k++;
		}
	Error[j]=sqrt(Correlation[0]*2*CorrelationTime[j]/(n+1));
	}
Results.push_back(Error);
Results.push_back(CorrelationTime);
	


fclose(ptr_file);
//printf("Correlation time %lf\n" ,results[1]);
//printf("Error %lf\n" ,results[2]);
}


vector<vector<double> > readIn2dData(const char* filename,int Skiprows)
{
/* Function takes a char* filename argument and returns a 
* 2d dynamic array containing the data
* Then after the matrix is created:



table.size() is the number of rows

table[i] is row i

table[i].size() is the number of cols. in row i

table[i][j] is the element at the j-th col. of row i

*/


	vector< vector<double> > table; 
	fstream ifs;

	//open file  
	ifs.open(filename);
	int i=0;  //Skiprow counter
	while (true)
	{
		string line;
		double buf;
		getline(ifs, line);
		//cout<<line<<"\n";
		i++;

		stringstream ss(line, ios_base::out|ios_base::in|ios_base::binary);

		if (!ifs)// mainly catch EOF
		    
		    break;
		if (line[0] == '#' || line.empty()||i <= Skiprows ) // catch empty lines or comment lines and skiprows
		    continue;

		vector<double> row;

		while (ss >> buf)
		    row.push_back(buf);
		    
		table.push_back(row);
		
	}

	ifs.close();
	n=table.size();
	m=table[0].size(); 

	return table;
}


