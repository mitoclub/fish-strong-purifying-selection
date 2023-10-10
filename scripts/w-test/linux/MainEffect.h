//snp, k, w, p-value
//output result<p threashold
#include"X2operation.h"
#include "chisqr.c"
#include "gamma.c"
void CalculateMainEffect()
{
	double p;
	cout<<"input p value threashold:\np = ";
	cin>>p;
	string outputFileName,tempFileName,tempExecuteCommand;
	cout<<"Result output file name:\nfile name = ";
	cin>>outputFileName;
	
	tempFileName=GetTempFileName();
	ofstream fout(tempFileName.c_str());
		
	int i,j,df;
	double t, n, p0, p1, N0r, N1r, a[3], b[3], w, ptemp, chisqrtemp; //a saves control; b saves case
	bool myflag;
	
	
	for (i=0;i<row;i++)
		{
			a[0]=a[1]=a[2]=b[0]=b[1]=b[2]=0; //initialization
			for (j=0;j<column;j++)
				if (y[j])
					++b[matrixOri[i][j]];
					else
						++a[matrixOri[i][j]];
			myflag=false;
			for (j=0;j<3;j++)
				if (((a[j]==0) && (b[j]!=0)) || ((b[j]==0) && (a[j]!=0)))
					{
						myflag=true;
						break;
					}
			n=0.0;
			if (myflag)
				{
					for (j=0;j<3;j++)
						if (!((a[j]==0)&&(b[j]==0)))
							{
								a[j]+=0.5;
								b[j]+=0.5;
								n+=0.5;
							}
				}
			N0r=1/(1.0*N0+n);
			N1r=1/(1.0*N1+n);
			w=0; df=2;
			for (j=0;j<3;j++)
				if (a[j])
					{
						p0=a[j]*N0r;
						p1=b[j]*N1r;
						t=log((1-p0)*p1/(p0*(1-p1)))/sqrt(1.0/(1.0*a[j])+1.0/(1.0*b[j])+1.0/(1.0*(N0+n-a[j]))+1.0/(1.0*(N1+n-b[j])));
						w+=t*t;
					}
					else
						--df;
			w=w*h1[df];
			chisqrtemp=	chisqr(df,w);
			if (chisqrtemp<=p)
				fout<<SNPname[i]<<"\t"<<df+1<<"\t"<<w<<"\t"<<chisqrtemp<<endl;
		}
	fout.close();
	cout<<"k\th\tf"<<endl;
	for (i=1;i<3;i++)
		cout<<i<<"\t"<<h1[i]<<"\t"<<f1[i]<<endl;

	tempExecuteCommand = "sort -k4g "+tempFileName+" | sed '1s/^/SNP\tk\tW\tp.value\n/' > "+outputFileName;
	system(tempExecuteCommand.c_str());
	
	tempFileName="rm "+ tempFileName;
	system(tempFileName.c_str());
}
