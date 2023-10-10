#include<stdio.h>
#include<string>
#include<fstream>
#include <time.h>
using namespace std;

void filterbyPvalue(double pthreshold, string inputFileName, string outputFileName)
{
	ifstream fin(inputFileName.c_str());
	ofstream fout(outputFileName.c_str());
	string s;
	getline(fin,s);
	fout<<s<<endl;
	string SNP1,SNP2;
	int k,k1,k2;
	double w,p_pair,wp_snp1,wp_snp2,p_snp1,p_snp2;
	while (fin>>SNP1)
		{
			fin>>SNP2>>w>>k>>p_pair>>wp_snp1>>k1>>p_snp1>>wp_snp2>>k2>>p_snp2;
			if (p_pair<=pthreshold)
				fout<<SNP1<<"\t"<<SNP2<<"\t"<<w<<"\t"<<k<<"\t"<<p_pair<<"\t"<<wp_snp1<<"\t"<<k1<<"\t"<<p_snp1<<"\t"<<wp_snp2<<"\t"<<k2<<"\t"<<p_snp2<<endl;
		}
	fin.close();
	fout.close();
}

void CalculatePairwiseInteractionsAll()
{
	//get single SNP info
	double t, n, p0, p1, N0r, N1r, a[9], b[9], p1value[maxMatrixRow],w1[maxMatrixRow];
	int i,j,k,df1[maxMatrixRow];
	bool myflag;
	
	for (i=0;i<row;i++)
		{
			a[0]=a[1]=a[2]=b[0]=b[1]=b[2]=0;
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
			w1[i]=0; df1[i]=2;
			for (j=0;j<3;j++)
				if (a[j])
					{
						p0=a[j]*N0r;
						p1=b[j]*N1r;
						t=log((1-p0)*p1/(p0*(1-p1)))/sqrt(1.0/(1.0*a[j])+1.0/(1.0*b[j])+1.0/(1.0*(N0+n-a[j]))+1.0/(1.0*(N1+n-b[j])));
						w1[i]+=t*t;
					}
					else
						--df1[i];
			w1[i]=w1[i]*h1[df1[i]];
			p1value[i]=chisqr(df1[i],w1[i]);
		}
	string outputFileName,tempFileName,tempFileName1,s;

	cout<<"Result output file name:\nfile name = ";
	cin>>outputFileName;
	double pthreshold;
	cout<<"Output p value threshold:\n p threshold = ";
	cin>>pthreshold;
	clock_t tStart = clock();
	int snp1,snp2;
	
	tempFileName=GetTempFileName();
	ofstream fout(tempFileName.c_str());
	double w;
	int df;
	for (i=0;i<row;i++)
		for(j=i+1;j<row;j++)
			{
				snp1=i;
				snp2=j;
				for (k=0;k<9;k++)
					a[k]=b[k]=0;
				for (k=0;k<N0reduced;k++)
					{
						a[0]+=__builtin_popcount(matrix[snp1][0][k] & matrix[snp2][0][k]);
						a[1]+=__builtin_popcount(matrix[snp1][0][k] & matrix[snp2][1][k]);
						a[2]+=__builtin_popcount(matrix[snp1][0][k] & matrix[snp2][2][k]);
						a[3]+=__builtin_popcount(matrix[snp1][1][k] & matrix[snp2][0][k]);
						a[4]+=__builtin_popcount(matrix[snp1][1][k] & matrix[snp2][1][k]);
						a[5]+=__builtin_popcount(matrix[snp1][1][k] & matrix[snp2][2][k]);
						a[6]+=__builtin_popcount(matrix[snp1][2][k] & matrix[snp2][0][k]);
						a[7]+=__builtin_popcount(matrix[snp1][2][k] & matrix[snp2][1][k]);
						a[8]+=__builtin_popcount(matrix[snp1][2][k] & matrix[snp2][2][k]);
					}
				for (k=N0reduced;k<N0reduced+N1reduced;k++)
					{
						b[0]+=__builtin_popcount(matrix[snp1][0][k] & matrix[snp2][0][k]);
						b[1]+=__builtin_popcount(matrix[snp1][0][k] & matrix[snp2][1][k]);
						b[2]+=__builtin_popcount(matrix[snp1][0][k] & matrix[snp2][2][k]);
						b[3]+=__builtin_popcount(matrix[snp1][1][k] & matrix[snp2][0][k]);
						b[4]+=__builtin_popcount(matrix[snp1][1][k] & matrix[snp2][1][k]);
						b[5]+=__builtin_popcount(matrix[snp1][1][k] & matrix[snp2][2][k]);
						b[6]+=__builtin_popcount(matrix[snp1][2][k] & matrix[snp2][0][k]);
						b[7]+=__builtin_popcount(matrix[snp1][2][k] & matrix[snp2][1][k]);
						b[8]+=__builtin_popcount(matrix[snp1][2][k] & matrix[snp2][2][k]);
						
					}
				myflag=false; //check if we need to increment 0.5
				for (k=0;k<9;k++)
					if (((a[k]==0) && (b[k]!=0)) || ((b[k]==0) && (a[k]!=0)))
						{
							myflag=true;
							break;
						}
				n=0.0;
				if (myflag)
					{
						for (k=0;k<9;k++)
							if (!((a[k]==0) && (b[k]==0)))
								{
									a[k]+=0.5;
									b[k]+=0.5;
									n+=0.5;
								}
					}
				N0r=1/(1.0*N0+n);
				N1r=1/(1.0*N1+n);
				w=0; df=8;
				for (k=0;k<9;k++)
					if (a[k])
						{
							p0=a[k]*N0r;
							p1=b[k]*N1r;
							t=log((1-p0)*p1/(p0*(1-p1)))/sqrt(1.0/(1.0*a[k])+1.0/(1.0*b[k])+1.0/(1.0*(N0+n-a[k]))+1.0/(1.0*(N1+n-b[k])));
							w+=t*t;
						}
						else
							--df;
				w=w*h2[df];
				//snp1, snp2, w, k, p-value, w1, k1, p1, w2,  k2,  p2
				fout<<SNPname[snp1]<<"\t"<<SNPname[snp2]<<"\t"<<w<<"\t"<<df+1<<"\t"<<0<<"\t"<<w1[snp1]<<"\t"<<df1[snp1]<<"\t"<<p1value[snp1]<<"\t"<<w1[snp2]<<"\t"<<df1[snp2]<<"\t"<<p1value[snp2]<<endl;
			}
	tempFileName1=GetTempFileName();
	CreatePFillInRProgram(tempFileName,tempFileName1);
	s="chmod a+x p.fill.in.pair.R";
	system(s.c_str());
	system("R --slave < p.fill.in.pair.R");
	string tempFileName2=GetTempFileName();
	s="sort -k5g "+tempFileName1+" | sed '1s/^/SNP1\tSNP2\tW\tk\tp.pair\tW.snp1\tk1\tp.snp1\tW.snp2\tk2\tp.snp2\n/' > "+tempFileName2;
	filterbyPvalue(pthreshold,tempFileName2,outputFileName);	
	system(s.c_str());
	cout<<"k\th\tf"<<endl;
	for (i=1;i<9;i++)
		cout<<i+1<<"\t"<<h2[i]<<" "<<f2[i]<<endl;
	system("rm p.fill.in.pair.R");
	tempFileName="rm "+ tempFileName;
	system(tempFileName.c_str());
	tempFileName1="rm "+ tempFileName1;
	system(tempFileName1.c_str());
	tempFileName2="rm "+ tempFileName2;
	system(tempFileName2.c_str());
	i=(int)(clock()-tStart)/CLOCKS_PER_SEC;
	cout<<"Time taken: "<<i/3600<<" hours "<<i%3600/60<<" minutes "<<i%60<<" seconds"<<endl;
}


void CalculatePairwiseInteractionsExhaustively()
{
	//get single SNP info
	double t, n, p0, p1, N0r, N1r, a[9], b[9], p1value[maxMatrixRow],w1[maxMatrixRow];
	int i,j,k,df1[maxMatrixRow];
	bool myflag;
	
	for (i=0;i<row;i++)
		{
			a[0]=a[1]=a[2]=b[0]=b[1]=b[2]=0;
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
			w1[i]=0; df1[i]=2;
			for (j=0;j<3;j++)
				if (a[j])
					{
						p0=a[j]*N0r;
						p1=b[j]*N1r;
						t=log((1-p0)*p1/(p0*(1-p1)))/sqrt(1.0/(1.0*a[j])+1.0/(1.0*b[j])+1.0/(1.0*(N0+n-a[j]))+1.0/(1.0*(N1+n-b[j])));
						w1[i]+=t*t;
					}
					else
						--df1[i];
			w1[i]=w1[i]*h1[df1[i]];
			p1value[i]=chisqr(df1[i],w1[i]);
		}
	string inputFileName,outputFileName,tempFileName,tempFileName1,s;
	cout<<"Input file cotains list of markers.\nfile name = ";
	cin>>inputFileName;
	cout<<"Result output file name:\nfile name = ";
	cin>>outputFileName;
	clock_t tStart = clock();
	int nmarker,snp1,snp2,markers[maxMatrixRow];
	
	ifstream fin(inputFileName.c_str());
	nmarker=0;
	getline(fin,s);//remove the header
	while (fin>>s)
		{
			for (i=0;i<row;i++)
				if (s==SNPname[i])
					{
						markers[nmarker++]=i;
						break;
					}
			getline(fin,s);
		}
		
	if (nmarker<2)
		{
			cout<<"Bad input file"<<endl;
			return;
		}
	
	tempFileName=GetTempFileName();
	ofstream fout(tempFileName.c_str());
	double w;
	int df;
	for (i=0;i<nmarker;i++)
		for(j=i+1;j<nmarker;j++)
			{
				snp1=markers[i];
				snp2=markers[j];
				for (k=0;k<9;k++)
					a[k]=b[k]=0;
				for (k=0;k<N0reduced;k++)
					{
						a[0]+=__builtin_popcount(matrix[snp1][0][k] & matrix[snp2][0][k]);
						a[1]+=__builtin_popcount(matrix[snp1][0][k] & matrix[snp2][1][k]);
						a[2]+=__builtin_popcount(matrix[snp1][0][k] & matrix[snp2][2][k]);
						a[3]+=__builtin_popcount(matrix[snp1][1][k] & matrix[snp2][0][k]);
						a[4]+=__builtin_popcount(matrix[snp1][1][k] & matrix[snp2][1][k]);
						a[5]+=__builtin_popcount(matrix[snp1][1][k] & matrix[snp2][2][k]);
						a[6]+=__builtin_popcount(matrix[snp1][2][k] & matrix[snp2][0][k]);
						a[7]+=__builtin_popcount(matrix[snp1][2][k] & matrix[snp2][1][k]);
						a[8]+=__builtin_popcount(matrix[snp1][2][k] & matrix[snp2][2][k]);
					}
				for (k=N0reduced;k<N0reduced+N1reduced;k++)
					{
						b[0]+=__builtin_popcount(matrix[snp1][0][k] & matrix[snp2][0][k]);
						b[1]+=__builtin_popcount(matrix[snp1][0][k] & matrix[snp2][1][k]);
						b[2]+=__builtin_popcount(matrix[snp1][0][k] & matrix[snp2][2][k]);
						b[3]+=__builtin_popcount(matrix[snp1][1][k] & matrix[snp2][0][k]);
						b[4]+=__builtin_popcount(matrix[snp1][1][k] & matrix[snp2][1][k]);
						b[5]+=__builtin_popcount(matrix[snp1][1][k] & matrix[snp2][2][k]);
						b[6]+=__builtin_popcount(matrix[snp1][2][k] & matrix[snp2][0][k]);
						b[7]+=__builtin_popcount(matrix[snp1][2][k] & matrix[snp2][1][k]);
						b[8]+=__builtin_popcount(matrix[snp1][2][k] & matrix[snp2][2][k]);
						
					}
				myflag=false; //check if we need to increment 0.5
				for (k=0;k<9;k++)
					if (((a[k]==0) && (b[k]!=0)) || ((b[k]==0) && (a[k]!=0)))
						{
							myflag=true;
							break;
						}
				n=0.0;
				if (myflag)
					{
						for (k=0;k<9;k++)
							if (!((a[k]==0) && (b[k]==0)))
								{
									a[k]+=0.5;
									b[k]+=0.5;
									n+=0.5;
								}
					}
				N0r=1/(1.0*N0+n);
				N1r=1/(1.0*N1+n);
				w=0; df=8;
				for (k=0;k<9;k++)
					if (a[k])
						{
							p0=a[k]*N0r;
							p1=b[k]*N1r;
							t=log((1-p0)*p1/(p0*(1-p1)))/sqrt(1.0/(1.0*a[k])+1.0/(1.0*b[k])+1.0/(1.0*(N0+n-a[k]))+1.0/(1.0*(N1+n-b[k])));
							w+=t*t;
						}
						else
							--df;
				w=w*h2[df];
				//snp1, snp2, w, k, p-value, w1, k1, p1, w2,  k2,  p2
				fout<<SNPname[snp1]<<"\t"<<SNPname[snp2]<<"\t"<<w<<"\t"<<df+1<<"\t"<<0<<"\t"<<w1[snp1]<<"\t"<<df1[snp1]<<"\t"<<p1value[snp1]<<"\t"<<w1[snp2]<<"\t"<<df1[snp2]<<"\t"<<p1value[snp2]<<endl;
			}
	fin.close();
	fout.close();
	tempFileName1=GetTempFileName();
	CreatePFillInRProgram(tempFileName,tempFileName1);
	s="chmod a+x p.fill.in.pair.R";
	system(s.c_str());
	system("R --slave < p.fill.in.pair.R");
	s="sort -k5g "+tempFileName1+" | sed '1s/^/SNP1\tSNP2\tW\tk\tp.pair\tW.snp1\tk1\tp.snp1\tW.snp2\tk2\tp.snp2\n/' > "+outputFileName;
	system(s.c_str());
	cout<<"k\th\tf"<<endl;
	for (i=1;i<9;i++)
		cout<<i+1<<"\t"<<h2[i]<<" "<<f2[i]<<endl;
	system("rm p.fill.in.pair.R");
	tempFileName="rm "+ tempFileName;
	system(tempFileName.c_str());
	tempFileName1="rm "+ tempFileName1;
	system(tempFileName1.c_str());
	i=(int)(clock()-tStart)/CLOCKS_PER_SEC;
	cout<<"Time taken: "<<i/3600<<" hours "<<i%3600/60<<" minutes "<<i%60<<" seconds"<<endl;
}

void CalculatePairwiseInteractionsOnGivenPairsOfMarkers()
{
	//get single SNP info
	double t, n, p0, p1, N0r, N1r, a[9], b[9], p1value[maxMatrixRow],w1[maxMatrixRow];
	int i,j,k,df1[maxMatrixRow];
	bool myflag;
	
	for (i=0;i<row;i++)
		{
			a[0]=a[1]=a[2]=b[0]=b[1]=b[2]=0;
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
			w1[i]=0; df1[i]=2;
			for (j=0;j<3;j++)
				if (a[j])
					{
						p0=a[j]*N0r;
						p1=b[j]*N1r;
						t=log((1-p0)*p1/(p0*(1-p1)))/sqrt(1.0/(1.0*a[j])+1.0/(1.0*b[j])+1.0/(1.0*(N0+n-a[j]))+1.0/(1.0*(N1+n-b[j])));
						w1[i]+=t*t;
					}
					else
						--df1[i];
			w1[i]=w1[i]*h1[df1[i]];
			p1value[i]=chisqr(df1[i],w1[i]);
		}

	string inputFileName,outputFileName,tempFileName,tempFileName1,s,s1,s2;
	cout<<"Input file cotains list of marker pairs.\nfile name = ";
	cin>>inputFileName;
	cout<<"Result output file name:\nfile name = ";
	cin>>outputFileName;
	clock_t tStart = clock();
	int nmarker,snp1,snp2;
	ifstream fin(inputFileName.c_str());
	nmarker=0;

	tempFileName=GetTempFileName();
	ofstream fout(tempFileName.c_str());
	double w;
	int df;
	getline(fin,s);//remove the header
	while (fin>>s1>>s2)	
		{
				++nmarker;
				for (k=0;k<row;k++)
					{
						if (s1==SNPname[k])
							snp1=k;
						if (s2==SNPname[k])
							snp2=k;
					}
				//snp1=snp1-1;
				//snp2=snp2-1;
				getline(fin,s);
				for (k=0;k<9;k++)
					a[k]=b[k]=0;
				for (k=0;k<N0reduced;k++)
					{
						a[0]+=__builtin_popcount(matrix[snp1][0][k] & matrix[snp2][0][k]);
						a[1]+=__builtin_popcount(matrix[snp1][0][k] & matrix[snp2][1][k]);
						a[2]+=__builtin_popcount(matrix[snp1][0][k] & matrix[snp2][2][k]);
						a[3]+=__builtin_popcount(matrix[snp1][1][k] & matrix[snp2][0][k]);
						a[4]+=__builtin_popcount(matrix[snp1][1][k] & matrix[snp2][1][k]);
						a[5]+=__builtin_popcount(matrix[snp1][1][k] & matrix[snp2][2][k]);
						a[6]+=__builtin_popcount(matrix[snp1][2][k] & matrix[snp2][0][k]);
						a[7]+=__builtin_popcount(matrix[snp1][2][k] & matrix[snp2][1][k]);
						a[8]+=__builtin_popcount(matrix[snp1][2][k] & matrix[snp2][2][k]);
					}
				for (k=N0reduced;k<N0reduced+N1reduced;k++)
					{
						b[0]+=__builtin_popcount(matrix[snp1][0][k] & matrix[snp2][0][k]);
						b[1]+=__builtin_popcount(matrix[snp1][0][k] & matrix[snp2][1][k]);
						b[2]+=__builtin_popcount(matrix[snp1][0][k] & matrix[snp2][2][k]);
						b[3]+=__builtin_popcount(matrix[snp1][1][k] & matrix[snp2][0][k]);
						b[4]+=__builtin_popcount(matrix[snp1][1][k] & matrix[snp2][1][k]);
						b[5]+=__builtin_popcount(matrix[snp1][1][k] & matrix[snp2][2][k]);
						b[6]+=__builtin_popcount(matrix[snp1][2][k] & matrix[snp2][0][k]);
						b[7]+=__builtin_popcount(matrix[snp1][2][k] & matrix[snp2][1][k]);
						b[8]+=__builtin_popcount(matrix[snp1][2][k] & matrix[snp2][2][k]);
						
					}
				myflag=false; //check if we need to increment 0.5
				for (k=0;k<9;k++)
					if (((a[k]==0) && (b[k]!=0)) || ((b[k]==0) && (a[k]!=0)))
						{
							myflag=true;
							break;
						}
				n=0.0;
				if (myflag)
					{
						for (k=0;k<9;k++)
							if (!((a[k]==0) && (b[k]==0)))
								{
									a[k]+=0.5;
									b[k]+=0.5;
									n+=0.5;
								}
					}
				N0r=1/(1.0*N0+n);
				N1r=1/(1.0*N1+n);
				w=0; df=8;
				for (k=0;k<9;k++)
					if (a[k])
						{
							p0=a[k]*N0r;
							p1=b[k]*N1r;
							t=log((1-p0)*p1/(p0*(1-p1)))/sqrt(1.0/(1.0*a[k])+1.0/(1.0*b[k])+1.0/(1.0*(N0+n-a[k]))+1.0/(1.0*(N1+n-b[k])));
							w+=t*t;
						}
						else
							--df;
				w=w*h2[df];
				//snp1, snp2, w, k, p-value, w1, k1, p1, w2,  k2,  p2
				fout<<SNPname[snp1]<<"\t"<<SNPname[snp2]<<"\t"<<w<<"\t"<<df+1<<"\t"<<0<<"\t"<<w1[snp1]<<"\t"<<df1[snp1]<<"\t"<<p1value[snp1]<<"\t"<<w1[snp2]<<"\t"<<df1[snp2]<<"\t"<<p1value[snp2]<<endl;
			}
	fin.close();
	fout.close();
	if (nmarker<1)
		{
			tempFileName="rm "+ tempFileName;
			system(tempFileName.c_str());
			cout<<"Bad input file"<<endl;
			return;
		}
	tempFileName1=GetTempFileName();
	CreatePFillInRProgram(tempFileName,tempFileName1);
	s="chmod a+x p.fill.in.pair.R";
	system(s.c_str());
	system("R --slave < p.fill.in.pair.R");
	s="sort -k5g "+tempFileName1+" | sed '1s/^/SNP1\tSNP2\tW\tk\tp.pair\tW.snp1\tk1\tp.snp1\tW.snp2\tk2\tp.snp2\n/' > "+outputFileName;
	system(s.c_str());
	cout<<"k\th\tf"<<endl;
	for (i=1;i<9;i++)
		cout<<i+1<<"\t"<<h2[i]<<" "<<f2[i]<<endl;
	system("rm p.fill.in.pair.R");
	tempFileName="rm "+ tempFileName;
	system(tempFileName.c_str());
	tempFileName1="rm "+ tempFileName1;
	system(tempFileName1.c_str());
	i=(int)(clock()-tStart)/CLOCKS_PER_SEC;
	cout<<"Time taken: "<<i/3600<<" hours "<<i%3600/60<<" minutes "<<i%60<<" seconds"<<endl;
}
