#include<string>
#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<time.h>
using namespace std;
const int maxMatrixRow=10001;//maximum number of genotype
const int maxMatrixColumn=1001;//number of phenotype

int row=0,column,N0,N1,N0reduced,N1reduced;

short y[maxMatrixColumn];
//boost Matrix
unsigned int matrix[maxMatrixRow][3][maxMatrixColumn/32+3];

short matrixOri[maxMatrixRow][maxMatrixColumn];

string SNPname[maxMatrixRow];

double h1[3],f1[3],h2[9],f2[9];

void PermuteY(short arr[],int n)
{
	short t,l;
	for (int i=n-1;i>0;i--)
		{
			l=rand()%(i+1);
			t=arr[i];
			arr[i]=arr[l];
			arr[l]=t;
		}
}

void calculateOrder2hf(int B, int p)
{
	int i,j,g1,g2;
	double Xbar[B][9],V[B][9],E[9],VarX[9];//modified
	double Xtemp[p*(p-1)];
	int nk[B][9],ktemp[p*(p-1)];//,nktemp[9];
	double n,N0r,N1r,p0,p1,t;
	short y_permuted[maxMatrixColumn];
	for (i=0;i<column;i++)
		y_permuted[i]=y[i];
	srand(time(NULL));
	double a[9],b[9];
	bool myflag;
	int selectedSNP;
	
	int snpSelect[p];
	bool selected[row];
	for (i=0;i<row;i++)
		{
			snpSelect[i]=i;
			selected[i]=false;
		}
	
	
	for (int bth=0;bth<B;bth++) //permute B times for Y
		{
			if (p!=row)
				{
					for (i=0;i<p;i++)
						{
							selectedSNP=rand()%row;
							while ((selected[selectedSNP]))
								{
									selectedSNP=rand()%row;
								}
							selected[selectedSNP]=true;
							snpSelect[i]=selectedSNP;
						}
					for (i=0;i<p;i++)
						selected[snpSelect[i]]=false;
				}
			PermuteY(y_permuted,column);
			i=0;//counter
			for (g1=0;g1<p;g1++)
				for (g2=g1+1;g2<p;g2++)
					{
						//calculate X2
						for (j=0;j<9;j++)
						a[j]=b[j]=0;
						for (j=0;j<column;j++)
							if (y_permuted[j])
								++b[matrixOri[snpSelect[g1]][j]*3+matrixOri[snpSelect[g2]][j]];
								else
									++a[matrixOri[snpSelect[g1]][j]*3+matrixOri[snpSelect[g2]][j]];
						myflag=false;
						for (j=0;j<9;j++)
							if (((a[j]==0) && (b[j]!=0)) || ((b[j]==0) && (a[j]!=0)))
								{
									myflag=true;
									break;
								}
						n=0.0;
						if (myflag)
							{
								for (j=0;j<9;j++)
									if (!((a[j]==0)&&(b[j]==0)))
										{
											a[j]+=0.5;
											b[j]+=0.5;
											n+=0.5;
										}
							}
						N0r=1/(1.0*N0+n);
						N1r=1/(1.0*N1+n);
						Xtemp[i]=0; ktemp[i]=8;
						for (j=0;j<9;j++)
							if (a[j])
								{
									p0=a[j]*N0r;
									p1=b[j]*N1r;
									t=log((1-p0)*p1/(p0*(1-p1)))/sqrt(1.0/(1.0*a[j])+1.0/(1.0*b[j])+1.0/(1.0*(N0+n-a[j]))+1.0/(1.0*(N1+n-b[j])));
									Xtemp[i]+=t*t;
								}
								else
									--ktemp[i];
						++i;
					}//end g2 loop
			j=i;
			for (i=1;i<9;i++)
				{
					Xbar[bth][i]=0;
					V[bth][i]=0;
					nk[bth][i]=0;
				}
			for (i=0;i<j;i++)
				{
					Xbar[bth][ktemp[i]]+=Xtemp[i];
					++nk[bth][ktemp[i]];
				}

			for (i=1;i<9;i++)
				if (nk[bth][i])
					Xbar[bth][i]=Xbar[bth][i]/nk[bth][i];
			for (i=0;i<j;i++)
				V[bth][ktemp[i]]+=(Xtemp[i]-Xbar[bth][ktemp[i]])*(Xtemp[i]-Xbar[bth][ktemp[i]]);
			for (i=1;i<9;i++)
				if (nk[bth][i]>1)
					V[bth][i]=V[bth][i]/((nk[bth][i]-1)*1.0);
		}//end bth loop
	h2[1]=0.5;h2[2]=0.667;h2[3]=0.75;h2[4]=0.8;h2[5]=0.833;h2[6]=0.857;h2[7]=0.875;h2[8]=0.889;
	f2[1]=1;f2[2]=2;f2[3]=3;f2[4]=4;f2[5]=5;f2[6]=6;f2[7]=7;f2[8]=8;
	for (i=1;i<9;i++)
		{
			E[i]=0;VarX[i]=0;
			for (j=0;j<B;j++)
				{
					E[i]+=Xbar[j][i];
					VarX[i]+=V[j][i];
				}
			E[i]=E[i]/(1.0*B);
			VarX[i]=VarX[i]/(1.0*B);
			if (VarX[i])
				{
					h2[i]=2*E[i]/VarX[i];
					f2[i]=2*E[i]*E[i]/VarX[i];
				}
		}
}


void calculateOrder1hf(int B, int p)
{
	int i,j,g1,g2;
	double X[B];
	double Xbar[B][3],V[B][3],E[3],VarX[3];//modified
	double Xtemp[row];
	int nk[B][3],ktemp[row];
	double n,N0r,N1r,p0,p1,t;
	short y_permuted[maxMatrixColumn];
	for (i=0;i<column;i++)
		y_permuted[i]=y[i];
	srand(time(NULL));
	double a[3],b[3];
	bool myflag;
	int selectedSNP;
	
	int snpSelect[p];
	bool selected[row];
	for (i=0;i<row;i++)
		{
			snpSelect[i]=i;
			selected[i]=false;
		}
	
	for (int bth=0;bth<B;bth++) //permute B times for Y
		{
			if (p!=row)
				{
					for (i=0;i<p;i++)
						{
							selectedSNP=rand()%row;
							while ((selected[selectedSNP]))
								{
									selectedSNP=rand()%row;
								}
							selected[selectedSNP]=true;
							snpSelect[i]=selectedSNP;
						}
					for (i=0;i<p;i++)
						selected[snpSelect[i]]=false;
				}
			PermuteY(y_permuted,column);
			
			for (i=0;i<p;i++) //random generate Row-number of record pairs
				{
					for (j=0;j<3;j++)
						a[j]=b[j]=0;
					for (j=0;j<column;j++)
						if (y_permuted[j])
							++b[matrixOri[snpSelect[i]][j]];
							else
								++a[matrixOri[snpSelect[i]][j]];
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
					Xtemp[i]=0; ktemp[i]=2;
					for (j=0;j<3;j++)
						if (a[j])
							{
								p0=a[j]*N0r;
								p1=b[j]*N1r;
								t=log((1-p0)*p1/(p0*(1-p1)))/sqrt(1.0/(1.0*a[j])+1.0/(1.0*b[j])+1.0/(1.0*(N0+n-a[j]))+1.0/(1.0*(N1+n-b[j])));
								Xtemp[i]+=t*t;
							}
							else
								--ktemp[i];
				}
			for (i=1;i<3;i++)
				{
					Xbar[bth][i]=0;
					V[bth][i]=0;
					nk[bth][i]=0;
				}
			for (i=0;i<p;i++)
				{
					Xbar[bth][ktemp[i]]+=Xtemp[i];
					++nk[bth][ktemp[i]];
				}
			for (i=1;i<3;i++)
				if (nk[bth][i])
					Xbar[bth][i]=Xbar[bth][i]/nk[bth][i];
			for (i=0;i<p;i++)
				V[bth][ktemp[i]]+=(Xtemp[i]-Xbar[bth][ktemp[i]])*(Xtemp[i]-Xbar[bth][ktemp[i]]);
			for (i=1;i<3;i++)
				if (nk[bth][i]>1)
					V[bth][i]=V[bth][i]/((nk[bth][i]-1)*1.0);
		}
	h1[1]=0.5;h1[2]=0.667;
	f1[1]=1;f1[2]=2;
	for (i=1;i<3;i++)
		{
			E[i]=0;VarX[i]=0;
			for (j=0;j<B;j++)
				{
					E[i]+=Xbar[j][i];
					VarX[i]+=V[j][i];
				}
			E[i]=E[i]/(1.0*B);
			VarX[i]=VarX[i]/(1.0*B);
			if (VarX[i])
				{
					h1[i]=2*E[i]/VarX[i];
					f1[i]=2*E[i]*E[i]/VarX[i];
				}
		}
}


void getMatrixColumn(string inputFileName)
{
 const char *chtemp=inputFileName.c_str(); 
 ifstream fin(chtemp);
 string s;
 getline(fin,s);
 column=0;
 s=s+" ";
 int i=0;
 while ((i<s.length()) && ((s[i]=='\t') || (s[i]==' ')))
       {
        ++i;
       }
 while ((i<s.length()) && ((s[i]!='\t') || (s[i]!=' ')))
       {
        ++i;
        while ((i<s.length()) && ((s[i]=='\t') || (s[i]==' ')))
              {
               ++column;
               ++i;
              }
       }
 fin.close();
}

void  getNumberSNP(string inputFileName)
{
 const char *chtemp=inputFileName.c_str(); 
 ifstream fin(chtemp);
 string s;
 getline(fin,s);
 row=0;
 s=s+" ";
 int i=0;
 while ((i<s.length()) && ((s[i]=='\t') || (s[i]==' ')))
       {
        ++i;
       }
 while ((i<s.length()) && ((s[i]!='\t') || (s[i]!=' ')))
       {
        ++i;
        while ((i<s.length()) && ((s[i]=='\t') || (s[i]==' ')))
              {
               ++row;
               ++i;
              }
       }
 fin.close();
}

void  getNumberPheno(string inputFileName)
{
	const char *chtemp=inputFileName.c_str(); 
	ifstream fin(chtemp);
	string s;
	getline(fin,s);
	getline(fin,s);
	column=0;
	while ((s!="\n") && (s!=""))
		{
			++column;
			getline(fin,s);
		}
	fin.close(); 
}


bool LoadData(){
	string inputFileName1,inputFileName2;
	cout<<"Input source genotype file name:\nfile name = ";
	cin>>inputFileName1;
	const char* chtempInput1=inputFileName1.c_str();
	ifstream fin(chtempInput1);
	if (!fin.is_open())
       {
        cout<<"Error: File not exist / File Corrupted / Incorrect Format"<<endl;
        fin.close();
        return false;
       }
	string s;
	int i,j;
	getNumberSNP(inputFileName1);
	getNumberPheno(inputFileName1);
	for (i=0;i<row;i++)
		fin>>SNPname[i];
	for (i=0;i<column;i++)
		{
			for (j=0;j<row;j++)
				{
					fin>>matrixOri[j][i];
					if ((matrixOri[j][i]>2)||(matrixOri[j][i]<0))
						{
							fin.close();
							cout<<"Error: data out of allowed range {0,1,2}"<<endl;
							row=0;
							return false;	
						}
				}
		}
	fin.close();
	if (row<2)
		{
			cout<<"Error: Input data only have one SNP info, the algorithm is unable to evaluate";
			row=0;
			return false;
		}
	//Reading case and control
	N0=N1=0;
	cout<<"Input source phenotype file name:\nfile name = ";
	cin>>inputFileName2;
	const char* chtempInput2=inputFileName2.c_str();
	ifstream fin2(chtempInput2);
	if (!fin2.is_open())
       {
        cout<<"Error: File not exist / File Corrupted / Incorrect Format"<<endl;
        fin2.close();
        return false;
       }
	for (i=0;i<column;i++)
		{
			fin2>>s>>y[i];
			if (y[i])
				++N1;
				else
					++N0;
		}
	fin2.close();
	//Initializing Boost Matrix
	N0reduced=N1reduced=0;
	int t,t0=0,t1=0;
	for (i=0;i<column;i++)
		{
			if (y[i]==0)//case Y=0
				{
					if (t0%32==0)
						{
							++N0reduced;++t0;
							for (j=0;j<row;j++)
								{
									t=matrixOri[j][i];
									if (t==2)
										{
											matrix[j][2][N0reduced-1]=1;
											matrix[j][1][N0reduced-1]=matrix[j][0][N0reduced-1]=0;
										}
										else if (t==1)
											{
												matrix[j][1][N0reduced-1]=1;
												matrix[j][2][N0reduced-1]=matrix[j][0][N0reduced-1]=0;
											}
											else
												{
													matrix[j][0][N0reduced-1]=1;
													matrix[j][1][N0reduced-1]=matrix[j][2][N0reduced-1]=0;
												}
								}
						}
						else
							{
								++t0;
								for (j=0;j<row;j++)
									{
										t=matrixOri[j][i];
										if (t==2)
											{
												matrix[j][2][N0reduced-1]=matrix[j][2][N0reduced-1]*2+1;
												matrix[j][1][N0reduced-1]=matrix[j][1][N0reduced-1]*2;
												matrix[j][0][N0reduced-1]=matrix[j][0][N0reduced-1]*2;
											}
											else if (t==1)
												{
													matrix[j][2][N0reduced-1]=matrix[j][2][N0reduced-1]*2;
													matrix[j][1][N0reduced-1]=matrix[j][1][N0reduced-1]*2+1;
													matrix[j][0][N0reduced-1]=matrix[j][0][N0reduced-1]*2;
												}
												else
													{
														matrix[j][2][N0reduced-1]=matrix[j][2][N0reduced-1]*2;
														matrix[j][1][N0reduced-1]=matrix[j][1][N0reduced-1]*2;
														matrix[j][0][N0reduced-1]=matrix[j][0][N0reduced-1]*2+1;
													}
									}
							}
				}
		}
		
	for (i=0;i<column;i++)
		{
			if (y[i]==1)//case Y=1
				{
					if (t1%32==0)
						{
							++N1reduced;++t1;
							for (j=0;j<row;j++)
								{
									t=matrixOri[j][i];
									if (t==2)
										{
											matrix[j][2][N0reduced+N1reduced-1]=1;
											matrix[j][1][N0reduced+N1reduced-1]=matrix[j][0][N0reduced+N1reduced-1]=0;
										}
										else if (t==1)
											{
												matrix[j][1][N0reduced+N1reduced-1]=1;
												matrix[j][2][N0reduced+N1reduced-1]=matrix[j][0][N0reduced+N1reduced-1]=0;
											}
											else
												{
													matrix[j][0][N0reduced+N1reduced-1]=1;
													matrix[j][1][N0reduced+N1reduced-1]=matrix[j][2][N0reduced+N1reduced-1]=0;
												}
								}
						}
						else
							{
								++t1;
								for (j=0;j<row;j++)
									{
										t=matrixOri[j][i];
										if (t==2)
											{
												matrix[j][2][N0reduced+N1reduced-1]=matrix[j][2][N0reduced+N1reduced-1]*2+1;
												matrix[j][1][N0reduced+N1reduced-1]=matrix[j][1][N0reduced+N1reduced-1]*2;
												matrix[j][0][N0reduced+N1reduced-1]=matrix[j][0][N0reduced+N1reduced-1]*2;
											}
											else if (t==1)
												{
													matrix[j][2][N0reduced+N1reduced-1]=matrix[j][2][N0reduced+N1reduced-1]*2;
													matrix[j][1][N0reduced+N1reduced-1]=matrix[j][1][N0reduced+N1reduced-1]*2+1;
													matrix[j][0][N0reduced+N1reduced-1]=matrix[j][0][N0reduced+N1reduced-1]*2;
												}
												else
													{
														matrix[j][2][N0reduced+N1reduced-1]=matrix[j][2][N0reduced+N1reduced-1]*2;
														matrix[j][1][N0reduced+N1reduced-1]=matrix[j][1][N0reduced+N1reduced-1]*2;
														matrix[j][0][N0reduced+N1reduced-1]=matrix[j][0][N0reduced+N1reduced-1]*2+1;
													}
									}
							}
				}
		}
	cout<<"Calculating h & f ..."<<endl;
	char ch;
	int B=200,p=50;
	cout<<"Would you like to change default value for hf evaluation? The default value are trial times B = 200, number of SNPs selected p = 50 (y/n)";
	cin>>ch;
	if ((ch=='Y') || (ch=='y'))
		{
			cout<<"B = ";
			cin>>B;
			cout<<"p = ";
			cin>>p;
		}
	if (p>row)
		{
			cout<<"p value exceeds the maximum of SNPs from input file, p value is set to be "<<row<<" instead"<<endl;
			p=row;
		}
	if (p<2)
		{
			cout<<"p value is too small for evalution"<<endl;
			row=0;
			return false;
		}
	calculateOrder1hf(B,p);
	calculateOrder2hf(B,p);
 	return true;
}

