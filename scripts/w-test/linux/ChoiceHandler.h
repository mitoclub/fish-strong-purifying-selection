#include<iostream>
#include<string>
#include<fstream>
#include<cmath>
#include"LoadData.h"
#include"MainEffect.h"
#include"Pairwise.h"
using namespace std;

string GetCommand(string s1)
{
 string s;
 cout<<"Enter your choice (0-"<<s1<<"): ";
 getline(cin,s);
 if (s=="")
    getline(cin,s);
 while (((s.length()!=1) && (s.length()!=s1.length()))  || (s<"0"))
       {
        cout<<"Invalid input, re-enter again (0-"<<s1<<"): ";
        getline(cin,s);
       }
 return s;
}


void ChoiceHandle(string choice, int argc, char *argv[])
{
 string secondChoice;
 if (choice=="1")
	{
		if (LoadData())//LoadData.h
			cout<<"Data Loaded successful"<<endl;
			else
				cout<<"Failed Loading Data"<<endl;
		cin.ignore(1024,'\n');
		cout<<"Press Enter to continue...";
		cin.get();
	}
 if (choice=="2")
 	{
 		if (row>1)
			CalculateMainEffect();//MainEffect.h
			else
				cout<<"No valid input data found"<<endl;
		cin.ignore(1024,'\n');
		cout<<"Press Enter to continue...";
		cin.get();
 	}
 if (choice=="3")
 	{
 		if (row>1)
		 	{
		 		while (true)
    				{
            			system("clear");
    					cout<<"Calculate Pairwise Interactions"<<endl;
 						cout<<"--------------------------------"<<endl<<endl;
						cout<<"0. Go Back"<<endl;
						cout<<"1. Calculates Pairwise Effect Automatically Using All The Input Data"<<endl;
    					cout<<"2. Calculate Pairwise Interactions Exhaustively From Main Effect Result (requires result from Main Effect)"<<endl;
            			cout<<"3. Calculate Pairwise Interactions On Given Pairs Of Markers"<<endl;
						secondChoice=GetCommand("3");
						if (secondChoice=="0")
    						return;
						system("clear");
						if (secondChoice=="1")
							{
								CalculatePairwiseInteractionsAll();//Pairwise.h
								cout<<"Calculation completed"<<endl;
								cin.ignore(1024,'\n');
								cout<<"Press Enter to continue...";
								cin.get();
							}
						if (secondChoice=="2")
							{
								CalculatePairwiseInteractionsExhaustively();//Pairwise.h
								cout<<"Calculation completed"<<endl;
								cin.ignore(1024,'\n');
								cout<<"Press Enter to continue...";
								cin.get();
							}
						if (secondChoice=="3")
							{
								CalculatePairwiseInteractionsOnGivenPairsOfMarkers();//Pairwise.h
								cout<<"Calculation completed"<<endl;
								cin.ignore(1024,'\n');
								cout<<"Press Enter to continue...";
								cin.get();
							}
					}
           }
           	else
           		{
           			cout<<"No valid input data found"<<endl;
					cin.ignore(1024,'\n');
					cout<<"Press Enter to continue...";
					cin.get();
           		}
    }
}
