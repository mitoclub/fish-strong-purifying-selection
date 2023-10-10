#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<cmath>
#include<stdlib.h>
/* 
 * Important: Global variable defined in LoadData.h
 */
#include"ChoiceHandler.h"
using namespace std;


int main(int argc, char *argv[]){
	string choice;
	while (true)
    	{
    		system("clear");
    		cout<<"System"<<endl;
    		cout<<"--------------------------------"<<endl;
    		cout<<"0. Exit"<<endl;
    		cout<<"1. Load Input Dataset"<<endl;
    		cout<<"2. Calculate Main Effect"<<endl;
        	cout<<"3. Calculate Pairwise Interactions"<<endl;
    		choice=GetCommand("3");
    		if (choice=="0")
        		{
        			cout<<"Thank you. Bye!"<<endl;
					break;	
        		}
    		system("clear");
    		ChoiceHandle(choice, argc, argv);//ChoiceHandler.h
       }
 return 0;
}
