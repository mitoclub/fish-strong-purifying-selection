#include<string>
string IntToString(int number)
{
 string st="";
 if (number==0)
    return "0";
 while (number!=0)
       {
        char ch=(number%10)+'0';
        number=number/10;
        st=ch+st;
       }
 return st;
}

string GetTempFileName()
{
	int i=1;
	string s;
	while (true)
		{
			s="temp"+IntToString(i);
			ifstream f(s.c_str());
			if (f.good())
				++i;
				else
					return s;
		}
}

void CreatePFillInRProgram(string s1,string s2)
{
	ofstream fout("p.fill.in.pair.R");
	fout<<"p.fill.in<-function(file.read){"<<endl;
	fout<<"data.test<-read.table(file.read)"<<endl;
	fout<<"p.value<-pchisq(data.test[,3],df=data.test[,4],lower.tail=F)  #calculate p-value"<<endl;
	fout<<"data.test[,5]<-p.value"<<endl;
	fout<<"file.out<-paste0(\""<<s2<<"\")"<<endl;
	fout<<"write.table(data.test,file=file.out,quote=F,col.name=F,row.name=F)"<<endl;
	fout<<"}"<<endl;
	fout<<"file.e3=c(\""<<s1<<"\")"<<endl;
	fout<<"p.fill.in(file.e3)"<<endl;
	fout.close();
}
