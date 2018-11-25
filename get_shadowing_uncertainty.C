{
	const int N = 8;
	ifstream fin0("data0");
	ifstream fin2("data1");
	ofstream fout("shadowing_uncertainty");
	double a,b,c,d,e,mean0,mean2,error0,error2;
	for(int i=0;i<N;i++){
		fin0>>a>>mean0>>error0>>a>>a;
		fin2>>a>>mean2>>error2>>a>>a;
		fout<<(mean2-mean0)/mean0<<"\t"<<sqrt(error0**2 + error2**2)/mean0<<"\n";
	}

}
