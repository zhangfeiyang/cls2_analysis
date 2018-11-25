{
	const int N = 8;
	ifstream fin0("data1");
	ifstream fin2("data2");
	ofstream fout("compton_uncertainty");
	double a,b,c,d,e,mean0,mean2;
	for(int i=0;i<N;i++){
		fin0>>a>>mean0>>a>>a>>a;
		fin2>>a>>mean2>>a>>a>>a;
		fout<<(mean2-mean0)/mean0<<"\n";
	}

}
