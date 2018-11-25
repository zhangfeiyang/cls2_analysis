{
	const int N = 8;
	ifstream fin("data2");
	ofstream fout("stat_uncertainty");
	double a,b,c,d,e;
	for(int i=0;i<N;i++){
		fin>>a>>b>>c>>d>>e;
		fout<<c/b<<"\n";
	}

}
