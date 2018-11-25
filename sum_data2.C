{
	
	ifstream fin("log0");
	string sources[8] = {"Ge68","Cs137","Mn54","Co60","K40","n-H","n-C","prompt"};
	double energy[8] = {0.5109989,0.6617,0.8348,1.2528685,1.4608,2.223,4.95,6.13};
	
	double tmp;
	for(int i=0;i<8;i++){
		string filename;
		cout<< energy[i]<<"\t";

		if(i<5)
			filename = "/junofs/production/public/users/zhangfy/non-uniform/offline_J17v1r1-Pre1/Examples/Tutorial/share/cls2/"+sources[i]+"/0_0/result_emc";
		else{
			if(i==5) filename = "/junofs/production/public/users/zhangfy/non-uniform/offline_J17v1r1-Pre1/Examples/Tutorial/share/cls2/n-H/0_0/result_nH"; 
			if(i==6) filename = "/junofs/production/public/users/zhangfy/non-uniform/offline_J17v1r1-Pre1/Examples/Tutorial/share/cls2/n-H/0_0/result_nC"; 
			if(i==7) filename = "/junofs/production/public/users/zhangfy/non-uniform/offline_J17v1r1-Pre1/Examples/Tutorial/share/cls2/n-H/0_0/result_prompt"; 
		}
		ifstream fin2(&filename[0]);
		for(int j=0;j<14;j++){
			fin2>>tmp;
			if(j>5 && j<10)	
				cout<<tmp<<"\t";
		}
		cout << "\n";
		
	}

}
