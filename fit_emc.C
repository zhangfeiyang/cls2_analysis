#include "rootheader.h"
double myfun(double *xx, double *par)
{
    double E = xx[0];
    double alpha = par[0];
    double mu = par[1];
    double sigma = par[2];
    double beta = par[3];
    double lambda = par[4];
    double gamma = par[5];
    double pi = 3.141592654;

    double gaus = alpha/sigma/sqrt(2*pi)*exp(-0.5*(E-mu)*(E-mu)/sigma/sigma);
    double exp_c1 = (beta)*lambda*exp((sigma*sigma*lambda*lambda+2*lambda*E)/2)/(exp(lambda*mu)-1);
    double exp1 = TMath::Erf((mu-E-sigma*sigma*lambda)/(TMath::Sqrt(2)*sigma))-TMath::Erf((-E-sigma*sigma*lambda)/sqrt(2)/sigma);
    double constant = (gamma)/mu*(TMath::Erf((mu-E)/sqrt(2)/sigma)-TMath::Erf(-E/sqrt(2)/sigma));

    return (gaus+exp_c1*exp1+constant)/(alpha+beta+gamma);

}

double myfun2(double *xx, double *par)
{
    double E = xx[0];
    double alpha = par[0];
    double mu = par[1];
    double sigma = par[2];
    double beta = par[3];
    double lambda = par[4];
    double gamma = par[5];
    double pi = 3.141592654;

    double gaus = alpha/sigma/sqrt(2*pi)*exp(-0.5*(E-mu)*(E-mu)/sigma/sigma);
    double exp_c1 = (beta)*lambda*exp((sigma*sigma*lambda*lambda+2*lambda*E)/2)/(exp(lambda*mu)-1);
    double exp1 = TMath::Erf((mu-E-sigma*sigma*lambda)/(TMath::Sqrt(2)*sigma))-TMath::Erf((-E-sigma*sigma*lambda)/sqrt(2)/sigma);
    double constant = (gamma)/mu*(TMath::Erf((mu-E)/sqrt(2)/sigma)-TMath::Erf(-E/sqrt(2)/sigma));

    return (gaus)/(alpha+beta+gamma);

}

double myfun3(double *xx, double *par)
{
    double E = xx[0];
    double alpha = par[0];
    double mu = par[1];
    double sigma = par[2];
    double beta = par[3];
    double lambda = par[4];
    double gamma = par[5];
    double pi = 3.141592654;
        double gaus = alpha/sigma/sqrt(2*pi)*exp(-0.5*(E-mu)*(E-mu)/sigma/sigma);    double exp_c1 = (beta)*lambda*exp((sigma*sigma*lambda*lambda+2*lambda*E)/2)/(exp(lambda*mu)-1);    double exp1 = TMath::Erf((mu-E-sigma*sigma*lambda)/(TMath::Sqrt(2)*sigma))-TMath::Erf((-E-sigma*sigma*lambda)/sqrt(2)/sigma);    double constant = (gamma)/mu*(TMath::Erf((mu-E)/sqrt(2)/sigma)-TMath::Erf(-E/sqrt(2)/sigma));
    
    return (exp_c1*exp1+constant)/(alpha+beta+gamma);

}
/*
double myfun4(double *xx, double *par)
{
    double E = xx[0];
    double c = par[0];
    double alpha = par[1];
    double mu = par[2];
    double sigma = par[3];
    double beta = par[4];
    double lambda = par[5];
    double gamma = par[6];
    double pi = 3.141592654;
    double gaus = alpha/sigma/sqrt(2*pi)*exp(-0.5*(E-mu)*(E-mu)/sigma/sigma);
    double exp_c1 = (beta)*lambda*exp((sigma*sigma*lambda*lambda+2*lambda*E)/2)/(exp(lambda*mu)-1);
    double exp1 = TMath::Erf((mu-E-sigma*sigma*lambda)/(TMath::Sqrt(2)*sigma))-TMath::Erf((-E-sigma*sigma*lambda)/sqrt(2)/sigma);
    double constant = (gamma)/mu*(TMath::Erf((mu-E)/sqrt(2)/sigma)-TMath::Erf(-E/sqrt(2)/sigma));
    return c*(exp_c1*exp1)/(alpha+beta+gamma);
}
*/
bool get_data(int argc,char **argv,double *pars,double *epars,double& chi2,double &ndf,string dirname,int index)
{

	string source = argv[1];
	string R = argv[2];	
	string Z = argv[3];	

	TCanvas *c1 = new TCanvas();
    gStyle->SetOptFit(1);
    gStyle->SetStatX(0.4);
    gStyle->SetStatY(0.9);
    gStyle->SetStatW(0.15);
    gStyle->SetStatH(0.15);
    TFile *file;
    TF1 *f = new TF1("f",myfun,0,20000,6);

	TChain *t = new TChain("evt");

//	t->Add("/junofs/production/public/users/zhangfy/non-uniform/offline_J17v1r1-Pre1/Examples/Tutorial/share/cls2/Co60/10000_10000/evt_*root");
	string filename = dirname+"evt_*root";
	t->Add(&filename[0]);
	
	int total_entries = t->GetEntries();	

	int max_entries;
	if(source=="Ge68" || source == "K40")
		max_entries = 400000;
	else
		max_entries = 40000;

	if(source=="Cs137" || source== "Mn54")
		max_entries = 82000;

	if(max_entries*(index+1) > total_entries) return false;
	
	cout << max_entries*(index+1) << "\t" << max_entries*index << "\n";
	t->Draw("totalPE>>h(200,0,0)","totalPE>0","",max_entries,max_entries*index);
	TH1F *h = (TH1F*)gDirectory->Get("h");

	int maxbin = h->GetMaximumBin();
	double maxbincenter = h->GetBinCenter(maxbin);

	t->Draw(Form("totalPE>>h(100,%i,%i)",int(maxbincenter-100),int(maxbincenter+100)),"","",max_entries,max_entries*index);	
	h = (TH1F*)gDirectory->Get("h");
	//h->Fit("gaus");
	
	
//	h->Fit("f");

/*  
    file = new TFile("1.022_AfterC.root","read");
  	t = (TTree*)file->Get("tt");
	TF1 *f = new TF1("f",myfun,0,20000,6);
	t->Draw("newPE>>h(200,1000,1600)","edep>1.0218");	
	TH1F *h = (TH1F*)gDirectory->Get("h");
*/	h->Fit("gaus","Q");
    TF1* fun = h->GetFunction("gaus");

    double C = fun->GetParameter(0);
    double mean = fun->GetParameter(1);
    double emean = fun->GetParError(2);
    double sigma = fun->GetParameter(2);

	//f->SetParameters(C*sqrt(2*3.14159)*sigma,0.99,mean,sigma,0.001,0.01,0.02);
    pars[0] = C*sqrt(2*3.14159)*sigma;
    pars[1] = mean;
    pars[2] = sigma;
    pars[3] = 0.01;
    pars[4] = 0.01;
    pars[5] = 0.01;

    f->SetParLimits(1,0,1000);
    f->SetParLimits(3,0,1000);
    f->SetParLimits(4,0,1);
    f->SetParLimits(5,0,1000);
	f->SetParameters(pars);

	f->SetParNames("#alpha","#mu","#sigma","#beta","#lambda","#gamma");
	t->Draw(Form("totalPE>>h(%i,%i,%i)",int(5*sigma),int(mean-10*sigma),int(mean-10*sigma)+int(5*sigma)*3),"","",max_entries,max_entries*index);	
	//t->Draw(Form("totalPE>>h(%i,%i,%i)",int(5*sigma),int(mean-5*sigma),int(mean-5*sigma)+int(5*sigma)*2),"","",max_entries*(index+1),max_entries*index);	
	//t->Draw(Form("totalPE>>h(%i,%i,%i)",int(5*sigma),int(mean-6*sigma),int(mean-6*sigma)+int(5*sigma)*2),"","",max_entries);	
	//t->Draw(Form("totalPE>>h(%i,%i,%i)",int(5*sigma),int(mean-6*sigma),int(mean-6*sigma)+int(5*sigma)*2),"","");	
	h = (TH1F*)gDirectory->Get("h");
	
	h->Fit(f,"M","");
	f->GetParameters(pars);

	if(source=="Ge68"){
		h->Fit(f,"M","",pars[2]-pars[3]*10,pars[2]+pars[3]);
		h->Fit(f,"M","",pars[2]-pars[3]*10,pars[2]+pars[3]);
		h->Fit(f,"M","",pars[2]-pars[3]*10,pars[2]+pars[3]);
	}
	else{
		h->Fit(f,"M","",int(mean-10*sigma),int(mean-5*sigma)+int(5*sigma)*2);
		h->Fit(f,"M","",int(mean-10*sigma),int(mean-5*sigma)+int(5*sigma)*2);
		h->Fit(f,"M","",int(mean-10*sigma),int(mean-5*sigma)+int(5*sigma)*2);
		h->Fit(f,"M","",int(mean-10*sigma),int(mean-5*sigma)+int(5*sigma)*2);
		h->Fit(f,"M","",int(mean-10*sigma),int(mean-5*sigma)+int(5*sigma)*2);
		h->Fit(f,"M","",int(mean-10*sigma),int(mean-5*sigma)+int(5*sigma)*2);
		h->Fit(f,"M","",int(mean-10*sigma),int(mean-5*sigma)+int(5*sigma)*2);
	}

	f->GetParameters(pars);
	double *tmps = f->GetParErrors();
	memcpy(epars,tmps,6*sizeof(double));

	ndf = f->GetNDF();
	chi2 = f->GetChisquare();

    TF1 *f2 = new TF1("f",myfun2,0,20000,7);
    TF1 *f3 = new TF1("f",myfun3,0,20000,7);
	f2->SetParameters(pars);
	f3->SetParameters(pars);

	//f2->SetParameter(0,pars[0]);	
	//f2->SetParameter(1,pars[1]);	
	//f2->SetParameter(2,pars[2]);	
	//f2->SetParameter(3,pars[3]);	
	//f2->SetParameter(4,pars[4]);	
	//f2->SetParameter(5,pars[5]);	
	//f2->SetParameter(6,pars[6]);	

	//f3->SetParameter(0,pars[0]);	
	//f3->SetParameter(1,pars[1]);	
	//f3->SetParameter(2,pars[2]);	
	//f3->SetParameter(3,pars[3]);	
	//f3->SetParameter(4,pars[4]);	
	//f3->SetParameter(5,pars[5]);	
	//f3->SetParameter(6,pars[6]);	

	f2->SetLineColor(kBlue);
	f3->SetLineColor(kBlack);
	f2->Draw("same");
	f3->Draw("same");

	filename = dirname+"result_emc.C";
	c1->SaveAs(&filename[0]);
	filename = dirname+"result_emc.png";
	c1->SaveAs(&filename[0]);
	
	for(int i=0;i<6;i++){
		cout << pars[i] <<"\t" << epars[i] <<"\t";
	}

	return true;

}
int main(int argc,char **argv){

	double *pars = new double[6];
	double *epars = new double[6];
	double chi2;
	double ndf;
	string source = argv[1];
	string R = argv[2];	
	string Z = argv[3];	

	string dirname = "/junofs/production/public/users/zhangfy/non-uniform/offline_J17v1r1-Pre1/Examples/Tutorial/share/cls2/"+source+"/"+R+"_"+Z+"/";
	string filename = dirname+"result3_emc";
	ofstream fout(&filename[0]);
	int index = 0;

	while(get_data(argc,argv,pars,epars,chi2,ndf,dirname,index)){
		
		fout<< chi2 <<"\t"<< ndf <<"\t";

		for(int i=0;i<6;i++){
			fout << pars[i] <<"\t" << epars[i] <<"\t";
		}

		fout << "\n";
		index++;
	}

	return 0;
}
