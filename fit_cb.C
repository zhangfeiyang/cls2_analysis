#include "rootheader.h"
/*
double myfun(double *xx, double *par)
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

    //return c*(gaus+exp_c1*exp1+constant)/(alpha+beta+gamma);
    return c*(gaus+exp_c1*exp1+constant);
}*/

double myfun(double *xx,double *par){
    double E = xx[0];
    double c = par[0];
	double mu = par[1];
	double sigma = par[2];
	double alpha = par[3];
	double n = par[4];

	double A = TMath::Power((n/TMath::Abs(alpha)),n)*exp(-0.5*alpha*alpha);
	
	double cut = (E-mu)/sigma;
	double result_cb;
	if(cut > -alpha){
		result_cb = exp(-0.5*(E-mu)*(E-mu)/sigma/sigma);
	}
	else{
		result_cb = A*TMath::Power((n/TMath::Abs(alpha)-TMath::Abs(alpha)-(E-mu)/sigma),-n);
	}

	return result_cb*c;
	
}

int main(int argc,char **argv)
{

	string source = argv[1];
	string R = argv[2];	
	string Z = argv[3];	
	string dirname = "/junofs/production/public/users/zhangfy/non-uniform/offline_J17v1r1-Pre1/Examples/Tutorial/share/cls2/"+source+"/"+R+"_"+Z+"/";

	TCanvas *c1 = new TCanvas();
    gStyle->SetOptFit(1);
    gStyle->SetStatX(0.4);
    gStyle->SetStatY(0.9);
    gStyle->SetStatW(0.15);
    gStyle->SetStatH(0.15);
    TFile *file;
	double *pars;
	pars = new double[7];

	TChain *t = new TChain("evt");

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
	
	total_entries = max_entries;
	
	t->Draw("totalPE>>h(200,0,0)","totalPE>0","",total_entries,total_entries-max_entries);
	TH1F *h = (TH1F*)gDirectory->Get("h");

	int maxbin = h->GetMaximumBin();
	double maxbincenter = h->GetBinCenter(maxbin);

	t->Draw(Form("totalPE>>h(100,%i,%i)",int(maxbincenter-100),int(maxbincenter+100)),"","",total_entries,total_entries-max_entries);	
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
	pars[3] = 2;
	pars[4] = 1;
    TF1 *f = new TF1("f",myfun,0,20000,5);
	f->SetParameters(pars);

	f->SetParNames("C","#mu","#sigma","#alpha","n");
	t->Draw(Form("totalPE>>h(%i,%i,%i)",int(5*sigma),int(mean-5*sigma),int(mean-5*sigma)+int(5*sigma)*2),"","",total_entries,total_entries-max_entries);	
	h = (TH1F*)gDirectory->Get("h");
	
	h->Fit(f,"M","");
	pars = f->GetParameters();

	if(source=="Ge68"){
		h->Fit(f,"M","",pars[1]-pars[2]*4,pars[1]+pars[2]);
		h->Fit(f,"M","",pars[1]-pars[2]*4,pars[1]+pars[2]);
		h->Fit(f,"M","",pars[1]-pars[2]*4,pars[1]+pars[2]);
	}
	else{
		h->Fit(f,"M","",int(mean-4*sigma),int(mean-5*sigma)+int(5*sigma)*2);
		h->Fit(f,"M","",int(mean-4*sigma),int(mean-5*sigma)+int(5*sigma)*2);
		h->Fit(f,"M","",int(mean-4*sigma),int(mean-5*sigma)+int(5*sigma)*2);
		h->Fit(f,"M","",int(mean-4*sigma),int(mean-5*sigma)+int(5*sigma)*2);
		h->Fit(f,"M","",int(mean-4*sigma),int(mean-5*sigma)+int(5*sigma)*2);
		h->Fit(f,"M","",int(mean-4*sigma),int(mean-5*sigma)+int(5*sigma)*2);
		h->Fit(f,"M","",int(mean-4*sigma),int(mean-5*sigma)+int(5*sigma)*2);
	}
	pars = f->GetParameters();
	double *epars = f->GetParErrors();
	int ndf = f->GetNDF();
	double chi2 = f->GetChisquare();
/*
    TF1 *f2 = new TF1("f",myfun,0,20000,5);
    TF1 *f3 = new TF1("f",myfun,0,20000,5);
	f2->SetParameters(pars);
	f3->SetParameters(pars);

	f2->SetParameter(0,pars[0]);	
	f2->SetParameter(1,0);	
	f2->SetParameter(2,pars[2]);	
	f2->SetParameter(3,pars[3]);	
	f2->SetParameter(4,pars[4]);	
	f2->SetParameter(5,pars[5]);	
	f2->SetParameter(6,pars[6]);	

	f3->SetParameter(0,pars[0]);	
	f3->SetParameter(1,pars[1]);	
	f3->SetParameter(2,pars[2]);	
	f3->SetParameter(3,pars[3]);	
	f3->SetParameter(4,0);	
	f3->SetParameter(5,pars[5]);	
	f3->SetParameter(6,0);	

	f2->SetLineColor(kBlue);
	f3->SetLineColor(kBlack);
	f2->Draw("same");
	f3->Draw("same");
*/
	filename = dirname+"result_cb.C";
	c1->SaveAs(&filename[0]);
	filename = dirname+"result_cb.png";
	c1->SaveAs(&filename[0]);

	filename = dirname+"result_cb";
	ofstream fout(&filename[0]);

	fout<< chi2 <<"\t"<< ndf <<"\t";

	for(int i=0;i<7;i++){
		if(i==1)	fout << 0 << "\t" << 0 <<"\t";
		fout << pars[i] <<"\t" << epars[i] <<"\t";
	}

	fout << "\n";

	return 0;
}
