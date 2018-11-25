#include "rootheader.h"
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

    return c*(gaus+exp_c1*exp1+constant)/(alpha+beta+gamma);
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
	double *pars;
	pars = new double[7];

	TChain *t = new TChain("t");

//	t->Add("/junofs/production/public/users/zhangfy/non-uniform/offline_J17v1r1-Pre1/Examples/Tutorial/share/cls2/Co60/10000_10000/evt_*root");
	string filename = dirname+"new_*root";
	t->Add(&filename[0]);
	int total_entries = t->GetEntries();
		
	int max_entries;
	if(source=="Ge68")
		max_entries = 4000000;
	else
		max_entries = 400000;
	
	t->Draw("delayPE>>h(200,5000,8000)","delayPE>0","",total_entries,total_entries-max_entries);
	TH1F *h = (TH1F*)gDirectory->Get("h");

    int maxbin = h->GetMaximumBin();
    double maxbincenter = h->GetBinCenter(maxbin);

    t->Draw(Form("delayPE>>h(100,%i,%i)",int(maxbincenter-100),int(maxbincenter+100)),"","",total_entries,total_entries-max_entries);
    h = (TH1F*)gDirectory->Get("h");
	h->Fit("gaus");
	TF1 *f = h->GetFunction("gaus"); 

    double C = f->GetParameter(0);
    double mean = f->GetParameter(1);
    double emean = f->GetParError(2);
    double sigma = f->GetParameter(2);
	t->Draw(Form("delayPE>>h(%i,%i,%i)",int(2.5*sigma),int(mean-5*sigma),int(mean-5*sigma)+int(2.5*sigma)*4),"","",total_entries,total_entries-max_entries);
    h = (TH1F*)gDirectory->Get("h");
    if(R=="16000" || R=="17000" || R=="17400"){
        f = new TF1("f",myfun,0,20000,7);
        pars[0] = C*sqrt(2*3.14159)*sigma;
        pars[1] = 0.99;
        pars[2] = mean;
        pars[3] = sigma;
        pars[4] = 0.001;
        pars[5] = 0.01;
        pars[6] = 0.01;
        f->SetParameters(pars);
        f->SetParNames("C","#alpha","#mu","#sigma","#beta","#lambda","#gamma");
        h->Fit(f);
    }
    else{
        h->Fit("gaus");
        f = h->GetFunction("gaus");
    }

	pars = f->GetParameters();
	double *epars = new double[7];
	epars = f->GetParErrors();
	
	filename = dirname+"result_nC.C";
	c1->SaveAs(&filename[0]);
	filename = dirname+"result_nC.png";
	c1->SaveAs(&filename[0]);

	filename = dirname+"result_nC";

	ofstream fout(&filename[0]);

    int ndf = f->GetNDF();
    double chi2 = f->GetChisquare();
	
	fout << chi2 << "\t"<< ndf <<"\t";

	for(int i=0;i<7;i++){
		if(i==1) fout << "0\t0\t";
		fout << pars[i] <<"\t" << epars[i] <<"\t";
	}

	fout << "\n";

	return 0;
}
