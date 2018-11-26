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

    //return c*(gaus+exp_c1*exp1+constant)/(alpha+beta+gamma);
    return c*(gaus+exp_c1*exp1+constant);
}

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
    TF1 *f = new TF1("f",myfun,0,20000,7);

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

	float edep;
	if(source == "Ge68") edep = 1.0219978;
	if(source == "Cs137") edep = 0.661657;
	if(source == "Mn54") edep = 0.834848;
	if(source == "Co60") edep = 2.5057385;
	if(source == "K40") edep = 1.4608;
	
	//t->Draw("totalPE>>h(200,0,0)",Form("totalPE>0 && TMath::Abs(edep-%f)<0.001",edep),"",max_entries);
	t->Draw("totalPE>>h(200,0,0)",Form("totalPE>0 && TMath::Abs(edep-%f)<0.001",edep),"",max_entries,max_entries*index);

	TH1F *h = (TH1F*)gDirectory->Get("h");

	int maxbin = h->GetMaximumBin();
	double maxbincenter = h->GetBinCenter(maxbin);

	t->Draw(Form("totalPE>>h(100,%i,%i)",int(maxbincenter-100),int(maxbincenter+100)),Form("totalPE>0 && TMath::Abs(edep-%f)<0.0001",edep),"",max_entries,max_entries*index);	
	//t->Draw(Form("totalPE>>h(100,%i,%i)",int(maxbincenter-100),int(maxbincenter+100)),Form("totalPE>0 && TMath::Abs(edep-%f)<0.001",edep),"",max_entries);	
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
    f = h->GetFunction("gaus");

    double C = f->GetParameter(0);
    double mean = f->GetParameter(1);
    double emean = f->GetParError(2);
    double sigma = f->GetParameter(2);

	//f->SetParameters(C*sqrt(2*3.14159)*sigma,0.99,mean,sigma,0.001,0.01,0.02);
	//t->Draw(Form("totalPE>>h(%i,%i,%i)",int(5*sigma),int(mean-5*sigma),int(mean-5*sigma)+int(5*sigma)*2),Form("totalPE>0 && TMath::Abs(edep-%f)<0.001",edep),"",max_entries);	
	t->Draw(Form("totalPE>>h(%i,%i,%i)",int(5*sigma),int(mean-5*sigma),int(mean-5*sigma)+int(5*sigma)*2),Form("totalPE>0 && TMath::Abs(edep-%f)<0.001",edep),"",max_entries,max_entries*index);	


	h = (TH1F*)gDirectory->Get("h");
	
	h->Fit("gaus","M","");
    f = h->GetFunction("gaus");
	
	f->GetParameters(pars);
	double *tmps = f->GetParErrors();
    memcpy(epars,tmps,7*sizeof(double));

    ndf = f->GetNDF();
    chi2 = f->GetChisquare();

	filename = dirname+"result_nocom.C";
	c1->SaveAs(&filename[0]);
	filename = dirname+"result_nocom.png";
	c1->SaveAs(&filename[0]);


	return true;
}
int main(int argc,char **argv){

    double *pars = new double[7];
    double *epars = new double[7];
    double chi2;
    double ndf;
    string source = argv[1];
    string R = argv[2];
    string Z = argv[3];

    string dirname = "/junofs/production/public/users/zhangfy/non-uniform/offline_J17v1r1-Pre1/Examples/Tutorial/share/cls2/"+source+"/"+R+"_"+Z+"/";
    string filename = dirname+"result_nocom";
    ofstream fout(&filename[0]);
    int index = 0;

    while(get_data(argc,argv,pars,epars,chi2,ndf,dirname,index) && index <1)
	{

        fout<< chi2 <<"\t"<< ndf <<"\t";

        for(int i=0;i<7;i++){
            fout << pars[i] <<"\t" << epars[i] <<"\t";
        }

        fout << "\n";
        index++;
    }

    return 0;
}
