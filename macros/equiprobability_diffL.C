//////////////////////////////////////////////////////////
///// Testing NovA numu2nue probability
//// check equiprobability different L
/////////////////////////////////////////////////////////
{

#include<TH1>
#include<TH2>
#include<TLatex>
#include<TString>
#include<TMath>
#include<TFile>    
    
    gROOT->ProcessLine(".L ../source/NuOscProb.C+");
    
    double x[1]=2.0; //Nova Energy
    
    TString fcont="nova_";
    TString fprob="model3flav_mat_kopp_2nd_MuToElec";//model3flav_mat_MuToElec
    TString fsubname="";
	//Oscillation parameter setup
    /* namespace OscPar
     {
     typedef enum EOscPar{
     kTh12 = 0,
     kTh23 = 1,
     kTh13 = 2,
     kDeltaM23 = 3,
     kDeltaM12 = 4,
     kDelta = 5,
     kDensity = 6,
     kL = 7,
     kNuAntiNu = 8,
     kUnknown = 9
     } OscPar_t;
     */
    
    double myDeltam32=2.42;//update
	double myDeltam21=7.54;
    double myTheta23=0.78;
	double myTheta12=0.61;
	double myTheta13=0.15;
    double mycp=4.0;
    double myL=810;//new baseline 1250
	double shifth13 =0.02;
    double problevel;
    
    double epsilon=0.00001;
    
    Int_t npointsx =50;//for theta13 steps
	double xmin =0;
	double xmax =0.20;
	Int_t npointsy=50;//for cp steps
	double ymin =0;
	double ymax = 2.0*TMath::Pi();
    
    //input L=810
    double parInput[9]={myTheta12,                       //Theta12
        myTheta23,                  //Theta23
        myTheta13,                       //Theta13
        myDeltam32*TMath::Power(10,-3.), //DeltaM32^2
        myDeltam21*TMath::Power(10,-5.),  //DeltaM21^2
        mycp,                       //DeltaCP (3Pi/2)
        2.71,                        //Rock Density
        myL,                        //L
        1};
    SetOscParam(parInput);
    x[0]=2.0+epsilon;
    double problevel=model3flav_mat_kopp_2nd_MuToElec(x);
    
    //input L=1250
    double parInputL2[9]={myTheta12,                       //Theta12
        myTheta23,                  //Theta23
        myTheta13,                       //Theta13
        myDeltam32*TMath::Power(10,-3.), //DeltaM32^2
        myDeltam21*TMath::Power(10,-5.),  //DeltaM21^2
        mycp,                       //DeltaCP (3Pi/2)
        2.71,                        //Rock Density
        1250,                        //L
        1};
    SetOscParam(parInputL2);
   // x[0]=2.0+epsilon;
    double problevelL2=model3flav_mat_kopp_2nd_MuToElec(x);
    
    
    TMarker *InputPoint= new TMarker(myTheta13,mycp,21);
    ////////////////////////////////
    
	double xcoord = 0.18;
    TString delMStr32=Form("%2.2f #times 10^{-3}",myDeltam32);
    TString legStr32="|#Delta m^{2}_{32}|=";
    TLatex* tlx=new TLatex(xcoord, 0.85,legStr32+delMStr32);
    tlx->SetNDC(kTRUE); // <- use NDC coordinate
    tlx->SetTextSize(0.04);
	
	TString delMStr21=Form("%2.2f #times 10^{-5}",myDeltam21);
    TString legStr21="|#Delta m^{2}_{21}|=";
    TLatex* tlx1=new TLatex(xcoord, 0.80,legStr21+delMStr21);
    tlx1->SetNDC(kTRUE); // <- use NDC coordinate
    tlx1->SetTextSize(0.04);
	
	
    TString Th23Str=Form("sin^{2}2#theta_{23}=%2.2f",TMath::Sin(2*myTheta23)*TMath::Sin(2*myTheta23));
    TLatex* tlx2=new TLatex(xcoord, 0.75,Th23Str);
    tlx2->SetNDC(kTRUE); // <- use NDC coordinate
    tlx2->SetTextSize(0.04);
	
	TString Th12Str=Form("sin^{2}2#theta_{12}=%2.2f",TMath::Sin(2*myTheta12)*TMath::Sin(2*myTheta12));
    TLatex* tlx3=new TLatex(xcoord, 0.70,Th12Str);
    tlx3->SetNDC(kTRUE); // <- use NDC coordinate
    tlx3->SetTextSize(0.04);
	
	TString Th13Str=Form("sin^{2}2#theta_{13}=%2.2f",TMath::Sin(2*myTheta13)*TMath::Sin(2*myTheta13));
    TLatex* tlx4=new TLatex(xcoord, 0.65,Th13Str);
    tlx4->SetNDC(kTRUE); // <- use NDC coordinate
    tlx4->SetTextSize(0.04);
	
	TString CPStr=Form("#delta_{CP}=%2.2f",mycp);
    TLatex* tlxcp=new TLatex(xcoord, 0.60,CPStr);
    tlxcp->SetNDC(kTRUE); // <- use NDC coordinate
    tlxcp->SetTextSize(0.04);
    
	
    TString LStr=Form("L=%2.0f km",myL);
    TLatex* tlx5=new TLatex(xcoord, 0.60,LStr);
    tlx5->SetNDC(kTRUE); // <- use NDC coordinate
    tlx5->SetTextSize(0.04);
    
    TString EnStr="E_{#nu}=2 GeV";
    TLatex* tlx6=new TLatex(xcoord, 0.55,EnStr);
    tlx6->SetNDC(kTRUE); // <- use NDC coordinate
    tlx6->SetTextSize(0.04);
    
    TString ProbStr=Form("P(#nu_{#mu}#rightarrow #nu_{e})=%2.2f",problevel);
    TLatex* tlxprob=new TLatex(xcoord, 0.50,ProbStr);
    tlxprob->SetNDC(kTRUE); // <- use NDC coordinate
    tlxprob->SetTextSize(0.04);
    
    TString NorHierStr="#color[4]{#Delta m^{2}_{32}>0}";
    TLatex* tlx7=new TLatex(0.63, 0.38,NorHierStr);
    tlx7->SetNDC(kTRUE); // <- use NDC coordinate
    tlx7->SetTextSize(0.04);
    
    TString InvHierStr="#color[2]{#Delta m^{2}_{32}<0}";
    TLatex* tlx8=new TLatex(0.2, 0.55,InvHierStr);
    tlx8->SetNDC(kTRUE); // <- use NDC coordinate
    tlx8->SetTextSize(0.04);
    
    
    TString LStrname=Form("%2.0f",myL);
    

    
    
    //neutrino normal hierarchy
    TH2F *hequiprob=0;
    TH1F *h1equiprob=new TH1F("h1","h1",npointsx,xmin,xmax);
    

    //-----------
    hequiprob = new TH2F("hcont_nu_norm","hcont_nu_norm",npointsx, xmin, xmax, npointsy, ymin, ymax);
    //hequiprob = new TH2F("hcont_nu_norm","hcont_nu_norm",npointsx, binx, npointsy, biny);
    hequiprob->GetXaxis()->SetTitle("#theta_{13}");
    hequiprob->GetYaxis()->SetTitle("CP phase");
    hequiprob->GetXaxis()->CenterTitle();
    hequiprob->GetYaxis()->CenterTitle();
	hequiprob->GetXaxis()->SetTitleSize(0.05);
    hequiprob->GetYaxis()->SetTitleSize(0.05);
	hequiprob->GetYaxis()->SetTitleOffset(1.0);
    hequiprob->GetXaxis()->SetTitleOffset(0.9);
    
    TH2F *hequiprobAnti=0;
    TH1F *h2equiprob=new TH1F("h2","h2",npointsx,xmin,xmax);
    hequiprobAnti = new TH2F("hcont_nu_normAnti","hcont_nu_normAnti",npointsx, xmin, xmax, npointsy, ymin, ymax);
    
    ///////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////
    //neutrino w L = 810km
	for (Int_t ix=0; ix < npointsx+1; ix++){
		double ptheta13 = ix*(xmax-xmin)/npointsx;
        for (Int_t iy=0; iy < npointsy+1; iy++){
			double pcp = iy*(ymax-ymin)/npointsy;
			
			double par[9]={myTheta12,                       //Theta12
				myTheta23,                  //Theta23
				ptheta13,                       //Theta13
				myDeltam32*TMath::Power(10,-3.), //DeltaM32^2
				myDeltam21*TMath::Power(10,-5.),  //DeltaM21^2
				pcp,                       //DeltaCP (3Pi/2)
				2.71,                        //Rock Density
				myL,                        //L
				1};
			SetOscParam(par);
			x[0]=2.0+epsilon;
			double pprob=model3flav_mat_kopp_2nd_MuToElec(x);
			cout<<"pmuex["<<ix<<","<<iy<<"] = "<<pprob<<endl;
			
			hequiprob->SetBinContent(ix,iy,pprob);
			
		}
	}
    
    //neutrino w L =1250km
	for (Int_t ix=0; ix < npointsx+1; ix++){
		double ptheta13 = ix*(xmax-xmin)/npointsx;
        for (Int_t iy=0; iy < npointsy+1; iy++){
			double pcp = iy*(ymax-ymin)/npointsy;
			
			double par[9]={myTheta12,                       //Theta12
				myTheta23,                  //Theta23
				ptheta13,                       //Theta13
				myDeltam32*TMath::Power(10,-3.), //DeltaM32^2
				myDeltam21*TMath::Power(10,-5.),  //DeltaM21^2
				pcp,                       //DeltaCP (3Pi/2)
				2.71,                        //Rock Density
				1250,                        //L
				1};
			SetOscParam(par);
			x[0]=2.0+epsilon;
			double pprob=model3flav_mat_kopp_2nd_MuToElec(x);
			cout<<"pmuex["<<ix<<","<<iy<<"] = "<<pprob<<endl;
			
			hequiprobAnti->SetBinContent(ix,iy,pprob);
			
		}
	}
    //fake histogram
    h1equiprob->SetLineColor(2);
    h2equiprob->SetLineColor(4);
    
    TLegend *leg = new TLegend(xcoord,0.2,xcoord+0.24,0.5);
    leg->SetFillStyle(0);
	leg->SetBorderSize(0);
    leg->SetTextSize(0.04);
    
    leg->AddEntry(InputPoint,"Input parameter","p");
	leg->AddEntry(h1equiprob,"Neutrino 810 km","l");
    leg->AddEntry(h2equiprob,"Neutrino 1250 km","l");
    
    //leg->AddEntry(MuElec_right,Th13StrR,"l");
    
	
			TCanvas *canvas1=0;
			canvas1 = new TCanvas("Can1","Can1",800,600);
			gStyle->SetOptStat(0);
			canvas1->cd();
			Int_t colors[]={2,4};
			gStyle->SetPalette((sizeof(colors)/sizeof(Int_t)), colors);
			Double_t levels[]={problevel,1.0};
			hequiprob->SetContour((sizeof(levels)/sizeof(Double_t)), levels);
			hequiprob->SetTitle("");
			hequiprob->Draw("cont1");
    
            Double_t levelsAnti[]={0.0,problevelL2};
            hequiprobAnti->SetContour((sizeof(levelsAnti)/sizeof(Double_t)), levelsAnti);
            hequiprobAnti->SetTitle("");
            hequiprobAnti->Draw("cont1 same");
    
            InputPoint->Draw("same");
			tlx->Draw("same");
			tlx1->Draw("same");
			tlx2->Draw("same");
			tlx3->Draw("same");
			tlx4->Draw("same");
            tlxcp->Draw("same");
			//tlx5->Draw("same");
			tlx6->Draw("same");
            //tlxprob->Draw("same");
            leg->Draw("same");
			canvas1->Update();
			canvas1->Print("../plots/equiprob_"+fprob+"_1cont_diffL.eps");
			canvas1->Print("../plots/equiprob_"+fprob+"_1cont_diffL.pdf");
    TFile *poutFile = new TFile("../plots/equiprob_"+fprob+"_1cont_diffL.root","RECREATE");
    hequiprob->Write("hequiprob");
    hequiprobAnti->Write("hequiprobAnti");
    poutFile->Close();
		

    exit(0);
    
}
