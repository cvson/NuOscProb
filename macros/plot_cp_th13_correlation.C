//////////////////////////////////////////////////////////
///// Testing NovA numu2nue probability
//// check equiprobability
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
    double myL=810;//new baseline 1250
	double shifth13 =0.02;
    
	double xcoord = 0.20;
    TString delMStr32=Form("%2.2f #times 10^{-3}",myDeltam32);
    TString legStr32="|#Delta m^{2}_{32}|=";
    TLatex* tlx=new TLatex(xcoord, 0.85,legStr32+delMStr32);
    tlx->SetNDC(kTRUE); // <- use NDC coordinate
    tlx->SetTextSize(0.04);
       tlx->SetTextAlign(12);
	
	TString delMStr21=Form("%2.2f #times 10^{-5}",myDeltam21);
    TString legStr21="|#Delta m^{2}_{21}|=";
    TLatex* tlx1=new TLatex(xcoord, 0.80,legStr21+delMStr21);
    tlx1->SetNDC(kTRUE); // <- use NDC coordinate
    tlx1->SetTextSize(0.04);
	tlx1->SetTextAlign(12);
	
    TString Th23Str=Form("sin^{2}2#theta_{23}=%2.2f",TMath::Sin(2*myTheta23)*TMath::Sin(2*myTheta23));
    TLatex* tlx2=new TLatex(xcoord, 0.75,Th23Str);
    tlx2->SetNDC(kTRUE); // <- use NDC coordinate
    tlx2->SetTextSize(0.04);
	tlx2->SetTextAlign(12);
	TString Th12Str=Form("sin^{2}2#theta_{12}=%2.2f",TMath::Sin(2*myTheta12)*TMath::Sin(2*myTheta12));
    TLatex* tlx3=new TLatex(xcoord, 0.70,Th12Str);
    tlx3->SetNDC(kTRUE); // <- use NDC coordinate
    tlx3->SetTextSize(0.04);
	tlx3->SetTextAlign(12);
	TString Th13Str=Form("sin^{2}2#theta_{13}=%2.2f",TMath::Sin(2*myTheta13)*TMath::Sin(2*myTheta13));
    TLatex* tlx4=new TLatex(xcoord, 0.65,Th13Str);
    tlx4->SetNDC(kTRUE); // <- use NDC coordinate
    tlx4->SetTextSize(0.04);
	tlx4->SetTextAlign(12);
	
	
    TString LStr=Form("L=%2.0f km",myL);
    TLatex* tlx5=new TLatex(xcoord, 0.65,LStr);
    tlx5->SetNDC(kTRUE); // <- use NDC coordinate
    tlx5->SetTextSize(0.04);
    tlx5->SetTextAlign(12);
    TString EnStr="E_{#nu}=2 GeV";
    TLatex* tlx6=new TLatex(xcoord, 0.60,EnStr);
    tlx6->SetNDC(kTRUE); // <- use NDC coordinate
    tlx6->SetTextSize(0.04);
    tlx6->SetTextAlign(12);
    TString NorHierStr="#color[4]{#Delta m^{2}_{32}>0}";
    TLatex* tlx7=new TLatex(0.63, 0.38,NorHierStr);
    tlx7->SetNDC(kTRUE); // <- use NDC coordinate
    tlx7->SetTextSize(0.04);
    tlx6->SetTextAlign(12);
    TString InvHierStr="#color[2]{#Delta m^{2}_{32}<0}";
    TLatex* tlx8=new TLatex(0.2, 0.55,InvHierStr);
    tlx8->SetNDC(kTRUE); // <- use NDC coordinate
    tlx8->SetTextSize(0.04);
    tlx8->SetTextAlign(12);
    

    TFile *poutFile = new TFile("../plots/equiprob_"+fprob+".root");
    TH2F *hequiprob=(TH2F*)poutFile->Get("hequiprob");
	
    TCanvas *canvas1=0;
    canvas1 = new TCanvas("Can1","Can1",800,600);
    gStyle->SetOptStat(0);
    canvas1->cd();
    gStyle->SetLineWidth(2);
    gPad->SetRightMargin(gPad->GetRightMargin()*2.);
    Int_t colors[]={0,41,42,45,95};
    gStyle->SetPalette((sizeof(colors)/sizeof(Int_t)), colors);
    //Double_t levels[]={0.001,0.01,0.05,0.1};
    //hequiprob->SetContour((sizeof(levels)/sizeof(Double_t)), levels);
    hequiprob->SetTitle("");
    hequiprob->GetYaxis()->CenterTitle();
    hequiprob->GetXaxis()->CenterTitle();
    hequiprob->GetXaxis()->SetLabelSize(hequiprob->GetXaxis()->GetTitleSize()*1.2);
    hequiprob->GetYaxis()->SetLabelSize(hequiprob->GetYaxis()->GetTitleSize()*1.2);
    hequiprob->GetXaxis()->SetTitleSize(hequiprob->GetXaxis()->GetLabelSize()*1.2);
    hequiprob->GetYaxis()->SetTitleSize(hequiprob->GetYaxis()->GetLabelSize()*1.2);
    hequiprob->GetZaxis()->SetTitleSize(hequiprob->GetZaxis()->GetLabelSize()*1.2);
    hequiprob->GetYaxis()->SetTitleOffset(0.9);
    hequiprob->Draw("colz");
    tlx->Draw("same");
    tlx1->Draw("same");
    tlx2->Draw("same");
    tlx3->Draw("same");
    //tlx4->Draw("same");
    tlx5->Draw("same");
    tlx6->Draw("same");
			
    canvas1->Update();
    canvas1->Print("../plots/equiprob_"+fprob+".eps");
    canvas1->Print("../plots/equiprob_"+fprob+".pdf");
		
 
    
}
