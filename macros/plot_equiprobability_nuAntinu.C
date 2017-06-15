//////////////////////////////////////////////////////////
///// Testing NovA numu2nue probability
//// check equiprobability using neutrino and antineutrino
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
    

    
    TMarker *InputPoint= new TMarker(myTheta13,mycp,21);
    ////////////////////////////////
    
	double xcoord = 0.21;
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
	
	TString CPStr=Form("#delta_{CP}=%2.2f",mycp);
    TLatex* tlxcp=new TLatex(xcoord, 0.60,CPStr);
    tlxcp->SetNDC(kTRUE); // <- use NDC coordinate
    tlxcp->SetTextSize(0.04);
    tlxcp->SetTextAlign(12);
	
    TString LStr=Form("L=%2.0f km",myL);
    TLatex* tlx5=new TLatex(xcoord, 0.60,LStr);
    tlx5->SetNDC(kTRUE); // <- use NDC coordinate
    tlx5->SetTextSize(0.04);
    tlx5->SetTextAlign(12);
    
    TString EnStr="E_{#nu}=2 GeV";
    TLatex* tlx6=new TLatex(xcoord, 0.55,EnStr);
    tlx6->SetNDC(kTRUE); // <- use NDC coordinate
    tlx6->SetTextSize(0.04);
    tlx6->SetTextAlign(12);
    
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
    
    
    TFile *poutFile = new TFile("../plots/equiprob_"+fprob+"_1cont_nuAntinu.root");
    TH2F* hequiprob = (TH2F*)poutFile->Get("hequiprob");
    TH2F* hequiprobAnti= (TH2F*)poutFile->Get("hequiprobAnti");
    //input neutrino L=810
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
    //input antineutrino L=810
    double parInputL2[9]={myTheta12,                       //Theta12
        myTheta23,                  //Theta23
        myTheta13,                       //Theta13
        myDeltam32*TMath::Power(10,-3.), //DeltaM32^2
        myDeltam21*TMath::Power(10,-5.),  //DeltaM21^2
        mycp,                       //DeltaCP (3Pi/2)
        2.71,                        //Rock Density
        myL,                        //L
        -1};
    SetOscParam(parInputL2);
    // x[0]=2.0+epsilon;
    double problevelL2=model3flav_mat_kopp_2nd_MuToElec(x);
    
    
  
    //fake histogram
    TH1* h1equiprob = new TH1F("h1","h1",npointsx,xmin,xmax);;
    TH1* h2equiprob = new TH1F("h2","h2",npointsx,xmin,xmax);;
    //UT orange
    Int_t ci1 =  TColor::GetColor("#B45F04");
    //navy blue
    Int_t ci2 = TColor::GetColor("#006699");
    
    h1equiprob->SetLineColor(ci1);
    h2equiprob->SetLineColor(ci2);
    
    TLegend *leg = new TLegend(xcoord,0.2,xcoord+0.24,0.5);
    leg->SetFillStyle(0);
	leg->SetBorderSize(0);
    leg->SetTextSize(0.04);
    
    leg->AddEntry(InputPoint,"Input parameter","p");
	leg->AddEntry(h1equiprob,"Neutrino 810 km","l");
    leg->AddEntry(h2equiprob,"AntiNeutrino 810 km","l");
    
    //leg->AddEntry(MuElec_right,Th13StrR,"l");
    
	
			TCanvas *canvas1=0;
        gStyle->SetLineWidth(2);

			canvas1 = new TCanvas("Can1","Can1",800,600);
			gStyle->SetOptStat(0);
			canvas1->cd();
			Int_t colors[]={ci1,ci2};
			gStyle->SetPalette((sizeof(colors)/sizeof(Int_t)), colors);
			Double_t levels[]={problevel,1.0};
			hequiprob->SetContour((sizeof(levels)/sizeof(Double_t)), levels);
			hequiprob->SetTitle("");
    hequiprob->GetXaxis()->SetLabelSize(hequiprob->GetXaxis()->GetTitleSize()*1.2);
    hequiprob->GetYaxis()->SetLabelSize(hequiprob->GetYaxis()->GetTitleSize()*1.2);
    hequiprob->GetXaxis()->SetTitleSize(hequiprob->GetXaxis()->GetLabelSize()*1.2);
    hequiprob->GetYaxis()->SetTitleSize(hequiprob->GetYaxis()->GetLabelSize()*1.2);
    hequiprob->GetZaxis()->SetTitleSize(hequiprob->GetZaxis()->GetLabelSize()*1.2);
    hequiprob->GetYaxis()->SetTitleOffset(0.9);
    hequiprob->GetXaxis()->SetTitle("Mixing angle #theta_{13}");
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
			canvas1->Print("../plots/equiprob_"+fprob+"_1cont_nuAntinu.eps");
			canvas1->Print("../plots/equiprob_"+fprob+"_1cont_nuAntinu.pdf");


    
}
