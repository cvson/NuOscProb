//////////////////////////////////////////////////////////
///// Testing NovA numu2nue probability
//// consider different octant
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
    double myTheta23=0.67264;//lower octant sin22theta23=0.95
	double myTheta23highOc=0.89815;//high octant sin22theta23=0.95
	double myTheta12=0.61;
	double myTheta13=0.15;
    double myL=1250;//new baseline 1250
    double xshift = 0.1;
    double yshift = 0.03;
    TString delMStr32=Form("%2.2f #times 10^{-3}",myDeltam32);
    TString legStr32="|#Delta m^{2}_{32}|=";
    TLatex* tlx=new TLatex(0.6+xshift, 0.85+yshift,legStr32+delMStr32);
    tlx->SetNDC(kTRUE); // <- use NDC coordinate
    tlx->SetTextSize(0.04);
    tlx->SetTextAlign(12);
	
	TString delMStr21=Form("%2.2f #times 10^{-5}",myDeltam21);
    TString legStr21="|#Delta m^{2}_{21}|=";
    TLatex* tlx1=new TLatex(0.6+xshift, 0.80+yshift,legStr21+delMStr21);
    tlx1->SetNDC(kTRUE); // <- use NDC coordinate
    tlx1->SetTextSize(0.04);
    tlx1->SetTextAlign(12);
	
	
    TString Th23Str=Form("sin^{2}2#theta_{23}=%2.2f",TMath::Sin(2*myTheta23)*TMath::Sin(2*myTheta23));
    TLatex* tlx2=new TLatex(0.6+xshift, 0.75+yshift,Th23Str);
    tlx2->SetNDC(kTRUE); // <- use NDC coordinate
    tlx2->SetTextSize(0.04);
    tlx2->SetTextAlign(12);
	
	TString Th12Str=Form("sin^{2}2#theta_{12}=%2.2f",TMath::Sin(2*myTheta12)*TMath::Sin(2*myTheta12));
    TLatex* tlx3=new TLatex(0.6+xshift, 0.70+yshift,Th12Str);
    tlx3->SetNDC(kTRUE); // <- use NDC coordinate
    tlx3->SetTextSize(0.04);
    tlx3->SetTextAlign(12);
	
	TString Th13Str=Form("sin^{2}2#theta_{13}=%2.2f",TMath::Sin(2*myTheta13)*TMath::Sin(2*myTheta13));
    TLatex* tlx4=new TLatex(0.6+xshift, 0.65+yshift,Th13Str);
    tlx4->SetNDC(kTRUE); // <- use NDC coordinate
    tlx4->SetTextSize(0.04);
    tlx4->SetTextAlign(12);
	
    TString LStr=Form("L=%2.0f km",myL);
    TLatex* tlx5=new TLatex(0.6+xshift, 0.60+yshift,LStr);
    tlx5->SetNDC(kTRUE); // <- use NDC coordinate
    tlx5->SetTextSize(0.04);
    tlx5->SetTextAlign(12);
    
    TString EnStr="E_{#nu}=2 GeV";
    TLatex* tlx6=new TLatex(0.6+xshift, 0.55+yshift,EnStr);
    tlx6->SetNDC(kTRUE); // <- use NDC coordinate
    tlx6->SetTextSize(0.04);
    tlx6->SetTextAlign(12);
    
    TString NorHierStr="#color[864]{#Delta m^{2}_{32}>0}";
    TLatex* tlx7=new TLatex(0.63, 0.40,NorHierStr);
    tlx7->SetNDC(kTRUE); // <- use NDC coordinate
    tlx7->SetTextSize(0.04);
    
    TString InvHierStr="#color[924]{#Delta m^{2}_{32}<0}";
    TLatex* tlx8=new TLatex(0.23, 0.53,InvHierStr);
    tlx8->SetNDC(kTRUE); // <- use NDC coordinate
    tlx8->SetTextSize(0.04);
    
    
    TString LStrname=Form("%2.0f",myL);
    
    double epsilon=0.00001;
    
    Int_t dcpstep =200;//for dcp steps
    
    // normal hierarchy
    TGraph *MuElec=0;
	Double_t pmuex[dcpstep+1], pmuey[dcpstep+1];
    
    // inverted hierarchy
	TGraph *MuElecInv=0;
    Double_t pmuexInv[dcpstep+1], pmueyInv[dcpstep+1];
    
    ///////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////
    ///neutrino
    // normal hierarchy
    for (Int_t idcp =0; idcp<dcpstep+1; idcp++) {
		//get dcp value
		double vdcp = idcp*2.0*TMath::Pi()/dcpstep;
        double par[9]={myTheta12,                       //Theta12
			myTheta23,                  //Theta23
			myTheta13,                       //Theta13
			myDeltam32*TMath::Power(10,-3.), //DeltaM32^2
			myDeltam21*TMath::Power(10,-5.),  //DeltaM21^2
			vdcp,                       //DeltaCP (3Pi/2)
			2.71,                        //Rock Density
			myL,                        //L
			1};
        SetOscParam(par);
        x[0]=2.0+epsilon;
        pmuex[idcp]=model3flav_mat_kopp_2nd_MuToElec(x);
        
        cout<<"pmuex["<<idcp<<"] = "<<pmuex[idcp]<<endl;
		

    }
    ///////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////
    ///Antineutrino
    // normal hierarchy
    for (Int_t idcp =0; idcp<dcpstep+1; idcp++) {
		//get dcp value
		double vdcp = idcp*2.0*TMath::Pi()/dcpstep;
        double par[9]={myTheta12,                       //Theta12
			myTheta23,                  //Theta23
			myTheta13,                       //Theta13
			myDeltam32*TMath::Power(10,-3.), //DeltaM32^2
			myDeltam21*TMath::Power(10,-5.),  //DeltaM21^2
			vdcp,                       //DeltaCP (3Pi/2)
			2.71,                        //Rock Density
			myL,                        //L
			-1};
        SetOscParam(par);
        x[0]=2.0+epsilon;
        pmuey[idcp]=model3flav_mat_kopp_2nd_MuToElec(x);
        
        cout<<"pmuex["<<idcp<<"] = "<<pmueyInv[idcp]<<endl;
		
        
    }
    
    ///////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////
    ///neutrino
    // inverted hierarchy
    for (Int_t idcp =0; idcp<dcpstep+1; idcp++) {
		//get dcp value
		double vdcp = idcp*2.0*TMath::Pi()/dcpstep;
        double par[9]={myTheta12,                       //Theta12
			myTheta23,                  //Theta23
			myTheta13,                       //Theta13
			-1.0*myDeltam32*TMath::Power(10,-3.), //DeltaM32^2
			myDeltam21*TMath::Power(10,-5.),  //DeltaM21^2
			vdcp,                       //DeltaCP (3Pi/2)
			2.71,                        //Rock Density
			myL,                        //L
			1};
        SetOscParam(par);
        x[0]=2.0+epsilon;
        pmuexInv[idcp]=model3flav_mat_kopp_2nd_MuToElec(x);
        
        cout<<"pmuex["<<idcp<<"] = "<<pmuexInv[idcp]<<endl;
		
        
    }
    
    ///////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////
    ///Antineutrino
    // inverted hierarchy
    for (Int_t idcp =0; idcp<dcpstep+1; idcp++) {
		//get dcp value
		double vdcp = idcp*2.0*TMath::Pi()/dcpstep;
        double par[9]={myTheta12,                       //Theta12
			myTheta23,                  //Theta23
			myTheta13,                       //Theta13
			-1*myDeltam32*TMath::Power(10,-3.), //DeltaM32^2
			myDeltam21*TMath::Power(10,-5.),  //DeltaM21^2
			vdcp,                       //DeltaCP (3Pi/2)
			2.71,                        //Rock Density
			myL,                        //L
			-1};
        SetOscParam(par);
        x[0]=2.0;
        pmueyInv[idcp]=model3flav_mat_kopp_2nd_MuToElec(x);
        
        cout<<"pmuex["<<idcp<<"] = "<<pmuey[idcp]<<endl;
		
        
    }

    // normal hierarchy
    TGraph *MuElechighOct=0;
    Double_t pmuehighx[dcpstep+1], pmuehighy[dcpstep+1];
		
    // inverted hierarchy
    TGraph *MuElechighOctInv=0;
    Double_t pmuehighxInv[dcpstep+1], pmuehighyInv[dcpstep+1];
		
    ///////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////
    ///neutrino
    // normal hierarchy
    for (Int_t idcp =0; idcp<dcpstep+1; idcp++) {
        //get dcp value
        double vdcp = idcp*2.0*TMath::Pi()/dcpstep;
        double par[9]={myTheta12,                       //Theta12
            myTheta23highOc,                  //Theta23
            myTheta13,                       //Theta13
            myDeltam32*TMath::Power(10,-3.), //DeltaM32^2
            myDeltam21*TMath::Power(10,-5.),  //DeltaM21^2
            vdcp,                       //DeltaCP (3Pi/2)
            2.71,                        //Rock Density
            myL,                        //L
            1};
        SetOscParam(par);
        x[0]=2.0+epsilon;
        pmuehighx[idcp]=model3flav_mat_kopp_2nd_MuToElec(x);
        
        cout<<"pmuex["<<idcp<<"] = "<<pmuehighx[idcp]<<endl;
			
			
    }
    ///////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////
    ///Antineutrino
    // normal hierarchy
    for (Int_t idcp =0; idcp<dcpstep+1; idcp++) {
        //get dcp value
        double vdcp = idcp*2.0*TMath::Pi()/dcpstep;
        double par[9]={myTheta12,                       //Theta12
            myTheta23highOc,                  //Theta23
            myTheta13,                       //Theta13
            myDeltam32*TMath::Power(10,-3.), //DeltaM32^2
            myDeltam21*TMath::Power(10,-5.),  //DeltaM21^2
            vdcp,                       //DeltaCP (3Pi/2)
            2.71,                        //Rock Density
            myL,                        //L
            -1};
        SetOscParam(par);
        x[0]=2.0+epsilon;
        pmuehighy[idcp]=model3flav_mat_kopp_2nd_MuToElec(x);
        cout<<"pmuex["<<idcp<<"] = "<<pmuehighyInv[idcp]<<endl;
			
			
    }
		
    ///////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////
    ///neutrino
    // inverted hierarchy
    for (Int_t idcp =0; idcp<dcpstep+1; idcp++) {
        //get dcp value
        double vdcp = idcp*2.0*TMath::Pi()/dcpstep;
        double par[9]={myTheta12,                       //Theta12
            myTheta23highOc,                  //Theta23
            myTheta13,                       //Theta13
            -1.0*myDeltam32*TMath::Power(10,-3.), //DeltaM32^2
            myDeltam21*TMath::Power(10,-5.),  //DeltaM21^2
            vdcp,                       //DeltaCP (3Pi/2)
            2.71,                        //Rock Density
            myL,                        //L
            1};
        SetOscParam(par);
        x[0]=2.0+epsilon;
        pmuehighxInv[idcp]=model3flav_mat_kopp_2nd_MuToElec(x);
			
        cout<<"pmuex["<<idcp<<"] = "<<pmuehighxInv[idcp]<<endl;
			
			
    }
		
    ///////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////
    ///Antineutrino
    // inverted hierarchy
    for (Int_t idcp =0; idcp<dcpstep+1; idcp++) {
        //get dcp value
        double vdcp = idcp*2.0*TMath::Pi()/dcpstep;
        double par[9]={myTheta12,                       //Theta12
            myTheta23highOc,                  //Theta23
            myTheta13,                       //Theta13
            -1*myDeltam32*TMath::Power(10,-3.), //DeltaM32^2
            myDeltam21*TMath::Power(10,-5.),  //DeltaM21^2
            vdcp,                       //DeltaCP (3Pi/2)
            2.71,                        //Rock Density
            myL,                        //L
            -1};
        SetOscParam(par);
        x[0]=2.0;
        pmuehighyInv[idcp]=model3flav_mat_kopp_2nd_MuToElec(x);
			
        cout<<"pmuex["<<idcp<<"] = "<<pmuehighy[idcp]<<endl;
			
			
    }
		
    //UT orange
    Int_t ci1 =  TColor::GetColor("#B45F04");
    //navy blue
    Int_t ci2 = TColor::GetColor("#006699");
    
    MuElec = new TGraph(dcpstep+1,pmuex,pmuey);
    MuElec->SetLineColor(ci2);
    MuElec->SetLineWidth(2);
    
    MuElecInv = new TGraph(dcpstep+1,pmuexInv,pmueyInv);
    MuElecInv->SetLineColor(ci1);
    MuElecInv->SetLineWidth(2);
	
    MuElechighOct = new TGraph(dcpstep+1,pmuehighx,pmuehighy);
    MuElechighOct->SetLineColor(ci2);
    MuElechighOct->SetLineStyle(7);
    MuElechighOct->SetLineWidth(2);
		
    MuElechighOctInv = new TGraph(dcpstep+1,pmuehighxInv,pmuehighyInv);
    MuElechighOctInv->SetLineColor(ci1);
    MuElechighOctInv->SetLineStyle(7);
    MuElechighOctInv->SetLineWidth(2);
    //UT orange
    Int_t ci4 =  TColor::GetColor("#CC3300");
    //navy blue
    Int_t ci3 = TColor::GetColor("#339900");
		
    //cp=0
    TMarker *MuEleccp0= new TMarker(pmuex[0],pmuey[0],21);
    TMarker *MuElecInvcp0= new TMarker(pmuexInv[0],pmueyInv[0],21);
	TMarker *MuElechighOctcp0= new TMarker(pmuehighx[0],pmuehighy[0],21);
	TMarker *MuElechighOctInvcp0= new TMarker(pmuehighxInv[0],pmuehighyInv[0],21);
		
    MuEleccp0->SetMarkerSize(1.3);
    MuElecInvcp0->SetMarkerSize(1.3);
    MuElechighOctcp0->SetMarkerSize(1.3);
    MuElechighOctInvcp0->SetMarkerSize(1.3);
    MuElechighOctcp0->SetMarkerColor(ci3);
    MuElechighOctInvcp0->SetMarkerColor(ci3);
    //cp=pi/2
    TMarker *MuEleccp12= new TMarker(pmuex[dcpstep/4],pmuey[dcpstep/4],8);
    TMarker *MuElecInvcp12= new TMarker(pmuexInv[dcpstep/4],pmueyInv[dcpstep/4],8);
    TMarker *MuElechighOctcp12= new TMarker(pmuehighx[dcpstep/4],pmuehighy[dcpstep/4],8);
    TMarker *MuElechighOctInvcp12= new TMarker(pmuehighxInv[dcpstep/4],pmuehighyInv[dcpstep/4],8);
		
    MuEleccp12->SetMarkerSize(1.3);
    MuElecInvcp12->SetMarkerSize(1.3);
    MuElechighOctcp12->SetMarkerSize(1.3);
    MuElechighOctInvcp12->SetMarkerSize(1.3);
    MuElechighOctcp12->SetMarkerColor(ci3);
    MuElechighOctInvcp12->SetMarkerColor(ci3);
    
    
    //cp=pi
    TMarker *MuEleccp1= new TMarker(pmuex[dcpstep/2],pmuey[dcpstep/2],25);
    TMarker *MuElecInvcp1= new TMarker(pmuexInv[dcpstep/2],pmueyInv[dcpstep/2],25);
    TMarker *MuElechighOctcp1= new TMarker(pmuehighx[dcpstep/2],pmuehighy[dcpstep/2],25);
    TMarker *MuElechighOctInvcp1= new TMarker(pmuehighxInv[dcpstep/2],pmuehighyInv[dcpstep/2],25);
		
    MuEleccp1->SetMarkerSize(1.3);
    MuElecInvcp1->SetMarkerSize(1.3);
    MuElechighOctcp1->SetMarkerSize(1.3);
    MuElechighOctInvcp1->SetMarkerSize(1.3);
    MuElechighOctcp1->SetMarkerColor(ci4);
    MuElechighOctInvcp1->SetMarkerColor(ci4);
    
    //cp=3pi/2
    TMarker *MuEleccp32= new TMarker(pmuex[dcpstep*3/4],pmuey[dcpstep*3/4],24);
    TMarker *MuElecInvcp32= new TMarker(pmuexInv[dcpstep*3/4],pmueyInv[dcpstep*3/4],24);
    TMarker *MuElechighOctcp32= new TMarker(pmuehighx[dcpstep*3/4],pmuehighy[dcpstep*3/4],24);
    TMarker *MuElechighOctInvcp32= new TMarker(pmuehighxInv[dcpstep*3/4],pmuehighyInv[dcpstep*3/4],24);
		
    MuEleccp32->SetMarkerSize(1.3);
    MuElecInvcp32->SetMarkerSize(1.3);
    MuElechighOctcp32->SetMarkerSize(1.3);
    MuElechighOctInvcp32->SetMarkerSize(1.3);
    MuElechighOctcp32->SetMarkerColor(ci4);
    MuElechighOctInvcp32->SetMarkerColor(ci4);
    
    
    
    //for  normal hierarchy only
    TCanvas *canvas1=0;
    gStyle->SetLineWidth(2);
    canvas1 = new TCanvas("Can1","Can1",800,600);
    
    TH2F *histo1=0;
    histo1 = new TH2F("NewOscProb","NewOscProb",100, 0., 0.09, 100, 0.0, 0.09);
	
	histo1->SetTitle("");
    //histo1->GetYaxis()->SetRangeUser(0.0,0.09);
    histo1->GetXaxis()->SetTitle("P(#nu_{#mu}#rightarrow #nu_{e})");
    histo1->GetYaxis()->SetTitle("P(#bar{#nu}_{#mu}#rightarrow #bar{#nu}_{e})");
    histo1->GetXaxis()->CenterTitle();
    histo1->GetYaxis()->CenterTitle();
    histo1->GetXaxis()->SetLabelSize(histo1->GetXaxis()->GetTitleSize()*1.2);
    histo1->GetYaxis()->SetLabelSize(histo1->GetYaxis()->GetTitleSize()*1.2);
    histo1->GetXaxis()->SetTitleSize(histo1->GetXaxis()->GetLabelSize()*1.2);
    histo1->GetYaxis()->SetTitleSize(histo1->GetYaxis()->GetLabelSize()*1.2);
    histo1->GetYaxis()->SetTitleOffset(1.0);
    
    canvas1->cd();
    gStyle->SetOptStat(0);
    gStyle->SetLabelSize(0.03,"xy");
    gStyle->SetTitleSize(0.04,"xy");
     //gStyle->SetTitleOffset(1.2,"xy");
    
    TLegend *leg = new TLegend(0.18,0.20,0.38,0.40);
    leg->SetFillStyle(0);
	leg->SetBorderSize(0);
    leg->AddEntry(MuEleccp0,"#delta_{CP}=0","p");
    leg->AddEntry(MuEleccp12,"#delta_{CP}=#pi/2","p");
    leg->AddEntry(MuEleccp1,"#delta_{CP}=#pi","p");
    leg->AddEntry(MuEleccp32,"#delta_{CP}=3#pi/2","p");
    
    histo1->Draw();
    MuElec->Draw("same");
    MuElecInv->Draw("same");
    MuElechighOct->Draw("same");
    MuElechighOctInv->Draw("same");
		
    MuEleccp0->Draw("same");
    MuElecInvcp0->Draw("same");
    MuEleccp12->Draw("same");
    MuElecInvcp12->Draw("same");
    MuEleccp1->Draw("same");
    MuElecInvcp1->Draw("same");
    MuEleccp32->Draw("same");
    MuElecInvcp32->Draw("same");
    MuElechighOctcp0->Draw("same");
    MuElechighOctInvcp0->Draw("same");
    MuElechighOctcp12->Draw("same");
    MuElechighOctInvcp12->Draw("same");
    MuElechighOctcp1->Draw("same");
    MuElechighOctInvcp1->Draw("same");
    MuElechighOctcp32->Draw("same");
    MuElechighOctInvcp32->Draw("same");
		
    tlx->Draw("same");
    tlx1->Draw("same");
    tlx2->Draw("same");
    tlx3->Draw("same");
    tlx4->Draw("same");
    tlx5->Draw("same");
    tlx6->Draw("same");
    tlx7->Draw("same");
    tlx8->Draw("same");
    leg->Draw("same");
    canvas1->Print("../plots/elec_appe"+fprob+LStrname+"octant.eps");
    canvas1->Print("../plots/elec_appe"+fprob+LStrname+"octant.pdf");
    canvas1->Print("../plots/elec_appe"+fprob+LStrname+"octant.png");
    
    exit(0);
    
}
