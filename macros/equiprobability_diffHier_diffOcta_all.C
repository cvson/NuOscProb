//////////////////////////////////////////////////////////
///// Testing NovA numu2nue probability
//// check equiprobability different hierarchy & different octant & nuAntinu
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
    //////////////////////////////////////////////////////////////////////////////////////
    //input neutrino normal hierarchy lower octant L=810
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
    
    //input neutrino inverted hierarchy lower octant L=810
    double parInputL2[9]={myTheta12,                       //Theta12
        myTheta23,                  //Theta23
        myTheta13,                       //Theta13
        -1.0*myDeltam32*TMath::Power(10,-3.), //DeltaM32^2
        myDeltam21*TMath::Power(10,-5.),  //DeltaM21^2
        mycp,                       //DeltaCP (3Pi/2)
        2.71,                        //Rock Density
        myL,                        //L
        1};
    SetOscParam(parInputL2);
   // x[0]=2.0+epsilon;
    double problevelL2=model3flav_mat_kopp_2nd_MuToElec(x);
    
    //input neutrino normal hierarchy higher octant L=810
    double parInputhigh[9]={myTheta12,                       //Theta12
        myTheta23highOc,                  //Theta23
        myTheta13,                       //Theta13
        myDeltam32*TMath::Power(10,-3.), //DeltaM32^2
        myDeltam21*TMath::Power(10,-5.),  //DeltaM21^2
        mycp,                       //DeltaCP (3Pi/2)
        2.71,                        //Rock Density
        myL,                        //L
        1};
    SetOscParam(parInputhigh);
    //x[0]=2.0+epsilon;
    double problevelhigh=model3flav_mat_kopp_2nd_MuToElec(x);
    
    //input neutrino inverted hierarchy higher octant L=810
    double parInputL2high[9]={myTheta12,                       //Theta12
        myTheta23highOc,                  //Theta23
        myTheta13,                       //Theta13
        -1.0*myDeltam32*TMath::Power(10,-3.), //DeltaM32^2
        myDeltam21*TMath::Power(10,-5.),  //DeltaM21^2
        mycp,                       //DeltaCP (3Pi/2)
        2.71,                        //Rock Density
        myL,                        //L
        1};
    SetOscParam(parInputL2high);
    // x[0]=2.0+epsilon;
    double problevelL2high=model3flav_mat_kopp_2nd_MuToElec(x);
    
    
    //////////////////////////////////////////////////////////////////////////////////////
    //input Antineutrino normal hierarchy lower octant L=810
    double parInputAnti[9]={myTheta12,                       //Theta12
        myTheta23,                  //Theta23
        myTheta13,                       //Theta13
        myDeltam32*TMath::Power(10,-3.), //DeltaM32^2
        myDeltam21*TMath::Power(10,-5.),  //DeltaM21^2
        mycp,                       //DeltaCP (3Pi/2)
        2.71,                        //Rock Density
        myL,                        //L
        -1};
    SetOscParam(parInputAnti);
    //x[0]=2.0+epsilon;
    double problevelAnti=model3flav_mat_kopp_2nd_MuToElec(x);
    
    //input Antineutrino inverted hierarchy lower octant L=810
    double parInputAntiL2[9]={myTheta12,                       //Theta12
        myTheta23,                  //Theta23
        myTheta13,                       //Theta13
        -1.0*myDeltam32*TMath::Power(10,-3.), //DeltaM32^2
        myDeltam21*TMath::Power(10,-5.),  //DeltaM21^2
        mycp,                       //DeltaCP (3Pi/2)
        2.71,                        //Rock Density
        myL,                        //L
        -1};
    SetOscParam(parInputAntiL2);
    // x[0]=2.0+epsilon;
    double problevelAntiL2=model3flav_mat_kopp_2nd_MuToElec(x);
    
    //input Antineutrino normal hierarchy higher octant L=810
    double parInputAntihigh[9]={myTheta12,                       //Theta12
        myTheta23highOc,                  //Theta23
        myTheta13,                       //Theta13
        myDeltam32*TMath::Power(10,-3.), //DeltaM32^2
        myDeltam21*TMath::Power(10,-5.),  //DeltaM21^2
        mycp,                       //DeltaCP (3Pi/2)
        2.71,                        //Rock Density
        myL,                        //L
        -1};
    SetOscParam(parInputAntihigh);
    //x[0]=2.0+epsilon;
    double problevelAntihigh=model3flav_mat_kopp_2nd_MuToElec(x);
    
    //input Antineutrino inverted hierarchy higher octant L=810
    double parInputAntiL2high[9]={myTheta12,                       //Theta12
        myTheta23highOc,                  //Theta23
        myTheta13,                       //Theta13
        -1.0*myDeltam32*TMath::Power(10,-3.), //DeltaM32^2
        myDeltam21*TMath::Power(10,-5.),  //DeltaM21^2
        mycp,                       //DeltaCP (3Pi/2)
        2.71,                        //Rock Density
        myL,                        //L
        -1};
    SetOscParam(parInputAntiL2high);
    // x[0]=2.0+epsilon;
    double problevelAntiL2high=model3flav_mat_kopp_2nd_MuToElec(x);
    //////////////////////////////////////////////////////////////////////////////////////
    
    
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
    
    /*TString ProbStr=Form("P(#nu_{#mu}#rightarrow #nu_{e})=%2.2f",problevel);
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
    tlx8->SetTextSize(0.04);*/
    
    
    //TString LStrname=Form("%2.0f",myL);
    

    
    
    //neutrino normal hierarchy lower octant
    TH2F *hequiprob=0;
    TH1F *h1equiprob=new TH1F("h1","h1",npointsx,xmin,xmax);
    hequiprob = new TH2F("hcont_nu_norm","hcont_nu_norm",npointsx, xmin, xmax, npointsy, ymin, ymax);
    hequiprob->GetXaxis()->SetTitle("Mixing angle #theta_{13}");
    hequiprob->GetYaxis()->SetTitle("CP phase");
    hequiprob->GetXaxis()->CenterTitle();
    hequiprob->GetYaxis()->CenterTitle();
	hequiprob->GetXaxis()->SetTitleSize(0.05);
    hequiprob->GetYaxis()->SetTitleSize(0.05);
	hequiprob->GetYaxis()->SetTitleOffset(1.0);
    hequiprob->GetXaxis()->SetTitleOffset(0.9);
    
    //neutrino inverted hierarchy lower octant
    TH2F *hequiprobInv=0;
    TH1F *h2equiprob=new TH1F("h2","h2",npointsx,xmin,xmax);
    hequiprobInv = new TH2F("hcont_nu_normInv","hcont_nu_normInv",npointsx, xmin, xmax, npointsy, ymin, ymax);
    
    //neutrino normal hierarchy higher octant
    TH2F *hequiprobhigh=0;
    TH1F *h3equiprob=new TH1F("h3","h3",npointsx,xmin,xmax);
    hequiprobhigh = new TH2F("hcont_nu_normhight","hcont_nu_normhigh",npointsx, xmin, xmax, npointsy, ymin, ymax);
    
    //neutrino inverted hierarchy lower octant
    TH2F *hequiprobInvhigh=0;
    TH1F *h4equiprob=new TH1F("h4","h4",npointsx,xmin,xmax);
    hequiprobInvhigh = new TH2F("hcont_nu_normInvhigh","hcont_nu_normInvhigh",npointsx, xmin, xmax, npointsy, ymin, ymax);
    
    /////////////////////////Antineutrino
    //Antineutrino normal hierarchy lower octant
    TH2F *hAntiequiprob=0;
    TH1F *hAnti1equiprob=new TH1F("hAnti1","hAnti1",npointsx,xmin,xmax);
    hAntiequiprob = new TH2F("hAnticont_nu_norm","hAnticont_nu_norm",npointsx, xmin, xmax, npointsy, ymin, ymax);
    
    //Antineutrino inverted hierarchy lower octant
    TH2F *hAntiequiprobInv=0;
    TH1F *hAnti2equiprob=new TH1F("hAnti2","hAnti2",npointsx,xmin,xmax);
    hAntiequiprobInv = new TH2F("hAnticont_nu_normInv","hAnticont_nu_normInv",npointsx, xmin, xmax, npointsy, ymin, ymax);
    
    //Antineutrino normal hierarchy higher octant
    TH2F *hAntiequiprobhigh=0;
    TH1F *hAnti3equiprob=new TH1F("hAnti3","hAnti3",npointsx,xmin,xmax);
    hAntiequiprobhigh = new TH2F("hAnticont_nu_normhight","hAnticont_nu_normhigh",npointsx, xmin, xmax, npointsy, ymin, ymax);
    
    //neutrino inverted hierarchy lower octant
    TH2F *hAntiequiprobInvhigh=0;
    TH1F *hAnti4equiprob=new TH1F("hAnti4","hAnti4",npointsx,xmin,xmax);
    hAntiequiprobInvhigh = new TH2F("hAnticont_nu_normInvhigh","hAnticont_nu_normInvhigh",npointsx, xmin, xmax, npointsy, ymin, ymax);
    
    ///////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////
    //neutrino w L = 810km normal hierarchy lower octant
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
    
    //neutrino w L =810km inverted hierarchy lower octant
	for (Int_t ix=0; ix < npointsx+1; ix++){
		double ptheta13 = ix*(xmax-xmin)/npointsx;
        for (Int_t iy=0; iy < npointsy+1; iy++){
			double pcp = iy*(ymax-ymin)/npointsy;
			
			double par[9]={myTheta12,                       //Theta12
				myTheta23,                  //Theta23
				ptheta13,                       //Theta13
				-1.0*myDeltam32*TMath::Power(10,-3.), //DeltaM32^2
				myDeltam21*TMath::Power(10,-5.),  //DeltaM21^2
				pcp,                       //DeltaCP (3Pi/2)
				2.71,                        //Rock Density
				myL,                        //L
				1};
			SetOscParam(par);
			x[0]=2.0+epsilon;
			double pprob=model3flav_mat_kopp_2nd_MuToElec(x);
			cout<<"pmuex["<<ix<<","<<iy<<"] = "<<pprob<<endl;
			
			hequiprobInv->SetBinContent(ix,iy,pprob);
			
		}
	}
    
    //neutrino w L = 810km normal hierarchy higher octant
	for (Int_t ix=0; ix < npointsx+1; ix++){
		double ptheta13 = ix*(xmax-xmin)/npointsx;
        for (Int_t iy=0; iy < npointsy+1; iy++){
			double pcp = iy*(ymax-ymin)/npointsy;
			
			double par[9]={myTheta12,                       //Theta12
				myTheta23highOc,                  //Theta23
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
			
			hequiprobhigh->SetBinContent(ix,iy,pprob);
			
		}
	}
    
    //neutrino w L =810km inverted hierarchy higher octant
	for (Int_t ix=0; ix < npointsx+1; ix++){
		double ptheta13 = ix*(xmax-xmin)/npointsx;
        for (Int_t iy=0; iy < npointsy+1; iy++){
			double pcp = iy*(ymax-ymin)/npointsy;
			
			double par[9]={myTheta12,                       //Theta12
				myTheta23highOc,                  //Theta23
				ptheta13,                       //Theta13
				-1.0*myDeltam32*TMath::Power(10,-3.), //DeltaM32^2
				myDeltam21*TMath::Power(10,-5.),  //DeltaM21^2
				pcp,                       //DeltaCP (3Pi/2)
				2.71,                        //Rock Density
				myL,                        //L
				1};
			SetOscParam(par);
			x[0]=2.0+epsilon;
			double pprob=model3flav_mat_kopp_2nd_MuToElec(x);
			cout<<"pmuex["<<ix<<","<<iy<<"] = "<<pprob<<endl;
			
			hequiprobInvhigh->SetBinContent(ix,iy,pprob);
			
		}
	}
    ///////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////
    //Antineutrino w L = 810km normal hierarchy lower octant
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
				-1};
			SetOscParam(par);
			x[0]=2.0+epsilon;
			double pprob=model3flav_mat_kopp_2nd_MuToElec(x);
			cout<<"pmuex["<<ix<<","<<iy<<"] = "<<pprob<<endl;
			
			hAntiequiprob->SetBinContent(ix,iy,pprob);
			
		}
	}
    
    //Antineutrino w L =810km inverted hierarchy lower octant
	for (Int_t ix=0; ix < npointsx+1; ix++){
		double ptheta13 = ix*(xmax-xmin)/npointsx;
        for (Int_t iy=0; iy < npointsy+1; iy++){
			double pcp = iy*(ymax-ymin)/npointsy;
			
			double par[9]={myTheta12,                       //Theta12
				myTheta23,                  //Theta23
				ptheta13,                       //Theta13
				-1.0*myDeltam32*TMath::Power(10,-3.), //DeltaM32^2
				myDeltam21*TMath::Power(10,-5.),  //DeltaM21^2
				pcp,                       //DeltaCP (3Pi/2)
				2.71,                        //Rock Density
				myL,                        //L
				-1};
			SetOscParam(par);
			x[0]=2.0+epsilon;
			double pprob=model3flav_mat_kopp_2nd_MuToElec(x);
			cout<<"pmuex["<<ix<<","<<iy<<"] = "<<pprob<<endl;
			
			hAntiequiprobInv->SetBinContent(ix,iy,pprob);
			
		}
	}
    
    //Antineutrino w L = 810km normal hierarchy higher octant
	for (Int_t ix=0; ix < npointsx+1; ix++){
		double ptheta13 = ix*(xmax-xmin)/npointsx;
        for (Int_t iy=0; iy < npointsy+1; iy++){
			double pcp = iy*(ymax-ymin)/npointsy;
			
			double par[9]={myTheta12,                       //Theta12
				myTheta23highOc,                  //Theta23
				ptheta13,                       //Theta13
				myDeltam32*TMath::Power(10,-3.), //DeltaM32^2
				myDeltam21*TMath::Power(10,-5.),  //DeltaM21^2
				pcp,                       //DeltaCP (3Pi/2)
				2.71,                        //Rock Density
				myL,                        //L
				-1};
			SetOscParam(par);
			x[0]=2.0+epsilon;
			double pprob=model3flav_mat_kopp_2nd_MuToElec(x);
			cout<<"pmuex["<<ix<<","<<iy<<"] = "<<pprob<<endl;
			
			hAntiequiprobhigh->SetBinContent(ix,iy,pprob);
			
		}
	}
    
    //Antineutrino w L =810km inverted hierarchy higher octant
	for (Int_t ix=0; ix < npointsx+1; ix++){
		double ptheta13 = ix*(xmax-xmin)/npointsx;
        for (Int_t iy=0; iy < npointsy+1; iy++){
			double pcp = iy*(ymax-ymin)/npointsy;
			
			double par[9]={myTheta12,                       //Theta12
				myTheta23highOc,                  //Theta23
				ptheta13,                       //Theta13
				-1.0*myDeltam32*TMath::Power(10,-3.), //DeltaM32^2
				myDeltam21*TMath::Power(10,-5.),  //DeltaM21^2
				pcp,                       //DeltaCP (3Pi/2)
				2.71,                        //Rock Density
				myL,                        //L
				-1};
			SetOscParam(par);
			x[0]=2.0+epsilon;
			double pprob=model3flav_mat_kopp_2nd_MuToElec(x);
			cout<<"pmuex["<<ix<<","<<iy<<"] = "<<pprob<<endl;
			
			hAntiequiprobInvhigh->SetBinContent(ix,iy,pprob);
			
		}
	}
    
    //fake histogram
    //UT orange
    Int_t ci1 =  TColor::GetColor("#B45F04");
    //navy blue
    Int_t ci2 = TColor::GetColor("#006699");
    h1equiprob->SetLineColor(ci1);
    h2equiprob->SetLineColor(ci2);
    
    h3equiprob->SetLineColor(ci1);
    h4equiprob->SetLineColor(ci2);
    h3equiprob->SetLineStyle(7);
    h4equiprob->SetLineStyle(7);
    
    //UT orange
    Int_t ci3 =  TColor::GetColor("#CC3300");
    //navy blue
    Int_t ci4 = TColor::GetColor("#339900");
    hAnti1equiprob->SetLineColor(ci3);
    hAnti2equiprob->SetLineColor(ci4);
    
    hAnti3equiprob->SetLineColor(ci3);
    hAnti4equiprob->SetLineColor(ci4);
    hAnti3equiprob->SetLineStyle(7);
    hAnti4equiprob->SetLineStyle(7);
    
    TLegend *leg = new TLegend(xcoord,0.20,xcoord+0.24,0.51);
    leg->SetFillStyle(0);
	leg->SetBorderSize(0);
    leg->SetTextSize(0.04);
    
    leg->AddEntry(InputPoint,"Input parameter","p");
	leg->AddEntry(h1equiprob,"#nu Normal & #theta_{23}<45^{o}","l");
    leg->AddEntry(h2equiprob,"#nu Inverted & #theta_{23}<45^{o}","l");
    leg->AddEntry(h3equiprob,"#nu Normal & #theta_{23}>45^{o}","l");
    leg->AddEntry(h4equiprob,"#nu Inverted & #theta_{23}>45^{o}","l");
    leg->AddEntry(hAnti1equiprob,"#bar{#nu} Normal & #theta_{23}<45^{o}","l");
    leg->AddEntry(hAnti2equiprob,"#bar{#nu} Inverted & #theta_{23}<45^{o}","l");
    leg->AddEntry(hAnti3equiprob,"#bar{#nu} Normal & #theta_{23}>45^{o}","l");
    leg->AddEntry(hAnti4equiprob,"#bar{#nu} Inverted & #theta_{23}>45^{o}","l");
    
    //leg->AddEntry(MuElec_right,Th13StrR,"l");
    
	
    TCanvas *canvas1=0;
    canvas1 = new TCanvas("Can1","Can1",800,600);
    gStyle->SetOptStat(0);
    gStyle->SetLineWidth(2);
    canvas1->cd();
    Int_t colors[]={ci1,ci2,ci3,ci4};
    gStyle->SetPalette((sizeof(colors)/sizeof(Int_t)), colors);
    
    Double_t levels[]={problevel,1.0,2.0,3.0};
    hequiprob->SetContour((sizeof(levels)/sizeof(Double_t)), levels);
    hequiprob->SetTitle("");
    hequiprob->GetXaxis()->SetLabelSize(hequiprob->GetXaxis()->GetTitleSize()*1.2);
    hequiprob->GetYaxis()->SetLabelSize(hequiprob->GetYaxis()->GetTitleSize()*1.2);
    hequiprob->GetXaxis()->SetTitleSize(hequiprob->GetXaxis()->GetLabelSize()*1.2);
    hequiprob->GetYaxis()->SetTitleSize(hequiprob->GetYaxis()->GetLabelSize()*1.2);
    hequiprob->GetZaxis()->SetTitleSize(hequiprob->GetZaxis()->GetLabelSize()*1.2);
    hequiprob->GetYaxis()->SetTitleOffset(0.9);
    hequiprob->Draw("cont1");
    
    Double_t levelsInv[]={0.0,problevelL2,1.0,2.0};
    hequiprobInv->SetContour((sizeof(levelsInv)/sizeof(Double_t)), levelsInv);
    hequiprobInv->SetTitle("");
    hequiprobInv->Draw("cont1 same");
    
    Double_t levelshigh[]={problevelhigh,1.0,2.0,3.0};
    hequiprobhigh->SetContour((sizeof(levelshigh)/sizeof(Double_t)), levelshigh);
    hequiprobhigh->SetTitle("");
    hequiprobhigh->SetLineStyle(7);
    hequiprobhigh->Draw("cont1 same");
    
    Double_t levelsInvhigh[]={0.0,problevelL2high,2.0,3.0};
    hequiprobInvhigh->SetContour((sizeof(levelsInvhigh)/sizeof(Double_t)), levelsInvhigh);
    hequiprobInvhigh->SetTitle("");
    hequiprobInvhigh->SetLineStyle(7);
    hequiprobInvhigh->Draw("cont1 same");
    
    Double_t levelsAnti[]={-1.0,0.0,problevelAnti,1.0};
    hAntiequiprob->SetContour((sizeof(levelsAnti)/sizeof(Double_t)), levelsAnti);
    hAntiequiprob->SetTitle("");
    hAntiequiprob->Draw("cont1 same");
    
    Double_t levelsAntiInv[]={-2.0,-1.0,0.0,problevelAntiL2};
    hAntiequiprobInv->SetContour((sizeof(levelsAntiInv)/sizeof(Double_t)), levelsAntiInv);
    hAntiequiprobInv->SetTitle("");
    hAntiequiprobInv->Draw("cont1 same");
    
    Double_t levelsAntihigh[]={-1.0,0.0,problevelAntihigh,1.0};
    hAntiequiprobhigh->SetContour((sizeof(levelsAntihigh)/sizeof(Double_t)), levelsAntihigh);
    hAntiequiprobhigh->SetTitle("");
    hAntiequiprobhigh->SetLineStyle(7);
    hAntiequiprobhigh->Draw("cont1 same");
    
    Double_t levelsAntiInvhigh[]={-2.0,-1.0,0.0,problevelAntiL2high};
    hAntiequiprobInvhigh->SetContour((sizeof(levelsAntiInvhigh)/sizeof(Double_t)), levelsAntiInvhigh);
    hAntiequiprobInvhigh->SetTitle("");
    hAntiequiprobInvhigh->SetLineStyle(7);
    hAntiequiprobInvhigh->Draw("cont1 same");
    
    
    
    
    InputPoint->Draw("same");
    tlx->Draw("same");
    tlx1->Draw("same");
    tlx2->Draw("same");
    tlx3->Draw("same");
    tlx4->Draw("same");
    //tlxcp->Draw("same");
    tlx5->Draw("same");
    tlx6->Draw("same");
    //tlxprob->Draw("same");
    leg->Draw("same");
    canvas1->Update();
    canvas1->Print("../plots/equiprob_"+fprob+"_1cont_diffHier_diffOcta_all.eps");
    canvas1->Print("../plots/equiprob_"+fprob+"_1cont_diffHier_diffOcta_all.pdf");
    
    TFile *poutFile = new TFile("../plots/equiprob_"+fprob+"_1cont_diffHier_diffOcta_all.root","RECREATE");
    hequiprob->Write("hequiprob");
    hequiprobInv->Write("hequiprobInv");
    hequiprobhigh->Write("hequiprobhigh");
    hequiprobInvhigh->Write("hequiprobInvhigh");
    hAntiequiprob->Write("hAntiequiprob");
    hAntiequiprobInv->Write("hAntiequiprobInv");
    hAntiequiprobhigh->Write("hAntiequiprobhigh");
    hAntiequiprobInvhigh->Write("hAntiequiprobInvhigh");
    poutFile->Close();

    exit(0);
    
}
