////////////////////////////////////////////////////////////////
//compare mu disappearance/appearance from variety of approximation
//in matter
////////////////////////////////////////////////////////////////
{

#include<TH1>
#include<TH2>
#include<TLatex>
#include<TString>
    
    gROOT->ProcessLine(".L ../source/NuOscProb.C+");
    //energy variable
    double x[1];
    
    //channel appearance (numu-->nue) or disappearance (numu-->numu)
    bool isAppearance=true;
    
    //change probability scale
    double ymax=isAppearance? 0.15:1.2;//for MutoE
    
    //savename
    TString fprob="comp_mat_MuToElec";

    /*
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

    //oscillation parameters
    double myDeltam32=2.42;
    double myDeltam21 = 7.54;
    
    double myTheta12 = 0.61;
    double myTheta23 = 0.78;
    double myTheta13 = 0.15;
    
    double myRockDensity = 2.71;
    double myCP = 0.0001;
    
    // baseline
    double myL=1250;
    TString flength=Form("L%2.0fkm",myL);
    
    double epsilon=0.00001;
    Int_t npoints=500;
    // fmodel3flav_mat = model3flav_mat_MuToElec(x)
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //Up to the second order
    double par[9]={myTheta12,                       //Theta12 0.61
        myTheta23,                  //Theta23
        myTheta13,                       //Theta13
        myDeltam32*TMath::Power(10,-3.), //DeltaM32^2
        myDeltam21*TMath::Power(10,-5.),  //DeltaM21^2 7.54
        myCP,                       //DeltaCP (3Pi/2)
        myRockDensity,                        //Rock Density
        myL,                        //L
        1};                        //NuAntiNu
        
    SetOscParam(par);
    GetOscParam(par);
    
    TGraph *fmodel3flav_mat=0;
    Double_t xmodel3flav_mat[npoints], ymodel3flav_mat[npoints];
    
    for (Int_t i=1; i < npoints+1; i++){
        x[0]=i/(npoints/10.) + epsilon;
        xmodel3flav_mat[i-1]=i/(npoints/10.) + epsilon;
        ymodel3flav_mat[i-1]=model3flav_mat_MuToElec(x);
    }
    
    fmodel3flav_mat = new TGraph(npoints,xmodel3flav_mat,ymodel3flav_mat);
    
    // fmodel3flav_mat_2nd_kopp = model3flav_mat_kopp_2nd_MuToElec(x)
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //This is the second order term only
    
    TGraph *fmodel3flav_mat_2nd_kopp=0;
    Double_t xmodel3flav_mat_2nd_kopp[npoints], ymodel3flav_mat_2nd_kopp[npoints];
    
    for (Int_t i=1; i < npoints+1; i++){
        x[0]=i/(npoints/10.) + epsilon;
        xmodel3flav_mat_2nd_kopp[i-1]=i/(npoints/10.) + epsilon;
        ymodel3flav_mat_2nd_kopp[i-1]=model3flav_mat_kopp_2nd_MuToElec(x);
        
    }
    
    fmodel3flav_mat_2nd_kopp = new TGraph(npoints,xmodel3flav_mat_2nd_kopp,ymodel3flav_mat_2nd_kopp);
    
    //fmodel3flav_mat_1stalpha = model3flav_vac_kopp_MuToElec(x)
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //First term in alpha

    TGraph *fmodel3flav_mat_1stalpha=0;
    Double_t xmodel3flav_mat_1stalpha[npoints], ymodel3flav_mat_1stalpha[npoints];
    
    for (Int_t i=1; i < npoints+1; i++){
        x[0]=i/(npoints/10.) + epsilon;
        xmodel3flav_mat_1stalpha[i-1]=i/(npoints/10.) + epsilon;
        ymodel3flav_mat_1stalpha[i-1]=model3flav_mat_1stalpha_MuToElec(x);
    }
    
    fmodel3flav_mat_1stalpha = new TGraph(npoints,xmodel3flav_mat_1stalpha,ymodel3flav_mat_1stalpha);
    
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    //First term in theta13
    
    TGraph *fmodel3flav_mat_1stth13=0;
    Double_t xmodel3flav_mat_1stth13[npoints], ymodel3flav_mat_1stth13[npoints];
    
    for (Int_t i=1; i < npoints+1; i++){
        x[0]=i/(npoints/10.) + epsilon;
        xmodel3flav_mat_1stth13[i-1]=i/(npoints/10.) + epsilon;
        ymodel3flav_mat_1stth13[i-1]=model3flav_mat_1stsin13_MuToElec(x);
    }
    
    fmodel3flav_mat_1stth13 = new TGraph(npoints,xmodel3flav_mat_1stth13,ymodel3flav_mat_1stth13);
    
    
    
    //Plot here
    TH2F *histo1=0;
    histo1 = new TH2F("NewOscProb","NewOscProb",100, 0., 10., 100, -.2., 1.);


    TCanvas *canvas1=0;
    gStyle->SetLineWidth(2);
    canvas1 = new TCanvas("Can1","Can1",1000,800);
       
    
    fmodel3flav_mat->SetLineColor(kBlue);
    fmodel3flav_mat->SetLineWidth(2);
    
    fmodel3flav_mat_2nd_kopp->SetLineColor(kRed);
    fmodel3flav_mat_2nd_kopp->SetLineWidth(2);
    
    fmodel3flav_mat_1stalpha->SetLineColor(kBlack);
    fmodel3flav_mat_1stalpha->SetLineWidth(2);
    
    fmodel3flav_mat_1stth13->SetLineColor(95);
    fmodel3flav_mat_1stth13->SetLineWidth(2);

    canvas1->cd();
    gStyle->SetOptStat(0);
    canvas1->SetLogx(1);
    TString delMStr=Form("%2.2f #times 10^{-3} eV^{2}",myDeltam32);

    TString legStr="#Delta m^{2}_{32}~";
    
    TLatex* tlx=new TLatex(0.5, 0.7,legStr+delMStr);
    tlx->SetNDC(kTRUE); // <- use NDC coordinate
    tlx->SetTextSize(0.04);

    TString ThStr=Form("sin^{2}2#theta_{23}=%2.2f",TMath::Sin(2*myTheta23)*TMath::Sin(2*myTheta23));

    TLatex* tlx1=new TLatex(0.5, 0.6,ThStr);
    tlx1->SetNDC(kTRUE); // <- use NDC coordinate
    tlx1->SetTextSize(0.04);

    TString LStr=Form("L=%2.0f",myL);

    TLatex* tlx2=new TLatex(0.5, 0.5,LStr);
    tlx2->SetNDC(kTRUE); // <- use NDC coordinate
    tlx2->SetTextSize(0.04);
    
    TLegend *leg = new TLegend(0.5,0.48,0.85,0.88);
    leg->SetFillStyle(0);
	leg->SetBorderSize(0);
    leg->SetTextSize(0.05);
    leg->SetTextFont(42);
    leg->SetHeader("Models");
    histo1->SetTitle("");
    histo1->GetYaxis()->SetRangeUser(0.0,ymax);
    //  histo1->GetYaxis()->SetMinimum(0.0);
    histo1->GetXaxis()->SetTitle("True Neutrino Energy (GeV)");
    histo1->GetYaxis()->SetTitle("Oscillation Prob.");
    histo1->GetXaxis()->CenterTitle();
    histo1->GetYaxis()->CenterTitle();
    histo1->GetYaxis()->SetTitleOffset(1.3);
    histo1->Draw();
    fmodel3flav_mat_1stalpha->Draw("same");
    fmodel3flav_mat->Draw("same");
    fmodel3flav_mat_2nd_kopp->Draw("same");
    fmodel3flav_mat_1stth13->Draw("same");
    leg->AddEntry(fmodel3flav_mat_1stth13,"3flav mat 1st th13","l");
    leg->AddEntry(fmodel3flav_mat_1stalpha,"3flav mat 1st alpha","l");
    leg->AddEntry(fmodel3flav_mat_2nd_kopp,"3flav mat 2nd only","l");
    leg->AddEntry(fmodel3flav_mat,"3flav mat up to 2nd","l");
    leg->Draw();


    canvas1->SaveAs("../plots/"+fprob+"_"+flength+".gif");
    canvas1->SaveAs("../plots/"+fprob+"_"+flength+".eps");
    
    exit(0);

  }
