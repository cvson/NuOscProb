{

    #include<TH1>
    #include<TH2>
    #include<TLatex>
    #include<TString>
    
    gROOT->ProcessLine(".L ../source/OscCalculatorPMNS.C+");
  
    double x=0.0; //Energy
    bool isNeutrino = false;
  
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
  
    //input parameters
    OscCalculatorPMNS myProb;
    double myDeltam2=2.42;
    double myDeltam21=7.54;
    double myTheta12=0.61;
    double myTheta23=0.78;
    double myTheta13=0.15;
    double myL=734.;
  
    //Set all parameters here
    myProb.SetL(myL);
  
    myProb.SetRho(2.71);
    myProb.SetDmsq21(myDeltam21*TMath::Power(10,-5.));
    myProb.SetDmsq32(myDeltam2*TMath::Power(10,-3.));

    myProb.SetTh12(myTheta12);
    myProb.SetTh23(myTheta23);
    myProb.SetTh13(myTheta13);
    myProb.SetdCP(0.0001);

    //Flag to trigger antineutrino probabilities +1=neutrinos, -1=antineutrinos
    int isAnti=1;

    //Note, probability flavors are defined by pdg code:
    //nue  =12, antinue  =-12
    //numu =14, antinue  =-14
    //nutau=16, antinutau=-16
    //For instance, to compute P(numu->nue) you would call:
    // myProb.P(14,12,E)

    //Read back parameters to check
    std::cout << myProb.GetL() << endl
        << myProb.GetRho() << endl
	    << myProb.GetDmsq21() << endl  
	    << myProb.GetDmsq32() << endl  
	    << myProb.GetTh12() << endl  
	    << myProb.GetTh23() << endl  
	    << myProb.GetTh13() << endl  
	    << myProb.GetdCP() << endl  
	    << endl;

    //This step is done automatically when computing probabilities,
    //but useful to do it here to check code works
    myProb.SetMix(myProb.GetTh12(),myProb.GetTh23(),myProb.GetTh13(),myProb.GetdCP());
    myProb.SetDeltaMsqrs(myProb.GetDmsq21(),myProb.GetDmsq32());

    //Print mixing matrix and delta_m_sqrs.
    myProb.PrintMix();
    myProb.PrintDeltaMsqrs();

    //Keep this here for reference
    /*
     double par[9]={0.61,                       //Theta12
                   myTheta23,                  //Theta23                   
                   0.086,                       //Theta13
                   myDeltam2*TMath::Power(10,-3.), //DeltaM32^2
                   7.59*TMath::Power(10,-5.),  //DeltaM21^2
                   0.0001,                       //DeltaCP (3Pi/2)
		   2.71,                        //Rock Density  
                   myL,                        //L
		   1};                        //NuAntiNu
     */
    double epsilon=0.00001;
  

    TGraph *MuSurvive=0;
    TGraph *Mu2e =0;
    TGraph *Mu2tau=0;
    

    
    Int_t npoints=50000;
    Double_t loveremax =38000;
    
    TH1 *h1MuSurvive= new TH1F("","",npoints,0,loveremax);
    TH1 *h1Mu2e =new TH1F("","",npoints,0,loveremax);
    TH1 *h1Mu2tau=new TH1F("","",npoints,0,loveremax);
    
    Double_t musx[npoints], musy[npoints];
    Double_t musynue[npoints];
    Double_t musynutau[npoints];
    
    for (Int_t i=1; i < npoints+1; i++){
        //x=i/(npoints/10.) + epsilon;
        x=2.0 + epsilon;
        myL = (i/(npoints*1.0/loveremax)+epsilon)*x;
      
        myProb.Reset();
        myProb.SetL(myL);
      
        myProb.SetRho(2.71);
        myProb.SetDmsq21(myDeltam21*TMath::Power(10,-5.));
        myProb.SetDmsq32(myDeltam2*TMath::Power(10,-3.));
      
        myProb.SetTh12(myTheta12);
        myProb.SetTh23(myTheta23);
        myProb.SetTh13(myTheta13);
        myProb.SetdCP(0.0001);
        myProb.SetMix(myProb.GetTh12(),myProb.GetTh23(),myProb.GetTh13(),myProb.GetdCP());
        myProb.SetDeltaMsqrs(myProb.GetDmsq21(),myProb.GetDmsq32());
      
        musx[i-1]=i/(npoints*1.0/loveremax) + epsilon;
        musy[i-1]=myProb.P(isAnti*14,isAnti*14,x);
        musynue[i-1]=myProb.P(isAnti*14,isAnti*12,x);
        musynutau[i-1]=myProb.P(isAnti*14,isAnti*16,x);
        h1MuSurvive->SetBinContent(i,musy[i-1]);
        /*
         cout << "E = " << x << " : "
         << "; MuSurvive =  " << myProb.P(isAnti*14,isAnti*14,x)
         << endl;
         */
    }
    
    MuSurvive = new TGraph(npoints,musx,musy);
    Mu2e = new TGraph(npoints,musx,musynue);
    Mu2tau = new TGraph(npoints,musx,musynutau);
    
  
  
    //Plot everything

    TH2F *histo1=0;
    histo1 = new TH2F("NewOscProb","NewOscProb",400000, 0., 40000, 500, 0, 1.);
  
    TCanvas *canvas1=0;
    canvas1 = new TCanvas("Can1","Can1",1000,800);
  

  
    
    canvas1->cd();
    gStyle->SetOptStat(0);
    gStyle->SetLineWidth(2);
    //gPad->SetLeftMargin(gPad->GetLeftMargin()*1.4);
    double xcoord = 0.58;
    TString delMStr32=Form("%2.2f #times 10^{-3} eV^{2}",myDeltam2);
    TString legStr32="|#Delta m^{2}_{32}|=";
    TLatex* tlx=new TLatex(xcoord, 0.85,legStr32+delMStr32);
    tlx->SetNDC(kTRUE); // <- use NDC coordinate
    tlx->SetTextSize(0.04);
    tlx->SetTextAlign(12);
	
	TString delMStr21=Form("%2.2f #times 10^{-5} eV^{2}",myDeltam21);
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
    TLatex* tlx5=new TLatex(xcoord, 0.60,LStr);
    tlx5->SetNDC(kTRUE); // <- use NDC coordinate
    tlx5->SetTextSize(0.04);
    tlx5->SetTextAlign(12);
  
  
    histo1->SetTitle("");
    //histo1->GetYaxis()->SetRangeUser(-0.2,0.2);
    //  histo1->GetYaxis()->SetMinimum(0.0);
    histo1->GetXaxis()->SetTitle("L/E (km/GeV)");
    histo1->GetYaxis()->SetTitle("Probability");
    histo1->GetYaxis()->CenterTitle();
    histo1->GetXaxis()->CenterTitle();
    histo1->GetXaxis()->SetLabelSize(histo1->GetXaxis()->GetTitleSize()*1.2);
    histo1->GetYaxis()->SetLabelSize(histo1->GetYaxis()->GetTitleSize()*1.2);
    histo1->GetXaxis()->SetTitleSize(histo1->GetXaxis()->GetLabelSize()*1.2);
    histo1->GetYaxis()->SetTitleSize(histo1->GetYaxis()->GetLabelSize()*1.2);
    histo1->GetYaxis()->SetTitleOffset(1.0);
    histo1->Draw();
    MuSurvive->SetLineWidth(3);
    Mu2e->SetLineWidth(3);
    Mu2tau->SetLineWidth(3);
    Int_t ci;
    //brow
    ci = TColor::GetColor("#CC3300");
    MuSurvive->SetLineColor(ci);
    MuSurvive->Draw("C same");
    TGraphSmooth *gs = new TGraphSmooth("normal");
    TGraph *grout;
    grout = gs->SmoothKern(MuSurvive,"normal",0.2);//SmoothSuper(MuSurvive,"",0,0);
    //grout->DrawClone("L");
    //TSpline3 *s3 = new TSpline3("s3",MuSurvive->GetX(),MuSurvive->GetY(),MuSurvive->GetN());
    //s3->SetLineColor(kRed);
    //s3->Draw("l same");
    // h1MuSurvive->Smooth();
    //h1MuSurvive->DrawClone("same");
    
    ci = TColor::GetColor("#0033CC");
    Mu2tau->SetLineColor(ci);
    Mu2tau->Draw("C same");
    
    ci = TColor::GetColor("#333333");
    Mu2e->SetLineColor(ci);
    Mu2e->Draw("C same");
    
    TLine *zero = new TLine(0,0,10,0);
    zero->SetLineWidth(2);
    zero->SetLineStyle(7);
    // dim gray / dim gre
    ci = TColor::GetColor("#696969");
    zero->SetLineColor(ci);
    zero->Draw("same");
    //fmodel2flav_vac->Draw("same");
    /*tlx->Draw("same");
    tlx1->Draw("same");
    tlx2->Draw("same");
    tlx3->Draw("same");
    tlx4->Draw("same");
    tlx5->Draw("same");*/
    Double_t xlegmin = 0.35;
    Double_t ylegmin = 0.65;
    TLegend* leg = new TLegend(xlegmin,ylegmin,xlegmin+0.3,ylegmin+0.23);
    leg->SetFillStyle(0);
	leg->SetBorderSize(0);
    leg->SetTextSize(0.06);
    leg->SetTextFont(42);
    leg->AddEntry(MuSurvive, "P(#nu_{#mu}#rightarrow #nu_{#mu})","l");
    leg->AddEntry(Mu2e, "P(#nu_{#mu}#rightarrow #nu_{e})","l");
    leg->AddEntry(Mu2tau, "P(#nu_{#mu}#rightarrow #nu_{#tau})","l");
    leg->Draw();
    canvas1->SaveAs("../plots/threeflav_exact_probability_matter.eps");
    
    //TCanvas *canvas2=0;
    //canvas2 = new TCanvas("Can1","Can1",800,600);
    
    

    
        

}
