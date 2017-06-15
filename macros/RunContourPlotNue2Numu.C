{

#include<TH1>
#include<TH2>
#include<TLatex>
#include<TString>
#include<TMath>
#include<TFile>    
    
    gROOT->ProcessLine(".L ../source/NuOscProb.C+");
    
    //energy variable
    double x[1];
    
    //channel appearance (numu-->nue) or disappearance (numu-->numu)
    //bool isAppearance=true;
    
    //change probability scale
    //double ymax=isAppearance? 0.15:1.2;
    
    //savename
    TString fcont="cont_";
    TString fprob="model3flav_mat_ElecToMu";
    TString fsubname="../plots/";


    //oscillation parameters
    double myDeltam32 = 2.42;
    double myDeltam21 = 7.54;
    
    double myTheta12 = 0.61;
    double myTheta23 = 0.78;
    double myTheta13 = 0.15;
    
    double myRockDensity = 2.71;
    double myCP = 0.0001;
    
    // baseline
    double myL;
    
    /*double par[9]={myTheta12,                       //Theta12
                   myTheta23,                  //Theta23
                   myTheta13,                       //Theta13
                   myDeltam32*TMath::Power(10,-3.), //DeltaM32^2
                   myDeltam21*TMath::Power(10,-5.),  //DeltaM21^2
                   myCP,                       //DeltaCP (3Pi/2)
		   myRockDensity,                        //Rock Density  
                   myL,                        //L
		   1};                        //NuAntiNu

    double epsilon=0.00001;

    SetOscParam(par);
    GetOscParam(par);*/
    
    //Neutrino Normal hierarchy
    Int_t npointsx=500;
    Int_t npointsy=500;
    //baseline length
    Double_t ymin =1e2;
    Double_t ymax = 1e5;
    
    //energy value
    Double_t xmin =1e-2;
    Double_t xmax =1e2;
    double epsilon=0.00001;
    Double_t binx[npointsx];
    Double_t biny[npointsy];
    for (Int_t i=0; i < npointsx+1; i++){
        binx[i]=TMath::Power(10.,(TMath::Log10(xmin)+((i+0.5)*1.0/npointsx)*(TMath::Log10(xmax)-TMath::Log10(xmin))));
        cout<<"binx "<<binx[i]<<endl;
    }
    
    for (Int_t i=0; i < npointsy+1; i++){
        biny[i] = TMath::Power(10.,(TMath::Log10(ymin)+((i+0.5)*1.0/npointsy)*(TMath::Log10(ymax)-TMath::Log10(ymin))));
        cout<<"biny "<<biny[i]<<endl;
    }
    
    //neutrino normal hierarchy
    TH2F *hcont_nu_norm=0;
    hcont_nu_norm = new TH2F("hcont_nu_norm","hcont_nu_norm",npointsx, binx, npointsy, biny);
    hcont_nu_norm->GetXaxis()->SetTitle("True Neutrino Energy (GeV)");
    hcont_nu_norm->GetYaxis()->SetTitle("Baseline (km)");
    hcont_nu_norm->GetXaxis()->CenterTitle();
    hcont_nu_norm->GetYaxis()->CenterTitle();

    for (Int_t ix=0; ix < npointsx+1; ix++){
        for (Int_t iy=0; iy < npointsy+1; iy++){
            x[0]= TMath::Power(10.,(TMath::Log10(xmin)+((ix+0.5)*1.0/npointsx)*(TMath::Log10(xmax)-TMath::Log10(xmin))));
            myL = TMath::Power(10.,(TMath::Log10(ymin)+((iy+0.5)*1.0/npointsy)*(TMath::Log10(ymax)-TMath::Log10(ymin))));
            double par[9]={myTheta12,                       //Theta12
                myTheta23,                  //Theta23
                myTheta13,                       //Theta13
                myDeltam32*TMath::Power(10,-3.), //DeltaM32^2
                myDeltam21*TMath::Power(10,-5.),  //DeltaM21^2
                myCP,                       //DeltaCP (3Pi/2)
                myRockDensity,                        //Rock Density
                myL,                        //L
                1};                        //NuAntiNu
            
            SetOscParam(par);
            hcont_nu_norm->SetBinContent(ix,iy, model3flav_mat_ElecToMu(x));
        }
    }
    
    //neutrino invert hierarchy
    TH2F *hcont_nu_invt=0;
    hcont_nu_invt = new TH2F("hcont_nu_invt","hcont_nu_invt",npointsx, binx, npointsy, biny);
    hcont_nu_invt->GetXaxis()->SetTitle("True Neutrino Energy (GeV)");
    hcont_nu_invt->GetYaxis()->SetTitle("Baseline (km)");
    hcont_nu_invt->GetXaxis()->CenterTitle();
    hcont_nu_invt->GetYaxis()->CenterTitle();
    
    for (Int_t ix=0; ix < npointsx+1; ix++){
        for (Int_t iy=0; iy < npointsy+1; iy++){
            x[0]= TMath::Power(10.,(TMath::Log10(xmin)+((ix+0.5)*1.0/npointsx)*(TMath::Log10(xmax)-TMath::Log10(xmin))));
            myL = TMath::Power(10.,(TMath::Log10(ymin)+((iy+0.5)*1.0/npointsy)*(TMath::Log10(ymax)-TMath::Log10(ymin))));
            double par[9]={myTheta12,                       //Theta12
                myTheta23,                  //Theta23
                myTheta13,                       //Theta13
                -1.0*myDeltam32*TMath::Power(10,-3.), //DeltaM32^2
                myDeltam21*TMath::Power(10,-5.),  //DeltaM21^2
                myCP,                       //DeltaCP (3Pi/2)
                myRockDensity,                        //Rock Density
                myL,                        //L
                1};                        //NuAntiNu
            
            SetOscParam(par);
            hcont_nu_invt->SetBinContent(ix,iy, model3flav_mat_ElecToMu(x));
        }
    }
    
    //antineutrino normal hierarchy
    TH2F *hcont_antinu_norm=0;
    hcont_antinu_norm = new TH2F("hcont_antinu_norm","hcont_antinu_norm",npointsx, binx, npointsy, biny);
    hcont_antinu_norm->GetXaxis()->SetTitle("True Neutrino Energy (GeV)");
    hcont_antinu_norm->GetYaxis()->SetTitle("Baseline (km)");
    hcont_antinu_norm->GetXaxis()->CenterTitle();
    hcont_antinu_norm->GetYaxis()->CenterTitle();
    
    for (Int_t ix=0; ix < npointsx+1; ix++){
        for (Int_t iy=0; iy < npointsy+1; iy++){
            x[0]= TMath::Power(10.,(TMath::Log10(xmin)+((ix+0.5)*1.0/npointsx)*(TMath::Log10(xmax)-TMath::Log10(xmin))));
            myL = TMath::Power(10.,(TMath::Log10(ymin)+((iy+0.5)*1.0/npointsy)*(TMath::Log10(ymax)-TMath::Log10(ymin))));
            double par[9]={myTheta12,                       //Theta12
                myTheta23,                  //Theta23
                myTheta13,                       //Theta13
                myDeltam32*TMath::Power(10,-3.), //DeltaM32^2
                myDeltam21*TMath::Power(10,-5.),  //DeltaM21^2
                myCP,                       //DeltaCP (3Pi/2)
                myRockDensity,                        //Rock Density
                myL,                        //L
                -1};                        //NuAntiNu
            
            SetOscParam(par);
            hcont_antinu_norm->SetBinContent(ix,iy, model3flav_mat_ElecToMu(x));
        }
    }
    
    //antineutrino invert hierarchy
    TH2F *hcont_antinu_invt=0;
    hcont_antinu_invt = new TH2F("hcont_antinu_invt","hcont_antinu_invt",npointsx, binx, npointsy, biny);
    hcont_antinu_invt->GetXaxis()->SetTitle("True Neutrino Energy (GeV)");
    hcont_antinu_invt->GetYaxis()->SetTitle("Baseline (km)");
    hcont_antinu_invt->GetXaxis()->CenterTitle();
    hcont_antinu_invt->GetYaxis()->CenterTitle();
    
    for (Int_t ix=0; ix < npointsx+1; ix++){
        for (Int_t iy=0; iy < npointsy+1; iy++){
            x[0]= TMath::Power(10.,(TMath::Log10(xmin)+((ix+0.5)*1.0/npointsx)*(TMath::Log10(xmax)-TMath::Log10(xmin))));
            myL = TMath::Power(10.,(TMath::Log10(ymin)+((iy+0.5)*1.0/npointsy)*(TMath::Log10(ymax)-TMath::Log10(ymin))));
            double par[9]={myTheta12,                       //Theta12
                myTheta23,                  //Theta23
                myTheta13,                       //Theta13
                -1.0*myDeltam32*TMath::Power(10,-3.), //DeltaM32^2
                myDeltam21*TMath::Power(10,-5.),  //DeltaM21^2
                myCP,                       //DeltaCP (3Pi/2)
                myRockDensity,                        //Rock Density
                myL,                        //L
                -1};                        //NuAntiNu
            
            SetOscParam(par);
            hcont_antinu_invt->SetBinContent(ix,iy, model3flav_mat_ElecToMu(x));
        }
    }
    
   
    TCanvas *canvas1=0;
    canvas1 = new TCanvas("Can1","Can1",800,600);
    gStyle->SetOptStat(0);
    canvas1->cd();
    Int_t colors[]={0,41,42,45,95};
    gStyle->SetPalette((sizeof(colors)/sizeof(Int_t)), colors);
    Double_t levels[]={0.001,0.01,0.05,0.1};
    hcont_nu_norm->SetContour((sizeof(levels)/sizeof(Double_t)), levels);
    canvas1->SetLogx(1);
    canvas1->SetLogy(1);
    canvas1->SetLogz(1);
    hcont_nu_norm->SetTitle("");
    hcont_nu_norm->Draw("colz");
    canvas1->SaveAs(fsubname+fcont+fprob+"_"+"nu_norm.gif");
    canvas1->SaveAs(fsubname+fcont+fprob+"_"+"nu_norm.eps");
    
    TCanvas *canvas2=0;
    canvas2 = new TCanvas("Can2","Can2",800,600);
    gStyle->SetOptStat(0);
    canvas2->cd();
    Int_t colors[]={0,41,42,45,95};
    gStyle->SetPalette((sizeof(colors)/sizeof(Int_t)), colors);
    Double_t levels[]={0.001,0.01,0.05,0.1};
    hcont_nu_invt->SetContour((sizeof(levels)/sizeof(Double_t)), levels);
    canvas2->SetLogx(1);
    canvas2->SetLogy(1);
    canvas2->SetLogz(1);
    hcont_nu_invt->SetTitle("");
    hcont_nu_invt->Draw("colz");
    canvas2->SaveAs(fsubname+fcont+fprob+"_"+"nu_invt.gif");
    canvas2->SaveAs(fsubname+fcont+fprob+"_"+"nu_invt.eps");
   
    TCanvas *canvas3=0;
    canvas3 = new TCanvas("Can3","Can3",800,600);
    gStyle->SetOptStat(0);
    canvas3->cd();
    Int_t colors[]={0,41,42,45,95};
    gStyle->SetPalette((sizeof(colors)/sizeof(Int_t)), colors);
    Double_t levels[]={0.001,0.01,0.05,0.1};
    hcont_antinu_norm->SetContour((sizeof(levels)/sizeof(Double_t)), levels);
    canvas3->SetLogx(1);
    canvas3->SetLogy(1);
    canvas3->SetLogz(1);
    hcont_antinu_norm->SetTitle("");
    hcont_antinu_norm->Draw("colz");
    canvas3->SaveAs(fsubname+fcont+fprob+"_"+"antinu_norm.gif");
    canvas3->SaveAs(fsubname+fcont+fprob+"_"+"antinu_norm.eps");
    
    TCanvas *canvas4=0;
    canvas4 = new TCanvas("Can4","Can4",800,600);
    gStyle->SetOptStat(0);
    canvas4->cd();
    Int_t colors[]={0,41,42,45,95};
    gStyle->SetPalette((sizeof(colors)/sizeof(Int_t)), colors);
    Double_t levels[]={0.001,0.01,0.05,0.1};
    hcont_antinu_invt->SetContour((sizeof(levels)/sizeof(Double_t)), levels);
    canvas4->SetLogx(1);
    canvas4->SetLogy(1);
    canvas4->SetLogz(1);
    hcont_antinu_invt->SetTitle("");
    hcont_antinu_invt->Draw("colz");
    canvas4->SaveAs(fsubname+fcont+fprob+"_"+"antinu_invt.gif");
    canvas4->SaveAs(fsubname+fcont+fprob+"_"+"antinu_invt.eps");
    
    TString rootOp =fsubname+fcont+fprob+".root";
    TFile *ffile = new TFile(rootOp,"RECREATE");
    hcont_nu_norm->Write("hcont_nu_norm");
    hcont_nu_invt->Write("hcont_nu_invt");
    hcont_antinu_norm->Write("hcont_antinu_norm");
    hcont_antinu_invt->Write("hcont_antinu_invt");
    
    ffile->Close();
    
    exit(0);
    
}
