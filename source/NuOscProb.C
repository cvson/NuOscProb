#include "OscWeight.h"
#include "TMath.h"
#include <cmath>
#include <cstdlib>
#include <iostream>

using namespace std;

void GetOscParam(double* par);
void SetOscParam(double * par, bool force = false);
void SetOscParam(OscPar::OscPar_t pos, double val, bool force = false);
//====================================================================
//====================================================================
//In vacuum
//Two-flavor vacuum function
double model2flav_vac_MuToElec(double *x);
double model2flav_vac_MuToTau(double *x);
double model2flav_vac_MuSurvive(double *x);
double model2flav_vac_ElecToElec(double *x);
double model2flav_vac_ElecToMu(double *x);
double model2flav_vac_ElecToTau(double *x);

//Three-flavor vacuum function
//No CP
double model3flav_vac_apx_nocp_MuToElec(double *x);
double model3flav_vac_apx_nocp_MuToTau(double *x);
double model3flav_vac_apx_nocp_MuSurvive(double *x);
double model3flav_vac_apx_nocp_ElecToElec(double *x);
double model3flav_vac_apx_nocp_ElecToMu(double *x);
double model3flav_vac_apx_nocp_ElecToTau(double *x);


//Three-flavor vacuum function
//First order in sin13
double model3flav_vac_1stsin13_MuToElec(double *x);
double model3flav_vac_1stsin13_MuToTau(double *x);
double model3flav_vac_1stsin13_MuSurvive(double *x);
double model3flav_vac_1stsin13_ElecToElec(double *x);
double model3flav_vac_1stsin13_ElecToMu(double *x);
double model3flav_vac_1stsin13_ElecToTau(double *x);

//First order in alpha
double model3flav_vac_1stalpha_MuToElec(double *x);
double model3flav_vac_1stalpha_MuToTau(double *x);
double model3flav_vac_1stalpha_MuSurvive(double *x);
double model3flav_vac_1stalpha_ElecToElec(double *x);
double model3flav_vac_1stalpha_ElecToMu(double *x);
double model3flav_vac_1stalpha_ElecToTau(double *x);

//second order in alpha and s13
double model3flav_vac_2nd_MuToElec(double *x);
double model3flav_vac_2nd_MuToTau(double *x);
double model3flav_vac_2nd_MuSurvive(double *x);
double model3flav_vac_2nd_ElecToElec(double *x);
double model3flav_vac_2nd_ElecToMu(double *x);
double model3flav_vac_2nd_ElecToTau(double *x);

//analytic prob from Joachim Kopp thesis
double model3flav_vac_kopp_MuToElec(double *x);
double model3flav_vac_kopp_ElecToElec(double *x);
double model2flav_mat_kopp_sol_MuToElec(double *x);
double model2flav_mat_kopp_atm_MuToElec(double *x);
double model3flav_mat_kopp_2nd_MuToElec(double *x);
double model3flav_vac_kopp_2nd_MuToElec(double *x);
//====================================================================
//Include matter effect
//====================================================================
//In matter code by Josua Havard
//First order in sin13
double model3flav_mat_1stsin13_MuToElec(double *x);
double model3flav_mat_1stsin13_MuToTau(double *x);
double model3flav_mat_1stsin13_MuSurvive(double *x);
double model3flav_mat_1stsin13_ElecToElec(double *x);
double model3flav_mat_1stsin13_ElecToMu(double *x);
double model3flav_mat_1stsin13_ElecToTau(double *x);

//First order in alpha
double model3flav_mat_1stalpha_MuToElec(double *x);
double model3flav_mat_1stalpha_MuToTau(double *x);
double model3flav_mat_1stalpha_MuSurvive(double *x);
double model3flav_mat_1stalpha_ElecToElec(double *x);
double model3flav_mat_1stalpha_ElecToMu(double *x);
double model3flav_mat_1stalpha_ElecToTau(double *x);

//Three-flavor function up to the second order of alpha and theta13
double model3flav_mat_MuToElec(double *x);
double model3flav_mat_MuToTau(double *x);
double model3flav_mat_MuSurvive(double *x);
double model3flav_mat_ElecToElec(double *x);
double model3flav_mat_ElecToMu(double *x);
double model3flav_mat_ElecToTau(double *x);

//Perdue formula
double perdue3flav_mat_MuToElec(double *x);


double fOscPar[10];

double fSinTh[3];
double fCosTh[3];
double fSin2Th[3];
double fCos2Th[3];
double fElecDensity;
double fV;
double fSinAL;
double fSinDCP;
double fCosDCP;

int main(void)
{}

void SetOscParam(OscPar::OscPar_t pos, double val, bool force)
{
  //cout<<"setting  "<<pos<<"  "<<val<<endl;
  if(fabs(fOscPar[pos] - val) > 1e-9 || force){
    fOscPar[pos] = val;
    if(pos < 3){
       fSinTh[pos] = sin(val);
       fCosTh[pos] = cos(val);
       fSin2Th[pos] = sin(2*val);
       fCos2Th[pos] = cos(2*val);
    }
    if(pos == OscPar::kDensity){
      double ne = OscPar::z_a*OscPar::A_av*val; //electron density #/cm^{3}
      fElecDensity = ne*OscPar::invCmToeV*OscPar::invCmToGeV*OscPar::invCmToGeV;
      //electron density with units Gev^{2} eV
      //Gev^{2} to cancel with GeV^{-2} in Gf
                                                                                
       fV = OscPar::root2*OscPar::Gf*fElecDensity; //eV
       fSinAL = sin(fV*fOscPar[OscPar::kL]/(OscPar::invKmToeV*2.));       
    }
    if(pos == OscPar::kDelta){
       fCosDCP = cos(val); 
       fSinDCP = sin(val);
    }
  }
}

void SetOscParam(double *par, bool force)
{
   for(int i = 0; i < int(OscPar::kUnknown); i++){
     SetOscParam(OscPar::OscPar_t(i), par[i], force);
     cout << "Setting osc param: " << i << " to: " << par[i] << endl;
   }
}

void GetOscParam(double *par)
{
   for(int i = 0; i < int(OscPar::kUnknown); i++){
     par[i] = fOscPar[i];
     cout << "Param: " << i << " is set to: " << par[i] << endl;
   }
}


//====================================================================
//====================================================================
// matter three-flavor
double model3flav_mat_MuToElec(double *x)
{
  double E = x[0]; //energy

  //Building the standard terms
  double L = fOscPar[OscPar::kL];  //baseline
  double dmsq_23 = fOscPar[OscPar::kDeltaM23];
  double dmsq_12 = fOscPar[OscPar::kDeltaM12];
  double dmsq_13 = dmsq_23+dmsq_12; //eV^{2}

  double sinsq_2th12 = fSin2Th[OscPar::kTh12]*fSin2Th[OscPar::kTh12];
  double sinsq_2th13 = fSin2Th[OscPar::kTh13]*fSin2Th[OscPar::kTh13];
                                                                                
  double cos_th23 = fCosTh[OscPar::kTh23];
  double cos_th12 = fCosTh[OscPar::kTh12];
  double sin_th13 = fSinTh[OscPar::kTh13];
                                                                                
  double sin_2th23 = fSin2Th[OscPar::kTh23];
  double sin_2th12 = fSin2Th[OscPar::kTh12];
  double cos_2th13 = fCos2Th[OscPar::kTh13];
  double cos_2th12 = fCos2Th[OscPar::kTh12];

  double sinsq_th23 = fSinTh[OscPar::kTh23]*fSinTh[OscPar::kTh23];
  double sinsq_th12 = fSinTh[OscPar::kTh12]*fSinTh[OscPar::kTh12];

  double d_cp = fOscPar[OscPar::kDelta];
  double cos_dcp = fCosDCP;
  double sin_dcp = fSinDCP;

  //Building the more complicated terms
  double Delta = dmsq_13*L/(4*E*1e9*OscPar::invKmToeV);
  double A = 2*fV*E*1e9/(dmsq_13);
  double alpha = dmsq_12/dmsq_13;

  // A and d_cp both change sign for antineutrinos
  double plusminus = int(fOscPar[OscPar::kNuAntiNu]);
  A *= plusminus;
  d_cp *= plusminus;
  sin_dcp *= plusminus;

  //Now calculate the resonance terms for alpha expansion (C13) and s13 expansion (C12)

  double C13 = TMath::Sqrt(sinsq_2th13+(A-cos_2th13)*(A-cos_2th13));

  double C12 = 1;  //really C12 -> infinity when alpha = 0 but not an option really
  if(fabs(alpha) > 1e-10){  //want to be careful here
    double temp = cos_2th12 - A/alpha;
    C12 = TMath::Sqrt(sinsq_2th12+(temp*temp));
  }

  //More complicated sin and cosine terms
  double cosC13Delta = cos(C13*Delta);
  double sinC13Delta = sin(C13*Delta);
  double sin1pADelta = sin((A+1)*Delta);
  double cos1pADelta = cos((A+1)*Delta);
  double sinADelta = sin(A*Delta);
  double sinAm1Delta = sin((A-1)*Delta);
  double cosdpDelta = cos(d_cp+Delta);
  double sinApam2Delta = sin((A+alpha-2)*Delta);
  double cosApam2Delta = cos((A+alpha-2)*Delta);

  double cosaC12Delta = 0;
  double sinaC12Delta = 0; 
  
  if(fabs(alpha) > 1e-10){
    cosaC12Delta = cos(alpha*C12*Delta);
    sinaC12Delta = sin(alpha*C12*Delta);
  } // otherwise not relevant

  //First we calculate the terms for the alpha expansion (good to all orders in th13)
  // this is the equivalent of Eq 47 & 48 corrected for Mu to E instead of E to Mu

  // Leading order term 
  double p1 = sinsq_th23*sinsq_2th13*sinC13Delta*sinC13Delta/(C13*C13);

  // terms that appear at order alpha
  double p2Inner = Delta*cosC13Delta*(1-A*cos_2th13)/C13 - 
                      A*sinC13Delta*(cos_2th13-A)/(C13*C13); 

  double p2 = -2*sinsq_th12*sinsq_th23*sinsq_2th13*sinC13Delta/(C13*C13)*p2Inner*alpha;

  double p3Inner = - sin_dcp*(cosC13Delta - cos1pADelta)*C13 
         + cos_dcp*(C13*sin1pADelta - (1-A*cos_2th13)*sinC13Delta);

  double p3 = sin_2th12*sin_2th23*sin_th13*sinC13Delta/(A*C13*C13)*p3Inner*alpha;

  //  p1 + p2 + p3 is the complete contribution for this expansion
  
  // Now for the expansion in orders of sin(th13) (good to all order alpha) 
  //  this is the equivalent of Eq 65 and 66

  // leading order term
  double pa1 = 0.0, pa2 = 0.0;

  if(fabs(alpha) > 1e-10){
    // leading order term
    pa1 = cos_th23*cos_th23*sinsq_2th12*sinaC12Delta*sinaC12Delta/(C12*C12);

    // and now to calculate the first order in s13 term
    double t1 = (cos_2th12 - A/alpha)/C12 
                  - alpha*A*C12*sinsq_2th12/(2*(1-alpha)*C12*C12);
    double t2 = -cos_dcp*(sinApam2Delta-sinaC12Delta*t1);
  
    double t3 = -(cosaC12Delta-cosApam2Delta)*sin_dcp;
 
    double denom = (1-A-alpha+A*alpha*cos_th12*cos_th12)*C12;
    double t4 = sin_2th12*sin_2th23*(1-alpha)*sinaC12Delta/denom;

    pa2 = t4*(t3+t2)*sin_th13;
  }
  //pa1+pa2 is the complete contribution from this expansion

  // In order to combine the information correctly we need to add the two
  //  expansions and subtract off the terms that are in both (alpha^1, s13^1) 
  //  these may be taken from the expansion to second order in both parameters
  //  Equation 31 

  double t1 = sinADelta*cosdpDelta*sinAm1Delta/(A*(A-1));
  double repeated = 2*alpha*sin_2th12*sin_2th23*sin_th13*t1;

  //  Calculate the total probability
  double totalP = p1+p2+p3 + (pa1+pa2) - repeated;
    if (totalP>1.0) return 1.0;
    else if (totalP<1e-4) return 1e-4;
    else return totalP;
}
//===================================
//MuToTau
double model3flav_mat_MuToTau(double *x)
{
  double E = x[0]; //energy

  double L = fOscPar[OscPar::kL];  //baseline
  double dmsq_23 = fOscPar[OscPar::kDeltaM23];
  double dmsq_12 = fOscPar[OscPar::kDeltaM12];
  double dmsq_13 = dmsq_23+dmsq_12; //eV^{2}

  double sinsq_2th12 = fSin2Th[OscPar::kTh12]*fSin2Th[OscPar::kTh12];
  double sinsq_2th13 = fSin2Th[OscPar::kTh13]*fSin2Th[OscPar::kTh13];
  double sinsq_2th23 = fSin2Th[OscPar::kTh23]*fSin2Th[OscPar::kTh23];
                                                                                
  double cos_th13 = fCosTh[OscPar::kTh13];
  double cos_th12 = fCosTh[OscPar::kTh12];
  double sin_th12 = fSinTh[OscPar::kTh12];
  double sin_th13 = fSinTh[OscPar::kTh13];
                                                                                
  double sin_2th23 = fSin2Th[OscPar::kTh23];
//  double sin_2th13 = fSin2Th[OscPar::kTh13];
  double sin_2th12 = fSin2Th[OscPar::kTh12];
                                                                                
  double cos_2th23 = fCos2Th[OscPar::kTh23];
  double cos_2th13 = fCos2Th[OscPar::kTh13];
  double cos_2th12 = fCos2Th[OscPar::kTh12];

  double d_cp = fOscPar[OscPar::kDelta];
  double cos_dcp = fCosDCP;
  double sin_dcp = fSinDCP;

  //Building the more complicated terms                                                                              
  double Delta = dmsq_13*L/(4*E*1e9*OscPar::invKmToeV);
  double A = 2*fV*E*1e9/(dmsq_13);
  double alpha = dmsq_12/dmsq_13;
  
  // A and d_cp both change sign for antineutrinos
  double plusminus = int(fOscPar[OscPar::kNuAntiNu]);
  A *= plusminus;
  d_cp *= plusminus;
  sin_dcp *= plusminus;

  //Now calculate the resonance terms for alpha expansion (C13) and s13 expansion (C12)
                                                                                                                       
  double C13 = TMath::Sqrt(sinsq_2th13+(A-cos_2th13)*(A-cos_2th13));
                                                                                                                       
  double C12 = 1;  //really C12 -> infinity when alpha = 0 but not an option really
  if(fabs(alpha) > 1e-10){  //want to be careful here
    double temp = cos_2th12 - A/alpha;
    C12 = TMath::Sqrt(sinsq_2th12+(temp*temp));
  }
                                                                                                                       
  //More complicated sin and cosine terms
  double cosC13Delta = cos(C13*Delta);
  double sinC13Delta = sin(C13*Delta);
  double sin1pADelta = sin((A+1)*Delta);
  double cos1pADelta = cos((A+1)*Delta);
  double sinADelta = sin((A)*Delta);

  double sin1pAmCDelta = sin(0.5*(A+1-C13)*Delta);
  double sin1pApCDelta = sin(0.5*(A+1+C13)*Delta);

  double cosaC12Delta = 0;
  double sinaC12Delta = 0;
                                                                                                                       
  if(fabs(alpha) > 1e-10){
    cosaC12Delta = cos(alpha*C12*Delta);
    sinaC12Delta = sin(alpha*C12*Delta);
  } // otherwise not relevant

  double sinApam2Delta = sin((A+alpha-2)*Delta);
  double cosApam2Delta = cos((A+alpha-2)*Delta);
  double sinAm1Delta = sin((A-1)*Delta);
  double cosAm1Delta = cos((A-1)*Delta);
  double sinDelta = sin(Delta);
  double sin2Delta = sin(2*Delta);

  double cosaC12pApam2Delta = 0;                                                                                                                       
  if(fabs(alpha) > 1e-10){
    cosaC12pApam2Delta = cos((alpha*C12+A+alpha-2)*Delta);
  }

  //First we calculate the terms for the alpha expansion (good to all orders in th13)
  // this is the equivalent of Eq 49 & 50 corrected for Mu to E instead of E to Mu

  // Leading order term
  double pmt_0 = 0.5*sinsq_2th23;
  pmt_0 *= (1 - (cos_2th13-A)/C13)*sin1pAmCDelta*sin1pAmCDelta 
             +  (1 + (cos_2th13-A)/C13)*sin1pApCDelta*sin1pApCDelta
             - 0.5*sinsq_2th13*sinC13Delta*sinC13Delta/(C13*C13);

  // terms that appear at order alpha
  double t0, t1, t2, t3;
  t0 = (cos_th12*cos_th12-sin_th12*sin_th12*sin_th13*sin_th13
          *(1+2*sin_th13*sin_th13*A+A*A)/(C13*C13))*cosC13Delta*sin1pADelta*2;
  t1 = 2*(cos_th12*cos_th12*cos_th13*cos_th13-cos_th12*cos_th12*sin_th13*sin_th13
         +sin_th12*sin_th12*sin_th13*sin_th13
         +(sin_th12*sin_th12*sin_th13*sin_th13-cos_th12*cos_th12)*A);
  t1 *= sinC13Delta*cos1pADelta/C13;

  t2 =  sin_th12*sin_th12*sinsq_2th13*sinC13Delta/(C13*C13*C13);
  t2 *= A/Delta*sin1pADelta+A/Delta*(cos_2th13-A)/C13*sinC13Delta
          - (1-A*cos_2th13)*cosC13Delta;

  double pmt_1 = -0.5*sinsq_2th23*Delta*(t0+t1+t2);   

  t0 = t1 = t2 = t3 = 0.0;

  t0 = cosC13Delta-cos1pADelta;
  t1 = 2*cos_th13*cos_th13*sin(d_cp)*sinC13Delta/C13*t0;
  t2 = -cos_2th23*cos_dcp*(1+A)*t0*t0;

  t3  = cos_2th23*cos_dcp*(sin1pADelta+(cos_2th13-A)/C13*sinC13Delta);
  t3 *= (1+2*sin_th13*sin_th13*A + A*A)*sinC13Delta/C13 - (1+A)*sin1pADelta;

//  cout<<t1<<"  "<<t2<<"  "<<t3<<endl;

  pmt_1 = pmt_1 + (t1+t2+t3)*sin_th13*sin_2th12*sin_2th23/(2*A*cos_th13*cos_th13);
  pmt_1 *= alpha;

  //  pmt_0 + pmt_1 is the complete contribution for this expansion
                                                                                                                       
  // Now for the expansion in orders of sin(th13) (good to all order alpha)
  //  this is the equivalent of Eq 67 and 68
                                                                                                                       
  // leading order term
  double pmt_a0 =  0.5*sinsq_2th23;

  pmt_a0 *= 1 - 0.5*sinsq_2th12*sinaC12Delta*sinaC12Delta/(C12*C12)
              - cosaC12pApam2Delta
              - (1 - (cos_2th12 - A/alpha)/C12)*sinaC12Delta*sinApam2Delta;
            
  double denom = (1-A-alpha+A*alpha*cos_th12*cos_th12)*C12;

  t0 = (cosaC12Delta-cosApam2Delta)*(cosaC12Delta-cosApam2Delta);
           
  t1 = (cos_2th12 - A/alpha)/C12*sinaC12Delta+sinApam2Delta;
            
  t2 = ((cos_2th12 - A/alpha)/C12+2*(1-alpha)/(alpha*A*C12))*sinaC12Delta
             + sinApam2Delta;

  t3 = (alpha*A*C12)/2.0*cos_2th23*cos_dcp*(t0 + t1*t2);
  
  t3 += sin_dcp*(1-alpha)*(cosaC12Delta-cosApam2Delta)*sinaC12Delta;

  double pmt_a1 = sin_th13*sin_2th12*sin_2th23/denom*t3;

  // pmt_a1+pmt_a2 is the complete contribution from this expansion
                                                                                                                       
  // In order to combine the information correctly we need to add the two
  //  expansions and subtract off the terms that are in both (alpha^1, s13^1)
  //  and lower order terms
  //  these may be taken from the expansion to second order in both parameters
  //  Equation 34


  // Now for the term of order alpha * s13 or lower order!
  t0 = t1 = t2 = t3 = 0.0;

  t1 = +sin_dcp*sinDelta*sinADelta*sinAm1Delta/(A*(A-1));
  t2 = -1/(A-1)*cos_dcp*sinDelta*(A*sinDelta-sinADelta*cosAm1Delta/A)*cos_2th23/denom;
  t0 =  2*alpha*sin_2th12*sin_2th23*sin_th13*(t1+t2);

  t1 = sinsq_2th23*sinDelta*sinDelta 
       - alpha*sinsq_2th23*cos_th12*cos_th12*Delta*sin2Delta;

  double repeated = t0+t1;

  //  Calculate the total probability
  double totalP = pmt_0 + pmt_1 + pmt_a0 + pmt_a1 - repeated;

    if (totalP>1.0) return 1.0;
    else if (totalP<1e-4) return 1e-4;
    else return totalP;
}

//===================================
//MuToMu
double model3flav_mat_MuSurvive(double *x)
{
  double p1 = 1. - model3flav_mat_MuToTau(x) - model3flav_mat_MuToElec(x);
  if(p1 < 0){ cout<<"Damnation! "<<x[0]<<" "<<p1<<endl; p1 = 0;}
    if (p1>1.0) return 1.0;
    else if (p1<1e-4) return 1e-4;
    else return p1;
}

double model3flav_mat_ElecToTau(double *x)
{
//  EtoTau is the same as E->Mu wich sinsq_th23 <-> cossq_th23, sin(2th23) <->-sin(2th23)
  double origCos = fCosTh[OscPar::kTh23];
  double origSin = fSinTh[OscPar::kTh23];
  double orig2Sin = fSin2Th[OscPar::kTh23];

  fCosTh[OscPar::kTh23] = origSin;
  fSinTh[OscPar::kTh23] = origCos;
  fSin2Th[OscPar::kTh23] = -orig2Sin;

  double prob = model3flav_mat_ElecToMu(x);

  //restore the world
  fCosTh[OscPar::kTh23] = origCos;
  fSinTh[OscPar::kTh23] = origSin;
  fSin2Th[OscPar::kTh23] = orig2Sin;

  return prob;

/*  This is an option as well but I think results in more cycles
  double p1 = 1. - OscCalc::model3flav_mat_ElecToMu(x) - OscCalc::model3flav_mat_ElecToElec(x);
  if(p1 < 0){ cout<<"Damnation! "<<x[0]<<" "<<p1<<endl; p1 = 0;}
  return p1;
*/ 
}

//===================================
//ElectoMu
double model3flav_mat_ElecToMu(double *x)
{
  // Flip delta to reverse direction
  double oldSinDelta = fSinDCP;
  double oldDelta =  fOscPar[OscPar::kDelta];

  fOscPar[OscPar::kDelta] = -oldDelta;
  fSinDCP = -oldSinDelta;
                                                                                                                       
  double prob = model3flav_mat_MuToElec(x);
                                                                                                                       
  //restore the world
  fOscPar[OscPar::kDelta] = oldDelta;
  fSinDCP = oldSinDelta;
 
  return prob;   
}
//===================================
//ElectoElec
double model3flav_mat_ElecToElec(double *x)
{
  double E = x[0]; //energy
                                                                                                                       
  //Building the standard terms
  double L = fOscPar[OscPar::kL];  //baseline
  double dmsq_23 = fOscPar[OscPar::kDeltaM23];
  double dmsq_12 = fOscPar[OscPar::kDeltaM12];
  double dmsq_13 = dmsq_23+dmsq_12; //eV^{2}
                                                                                                                       
  double sinsq_2th12 = fSin2Th[OscPar::kTh12]*fSin2Th[OscPar::kTh12];
  double sinsq_2th13 = fSin2Th[OscPar::kTh13]*fSin2Th[OscPar::kTh13];
                                                                                                                       
  double sin_th12 = fSinTh[OscPar::kTh12];
 
//  double cos_2th23 = fCos2Th[OscPar::kTh23];
  double cos_2th13 = fCos2Th[OscPar::kTh13];
  double cos_2th12 = fCos2Th[OscPar::kTh12];
                                                                                                                       
  double d_cp = fOscPar[OscPar::kDelta];
  double sin_dcp = fSinDCP;
                                                                                                                       
  //Building the more complicated terms
  double Delta = dmsq_13*L/(4*E*1e9*OscPar::invKmToeV);
  double A = 2*fV*E*1e9/(dmsq_13);
  double alpha = dmsq_12/dmsq_13;
                                                                                                                       
  // A and d_cp both change sign for antineutrinos
  double plusminus = int(fOscPar[OscPar::kNuAntiNu]);
  A *= plusminus;
  d_cp *= plusminus;
  sin_dcp *= plusminus;
                                                                                                                       
  //Now calculate the resonance terms for alpha expansion (C13) and s13 expansion (C12)
  double C13 = TMath::Sqrt(sinsq_2th13+(A-cos_2th13)*(A-cos_2th13));
                                                                                                                       
  double C12 = 1;  //really C12 -> infinity when alpha = 0 but not an option really
  if(fabs(alpha) > 1e-10){  //want to be careful here
    double temp = cos_2th12 - A/alpha;
    C12 = TMath::Sqrt(sinsq_2th12+(temp*temp));
  }
                                                                                                                       
  //More complicated sin and cosine terms
  double cosC13Delta = cos(C13*Delta);
  double sinC13Delta = sin(C13*Delta);
                                                                                                                       
  double cosaC12Delta = 0;
  double sinaC12Delta = 0;
                                                                                                                       
  if(fabs(alpha) > 1e-10){
    cosaC12Delta = cos(alpha*C12*Delta);
    sinaC12Delta = sin(alpha*C12*Delta);
  } // otherwise not relevant
                                                                                                                       
  //First we calculate the terms for the alpha expansion (good to all orders in th13)
  // this is the equivalent of Eq 45 & 46 corrected for Mu to E instead of E to Mu
                                                                                                                       
  // Leading order term
  double p1 = 1 - sinsq_2th13*sinC13Delta*sinC13Delta/(C13*C13);
                                                                                                                       
  // terms that appear at order alpha
  double p2Inner = Delta*cosC13Delta*(1-A*cos_2th13)/C13 -
                      A*sinC13Delta*(cos_2th13-A)/(C13*C13);
                                                                                                                       
  double p2 = +2*sin_th12*sin_th12*sinsq_2th13*sinC13Delta/(C13*C13)*p2Inner*alpha;
  //  p1 + p2 is the complete contribution for this expansion
                                                                                                                       
  // Now for the expansion in orders of sin(th13) (good to all order alpha)
  //  this is the equivalent of Eq 63 and 64
                                                                                                                       
  // leading order term
  double pa1 = 1.0, pa2 = 0.0;
                                                                                                                       
  if(fabs(alpha) > 1e-10){
    // leading order term
    pa1 = 1.0 - sinsq_2th12*sinaC12Delta*sinaC12Delta/(C12*C12);
  }
  //pa1 is the complete contribution from this expansion, there is no order s13^1 term
                                                                                                                       
  // In order to combine the information correctly we need to add the two
  //  expansions and subtract off the terms that are in both (alpha^1, s13^1)
  //  these may be taken from the expansion to second order in both parameters
  //  Equation 30
                                                                                                                       
  double repeated = 1;
                                                                                                                       
  //  Calculate the total probability
  double totalP = p1+p2 + (pa1+pa2) - repeated;
                                                                                                                       
  return totalP;
}

//====================================================================================
//====================================================================
//Two flavor formular
//===================================
//MuToElec
double model2flav_vac_MuToElec(double *x)
{
    return 0.0;
}
//===================================
//MuToTau
double model2flav_vac_MuToTau(double *x)
{
    double E = x[0]; //energy
    
    double L = fOscPar[OscPar::kL];  //baseline
    double dmsq_23 = fOscPar[OscPar::kDeltaM23];
    
    double sinsq_2th23 = fSin2Th[OscPar::kTh23]*fSin2Th[OscPar::kTh23];
    
    double sin_BL = sin(1.267*L*dmsq_23/E);
    double p1 = sinsq_2th23*sinsq_2th23*sin_BL*sin_BL;
    return p1;
}
//===================================
//MuToMu

double model2flav_vac_MuSurvive(double *x)
{
    double p1 = model2flav_vac_MuToTau(x);
    return 1.-p1;
}
//===================================
//ElecToTau
double model2flav_vac_ElecToTau(double *x)
{
    
    
    return 0.0;
    
    
}
//===================================
//ElecToMu
double model2flav_vac_ElecToMu(double *x)
{
    double E = x[0]; //energy
    
    double L = fOscPar[OscPar::kL];  //baseline
    double dmsq_12 = fOscPar[OscPar::kDeltaM12];
    
    double sinsq_2th12 = fSin2Th[OscPar::kTh12]*fSin2Th[OscPar::kTh12];
    
    double sin_BL = sin(1.267*L*dmsq_12/E);
    double p1 = sinsq_2th12*sinsq_2th12*sin_BL*sin_BL;
    return p1;
}

//===================================
//ElecToElec
double model2flav_vac_ElecToElec(double *x)
{

    double p1 = model2flav_vac_ElecToMu(x);
    return 1.-p1;
}
//====================================================================================
//Three-flavor vacuum function
//http://iopscience.iop.org/1126-6708/2004/04/078
//First order in sin13
double model3flav_vac_1stsin13_MuToElec(double *x){
    double E = x[0]; //energy
    double L = fOscPar[OscPar::kL];  //baseline
    double dmsq_23 = fOscPar[OscPar::kDeltaM23];
    double dmsq_12 = fOscPar[OscPar::kDeltaM12];
    double dmsq_13 = dmsq_23+dmsq_12; //eV^{2}
    
    double sinsq_2th12 = fSin2Th[OscPar::kTh12]*fSin2Th[OscPar::kTh12];
    double sinsq_2th13 = fSin2Th[OscPar::kTh13]*fSin2Th[OscPar::kTh13];
    double sinsq_2th23 = fSin2Th[OscPar::kTh23]*fSin2Th[OscPar::kTh23];
    
    double cos_th13 = fCosTh[OscPar::kTh13];
    double cos_th23 = fCosTh[OscPar::kTh23];
    double cos_th12 = fCosTh[OscPar::kTh12];
    double sin_th13 = fSinTh[OscPar::kTh13];
    
    double sin_2th23 = fSin2Th[OscPar::kTh23];
    double sin_2th12 = fSin2Th[OscPar::kTh12];
    double sin_2th13 = fSin2Th[OscPar::kTh13];
    
    double cos_2th13 = fCos2Th[OscPar::kTh13];
    double cos_2th12 = fCos2Th[OscPar::kTh12];
    double cos_2th23 = fCos2Th[OscPar::kTh23];
    
    double sinsq_th23 = fSinTh[OscPar::kTh23]*fSinTh[OscPar::kTh23];
    double sinsq_th12 = fSinTh[OscPar::kTh12]*fSinTh[OscPar::kTh12];
    double sinsq_th13 = fSinTh[OscPar::kTh13]*fSinTh[OscPar::kTh13];
    
    double d_cp = fOscPar[OscPar::kDelta];
    double cos_dcp = fCosDCP;
    double sin_dcp = fSinDCP;
    
    //Building the more complicated terms
    double Delta = dmsq_13*L/(4*E*1e9*OscPar::invKmToeV);
    double A = 2*fV*E*1e9/(dmsq_13);
    double alpha = dmsq_12/dmsq_13;
    
    // A and d_cp both change sign for antineutrinos
    double plusminus = int(fOscPar[OscPar::kNuAntiNu]);
    A *= plusminus;
    d_cp *= plusminus;
    sin_dcp *= plusminus;
    
    //Now calculate the resonance terms for alpha expansion (C13) and s13 expansion (C12)
    
    //first term
    double tmue0 = sinsq_2th12*pow(cos_th23,2)*pow(sin(alpha*Delta),2);
    //double tmue0 = sinsq_2th12*pow(cos_th23,2)*pow(sin(alpha),2);
    //double tmue0 = sinsq_2th12*pow(cos_th23,2)*pow(sin(alpha),2)*Delta;
    //second term
    double tmue11 = cos_dcp*sin_2th12*sin_2th23*(pow(sin(Delta),2)-pow(sin((1.-alpha)*Delta),2)+cos_2th12*pow(sin(alpha*Delta),2));
    //change the sign here for CP, comparing with the paper, since we calculate for mu2e not e2mu
    double tmue12 = -0.5*sin_dcp*sin_2th12*sin_2th23*(-sin(2*Delta)+sin(2*(1-alpha)*Delta)+sin(2*alpha*Delta));
    
    double mu2e = tmue0+sin_th13*(tmue11+tmue12);
    return mu2e;
}
//-----------------------------------------------
double model3flav_vac_1stsin13_MuToTau(double *x){
    double E = x[0]; //energy
    double L = fOscPar[OscPar::kL];  //baseline
    double dmsq_23 = fOscPar[OscPar::kDeltaM23];
    double dmsq_12 = fOscPar[OscPar::kDeltaM12];
    double dmsq_13 = dmsq_23+dmsq_12; //eV^{2}
    
    double sinsq_2th12 = fSin2Th[OscPar::kTh12]*fSin2Th[OscPar::kTh12];
    double sinsq_2th13 = fSin2Th[OscPar::kTh13]*fSin2Th[OscPar::kTh13];
    double sinsq_2th23 = fSin2Th[OscPar::kTh23]*fSin2Th[OscPar::kTh23];
    
    double cos_th13 = fCosTh[OscPar::kTh13];
    double cos_th23 = fCosTh[OscPar::kTh23];
    double cos_th12 = fCosTh[OscPar::kTh12];
    double sin_th13 = fSinTh[OscPar::kTh13];
    
    double sin_2th23 = fSin2Th[OscPar::kTh23];
    double sin_2th12 = fSin2Th[OscPar::kTh12];
    double sin_2th13 = fSin2Th[OscPar::kTh13];
    
    double cos_2th13 = fCos2Th[OscPar::kTh13];
    double cos_2th12 = fCos2Th[OscPar::kTh12];
    double cos_2th23 = fCos2Th[OscPar::kTh23];
    
    double sinsq_th23 = fSinTh[OscPar::kTh23]*fSinTh[OscPar::kTh23];
    double sinsq_th12 = fSinTh[OscPar::kTh12]*fSinTh[OscPar::kTh12];
    double sinsq_th13 = fSinTh[OscPar::kTh13]*fSinTh[OscPar::kTh13];
    
    double d_cp = fOscPar[OscPar::kDelta];
    double cos_dcp = fCosDCP;
    double sin_dcp = fSinDCP;
    
    //Building the more complicated terms
    double Delta = dmsq_13*L/(4*E*1e9*OscPar::invKmToeV);
    double A = 2*fV*E*1e9/(dmsq_13);
    double alpha = dmsq_12/dmsq_13;
    
    // A and d_cp both change sign for antineutrinos
    double plusminus = int(fOscPar[OscPar::kNuAntiNu]);
    A *= plusminus;
    d_cp *= plusminus;
    sin_dcp *= plusminus;
    
    double t0 = pow(cos_th13,4)*sinsq_2th23*sin(Delta)*sin(Delta);
    double t11 = -1.0*Delta*cos_th13*cos_th13*sinsq_2th23*(cos_th12*cos_th12-sinsq_th12*sinsq_th13)*sin(2*Delta);
    double t12 = sin_th13*sin_2th12*sin_2th23*Delta*cos_th13*cos_th13*(2*sin_dcp*pow(sin(Delta),2)+cos_dcp*cos_2th23*sin(2*Delta));
    return t0+alpha*(t11+t12);
}
//----------------------------------------------------------------------------
double model3flav_vac_1stsin13_MuSurvive(double *x){

    
    double pmue = model3flav_vac_1stsin13_MuToElec(x);
    
    double pmutau = model3flav_vac_1stsin13_MuToElec(x);
    
    return 1.-pmue-pmutau;
}
//ElecToElec
double model3flav_vac_1stsin13_ElecToElec(double *x){
    double E = x[0]; //energy
    double L = fOscPar[OscPar::kL];  //baseline
    double dmsq_23 = fOscPar[OscPar::kDeltaM23];
    double dmsq_12 = fOscPar[OscPar::kDeltaM12];
    double dmsq_13 = dmsq_23+dmsq_12; //eV^{2}
    
    double sinsq_2th12 = fSin2Th[OscPar::kTh12]*fSin2Th[OscPar::kTh12];
    double sinsq_2th13 = fSin2Th[OscPar::kTh13]*fSin2Th[OscPar::kTh13];
    double sinsq_2th23 = fSin2Th[OscPar::kTh23]*fSin2Th[OscPar::kTh23];
    
    double cos_th13 = fCosTh[OscPar::kTh13];
    double cos_th23 = fCosTh[OscPar::kTh23];
    double cos_th12 = fCosTh[OscPar::kTh12];
    double sin_th13 = fSinTh[OscPar::kTh13];
    
    double sin_2th23 = fSin2Th[OscPar::kTh23];
    double sin_2th12 = fSin2Th[OscPar::kTh12];
    double sin_2th13 = fSin2Th[OscPar::kTh13];
    
    double cos_2th13 = fCos2Th[OscPar::kTh13];
    double cos_2th12 = fCos2Th[OscPar::kTh12];
    double cos_2th23 = fCos2Th[OscPar::kTh23];
    
    double sinsq_th23 = fSinTh[OscPar::kTh23]*fSinTh[OscPar::kTh23];
    double sinsq_th12 = fSinTh[OscPar::kTh12]*fSinTh[OscPar::kTh12];
    double sinsq_th13 = fSinTh[OscPar::kTh13]*fSinTh[OscPar::kTh13];
    
    double d_cp = fOscPar[OscPar::kDelta];
    double cos_dcp = fCosDCP;
    double sin_dcp = fSinDCP;
    
    //Building the more complicated terms
    double Delta = dmsq_13*L/(4*E*1e9*OscPar::invKmToeV);
    double A = 2*fV*E*1e9/(dmsq_13);
    double alpha = dmsq_12/dmsq_13;
    
    // A and d_cp both change sign for antineutrinos
    double plusminus = int(fOscPar[OscPar::kNuAntiNu]);
    A *= plusminus;
    d_cp *= plusminus;
    sin_dcp *= plusminus;

    double tee0 = 1-sinsq_2th13*pow(sin(Delta),2);
    double tee1 = Delta* sinsq_2th12*sinsq_2th13*sin(2*Delta);
    return tee0+alpha*tee1;
}
//ElecToMu
double model3flav_vac_1stsin13_ElecToMu(double *x){
    double E = x[0]; //energy
    double L = fOscPar[OscPar::kL];  //baseline
    double dmsq_23 = fOscPar[OscPar::kDeltaM23];
    double dmsq_12 = fOscPar[OscPar::kDeltaM12];
    double dmsq_13 = dmsq_23+dmsq_12; //eV^{2}
    
    double sinsq_2th12 = fSin2Th[OscPar::kTh12]*fSin2Th[OscPar::kTh12];
    double sinsq_2th13 = fSin2Th[OscPar::kTh13]*fSin2Th[OscPar::kTh13];
    double sinsq_2th23 = fSin2Th[OscPar::kTh23]*fSin2Th[OscPar::kTh23];
    
    double cos_th13 = fCosTh[OscPar::kTh13];
    double cos_th23 = fCosTh[OscPar::kTh23];
    double cos_th12 = fCosTh[OscPar::kTh12];
    double sin_th13 = fSinTh[OscPar::kTh13];
    
    double sin_2th23 = fSin2Th[OscPar::kTh23];
    double sin_2th12 = fSin2Th[OscPar::kTh12];
    double sin_2th13 = fSin2Th[OscPar::kTh13];
    
    double cos_2th13 = fCos2Th[OscPar::kTh13];
    double cos_2th12 = fCos2Th[OscPar::kTh12];
    double cos_2th23 = fCos2Th[OscPar::kTh23];
    
    double sinsq_th23 = fSinTh[OscPar::kTh23]*fSinTh[OscPar::kTh23];
    double sinsq_th12 = fSinTh[OscPar::kTh12]*fSinTh[OscPar::kTh12];
    double sinsq_th13 = fSinTh[OscPar::kTh13]*fSinTh[OscPar::kTh13];
    
    double d_cp = fOscPar[OscPar::kDelta];
    double cos_dcp = fCosDCP;
    double sin_dcp = fSinDCP;
    
    //Building the more complicated terms
    double Delta = dmsq_13*L/(4*E*1e9*OscPar::invKmToeV);
    double A = 2*fV*E*1e9/(dmsq_13);
    double alpha = dmsq_12/dmsq_13;
    
    // A and d_cp both change sign for antineutrinos
    double plusminus = int(fOscPar[OscPar::kNuAntiNu]);
    A *= plusminus;
    d_cp *= plusminus;
    sin_dcp *= plusminus;

    double temu0 = sinsq_2th13*sinsq_th23*pow(sin(Delta),2);
    double temu11 = -Delta*sinsq_th12*sinsq_2th13*sinsq_th23*sin(2*Delta);
    double temu12 = Delta*sin_2th12*sin_th13*pow(cos_th13,2)*sin_2th23*(2*sin_dcp*pow(sin(Delta),2)+cos_dcp*sin(2*Delta));
    return temu0+alpha*(temu11+temu12);
}
//ElecToTau
double model3flav_vac_1stsin13_ElecToTau(double *x){

    double etoe= model3flav_vac_1stsin13_ElecToElec(x);
    double etomu= model3flav_vac_1stsin13_ElecToMu(x);
    
    return 1.-(etoe+etomu);
}
//First order in alpha
//http://iopscience.iop.org/1126-6708/2004/04/078
//MuToElec
double model3flav_vac_1stalpha_MuToElec(double *x){
    double E = x[0]; //energy
    double L = fOscPar[OscPar::kL];  //baseline
    double dmsq_23 = fOscPar[OscPar::kDeltaM23];
    double dmsq_12 = fOscPar[OscPar::kDeltaM12];
    double dmsq_13 = dmsq_23+dmsq_12; //eV^{2}
    
    double sinsq_2th12 = fSin2Th[OscPar::kTh12]*fSin2Th[OscPar::kTh12];
    double sinsq_2th13 = fSin2Th[OscPar::kTh13]*fSin2Th[OscPar::kTh13];
    double sinsq_2th23 = fSin2Th[OscPar::kTh23]*fSin2Th[OscPar::kTh23];
    
    double cos_th13 = fCosTh[OscPar::kTh13];
    double cos_th23 = fCosTh[OscPar::kTh23];
    double cos_th12 = fCosTh[OscPar::kTh12];
    double sin_th13 = fSinTh[OscPar::kTh13];
    
    double sin_2th23 = fSin2Th[OscPar::kTh23];
    double sin_2th12 = fSin2Th[OscPar::kTh12];
    double sin_2th13 = fSin2Th[OscPar::kTh13];
    
    double cos_2th13 = fCos2Th[OscPar::kTh13];
    double cos_2th12 = fCos2Th[OscPar::kTh12];
    double cos_2th23 = fCos2Th[OscPar::kTh23];
    
    double sinsq_th23 = fSinTh[OscPar::kTh23]*fSinTh[OscPar::kTh23];
    double sinsq_th12 = fSinTh[OscPar::kTh12]*fSinTh[OscPar::kTh12];
    double sinsq_th13 = fSinTh[OscPar::kTh13]*fSinTh[OscPar::kTh13];
    
    double d_cp = fOscPar[OscPar::kDelta];
    double cos_dcp = fCosDCP;
    double sin_dcp = fSinDCP;
    
    //Building the more complicated terms
    double Delta = dmsq_13*L/(4*E*1e9*OscPar::invKmToeV);
    double A = 2*fV*E*1e9/(dmsq_13);
    double alpha = dmsq_12/dmsq_13;
    
    // A and d_cp both change sign for antineutrinos
    double plusminus = int(fOscPar[OscPar::kNuAntiNu]);
    A *= plusminus;
    d_cp *= plusminus;
    sin_dcp *= plusminus;
    
    //first term
    double tmue0 = sinsq_2th13*sinsq_th23*pow(sin(Delta),2);
    //second term
    double tmue11=-Delta*sinsq_th12*sinsq_2th13*sinsq_th23*sin(2*Delta);
    //change the side here for deltaCP for mue, in the paper, that is emu probability
    double tmue12=Delta*sin_2th12*sin_th13*pow(cos_th13,2)*sin_2th23*(-2*sin_dcp*pow(sin(Delta),2)+cos_dcp*sin(2*Delta));
    
    double mu2e = tmue0+alpha*(tmue11+tmue12);
    return mu2e;
}
//MuToTau
double model3flav_vac_1stalpha_MuToTau(double *x){
    return 0.0;
}
//MuSurvive
double model3flav_vac_1stalpha_MuSurvive(double *x){
    return 0.0;
}
//ElecToElec
double model3flav_vac_1stalpha_ElecToElec(double *x){
    return 0.0;
}
//ElecToMu
double model3flav_vac_1stalpha_ElecToMu(double *x){
    return 0.0;
}
//ElecToTau
double model3flav_vac_1stalpha_ElecToTau(double *x){
    return 0.0;
}
//====================================================================================
//Three-flavor vacuum function
//This is used for fiting
double model3flav_vac_apx_nocp_MuToElec(double *x){
    double E = x[0]; //energy
    double L = fOscPar[OscPar::kL];  //baseline
    double dmsq_23 = fOscPar[OscPar::kDeltaM23];
    double dmsq_12 = fOscPar[OscPar::kDeltaM12];
    double dmsq_13 = dmsq_23+dmsq_12; //eV^{2}
    
    double sinsq_2th13 = fSin2Th[OscPar::kTh13]*fSin2Th[OscPar::kTh13];
    
    double sinsq_th23 = fSinTh[OscPar::kTh23]*fSinTh[OscPar::kTh23];
 
    
    //Building the more complicated terms
    double Delta = dmsq_13*L/(4*E*1e9*OscPar::invKmToeV);
    double mu2e = sinsq_2th13*sinsq_th23*pow(sin(Delta),2);
    return mu2e;
    
}

double model3flav_vac_apx_nocp_MuSurvive(double *x){
    double E = x[0]; //energy
    double L = fOscPar[OscPar::kL];  //baseline
    double dmsq_23 = fOscPar[OscPar::kDeltaM23];
    double dmsq_12 = fOscPar[OscPar::kDeltaM12];
    double dmsq_13 = dmsq_23+dmsq_12; //eV^{2}
    
    double cos_th13 = fCosTh[OscPar::kTh13];
    double sinsq_2th23 = fSin2Th[OscPar::kTh23]*fSin2Th[OscPar::kTh23];
    
    
    //Building the more complicated terms
    double Delta = dmsq_13*L/(4*E*1e9*OscPar::invKmToeV);
    double mu2mu = 1. -pow(cos_th13,2)*sinsq_2th23*pow(sin(Delta),2);
    return mu2mu;

    
}
double model3flav_vac_apx_nocp_MuToTau(double *x){

    
    double mu2e = model3flav_vac_apx_nocp_MuToElec(x);
    double mu2mu = model3flav_vac_apx_nocp_MuSurvive(x);
    return 1-mu2mu-mu2e;
    
    
}

double model3flav_vac_apx_nocp_ElecToElec(double *x){
    double E = x[0]; //energy
    double L = fOscPar[OscPar::kL];  //baseline
    double dmsq_23 = fOscPar[OscPar::kDeltaM23];
    double dmsq_12 = fOscPar[OscPar::kDeltaM12];
    double dmsq_13 = dmsq_23+dmsq_12; //eV^{2}
    
    double sinsq_2th13 = fSin2Th[OscPar::kTh13]*fSin2Th[OscPar::kTh13];
        
    //Building the more complicated terms
    double Delta = dmsq_13*L/(4*E*1e9*OscPar::invKmToeV);
    double e2e = 1-sinsq_2th13*pow(sin(Delta),2);
    return e2e;

    
}
double model3flav_vac_apx_nocp_ElecToMu(double *x){
   //CP equals to zeo, e2mu = mu2e
    double E = x[0]; //energy
    double L = fOscPar[OscPar::kL];  //baseline
    double dmsq_23 = fOscPar[OscPar::kDeltaM23];
    double dmsq_12 = fOscPar[OscPar::kDeltaM12];
    double dmsq_13 = dmsq_23+dmsq_12; //eV^{2}
    
    double sinsq_2th13 = fSin2Th[OscPar::kTh13]*fSin2Th[OscPar::kTh13];
    
    double sinsq_th23 = fSinTh[OscPar::kTh23]*fSinTh[OscPar::kTh23];
    
    
    //Building the more complicated terms
    double Delta = dmsq_13*L/(4*E*1e9*OscPar::invKmToeV);
    double mu2e = sinsq_2th13*sinsq_th23*pow(sin(Delta),2);
    return mu2e;
    
}
double model3flav_vac_apx_nocp_ElecToTau(double *x){
    double e2e = model3flav_vac_apx_nocp_ElecToElec(x);
    double mu2e = model3flav_vac_apx_nocp_ElecToMu(x);
    return 1. -e2e-mu2e;

    
}
//====================================================================================
//====================================================================================
//second order in alpha and s13
double model3flav_vac_2nd_MuToElec(double *x){
    double E = x[0]; //energy
    double L = fOscPar[OscPar::kL];  //baseline
    double dmsq_23 = fOscPar[OscPar::kDeltaM23];
    double dmsq_12 = fOscPar[OscPar::kDeltaM12];
    double dmsq_13 = dmsq_23+dmsq_12; //eV^{2}
    
    double sinsq_2th12 = fSin2Th[OscPar::kTh12]*fSin2Th[OscPar::kTh12];
    double sinsq_2th13 = fSin2Th[OscPar::kTh13]*fSin2Th[OscPar::kTh13];
    double sinsq_2th23 = fSin2Th[OscPar::kTh23]*fSin2Th[OscPar::kTh23];
    
    double cos_th13 = fCosTh[OscPar::kTh13];
    double cos_th23 = fCosTh[OscPar::kTh23];
    double cos_th12 = fCosTh[OscPar::kTh12];
    double sin_th13 = fSinTh[OscPar::kTh13];
    
    double sin_2th23 = fSin2Th[OscPar::kTh23];
    double sin_2th12 = fSin2Th[OscPar::kTh12];
    double sin_2th13 = fSin2Th[OscPar::kTh13];
    
    double cos_2th13 = fCos2Th[OscPar::kTh13];
    double cos_2th12 = fCos2Th[OscPar::kTh12];
    double cos_2th23 = fCos2Th[OscPar::kTh23];
    
    double sinsq_th23 = fSinTh[OscPar::kTh23]*fSinTh[OscPar::kTh23];
    double sinsq_th12 = fSinTh[OscPar::kTh12]*fSinTh[OscPar::kTh12];
    double sinsq_th13 = fSinTh[OscPar::kTh13]*fSinTh[OscPar::kTh13];
    
    double d_cp = fOscPar[OscPar::kDelta];
    double cos_dcp = fCosDCP;
    double sin_dcp = fSinDCP;
    
    //Building the more complicated terms
    double Delta = dmsq_13*L/(4*E*1e9*OscPar::invKmToeV);
    double A = 2*fV*E*1e9/(dmsq_13);
    double alpha = dmsq_12/dmsq_13;
    
    // A and d_cp both change sign for antineutrinos
    double plusminus = int(fOscPar[OscPar::kNuAntiNu]);
    A *= plusminus;
    d_cp *= plusminus;
    sin_dcp *= plusminus;
    //change the sign of d_cp here
    double mu2e = pow(alpha,2)*sinsq_2th12*pow(cos_th23,2)*pow(Delta,2)+4*sinsq_th13*sinsq_th23*pow(sin(Delta),2) + 2*alpha*sin_th13*sin_2th12*sin_2th23*cos(Delta+d_cp)*Delta*sin(Delta);
    return mu2e;
}
double model3flav_vac_2nd_MuToTau(double *x){
    double E = x[0]; //energy
    double L = fOscPar[OscPar::kL];  //baseline
    double dmsq_23 = fOscPar[OscPar::kDeltaM23];
    double dmsq_12 = fOscPar[OscPar::kDeltaM12];
    double dmsq_13 = dmsq_23+dmsq_12; //eV^{2}
    
    double sinsq_2th12 = fSin2Th[OscPar::kTh12]*fSin2Th[OscPar::kTh12];
    double sinsq_2th13 = fSin2Th[OscPar::kTh13]*fSin2Th[OscPar::kTh13];
    double sinsq_2th23 = fSin2Th[OscPar::kTh23]*fSin2Th[OscPar::kTh23];
    
    double cos_th13 = fCosTh[OscPar::kTh13];
    double cos_th23 = fCosTh[OscPar::kTh23];
    double cos_th12 = fCosTh[OscPar::kTh12];
    double sin_th13 = fSinTh[OscPar::kTh13];
    
    double sin_2th23 = fSin2Th[OscPar::kTh23];
    double sin_2th12 = fSin2Th[OscPar::kTh12];
    double sin_2th13 = fSin2Th[OscPar::kTh13];
    
    double cos_2th13 = fCos2Th[OscPar::kTh13];
    double cos_2th12 = fCos2Th[OscPar::kTh12];
    double cos_2th23 = fCos2Th[OscPar::kTh23];
    
    double sinsq_th23 = fSinTh[OscPar::kTh23]*fSinTh[OscPar::kTh23];
    double sinsq_th12 = fSinTh[OscPar::kTh12]*fSinTh[OscPar::kTh12];
    double sinsq_th13 = fSinTh[OscPar::kTh13]*fSinTh[OscPar::kTh13];
    
    double d_cp = fOscPar[OscPar::kDelta];
    double cos_dcp = fCosDCP;
    double sin_dcp = fSinDCP;
    
    //Building the more complicated terms
    double Delta = dmsq_13*L/(4*E*1e9*OscPar::invKmToeV);
    double A = 2*fV*E*1e9/(dmsq_13);
    double alpha = dmsq_12/dmsq_13;
    
    // A and d_cp both change sign for antineutrinos
    double plusminus = int(fOscPar[OscPar::kNuAntiNu]);
    A *= plusminus;
    d_cp *= plusminus;
    sin_dcp *= plusminus;
    
    double mu2tau0 = sinsq_2th23*pow(sin(Delta),2)-alpha*pow(cos_th12,2)*sinsq_2th23*Delta*sin(2*Delta);
    double mu2tau1 = pow(alpha,2)*sinsq_2th23*pow(Delta,2)*(pow(cos_th12,4)*cos(2*Delta)-0.5*sinsq_2th12*pow(sin(Delta),2));
    double mu2tau2 = 2*sinsq_th13*sinsq_2th23*pow(sin(Delta),2);
    double mu2tau3 = 2*alpha*sin_th13*sin_2th12*sin_2th23*(sin_dcp*sin(Delta)-cos_2th23*cos_dcp*cos(Delta))*Delta*sin(Delta);
    double mu2tau = mu2tau0+mu2tau1+mu2tau2+mu2tau3;
    
    return mu2tau;
}
//can make double check here for mu survive
double model3flav_vac_2nd_MuSurvive(double *x){
    double E = x[0]; //energy
    double L = fOscPar[OscPar::kL];  //baseline
    double dmsq_23 = fOscPar[OscPar::kDeltaM23];
    double dmsq_12 = fOscPar[OscPar::kDeltaM12];
    double dmsq_13 = dmsq_23+dmsq_12; //eV^{2}
    
    double sinsq_2th12 = fSin2Th[OscPar::kTh12]*fSin2Th[OscPar::kTh12];
    double sinsq_2th13 = fSin2Th[OscPar::kTh13]*fSin2Th[OscPar::kTh13];
    double sinsq_2th23 = fSin2Th[OscPar::kTh23]*fSin2Th[OscPar::kTh23];
    
    double cos_th13 = fCosTh[OscPar::kTh13];
    double cos_th23 = fCosTh[OscPar::kTh23];
    double cos_th12 = fCosTh[OscPar::kTh12];
    double sin_th13 = fSinTh[OscPar::kTh13];
    
    double sin_2th23 = fSin2Th[OscPar::kTh23];
    double sin_2th12 = fSin2Th[OscPar::kTh12];
    double sin_2th13 = fSin2Th[OscPar::kTh13];
    
    double cos_2th13 = fCos2Th[OscPar::kTh13];
    double cos_2th12 = fCos2Th[OscPar::kTh12];
    double cos_2th23 = fCos2Th[OscPar::kTh23];
    
    double sinsq_th23 = fSinTh[OscPar::kTh23]*fSinTh[OscPar::kTh23];
    double sinsq_th12 = fSinTh[OscPar::kTh12]*fSinTh[OscPar::kTh12];
    double sinsq_th13 = fSinTh[OscPar::kTh13]*fSinTh[OscPar::kTh13];
    
    double d_cp = fOscPar[OscPar::kDelta];
    double cos_dcp = fCosDCP;
    double sin_dcp = fSinDCP;
    
    //Building the more complicated terms
    double Delta = dmsq_13*L/(4*E*1e9*OscPar::invKmToeV);
    double A = 2*fV*E*1e9/(dmsq_13);
    double alpha = dmsq_12/dmsq_13;
    
    // A and d_cp both change sign for antineutrinos
    double plusminus = int(fOscPar[OscPar::kNuAntiNu]);
    A *= plusminus;
    d_cp *= plusminus;
    sin_dcp *= plusminus;
    
    double mu2mu0 = 1-sinsq_2th23*pow(sin(Delta),2)+alpha*pow(cos_th12,2)*sinsq_2th23*Delta*sin(2*Delta);
    double mu2mu1 = -pow(alpha,2)*pow(Delta,2)*(sinsq_2th12*pow(cos_th23,2)+pow(cos_th12,2)*sinsq_2th23*(cos(2*Delta)-sinsq_th12));
    double mu2mu2 = 4*sinsq_th13*sinsq_th23*cos_2th23*pow(sin(Delta),2);
    double mu2mu3 = -2*alpha*sin_th13*sin_2th12*sinsq_th23*sin_2th23*cos_dcp*Delta*sin(2*Delta);
    double mu2mu = mu2mu0+mu2mu1+mu2mu2+mu2mu3;
    return mu2mu;
}
double model3flav_vac_2nd_ElecToElec(double *x){
    double E = x[0]; //energy
    double L = fOscPar[OscPar::kL];  //baseline
    double dmsq_23 = fOscPar[OscPar::kDeltaM23];
    double dmsq_12 = fOscPar[OscPar::kDeltaM12];
    double dmsq_13 = dmsq_23+dmsq_12; //eV^{2}
    
    double sinsq_2th12 = fSin2Th[OscPar::kTh12]*fSin2Th[OscPar::kTh12];
    double sinsq_2th13 = fSin2Th[OscPar::kTh13]*fSin2Th[OscPar::kTh13];
    double sinsq_2th23 = fSin2Th[OscPar::kTh23]*fSin2Th[OscPar::kTh23];
    
    double cos_th13 = fCosTh[OscPar::kTh13];
    double cos_th23 = fCosTh[OscPar::kTh23];
    double cos_th12 = fCosTh[OscPar::kTh12];
    double sin_th13 = fSinTh[OscPar::kTh13];
    
    double sin_2th23 = fSin2Th[OscPar::kTh23];
    double sin_2th12 = fSin2Th[OscPar::kTh12];
    double sin_2th13 = fSin2Th[OscPar::kTh13];
    
    double cos_2th13 = fCos2Th[OscPar::kTh13];
    double cos_2th12 = fCos2Th[OscPar::kTh12];
    double cos_2th23 = fCos2Th[OscPar::kTh23];
    
    double sinsq_th23 = fSinTh[OscPar::kTh23]*fSinTh[OscPar::kTh23];
    double sinsq_th12 = fSinTh[OscPar::kTh12]*fSinTh[OscPar::kTh12];
    double sinsq_th13 = fSinTh[OscPar::kTh13]*fSinTh[OscPar::kTh13];
    
    double d_cp = fOscPar[OscPar::kDelta];
    double cos_dcp = fCosDCP;
    double sin_dcp = fSinDCP;
    
    //Building the more complicated terms
    double Delta = dmsq_13*L/(4*E*1e9*OscPar::invKmToeV);
    double A = 2*fV*E*1e9/(dmsq_13);
    double alpha = dmsq_12/dmsq_13;
    
    // A and d_cp both change sign for antineutrinos
    double plusminus = int(fOscPar[OscPar::kNuAntiNu]);
    A *= plusminus;
    d_cp *= plusminus;
    sin_dcp *= plusminus;
    double e2e = 1-pow(alpha,2)*sinsq_2th12*pow(Delta,2)-4*sinsq_th13*pow(sin(Delta),2);
    return e2e;
}
double model3flav_vac_2nd_ElecToMu(double *x){
    //only change the sign of delta CP comparing with mu2e
    double E = x[0]; //energy
    double L = fOscPar[OscPar::kL];  //baseline
    double dmsq_23 = fOscPar[OscPar::kDeltaM23];
    double dmsq_12 = fOscPar[OscPar::kDeltaM12];
    double dmsq_13 = dmsq_23+dmsq_12; //eV^{2}
    
    double sinsq_2th12 = fSin2Th[OscPar::kTh12]*fSin2Th[OscPar::kTh12];
    double sinsq_2th13 = fSin2Th[OscPar::kTh13]*fSin2Th[OscPar::kTh13];
    double sinsq_2th23 = fSin2Th[OscPar::kTh23]*fSin2Th[OscPar::kTh23];
    
    double cos_th13 = fCosTh[OscPar::kTh13];
    double cos_th23 = fCosTh[OscPar::kTh23];
    double cos_th12 = fCosTh[OscPar::kTh12];
    double sin_th13 = fSinTh[OscPar::kTh13];
    
    double sin_2th23 = fSin2Th[OscPar::kTh23];
    double sin_2th12 = fSin2Th[OscPar::kTh12];
    double sin_2th13 = fSin2Th[OscPar::kTh13];
    
    double cos_2th13 = fCos2Th[OscPar::kTh13];
    double cos_2th12 = fCos2Th[OscPar::kTh12];
    double cos_2th23 = fCos2Th[OscPar::kTh23];
    
    double sinsq_th23 = fSinTh[OscPar::kTh23]*fSinTh[OscPar::kTh23];
    double sinsq_th12 = fSinTh[OscPar::kTh12]*fSinTh[OscPar::kTh12];
    double sinsq_th13 = fSinTh[OscPar::kTh13]*fSinTh[OscPar::kTh13];
    
    double d_cp = fOscPar[OscPar::kDelta];
    double cos_dcp = fCosDCP;
    double sin_dcp = fSinDCP;
    
    //Building the more complicated terms
    double Delta = dmsq_13*L/(4*E*1e9*OscPar::invKmToeV);
    double A = 2*fV*E*1e9/(dmsq_13);
    double alpha = dmsq_12/dmsq_13;
    
    // A and d_cp both change sign for antineutrinos
    double plusminus = int(fOscPar[OscPar::kNuAntiNu]);
    A *= plusminus;
    d_cp *= plusminus;
    sin_dcp *= plusminus;
    //change the sign of d_cp here
    double e2mu = pow(alpha,2)*sinsq_2th12*pow(cos_th23,2)*pow(Delta,2)+4*sinsq_th13*sinsq_th23*pow(sin(Delta),2) + 2*alpha*sin_th13*sin_2th12*sin_2th23*cos(Delta-d_cp)*Delta*sin(Delta);
    return e2mu;
}
double model3flav_vac_2nd_ElecToTau(double *x){
    //only change the sign of delta CP comparing with mu2e
    double E = x[0]; //energy
    double L = fOscPar[OscPar::kL];  //baseline
    double dmsq_23 = fOscPar[OscPar::kDeltaM23];
    double dmsq_12 = fOscPar[OscPar::kDeltaM12];
    double dmsq_13 = dmsq_23+dmsq_12; //eV^{2}
    
    double sinsq_2th12 = fSin2Th[OscPar::kTh12]*fSin2Th[OscPar::kTh12];
    double sinsq_2th13 = fSin2Th[OscPar::kTh13]*fSin2Th[OscPar::kTh13];
    double sinsq_2th23 = fSin2Th[OscPar::kTh23]*fSin2Th[OscPar::kTh23];
    
    double cos_th13 = fCosTh[OscPar::kTh13];
    double cos_th23 = fCosTh[OscPar::kTh23];
    double cos_th12 = fCosTh[OscPar::kTh12];
    double sin_th13 = fSinTh[OscPar::kTh13];
    
    double sin_2th23 = fSin2Th[OscPar::kTh23];
    double sin_2th12 = fSin2Th[OscPar::kTh12];
    double sin_2th13 = fSin2Th[OscPar::kTh13];
    
    double cos_2th13 = fCos2Th[OscPar::kTh13];
    double cos_2th12 = fCos2Th[OscPar::kTh12];
    double cos_2th23 = fCos2Th[OscPar::kTh23];
    
    double sinsq_th23 = fSinTh[OscPar::kTh23]*fSinTh[OscPar::kTh23];
    double sinsq_th12 = fSinTh[OscPar::kTh12]*fSinTh[OscPar::kTh12];
    double sinsq_th13 = fSinTh[OscPar::kTh13]*fSinTh[OscPar::kTh13];
    
    double d_cp = fOscPar[OscPar::kDelta];
    double cos_dcp = fCosDCP;
    double sin_dcp = fSinDCP;
    
    //Building the more complicated terms
    double Delta = dmsq_13*L/(4*E*1e9*OscPar::invKmToeV);
    double A = 2*fV*E*1e9/(dmsq_13);
    double alpha = dmsq_12/dmsq_13;
    
    // A and d_cp both change sign for antineutrinos
    double plusminus = int(fOscPar[OscPar::kNuAntiNu]);
    A *= plusminus;
    d_cp *= plusminus;
    sin_dcp *= plusminus;
    
    double e2tau0 = pow(alpha,2)*sinsq_2th12*sinsq_th23*pow(Delta,2)+4*sinsq_th13*pow(cos_th23,2)*pow(sin(Delta),2);
    double e2tau1 = -2*alpha*sin_th13*sin_2th12*sin_2th23*cos(Delta-d_cp)*Delta*sin(Delta);
    double e2tau = e2tau0+e2tau1;
    return e2tau;
}
//====================================================================================
//====================================================================================
//analytic prob from Joachim Kopp thesis
double model3flav_vac_kopp_MuToElec(double *x){
    double E = x[0]; //energy
    double L = fOscPar[OscPar::kL];  //baseline
    double dmsq_23 = fOscPar[OscPar::kDeltaM23];
    double dmsq_12 = fOscPar[OscPar::kDeltaM12];
    double dmsq_13 = dmsq_23+dmsq_12; //eV^{2}
    
    double sinsq_2th12 = fSin2Th[OscPar::kTh12]*fSin2Th[OscPar::kTh12];
    double sinsq_2th13 = fSin2Th[OscPar::kTh13]*fSin2Th[OscPar::kTh13];
    double sinsq_2th23 = fSin2Th[OscPar::kTh23]*fSin2Th[OscPar::kTh23];
    
    double cos_th13 = fCosTh[OscPar::kTh13];
    double cos_th23 = fCosTh[OscPar::kTh23];
    double cos_th12 = fCosTh[OscPar::kTh12];
    double sin_th13 = fSinTh[OscPar::kTh13];
    
    double sin_2th23 = fSin2Th[OscPar::kTh23];
    double sin_2th12 = fSin2Th[OscPar::kTh12];
    double sin_2th13 = fSin2Th[OscPar::kTh13];
    
    double cos_2th13 = fCos2Th[OscPar::kTh13];
    double cos_2th12 = fCos2Th[OscPar::kTh12];
    double cos_2th23 = fCos2Th[OscPar::kTh23];
    
    double sinsq_th23 = fSinTh[OscPar::kTh23]*fSinTh[OscPar::kTh23];
    double sinsq_th12 = fSinTh[OscPar::kTh12]*fSinTh[OscPar::kTh12];
    double sinsq_th13 = fSinTh[OscPar::kTh13]*fSinTh[OscPar::kTh13];
    
    double d_cp = fOscPar[OscPar::kDelta];
    double cos_dcp = fCosDCP;
    double sin_dcp = fSinDCP;
    
    //Building the more complicated terms
    double Delta = dmsq_13*L/(4*E*1e9*OscPar::invKmToeV);
    double A = 2*fV*E*1e9/(dmsq_13);
    double alpha = dmsq_12/dmsq_13;
    
    // A and d_cp both change sign for antineutrinos
    double plusminus = int(fOscPar[OscPar::kNuAntiNu]);
    A *= plusminus;
    d_cp *= plusminus;
    sin_dcp *= plusminus;
    double mu2e0 = sinsq_2th12*pow(cos_th23,2)*pow(cos_th13,2)*pow(sin(alpha*Delta),2);
    double mu2e1 = -0.5*sin_2th12*sin_2th13*sin_2th23*cos_th13*sin(alpha*Delta);
    double mu2e1fact = sin((alpha-2)*Delta-d_cp)+sin_dcp*cos(alpha*Delta)-cos_2th12*cos_dcp*sin(alpha*Delta);
    double mu2e2 = 0.25*sinsq_2th13*sinsq_th23;
    double mu2e2fact = 2-sinsq_2th12*pow(sin(alpha*Delta),2)-2*pow(cos_th12,2)*cos(2*Delta)-2*sinsq_th12*cos(2*(alpha-1)*Delta);
    double mu2e = mu2e0+ mu2e1*mu2e1fact+mu2e2*mu2e2fact;
    if (mu2e>1.0) {
        return 1.0;
    }
    else if (mu2e<1e-4) return 1e-4;
    else return mu2e;
}
//page 32 in Kopp thesis
//modify cos2 to sinsq
double model3flav_vac_kopp_ElecToElec(double *x){
    double E = x[0]; //energy
    double L = fOscPar[OscPar::kL];  //baseline
    double dmsq_23 = fOscPar[OscPar::kDeltaM23];
    double dmsq_12 = fOscPar[OscPar::kDeltaM12];
    double dmsq_13 = dmsq_23+dmsq_12; //eV^{2}
    
    double sinsq_2th12 = fSin2Th[OscPar::kTh12]*fSin2Th[OscPar::kTh12];
    double sinsq_2th13 = fSin2Th[OscPar::kTh13]*fSin2Th[OscPar::kTh13];
    double sinsq_2th23 = fSin2Th[OscPar::kTh23]*fSin2Th[OscPar::kTh23];
    
    double cos_th13 = fCosTh[OscPar::kTh13];
    double cos_th23 = fCosTh[OscPar::kTh23];
    double cos_th12 = fCosTh[OscPar::kTh12];
    double sin_th13 = fSinTh[OscPar::kTh13];
    
    double sin_2th23 = fSin2Th[OscPar::kTh23];
    double sin_2th12 = fSin2Th[OscPar::kTh12];
    double sin_2th13 = fSin2Th[OscPar::kTh13];
    
    double cos_2th13 = fCos2Th[OscPar::kTh13];
    double cos_2th12 = fCos2Th[OscPar::kTh12];
    double cos_2th23 = fCos2Th[OscPar::kTh23];
    
    double sinsq_th23 = fSinTh[OscPar::kTh23]*fSinTh[OscPar::kTh23];
    double sinsq_th12 = fSinTh[OscPar::kTh12]*fSinTh[OscPar::kTh12];
    double sinsq_th13 = fSinTh[OscPar::kTh13]*fSinTh[OscPar::kTh13];
    
    double d_cp = fOscPar[OscPar::kDelta];
    double cos_dcp = fCosDCP;
    double sin_dcp = fSinDCP;
    
    //Building the more complicated terms
    double Delta = dmsq_13*L/(4*E*1e9*OscPar::invKmToeV);
    double A = 2*fV*E*1e9/(dmsq_13);
    double alpha = dmsq_12/dmsq_13;
    
    // A and d_cp both change sign for antineutrinos
    double plusminus = int(fOscPar[OscPar::kNuAntiNu]);
    A *= plusminus;
    d_cp *= plusminus;
    sin_dcp *= plusminus;
    double e2e0 = -1.0*pow(cos_th13,4)*sinsq_2th12*pow(sin(alpha*Delta),2);
    double e2e1 = -1.0*pow(cos_th12,2)*sinsq_2th13*pow(sin(Delta),2);
    double e2e2 = -1.0*sinsq_th12*sinsq_2th13*pow(sin((alpha-1)*Delta),2);
    
    double mu2e = 1+e2e0+e2e1+e2e2;
    if (mu2e>1.0) {
        return 1.0;
    }
    else if (mu2e<1e-4) return 1e-4;
    else return mu2e;
}

double model2flav_mat_kopp_sol_MuToElec(double *x){
    double E = x[0]; //energy
    double L = fOscPar[OscPar::kL];  //baseline
    double dmsq_23 = fOscPar[OscPar::kDeltaM23];
    double dmsq_12 = fOscPar[OscPar::kDeltaM12];
    double dmsq_13 = dmsq_23+dmsq_12; //eV^{2}
    
    double sinsq_2th12 = fSin2Th[OscPar::kTh12]*fSin2Th[OscPar::kTh12];
    double sinsq_2th13 = fSin2Th[OscPar::kTh13]*fSin2Th[OscPar::kTh13];
    double sinsq_2th23 = fSin2Th[OscPar::kTh23]*fSin2Th[OscPar::kTh23];
    
    double cos_th13 = fCosTh[OscPar::kTh13];
    double cos_th23 = fCosTh[OscPar::kTh23];
    double cos_th12 = fCosTh[OscPar::kTh12];
    double sin_th13 = fSinTh[OscPar::kTh13];
    
    double sin_2th23 = fSin2Th[OscPar::kTh23];
    double sin_2th12 = fSin2Th[OscPar::kTh12];
    double sin_2th13 = fSin2Th[OscPar::kTh13];
    
    double cos_2th13 = fCos2Th[OscPar::kTh13];
    double cos_2th12 = fCos2Th[OscPar::kTh12];
    double cos_2th23 = fCos2Th[OscPar::kTh23];
    
    double sinsq_th23 = fSinTh[OscPar::kTh23]*fSinTh[OscPar::kTh23];
    double sinsq_th12 = fSinTh[OscPar::kTh12]*fSinTh[OscPar::kTh12];
    double sinsq_th13 = fSinTh[OscPar::kTh13]*fSinTh[OscPar::kTh13];
    
    double d_cp = fOscPar[OscPar::kDelta];
    double cos_dcp = fCosDCP;
    double sin_dcp = fSinDCP;
    
    //Building the more complicated terms
    double Delta = dmsq_13*L/(4*E*1e9*OscPar::invKmToeV);
    double A = 2*fV*E*1e9/(dmsq_13);
    double alpha = dmsq_12/dmsq_13;
    
    // A and d_cp both change sign for antineutrinos
    double plusminus = int(fOscPar[OscPar::kNuAntiNu]);
    A *= plusminus;
    d_cp *= plusminus;
    sin_dcp *= plusminus;
    //Now calculate the resonance terms for alpha expansion (C13) and s13 expansion (C12)
    double C13 = TMath::Sqrt(sinsq_2th13+(A-cos_2th13)*(A-cos_2th13));
    
    double C12 = 1;  //really C12 -> infinity when alpha = 0 but not an option really
    if(fabs(alpha) > 1e-10){  //want to be careful here
        double temp = cos_2th12 - A/alpha;
        C12 = TMath::Sqrt(sinsq_2th12+(temp*temp));
    }
    
    //More complicated sin and cosine terms
    double cosC13Delta = cos(C13*Delta);
    double sinC13Delta = sin(C13*Delta);
    
    double cosaC12Delta = 0;
    double sinaC12Delta = 0;
    
    if(fabs(alpha) > 1e-10){
        cosaC12Delta = cos(alpha*C12*Delta);
        sinaC12Delta = sin(alpha*C12*Delta);
    }
    double mu2e =pow(cos_th23,2)*sinsq_2th12*pow(sinaC12Delta,2)/pow(C12,2);
    if (mu2e>1.) {
        return 1.0;
    }
    else if (mu2e<1e-8) return 1e-4;
    return mu2e;
                           
}
//dominant by the atmospheric term (neglec solar term)
//
double model2flav_mat_kopp_atm_MuToElec(double *x){
    
    double E = x[0]; //energy
    double L = fOscPar[OscPar::kL];  //baseline
    double dmsq_23 = fOscPar[OscPar::kDeltaM23];
    double dmsq_12 = fOscPar[OscPar::kDeltaM12];
    double dmsq_13 = dmsq_23+dmsq_12; //eV^{2}
    
    double sinsq_2th12 = fSin2Th[OscPar::kTh12]*fSin2Th[OscPar::kTh12];
    double sinsq_2th13 = fSin2Th[OscPar::kTh13]*fSin2Th[OscPar::kTh13];
    double sinsq_2th23 = fSin2Th[OscPar::kTh23]*fSin2Th[OscPar::kTh23];
    
    double cos_th13 = fCosTh[OscPar::kTh13];
    double cos_th23 = fCosTh[OscPar::kTh23];
    double cos_th12 = fCosTh[OscPar::kTh12];
    double sin_th13 = fSinTh[OscPar::kTh13];
    
    double sin_2th23 = fSin2Th[OscPar::kTh23];
    double sin_2th12 = fSin2Th[OscPar::kTh12];
    double sin_2th13 = fSin2Th[OscPar::kTh13];
    
    double cos_2th13 = fCos2Th[OscPar::kTh13];
    double cos_2th12 = fCos2Th[OscPar::kTh12];
    double cos_2th23 = fCos2Th[OscPar::kTh23];
    
    double sinsq_th23 = fSinTh[OscPar::kTh23]*fSinTh[OscPar::kTh23];
    double sinsq_th12 = fSinTh[OscPar::kTh12]*fSinTh[OscPar::kTh12];
    double sinsq_th13 = fSinTh[OscPar::kTh13]*fSinTh[OscPar::kTh13];
    
    double d_cp = fOscPar[OscPar::kDelta];
    double cos_dcp = fCosDCP;
    double sin_dcp = fSinDCP;
    
    //Building the more complicated terms
    double Delta = dmsq_13*L/(4*E*1e9*OscPar::invKmToeV);
    double A = 2*fV*E*1e9/(dmsq_13);
    double alpha = dmsq_12/dmsq_13;
    
    // A and d_cp both change sign for antineutrinos
    double plusminus = int(fOscPar[OscPar::kNuAntiNu]);
    A *= plusminus;
    d_cp *= plusminus;
    sin_dcp *= plusminus;
    //Now calculate the resonance terms for alpha expansion (C13) and s13 expansion (C12)
    double C13 = TMath::Sqrt(sinsq_2th13+(A-cos_2th13)*(A-cos_2th13));
    
    double C12 = 1;  //really C12 -> infinity when alpha = 0 but not an option really
    if(fabs(alpha) > 1e-10){  //want to be careful here
        double temp = cos_2th12 - A/alpha;
        C12 = TMath::Sqrt(sinsq_2th12+(temp*temp));
    }
    
    //More complicated sin and cosine terms
    double cosC13Delta = cos(C13*Delta);
    double sinC13Delta = sin(C13*Delta);
    double sinaC13Delta =0;
    if(fabs(alpha) > 1e-10) sinaC13Delta = sin(alpha*C13*Delta);
    
    double cosaC12Delta = 0;
    double sinaC12Delta = 0;
    
    if(fabs(alpha) > 1e-10){
        cosaC12Delta = cos(alpha*C12*Delta);
        sinaC12Delta = sin(alpha*C12*Delta);
    }
    double mu2e =sinsq_th23*sinsq_2th13*pow(sinaC13Delta,2)/pow(C13,2);
    //double mu2e =sinsq_th23*sinsq_2th13*pow(sinC13Delta,2)/pow(C13,2);
    //return mu2e;
    if (mu2e>1.) {
        return 1.0;
    }
    else if (mu2e<1e-4) return 1e-4;
    else return mu2e;

    
}
double model3flav_mat_kopp_2nd_MuToElec(double *x){
    
    double E = x[0]; //energy
    double L = fOscPar[OscPar::kL];  //baseline
    double dmsq_23 = fOscPar[OscPar::kDeltaM23];
    double dmsq_12 = fOscPar[OscPar::kDeltaM12];
    double dmsq_13 = dmsq_23+dmsq_12; //eV^{2}
    
    double sinsq_2th12 = fSin2Th[OscPar::kTh12]*fSin2Th[OscPar::kTh12];
    double sinsq_2th13 = fSin2Th[OscPar::kTh13]*fSin2Th[OscPar::kTh13];
    double sinsq_2th23 = fSin2Th[OscPar::kTh23]*fSin2Th[OscPar::kTh23];
    
    double cos_th13 = fCosTh[OscPar::kTh13];
    double cos_th23 = fCosTh[OscPar::kTh23];
    double cos_th12 = fCosTh[OscPar::kTh12];
    double sin_th13 = fSinTh[OscPar::kTh13];
    
    double sin_2th23 = fSin2Th[OscPar::kTh23];
    double sin_2th12 = fSin2Th[OscPar::kTh12];
    double sin_2th13 = fSin2Th[OscPar::kTh13];
    
    double cos_2th13 = fCos2Th[OscPar::kTh13];
    double cos_2th12 = fCos2Th[OscPar::kTh12];
    double cos_2th23 = fCos2Th[OscPar::kTh23];
    
    double sinsq_th23 = fSinTh[OscPar::kTh23]*fSinTh[OscPar::kTh23];
    double sinsq_th12 = fSinTh[OscPar::kTh12]*fSinTh[OscPar::kTh12];
    double sinsq_th13 = fSinTh[OscPar::kTh13]*fSinTh[OscPar::kTh13];
    
    double d_cp = fOscPar[OscPar::kDelta];
    double cos_dcp = fCosDCP;
    double sin_dcp = fSinDCP;
    
    //Building the more complicated terms
    double Delta = dmsq_13*L/(4*E*1e9*OscPar::invKmToeV);
    double A = 2*fV*E*1e9/(dmsq_13);
    double alpha = dmsq_12/dmsq_13;
    
    // A and d_cp both change sign for antineutrinos
    double plusminus = int(fOscPar[OscPar::kNuAntiNu]);
    A *= plusminus;
    d_cp *= plusminus;
    sin_dcp *= plusminus;
    double mu2e0 = pow(alpha,2)*sinsq_2th12*pow(cos_th23,2)*pow(sin(A*Delta),2)/(pow(A,2));
    double mu2e1 = 4*sinsq_th13*sinsq_th23*pow(sin((A-1)*Delta),2)/pow((A-1),2);
    double mu2e2 = 2*alpha*sin_th13*sin_2th12*sin_2th23*cos(Delta+d_cp)*sin(A*Delta)*sin((A-1)*Delta)/(A*(A-1));
    double mu2e = mu2e0+mu2e1+mu2e2;
    if (mu2e>1) {
        cout<<"warning "<<mu2e<<" at E"<<E<<" and L"<<L<<endl;
    }
    //return mu2e;
    if (mu2e>1.0) {
        return 1.0;
    }
    else if (mu2e<1e-4) return 1e-4;
        else return mu2e;
}
double model3flav_vac_kopp_2nd_MuToElec(double *x){
    double E = x[0]; //energy
    double L = fOscPar[OscPar::kL];  //baseline
    double dmsq_23 = fOscPar[OscPar::kDeltaM23];
    double dmsq_12 = fOscPar[OscPar::kDeltaM12];
    double dmsq_13 = dmsq_23+dmsq_12; //eV^{2}
    
    double sinsq_2th12 = fSin2Th[OscPar::kTh12]*fSin2Th[OscPar::kTh12];
    double sinsq_2th13 = fSin2Th[OscPar::kTh13]*fSin2Th[OscPar::kTh13];
    double sinsq_2th23 = fSin2Th[OscPar::kTh23]*fSin2Th[OscPar::kTh23];
    
    double cos_th13 = fCosTh[OscPar::kTh13];
    double cos_th23 = fCosTh[OscPar::kTh23];
    double cos_th12 = fCosTh[OscPar::kTh12];
    double sin_th13 = fSinTh[OscPar::kTh13];
    
    double sin_2th23 = fSin2Th[OscPar::kTh23];
    double sin_2th12 = fSin2Th[OscPar::kTh12];
    double sin_2th13 = fSin2Th[OscPar::kTh13];
    
    double cos_2th13 = fCos2Th[OscPar::kTh13];
    double cos_2th12 = fCos2Th[OscPar::kTh12];
    double cos_2th23 = fCos2Th[OscPar::kTh23];
    
    double sinsq_th23 = fSinTh[OscPar::kTh23]*fSinTh[OscPar::kTh23];
    double sinsq_th12 = fSinTh[OscPar::kTh12]*fSinTh[OscPar::kTh12];
    double sinsq_th13 = fSinTh[OscPar::kTh13]*fSinTh[OscPar::kTh13];
    
    double d_cp = fOscPar[OscPar::kDelta];
    double cos_dcp = fCosDCP;
    double sin_dcp = fSinDCP;
    
    //Building the more complicated terms
    double Delta = dmsq_13*L/(4*E*1e9*OscPar::invKmToeV);
    double A = 2*fV*E*1e9/(dmsq_13);
    double alpha = dmsq_12/dmsq_13;
    
    // A and d_cp both change sign for antineutrinos
    double plusminus = int(fOscPar[OscPar::kNuAntiNu]);
    A *= plusminus;
    d_cp *= plusminus;
    sin_dcp *= plusminus;
    double mu2e0 = pow(alpha,2)*sinsq_2th12*pow(cos_th23,2);
    double mu2e1 = 4*sinsq_th13*sinsq_th23*pow(sin(Delta),2);
    double mu2e2 = 2*alpha*sin_th13*sin_2th12*sin_2th23*cos(Delta+d_cp)*sin(Delta);
    double mu2e = mu2e0+mu2e1+mu2e2;
    if (mu2e>1.0) {
        return 1.0;
    }
    else if (mu2e<1e-4) return 1e-4;
    else return mu2e;
   
}
//====================================================================================
//====================================================================================
//First order in sin13
double model3flav_mat_1stsin13_MuToElec(double *x){
    double E = x[0]; //energy
    
    //Building the standard terms
    double L = fOscPar[OscPar::kL];  //baseline
    double dmsq_23 = fOscPar[OscPar::kDeltaM23];
    double dmsq_12 = fOscPar[OscPar::kDeltaM12];
    double dmsq_13 = dmsq_23+dmsq_12; //eV^{2}
    
    double sinsq_2th12 = fSin2Th[OscPar::kTh12]*fSin2Th[OscPar::kTh12];
    double sinsq_2th13 = fSin2Th[OscPar::kTh13]*fSin2Th[OscPar::kTh13];
    
    double cos_th23 = fCosTh[OscPar::kTh23];
    double cos_th12 = fCosTh[OscPar::kTh12];
    double sin_th13 = fSinTh[OscPar::kTh13];
    
    double sin_2th23 = fSin2Th[OscPar::kTh23];
    double sin_2th12 = fSin2Th[OscPar::kTh12];
    double cos_2th13 = fCos2Th[OscPar::kTh13];
    double cos_2th12 = fCos2Th[OscPar::kTh12];
    
    double sinsq_th23 = fSinTh[OscPar::kTh23]*fSinTh[OscPar::kTh23];
    double sinsq_th12 = fSinTh[OscPar::kTh12]*fSinTh[OscPar::kTh12];
    
    double d_cp = fOscPar[OscPar::kDelta];
    double cos_dcp = fCosDCP;
    double sin_dcp = fSinDCP;
    
    //Building the more complicated terms
    double Delta = dmsq_13*L/(4*E*1e9*OscPar::invKmToeV);
    double A = 2*fV*E*1e9/(dmsq_13);
    double alpha = dmsq_12/dmsq_13;
    
    // A and d_cp both change sign for antineutrinos
    double plusminus = int(fOscPar[OscPar::kNuAntiNu]);
    A *= plusminus;
    d_cp *= plusminus;
    sin_dcp *= plusminus;
    
    //Now calculate the resonance terms for alpha expansion (C13) and s13 expansion (C12)
    
    double C13 = TMath::Sqrt(sinsq_2th13+(A-cos_2th13)*(A-cos_2th13));
    
    double C12 = 1;  //really C12 -> infinity when alpha = 0 but not an option really
    if(fabs(alpha) > 1e-10){  //want to be careful here
        double temp = cos_2th12 - A/alpha;
        C12 = TMath::Sqrt(sinsq_2th12+(temp*temp));
    }
    
    //More complicated sin and cosine terms
    double cosC13Delta = cos(C13*Delta);
    double sinC13Delta = sin(C13*Delta);
    double sin1pADelta = sin((A+1)*Delta);
    double cos1pADelta = cos((A+1)*Delta);
    double sinADelta = sin(A*Delta);
    double sinAm1Delta = sin((A-1)*Delta);
    double cosdpDelta = cos(d_cp+Delta);
    double sinApam2Delta = sin((A+alpha-2)*Delta);
    double cosApam2Delta = cos((A+alpha-2)*Delta);
    
    double cosaC12Delta = 0;
    double sinaC12Delta = 0;
    
    if(fabs(alpha) > 1e-10){
        cosaC12Delta = cos(alpha*C12*Delta);
        sinaC12Delta = sin(alpha*C12*Delta);
    } // otherwise not relevant
    
    
    // Now for the expansion in orders of sin(th13) (good to all order alpha)
    //  this is the equivalent of Eq 65 and 66
    
    // leading order term
    double pa1 = 0.0, pa2 = 0.0;
    
    if(fabs(alpha) > 1e-10){
        // leading order term
        pa1 = cos_th23*cos_th23*sinsq_2th12*sinaC12Delta*sinaC12Delta/(C12*C12);
        
        // and now to calculate the first order in s13 term
        double t1 = (cos_2th12 - A/alpha)/C12
        - alpha*A*C12*sinsq_2th12/(2*(1-alpha)*C12*C12);
        double t2 = -cos_dcp*(sinApam2Delta-sinaC12Delta*t1);
        
        double t3 = -(cosaC12Delta-cosApam2Delta)*sin_dcp;
        
        double denom = (1-A-alpha+A*alpha*cos_th12*cos_th12)*C12;
        double t4 = sin_2th12*sin_2th23*(1-alpha)*sinaC12Delta/denom;
        
        pa2 = t4*(t3+t2)*sin_th13;
    }
   
    
    //  Calculate the total probability
    double mu2e1stth13 = (pa1+pa2);
    if (mu2e1stth13>1.0) return 1.0;
    else if (mu2e1stth13<1e-4) return 1e-4;
    else return mu2e1stth13;
}
double model3flav_mat_1stsin13_MuToTau(double *x){
    return 0.0;
}
double model3flav_mat_1stsin13_MuSurvive(double *x){
    return 0.0;
}
double model3flav_mat_1stsin13_ElecToElec(double *x){
    return 0.0;
}
double model3flav_mat_1stsin13_ElecToMu(double *x){
    return 0.0;
}
double model3flav_mat_1stsin13_ElecToTau(double *x){
    return 0.0;
}
//====================================================================================
//====================================================================================
//First order in alpha
double model3flav_mat_1stalpha_MuToElec(double *x){
    double E = x[0]; //energy
    
    //Building the standard terms
    double L = fOscPar[OscPar::kL];  //baseline
    double dmsq_23 = fOscPar[OscPar::kDeltaM23];
    double dmsq_12 = fOscPar[OscPar::kDeltaM12];
    double dmsq_13 = dmsq_23+dmsq_12; //eV^{2}
    
    double sinsq_2th12 = fSin2Th[OscPar::kTh12]*fSin2Th[OscPar::kTh12];
    double sinsq_2th13 = fSin2Th[OscPar::kTh13]*fSin2Th[OscPar::kTh13];
    
    double cos_th23 = fCosTh[OscPar::kTh23];
    double cos_th12 = fCosTh[OscPar::kTh12];
    double sin_th13 = fSinTh[OscPar::kTh13];
    
    double sin_2th23 = fSin2Th[OscPar::kTh23];
    double sin_2th12 = fSin2Th[OscPar::kTh12];
    double cos_2th13 = fCos2Th[OscPar::kTh13];
    double cos_2th12 = fCos2Th[OscPar::kTh12];
    
    double sinsq_th23 = fSinTh[OscPar::kTh23]*fSinTh[OscPar::kTh23];
    double sinsq_th12 = fSinTh[OscPar::kTh12]*fSinTh[OscPar::kTh12];
    
    double d_cp = fOscPar[OscPar::kDelta];
    double cos_dcp = fCosDCP;
    double sin_dcp = fSinDCP;
    
    //Building the more complicated terms
    double Delta = dmsq_13*L/(4*E*1e9*OscPar::invKmToeV);
    double A = 2*fV*E*1e9/(dmsq_13);
    double alpha = dmsq_12/dmsq_13;
    
    // A and d_cp both change sign for antineutrinos
    double plusminus = int(fOscPar[OscPar::kNuAntiNu]);
    A *= plusminus;
    d_cp *= plusminus;
    sin_dcp *= plusminus;
    
    //Now calculate the resonance terms for alpha expansion (C13) and s13 expansion (C12)
    
    double C13 = TMath::Sqrt(sinsq_2th13+(A-cos_2th13)*(A-cos_2th13));
    
    double C12 = 1;  //really C12 -> infinity when alpha = 0 but not an option really
    if(fabs(alpha) > 1e-10){  //want to be careful here
        double temp = cos_2th12 - A/alpha;
        C12 = TMath::Sqrt(sinsq_2th12+(temp*temp));
    }
    
    //More complicated sin and cosine terms
    double cosC13Delta = cos(C13*Delta);
    double sinC13Delta = sin(C13*Delta);
    double sin1pADelta = sin((A+1)*Delta);
    double cos1pADelta = cos((A+1)*Delta);
    double sinADelta = sin(A*Delta);
    double sinAm1Delta = sin((A-1)*Delta);
    double cosdpDelta = cos(d_cp+Delta);
    double sinApam2Delta = sin((A+alpha-2)*Delta);
    double cosApam2Delta = cos((A+alpha-2)*Delta);
    
    double cosaC12Delta = 0;
    double sinaC12Delta = 0;
    
    if(fabs(alpha) > 1e-10){
        cosaC12Delta = cos(alpha*C12*Delta);
        sinaC12Delta = sin(alpha*C12*Delta);
    } // otherwise not relevant
    
    //First we calculate the terms for the alpha expansion (good to all orders in th13)
    // this is the equivalent of Eq 47 & 48 corrected for Mu to E instead of E to Mu
    
    // Leading order term
    double p1 = sinsq_th23*sinsq_2th13*sinC13Delta*sinC13Delta/(C13*C13);
    
    // terms that appear at order alpha
    double p2Inner = Delta*cosC13Delta*(1-A*cos_2th13)/C13 -
    A*sinC13Delta*(cos_2th13-A)/(C13*C13);
    
    double p2 = -2*sinsq_th12*sinsq_th23*sinsq_2th13*sinC13Delta/(C13*C13)*p2Inner*alpha;
    
    double p3Inner = - sin_dcp*(cosC13Delta - cos1pADelta)*C13
    + cos_dcp*(C13*sin1pADelta - (1-A*cos_2th13)*sinC13Delta);
    
    double p3 = sin_2th12*sin_2th23*sin_th13*sinC13Delta/(A*C13*C13)*p3Inner*alpha;
    
    double mu2e1stalpha = p1+p2+p3;
    
    if (mu2e1stalpha>1.0) return 1.0;
    else if (mu2e1stalpha<1e-4) return 1e-4;
    else return mu2e1stalpha;

}
double model3flav_mat_1stalpha_MuToTau(double *x){
    return 0.0;
}
double model3flav_mat_1stalpha_MuSurvive(double *x){
    return 0.0;
}
double model3flav_mat_1stalpha_ElecToElec(double *x){
    return 0.0;
}
double model3flav_mat_1stalpha_ElecToMu(double *x){
    return 0.0;
}
double model3flav_mat_1stalpha_ElecToTau(double *x){
    return 0.0;
}
//====================================================================================
//====================================================================================
// From Perdue presentation

double perdue3flav_mat_MuToElec(double *x){
    double E = x[0]; //energy
    
    //Building the standard terms
    double L = fOscPar[OscPar::kL];  //baseline
    double dmsq_23 = fOscPar[OscPar::kDeltaM23];
    double dmsq_12 = fOscPar[OscPar::kDeltaM12];
    double dmsq_13 = dmsq_23+dmsq_12; //eV^{2}
    
    double sinsq_2th12 = fSin2Th[OscPar::kTh12]*fSin2Th[OscPar::kTh12];
    double sinsq_2th13 = fSin2Th[OscPar::kTh13]*fSin2Th[OscPar::kTh13];
    
    double cos_th23 = fCosTh[OscPar::kTh23];
    double cos_th12 = fCosTh[OscPar::kTh12];
    double sin_th13 = fSinTh[OscPar::kTh13];
    
    double sin_2th23 = fSin2Th[OscPar::kTh23];
    double sin_2th12 = fSin2Th[OscPar::kTh12];
    double cos_2th13 = fCos2Th[OscPar::kTh13];
    double cos_2th12 = fCos2Th[OscPar::kTh12];
    
    double sinsq_th23 = fSinTh[OscPar::kTh23]*fSinTh[OscPar::kTh23];
    double sinsq_th12 = fSinTh[OscPar::kTh12]*fSinTh[OscPar::kTh12];
    
    double d_cp = fOscPar[OscPar::kDelta];
    double cos_dcp = fCosDCP;
    double sin_dcp = fSinDCP;
    
    //Building the more complicated terms
    double Delta = dmsq_13*L/(4*E*1e9*OscPar::invKmToeV);
    double A = 2*fV*E*1e9/(dmsq_13);
    double alpha = dmsq_12/dmsq_13;
    
    // A and d_cp both change sign for antineutrinos
    double plusminus = int(fOscPar[OscPar::kNuAntiNu]);
    A *= plusminus;
    d_cp *= plusminus;
    sin_dcp *= plusminus;
    
    //matter effect
    double a = 4000.;  // km
    a *= plusminus;
    
    
    double D31=1.27*dmsq_13*L/E;
    double D21=1.27*dmsq_12*L/E;
    double D32=1.27*dmsq_23*L/E;
    double matter =1.0;
    
    double atm =  sinsq_th23 *  sinsq_2th13   * pow(sin(D31  -  matter*L/a) , 2) *  pow( 1./(1. - matter*L/(a*D31))  ,  2);//32
    //should be possitive
    double sol = pow(cos_th23,2)   * sinsq_2th12  * pow(sin(D21 + matter*(L/a - D21)) , 2) * pow(1. + matter*(D21/(L/a) - 1.), 2);
    //correlation term
    double interf = 2.*sqrt(atm*sol)*(cos(D32) * cos_dcp - plusminus*sin(D32) * sin_dcp);
    return atm +sol + interf;
    
}

