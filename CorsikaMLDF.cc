//Provided by Alex Tapia
#include <cmath>
#include <vector>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <iostream>
using namespace std;

#include "TF1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TProfile.h"

#include "corsika-tape.h"
#include "corsika-particle.h"


#define QUIET  1

#ifdef QUIET
 #define Q  0
#else
 #define Q  1
#endif

const unsigned int kBINS             = 140;

const double       kUndergroundDepth = 1.3;//2.5;//2.25;  /// [ m ]

const double       kSaveSeeds        = false;

struct mVector {
  /// Surface level
  std::vector<double> vX1;  
  std::vector<double> vY1;  
  /// Underground level
  std::vector<double> vX2;  
  std::vector<double> vY2;  
};

/// Perpendicular distance from position "x,y" to shower axis
double RPerp(const double x, const double y, const double theta, const double phi) {

  /// Axis versor
  const double nX   = std::sin(theta) * std::cos(phi); 
  const double nY   = std::sin(theta) * std::sin(phi);
  
  /// Scalar dot    
  const double scal = nX * x + nY * y;

  /// Squared position
  const double pos2 = x*x + y*y;

  return std::sqrt( pos2 - scal*scal );

} /// End RPerp()

double fitf(double *x, double *par) {

  /*
    The AMIGA muon detectors of the Pierre Auger Observatory: overview and status
    F. Suarez for the Pierre Auger Collaboration. 33RD INTERNATIONAL COSMIC RAY 
    CONFERENCE, RIO DE JANEIRO 2013
  */

  ///const double x0    = 150.;  /// [ m ]
  ///const double alpha = 1.00;
  ///const double gamma = 1.85;
 
  /*
     KASCADE-Grande Collaboration, Proc. 29th ICRC, Pune, India, 6 (2005) 301.
  */

  /// const double x0    = 320.;  /// [ m ]
  /// const double alpha = 0.69;
  /// const double gamma = 1.00;

  /*
    KASCADE-Grande Collaboration. See Supanitsky PhD Thesis. Detectores de 
    Superficie y la Composicion Quimica de los Rayos Cosmicos (2007).
  */

  const double x0    = 320.;  /// [ m ]
  const double alpha = 0.75;
  const double gamma = 2.93;

  return ( par[0] * 
           std::pow(x[0]/x0,-alpha) * 
           std::pow(1.+x[0]/x0,-par[1]) *
           std::pow(1. + std::pow(x[0]/(10.*x0),2),-gamma) );

  /*
    First showers seen by the AMIGA Muon Detector. F. Sanchez et al. See
    GAP-Note 2012-120.
  */

  /// const double x0 = 450.;  /// [ m ]

  /// return ( par[0] * std::pow(x[0]/x0,-par[1]) );

} /// End fitf() 

bool IsUnderground(const double energy, const double cosZenith) {

  /// Energy loss
  const double dEdX = 1.808e-7;     /// [ GeV m^2/g ]
  /// Rock density
  const double rho  = 1.8e+6;     /// [ g/m^3 ]
  
  double eUnderground = energy - dEdX * rho * kUndergroundDepth / cosZenith;
  
  if(eUnderground > 0.)
   return true;
  else
    return false;

} /// End IsUnderground()

int main(int argc, char **argv) {

  if(argc == 1) {
   std::cout << "\n Usage: ./CorsikaMLDF.exe CORSIKA file \n" << endl;
   std::exit(1);
  } /// End if()
    
  const double dr = 25.;  /// [ m ]

  double mOn[kBINS], mUn[kBINS], r_j[kBINS];
  for(unsigned int jBin = 0; jBin <= kBINS; ++jBin) { 
   mOn[jBin] = 0.;
   mUn[jBin] = 0.;
   r_j[jBin] = ((double) jBin) * dr; 
  } /// End for()
 
  std::vector<mVector> mVec;

  std::fstream outFile;
  if(kSaveSeeds) {
   outFile.open("seeds.dat",ios::out);
   outFile << "#" << endl;
  }

  unsigned int nFile = 0;
  for(unsigned int jIter = 1; jIter < (unsigned int) argc; ++jIter) {
   corsika_file  *tape;
   sub_block     *sb;
   particle_data *pd;

   RWORD desc;
     
   if((tape = corsika_fopen(argv[jIter], "r")) == NULL) {
    std::fprintf(stderr, "\n Cannot open file %s. Skiping file... \n", argv[jIter]);
    continue;
   } /// End if()

   if(std::strstr(argv[jIter],".tab") || 
      std::strstr(argv[jIter],".long") || 
      std::strstr(argv[jIter],".dbase") ||
      std::strstr(argv[jIter],".db") || 
      std::strstr(argv[jIter],".lst"))
      
    continue;

   nFile++;

   std::cout << "\n Processing file " << nFile << "\n" << endl;

   double zenith     = 0.;
   double azimuth    = 0.;
   double altitude   = 0.;  

   unsigned int jMax = 0;
   unsigned int seed[10];
   for(unsigned int j = 0; j < 10; ++j)
    seed[j] = 0;

   while((sb = corsika_get(tape)) != NULL) {
     if(is_event_header(sb->eh.id)) {
      zenith   = sb->eh.theta;                    /// [ rad ]
      azimuth  = sb->eh.phi;                      /// [ rad ]
      altitude = sb->eh.observation_height[ 0 ];  /// [ cm ]

      std::cout << "   zenith   = " << zenith  * 180. / 3.14159 << " deg " << endl;
      std::cout << "   azimuth  = " << azimuth * 180. / 3.14159 << " deg " << endl;
      std::cout << "   altitude = " << altitude / 100           << " m "   << endl;

      jMax = (unsigned int) sb->eh.random_sequences;
      for(unsigned int j = 0; j < jMax; ++j)
       seed[j] = sb->eh.seeds[j].seed; 
     } /// End if()

     if(is_event_end(sb->eh.id)) {
      std::cout << "   event id = " << sb->ee.event_number << endl; 
      std::cout << "   seeds    = "; 
      for(unsigned int j = 0; j < jMax; ++j)
       std::cout << " " << seed[j];
      std::cout << endl;

      outFile << " " << sb->ee.event_number;
      for(unsigned int j = 0; j < jMax; ++j)
       outFile << " " << seed[j];
      outFile << endl; 
    } /// End if()

     if(is_control(sb->rh.id)) {
      Q && std::printf("%4.4s\n", sb->rh.id);
     } 
     else {
      for(int i = 0; i < PARTICLES_IN_SUB_BLOCK; ++i) {
       pd   = sb->pb.particle + i;
       desc = pd->description;
       
       if(is_particle_r(desc)) {
        std::string s = particle[(int) desc/1000].name;

        /// Only mu+/mu-
        if(s != "mu+" && s != "mu-") 
         continue; 
       }
       else if(is_nucleus_r(desc))
         continue;
       else
         continue;
    
       if(desc != 0) {
        /* Surface level */
        double rPart = RPerp(pd->x, pd->y, zenith, azimuth) / 100;  /// [ m ]

        for(unsigned int jBin = 0; jBin < kBINS; ++jBin) {
         if(rPart >= r_j[jBin] && rPart < r_j[jBin+1]) {
#ifdef THINNING
	  mOn[jBin] += pd->weight;
#else
	  mOn[jBin] += 1.;
#endif
	  break;
         } /// End if()
        } /// End for()

        /* Underground level */
        /// Energy of the particle on the ground level
        const double energy   = std::sqrt(pd->px * pd->px + pd->py * pd->py + pd->pz * pd->pz);  /// [ GeV ]

        /// cos(theta) of the particle on the ground level
        const double cosTheta = pd->pz / energy;

        /// Here from before condition the energy > 0.
        if(IsUnderground(energy, cosTheta)) { 
         /// If pz = 0 then RPerp is infinity.
         if(!pd->pz)
          continue;

         const double xUnGrd = pd->x - 100. * kUndergroundDepth * pd->px / pd->pz;  /// [ cm ]
         const double yUnGrd = pd->y - 100. * kUndergroundDepth * pd->py / pd->pz;  /// [ cm ]

         const double rUnGrd = RPerp(xUnGrd, yUnGrd, zenith, azimuth) / 100.; 
         
         for(unsigned int jBin = 0; jBin < kBINS; ++jBin) { 
          if(rUnGrd >= r_j[jBin] && rUnGrd < r_j[jBin+1]) {
#ifdef THINNING
	   mUn[jBin] += pd->weight;
#else
	   mUn[jBin] += 1.;
#endif
	   break;
          } /// End if()
         } /// End for()
        } /// End if()
       } /// End if()
      } /// End for(i)
     } /// End else{}
   } /// End while()

   struct mVector vDummy;
   for(unsigned int jBin = 0; jBin < kBINS; ++jBin) {
    double area = 3.14159 * (r_j[jBin+1]*r_j[jBin+1] - r_j[jBin]*r_j[jBin]);

    vDummy.vX1.push_back( r_j[jBin] + dr/2. );
    vDummy.vY1.push_back( mOn[jBin] / area ); 

    vDummy.vX2.push_back( r_j[jBin] + dr/2. );
    vDummy.vY2.push_back( mUn[jBin] / area ); 
   } /// End for()
 
   mVec.push_back(vDummy);

   for(unsigned int jBin = 0; jBin <= kBINS; ++jBin) {
    mOn[jBin] = 0.;
    mUn[jBin] = 0.;
   } /// End for()
  } /// End for(jIter)

  if(kSaveSeeds) {
   outFile << "#" << endl;
   outFile.close();
  } 

  ///-----

  TLegend *lg = new TLegend(0.6611,0.6828,0.8389,0.8387,NULL,"brNDC");
  lg->SetBorderSize(1);
  lg->SetTextFont(62);
  lg->SetTextSize(0.03146);
  lg->SetLineColor(1);
  lg->SetLineStyle(1);
  lg->SetLineWidth(1);
  lg->SetFillColor(10);
  lg->SetFillStyle(1001);

  std::fstream outFile2;
  outFile2.open("fit.dat",ios::out);

  TFile *rootFile = new TFile("rootFile.root","RECREATE");

  /* PLOTTER */

  TCanvas  *cA[2]; 
  TProfile *prof[2];
  
  for(unsigned int nIter = 0; nIter < 2; ++nIter) {
   char c1[25];
   std::sprintf(c1,"cA%d", nIter + 1);

   cA[nIter] = new TCanvas(c1,c1,600,400);
   cA[nIter]->Clear();
   cA[nIter]->SetFillColor(10);

    char c2[25];
    std::sprintf(c2,"prof%d", nIter + 1);

    prof[nIter] = new TProfile(c2,c2,120,0.,3000.,0.,6000);
    prof[nIter]->Clear();

    prof[nIter]->SetMarkerStyle(8);
    prof[nIter]->SetMarkerColor(1); 
    prof[nIter]->SetMarkerSize(0.7);
   } /// End for()

   TGraph   *gr[2];

   for(unsigned int jIter = 0; jIter < (unsigned int) mVec.size(); ++jIter) {
    for(unsigned int nIter = 0; nIter < 2; ++nIter) {
     char c3[25];

     if(!nIter) {
      if(jIter < 9)
       std::sprintf(c3,"shower_0%d_surface", jIter + 1);
      else
       std::sprintf(c3,"shower_%d_surface", jIter + 1);

      for(unsigned int i = 0; i < (unsigned int) mVec[jIter].vX1.size(); ++i) {
       if(mVec[jIter].vY1[i] <= 0.) {
        mVec[jIter].vX1.erase(mVec[jIter].vX1.begin() + i);
        mVec[jIter].vY1.erase(mVec[jIter].vY1.begin() + i);
       } /// End if()
      } /// End for()
     }
     else {
       if(jIter < 9)
        std::sprintf(c3,"shower_0%d_underground", jIter + 1);
       else
        std::sprintf(c3,"shower_%d_underground", jIter + 1);
 
       for(unsigned int i = 0; i < (unsigned int) mVec[jIter].vX2.size(); ++i) {
        if(mVec[jIter].vY2[i] <= 0.) {
         mVec[jIter].vX2.erase(mVec[jIter].vX2.begin() + i);
         mVec[jIter].vY2.erase(mVec[jIter].vY2.begin() + i);
        } /// End if() 
       } /// End for()
     } /// End else{}

     if(!nIter)
      gr[0] = new TGraph(mVec[jIter].vX1.size(), &mVec[jIter].vX1.front(), &mVec[jIter].vY1.front());
     else
      gr[1] = new TGraph(mVec[jIter].vX2.size(), &mVec[jIter].vX2.front(), &mVec[jIter].vY2.front());

     gr[nIter]->SetMarkerColor(2);
     gr[nIter]->SetLineColor(2);


     TCanvas *cB = new TCanvas(c3,c3,600,400);
     cB->Clear();
     cB->SetFillColor(10);

     gPad->SetGrid();
     gPad->SetLogy();

     gr[nIter]->Draw("APL");

     char title[54];
     if(!nIter)
      std::sprintf(title,"Lateral distribution of #mu+/#mu- (surface level) - Shower = %d",jIter+1);
     else
      std::sprintf(title,"Lateral distribution of #mu+/#mu- (underground level) - Shower = %d",jIter+1);

     gr[nIter]->SetTitle(title);
     gr[nIter]->GetXaxis()->SetTitle(" r [ m ] ");
     gr[nIter]->GetYaxis()->SetTitle(" #rho_{#mu} [ 1/m^{2} ] ");
     gr[nIter]->GetXaxis()->CenterTitle();
     gr[nIter]->GetYaxis()->CenterTitle();
  
     gr[nIter]->SetMinimum(1.e-2);

     TF1 *func = new TF1("func",fitf,150.,2500.,2);
   
     gr[nIter]->Fit("func","NR");

     func->SetLineColor(1);
     func->Draw("SAME");
 
     if(!nIter)
       outFile2 << " SURFACE ";
     else
       outFile2 << " UNDERGROUND ";

     outFile2 << jIter+1 
              << "  " << func->GetParameter(0) 
              << "  " << func->GetParError(0)
              << "  " << func->GetParameter(1) 
              << "  " << func->GetParError(1)
              << endl;
 
     cB->Update();
     cB->Write();
  
     delete func;
     delete cB;
   
     ///-----

     cA[nIter]->cd();
   
     gPad->SetGrid();
     gPad->SetLogy();

     if(!jIter) {
      //gr[nIter]->Draw("APL");//ATC:10/11/2015

      gr[nIter]->SetTitle("Lateral distribution of #mu+/#mu-: Average distribution and fit");

      gr[nIter]->GetXaxis()->SetTitle(" r [ m ] ");
      gr[nIter]->GetYaxis()->SetTitle(" #rho_{#mu} [ 1/m^{2} ] ");
      gr[nIter]->GetXaxis()->CenterTitle();
      gr[nIter]->GetYaxis()->CenterTitle();  

      gr[nIter]->SetMinimum(1.e-2);
    
      //if(!nIter)
       //lg->AddEntry(gr[nIter],"showers","lp");//ATC:10/11/2015
     }
     else
       gr[nIter]->Draw("PL");

     cA[nIter]->Update();

     ///-----

     if(!nIter) {
      for(unsigned int i = 0; i < (unsigned int) mVec[jIter].vX1.size(); ++i)
       prof[0]->Fill(mVec[jIter].vX1[i], mVec[jIter].vY1[i]);
     }
     else 
       for(unsigned int i = 0; i < (unsigned int) mVec[jIter].vX2.size(); ++i)
        prof[1]->Fill(mVec[jIter].vX2[i], mVec[jIter].vY2[i]);
   } /// End for(nIter)
  } /// End for(jIter)

  ///-----

  for(unsigned int nIter = 0; nIter < 2; ++nIter) {
   cA[nIter]->cd();
   //prof[nIter]->Draw("SAME");//ATC:10/11/2015 
   prof[nIter]->Draw(""); 

   prof[nIter]->SetTitle("Lateral distribution of #mu+/#mu-: Average distribution and fit");

   prof[nIter]->GetXaxis()->SetTitle(" r [ m ] ");
   prof[nIter]->GetYaxis()->SetTitle(" #rho_{#mu} [ 1/m^{2} ] ");
   prof[nIter]->GetXaxis()->CenterTitle();
   prof[nIter]->GetYaxis()->CenterTitle();
   prof[nIter]->SetMarkerColor(kBlack);
   prof[nIter]->SetLineColor(kBlack);


   TF1 *func = new TF1("func",fitf,250.,2500.,2);
   prof[nIter]->Fit("func","NR");

   func->SetLineColor(4);
   func->Draw("SAME");

   if(!nIter) {
    lg->AddEntry(prof[nIter],"average","p");
    lg->AddEntry(func,"fit","lp");
   }

   lg->Draw("SAME");

   cA[nIter]->Update();
   cA[nIter]->Write();
 
   /// outFile << "## AVERAGE " << endl;
   /// outFile << "#" << endl;

   for(unsigned int jBin = 1; jBin <= (unsigned int) prof[nIter]->GetNbinsX(); ++jBin) {
    double x = (double) prof[nIter]->GetBinCenter(jBin);
    double y = (double) prof[nIter]->GetBinContent(jBin);
    double w = (double) prof[nIter]->GetBinError(jBin);

    /// if(y)
    ///  outFile << "  " << x << "  " << y << "  " << w << endl;
   } /// End for()

   /// outFile << "#" << endl;

   outFile2 << "#" << endl;

   if(!nIter)
     outFile2 << " AVERAGE-SURFACE ";
   else
     outFile2 << " AVERAGE-UNDERGROUND ";

   outFile2 << "  " << func->GetParameter(0) 
            << "  " << func->GetParError(0)
            << "  " << func->GetParameter(1) 
            << "  " << func->GetParError(1)
            << endl;

   outFile2 << "#" << endl;

   delete func; 
  } /// End for(nIter)

  outFile2.close();
 
  delete lg;

  for(unsigned int nIter = 0; nIter < 2; ++nIter) {
   delete prof[nIter]; prof[nIter] = 0;
   delete gr[nIter]; gr[nIter] = 0;
   delete cA[nIter]; cA[nIter] = 0;
  } /// End for()

  ///-----

  rootFile->Close();

  return 0;

} /// End main()
