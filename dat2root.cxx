#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

void dat2root(const char* inFilename  = "MuC_test.txt",
                        const char* outFilename = "simc_out.root",
                        Double_t Ebeam         = 7000,    // [GeV]
                        Double_t mProton       = 0.938272 // [GeV]
                       )
{
   // Input text
   std::ifstream inFile(inFilename);
   if (!inFile.is_open()) {
      std::cerr << "Error: Could not open input file " << inFilename << std::endl;
      return;
   }

   // Output ROOT file, TTree
   TFile* outFile = new TFile(outFilename, "RECREATE");
   if (!outFile || outFile->IsZombie()) {
      std::cerr << "Error creating ROOT file " << outFilename << std::endl;
      return;
   }
   TTree* tree = new TTree("simc", "SIMC proton/electron events");

   // We'll store TLorentzVectors for the final-state proton/electron
   TLorentzVector lv_p; // (px_p, py_p, pz_p, E_p)
   TLorentzVector lv_e; // (px_e, py_e, pz_e, E_e)

   // Also store the vertex positions (vz_p, vz_e), and weight
   Double_t vz_p, vz_e, weight;

   // We'll store kinematic variables:
   Double_t Q2, W, yy, tAbs, xbj, mx; // Q^2, W, y, |t|, Bjorken x

   // Create branches
   tree->Branch("proton",    &lv_p);
   tree->Branch("electron",  &lv_e);
   tree->Branch("vz_p",      &vz_p,    "vz_p/D");
   tree->Branch("vz_e",      &vz_e,    "vz_e/D");
   tree->Branch("weight",    &weight,  "weight/D");

   tree->Branch("Q2",  &Q2,  "Q2/D");
   tree->Branch("W",   &W,   "W/D");
   tree->Branch("y",   &yy,  "y/D");
   tree->Branch("tAbs",&tAbs,"tAbs/D");
   tree->Branch("x",   &xbj, "x/D");
   tree->Branch("mx",   &mx, "mx/D");

   // We'll define the beam and target in the lab frame:
   // beam along +z with energy Ebeam
   TLorentzVector e_in(0.0, 0.0, Ebeam, Ebeam);
   // stationary proton target
   TLorentzVector p_target(0.0, 0.0, 0.0, mProton);

   // Read each line
   std::string line;
   while (true) {
      if (!std::getline(inFile, line)) break; // end-of-file
      if (line.size() < 2) continue;         // skip empty or short lines

      std::istringstream iss(line);
      // Expect 11 columns:
      // px_p  py_p  pz_p  E_p  vz_p  px_e  py_e  pz_e  E_e  vz_e  weight
      Double_t px_p_d, py_p_d, pz_p_d, E_p_d;
      Double_t px_e_d, py_e_d, pz_e_d, E_e_d;
      if (!(iss >> px_p_d >> py_p_d >> pz_p_d >> E_p_d >> vz_p
                >> px_e_d >> py_e_d >> pz_e_d >> E_e_d >> vz_e >> weight)) {
         // line not in expected format
         continue;
      }

      // Fill TLorentzVectors
      lv_p.SetPxPyPzE(px_p_d, py_p_d, pz_p_d, E_p_d);
      lv_e.SetPxPyPzE(px_e_d, py_e_d, pz_e_d, E_e_d);

      // Let's compute the 4-vector of momentum transfer: q = e_in - lv_e
      TLorentzVector q = e_in - lv_e;

      // Q^2 = - (q^2)
      Q2 = -(q.Mag2()); // prefer positive

      // W^2 = (p_target + q)^2
      TLorentzVector hadronic = p_target + q;
      Double_t W2 = hadronic.Mag2();
      if (W2 > 0.0) W = std::sqrt(W2);
      else          W = 0.0;

      // y = (p_target · q) / (p_target · e_in)
      //   = (E_p * nu - p_vec dot q_vec) / (E_p * Ebeam - p_vec dot e_in_vec)
      // but p_target is at rest => p_target dot e_in = mProton * Ebeam
      // and p_target dot q = mProton * nu, where nu = q.E()
      Double_t pDotEIn = mProton* e_in.E();
      Double_t pDotQ   = mProton* q.E();
      yy = pDotQ / pDotEIn;

      // x = Q^2 / (2 * p_target dot q)
      //   = Q^2 / (2 * mProton * nu)
      // be sure to check for nonzero
      Double_t nu = q.E();
      if (std::fabs(nu) < 1e-9) xbj = 0.0;
      else xbj = Q2 / (2.0*mProton*nu);

      // t = (p_target - p_final)^2 if the final proton is the p we see
      // assume p_final = lv_p, so t = (p_target - lv_p)^2
      TLorentzVector tvec = p_target - lv_p;
      Double_t tVal = tvec.Mag2(); // might be negative
      tAbs = std::fabs(tVal);

      mx = ((e_in+p_target) - (lv_p+lv_e)).M();

      // Fill the TTree
      tree->Fill();
   }

   // Write
   outFile->cd();
   tree->Write();
   outFile->Close();
   inFile.close();

   std::cout << "Converted " << inFilename << " => " << outFilename
             << ", TTree 'simc', storing proton/electron TLorentzVectors"
             << " plus kinematics (Q2, W, y, tAbs, x)." << std::endl;
}
