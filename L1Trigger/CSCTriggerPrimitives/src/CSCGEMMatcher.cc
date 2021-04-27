#include "L1Trigger/CSCTriggerPrimitives/interface/CSCGEMMatcher.h"
#include "L1Trigger/CSCTriggerPrimitives/interface/CSCLUTReader.h"
#include "L1Trigger/CSCTriggerPrimitives/interface/GEMInternalCluster.h"
#include "DataFormats/CSCDigi/interface/CSCConstants.h"
#include "DataFormats/CSCDigi/interface/CSCALCTDigi.h"
#include "DataFormats/CSCDigi/interface/CSCCLCTDigi.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <algorithm>
#include <math.h>

CSCGEMMatcher::CSCGEMMatcher(int endcap, unsigned station, unsigned chamber, const edm::ParameterSet& luts)
    : endcap_(endcap), station_(station), chamber_(chamber)
{

}

int CSCGEMMatcher::CSCGEMSlopeCorrector(bool isFacing, bool isL1orCopad, int cscSlope) const {

  //Slope correction fit values for GEM to any CSC combinations
  float Erf[4][2]={{63., 9.6},{118., 7.4},{52., 9.5},{106., 7.1}}; //constants for facing layers
  float Lin[4][2]={{1.1, 5.8},{0.7, 16.7},{-1.3, 4.2},{-7.4, 14.}};//constants for off-side layers

  //determine shift by slope correction
  float SlopeShift=0.;
  int ChooseCorr = (chamber_ % 2) + (isL1orCopad ? 0 : 2); //chooses applicable constants
  if(isFacing) SlopeShift = Erf[ChooseCorr][0] * std::erf(cscSlope / Erf[ChooseCorr][1]);
  else         SlopeShift = Lin[ChooseCorr + 1][0] + Lin[ChooseCorr + 1][1] * cscSlope;
  
  return round(SlopeShift * endcap_);
}

// match an ALCT to GEMInternalCluster by bunch-crossing
CSCGEMMatcher::GEMInternalClusters CSCGEMMatcher::matchingClustersBX(const CSCALCTDigi& alct,
                                                                     const GEMInternalClusters& clusters) const {
  GEMInternalClusters output;

  if (!alct.isValid()) return output;

  // select clusters matched in time
  for (const auto& cl : clusters) {
    if (std::abs(alct.getBX() - cl.bx()) <= 1)
      output.push_back(cl);
  }

  return output;
}

// match a CLCT to GEMInternalCluster by bunch-crossing
CSCGEMMatcher::GEMInternalClusters CSCGEMMatcher::matchingClustersBX(const CSCCLCTDigi& clct,
                                                                     const GEMInternalClusters& clusters) const {
  GEMInternalClusters output;

  if (!clct.isValid()) return output;

  // select clusters matched in time
  for (const auto& cl : clusters) {
    if (std::abs(clct.getBX() - cl.bx()) <= 3)
      output.push_back(cl);
  }

  return output;
}

// match an ALCT and CLCT to GEMInternalCluster by bunch-crossing
CSCGEMMatcher::GEMInternalClusters CSCGEMMatcher::matchingClustersBX(const CSCALCTDigi& alct,
                                                                     const CSCCLCTDigi& clct,
                                                                     const GEMInternalClusters& clusters) const {
  GEMInternalClusters output;

  // both need to be valid
  if (!alct.isValid() or !clct.isValid()) return output;

  // get the single matches
  const auto& alctClusters = matchingClustersBX(alct, clusters);
  const auto& clctClusters = matchingClustersBX(clct, clusters);

  // get the intersection
  for (const auto& p : alctClusters) {
    for (const auto& q : clctClusters) {
      if (p == q) {
        output.push_back(p);
      }
    }
  }

  return output;
}

CSCGEMMatcher::GEMInternalClusters CSCGEMMatcher::matchingClustersLoc(const CSCALCTDigi& alct, const GEMInternalClusters& clusters) const{

  GEMInternalClusters output;

  if (!alct.isValid()) return output;

  // select clusters matched in wiregroup
  for (const auto& cl : clusters) {
    if (cl.min_wg() <= alct.getKeyWG() and alct.getKeyWG() <= cl.max_wg())
      output.push_back(cl);
  }

  return output;
}

CSCGEMMatcher::GEMInternalClusters CSCGEMMatcher::matchingClustersLoc(const CSCCLCTDigi& clct, const GEMInternalClusters& clusters) const{

  GEMInternalClusters output;

  if (!clct.isValid()) return output;

  // select clusters matched in wiregroup
  for (const auto& cl : clusters) {

    // key 1/8-strip
    int key_es = -1;
    
    //modification of DeltaStrip by CLCT slope
    int SlopeShift = 0;
    float clctSlope = pow(-1, clct.getBend()) * clct.getSlope();

    // for coincidences or single clusters in L1
    if (cl.isCoincidence() or cl.id().layer() == 1) {
      key_es = cl.layer1_middle_es();
      if (station_ == 1 and clct.getKeyStrip() > CSCConstants::MAX_HALF_STRIP_ME1B)
        key_es = cl.layer1_middle_es_me1a();

      //set SlopeShift for L1 or Copad case
      SlopeShift = CSCGEMSlopeCorrector(true, true, clctSlope); // currently fixed to facing detectors, must be determined at motherboard level
    }

    // for single clusters in L2
    else if (cl.id().layer() == 2) {
      key_es = cl.layer2_middle_es();
      if (station_ == 1 and clct.getKeyStrip() > CSCConstants::MAX_HALF_STRIP_ME1B)
        key_es = cl.layer2_middle_es_me1a();

      //set SlopeShift for L2 case
      SlopeShift = CSCGEMSlopeCorrector(true, false, clctSlope); // currently fixed to facing detectors, must be determined at motherboard level
      
    }
    
    else edm::LogWarning("CSCGEMMatcher") << "cluster.id().layer =" << cl.id().layer() << " out of acceptable range 1-2!";

    // matching by 1/8-strip
    // determine matching window by chamber, assuming facing chambers only are processed
    int window = chamber_ % 2 == 0 ? 20 : 40;

    if (abs(clct.getKeyStrip(8) - key_es + SlopeShift) < window)
      output.push_back(cl); // currently pushes back all clusters in window, might want to only push back closest?
  }



  return output;
}

CSCGEMMatcher::GEMInternalClusters CSCGEMMatcher::matchingClustersLoc(const CSCALCTDigi& alct,
                                                                      const CSCCLCTDigi& clct,
                                                                      const GEMInternalClusters& clusters) const
{
  GEMInternalClusters output;

  // both need to be valid
  if (!alct.isValid() or !clct.isValid()) return output;

  // get the single matches
  const auto& alctClusters = matchingClustersLoc(alct, clusters);
  const auto& clctClusters = matchingClustersLoc(clct, clusters);

  // get the intersection
  for (const auto& p : alctClusters) {
    for (const auto& q : clctClusters) {
      if (p == q) {
        output.push_back(p);
      }
    }
  }

  return output;
}

CSCGEMMatcher::GEMInternalClusters CSCGEMMatcher::matchingClustersBXLoc(const CSCALCTDigi& alct, const GEMInternalClusters& clusters) const
{
  // match by BX
 const auto& alctClustersBX = matchingClustersBX(alct, clusters);

 // match spatially
 const auto& alctClustersBXLoc = matchingClustersBX(alct, alctClustersBX);

 return alctClustersBXLoc;
}

CSCGEMMatcher::GEMInternalClusters CSCGEMMatcher::matchingClustersBXLoc(const CSCCLCTDigi& clct, const GEMInternalClusters& clusters) const
{
  // match by BX
 const auto& clctClustersBX = matchingClustersBX(clct, clusters);

 // match spatially
 const auto& clctClustersBXLoc = matchingClustersBX(clct, clctClustersBX);

 return clctClustersBXLoc;
}

CSCGEMMatcher::GEMInternalClusters CSCGEMMatcher::matchingClustersBXLoc(const CSCALCTDigi& alct,
                                            const CSCCLCTDigi& clct,
                                            const GEMInternalClusters& clusters) const
{
  // match by BX
  const auto& clustersBX = matchingClustersBX(alct, clct, clusters);

  // match spatially
  const auto& clustersBXLoc = matchingClustersBX(alct, clct, clustersBX);

  return clustersBXLoc;
}
