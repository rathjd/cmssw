#include "L1Trigger/CSCTriggerPrimitives/interface/CSCGEMMatcher.h"
#include "L1Trigger/CSCTriggerPrimitives/interface/CSCLUTReader.h"
#include "L1Trigger/CSCTriggerPrimitives/interface/GEMInternalCluster.h"
#include "DataFormats/CSCDigi/interface/CSCConstants.h"
#include "DataFormats/CSCDigi/interface/CSCALCTDigi.h"
#include "DataFormats/CSCDigi/interface/CSCCLCTDigi.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <algorithm>

CSCGEMMatcher::CSCGEMMatcher(unsigned endcap, unsigned station, unsigned chamber, const edm::ParameterSet& luts)
    : endcap_(endcap), station_(station), chamber_(chamber)
{

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

    // for coincidences or single clusters in L1
    if (cl.isCoincidence() or cl.id().layer() == 1) {
      key_es = cl.layer1_middle_es();
      if (station_ == 1 and clct.getKeyStrip() > CSCConstants::MAX_HALF_STRIP_ME1B)
        key_es = cl.layer1_middle_es_me1a();
    }

    // for single clusters in L2
    else if (cl.id().layer() == 2) {
      key_es = cl.layer2_middle_es();
      if (station_ == 1 and clct.getKeyStrip() > CSCConstants::MAX_HALF_STRIP_ME1B)
        key_es = cl.layer2_middle_es_me1a();
    }
    
    else edm::LogWarning("CSCGEMMatcher") << "cluster.id().layer =" << cl.id().layer() << " out of acceptable range 1-2!";

    // matching by 1/8-strip
    // need new function from Denis here!!!
    if (key_es <= clct.getKeyStrip(8) and clct.getKeyStrip(8) <= key_es)
      output.push_back(cl);
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
