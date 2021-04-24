#ifndef L1Trigger_CSCTriggerPrimitives_CSCGEMMatcher
#define L1Trigger_CSCTriggerPrimitives_CSCGEMMatcher

/** \class CSCGEMMatcher
 *
 * Helper class to check if an ALCT or CLCT matches with a GEMInternalCluster
 *
 * \author Sven Dildick (TAMU)
 *
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <string>
#include <vector>

class CSCLUTReader;
class CSCALCTDigi;
class CSCCLCTDigi;
class GEMInternalCluster;

class CSCGEMMatcher {
public:
  typedef std::vector<GEMInternalCluster> GEMInternalClusters;

  CSCGEMMatcher(unsigned endcap, unsigned station, unsigned chamber, const edm::ParameterSet& luts);

  // match by BX

  // coincidences
  GEMInternalClusters matchingClustersBX(const CSCALCTDigi& alct, const GEMInternalClusters& clusters) const;

  // coincidences
  GEMInternalClusters matchingClustersBX(const CSCCLCTDigi& clct, const GEMInternalClusters& clusters) const;

  // coincidences or single clusters
  GEMInternalClusters matchingClustersBX(const CSCALCTDigi& alct,
                                         const CSCCLCTDigi& clct,
                                         const GEMInternalClusters& clusters) const;

  // match by location

  // coincidences
  GEMInternalClusters matchingClustersLoc(const CSCALCTDigi& alct, const GEMInternalClusters& clusters) const;

  // coincidences
  GEMInternalClusters matchingClustersLoc(const CSCCLCTDigi& clct, const GEMInternalClusters& clusters) const;

  // coincidences or single clusters
  GEMInternalClusters matchingClustersLoc(const CSCALCTDigi& alct,
                                          const CSCCLCTDigi& clct,
                                          const GEMInternalClusters& clusters) const;

  // match by BX and location

  // coincidences
  GEMInternalClusters matchingClustersBXLoc(const CSCALCTDigi& alct, const GEMInternalClusters& clusters) const;

  // coincidences
  GEMInternalClusters matchingClustersBXLoc(const CSCCLCTDigi& clct, const GEMInternalClusters& clusters) const;

  // coincidences or single clusters
  GEMInternalClusters matchingClustersBXLoc(const CSCALCTDigi& alct,
                                            const CSCCLCTDigi& clct,
                                            const GEMInternalClusters& clusters) const;

private:
  unsigned endcap_;
  unsigned station_;
  unsigned ring_;
  unsigned chamber_;
};

#endif
