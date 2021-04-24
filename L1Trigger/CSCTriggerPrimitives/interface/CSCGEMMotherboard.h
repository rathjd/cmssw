#ifndef L1Trigger_CSCTriggerPrimitives_CSCGEMMotherboard_h
#define L1Trigger_CSCTriggerPrimitives_CSCGEMMotherboard_h

/** \class CSCGEMMotherboard
 *
 * Base class for TMBs for the GEM-CSC integrated local trigger. Inherits
 * from CSCUpgradeMotherboard. Provides common functionality to match
 * ALCT/CLCT to GEM single clusters or coincidences of clusters
 *
 * \author Sven Dildick (TAMU)
 *
 */

#include "L1Trigger/CSCTriggerPrimitives/interface/CSCUpgradeMotherboard.h"
#include "L1Trigger/CSCTriggerPrimitives/interface/GEMCoPadProcessor.h"
#include "L1Trigger/CSCTriggerPrimitives/interface/CSCGEMMatcher.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "DataFormats/GEMDigi/interface/GEMPadDigiCollection.h"
#include "DataFormats/GEMDigi/interface/GEMCoPadDigiCollection.h"

class CSCGEMMotherboard : public CSCUpgradeMotherboard {
public:
  typedef std::vector<GEMInternalCluster> GEMInternalClusters;

  // standard constructor
  CSCGEMMotherboard(unsigned endcap,
                    unsigned station,
                    unsigned sector,
                    unsigned subsector,
                    unsigned chamber,
                    const edm::ParameterSet& conf);

  ~CSCGEMMotherboard() override;

  // clear stored pads and copads
  void clear();

  // run TMB with GEM pad clusters as input
  virtual void run(const CSCWireDigiCollection* wiredc,
                   const CSCComparatorDigiCollection* compdc,
                   const GEMPadDigiClusterCollection* gemPads);

  /** additional processor for GEMs */
  std::unique_ptr<GEMCoPadProcessor> coPadProcessor;

  /// set CSC and GEM geometries for the matching needs
  void setGEMGeometry(const GEMGeometry* g) { gem_g = g; }

protected:

  /* correlate a pair of ALCTs with matched clusters or coclusters
     the output is up to two LCTs */
  void correlateLCTsGEM(const CSCALCTDigi& bestALCT,
                        const CSCALCTDigi& secondALCT,
                        const GEMInternalClusters& clusters,
                        CSCCorrelatedLCTDigi& lct1,
                        CSCCorrelatedLCTDigi& lct2) const;

  /* correlate a pair of CLCTs with matched clusters or coclusters
     the output is up to two LCTs */
  void correlateLCTsGEM(const CSCCLCTDigi& bestCLCT,
                        const CSCCLCTDigi& secondCLCT,
                        const GEMInternalClusters& clusters,
                        CSCCorrelatedLCTDigi& lct1,
                        CSCCorrelatedLCTDigi& lct2) const;

  /* correlate a pair of ALCTs and a pair of CLCTs with matched clusters or coclusters
     the output is up to two LCTs */
  void correlateLCTsGEM(const CSCALCTDigi& bestALCT,
                        const CSCALCTDigi& secondALCT,
                        const CSCCLCTDigi& bestCLCT,
                        const CSCCLCTDigi& secondCLCT,
                        const GEMInternalClusters& clusters,
                        CSCCorrelatedLCTDigi& lct1,
                        CSCCorrelatedLCTDigi& lct2) const;

  /*
   * General function to construct integrated stubs from CSC and GEM information.
   * Options are:
   * 1. ALCT-CLCT-GEMPad
   * 2. ALCT-CLCT-GEMCoPad
   * 3. ALCT-GEMCoPad
   * 4. CLCT-GEMCoPad
   */
  // last argument is the LCT number (1 or 2)
  CSCCorrelatedLCTDigi constructLCTsGEM(const CSCALCTDigi& alct,
                                        const CSCCLCTDigi& clct,
                                        const GEMInternalCluster& gem,
                                        int i) const;

  // get the clusters in handy container
  void processGEMClusters(const GEMPadDigiClusterCollection* pads);

  /** Chamber id (trigger-type labels). */
  unsigned gemId;
  const GEMGeometry* gem_g;

  // for convenience...
  std::map<int, GEMInternalClusters > clusters_;

  std::unique_ptr<CSCGEMMatcher> cscGEMMatcher_;

private:
  // promote ALCT-GEM pattern
  bool promoteALCTGEMpattern_;

  // promote ALCT-GEM quality
  bool promoteALCTGEMquality_;
  bool promoteCLCTGEMquality_;
  bool promoteCLCTGEMquality_ME1a_;
  bool promoteCLCTGEMquality_ME1b_;

  // Drop low quality stubs if they don't have GEMs
  bool dropLowQualityCLCTsNoGEMs_ME1a_;
  bool dropLowQualityCLCTsNoGEMs_ME1b_;
  bool dropLowQualityALCTsNoGEMs_ME1a_;
  bool dropLowQualityALCTsNoGEMs_ME1b_;
  // drop low quality stubs if they don't have GEMs
  bool dropLowQualityCLCTsNoGEMs_;
  bool dropLowQualityALCTsNoGEMs_;

  // build LCT from ALCT and GEM
  bool buildLCTfromALCTandGEM_ME1a_;
  bool buildLCTfromALCTandGEM_ME1b_;
  bool buildLCTfromCLCTandGEM_ME1a_;
  bool buildLCTfromCLCTandGEM_ME1b_;

  // build LCT from ALCT and GEM
  bool buildLCTfromALCTandGEM_;
  bool buildLCTfromCLCTandGEM_;
};

#endif
