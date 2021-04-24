#ifndef L1Trigger_CSCTriggerPrimitives_CSCUpgradeMotherboard_h
#define L1Trigger_CSCTriggerPrimitives_CSCUpgradeMotherboard_h

/** \class CSCUpgradeMotherboard
 *
 * Base class for upgrade TMBs (MEX/1) chambers, that either run the
 * upgrade CSC-only TMB algorithm or the CSC-GEM algorithm
 *
 * \author Sven Dildick (TAMU)
 *
 */

#include "L1Trigger/CSCTriggerPrimitives/interface/CSCMotherboard.h"
#include "L1Trigger/CSCTriggerPrimitives/interface/CSCUpgradeAnodeLCTProcessor.h"
#include "L1Trigger/CSCTriggerPrimitives/interface/CSCUpgradeCathodeLCTProcessor.h"
#include "L1Trigger/CSCTriggerPrimitives/interface/LCTContainer.h"
#include "L1Trigger/CSCTriggerPrimitives/interface/CSCOverlap.h"

class CSCUpgradeMotherboard : public CSCMotherboard {
public:
  // standard constructor
  CSCUpgradeMotherboard(unsigned endcap,
                        unsigned station,
                        unsigned sector,
                        unsigned subsector,
                        unsigned chamber,
                        const edm::ParameterSet& conf);

  ~CSCUpgradeMotherboard() override {}

  // Empty the LCT container
  void clear();

  // run TMB with GEM pad clusters as input
  void run(const CSCWireDigiCollection* wiredc, const CSCComparatorDigiCollection* compdc) override;

  /* Returns vector of read-out correlated LCTs, if any.  Starts with
     the vector of all found LCTs and selects the ones in the read-out
     time window.
     When ME1/a is disabled, LCTs from ME1/a are rejected
  */
  std::vector<CSCCorrelatedLCTDigi> readoutLCTs() const override;

protected:
  void correlateLCTs(const CSCALCTDigi& bestALCT,
                     const CSCALCTDigi& secondALCT,
                     const CSCCLCTDigi& bestCLCT,
                     const CSCCLCTDigi& secondCLCT,
                     CSCCorrelatedLCTDigi& lct1,
                     CSCCorrelatedLCTDigi& lct2) const;

  // special cases for ME1/1 (when GEMs are not used)
  bool doesALCTCrossCLCT(const CSCALCTDigi& a, const CSCCLCTDigi& c) const;

  // special correlation function for ME1/1
  void correlateLCTsME11(const CSCALCTDigi& bALCT,
                         const CSCALCTDigi& sALCT,
                         const CSCCLCTDigi& bCLCT,
                         const CSCCLCTDigi& sCLCT,
                         CSCCorrelatedLCTDigi& lct1,
                         CSCCorrelatedLCTDigi& lct2) const;

  // Cross-BX sorting algorithm
  // where the central match BX goes first,
  // then the closest early, the closest late, etc.
  void crossBxSorting();

  // Compare two matches of type <ID,DIGI>
  // The template is match<GEMPadDigi> or match<GEMCoPadDigi>
  template <class S>
  bool compare(const S& p, const S& q) const;

  // Get the common matches of type <ID,DIGI>. Could be more than 1
  // The template is matches<GEMPadDigi> or matches<GEMCoPadDigi>
  template <class S>
  void intersection(const S& d1, const S& d2, S& result) const;

  /** for the case when more than 2 LCTs/BX are allowed;
      maximum match window = 15 */
  LCTContainer allLCTs;

  /** "preferential" index array in matching window for cross-BX sorting */
  int pref[CSCConstants::MAX_LCT_TBINS];

  bool match_earliest_alct_only;
  bool match_earliest_clct_only;

  /* type of algorithm to sort the stubs */
  bool tmb_cross_bx_sorting_;

  // debug gem matching
  bool debug_matching;

  // ignore unphysical ALCT-CLCT matches
  bool ignoreAlctCrossClct;

  std::unique_ptr<CSCOverlap> cscOverlap_;
};

template <class S>
bool CSCUpgradeMotherboard::compare(const S& p, const S& q) const {
  return (p.first == q.first) and (p.second == q.second);
}

template <class S>
void CSCUpgradeMotherboard::intersection(const S& d1, const S& d2, S& result) const {
  for (const auto& p : d1) {
    for (const auto& q : d2) {
      if (compare(p, q)) {
        result.push_back(p);
      }
    }
  }
}

#endif
