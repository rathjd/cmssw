#ifndef L1Trigger_CSCTriggerPrimitives_CSCOverlap
#define L1Trigger_CSCTriggerPrimitives_CSCOverlap

/** \class CSCOverlap
 *
 * Helper class to check if an ALCT overlaps with a CLCT in ME1/1
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

class CSCOverlap {
public:
  CSCOverlap(unsigned endcap, unsigned station, unsigned ring, bool isganged, const edm::ParameterSet& luts);

  // check if an ALCT can cross a CLCT. Not always the case for ME1/1
  bool doesALCTCrossCLCT(const CSCALCTDigi& a, const CSCCLCTDigi& c) const;

private:
  bool doesWiregroupCrossHalfStrip(int wg, int keystrip) const;

  unsigned endcap_;
  unsigned station_;
  unsigned ring_;
  bool isganged_;

  // strings to paths of LUTs
  std::vector<std::string> wgCrossHsME1aFiles_;
  std::vector<std::string> wgCrossHsME1aGangedFiles_;
  std::vector<std::string> wgCrossHsME1bFiles_;

  // unique pointers to the luts
  std::unique_ptr<CSCLUTReader> wg_cross_min_hs_ME1a_;
  std::unique_ptr<CSCLUTReader> wg_cross_max_hs_ME1a_;
  std::unique_ptr<CSCLUTReader> wg_cross_min_hs_ME1a_ganged_;
  std::unique_ptr<CSCLUTReader> wg_cross_max_hs_ME1a_ganged_;
  std::unique_ptr<CSCLUTReader> wg_cross_min_hs_ME1b_;
  std::unique_ptr<CSCLUTReader> wg_cross_max_hs_ME1b_;
};

#endif
