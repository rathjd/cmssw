#ifndef L1Trigger_CSCTriggerPrimitives_GEMCoPadProcessor_h
#define L1Trigger_CSCTriggerPrimitives_GEMCoPadProcessor_h

/** \class GEMCoPadProcessor
 *
 * \author Sven Dildick (TAMU)
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/GEMDigi/interface/GEMPadDigiCollection.h"
#include "DataFormats/GEMDigi/interface/GEMPadDigiClusterCollection.h"
#include "DataFormats/GEMDigi/interface/GEMCoPadDigi.h"
#include "L1Trigger/CSCTriggerPrimitives/interface/GEMInternalCluster.h"
#include "L1Trigger/CSCTriggerPrimitives/interface/CSCLUTReader.h"

#include <vector>

class GEMCoPadProcessor {
public:
  /** Normal constructor. */
  GEMCoPadProcessor(unsigned region,
                    unsigned station,
                    unsigned chamber,
                    const edm::ParameterSet& copad,
                    const edm::ParameterSet& luts);

  /** Clear copad vector */
  void clear();

  /** Runs the CoPad processor code. */
  std::vector<GEMInternalCluster> run(const GEMPadDigiClusterCollection*);

private:
  void addSingleClusters(const GEMPadDigiClusterCollection*);

  void addCoincidenceClusters(const GEMPadDigiClusterCollection*);

  // translate the cluster central pad numbers into 1/8-strip number,
  // and roll numbers into min and max wiregroup numbers
  // for matching with CSC trigger primitives
  void doCoordinateConversion();

  /** Chamber id (trigger-type labels). */
  int theRegion;
  int theStation;
  int theChamber;
  bool isEven_;

  unsigned int maxDeltaPad_;
  unsigned int maxDeltaBX_;
  unsigned int maxDeltaRoll_;

  // output collection
  std::vector<GEMInternalCluster> clusters_;

  // strings to paths of LUTs
  std::vector<std::string> padToEsME1aFiles_;
  std::vector<std::string> padToEsME1bFiles_;
  std::vector<std::string> padToEsME21Files_;
  std::vector<std::string> rollToMaxWgME11Files_;
  std::vector<std::string> rollToMinWgME11Files_;
  std::vector<std::string> rollToMaxWgME21Files_;
  std::vector<std::string> rollToMinWgME21Files_;

  // unique pointers to the luts
  std::unique_ptr<CSCLUTReader> GEMCSCLUT_pad_es_ME1a_even_;
  std::unique_ptr<CSCLUTReader> GEMCSCLUT_pad_es_ME1a_odd_;
  std::unique_ptr<CSCLUTReader> GEMCSCLUT_pad_es_ME1b_even_;
  std::unique_ptr<CSCLUTReader> GEMCSCLUT_pad_es_ME1b_odd_;
  std::unique_ptr<CSCLUTReader> GEMCSCLUT_pad_es_ME21_even_;
  std::unique_ptr<CSCLUTReader> GEMCSCLUT_pad_es_ME21_odd_;

  std::unique_ptr<CSCLUTReader> GEMCSCLUT_roll_max_wg_ME11_even_;
  std::unique_ptr<CSCLUTReader> GEMCSCLUT_roll_max_wg_ME11_odd_;
  std::unique_ptr<CSCLUTReader> GEMCSCLUT_roll_min_wg_ME11_even_;
  std::unique_ptr<CSCLUTReader> GEMCSCLUT_roll_min_wg_ME11_odd_;
  std::unique_ptr<CSCLUTReader> GEMCSCLUT_roll_max_wg_ME21_even_;
  std::unique_ptr<CSCLUTReader> GEMCSCLUT_roll_max_wg_ME21_odd_;
  std::unique_ptr<CSCLUTReader> GEMCSCLUT_roll_min_wg_ME21_even_;
  std::unique_ptr<CSCLUTReader> GEMCSCLUT_roll_min_wg_ME21_odd_;
};

#endif
