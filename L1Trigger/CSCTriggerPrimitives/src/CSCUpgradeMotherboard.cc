#include <memory>

#include "L1Trigger/CSCTriggerPrimitives/interface/CSCUpgradeMotherboard.h"

CSCUpgradeMotherboard::CSCUpgradeMotherboard(unsigned endcap,
                                             unsigned station,
                                             unsigned sector,
                                             unsigned subsector,
                                             unsigned chamber,
                                             const edm::ParameterSet& conf)
    :  // special configuration parameters for ME11 treatment
      CSCMotherboard(endcap, station, sector, subsector, chamber, conf),
      allLCTs(match_trig_window_size) {
  if (!runPhase2_)
    edm::LogError("CSCUpgradeMotherboard|SetupError") << "+++ TMB constructed while runPhase2_ is not set! +++\n";

  if (theRing == 1) {
    if (theStation == 1 and !runME11Up_)
      edm::LogError("CSCUpgradeMotherboard|SetupError") << "+++ TMB constructed while runME11Up_ is not set! +++\n";
    if (theStation == 2 and !runME21Up_)
      edm::LogError("CSCUpgradeMotherboard|SetupError") << "+++ TMB constructed while runME21Up_ is not set! +++\n";
    if (theStation == 3 and !runME31Up_)
      edm::LogError("CSCUpgradeMotherboard|SetupError") << "+++ TMB constructed while runME31Up_ is not set! +++\n";
    if (theStation == 4 and !runME41Up_)
      edm::LogError("CSCUpgradeMotherboard|SetupError") << "+++ TMB constructed while runME41Up_ is not set! +++\n";
  }

  // enable the upgrade processors
  if (runPhase2_ and theRing == 1) {
    clctProc = std::make_unique<CSCUpgradeCathodeLCTProcessor>(endcap, station, sector, subsector, chamber, conf);
    if (enableAlctPhase2_) {
      alctProc = std::make_unique<CSCUpgradeAnodeLCTProcessor>(endcap, station, sector, subsector, chamber, conf);
    }
  }

  match_earliest_alct_only = tmbParams_.getParameter<bool>("matchEarliestAlctOnly");
  match_earliest_clct_only = tmbParams_.getParameter<bool>("matchEarliestClctOnly");
  drop_used_clcts = tmbParams_.getParameter<bool>("tmbDropUsedClcts");
  tmb_cross_bx_sorting_ = tmbParams_.getParameter<bool>("tmbCrossBxSorting");
  debug_matching = tmbParams_.getParameter<bool>("debugMatching");
  // ignore unphysical ALCT-CLCT matches
  ignoreAlctCrossClct = tmbParams_.getParameter<bool>("ignoreAlctCrossClct");

  const edm::ParameterSet me11luts(tmbParams_.getParameter<edm::ParameterSet>("wgCrossHsME11Params"));
  cscOverlap_ = std::make_unique<CSCOverlap>(endcap, station, theRing, gangedME1a_, me11luts);

  // set prefential index for cross-bx matching
  pref[0] = match_trig_window_size / 2;
  for (unsigned int m = 2; m < match_trig_window_size; m += 2) {
    pref[m - 1] = pref[0] - m / 2;
    pref[m] = pref[0] + m / 2;
  }
}

void CSCUpgradeMotherboard::clear() {
  CSCMotherboard::clear();
  allLCTs.clear();
}

void CSCUpgradeMotherboard::run(const CSCWireDigiCollection* wiredc, const CSCComparatorDigiCollection* compdc) {
  clear();

  // Step 1: Setup
  if (!(alctProc and clctProc)) {
    edm::LogError("CSCUpgradeMotherboard|SetupError")
        << "+++ run() called for non-existing ALCT/CLCT processor! +++ \n";
    return;
  }

  alctProc->setCSCGeometry(cscGeometry_);
  clctProc->setCSCGeometry(cscGeometry_);

  alctV = alctProc->run(wiredc);  // run anodeLCT
  clctV = clctProc->run(compdc);  // run cathodeLCT

  // if there are no ALCTs and no CLCTs, it does not make sense to run this TMB
  if (alctV.empty() and clctV.empty())
    return;

  // encode high multiplicity bits
  encodeHighMultiplicityBits();

  int used_clct_mask[20];
  for (int c = 0; c < 20; ++c)
    used_clct_mask[c] = 0;

  // Step 2: ALCT-centric matching
  for (int bx_alct = 0; bx_alct < CSCConstants::MAX_ALCT_TBINS; bx_alct++) {
    if (alctProc->getBestALCT(bx_alct).isValid()) {
      const int bx_clct_start(bx_alct - match_trig_window_size / 2 - CSCConstants::ALCT_CLCT_OFFSET);
      const int bx_clct_stop(bx_alct + match_trig_window_size / 2 - CSCConstants::ALCT_CLCT_OFFSET);

      // ALCT-to-CLCT
      for (int bx_clct = bx_clct_start; bx_clct <= bx_clct_stop; bx_clct++) {
        if (bx_clct < 0 or bx_clct >= CSCConstants::MAX_CLCT_TBINS)
          continue;
        if (drop_used_clcts and used_clct_mask[bx_clct])
          continue;
        if (clctProc->getBestCLCT(bx_clct).isValid()) {
          if (debug_matching)
            LogTrace("CSCUpgradeMotherboard") << "++Valid ME21 CLCT: " << clctProc->getBestCLCT(bx_clct) << std::endl;

          int mbx = bx_clct - bx_clct_start;
          CSCUpgradeMotherboard::correlateLCTs(alctProc->getBestALCT(bx_alct),
                                               alctProc->getSecondALCT(bx_alct),
                                               clctProc->getBestCLCT(bx_clct),
                                               clctProc->getSecondCLCT(bx_clct),
                                               allLCTs(bx_alct, mbx, 0),
                                               allLCTs(bx_alct, mbx, 1));
          if (infoV > 1)
            LogTrace("CSCUpgradeMotherboard")
                << "Successful ALCT-CLCT match in ME21: bx_alct = " << bx_alct << "; match window: [" << bx_clct_start
                << "; " << bx_clct_stop << "]; bx_clct = " << bx_clct << std::endl;
          LogTrace("CSCUpgradeMotherboard") << "+++ Best CLCT Details: ";
          clctProc->getBestCLCT(bx_clct).print();
          LogTrace("CSCUpgradeMotherboard") << "+++ Second CLCT Details: ";
          clctProc->getSecondCLCT(bx_clct).print();
          if (allLCTs(bx_alct, mbx, 0).isValid()) {
            used_clct_mask[bx_clct] += 1;
            if (match_earliest_clct_only)
              break;
          }
        }
      }
    }
  }

  // Step 3: Counting and sorting of LCTs
  crossBxSorting();
}

void CSCUpgradeMotherboard::correlateLCTs(const CSCALCTDigi& bALCT,
                                          const CSCALCTDigi& sALCT,
                                          const CSCCLCTDigi& bCLCT,
                                          const CSCCLCTDigi& sCLCT,
                                          CSCCorrelatedLCTDigi& lct1,
                                          CSCCorrelatedLCTDigi& lct2) const {
  CSCALCTDigi bestALCT = bALCT;
  CSCALCTDigi secondALCT = sALCT;
  CSCCLCTDigi bestCLCT = bCLCT;
  CSCCLCTDigi secondCLCT = sCLCT;

  const bool anodeBestValid = bestALCT.isValid();
  const bool anodeSecondValid = secondALCT.isValid();
  const bool cathodeBestValid = bestCLCT.isValid();
  const bool cathodeSecondValid = secondCLCT.isValid();

  if (anodeBestValid and !anodeSecondValid)
    secondALCT = bestALCT;
  if (!anodeBestValid and anodeSecondValid)
    bestALCT = secondALCT;
  if (cathodeBestValid and !cathodeSecondValid)
    secondCLCT = bestCLCT;
  if (!cathodeBestValid and cathodeSecondValid)
    bestCLCT = secondCLCT;

  // ALCT-CLCT matching conditions are defined by "trig_enable" configuration
  // parameters.
  if ((alct_trig_enable and bestALCT.isValid()) or (clct_trig_enable and bestCLCT.isValid()) or
      (match_trig_enable and bestALCT.isValid() and bestCLCT.isValid())) {
    lct1 = constructLCTs(bestALCT, bestCLCT, CSCCorrelatedLCTDigi::ALCTCLCT, 1);
  }

  if (((secondALCT != bestALCT) or (secondCLCT != bestCLCT)) and
      ((alct_trig_enable and secondALCT.isValid()) or (clct_trig_enable and secondCLCT.isValid()) or
       (match_trig_enable and secondALCT.isValid() and secondCLCT.isValid()))) {
    lct2 = constructLCTs(secondALCT, secondCLCT, CSCCorrelatedLCTDigi::ALCTCLCT, 2);
  }
}

bool CSCUpgradeMotherboard::doesALCTCrossCLCT(const CSCALCTDigi& a, const CSCCLCTDigi& c) const {
  return cscOverlap_->doesALCTCrossCLCT(a, c);
}

void CSCUpgradeMotherboard::correlateLCTsME11(const CSCALCTDigi& bALCT,
                                              const CSCALCTDigi& sALCT,
                                              const CSCCLCTDigi& bCLCT,
                                              const CSCCLCTDigi& sCLCT,
                                              CSCCorrelatedLCTDigi& lct1,
                                              CSCCorrelatedLCTDigi& lct2) const {
  // assume that always anodeBestValid && cathodeBestValid
  CSCALCTDigi bestALCT = bALCT;
  CSCALCTDigi secondALCT = sALCT;
  CSCCLCTDigi bestCLCT = bCLCT;
  CSCCLCTDigi secondCLCT = sCLCT;

  if (ignoreAlctCrossClct) {
    const bool anodeBestValid = bestALCT.isValid();
    const bool anodeSecondValid = secondALCT.isValid();
    const bool cathodeBestValid = bestCLCT.isValid();
    const bool cathodeSecondValid = secondCLCT.isValid();
    if (anodeBestValid and !anodeSecondValid)
      secondALCT = bestALCT;
    if (!anodeBestValid and anodeSecondValid)
      bestALCT = secondALCT;
    if (cathodeBestValid and !cathodeSecondValid)
      secondCLCT = bestCLCT;
    if (!cathodeBestValid and cathodeSecondValid)
      bestCLCT = secondCLCT;
    // ALCT-CLCT matching conditions are defined by "trig_enable" configuration
    // parameters.
    if ((alct_trig_enable and bestALCT.isValid()) or (clct_trig_enable and bestCLCT.isValid()) or
        (match_trig_enable and bestALCT.isValid() and bestCLCT.isValid())) {
      lct1 = constructLCTs(bestALCT, bestCLCT, CSCCorrelatedLCTDigi::ALCTCLCT, 1);
    }
    if (((secondALCT != bestALCT) or (secondCLCT != bestCLCT)) and
        ((alct_trig_enable and secondALCT.isValid()) or (clct_trig_enable and secondCLCT.isValid()) or
         (match_trig_enable and secondALCT.isValid() and secondCLCT.isValid()))) {
      lct2 = constructLCTs(secondALCT, secondCLCT, CSCCorrelatedLCTDigi::ALCTCLCT, 2);
    }
    return;
  } else {
    if (secondALCT == bestALCT)
      secondALCT.clear();
    if (secondCLCT == bestCLCT)
      secondCLCT.clear();

    const int ok11 = doesALCTCrossCLCT(bestALCT, bestCLCT);
    const int ok12 = doesALCTCrossCLCT(bestALCT, secondCLCT);
    const int ok21 = doesALCTCrossCLCT(secondALCT, bestCLCT);
    const int ok22 = doesALCTCrossCLCT(secondALCT, secondCLCT);
    const int code = (ok11 << 3) | (ok12 << 2) | (ok21 << 1) | (ok22);

    int dbg = 0;
    if (dbg)
      LogTrace("CSCMotherboardME11") << "debug correlateLCTs in ME11 " << cscId_ << std::endl
                                     << "ALCT1: " << bestALCT << std::endl
                                     << "ALCT2: " << secondALCT << std::endl
                                     << "CLCT1: " << bestCLCT << std::endl
                                     << "CLCT2: " << secondCLCT << std::endl
                                     << "ok 11 12 21 22 code = " << ok11 << " " << ok12 << " " << ok21 << " " << ok22
                                     << " " << code << std::endl;

    if (code == 0)
      return;

    // LUT defines correspondence between possible ok## combinations
    // and resulting lct1 and lct2
    int lut[16][2] = {
        //ok: 11 12 21 22
        {0, 0},    // 0  0  0  0
        {22, 0},   // 0  0  0  1
        {21, 0},   // 0  0  1  0
        {21, 22},  // 0  0  1  1
        {12, 0},   // 0  1  0  0
        {12, 22},  // 0  1  0  1
        {12, 21},  // 0  1  1  0
        {12, 21},  // 0  1  1  1
        {11, 0},   // 1  0  0  0
        {11, 22},  // 1  0  0  1
        {11, 21},  // 1  0  1  0
        {11, 22},  // 1  0  1  1
        {11, 12},  // 1  1  0  0
        {11, 22},  // 1  1  0  1
        {11, 12},  // 1  1  1  0
        {11, 22},  // 1  1  1  1
    };

    if (dbg)
      LogTrace("CSCMotherboardME11") << "lut 0 1 = " << lut[code][0] << " " << lut[code][1] << std::endl;

    switch (lut[code][0]) {
      case 11:
        lct1 = constructLCTs(bestALCT, bestCLCT, CSCCorrelatedLCTDigi::ALCTCLCT, 1);
        break;
      case 12:
        lct1 = constructLCTs(bestALCT, secondCLCT, CSCCorrelatedLCTDigi::ALCTCLCT, 1);
        break;
      case 21:
        lct1 = constructLCTs(secondALCT, bestCLCT, CSCCorrelatedLCTDigi::ALCTCLCT, 1);
        break;
      case 22:
        lct1 = constructLCTs(secondALCT, secondCLCT, CSCCorrelatedLCTDigi::ALCTCLCT, 1);
        break;
      default:
        return;
    }

    if (dbg)
      LogTrace("CSCMotherboardME11") << "lct1: " << lct1 << std::endl;

    switch (lut[code][1]) {
      case 12:
        lct2 = constructLCTs(bestALCT, secondCLCT, CSCCorrelatedLCTDigi::ALCTCLCT, 2);
        if (dbg)
          LogTrace("CSCMotherboardME11") << "lct2: " << lct2 << std::endl;
        return;
      case 21:
        lct2 = constructLCTs(secondALCT, bestCLCT, CSCCorrelatedLCTDigi::ALCTCLCT, 2);
        if (dbg)
          LogTrace("CSCMotherboardME11") << "lct2: " << lct2 << std::endl;
        return;
      case 22:
        lct2 = constructLCTs(secondALCT, secondCLCT, CSCCorrelatedLCTDigi::ALCTCLCT, 2);
        if (dbg)
          LogTrace("CSCMotherboardME11") << "lct2: " << lct2 << std::endl;
        return;
      default:
        return;
    }
    if (dbg)
      LogTrace("CSCMotherboardME11") << "out of correlateLCTsME11" << std::endl;

    return;
  }
}

void CSCUpgradeMotherboard::crossBxSorting() {
  // reduction of nLCTs per each BX
  for (int bx = 0; bx < CSCConstants::MAX_LCT_TBINS; bx++) {
    // counting across the trigger window
    unsigned int n = 0;
    for (unsigned int mbx = 0; mbx < match_trig_window_size; mbx++) {
      for (int i = 0; i < CSCConstants::MAX_LCTS_PER_CSC; i++) {
        if (allLCTs(bx, mbx, i).isValid()) {
          ++n;
          if (infoV > 0) {
            LogDebug("CSCUpgradeMotherboard")
                << "Candidate LCT" << i + 1 << " " << bx << "/" << bx + mbx - match_trig_window_size / 2 << ": "
                << allLCTs(bx, mbx, i) << std::endl;
          }
        }
      }
    }

    // some simple cross-bx sorting algorithms
    if (tmb_cross_bx_sorting_ and (n > CSCConstants::MAX_LCTS_PER_CSC)) {
      n = 0;
      for (unsigned int mbx = 0; mbx < match_trig_window_size; mbx++) {
        for (int i = 0; i < CSCConstants::MAX_LCTS_PER_CSC; i++) {
          if (allLCTs(bx, pref[mbx], i).isValid()) {
            n++;
            if (n > CSCConstants::MAX_LCTS_PER_CSC)
              allLCTs(bx, pref[mbx], i).clear();
          }
        }
      }

      // show the final LCTs
      n = 0;
      for (unsigned int mbx = 0; mbx < match_trig_window_size; mbx++) {
        for (int i = 0; i < CSCConstants::MAX_LCTS_PER_CSC; i++) {
          if (allLCTs(bx, mbx, i).isValid()) {
            n++;
            if (infoV > 0) {
              LogDebug("CSCUpgradeMotherboard")
                  << "Selected LCT" << i + 1 << " " << bx << "/" << bx + mbx - match_trig_window_size / 2 << ": "
                  << allLCTs(bx, mbx, i) << std::endl;
            }
          }
        }
      }
      if (infoV > 0 and n > 0)
        LogDebug("CSCUpgradeMotherboard") << "bx " << bx << " nnLCT:" << n << " " << n << std::endl;
    }  // x-bx sorting
  }
}

std::vector<CSCCorrelatedLCTDigi> CSCUpgradeMotherboard::readoutLCTs() const {
  // get the time matched LCTs (should be done after sorting and selecting)
  std::vector<CSCCorrelatedLCTDigi> result;
  allLCTs.getMatched(result);

  // do a final check on the LCTs in readout
  qualityControl_->checkMultiplicityBX(result);
  for (const auto& lct : result) {
    qualityControl_->checkValid(lct);
  }
  return result;
}
