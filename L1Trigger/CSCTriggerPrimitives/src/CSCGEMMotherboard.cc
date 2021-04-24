#include <memory>

#include "L1Trigger/CSCTriggerPrimitives/interface/CSCGEMMotherboard.h"

CSCGEMMotherboard::CSCGEMMotherboard(unsigned endcap,
                                     unsigned station,
                                     unsigned sector,
                                     unsigned subsector,
                                     unsigned chamber,
                                     const edm::ParameterSet& conf)
    : CSCUpgradeMotherboard(endcap, station, sector, subsector, chamber, conf),

  promoteALCTGEMpattern_(tmbParams_.getParameter<bool>("promoteALCTGEMpattern")),

  promoteALCTGEMquality_(tmbParams_.getParameter<bool>("promoteALCTGEMquality")),
  promoteCLCTGEMquality_(tmbParams_.getParameter<bool>("promoteCLCTGEMquality")),
  promoteCLCTGEMquality_ME1a_(tmbParams_.getParameter<bool>("promoteCLCTGEMquality_ME1a")),
  promoteCLCTGEMquality_ME1b_(tmbParams_.getParameter<bool>("promoteCLCTGEMquality_ME1b")),

      dropLowQualityCLCTsNoGEMs_ME1a_(tmbParams_.getParameter<bool>("dropLowQualityCLCTsNoGEMs_ME1a")),
      dropLowQualityCLCTsNoGEMs_ME1b_(tmbParams_.getParameter<bool>("dropLowQualityCLCTsNoGEMs_ME1b")),
  dropLowQualityALCTsNoGEMs_ME1a_(tmbParams_.getParameter<bool>("dropLowQualityALCTsNoGEMs_ME1a")),
  dropLowQualityALCTsNoGEMs_ME1b_(tmbParams_.getParameter<bool>("dropLowQualityALCTsNoGEMs_ME1b")),
  dropLowQualityCLCTsNoGEMs_(tmbParams_.getParameter<bool>("dropLowQualityCLCTsNoGEMs")),
  dropLowQualityALCTsNoGEMs_(tmbParams_.getParameter<bool>("dropLowQualityALCTsNoGEMs")),

  buildLCTfromALCTandGEM_ME1a_(tmbParams_.getParameter<bool>("buildLCTfromALCTandGEM_ME1a")),
  buildLCTfromALCTandGEM_ME1b_(tmbParams_.getParameter<bool>("buildLCTfromALCTandGEM_ME1b")),
  buildLCTfromCLCTandGEM_ME1a_(tmbParams_.getParameter<bool>("buildLCTfromCLCTandGEM_ME1a")),
  buildLCTfromCLCTandGEM_ME1b_(tmbParams_.getParameter<bool>("buildLCTfromCLCTandGEM_ME1b")),
  buildLCTfromALCTandGEM_(tmbParams_.getParameter<bool>("buildLCTfromALCTandGEM")),
  buildLCTfromCLCTandGEM_(tmbParams_.getParameter<bool>("buildLCTfromCLCTandGEM"))
{

  if (!runPhase2_) {
    edm::LogError("CSCGEMMotherboard|SetupError") << "+++ TMB constructed while runPhase2 is not set! +++\n";
  }

  if (!runME11ILT_) {
    edm::LogError("CSCGEMMotherboard|SetupError") << "+++ TMB constructed while runME11ILT_ is not set! +++\n";
  };

  // super chamber has layer=0!
  gemId = GEMDetId(theRegion, 1, theStation, 0, theChamber, 0).rawId();

  const edm::ParameterSet coPadParams(station == 1 ? conf.getParameter<edm::ParameterSet>("copadParamGE11")
                                                   : conf.getParameter<edm::ParameterSet>("copadParamGE21"));

  const edm::ParameterSet gemcscluts(conf.getParameter<edm::ParameterSet>("gemcscParams"));
  coPadProcessor = std::make_unique<GEMCoPadProcessor>(theRegion, theStation, theChamber, coPadParams, gemcscluts);

  cscGEMMatcher_ = std::make_unique<CSCGEMMatcher>(theRegion, theStation, theChamber, gemcscluts);
}

CSCGEMMotherboard::~CSCGEMMotherboard() {}

void CSCGEMMotherboard::clear() {
  CSCUpgradeMotherboard::clear();
  coPadProcessor->clear();
}

//===============================================
//use ALCTs, CLCTs, GEMs to build LCTs
//loop over each BX to find valid ALCT in ME1b, try ALCT-CLCT match, if ALCT-CLCT match failed, try ALCT-GEM match
//do the same in ME1a
//sort LCTs according to different algorithm, and send out best 2 LCTs
//===============================================

void CSCGEMMotherboard::run(const CSCWireDigiCollection* wiredc,
                                const CSCComparatorDigiCollection* compdc,
                                const GEMPadDigiClusterCollection* gemClusters) {
  // Step 1: Setup
  clear();

  // encode high multiplicity bits
  encodeHighMultiplicityBits();

  // check for GEM geometry
  if (gem_g != nullptr) {
    edm::LogError("CSCGEMMotherboard|SetupError")
        << "+++ run() called for GEM-CSC integrated trigger without valid GEM geometry! +++ \n";
    return;
  }

  if (!(alctProc and clctProc)) {
    edm::LogError("CSCGEMMotherboard|SetupError")
        << "+++ run() called for non-existing ALCT/CLCT processor! +++ \n";
    return;
  }

  alctProc->setCSCGeometry(cscGeometry_);
  clctProc->setCSCGeometry(cscGeometry_);

  alctV = alctProc->run(wiredc);  // run anodeLCT
  clctV = clctProc->run(compdc);  // run cathodeLCT

  processGEMClusters(gemClusters);

  // if there are no ALCTs and no CLCTs, it does not make sense to run this TMB
  if (alctV.empty() and clctV.empty())
    return;

  int used_clct_mask[20];
  for (int b = 0; b < 20; b++)
    used_clct_mask[b] = 0;

  // Step 2: ALCT-centric matching
  for (int bx_alct = 0; bx_alct < CSCConstants::MAX_ALCT_TBINS; bx_alct++) {
    if (alctProc->getBestALCT(bx_alct).isValid()) {
      const int bx_clct_start(bx_alct - match_trig_window_size / 2 - CSCConstants::ALCT_CLCT_OFFSET);
      const int bx_clct_stop(bx_alct + match_trig_window_size / 2 - CSCConstants::ALCT_CLCT_OFFSET);
      const int bx_copad_start(bx_alct - 3);
      const int bx_copad_stop(bx_alct + 3);

      if (debug_matching) {
        LogTrace("CSCGEMMotherboard")
            << "========================================================================\n"
            << "ALCT-CLCT matching in ME1/1 chamber: " << cscId_ << "\n"
            << "------------------------------------------------------------------------\n"
            << "+++ Best ALCT Details: " << alctProc->getBestALCT(bx_alct) << "\n"
            << "+++ Second ALCT Details: " << alctProc->getSecondALCT(bx_alct) << std::endl;

        LogTrace("CSCGEMMotherboard")
            << "------------------------------------------------------------------------ \n"
            << "Attempt ALCT-CLCT matching in ME1/b in bx range: [" << bx_clct_start << "," << bx_clct_stop << "]"
            << std::endl;
      }

      // ALCT-to-CLCT matching in ME1b
      int nSuccesFulMatches = 0;
      for (int bx_clct = bx_clct_start; bx_clct <= bx_clct_stop; bx_clct++) {
        if (bx_clct < 0 or bx_clct >= CSCConstants::MAX_CLCT_TBINS)
          continue;
        if (drop_used_clcts and used_clct_mask[bx_clct])
          continue;
        if (clctProc->getBestCLCT(bx_clct).isValid()) {
          if (debug_matching)
            LogTrace("CSCGEMMotherboard") << "++Valid ME1b CLCT: " << clctProc->getBestCLCT(bx_clct) << std::endl;

          int mbx = bx_clct - bx_clct_start;
          correlateLCTsGEM(alctProc->getBestALCT(bx_alct),
                           alctProc->getSecondALCT(bx_alct),
                           clctProc->getBestCLCT(bx_clct),
                           clctProc->getSecondCLCT(bx_clct),
                           clusters_[bx_alct],
                           allLCTs(bx_alct, mbx, 0),
                           allLCTs(bx_alct, mbx, 1));

          if (allLCTs(bx_alct, mbx, 0).isValid()) {
            ++nSuccesFulMatches;
            used_clct_mask[bx_clct] += 1;

            if (match_earliest_clct_only)
              break;
          } else
            LogTrace("CSCGEMMotherboard") << "No valid LCT is built from ALCT-CLCT matching in ME1b" << std::endl;
        }
      }

      // ALCT-to-GEM matching in ME1b
      int nSuccesFulGEMMatches = 0;
      if (nSuccesFulMatches == 0) {
        if (debug_matching)
          LogTrace("CSCGEMMotherboard") << "++No valid ALCT-CLCT matches in ME1b" << std::endl;
        for (int bx_gem = bx_copad_start; bx_gem <= bx_copad_stop; bx_gem++) {

          correlateLCTsGEM(alctProc->getBestALCT(bx_alct),
                           alctProc->getSecondALCT(bx_alct),
                           clusters_[bx_alct],
                           allLCTs(bx_alct, 0, 0),
                           allLCTs(bx_alct, 0, 1));

          if (allLCTs(bx_alct, 0, 0).isValid()) {
            ++nSuccesFulGEMMatches;

            if (match_earliest_clct_only)
              break;
          }
        }
      }
    }  // end of ALCT valid block
    else {
      if (true) {
        // keep it simple for the time being, only consider the first copad
        const int bx_clct_start(bx_alct - match_trig_window_size / 2 - CSCConstants::ALCT_CLCT_OFFSET);
        const int bx_clct_stop(bx_alct + match_trig_window_size / 2 - CSCConstants::ALCT_CLCT_OFFSET);

        for (int bx_clct = bx_clct_start; bx_clct <= bx_clct_stop; bx_clct++) {
          if (bx_clct < 0 or bx_clct >= CSCConstants::MAX_CLCT_TBINS)
            continue;
          if (drop_used_clcts and used_clct_mask[bx_clct])
            continue;
          if (clctProc->getBestCLCT(bx_clct).isValid()) {
            const int quality(clctProc->getBestCLCT(bx_clct).getQuality());

            int mbx = bx_clct - bx_clct_start;

            correlateLCTsGEM(clctProc->getBestCLCT(bx_clct),
                             clctProc->getSecondCLCT(bx_clct),
                             clusters_[bx_clct],
                             allLCTs(bx_alct, mbx, 0),
                             allLCTs(bx_alct, mbx, 1));

            if (allLCTs(bx_alct, mbx, 0).isValid()) {
              used_clct_mask[bx_clct] += 1;

              if (match_earliest_clct_only)
                break;
            }
          }
        }  //end of clct loop
      }
    }
  }  // end of ALCT-centric matching

  // Step 3: Counting and sorting of LCTs
  CSCUpgradeMotherboard::crossBxSorting();
}

void CSCGEMMotherboard::processGEMClusters(const GEMPadDigiClusterCollection* gemClusters) {
  auto tempClusters = coPadProcessor->run(gemClusters);

  // sort clusters by bx
  for (const auto& cl : tempClusters) {
    // add the central BX for convenience
    clusters_[cl.bx()].push_back(cl);
  }
}

void CSCGEMMotherboard::correlateLCTsGEM(const CSCALCTDigi& bestALCT,
                                         const CSCALCTDigi& secondALCT,
                                         const GEMInternalClusters& clusters,
                                         CSCCorrelatedLCTDigi& lct1,
                                         CSCCorrelatedLCTDigi& lct2) const
{
  // output candidates
  std::vector<CSCCorrelatedLCTDigi> lcts;

  const auto& bestClusters = cscGEMMatcher_->matchingClustersBXLoc(bestALCT, clusters);
  const auto& secondClusters = cscGEMMatcher_->matchingClustersBXLoc(secondALCT, clusters);

  // check all possibilities
  // push all to the output
  for (const auto& cl : bestClusters) {
    lcts.push_back(constructLCTsGEM(bestALCT, CSCCLCTDigi(), cl, 0));
  }

  for (const auto& cl : secondClusters) {
    lcts.push_back(constructLCTsGEM(secondALCT, CSCCLCTDigi(), cl, 0));
  }

  // retain first 2

  // assign track number
}

void CSCGEMMotherboard::correlateLCTsGEM(const CSCCLCTDigi& bestCLCT,
                                         const CSCCLCTDigi& secondCLCT,
                                         const GEMInternalClusters& clusters,
                                         CSCCorrelatedLCTDigi& lct1,
                                         CSCCorrelatedLCTDigi& lct2) const
{
  // output candidates
  std::vector<CSCCorrelatedLCTDigi> lcts;

  const auto& bestClusters = cscGEMMatcher_->matchingClustersBXLoc(bestCLCT, clusters);
  const auto& secondClusters = cscGEMMatcher_->matchingClustersBXLoc(secondCLCT, clusters);

  // check all possibilities
  // push all to the output
  for (const auto& cl : bestClusters) {
    lcts.push_back(constructLCTsGEM(CSCALCTDigi(), bestCLCT, cl, 0));
  }

  for (const auto& cl : secondClusters) {
    lcts.push_back(constructLCTsGEM(CSCALCTDigi(), secondCLCT, cl, 0));
  }

  // sort by bending angle

  // retain 2 with smallest bending angle

  // assign track number
}

void CSCGEMMotherboard::correlateLCTsGEM(const CSCALCTDigi& bestALCT,
                                         const CSCALCTDigi& secondALCT,
                                         const CSCCLCTDigi& bestCLCT,
                                         const CSCCLCTDigi& secondCLCT,
                                         const GEMInternalClusters& clusters,
                                         CSCCorrelatedLCTDigi& lct1,
                                         CSCCorrelatedLCTDigi& lct2) const
{
  // output candidates
  std::vector<CSCCorrelatedLCTDigi> lcts;

  const auto& bbClusters = cscGEMMatcher_->matchingClustersBXLoc(bestALCT, bestCLCT, clusters);
  const auto& bsClusters = cscGEMMatcher_->matchingClustersBXLoc(bestALCT, secondCLCT, clusters);
  const auto& sbClusters = cscGEMMatcher_->matchingClustersBXLoc(secondALCT, bestCLCT, clusters);
  const auto& ssClusters = cscGEMMatcher_->matchingClustersBXLoc(secondALCT, secondCLCT, clusters);

  // check all possibilities
  // push all to the output
  for (const auto& cl : bbClusters) {
    lcts.push_back(constructLCTsGEM(bestALCT, bestCLCT, cl, 0));
  }

  for (const auto& cl : bsClusters) {
    lcts.push_back(constructLCTsGEM(bestALCT, secondCLCT, cl, 0));
  }

  for (const auto& cl : sbClusters) {
    lcts.push_back(constructLCTsGEM(secondALCT, bestCLCT, cl, 0));
  }

  for (const auto& cl : ssClusters) {
    lcts.push_back(constructLCTsGEM(secondALCT, secondCLCT, cl, 0));
  }

  // sort by bending angle

  // retain 2 with smallest bending angle

  // assign track number
}

CSCCorrelatedLCTDigi CSCGEMMotherboard::constructLCTsGEM(const CSCALCTDigi& alct,
                                                         const CSCCLCTDigi& clct,
                                                         const GEMInternalCluster& gem1,
                                                         int trknmb) const {
  int pattern = 0, quality = 0, bx = 0, keyStrip = 0, keyWG = 0, bend = 0, valid = 0;

  // make a new LCT
  CSCCorrelatedLCTDigi thisLCT;

  if (!alct.isValid() and !clct.isValid()) {
    edm::LogError("CSCGEMCMotherboard") << "Warning!!! neither ALCT nor CLCT valid, return invalid LCT";
    return thisLCT;
  }


  /*

    // Determine the case and assign properties depending on the LCT dataformat (old/new)
    if (alct.isValid() and clct.isValid() and gem1.isValid() and not gem2.isValid()) {
    pattern = encodePattern(clct.getPattern());
    if (runCCLUT_) {
    quality = qualityAssignment_->findQualityGEMv2(alct, clct, 1);
    } else {
    quality = qualityAssignment_->findQualityGEMv1(alct, clct, 1);
    }
    bx = alct.getBX();
    keyStrip = clct.getKeyStrip();
    keyWG = alct.getKeyWG();
    bend = clct.getBend();
    thisLCT.setALCT(getBXShiftedALCT(alct));
    thisLCT.setCLCT(getBXShiftedCLCT(clct));
    thisLCT.setGEM1(gem1);
    thisLCT.setType(CSCCorrelatedLCTDigi::ALCTCLCTGEM);
    valid = doesALCTCrossCLCT(alct, clct);
    if (runCCLUT_) {
    thisLCT.setRun3(true);
    // 4-bit slope value derived with the CCLUT algorithm
    thisLCT.setSlope(clct.getSlope());
    thisLCT.setQuartStrip(clct.getQuartStrip());
    thisLCT.setEighthStrip(clct.getEighthStrip());
    thisLCT.setRun3Pattern(clct.getRun3Pattern());
    }
    } else if (alct.isValid() and clct.isValid() and not gem1.isValid() and gem2.isValid()) {
    pattern = encodePattern(clct.getPattern());
    if (runCCLUT_) {
    quality = qualityAssignment_->findQualityGEMv2(alct, clct, 2);
    } else {
    quality = qualityAssignment_->findQualityGEMv1(alct, clct, 2);
    }
    bx = alct.getBX();
    keyStrip = clct.getKeyStrip();
    keyWG = alct.getKeyWG();
    bend = clct.getBend();
    thisLCT.setALCT(getBXShiftedALCT(alct));
    thisLCT.setCLCT(getBXShiftedCLCT(clct));
    thisLCT.setGEM1(gem2.first());
    thisLCT.setGEM2(gem2.second());
    thisLCT.setType(CSCCorrelatedLCTDigi::ALCTCLCT2GEM);
    valid = doesALCTCrossCLCT(alct, clct);
    if (runCCLUT_) {
    thisLCT.setRun3(true);
    // 4-bit slope value derived with the CCLUT algorithm
    thisLCT.setSlope(clct.getSlope());
    thisLCT.setQuartStrip(clct.getQuartStrip());
    thisLCT.setEighthStrip(clct.getEighthStrip());
    thisLCT.setRun3Pattern(clct.getRun3Pattern());
    }
    } else if (alct.isValid() and gem2.isValid() and not clct.isValid()) {
    //in ME11
    //ME1b: keyWG >15,
    //ME1a and ME1b overlap:  10<=keyWG<=15
    //ME1a: keyWG < 10
    //in overlap region, firstly try a match in ME1b

    keyWG = alct.getKeyWG();

    pattern = promoteALCTGEMpattern_ ? 10 : 0;
    quality = promoteALCTGEMquality_ ? 15 : 11;
    bx = alct.getBX();
    thisLCT.setALCT(getBXShiftedALCT(alct));
    thisLCT.setGEM1(gem2.first());
    thisLCT.setGEM2(gem2.second());
    thisLCT.setType(CSCCorrelatedLCTDigi::ALCT2GEM);
    valid = true;
    } else if (clct.isValid() and gem2.isValid() and not alct.isValid()) {
    // min roll number is always 1
    // max roll number is 8 or 16, depending on the station

    pattern = encodePattern(clct.getPattern());
    quality = promoteCLCTGEMquality_ ? 15 : 11;
    bx = gem2.bx(1) + CSCConstants::LCT_CENTRAL_BX;
    keyStrip = clct.getKeyStrip();
    // choose the corresponding wire-group in the middle of the partition
    // keyWG = mymap2.at(gem2.roll() - 1);
    bend = clct.getBend();
    thisLCT.setCLCT(clct);
    thisLCT.setGEM1(gem2.first());
    thisLCT.setGEM2(gem2.second());
    thisLCT.setType(CSCCorrelatedLCTDigi::CLCT2GEM);
    valid = true;
    if (runCCLUT_) {
    thisLCT.setRun3(true);
    // 4-bit slope value derived with the CCLUT algorithm
    thisLCT.setSlope(clct.getSlope());
    thisLCT.setQuartStrip(clct.getQuartStrip());
    thisLCT.setEighthStrip(clct.getEighthStrip());
    thisLCT.setRun3Pattern(clct.getRun3Pattern());
    }
    }

    if (valid == 0)
    LogTrace("CSCGEMCMotherboard") << "Warning!!! wiregroup and strip pair are not crossing each other"
    << " detid " << cscId_ << " with wiregroup " << keyWG << "keyStrip " << keyStrip
    << " \n";

    // fill the rest of the properties
    thisLCT.setTrknmb(trknmb);
    thisLCT.setValid(valid);
    thisLCT.setQuality(quality);
    thisLCT.setWireGroup(keyWG);
    thisLCT.setStrip(keyStrip);
    thisLCT.setPattern(pattern);
    thisLCT.setBend(bend);
    thisLCT.setBX(bx);
    thisLCT.setMPCLink(0);
    thisLCT.setBX0(0);
    // Not used in Run-2. Will not be assigned in Run-3
    thisLCT.setSyncErr(0);
    thisLCT.setCSCID(theTrigChamber);

    // future work: add a section that produces LCTs according
    // to the new LCT dataformat (not yet defined)
    */

  // return new LCT
  return thisLCT;
}
