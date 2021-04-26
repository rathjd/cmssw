import FWCore.ParameterSet.Config as cms

tmbPhase1 = cms.PSet(
    mpcBlockMe1a    = cms.uint32(0),
    alctTrigEnable  = cms.uint32(0),
    clctTrigEnable  = cms.uint32(0),
    matchTrigEnable = cms.uint32(1),
    matchTrigWindowSize = cms.uint32(7),
    tmbL1aWindowSize = cms.uint32(7),
    # Debug
    verbosity = cms.int32(0),

    # Configure early_tbins instead of hardcoding it
    tmbEarlyTbins = cms.int32(4),

    # Flag for whether to readout only the earliest max two LCTs in a
    # L1A readout window, as there is only room just for two in the TMB header.
    # If false, all LCTs would be readout in L1A window.
    tmbReadoutEarliest2 = cms.bool(True),

    # For CLCT-centric matching, whether to drop ALCTs that were matched
    # to CLCTs in this BX, and not use them in the following BX
    tmbDropUsedAlcts = cms.bool(True),

    # For ALCT-centric matching, whether to drop CLCTs that were matched
    # to ALCTs in this BX, and not use them in the following BX
    tmbDropUsedClcts = cms.bool(False),

    ## bits for high-multiplicity triggers
    useHighMultiplicityBits = cms.bool(False),
)

# to be used by ME11 chambers with upgraded TMB and ALCT
tmbPhase2 = tmbPhase1.clone(
    # reduce ALCT-CLCT matching window size from 7 to 3
    matchTrigWindowSize = 3,
    tmbReadoutEarliest2 = False,
    tmbDropUsedAlcts = False,
    tmbDropUsedClcts = False,

    # For CLCT-centric matching in ME11, break after finding
    # the first BX with matching ALCT
    matchEarliestAlctOnly = cms.bool(False),

    # For ALCT-centric matching in ME11, break after finding
    # the first BX with matching CLCT
    matchEarliestClctOnly = cms.bool(False),

    # False = default "non-X-BX" sorting algorithm,
    #     where the first BX with match goes first
    # True = simple X-BX sorting algorithm,
    #     where the central match BX goes first,
    #     then the closest early, the closest late, etc.
    tmbCrossBxSorting = cms.bool(True),

    # True: allow construction of unphysical LCTs
    # in ME11 for which WG and HS do not intersect
    # False: do not build unphysical LCTs
    ignoreAlctCrossClct = cms.bool(True),

    ## run in debug mode
    debugLUTs = cms.bool(False),
    debugMatching = cms.bool(False),
)

# to be used by ME11 chambers with GEM-CSC ILT
tmbPhase2GE11 = tmbPhase2.clone(
    ## matching to pads
    maxDeltaBXPad = cms.int32(1),
    maxDeltaBXCoPad = cms.int32(1),
    maxDeltaPadL1Even = cms.int32(12),
    maxDeltaPadL1Odd = cms.int32(24),
    maxDeltaPadL2Even = cms.int32(12),
    maxDeltaPadL2Odd = cms.int32(24),

    ## efficiency recovery switches
    dropLowQualityCLCTsNoGEMs_ME1a = cms.bool(False),
    dropLowQualityCLCTsNoGEMs_ME1b = cms.bool(True),
    dropLowQualityALCTsNoGEMs_ME1a = cms.bool(False),
    dropLowQualityALCTsNoGEMs_ME1b = cms.bool(False),
    buildLCTfromALCTandGEM_ME1a = cms.bool(False),
    buildLCTfromALCTandGEM_ME1b = cms.bool(True),
    buildLCTfromCLCTandGEM_ME1a = cms.bool(False),
    buildLCTfromCLCTandGEM_ME1b = cms.bool(True),
    promoteALCTGEMpattern = cms.bool(True),
    promoteALCTGEMquality = cms.bool(True),
    promoteCLCTGEMquality_ME1a = cms.bool(True),
    promoteCLCTGEMquality_ME1b = cms.bool(True),
)

# to be used by ME21 chambers with GEM-CSC ILT
tmbPhase2GE21 = tmbPhase2.clone(
    ## matching to pads
    maxDeltaBXPad = cms.int32(1),
    maxDeltaBXCoPad = cms.int32(1),
    maxDeltaPadL1Even = cms.int32(12),
    maxDeltaPadL1Odd = cms.int32(24),
    maxDeltaPadL2Even = cms.int32(12),
    maxDeltaPadL2Odd = cms.int32(24),

    ## efficiency recovery switches
    dropLowQualityALCTsNoGEMs = cms.bool(True),
    dropLowQualityCLCTsNoGEMs = cms.bool(True),
    buildLCTfromALCTandGEM = cms.bool(True),
    buildLCTfromCLCTandGEM = cms.bool(True),
    promoteALCTGEMpattern = cms.bool(True),
    promoteALCTGEMquality = cms.bool(True),
    promoteCLCTGEMquality = cms.bool(True),
)

## LUTs to map wiregroup onto min and max half-strip number that it crosses in ME1/1
wgCrossHsME11Params = cms.PSet(
    wgCrossHsME1aFiles = cms.vstring(
        "L1Trigger/CSCTriggerPrimitives/data/ME11/CSCLUT_wg_min_hs_ME1a.txt",
        "L1Trigger/CSCTriggerPrimitives/data/ME11/CSCLUT_wg_max_hs_ME1a.txt",
    ),
    wgCrossHsME1aGangedFiles = cms.vstring(
        "L1Trigger/CSCTriggerPrimitives/data/ME11/CSCLUT_wg_min_hs_ME1a_ganged.txt",
        "L1Trigger/CSCTriggerPrimitives/data/ME11/CSCLUT_wg_max_hs_ME1a_ganged.txt",
    ),
    wgCrossHsME1bFiles = cms.vstring(
        "L1Trigger/CSCTriggerPrimitives/data/ME11/CSCLUT_wg_min_hs_ME1b.txt",
        "L1Trigger/CSCTriggerPrimitives/data/ME11/CSCLUT_wg_max_hs_ME1b.txt",
    )
)

tmbPSets = cms.PSet(
    tmbPhase1 = tmbPhase1.clone(),
    tmbPhase2 = tmbPhase2.clone(),
    wgCrossHsME11Params = wgCrossHsME11Params.clone(),
    tmbPhase2GE11 = tmbPhase2GE11.clone(),
    tmbPhase2GE21 = tmbPhase2GE21.clone(),
)
