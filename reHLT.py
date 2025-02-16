# hltGetConfiguration /users/sharper/2022/scouting/FullAndPixelScouting/V7 --globaltag auto:phase1_2021_realistic --mc --unprescale --eras Run3 --output minimal --max-events 100 --l1-emulator uGT --l1 L1Menu_Collisions2022_v1_0_0_xml --input /store/mc/Run3Summer21DRPremix/TT_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW/120X_mcRun3_2021_realistic_v6-v2/2540000/b354245e-d8bc-424d-b527-58815586a6a5.root --customise HLTrigger/Configuration/customizeHLTforPatatrack.customizeHLTforPatatrackTriplets

# /users/sharper/2022/scouting/FullAndPixelScouting/V7 (CMSSW_12_3_0)

import FWCore.ParameterSet.Config as cms

from HeterogeneousCore.CUDACore.SwitchProducerCUDA import SwitchProducerCUDA

from Configuration.Eras.Era_Run3_cff import Run3
process = cms.Process( "HLTX", Run3 )

process.HLTConfigVersion = cms.PSet(
  tableName = cms.string('/users/sharper/2022/scouting/FullAndPixelScouting/V7')
)

process.transferSystem = cms.PSet( 
  destinations = cms.vstring( 'Tier0',
    'DQM',
    'ECAL',
    'EventDisplay',
    'Lustre',
    'None' ),
  transferModes = cms.vstring( 'default',
    'test',
    'emulator' ),
  streamA = cms.PSet( 
    default = cms.vstring( 'Tier0' ),
    test = cms.vstring( 'Lustre' ),
    emulator = cms.vstring( 'Lustre' )
  ),
  streamCalibration = cms.PSet( 
    default = cms.vstring( 'Tier0' ),
    test = cms.vstring( 'Lustre' ),
    emulator = cms.vstring( 'None' )
  ),
  streamDQM = cms.PSet( 
    default = cms.vstring( 'DQM' ),
    test = cms.vstring( 'DQM',
      'Lustre' ),
    emulator = cms.vstring( 'None' )
  ),
  streamDQMCalibration = cms.PSet( 
    default = cms.vstring( 'DQM' ),
    test = cms.vstring( 'DQM',
      'Lustre' ),
    emulator = cms.vstring( 'None' )
  ),
  streamEcalCalibration = cms.PSet( 
    default = cms.vstring( 'ECAL' ),
    test = cms.vstring( 'ECAL' ),
    emulator = cms.vstring( 'None' )
  ),
  streamEventDisplay = cms.PSet( 
    default = cms.vstring( 'EventDisplay',
      'Tier0' ),
    test = cms.vstring( 'EventDisplay',
      'Lustre' ),
    emulator = cms.vstring( 'None' )
  ),
  streamExpressCosmics = cms.PSet( 
    default = cms.vstring( 'Tier0' ),
    test = cms.vstring( 'Lustre' ),
    emulator = cms.vstring( 'Lustre' )
  ),
  streamNanoDST = cms.PSet( 
    default = cms.vstring( 'Tier0' ),
    test = cms.vstring( 'Lustre' ),
    emulator = cms.vstring( 'None' )
  ),
  streamRPCMON = cms.PSet( 
    default = cms.vstring( 'Tier0' ),
    test = cms.vstring( 'Lustre' ),
    emulator = cms.vstring( 'None' )
  ),
  streamTrackerCalibration = cms.PSet( 
    default = cms.vstring( 'Tier0' ),
    test = cms.vstring( 'Lustre' ),
    emulator = cms.vstring( 'None' )
  ),
  default = cms.PSet( 
    default = cms.vstring( 'Tier0' ),
    test = cms.vstring( 'Lustre' ),
    emulator = cms.vstring( 'Lustre' ),
    streamLookArea = cms.PSet(  )
  ),
  streamLookArea = cms.PSet( 
    default = cms.vstring( 'DQM' ),
    test = cms.vstring( 'DQM',
      'Lustre' ),
    emulator = cms.vstring( 'None' )
  )
)
process.HLTIter0PSetTrajectoryBuilderIT = cms.PSet( 
  ComponentType = cms.string( "CkfTrajectoryBuilder" ),
  lostHitPenalty = cms.double( 30.0 ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTIter0PSetTrajectoryFilterIT" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  maxCand = cms.int32( 2 ),
  alwaysUseInvalidHits = cms.bool( False ),
  estimator = cms.string( "hltESPChi2ChargeMeasurementEstimator9" ),
  intermediateCleaning = cms.bool( True ),
  updator = cms.string( "hltESPKFUpdator" ),
  seedAs5DHit = cms.bool( False )
)
process.HLTIter4PSetTrajectoryBuilderIT = cms.PSet( 
  ComponentType = cms.string( "CkfTrajectoryBuilder" ),
  lostHitPenalty = cms.double( 30.0 ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTIter4PSetTrajectoryFilterIT" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  maxCand = cms.int32( 1 ),
  alwaysUseInvalidHits = cms.bool( False ),
  estimator = cms.string( "hltESPChi2ChargeMeasurementEstimator16" ),
  intermediateCleaning = cms.bool( True ),
  updator = cms.string( "hltESPKFUpdator" ),
  seedAs5DHit = cms.bool( False )
)
process.HLTIter0GroupedCkfTrajectoryBuilderIT = cms.PSet( 
  keepOriginalIfRebuildFails = cms.bool( False ),
  lockHits = cms.bool( True ),
  maxDPhiForLooperReconstruction = cms.double( 2.0 ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTIter0PSetTrajectoryFilterIT" ) ),
  maxCand = cms.int32( 2 ),
  estimator = cms.string( "hltESPChi2ChargeMeasurementEstimator9" ),
  intermediateCleaning = cms.bool( True ),
  bestHitOnly = cms.bool( True ),
  useSameTrajFilter = cms.bool( True ),
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  lostHitPenalty = cms.double( 30.0 ),
  requireSeedHitsInRebuild = cms.bool( True ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  maxPtForLooperReconstruction = cms.double( 0.7 ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  minNrOfHitsForRebuild = cms.int32( 5 ),
  alwaysUseInvalidHits = cms.bool( False ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTIter0PSetTrajectoryFilterIT" ) ),
  foundHitBonus = cms.double( 5.0 ),
  updator = cms.string( "hltESPKFUpdator" ),
  seedAs5DHit = cms.bool( False )
)
process.HLTIter4PSetTrajectoryFilterIT = cms.PSet( 
  minimumNumberOfHits = cms.int32( 6 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 9999 ),
  minPt = cms.double( 0.3 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 1.0 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 999.0 ),
  maxLostHits = cms.int32( 0 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetPvClusterComparerForIT = cms.PSet( 
  track_chi2_max = cms.double( 20.0 ),
  track_pt_max = cms.double( 20.0 ),
  track_prob_min = cms.double( -1.0 ),
  track_pt_min = cms.double( 1.0 )
)
process.HLTPSetMixedStepTrajectoryFilter = cms.PSet( 
  minimumNumberOfHits = cms.int32( 3 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 9999 ),
  minPt = cms.double( 0.1 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 1.4 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  maxLostHits = cms.int32( 999 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetMuonCkfTrajectoryBuilder = cms.PSet( 
  rescaleErrorIfFail = cms.double( 1.0 ),
  ComponentType = cms.string( "MuonCkfTrajectoryBuilder" ),
  lostHitPenalty = cms.double( 30.0 ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetMuonCkfTrajectoryFilter" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterial" ),
  maxCand = cms.int32( 5 ),
  alwaysUseInvalidHits = cms.bool( True ),
  estimator = cms.string( "hltESPChi2ChargeMeasurementEstimator30" ),
  intermediateCleaning = cms.bool( False ),
  propagatorProximity = cms.string( "SteppingHelixPropagatorAny" ),
  updator = cms.string( "hltESPKFUpdator" ),
  deltaEta = cms.double( -1.0 ),
  useSeedLayer = cms.bool( False ),
  deltaPhi = cms.double( -1.0 ),
  seedAs5DHit = cms.bool( False )
)
process.HLTIter0HighPtTkMuPSetTrajectoryFilterIT = cms.PSet( 
  minimumNumberOfHits = cms.int32( 3 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 9999 ),
  minPt = cms.double( 0.3 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 1.0 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 999.0 ),
  maxLostHits = cms.int32( 1 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetPvClusterComparerForBTag = cms.PSet( 
  track_chi2_max = cms.double( 20.0 ),
  track_pt_max = cms.double( 20.0 ),
  track_prob_min = cms.double( -1.0 ),
  track_pt_min = cms.double( 0.1 )
)
process.HLTIter2GroupedCkfTrajectoryBuilderIT = cms.PSet( 
  keepOriginalIfRebuildFails = cms.bool( False ),
  lockHits = cms.bool( True ),
  maxDPhiForLooperReconstruction = cms.double( 2.0 ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTIter2PSetTrajectoryFilterIT" ) ),
  maxCand = cms.int32( 2 ),
  estimator = cms.string( "hltESPChi2ChargeMeasurementEstimator16" ),
  intermediateCleaning = cms.bool( True ),
  bestHitOnly = cms.bool( True ),
  useSameTrajFilter = cms.bool( True ),
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  lostHitPenalty = cms.double( 30.0 ),
  requireSeedHitsInRebuild = cms.bool( True ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  maxPtForLooperReconstruction = cms.double( 0.7 ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  minNrOfHitsForRebuild = cms.int32( 5 ),
  alwaysUseInvalidHits = cms.bool( False ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTIter2PSetTrajectoryFilterIT" ) ),
  foundHitBonus = cms.double( 5.0 ),
  updator = cms.string( "hltESPKFUpdator" ),
  seedAs5DHit = cms.bool( False )
)
process.HLTSiStripClusterChargeCutTight = cms.PSet(  value = cms.double( 1945.0 ) )
process.HLTPSetDetachedStepTrajectoryFilterBase = cms.PSet( 
  minimumNumberOfHits = cms.int32( 3 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 2 ),
  minPt = cms.double( 0.075 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 2.0 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutLoose" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  maxLostHits = cms.int32( 999 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetMuonTrackingRegionBuilder8356 = cms.PSet( 
  Rescale_Dz = cms.double( 3.0 ),
  Pt_fixed = cms.bool( False ),
  Eta_fixed = cms.bool( False ),
  Eta_min = cms.double( 0.1 ),
  DeltaZ = cms.double( 15.9 ),
  maxRegions = cms.int32( 2 ),
  EtaR_UpperLimit_Par1 = cms.double( 0.25 ),
  UseVertex = cms.bool( False ),
  Z_fixed = cms.bool( True ),
  PhiR_UpperLimit_Par1 = cms.double( 0.6 ),
  PhiR_UpperLimit_Par2 = cms.double( 0.2 ),
  Rescale_phi = cms.double( 3.0 ),
  DeltaEta = cms.double( 0.2 ),
  precise = cms.bool( True ),
  OnDemand = cms.int32( -1 ),
  EtaR_UpperLimit_Par2 = cms.double( 0.15 ),
  MeasurementTrackerName = cms.InputTag( "hltESPMeasurementTracker" ),
  vertexCollection = cms.InputTag( "pixelVertices" ),
  Pt_min = cms.double( 1.5 ),
  beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
  Phi_fixed = cms.bool( False ),
  DeltaR = cms.double( 0.2 ),
  input = cms.InputTag( 'hltL2Muons','UpdatedAtVtx' ),
  DeltaPhi = cms.double( 0.2 ),
  Phi_min = cms.double( 0.1 ),
  Rescale_eta = cms.double( 3.0 )
)
process.HLTPSetJetCoreStepTrajectoryFilter = cms.PSet( 
  minimumNumberOfHits = cms.int32( 4 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 9999 ),
  minPt = cms.double( 0.1 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 2.0 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  maxLostHits = cms.int32( 999 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTIter2PSetTrajectoryFilterIT = cms.PSet( 
  minimumNumberOfHits = cms.int32( 3 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 1 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 0 ),
  minPt = cms.double( 0.3 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 1.0 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 999.0 ),
  maxLostHits = cms.int32( 1 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetMuTrackJpsiTrajectoryBuilder = cms.PSet( 
  ComponentType = cms.string( "CkfTrajectoryBuilder" ),
  lostHitPenalty = cms.double( 30.0 ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetMuTrackJpsiTrajectoryFilter" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterial" ),
  maxCand = cms.int32( 1 ),
  alwaysUseInvalidHits = cms.bool( False ),
  estimator = cms.string( "hltESPChi2ChargeMeasurementEstimator30" ),
  intermediateCleaning = cms.bool( True ),
  updator = cms.string( "hltESPKFUpdator" ),
  seedAs5DHit = cms.bool( False )
)
process.HLTPSetTrajectoryBuilderForGsfElectrons = cms.PSet( 
  ComponentType = cms.string( "CkfTrajectoryBuilder" ),
  lostHitPenalty = cms.double( 90.0 ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  propagatorOpposite = cms.string( "hltESPBwdElectronPropagator" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetTrajectoryFilterForElectrons" ) ),
  propagatorAlong = cms.string( "hltESPFwdElectronPropagator" ),
  maxCand = cms.int32( 5 ),
  alwaysUseInvalidHits = cms.bool( True ),
  estimator = cms.string( "hltESPChi2ChargeMeasurementEstimator2000" ),
  intermediateCleaning = cms.bool( False ),
  updator = cms.string( "hltESPKFUpdator" ),
  seedAs5DHit = cms.bool( False )
)
process.HLTSiStripClusterChargeCutNone = cms.PSet(  value = cms.double( -1.0 ) )
process.HLTPSetMuonCkfTrajectoryFilter = cms.PSet( 
  minimumNumberOfHits = cms.int32( 5 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 9999 ),
  minPt = cms.double( 0.9 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 1.0 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( -1 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 999.0 ),
  maxLostHits = cms.int32( 1 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetDetachedStepTrajectoryFilter = cms.PSet( 
  ComponentType = cms.string( "CompositeTrajectoryFilter" ),
  filters = cms.VPSet( 
    cms.PSet(  refToPSet_ = cms.string( "HLTPSetDetachedStepTrajectoryFilterBase" )    )
  )
)
process.HLTIter1PSetTrajectoryFilterIT = cms.PSet( 
  minimumNumberOfHits = cms.int32( 3 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 0 ),
  minPt = cms.double( 0.2 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 1.0 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 999.0 ),
  maxLostHits = cms.int32( 1 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetMixedStepTrajectoryBuilder = cms.PSet( 
  useSameTrajFilter = cms.bool( True ),
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  keepOriginalIfRebuildFails = cms.bool( False ),
  lostHitPenalty = cms.double( 30.0 ),
  lockHits = cms.bool( True ),
  requireSeedHitsInRebuild = cms.bool( True ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  maxDPhiForLooperReconstruction = cms.double( 2.0 ),
  maxPtForLooperReconstruction = cms.double( 0.7 ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialForMixedStepOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetMixedStepTrajectoryFilter" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterialForMixedStep" ),
  minNrOfHitsForRebuild = cms.int32( 5 ),
  maxCand = cms.int32( 2 ),
  alwaysUseInvalidHits = cms.bool( True ),
  estimator = cms.string( "hltESPChi2ChargeTightMeasurementEstimator16" ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetMixedStepTrajectoryFilter" ) ),
  intermediateCleaning = cms.bool( True ),
  foundHitBonus = cms.double( 5.0 ),
  updator = cms.string( "hltESPKFUpdator" ),
  bestHitOnly = cms.bool( True ),
  seedAs5DHit = cms.bool( False )
)
process.HLTSeedFromProtoTracks = cms.PSet( 
  TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
  SeedMomentumForBOFF = cms.double( 5.0 ),
  propagator = cms.string( "PropagatorWithMaterialParabolicMf" ),
  forceKinematicWithRegionDirection = cms.bool( False ),
  magneticField = cms.string( "ParabolicMf" ),
  OriginTransverseErrorMultiplier = cms.double( 1.0 ),
  ComponentName = cms.string( "SeedFromConsecutiveHitsCreator" ),
  MinOneOverPtError = cms.double( 1.0 )
)
process.HLTIter2PSetTrajectoryBuilderIT = cms.PSet( 
  ComponentType = cms.string( "CkfTrajectoryBuilder" ),
  lostHitPenalty = cms.double( 30.0 ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTIter2PSetTrajectoryFilterIT" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  maxCand = cms.int32( 2 ),
  alwaysUseInvalidHits = cms.bool( False ),
  estimator = cms.string( "hltESPChi2ChargeMeasurementEstimator16" ),
  intermediateCleaning = cms.bool( True ),
  updator = cms.string( "hltESPKFUpdator" ),
  seedAs5DHit = cms.bool( False )
)
process.HLTPSetMuTrackJpsiTrajectoryFilter = cms.PSet( 
  minimumNumberOfHits = cms.int32( 5 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 9999 ),
  minPt = cms.double( 10.0 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 1.0 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( 8 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 999.0 ),
  maxLostHits = cms.int32( 1 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetDetachedStepTrajectoryBuilder = cms.PSet( 
  useSameTrajFilter = cms.bool( True ),
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  keepOriginalIfRebuildFails = cms.bool( False ),
  lostHitPenalty = cms.double( 30.0 ),
  lockHits = cms.bool( True ),
  requireSeedHitsInRebuild = cms.bool( True ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  maxDPhiForLooperReconstruction = cms.double( 2.0 ),
  maxPtForLooperReconstruction = cms.double( 0.7 ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetDetachedStepTrajectoryFilter" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  minNrOfHitsForRebuild = cms.int32( 5 ),
  maxCand = cms.int32( 3 ),
  alwaysUseInvalidHits = cms.bool( False ),
  estimator = cms.string( "hltESPChi2ChargeMeasurementEstimator9" ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetDetachedStepTrajectoryFilter" ) ),
  intermediateCleaning = cms.bool( True ),
  foundHitBonus = cms.double( 5.0 ),
  updator = cms.string( "hltESPKFUpdator" ),
  bestHitOnly = cms.bool( True ),
  seedAs5DHit = cms.bool( False )
)
process.HLTIter0PSetTrajectoryFilterIT = cms.PSet( 
  minimumNumberOfHits = cms.int32( 3 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 0 ),
  minPt = cms.double( 0.3 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 1.0 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 999.0 ),
  maxLostHits = cms.int32( 1 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetLowPtStepTrajectoryFilter = cms.PSet( 
  minimumNumberOfHits = cms.int32( 3 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 1 ),
  minPt = cms.double( 0.075 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 2.0 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutLoose" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  maxLostHits = cms.int32( 999 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTSeedFromConsecutiveHitsCreator = cms.PSet( 
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  SeedMomentumForBOFF = cms.double( 5.0 ),
  propagator = cms.string( "PropagatorWithMaterial" ),
  forceKinematicWithRegionDirection = cms.bool( False ),
  magneticField = cms.string( "" ),
  OriginTransverseErrorMultiplier = cms.double( 1.0 ),
  ComponentName = cms.string( "SeedFromConsecutiveHitsCreator" ),
  MinOneOverPtError = cms.double( 1.0 )
)
process.HLTIter1PSetTrajectoryBuilderIT = cms.PSet( 
  ComponentType = cms.string( "CkfTrajectoryBuilder" ),
  lostHitPenalty = cms.double( 30.0 ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTIter1PSetTrajectoryFilterIT" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  maxCand = cms.int32( 2 ),
  alwaysUseInvalidHits = cms.bool( False ),
  estimator = cms.string( "hltESPChi2ChargeMeasurementEstimator16" ),
  intermediateCleaning = cms.bool( True ),
  updator = cms.string( "hltESPKFUpdator" ),
  seedAs5DHit = cms.bool( False )
)
process.HLTSiStripClusterChargeCutForHI = cms.PSet(  value = cms.double( 2069.0 ) )
process.HLTPSetLowPtStepTrajectoryBuilder = cms.PSet( 
  useSameTrajFilter = cms.bool( True ),
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  keepOriginalIfRebuildFails = cms.bool( False ),
  lostHitPenalty = cms.double( 30.0 ),
  lockHits = cms.bool( True ),
  requireSeedHitsInRebuild = cms.bool( True ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  maxDPhiForLooperReconstruction = cms.double( 2.0 ),
  maxPtForLooperReconstruction = cms.double( 0.7 ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetLowPtStepTrajectoryFilter" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  minNrOfHitsForRebuild = cms.int32( 5 ),
  maxCand = cms.int32( 4 ),
  alwaysUseInvalidHits = cms.bool( True ),
  estimator = cms.string( "hltESPChi2ChargeMeasurementEstimator9" ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetLowPtStepTrajectoryFilter" ) ),
  intermediateCleaning = cms.bool( True ),
  foundHitBonus = cms.double( 5.0 ),
  updator = cms.string( "hltESPKFUpdator" ),
  bestHitOnly = cms.bool( True ),
  seedAs5DHit = cms.bool( False )
)
process.HLTPSetTrajectoryFilterForElectrons = cms.PSet( 
  minimumNumberOfHits = cms.int32( 5 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 9999 ),
  minPt = cms.double( 2.0 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 1.0 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( -1 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( -1 ),
  maxLostHitsFraction = cms.double( 999.0 ),
  maxLostHits = cms.int32( 1 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetJetCoreStepTrajectoryBuilder = cms.PSet( 
  useSameTrajFilter = cms.bool( True ),
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  keepOriginalIfRebuildFails = cms.bool( False ),
  lostHitPenalty = cms.double( 30.0 ),
  lockHits = cms.bool( True ),
  requireSeedHitsInRebuild = cms.bool( True ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  maxDPhiForLooperReconstruction = cms.double( 2.0 ),
  maxPtForLooperReconstruction = cms.double( 0.7 ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetJetCoreStepTrajectoryFilter" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  minNrOfHitsForRebuild = cms.int32( 5 ),
  maxCand = cms.int32( 50 ),
  alwaysUseInvalidHits = cms.bool( True ),
  estimator = cms.string( "hltESPChi2MeasurementEstimator30" ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetJetCoreStepTrajectoryFilter" ) ),
  intermediateCleaning = cms.bool( True ),
  foundHitBonus = cms.double( 5.0 ),
  updator = cms.string( "hltESPKFUpdator" ),
  bestHitOnly = cms.bool( True ),
  seedAs5DHit = cms.bool( False )
)
process.HLTPSetPvClusterComparer = cms.PSet( 
  track_chi2_max = cms.double( 9999999.0 ),
  track_pt_max = cms.double( 10.0 ),
  track_prob_min = cms.double( -1.0 ),
  track_pt_min = cms.double( 2.5 )
)
process.HLTIter0HighPtTkMuPSetTrajectoryBuilderIT = cms.PSet( 
  ComponentType = cms.string( "CkfTrajectoryBuilder" ),
  lostHitPenalty = cms.double( 30.0 ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTIter0HighPtTkMuPSetTrajectoryFilterIT" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  maxCand = cms.int32( 4 ),
  alwaysUseInvalidHits = cms.bool( True ),
  estimator = cms.string( "hltESPChi2ChargeMeasurementEstimator30" ),
  intermediateCleaning = cms.bool( True ),
  updator = cms.string( "hltESPKFUpdator" ),
  seedAs5DHit = cms.bool( False )
)
process.HLTIter1GroupedCkfTrajectoryBuilderIT = cms.PSet( 
  useSameTrajFilter = cms.bool( True ),
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  keepOriginalIfRebuildFails = cms.bool( False ),
  lostHitPenalty = cms.double( 30.0 ),
  lockHits = cms.bool( True ),
  requireSeedHitsInRebuild = cms.bool( True ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTIter1PSetTrajectoryFilterIT" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  minNrOfHitsForRebuild = cms.int32( 5 ),
  maxCand = cms.int32( 2 ),
  alwaysUseInvalidHits = cms.bool( False ),
  estimator = cms.string( "hltESPChi2ChargeMeasurementEstimator16" ),
  intermediateCleaning = cms.bool( True ),
  foundHitBonus = cms.double( 5.0 ),
  updator = cms.string( "hltESPKFUpdator" ),
  bestHitOnly = cms.bool( True ),
  seedAs5DHit = cms.bool( False ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTIter1PSetTrajectoryFilterIT" ) )
)
process.HLTPSetInitialStepTrajectoryBuilder = cms.PSet( 
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  bestHitOnly = cms.bool( True ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetInitialStepTrajectoryFilter" ) ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetInitialStepTrajectoryFilter" ) ),
  useSameTrajFilter = cms.bool( True ),
  maxCand = cms.int32( 3 ),
  intermediateCleaning = cms.bool( True ),
  lostHitPenalty = cms.double( 30.0 ),
  foundHitBonus = cms.double( 10.0 ),
  lockHits = cms.bool( True ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  updator = cms.string( "hltESPKFUpdator" ),
  alwaysUseInvalidHits = cms.bool( True ),
  requireSeedHitsInRebuild = cms.bool( True ),
  keepOriginalIfRebuildFails = cms.bool( True ),
  estimator = cms.string( "hltESPInitialStepChi2ChargeMeasurementEstimator30" ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  minNrOfHitsForRebuild = cms.int32( 1 ),
  maxDPhiForLooperReconstruction = cms.double( 2.0 ),
  maxPtForLooperReconstruction = cms.double( 0.7 ),
  seedAs5DHit = cms.bool( False )
)
process.HLTPSetInitialStepTrajectoryFilter = cms.PSet( 
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  minimumNumberOfHits = cms.int32( 3 ),
  seedPairPenalty = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  minPt = cms.double( 0.2 ),
  nSigmaMinPt = cms.double( 5.0 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHits = cms.int32( 999 ),
  maxConsecLostHits = cms.int32( 1 ),
  maxNumberOfHits = cms.int32( 100 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  constantValueForLostHitsFractionFilter = cms.double( 2.0 ),
  seedExtension = cms.int32( 0 ),
  strictSeedExtension = cms.bool( False ),
  pixelSeedExtension = cms.bool( False ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  maxCCCLostHits = cms.int32( 0 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutLoose" ) ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetPixelLessStepTrajectoryBuilder = cms.PSet( 
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  bestHitOnly = cms.bool( True ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetPixelLessStepTrajectoryFilter" ) ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetPixelLessStepTrajectoryFilter" ) ),
  useSameTrajFilter = cms.bool( True ),
  maxCand = cms.int32( 2 ),
  intermediateCleaning = cms.bool( True ),
  lostHitPenalty = cms.double( 30.0 ),
  foundHitBonus = cms.double( 10.0 ),
  lockHits = cms.bool( True ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  updator = cms.string( "hltESPKFUpdator" ),
  alwaysUseInvalidHits = cms.bool( False ),
  requireSeedHitsInRebuild = cms.bool( True ),
  keepOriginalIfRebuildFails = cms.bool( False ),
  estimator = cms.string( "hltESPPixelLessStepChi2ChargeMeasurementEstimator16" ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  minNrOfHitsForRebuild = cms.int32( 4 ),
  maxDPhiForLooperReconstruction = cms.double( 2.0 ),
  maxPtForLooperReconstruction = cms.double( 0.7 ),
  seedAs5DHit = cms.bool( False )
)
process.HLTPSetPixelLessStepTrajectoryFilter = cms.PSet( 
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  minimumNumberOfHits = cms.int32( 4 ),
  seedPairPenalty = cms.int32( 1 ),
  chargeSignificance = cms.double( -1.0 ),
  minPt = cms.double( 0.1 ),
  nSigmaMinPt = cms.double( 5.0 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHits = cms.int32( 0 ),
  maxConsecLostHits = cms.int32( 1 ),
  maxNumberOfHits = cms.int32( 100 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  constantValueForLostHitsFractionFilter = cms.double( 2.0 ),
  seedExtension = cms.int32( 0 ),
  strictSeedExtension = cms.bool( False ),
  pixelSeedExtension = cms.bool( False ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  maxCCCLostHits = cms.int32( 9999 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetTobTecStepTrajectoryFilter = cms.PSet( 
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  minimumNumberOfHits = cms.int32( 5 ),
  seedPairPenalty = cms.int32( 1 ),
  chargeSignificance = cms.double( -1.0 ),
  minPt = cms.double( 0.1 ),
  nSigmaMinPt = cms.double( 5.0 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHits = cms.int32( 0 ),
  maxConsecLostHits = cms.int32( 1 ),
  maxNumberOfHits = cms.int32( 100 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  constantValueForLostHitsFractionFilter = cms.double( 2.0 ),
  seedExtension = cms.int32( 0 ),
  strictSeedExtension = cms.bool( False ),
  pixelSeedExtension = cms.bool( False ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  maxCCCLostHits = cms.int32( 9999 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetTobTecStepInOutTrajectoryFilter = cms.PSet( 
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  minimumNumberOfHits = cms.int32( 4 ),
  seedPairPenalty = cms.int32( 1 ),
  chargeSignificance = cms.double( -1.0 ),
  minPt = cms.double( 0.1 ),
  nSigmaMinPt = cms.double( 5.0 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHits = cms.int32( 0 ),
  maxConsecLostHits = cms.int32( 1 ),
  maxNumberOfHits = cms.int32( 100 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  constantValueForLostHitsFractionFilter = cms.double( 2.0 ),
  seedExtension = cms.int32( 0 ),
  strictSeedExtension = cms.bool( False ),
  pixelSeedExtension = cms.bool( False ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  maxCCCLostHits = cms.int32( 9999 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetTobTecStepTrajectoryBuilder = cms.PSet( 
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  bestHitOnly = cms.bool( True ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetTobTecStepTrajectoryFilter" ) ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetTobTecStepInOutTrajectoryFilter" ) ),
  useSameTrajFilter = cms.bool( False ),
  maxCand = cms.int32( 2 ),
  intermediateCleaning = cms.bool( True ),
  lostHitPenalty = cms.double( 30.0 ),
  foundHitBonus = cms.double( 10.0 ),
  lockHits = cms.bool( True ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  updator = cms.string( "hltESPKFUpdator" ),
  alwaysUseInvalidHits = cms.bool( False ),
  requireSeedHitsInRebuild = cms.bool( True ),
  keepOriginalIfRebuildFails = cms.bool( False ),
  estimator = cms.string( "hltESPTobTecStepChi2ChargeMeasurementEstimator16" ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  minNrOfHitsForRebuild = cms.int32( 4 ),
  maxDPhiForLooperReconstruction = cms.double( 2.0 ),
  maxPtForLooperReconstruction = cms.double( 0.7 ),
  seedAs5DHit = cms.bool( False )
)
process.HLTIter0IterL3MuonPSetGroupedCkfTrajectoryBuilderIT = cms.PSet( 
  useSameTrajFilter = cms.bool( True ),
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  keepOriginalIfRebuildFails = cms.bool( True ),
  lostHitPenalty = cms.double( 1.0 ),
  lockHits = cms.bool( True ),
  requireSeedHitsInRebuild = cms.bool( True ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTIter0IterL3MuonGroupedCkfTrajectoryFilterIT" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterial" ),
  minNrOfHitsForRebuild = cms.int32( 2 ),
  maxCand = cms.int32( 5 ),
  alwaysUseInvalidHits = cms.bool( True ),
  estimator = cms.string( "hltESPChi2ChargeMeasurementEstimator30" ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTIter0IterL3MuonGroupedCkfTrajectoryFilterIT" ) ),
  intermediateCleaning = cms.bool( True ),
  foundHitBonus = cms.double( 1000.0 ),
  updator = cms.string( "hltESPKFUpdator" ),
  bestHitOnly = cms.bool( True ),
  seedAs5DHit = cms.bool( False )
)
process.HLTIter0IterL3FromL1MuonGroupedCkfTrajectoryFilterIT = cms.PSet( 
  minimumNumberOfHits = cms.int32( 3 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 9999 ),
  minPt = cms.double( 0.9 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 10.0 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  maxLostHits = cms.int32( 999 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTIter0IterL3FromL1MuonPSetGroupedCkfTrajectoryBuilderIT = cms.PSet( 
  useSameTrajFilter = cms.bool( True ),
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  keepOriginalIfRebuildFails = cms.bool( True ),
  lostHitPenalty = cms.double( 1.0 ),
  lockHits = cms.bool( True ),
  requireSeedHitsInRebuild = cms.bool( True ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTIter0IterL3FromL1MuonGroupedCkfTrajectoryFilterIT" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterial" ),
  minNrOfHitsForRebuild = cms.int32( 2 ),
  maxCand = cms.int32( 5 ),
  alwaysUseInvalidHits = cms.bool( True ),
  estimator = cms.string( "hltESPChi2ChargeMeasurementEstimator30" ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTIter0IterL3FromL1MuonGroupedCkfTrajectoryFilterIT" ) ),
  intermediateCleaning = cms.bool( True ),
  foundHitBonus = cms.double( 1000.0 ),
  updator = cms.string( "hltESPKFUpdator" ),
  bestHitOnly = cms.bool( True ),
  seedAs5DHit = cms.bool( False )
)
process.HLTIter0IterL3MuonGroupedCkfTrajectoryFilterIT = cms.PSet( 
  minimumNumberOfHits = cms.int32( 3 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 9999 ),
  minPt = cms.double( 0.9 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 10.0 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  maxLostHits = cms.int32( 999 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTIter2IterL3FromL1MuonPSetTrajectoryFilterIT = cms.PSet( 
  minimumNumberOfHits = cms.int32( 5 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 9999 ),
  minPt = cms.double( 0.3 ),
  maxConsecLostHits = cms.int32( 3 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 1.0 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 999.0 ),
  maxLostHits = cms.int32( 1 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTIter2IterL3FromL1MuonPSetGroupedCkfTrajectoryBuilderIT = cms.PSet( 
  useSameTrajFilter = cms.bool( True ),
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  keepOriginalIfRebuildFails = cms.bool( False ),
  lostHitPenalty = cms.double( 30.0 ),
  lockHits = cms.bool( True ),
  requireSeedHitsInRebuild = cms.bool( False ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTIter2IterL3FromL1MuonPSetTrajectoryFilterIT" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  minNrOfHitsForRebuild = cms.int32( 5 ),
  maxCand = cms.int32( 2 ),
  alwaysUseInvalidHits = cms.bool( False ),
  estimator = cms.string( "hltESPChi2ChargeMeasurementEstimator30" ),
  intermediateCleaning = cms.bool( True ),
  foundHitBonus = cms.double( 1000.0 ),
  updator = cms.string( "hltESPKFUpdator" ),
  bestHitOnly = cms.bool( True ),
  seedAs5DHit = cms.bool( False ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTIter2IterL3FromL1MuonPSetTrajectoryFilterIT" ) )
)
process.HLTIter2IterL3MuonPSetTrajectoryFilterIT = cms.PSet( 
  minimumNumberOfHits = cms.int32( 5 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 9999 ),
  minPt = cms.double( 0.3 ),
  maxConsecLostHits = cms.int32( 3 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 1.0 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 999.0 ),
  maxLostHits = cms.int32( 1 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTIter2IterL3MuonPSetGroupedCkfTrajectoryBuilderIT = cms.PSet( 
  useSameTrajFilter = cms.bool( True ),
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  keepOriginalIfRebuildFails = cms.bool( False ),
  lostHitPenalty = cms.double( 30.0 ),
  lockHits = cms.bool( True ),
  requireSeedHitsInRebuild = cms.bool( False ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTIter2IterL3MuonPSetTrajectoryFilterIT" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  minNrOfHitsForRebuild = cms.int32( 5 ),
  maxCand = cms.int32( 2 ),
  alwaysUseInvalidHits = cms.bool( False ),
  estimator = cms.string( "hltESPChi2ChargeMeasurementEstimator30" ),
  intermediateCleaning = cms.bool( True ),
  foundHitBonus = cms.double( 1000.0 ),
  updator = cms.string( "hltESPKFUpdator" ),
  bestHitOnly = cms.bool( True ),
  seedAs5DHit = cms.bool( False ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTIter2IterL3MuonPSetTrajectoryFilterIT" ) )
)
process.HLTPSetPixelPairStepTrajectoryFilter = cms.PSet( 
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  minimumNumberOfHits = cms.int32( 4 ),
  seedPairPenalty = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  minPt = cms.double( 0.1 ),
  nSigmaMinPt = cms.double( 5.0 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHits = cms.int32( 999 ),
  maxConsecLostHits = cms.int32( 1 ),
  maxNumberOfHits = cms.int32( 100 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  constantValueForLostHitsFractionFilter = cms.double( 2.0 ),
  seedExtension = cms.int32( 0 ),
  strictSeedExtension = cms.bool( False ),
  pixelSeedExtension = cms.bool( False ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  maxCCCLostHits = cms.int32( 0 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutLoose" ) ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetPixelPairStepTrajectoryFilterInOut = cms.PSet( 
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  minimumNumberOfHits = cms.int32( 4 ),
  seedPairPenalty = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  minPt = cms.double( 0.1 ),
  nSigmaMinPt = cms.double( 5.0 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHits = cms.int32( 999 ),
  maxConsecLostHits = cms.int32( 1 ),
  maxNumberOfHits = cms.int32( 100 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  constantValueForLostHitsFractionFilter = cms.double( 2.0 ),
  seedExtension = cms.int32( 1 ),
  strictSeedExtension = cms.bool( False ),
  pixelSeedExtension = cms.bool( False ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  maxCCCLostHits = cms.int32( 0 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutLoose" ) ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetPixelPairStepTrajectoryBuilder = cms.PSet( 
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  bestHitOnly = cms.bool( True ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetPixelPairStepTrajectoryFilter" ) ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetPixelPairStepTrajectoryFilterInOut" ) ),
  useSameTrajFilter = cms.bool( False ),
  maxCand = cms.int32( 3 ),
  intermediateCleaning = cms.bool( True ),
  lostHitPenalty = cms.double( 30.0 ),
  foundHitBonus = cms.double( 10.0 ),
  lockHits = cms.bool( True ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  updator = cms.string( "hltESPKFUpdator" ),
  alwaysUseInvalidHits = cms.bool( True ),
  requireSeedHitsInRebuild = cms.bool( True ),
  keepOriginalIfRebuildFails = cms.bool( False ),
  estimator = cms.string( "hltESPPixelPairStepChi2ChargeMeasurementEstimator9" ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  minNrOfHitsForRebuild = cms.int32( 5 ),
  maxDPhiForLooperReconstruction = cms.double( 2.0 ),
  maxPtForLooperReconstruction = cms.double( 0.7 ),
  seedAs5DHit = cms.bool( False )
)
process.HLTPSetTobTecStepTrajectoryBuilderPPOnAA = cms.PSet( 
  useSameTrajFilter = cms.bool( False ),
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  keepOriginalIfRebuildFails = cms.bool( False ),
  lostHitPenalty = cms.double( 30.0 ),
  lockHits = cms.bool( True ),
  requireSeedHitsInRebuild = cms.bool( True ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  maxDPhiForLooperReconstruction = cms.double( 2.0 ),
  maxPtForLooperReconstruction = cms.double( 0.7 ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetTobTecStepTrajectoryFilterPPOnAA" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  minNrOfHitsForRebuild = cms.int32( 4 ),
  maxCand = cms.int32( 2 ),
  alwaysUseInvalidHits = cms.bool( False ),
  estimator = cms.string( "hltESPTobTecStepChi2ChargeMeasurementEstimator16" ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetTobTecStepInOutTrajectoryFilterPPOnAA" ) ),
  intermediateCleaning = cms.bool( True ),
  foundHitBonus = cms.double( 10.0 ),
  updator = cms.string( "hltESPKFUpdator" ),
  bestHitOnly = cms.bool( True ),
  seedAs5DHit = cms.bool( False )
)
process.HLTPSetLowPtTripletStepTrajectoryFilterPPOnAA = cms.PSet( 
  minimumNumberOfHits = cms.int32( 3 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 0 ),
  minPt = cms.double( 0.49 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 2.0 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutLoose" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  maxLostHits = cms.int32( 999 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetDetachedQuadStepTrajectoryFilterPPOnAA = cms.PSet( 
  minimumNumberOfHits = cms.int32( 3 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 0 ),
  minPt = cms.double( 0.9 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 2.0 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutLoose" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  maxLostHits = cms.int32( 999 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetCkfBaseTrajectoryFilter_block = cms.PSet( 
  minimumNumberOfHits = cms.int32( 5 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 9999 ),
  minPt = cms.double( 0.9 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 2.0 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  maxLostHits = cms.int32( 999 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetInitialStepTrajectoryBuilderPPOnAA = cms.PSet( 
  useSameTrajFilter = cms.bool( True ),
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  keepOriginalIfRebuildFails = cms.bool( True ),
  lostHitPenalty = cms.double( 30.0 ),
  lockHits = cms.bool( True ),
  requireSeedHitsInRebuild = cms.bool( True ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  maxDPhiForLooperReconstruction = cms.double( 2.0 ),
  maxPtForLooperReconstruction = cms.double( 0.7 ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetInitialStepTrajectoryFilterPPOnAA" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  minNrOfHitsForRebuild = cms.int32( 1 ),
  maxCand = cms.int32( 3 ),
  alwaysUseInvalidHits = cms.bool( True ),
  estimator = cms.string( "hltESPInitialStepChi2ChargeMeasurementEstimator30" ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetCkfBaseTrajectoryFilter_block" ) ),
  intermediateCleaning = cms.bool( True ),
  foundHitBonus = cms.double( 10.0 ),
  updator = cms.string( "hltESPKFUpdator" ),
  bestHitOnly = cms.bool( True ),
  seedAs5DHit = cms.bool( False )
)
process.HLTPSetLowPtTripletStepTrajectoryBuilderPPOnAA = cms.PSet( 
  useSameTrajFilter = cms.bool( True ),
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  keepOriginalIfRebuildFails = cms.bool( False ),
  lostHitPenalty = cms.double( 30.0 ),
  lockHits = cms.bool( True ),
  requireSeedHitsInRebuild = cms.bool( True ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  maxDPhiForLooperReconstruction = cms.double( 2.0 ),
  maxPtForLooperReconstruction = cms.double( 0.7 ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetLowPtTripletStepTrajectoryFilterPPOnAA" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  minNrOfHitsForRebuild = cms.int32( 5 ),
  maxCand = cms.int32( 4 ),
  alwaysUseInvalidHits = cms.bool( True ),
  estimator = cms.string( "hltESPLowPtTripletStepChi2ChargeMeasurementEstimator9" ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetCkfBaseTrajectoryFilter_block" ) ),
  intermediateCleaning = cms.bool( True ),
  foundHitBonus = cms.double( 10.0 ),
  updator = cms.string( "hltESPKFUpdator" ),
  bestHitOnly = cms.bool( True ),
  seedAs5DHit = cms.bool( False )
)
process.HLTPSetTobTecStepInOutTrajectoryFilterPPOnAA = cms.PSet( 
  minimumNumberOfHits = cms.int32( 4 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 9999 ),
  minPt = cms.double( 2.0 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 2.0 ),
  seedPairPenalty = cms.int32( 1 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  maxLostHits = cms.int32( 0 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetInitialStepTrajectoryFilterBasePreSplittingPPOnAA = cms.PSet( 
  minimumNumberOfHits = cms.int32( 4 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 0 ),
  minPt = cms.double( 0.2 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 2.0 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutLoose" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  maxLostHits = cms.int32( 999 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetDetachedTripletStepTrajectoryBuilderPPOnAA = cms.PSet( 
  useSameTrajFilter = cms.bool( True ),
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  keepOriginalIfRebuildFails = cms.bool( False ),
  lostHitPenalty = cms.double( 30.0 ),
  lockHits = cms.bool( True ),
  requireSeedHitsInRebuild = cms.bool( True ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  maxDPhiForLooperReconstruction = cms.double( 2.0 ),
  maxPtForLooperReconstruction = cms.double( 0.7 ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetDetachedTripletStepTrajectoryFilterPPOnAA" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  minNrOfHitsForRebuild = cms.int32( 5 ),
  maxCand = cms.int32( 3 ),
  alwaysUseInvalidHits = cms.bool( True ),
  estimator = cms.string( "hltESPDetachedTripletStepChi2ChargeMeasurementEstimator9" ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetCkfBaseTrajectoryFilter_block" ) ),
  intermediateCleaning = cms.bool( True ),
  foundHitBonus = cms.double( 10.0 ),
  updator = cms.string( "hltESPKFUpdator" ),
  bestHitOnly = cms.bool( True ),
  seedAs5DHit = cms.bool( False )
)
process.HLTPSetHighPtTripletStepTrajectoryBuilderPPOnAA = cms.PSet( 
  useSameTrajFilter = cms.bool( True ),
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  keepOriginalIfRebuildFails = cms.bool( False ),
  lostHitPenalty = cms.double( 30.0 ),
  lockHits = cms.bool( True ),
  requireSeedHitsInRebuild = cms.bool( True ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  maxDPhiForLooperReconstruction = cms.double( 2.0 ),
  maxPtForLooperReconstruction = cms.double( 0.7 ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetHighPtTripletStepTrajectoryFilterPPOnAA" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  minNrOfHitsForRebuild = cms.int32( 5 ),
  maxCand = cms.int32( 3 ),
  alwaysUseInvalidHits = cms.bool( True ),
  estimator = cms.string( "hltESPHighPtTripletStepChi2ChargeMeasurementEstimator30" ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetCkfBaseTrajectoryFilter_block" ) ),
  intermediateCleaning = cms.bool( True ),
  foundHitBonus = cms.double( 10.0 ),
  updator = cms.string( "hltESPKFUpdator" ),
  bestHitOnly = cms.bool( True ),
  seedAs5DHit = cms.bool( False )
)
process.HLTPSetDetachedTripletStepTrajectoryFilterPPOnAA = cms.PSet( 
  minimumNumberOfHits = cms.int32( 3 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 0 ),
  minPt = cms.double( 0.9 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 2.0 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutLoose" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  maxLostHits = cms.int32( 999 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetInitialStepTrajectoryBuilderPreSplittingPPOnAA = cms.PSet( 
  useSameTrajFilter = cms.bool( True ),
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  keepOriginalIfRebuildFails = cms.bool( False ),
  lostHitPenalty = cms.double( 30.0 ),
  lockHits = cms.bool( True ),
  requireSeedHitsInRebuild = cms.bool( True ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  maxDPhiForLooperReconstruction = cms.double( 2.0 ),
  maxPtForLooperReconstruction = cms.double( 0.7 ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetInitialStepTrajectoryFilterPreSplittingPPOnAA" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  minNrOfHitsForRebuild = cms.int32( 5 ),
  maxCand = cms.int32( 3 ),
  alwaysUseInvalidHits = cms.bool( True ),
  estimator = cms.string( "hltESPInitialStepChi2ChargeMeasurementEstimator30" ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetCkfBaseTrajectoryFilter_block" ) ),
  intermediateCleaning = cms.bool( True ),
  foundHitBonus = cms.double( 10.0 ),
  updator = cms.string( "hltESPKFUpdator" ),
  bestHitOnly = cms.bool( True ),
  seedAs5DHit = cms.bool( False )
)
process.HLTPSetLowPtQuadStepTrajectoryFilterPPOnAA = cms.PSet( 
  minimumNumberOfHits = cms.int32( 3 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 0 ),
  minPt = cms.double( 0.49 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 2.0 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutLoose" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  maxLostHits = cms.int32( 999 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetMixedTripletStepTrajectoryFilterPPOnAA = cms.PSet( 
  minimumNumberOfHits = cms.int32( 3 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 9999 ),
  minPt = cms.double( 0.4 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 1.4 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  maxLostHits = cms.int32( 999 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetInitialStepTrajectoryFilterPreSplittingPPOnAA = cms.PSet( 
  ComponentType = cms.string( "CompositeTrajectoryFilter" ),
  filters = cms.VPSet( 
    cms.PSet(  refToPSet_ = cms.string( "HLTPSetInitialStepTrajectoryFilterBasePreSplittingPPOnAA" )    ),
    cms.PSet(  refToPSet_ = cms.string( "HLTPSetInitialStepTrajectoryFilterShapePreSplittingPPOnAA" )    )
  )
)
process.HLTPSetMixedTripletStepTrajectoryBuilderPPOnAA = cms.PSet( 
  useSameTrajFilter = cms.bool( True ),
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  keepOriginalIfRebuildFails = cms.bool( False ),
  lostHitPenalty = cms.double( 30.0 ),
  lockHits = cms.bool( True ),
  requireSeedHitsInRebuild = cms.bool( True ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  maxDPhiForLooperReconstruction = cms.double( 2.0 ),
  maxPtForLooperReconstruction = cms.double( 0.7 ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialForMixedStepOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetMixedTripletStepTrajectoryFilterPPOnAA" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterialForMixedStep" ),
  minNrOfHitsForRebuild = cms.int32( 5 ),
  maxCand = cms.int32( 2 ),
  alwaysUseInvalidHits = cms.bool( True ),
  estimator = cms.string( "hltESPMixedTripletStepChi2ChargeMeasurementEstimator16" ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetCkfBaseTrajectoryFilter_block" ) ),
  intermediateCleaning = cms.bool( True ),
  foundHitBonus = cms.double( 10.0 ),
  updator = cms.string( "hltESPKFUpdator" ),
  bestHitOnly = cms.bool( True ),
  seedAs5DHit = cms.bool( False )
)
process.HLTPSetHighPtTripletStepTrajectoryFilterPPOnAA = cms.PSet( 
  minimumNumberOfHits = cms.int32( 3 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 0 ),
  minPt = cms.double( 0.7 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 2.0 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutLoose" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  maxLostHits = cms.int32( 999 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetLowPtQuadStepTrajectoryBuilderPPOnAA = cms.PSet( 
  useSameTrajFilter = cms.bool( True ),
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  keepOriginalIfRebuildFails = cms.bool( False ),
  lostHitPenalty = cms.double( 30.0 ),
  lockHits = cms.bool( True ),
  requireSeedHitsInRebuild = cms.bool( True ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  maxDPhiForLooperReconstruction = cms.double( 2.0 ),
  maxPtForLooperReconstruction = cms.double( 0.7 ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetLowPtQuadStepTrajectoryFilterPPOnAA" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  minNrOfHitsForRebuild = cms.int32( 5 ),
  maxCand = cms.int32( 4 ),
  alwaysUseInvalidHits = cms.bool( True ),
  estimator = cms.string( "hltESPLowPtQuadStepChi2ChargeMeasurementEstimator9" ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetCkfBaseTrajectoryFilter_block" ) ),
  intermediateCleaning = cms.bool( True ),
  foundHitBonus = cms.double( 10.0 ),
  updator = cms.string( "hltESPKFUpdator" ),
  bestHitOnly = cms.bool( True ),
  seedAs5DHit = cms.bool( False )
)
process.HLTPSetPixelLessStepTrajectoryBuilderPPOnAA = cms.PSet( 
  useSameTrajFilter = cms.bool( True ),
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  keepOriginalIfRebuildFails = cms.bool( False ),
  lostHitPenalty = cms.double( 30.0 ),
  lockHits = cms.bool( True ),
  requireSeedHitsInRebuild = cms.bool( True ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  maxDPhiForLooperReconstruction = cms.double( 2.0 ),
  maxPtForLooperReconstruction = cms.double( 0.7 ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetPixelLessStepTrajectoryFilterPPOnAA" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  minNrOfHitsForRebuild = cms.int32( 4 ),
  maxCand = cms.int32( 2 ),
  alwaysUseInvalidHits = cms.bool( False ),
  estimator = cms.string( "hltESPPixelLessStepChi2ChargeMeasurementEstimator16" ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetCkfBaseTrajectoryFilter_block" ) ),
  intermediateCleaning = cms.bool( True ),
  foundHitBonus = cms.double( 10.0 ),
  updator = cms.string( "hltESPKFUpdator" ),
  bestHitOnly = cms.bool( True ),
  seedAs5DHit = cms.bool( False )
)
process.HLTSiStripClusterChargeCutLoose = cms.PSet(  value = cms.double( 1620.0 ) )
process.HLTPSetDetachedQuadStepTrajectoryBuilderPPOnAA = cms.PSet( 
  useSameTrajFilter = cms.bool( True ),
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  keepOriginalIfRebuildFails = cms.bool( False ),
  lostHitPenalty = cms.double( 30.0 ),
  lockHits = cms.bool( True ),
  requireSeedHitsInRebuild = cms.bool( True ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  maxDPhiForLooperReconstruction = cms.double( 2.0 ),
  maxPtForLooperReconstruction = cms.double( 0.7 ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetDetachedQuadStepTrajectoryFilterPPOnAA" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  minNrOfHitsForRebuild = cms.int32( 5 ),
  maxCand = cms.int32( 3 ),
  alwaysUseInvalidHits = cms.bool( True ),
  estimator = cms.string( "hltESPDetachedQuadStepChi2ChargeMeasurementEstimator9" ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetCkfBaseTrajectoryFilter_block" ) ),
  intermediateCleaning = cms.bool( True ),
  foundHitBonus = cms.double( 10.0 ),
  updator = cms.string( "hltESPKFUpdator" ),
  bestHitOnly = cms.bool( True ),
  seedAs5DHit = cms.bool( False )
)
process.HLTPSetInitialStepTrajectoryFilterShapePreSplittingPPOnAA = cms.PSet( 
  ComponentType = cms.string( "StripSubClusterShapeTrajectoryFilter" ),
  subclusterCutSN = cms.double( 12.0 ),
  trimMaxADC = cms.double( 30.0 ),
  seedCutMIPs = cms.double( 0.35 ),
  subclusterCutMIPs = cms.double( 0.45 ),
  subclusterWindow = cms.double( 0.7 ),
  maxNSat = cms.uint32( 3 ),
  trimMaxFracNeigh = cms.double( 0.25 ),
  maxTrimmedSizeDiffNeg = cms.double( 1.0 ),
  seedCutSN = cms.double( 7.0 ),
  layerMask = cms.PSet( 
    TOB = cms.bool( False ),
    TIB = cms.vuint32( 1, 2 ),
    TID = cms.vuint32( 1, 2 ),
    TEC = cms.bool( False )
  ),
  maxTrimmedSizeDiffPos = cms.double( 0.7 ),
  trimMaxFracTotal = cms.double( 0.15 )
)
process.HLTPSetInitialStepTrajectoryFilterPPOnAA = cms.PSet( 
  minimumNumberOfHits = cms.int32( 3 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 0 ),
  minPt = cms.double( 0.6 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 2.0 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutLoose" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  maxLostHits = cms.int32( 999 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetTobTecStepTrajectoryFilterPPOnAA = cms.PSet( 
  minimumNumberOfHits = cms.int32( 5 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 9999 ),
  minPt = cms.double( 2.0 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 2.0 ),
  seedPairPenalty = cms.int32( 1 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  maxLostHits = cms.int32( 0 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetPixelLessStepTrajectoryFilterPPOnAA = cms.PSet( 
  minimumNumberOfHits = cms.int32( 4 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 9999 ),
  minPt = cms.double( 2.0 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 2.0 ),
  seedPairPenalty = cms.int32( 1 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  maxLostHits = cms.int32( 0 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetPixelPairStepTrajectoryFilterPPOnAA = cms.PSet( 
  minimumNumberOfHits = cms.int32( 4 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 0 ),
  minPt = cms.double( 0.1 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 2.0 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutLoose" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  maxLostHits = cms.int32( 999 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetPixelPairStepTrajectoryFilterInOutPPOnAA = cms.PSet( 
  minimumNumberOfHits = cms.int32( 4 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 1 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 0 ),
  minPt = cms.double( 0.1 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 2.0 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutLoose" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  maxLostHits = cms.int32( 999 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetPixelPairStepTrajectoryBuilderPPOnAA = cms.PSet( 
  useSameTrajFilter = cms.bool( False ),
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  keepOriginalIfRebuildFails = cms.bool( False ),
  lostHitPenalty = cms.double( 30.0 ),
  lockHits = cms.bool( True ),
  requireSeedHitsInRebuild = cms.bool( True ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  maxDPhiForLooperReconstruction = cms.double( 2.0 ),
  maxPtForLooperReconstruction = cms.double( 0.7 ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetPixelPairStepTrajectoryFilterPPOnAA" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  minNrOfHitsForRebuild = cms.int32( 5 ),
  maxCand = cms.int32( 3 ),
  alwaysUseInvalidHits = cms.bool( True ),
  estimator = cms.string( "hltESPPixelPairStepChi2ChargeMeasurementEstimator9" ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetPixelPairStepTrajectoryFilterInOutPPOnAA" ) ),
  intermediateCleaning = cms.bool( True ),
  foundHitBonus = cms.double( 10.0 ),
  updator = cms.string( "hltESPKFUpdator" ),
  bestHitOnly = cms.bool( True ),
  seedAs5DHit = cms.bool( False )
)
process.HLTPSetInitialStepTrajectoryFilterBasePreSplittingForFullTrackingPPOnAA = cms.PSet( 
  minimumNumberOfHits = cms.int32( 4 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  maxCCCLostHits = cms.int32( 0 ),
  nSigmaMinPt = cms.double( 5.0 ),
  minPt = cms.double( 1.0 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 2.0 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutLoose" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  maxLostHits = cms.int32( 999 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetInitialStepTrajectoryBuilderPreSplittingForFullTrackingPPOnAA = cms.PSet( 
  useSameTrajFilter = cms.bool( True ),
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  keepOriginalIfRebuildFails = cms.bool( False ),
  lostHitPenalty = cms.double( 30.0 ),
  lockHits = cms.bool( True ),
  requireSeedHitsInRebuild = cms.bool( True ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  maxDPhiForLooperReconstruction = cms.double( 2.0 ),
  maxPtForLooperReconstruction = cms.double( 0.7 ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetInitialStepTrajectoryFilterPreSplittingForFullTrackingPPOnAA" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  minNrOfHitsForRebuild = cms.int32( 5 ),
  maxCand = cms.int32( 3 ),
  alwaysUseInvalidHits = cms.bool( True ),
  estimator = cms.string( "hltESPInitialStepChi2ChargeMeasurementEstimator30" ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetCkfBaseTrajectoryFilter_block" ) ),
  intermediateCleaning = cms.bool( True ),
  foundHitBonus = cms.double( 10.0 ),
  updator = cms.string( "hltESPKFUpdator" ),
  bestHitOnly = cms.bool( True ),
  seedAs5DHit = cms.bool( False )
)
process.HLTPSetInitialStepTrajectoryFilterPreSplittingForFullTrackingPPOnAA = cms.PSet( 
  ComponentType = cms.string( "CompositeTrajectoryFilter" ),
  filters = cms.VPSet( 
    cms.PSet(  refToPSet_ = cms.string( "HLTPSetInitialStepTrajectoryFilterBasePreSplittingForFullTrackingPPOnAA" )    ),
    cms.PSet(  refToPSet_ = cms.string( "HLTPSetInitialStepTrajectoryFilterShapePreSplittingPPOnAA" )    )
  )
)
process.HLTPSetInitialStepTrajectoryFilterForFullTrackingPPOnAA = cms.PSet( 
  minimumNumberOfHits = cms.int32( 4 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  maxCCCLostHits = cms.int32( 0 ),
  nSigmaMinPt = cms.double( 5.0 ),
  minPt = cms.double( 1.0 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 2.0 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutLoose" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  maxLostHits = cms.int32( 999 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetInitialStepTrajectoryBuilderForFullTrackingPPOnAA = cms.PSet( 
  useSameTrajFilter = cms.bool( True ),
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  keepOriginalIfRebuildFails = cms.bool( True ),
  lostHitPenalty = cms.double( 30.0 ),
  lockHits = cms.bool( True ),
  requireSeedHitsInRebuild = cms.bool( True ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  maxDPhiForLooperReconstruction = cms.double( 2.0 ),
  maxPtForLooperReconstruction = cms.double( 0.7 ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetInitialStepTrajectoryFilterForFullTrackingPPOnAA" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  minNrOfHitsForRebuild = cms.int32( 1 ),
  maxCand = cms.int32( 3 ),
  alwaysUseInvalidHits = cms.bool( True ),
  estimator = cms.string( "hltESPInitialStepChi2ChargeMeasurementEstimator30" ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetCkfBaseTrajectoryFilter_block" ) ),
  intermediateCleaning = cms.bool( True ),
  foundHitBonus = cms.double( 10.0 ),
  updator = cms.string( "hltESPKFUpdator" ),
  bestHitOnly = cms.bool( True ),
  seedAs5DHit = cms.bool( False )
)
process.HLTPSetLowPtQuadStepTrajectoryFilterForFullTrackingPPOnAA = cms.PSet( 
  minimumNumberOfHits = cms.int32( 3 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 0 ),
  minPt = cms.double( 1.0 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 2.0 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutLoose" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  maxLostHits = cms.int32( 999 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetLowPtQuadStepTrajectoryBuilderForFullTrackingPPOnAA = cms.PSet( 
  useSameTrajFilter = cms.bool( True ),
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  keepOriginalIfRebuildFails = cms.bool( False ),
  lostHitPenalty = cms.double( 30.0 ),
  lockHits = cms.bool( True ),
  requireSeedHitsInRebuild = cms.bool( True ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  maxDPhiForLooperReconstruction = cms.double( 2.0 ),
  maxPtForLooperReconstruction = cms.double( 0.7 ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetLowPtQuadStepTrajectoryFilterForFullTrackingPPOnAA" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  minNrOfHitsForRebuild = cms.int32( 5 ),
  maxCand = cms.int32( 4 ),
  alwaysUseInvalidHits = cms.bool( True ),
  estimator = cms.string( "hltESPLowPtQuadStepChi2ChargeMeasurementEstimator9" ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetCkfBaseTrajectoryFilter_block" ) ),
  intermediateCleaning = cms.bool( True ),
  foundHitBonus = cms.double( 10.0 ),
  updator = cms.string( "hltESPKFUpdator" ),
  bestHitOnly = cms.bool( True ),
  seedAs5DHit = cms.bool( False )
)
process.HLTPSetHighPtTripletStepTrajectoryFilterForFullTrackingPPOnAA = cms.PSet( 
  minimumNumberOfHits = cms.int32( 3 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 0 ),
  minPt = cms.double( 1.0 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 2.0 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutLoose" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  maxLostHits = cms.int32( 999 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetHighPtTripletStepTrajectoryBuilderForFullTrackingPPOnAA = cms.PSet( 
  useSameTrajFilter = cms.bool( True ),
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  keepOriginalIfRebuildFails = cms.bool( False ),
  lostHitPenalty = cms.double( 30.0 ),
  lockHits = cms.bool( True ),
  requireSeedHitsInRebuild = cms.bool( True ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  maxDPhiForLooperReconstruction = cms.double( 2.0 ),
  maxPtForLooperReconstruction = cms.double( 0.7 ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetHighPtTripletStepTrajectoryFilterForFullTrackingPPOnAA" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  minNrOfHitsForRebuild = cms.int32( 5 ),
  maxCand = cms.int32( 3 ),
  alwaysUseInvalidHits = cms.bool( True ),
  estimator = cms.string( "hltESPHighPtTripletStepChi2ChargeMeasurementEstimator30" ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetCkfBaseTrajectoryFilter_block" ) ),
  intermediateCleaning = cms.bool( True ),
  foundHitBonus = cms.double( 10.0 ),
  updator = cms.string( "hltESPKFUpdator" ),
  bestHitOnly = cms.bool( True ),
  seedAs5DHit = cms.bool( False )
)
process.HLTPSetLowPtTripletStepTrajectoryFilterForFullTrackingPPOnAA = cms.PSet( 
  minimumNumberOfHits = cms.int32( 3 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 0 ),
  minPt = cms.double( 2.8 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 2.0 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutLoose" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  maxLostHits = cms.int32( 999 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetLowPtTripletStepTrajectoryBuilderForFullTrackingPPOnAA = cms.PSet( 
  useSameTrajFilter = cms.bool( True ),
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  keepOriginalIfRebuildFails = cms.bool( False ),
  lostHitPenalty = cms.double( 30.0 ),
  lockHits = cms.bool( True ),
  requireSeedHitsInRebuild = cms.bool( True ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  maxDPhiForLooperReconstruction = cms.double( 2.0 ),
  maxPtForLooperReconstruction = cms.double( 0.7 ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetLowPtTripletStepTrajectoryFilterForFullTrackingPPOnAA" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  minNrOfHitsForRebuild = cms.int32( 5 ),
  maxCand = cms.int32( 4 ),
  alwaysUseInvalidHits = cms.bool( True ),
  estimator = cms.string( "hltESPLowPtTripletStepChi2ChargeMeasurementEstimator9" ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetCkfBaseTrajectoryFilter_block" ) ),
  intermediateCleaning = cms.bool( True ),
  foundHitBonus = cms.double( 10.0 ),
  updator = cms.string( "hltESPKFUpdator" ),
  bestHitOnly = cms.bool( True ),
  seedAs5DHit = cms.bool( False )
)
process.HLTPSetDetachedQuadStepTrajectoryFilterForFullTrackingPPOnAA = cms.PSet( 
  minimumNumberOfHits = cms.int32( 3 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 0 ),
  minPt = cms.double( 5.0 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 2.0 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutLoose" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  maxLostHits = cms.int32( 999 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetDetachedTripletStepTrajectoryFilterForFullTrackingPPOnAA = cms.PSet( 
  minimumNumberOfHits = cms.int32( 3 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 0 ),
  minPt = cms.double( 5.0 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 2.0 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutLoose" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  maxLostHits = cms.int32( 999 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetPixelPairStepTrajectoryFilterForFullTrackingPPOnAA = cms.PSet( 
  minimumNumberOfHits = cms.int32( 4 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 0 ),
  minPt = cms.double( 5.0 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 2.0 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutLoose" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  maxLostHits = cms.int32( 999 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetPixelPairStepTrajectoryBuilderForFullTrackingPPOnAA = cms.PSet( 
  useSameTrajFilter = cms.bool( False ),
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  keepOriginalIfRebuildFails = cms.bool( False ),
  lostHitPenalty = cms.double( 30.0 ),
  lockHits = cms.bool( True ),
  requireSeedHitsInRebuild = cms.bool( True ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  maxDPhiForLooperReconstruction = cms.double( 2.0 ),
  maxPtForLooperReconstruction = cms.double( 0.7 ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetPixelPairStepTrajectoryFilterForFullTrackingPPOnAA" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  minNrOfHitsForRebuild = cms.int32( 5 ),
  maxCand = cms.int32( 3 ),
  alwaysUseInvalidHits = cms.bool( True ),
  estimator = cms.string( "hltESPPixelPairStepChi2ChargeMeasurementEstimator9" ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetPixelPairStepTrajectoryFilterInOutForFullTrackingPPOnAA" ) ),
  intermediateCleaning = cms.bool( True ),
  foundHitBonus = cms.double( 10.0 ),
  updator = cms.string( "hltESPKFUpdator" ),
  bestHitOnly = cms.bool( True ),
  seedAs5DHit = cms.bool( False )
)
process.HLTPSetMixedTripletStepTrajectoryFilterForFullTrackingPPOnAA = cms.PSet( 
  minimumNumberOfHits = cms.int32( 3 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 9999 ),
  minPt = cms.double( 5.0 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 1.4 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  maxLostHits = cms.int32( 999 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetPixelLessStepTrajectoryFilterForFullTrackingPPOnAA = cms.PSet( 
  minimumNumberOfHits = cms.int32( 4 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 9999 ),
  minPt = cms.double( 5.0 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 2.0 ),
  seedPairPenalty = cms.int32( 1 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  maxLostHits = cms.int32( 0 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetPixelLessStepTrajectoryBuilderForFullTrackingPPOnAA = cms.PSet( 
  useSameTrajFilter = cms.bool( True ),
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  keepOriginalIfRebuildFails = cms.bool( False ),
  lostHitPenalty = cms.double( 30.0 ),
  lockHits = cms.bool( True ),
  requireSeedHitsInRebuild = cms.bool( True ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  maxDPhiForLooperReconstruction = cms.double( 2.0 ),
  maxPtForLooperReconstruction = cms.double( 0.7 ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetPixelLessStepTrajectoryFilterForFullTrackingPPOnAA" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  minNrOfHitsForRebuild = cms.int32( 4 ),
  maxCand = cms.int32( 2 ),
  alwaysUseInvalidHits = cms.bool( False ),
  estimator = cms.string( "hltESPPixelLessStepChi2ChargeMeasurementEstimator16" ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetCkfBaseTrajectoryFilter_block" ) ),
  intermediateCleaning = cms.bool( True ),
  foundHitBonus = cms.double( 10.0 ),
  updator = cms.string( "hltESPKFUpdator" ),
  bestHitOnly = cms.bool( True ),
  seedAs5DHit = cms.bool( False )
)
process.HLTPSetTobTecStepTrajectoryFilterForFullTrackingPPOnAA = cms.PSet( 
  minimumNumberOfHits = cms.int32( 5 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 9999 ),
  minPt = cms.double( 5.0 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 2.0 ),
  seedPairPenalty = cms.int32( 1 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  maxLostHits = cms.int32( 0 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetTobTecStepInOutTrajectoryFilterForFullTrackingPPOnAA = cms.PSet( 
  minimumNumberOfHits = cms.int32( 4 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 9999 ),
  minPt = cms.double( 5.0 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 2.0 ),
  seedPairPenalty = cms.int32( 1 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  maxLostHits = cms.int32( 0 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetTobTecStepTrajectoryBuilderForFullTrackingPPOnAA = cms.PSet( 
  useSameTrajFilter = cms.bool( False ),
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  keepOriginalIfRebuildFails = cms.bool( False ),
  lostHitPenalty = cms.double( 30.0 ),
  lockHits = cms.bool( True ),
  requireSeedHitsInRebuild = cms.bool( True ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  maxDPhiForLooperReconstruction = cms.double( 2.0 ),
  maxPtForLooperReconstruction = cms.double( 0.7 ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetTobTecStepTrajectoryFilterForFullTrackingPPOnAA" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  minNrOfHitsForRebuild = cms.int32( 4 ),
  maxCand = cms.int32( 2 ),
  alwaysUseInvalidHits = cms.bool( False ),
  estimator = cms.string( "hltESPTobTecStepChi2ChargeMeasurementEstimator16" ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetTobTecStepInOutTrajectoryFilterForFullTrackingPPOnAA" ) ),
  intermediateCleaning = cms.bool( True ),
  foundHitBonus = cms.double( 10.0 ),
  updator = cms.string( "hltESPKFUpdator" ),
  bestHitOnly = cms.bool( True ),
  seedAs5DHit = cms.bool( False )
)
process.HLTPSetJetCoreStepTrajectoryFilterForFullTrackingPPOnAA = cms.PSet( 
  minimumNumberOfHits = cms.int32( 4 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 9999 ),
  minPt = cms.double( 5.0 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 2.0 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  maxLostHits = cms.int32( 999 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetJetCoreStepTrajectoryBuilderForFullTrackingPPOnAA = cms.PSet( 
  useSameTrajFilter = cms.bool( True ),
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  keepOriginalIfRebuildFails = cms.bool( False ),
  lostHitPenalty = cms.double( 30.0 ),
  lockHits = cms.bool( True ),
  requireSeedHitsInRebuild = cms.bool( True ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  maxDPhiForLooperReconstruction = cms.double( 2.0 ),
  maxPtForLooperReconstruction = cms.double( 0.7 ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetJetCoreStepTrajectoryFilterForFullTrackingPPOnAA" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  minNrOfHitsForRebuild = cms.int32( 5 ),
  maxCand = cms.int32( 50 ),
  alwaysUseInvalidHits = cms.bool( True ),
  estimator = cms.string( "hltESPChi2MeasurementEstimator30" ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetCkfBaseTrajectoryFilter_block" ) ),
  intermediateCleaning = cms.bool( True ),
  foundHitBonus = cms.double( 10.0 ),
  updator = cms.string( "hltESPKFUpdator" ),
  bestHitOnly = cms.bool( True ),
  seedAs5DHit = cms.bool( False )
)
process.HLTPSetPixelPairStepTrajectoryFilterInOutForFullTrackingPPOnAA = cms.PSet( 
  minimumNumberOfHits = cms.int32( 4 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 1 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 0 ),
  minPt = cms.double( 5.0 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 2.0 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutLoose" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  maxLostHits = cms.int32( 999 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetMixedTripletStepTrajectoryBuilderForFullTrackingPPOnAA = cms.PSet( 
  useSameTrajFilter = cms.bool( True ),
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  keepOriginalIfRebuildFails = cms.bool( False ),
  lostHitPenalty = cms.double( 30.0 ),
  lockHits = cms.bool( True ),
  requireSeedHitsInRebuild = cms.bool( True ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  maxDPhiForLooperReconstruction = cms.double( 2.0 ),
  maxPtForLooperReconstruction = cms.double( 0.7 ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialForMixedStepOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetMixedTripletStepTrajectoryFilterForFullTrackingPPOnAA" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterialForMixedStep" ),
  minNrOfHitsForRebuild = cms.int32( 5 ),
  maxCand = cms.int32( 2 ),
  alwaysUseInvalidHits = cms.bool( True ),
  estimator = cms.string( "hltESPMixedTripletStepChi2ChargeMeasurementEstimator16" ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetCkfBaseTrajectoryFilter_block" ) ),
  intermediateCleaning = cms.bool( True ),
  foundHitBonus = cms.double( 10.0 ),
  updator = cms.string( "hltESPKFUpdator" ),
  bestHitOnly = cms.bool( True ),
  seedAs5DHit = cms.bool( False )
)
process.HLTPSetDetachedQuadStepTrajectoryBuilderForFullTrackingPPOnAA = cms.PSet( 
  useSameTrajFilter = cms.bool( True ),
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  keepOriginalIfRebuildFails = cms.bool( False ),
  lostHitPenalty = cms.double( 30.0 ),
  lockHits = cms.bool( True ),
  requireSeedHitsInRebuild = cms.bool( True ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  maxDPhiForLooperReconstruction = cms.double( 2.0 ),
  maxPtForLooperReconstruction = cms.double( 0.7 ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetDetachedQuadStepTrajectoryFilterForFullTrackingPPOnAA" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  minNrOfHitsForRebuild = cms.int32( 5 ),
  maxCand = cms.int32( 3 ),
  alwaysUseInvalidHits = cms.bool( True ),
  estimator = cms.string( "hltESPDetachedQuadStepChi2ChargeMeasurementEstimator9" ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetCkfBaseTrajectoryFilter_block" ) ),
  intermediateCleaning = cms.bool( True ),
  foundHitBonus = cms.double( 10.0 ),
  updator = cms.string( "hltESPKFUpdator" ),
  bestHitOnly = cms.bool( True ),
  seedAs5DHit = cms.bool( False )
)
process.HLTPSetDetachedTripletStepTrajectoryBuilderForFullTrackingPPOnAA = cms.PSet( 
  useSameTrajFilter = cms.bool( True ),
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  keepOriginalIfRebuildFails = cms.bool( False ),
  lostHitPenalty = cms.double( 30.0 ),
  lockHits = cms.bool( True ),
  requireSeedHitsInRebuild = cms.bool( True ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  maxDPhiForLooperReconstruction = cms.double( 2.0 ),
  maxPtForLooperReconstruction = cms.double( 0.7 ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetDetachedTripletStepTrajectoryFilterForFullTrackingPPOnAA" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  minNrOfHitsForRebuild = cms.int32( 5 ),
  maxCand = cms.int32( 3 ),
  alwaysUseInvalidHits = cms.bool( True ),
  estimator = cms.string( "hltESPDetachedTripletStepChi2ChargeMeasurementEstimator9" ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetCkfBaseTrajectoryFilter_block" ) ),
  intermediateCleaning = cms.bool( True ),
  foundHitBonus = cms.double( 10.0 ),
  updator = cms.string( "hltESPKFUpdator" ),
  bestHitOnly = cms.bool( True ),
  seedAs5DHit = cms.bool( False )
)
process.HLTPSetInitialStepTrajectoryFilterForDmesonPPOnAA = cms.PSet( 
  minimumNumberOfHits = cms.int32( 4 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  maxCCCLostHits = cms.int32( 0 ),
  nSigmaMinPt = cms.double( 5.0 ),
  minPt = cms.double( 3.0 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 2.0 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutLoose" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  maxLostHits = cms.int32( 999 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetInitialStepTrajectoryBuilderForDmesonPPOnAA = cms.PSet( 
  useSameTrajFilter = cms.bool( True ),
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  keepOriginalIfRebuildFails = cms.bool( True ),
  lostHitPenalty = cms.double( 30.0 ),
  lockHits = cms.bool( True ),
  requireSeedHitsInRebuild = cms.bool( True ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  maxDPhiForLooperReconstruction = cms.double( 2.0 ),
  maxPtForLooperReconstruction = cms.double( 0.7 ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetInitialStepTrajectoryFilterForDmesonPPOnAA" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  minNrOfHitsForRebuild = cms.int32( 1 ),
  maxCand = cms.int32( 3 ),
  alwaysUseInvalidHits = cms.bool( True ),
  estimator = cms.string( "hltESPInitialStepChi2ChargeMeasurementEstimator30" ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetCkfBaseTrajectoryFilter_block" ) ),
  intermediateCleaning = cms.bool( True ),
  foundHitBonus = cms.double( 10.0 ),
  updator = cms.string( "hltESPKFUpdator" ),
  bestHitOnly = cms.bool( True ),
  seedAs5DHit = cms.bool( False )
)
process.HLTPSetLowPtQuadStepTrajectoryFilterForDmesonPPOnAA = cms.PSet( 
  minimumNumberOfHits = cms.int32( 3 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 0 ),
  minPt = cms.double( 2.8 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 2.0 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutLoose" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  maxLostHits = cms.int32( 999 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetLowPtQuadStepTrajectoryBuilderForDmesonPPOnAA = cms.PSet( 
  useSameTrajFilter = cms.bool( True ),
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  keepOriginalIfRebuildFails = cms.bool( False ),
  lostHitPenalty = cms.double( 30.0 ),
  lockHits = cms.bool( True ),
  requireSeedHitsInRebuild = cms.bool( True ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  maxDPhiForLooperReconstruction = cms.double( 2.0 ),
  maxPtForLooperReconstruction = cms.double( 0.7 ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetLowPtQuadStepTrajectoryFilterForDmesonPPOnAA" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  minNrOfHitsForRebuild = cms.int32( 5 ),
  maxCand = cms.int32( 4 ),
  alwaysUseInvalidHits = cms.bool( True ),
  estimator = cms.string( "hltESPLowPtQuadStepChi2ChargeMeasurementEstimator9" ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetCkfBaseTrajectoryFilter_block" ) ),
  intermediateCleaning = cms.bool( True ),
  foundHitBonus = cms.double( 10.0 ),
  updator = cms.string( "hltESPKFUpdator" ),
  bestHitOnly = cms.bool( True ),
  seedAs5DHit = cms.bool( False )
)
process.HLTPSetHighPtTripletStepTrajectoryFilterForDmesonPPOnAA = cms.PSet( 
  minimumNumberOfHits = cms.int32( 3 ),
  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  seedExtension = cms.int32( 0 ),
  chargeSignificance = cms.double( -1.0 ),
  pixelSeedExtension = cms.bool( False ),
  strictSeedExtension = cms.bool( False ),
  nSigmaMinPt = cms.double( 5.0 ),
  maxCCCLostHits = cms.int32( 0 ),
  minPt = cms.double( 3.5 ),
  maxConsecLostHits = cms.int32( 1 ),
  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  constantValueForLostHitsFractionFilter = cms.double( 2.0 ),
  seedPairPenalty = cms.int32( 0 ),
  maxNumberOfHits = cms.int32( 100 ),
  minNumberOfHitsForLoopers = cms.int32( 13 ),
  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutLoose" ) ),
  minNumberOfHitsPerLoop = cms.int32( 4 ),
  minHitsMinPt = cms.int32( 3 ),
  maxLostHitsFraction = cms.double( 0.1 ),
  maxLostHits = cms.int32( 999 ),
  highEtaSwitch = cms.double( 5.0 ),
  minHitsAtHighEta = cms.int32( 5 )
)
process.HLTPSetHighPtTripletStepTrajectoryBuilderForDmesonPPOnAA = cms.PSet( 
  useSameTrajFilter = cms.bool( True ),
  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  keepOriginalIfRebuildFails = cms.bool( False ),
  lostHitPenalty = cms.double( 30.0 ),
  lockHits = cms.bool( True ),
  requireSeedHitsInRebuild = cms.bool( True ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  maxDPhiForLooperReconstruction = cms.double( 2.0 ),
  maxPtForLooperReconstruction = cms.double( 0.7 ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetHighPtTripletStepTrajectoryFilterForDmesonPPOnAA" ) ),
  propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  minNrOfHitsForRebuild = cms.int32( 5 ),
  maxCand = cms.int32( 3 ),
  alwaysUseInvalidHits = cms.bool( True ),
  estimator = cms.string( "hltESPHighPtTripletStepChi2ChargeMeasurementEstimator30" ),
  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetCkfBaseTrajectoryFilter_block" ) ),
  intermediateCleaning = cms.bool( True ),
  foundHitBonus = cms.double( 10.0 ),
  updator = cms.string( "hltESPKFUpdator" ),
  bestHitOnly = cms.bool( True ),
  seedAs5DHit = cms.bool( False )
)

process.ppsPixelTopologyESSource = cms.ESSource( "PPSPixelTopologyESSource",
    RunType = cms.string( "Run3" ),
    PitchSimY = cms.double( 0.15 ),
    PitchSimX = cms.double( 0.1 ),
    thickness = cms.double( 0.23 ),
    noOfPixelSimX = cms.int32( 160 ),
    noOfPixelSimY = cms.int32( 104 ),
    noOfPixels = cms.int32( 16640 ),
    simXWidth = cms.double( 16.6 ),
    simYWidth = cms.double( 16.2 ),
    deadEdgeWidth = cms.double( 0.2 ),
    activeEdgeSigma = cms.double( 0.02 ),
    physActiveEdgeDist = cms.double( 0.15 ),
    appendToDataLabel = cms.string( "" )
)
process.hltESSHcalSeverityLevel = cms.ESSource( "EmptyESSource",
    recordName = cms.string( "HcalSeverityLevelComputerRcd" ),
    iovIsRunNotTime = cms.bool( True ),
    firstValid = cms.vuint32( 1 )
)
process.hltESSEcalSeverityLevel = cms.ESSource( "EmptyESSource",
    recordName = cms.string( "EcalSeverityLevelAlgoRcd" ),
    iovIsRunNotTime = cms.bool( True ),
    firstValid = cms.vuint32( 1 )
)
process.hltESSBTagRecord = cms.ESSource( "EmptyESSource",
    recordName = cms.string( "JetTagComputerRecord" ),
    iovIsRunNotTime = cms.bool( True ),
    firstValid = cms.vuint32( 1 )
)
process.es_hardcode = cms.ESSource( "HcalHardcodeCalibrations",
    fromDDD = cms.untracked.bool( False ),
    toGet = cms.untracked.vstring( 'GainWidths' )
)
process.eegeom = cms.ESSource( "EmptyESSource",
    recordName = cms.string( "EcalMappingRcd" ),
    iovIsRunNotTime = cms.bool( True ),
    firstValid = cms.vuint32( 1 )
)
process.HepPDTESSource = cms.ESSource( "HepPDTESSource",
    pdtFileName = cms.FileInPath( "SimGeneral/HepPDTESSource/data/pythiaparticle.tbl" )
)
process.HcalTimeSlewEP = cms.ESSource( "HcalTimeSlewEP",
    appendToDataLabel = cms.string( "HBHE" ),
    timeSlewParametersM2 = cms.VPSet( 
      cms.PSet(  slope = cms.double( -3.178648 ),
        tmax = cms.double( 16.0 ),
        tzero = cms.double( 23.960177 )
      ),
      cms.PSet(  slope = cms.double( -1.5610227 ),
        tmax = cms.double( 10.0 ),
        tzero = cms.double( 11.977461 )
      ),
      cms.PSet(  slope = cms.double( -1.075824 ),
        tmax = cms.double( 6.25 ),
        tzero = cms.double( 9.109694 )
      )
    ),
    timeSlewParametersM3 = cms.VPSet( 
      cms.PSet(  tspar2 = cms.double( 0.0 ),
        tspar0 = cms.double( 12.2999 ),
        tspar1 = cms.double( -2.19142 ),
        tspar1_siPM = cms.double( 0.0 ),
        cap = cms.double( 6.0 ),
        tspar0_siPM = cms.double( 0.0 ),
        tspar2_siPM = cms.double( 0.0 )
      ),
      cms.PSet(  tspar2 = cms.double( 32.0 ),
        tspar0 = cms.double( 15.5 ),
        tspar1 = cms.double( -3.2 ),
        tspar1_siPM = cms.double( 0.0 ),
        cap = cms.double( 6.0 ),
        tspar0_siPM = cms.double( 0.0 ),
        tspar2_siPM = cms.double( 0.0 )
      ),
      cms.PSet(  tspar2 = cms.double( 0.0 ),
        tspar0 = cms.double( 12.2999 ),
        tspar1 = cms.double( -2.19142 ),
        tspar1_siPM = cms.double( 0.0 ),
        cap = cms.double( 6.0 ),
        tspar0_siPM = cms.double( 0.0 ),
        tspar2_siPM = cms.double( 0.0 )
      ),
      cms.PSet(  tspar2 = cms.double( 0.0 ),
        tspar0 = cms.double( 12.2999 ),
        tspar1 = cms.double( -2.19142 ),
        tspar1_siPM = cms.double( 0.0 ),
        cap = cms.double( 6.0 ),
        tspar0_siPM = cms.double( 0.0 ),
        tspar2_siPM = cms.double( 0.0 )
      )
    )
)
process.GlobalTag = cms.ESSource( "PoolDBESSource",
    DBParameters = cms.PSet( 
      connectionRetrialTimeOut = cms.untracked.int32( 60 ),
      idleConnectionCleanupPeriod = cms.untracked.int32( 10 ),
      enableReadOnlySessionOnUpdateConnection = cms.untracked.bool( False ),
      enablePoolAutomaticCleanUp = cms.untracked.bool( False ),
      messageLevel = cms.untracked.int32( 0 ),
      authenticationPath = cms.untracked.string( "." ),
      connectionRetrialPeriod = cms.untracked.int32( 10 ),
      connectionTimeOut = cms.untracked.int32( 0 ),
      enableConnectionSharing = cms.untracked.bool( True )
    ),
    connect = cms.string( "frontier://FrontierProd/CMS_CONDITIONS" ),
    globaltag = cms.string( "None" ),
    snapshotTime = cms.string( "" ),
    toGet = cms.VPSet( 
      cms.PSet(  record = cms.string( "BeamSpotOnlineLegacyObjectsRcd" ),
        refreshTime = cms.uint64( 2 )
      ),
      cms.PSet(  record = cms.string( "BeamSpotOnlineHLTObjectsRcd" ),
        refreshTime = cms.uint64( 2 )
      )
    ),
    DumpStat = cms.untracked.bool( False ),
    ReconnectEachRun = cms.untracked.bool( False ),
    RefreshAlways = cms.untracked.bool( False ),
    RefreshEachRun = cms.untracked.bool( False ),
    RefreshOpenIOVs = cms.untracked.bool( False ),
    pfnPostfix = cms.untracked.string( "None" )
)
process.GlobalParametersRcdSource = cms.ESSource( "EmptyESSource",
    recordName = cms.string( "L1TGlobalParametersRcd" ),
    iovIsRunNotTime = cms.bool( True ),
    firstValid = cms.vuint32( 1 )
)
process.CSCINdexerESSource = cms.ESSource( "EmptyESSource",
    recordName = cms.string( "CSCIndexerRecord" ),
    iovIsRunNotTime = cms.bool( True ),
    firstValid = cms.vuint32( 1 )
)
process.CSCChannelMapperESSource = cms.ESSource( "EmptyESSource",
    recordName = cms.string( "CSCChannelMapperRecord" ),
    iovIsRunNotTime = cms.bool( True ),
    firstValid = cms.vuint32( 1 )
)

process.trackerTopology = cms.ESProducer( "TrackerTopologyEP",
  appendToDataLabel = cms.string( "" )
)
process.sistripconn = cms.ESProducer( "SiStripConnectivity" )
process.siStripLorentzAngleDepESProducer = cms.ESProducer( "SiStripLorentzAngleDepESProducer",
  LatencyRecord = cms.PSet( 
    label = cms.untracked.string( "" ),
    record = cms.string( "SiStripLatencyRcd" )
  ),
  LorentzAnglePeakMode = cms.PSet( 
    label = cms.untracked.string( "peak" ),
    record = cms.string( "SiStripLorentzAngleRcd" )
  ),
  LorentzAngleDeconvMode = cms.PSet( 
    label = cms.untracked.string( "deconvolution" ),
    record = cms.string( "SiStripLorentzAngleRcd" )
  )
)
process.siStripBackPlaneCorrectionDepESProducer = cms.ESProducer( "SiStripBackPlaneCorrectionDepESProducer",
  LatencyRecord = cms.PSet( 
    label = cms.untracked.string( "" ),
    record = cms.string( "SiStripLatencyRcd" )
  ),
  BackPlaneCorrectionPeakMode = cms.PSet( 
    label = cms.untracked.string( "peak" ),
    record = cms.string( "SiStripBackPlaneCorrectionRcd" )
  ),
  BackPlaneCorrectionDeconvMode = cms.PSet( 
    label = cms.untracked.string( "deconvolution" ),
    record = cms.string( "SiStripBackPlaneCorrectionRcd" )
  )
)
process.siPixelTemplateDBObjectESProducer = cms.ESProducer( "SiPixelTemplateDBObjectESProducer" )
process.siPixelQualityESProducer = cms.ESProducer( "SiPixelQualityESProducer",
  siPixelQualityLabel = cms.string( "" )
)
process.preshowerDetIdAssociator = cms.ESProducer( "DetIdAssociatorESProducer",
  ComponentName = cms.string( "PreshowerDetIdAssociator" ),
  etaBinSize = cms.double( 0.1 ),
  nEta = cms.int32( 60 ),
  nPhi = cms.int32( 30 ),
  hcalRegion = cms.int32( 2 ),
  includeBadChambers = cms.bool( False ),
  includeGEM = cms.bool( False ),
  includeME0 = cms.bool( False )
)
process.navigationSchoolESProducer = cms.ESProducer( "NavigationSchoolESProducer",
  ComponentName = cms.string( "SimpleNavigationSchool" ),
  SimpleMagneticField = cms.string( "ParabolicMf" )
)
process.muonSeededTrajectoryCleanerBySharedHits = cms.ESProducer( "TrajectoryCleanerESProducer",
  ComponentName = cms.string( "muonSeededTrajectoryCleanerBySharedHits" ),
  ComponentType = cms.string( "TrajectoryCleanerBySharedHits" ),
  fractionShared = cms.double( 0.1 ),
  ValidHitBonus = cms.double( 1000.0 ),
  MissingHitPenalty = cms.double( 1.0 ),
  allowSharedFirstHit = cms.bool( True )
)
process.muonDetIdAssociator = cms.ESProducer( "DetIdAssociatorESProducer",
  ComponentName = cms.string( "MuonDetIdAssociator" ),
  etaBinSize = cms.double( 0.125 ),
  nEta = cms.int32( 48 ),
  nPhi = cms.int32( 48 ),
  hcalRegion = cms.int32( 2 ),
  includeBadChambers = cms.bool( False ),
  includeGEM = cms.bool( False ),
  includeME0 = cms.bool( False )
)
process.multipleScatteringParametrisationMakerESProducer = cms.ESProducer( "MultipleScatteringParametrisationMakerESProducer",
  appendToDataLabel = cms.string( "" )
)
process.hoDetIdAssociator = cms.ESProducer( "DetIdAssociatorESProducer",
  ComponentName = cms.string( "HODetIdAssociator" ),
  etaBinSize = cms.double( 0.087 ),
  nEta = cms.int32( 30 ),
  nPhi = cms.int32( 72 ),
  hcalRegion = cms.int32( 2 ),
  includeBadChambers = cms.bool( False ),
  includeGEM = cms.bool( False ),
  includeME0 = cms.bool( False )
)
process.hltTrackCleaner = cms.ESProducer( "TrackCleanerESProducer",
  ComponentName = cms.string( "hltTrackCleaner" ),
  appendToDataLabel = cms.string( "" )
)
process.hltPixelTracksCleanerBySharedHits = cms.ESProducer( "PixelTrackCleanerBySharedHitsESProducer",
  ComponentName = cms.string( "hltPixelTracksCleanerBySharedHits" ),
  useQuadrupletAlgo = cms.bool( False ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPTrajectorySmootherRK = cms.ESProducer( "KFTrajectorySmootherESProducer",
  ComponentName = cms.string( "hltESPTrajectorySmootherRK" ),
  Propagator = cms.string( "hltESPRungeKuttaTrackerPropagator" ),
  Updator = cms.string( "hltESPKFUpdator" ),
  Estimator = cms.string( "hltESPChi2MeasurementEstimator30" ),
  RecoGeometry = cms.string( "hltESPDummyDetLayerGeometry" ),
  errorRescaling = cms.double( 100.0 ),
  minHits = cms.int32( 3 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPTrajectoryFitterRK = cms.ESProducer( "KFTrajectoryFitterESProducer",
  ComponentName = cms.string( "hltESPTrajectoryFitterRK" ),
  Propagator = cms.string( "hltESPRungeKuttaTrackerPropagator" ),
  Updator = cms.string( "hltESPKFUpdator" ),
  Estimator = cms.string( "hltESPChi2MeasurementEstimator30" ),
  RecoGeometry = cms.string( "hltESPDummyDetLayerGeometry" ),
  minHits = cms.int32( 3 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPTrajectoryCleanerBySharedHits = cms.ESProducer( "TrajectoryCleanerESProducer",
  ComponentName = cms.string( "hltESPTrajectoryCleanerBySharedHits" ),
  ComponentType = cms.string( "TrajectoryCleanerBySharedHits" ),
  fractionShared = cms.double( 0.5 ),
  ValidHitBonus = cms.double( 100.0 ),
  MissingHitPenalty = cms.double( 0.0 ),
  allowSharedFirstHit = cms.bool( False )
)
process.hltESPTrackerRecoGeometryESProducer = cms.ESProducer( "TrackerRecoGeometryESProducer",
  usePhase2Stacks = cms.bool( False ),
  trackerGeometryLabel = cms.untracked.string( "" ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPTrackAlgoPriorityOrder = cms.ESProducer( "TrackAlgoPriorityOrderESProducer",
  ComponentName = cms.string( "hltESPTrackAlgoPriorityOrder" ),
  algoOrder = cms.vstring(  ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPTobTecStepTrajectoryCleanerBySharedHits = cms.ESProducer( "TrajectoryCleanerESProducer",
  ComponentName = cms.string( "hltESPTobTecStepTrajectoryCleanerBySharedHits" ),
  ComponentType = cms.string( "TrajectoryCleanerBySharedHits" ),
  fractionShared = cms.double( 0.09 ),
  ValidHitBonus = cms.double( 5.0 ),
  MissingHitPenalty = cms.double( 20.0 ),
  allowSharedFirstHit = cms.bool( True )
)
process.hltESPTobTecStepRKTrajectorySmootherForLoopers = cms.ESProducer( "KFTrajectorySmootherESProducer",
  ComponentName = cms.string( "hltESPTobTecStepRKSmootherForLoopers" ),
  Propagator = cms.string( "PropagatorWithMaterialForLoopers" ),
  Updator = cms.string( "hltESPKFUpdator" ),
  Estimator = cms.string( "hltESPChi2MeasurementEstimator30" ),
  RecoGeometry = cms.string( "hltESPDummyDetLayerGeometry" ),
  errorRescaling = cms.double( 10.0 ),
  minHits = cms.int32( 7 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPTobTecStepRKTrajectorySmoother = cms.ESProducer( "KFTrajectorySmootherESProducer",
  ComponentName = cms.string( "hltESPTobTecStepRKSmoother" ),
  Propagator = cms.string( "PropagatorWithMaterialParabolicMf" ),
  Updator = cms.string( "hltESPKFUpdator" ),
  Estimator = cms.string( "hltESPChi2MeasurementEstimator30" ),
  RecoGeometry = cms.string( "hltESPDummyDetLayerGeometry" ),
  errorRescaling = cms.double( 10.0 ),
  minHits = cms.int32( 7 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPTobTecStepRKTrajectoryFitterForLoopers = cms.ESProducer( "KFTrajectoryFitterESProducer",
  ComponentName = cms.string( "hltESPTobTecStepRKFitterForLoopers" ),
  Propagator = cms.string( "PropagatorWithMaterialForLoopers" ),
  Updator = cms.string( "hltESPKFUpdator" ),
  Estimator = cms.string( "hltESPChi2MeasurementEstimator30" ),
  RecoGeometry = cms.string( "hltESPDummyDetLayerGeometry" ),
  minHits = cms.int32( 7 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPTobTecStepRKTrajectoryFitter = cms.ESProducer( "KFTrajectoryFitterESProducer",
  ComponentName = cms.string( "hltESPTobTecStepRKFitter" ),
  Propagator = cms.string( "PropagatorWithMaterialParabolicMf" ),
  Updator = cms.string( "hltESPKFUpdator" ),
  Estimator = cms.string( "hltESPChi2MeasurementEstimator30" ),
  RecoGeometry = cms.string( "hltESPDummyDetLayerGeometry" ),
  minHits = cms.int32( 7 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPTobTecStepFlexibleKFFittingSmoother = cms.ESProducer( "FlexibleKFFittingSmootherESProducer",
  ComponentName = cms.string( "hltESPTobTecStepFlexibleKFFittingSmoother" ),
  standardFitter = cms.string( "hltESPTobTecStepFitterSmoother" ),
  looperFitter = cms.string( "hltESPTobTecStepFitterSmootherForLoopers" ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPTobTecStepFittingSmootherForLoopers = cms.ESProducer( "KFFittingSmootherESProducer",
  ComponentName = cms.string( "hltESPTobTecStepFitterSmootherForLoopers" ),
  Fitter = cms.string( "hltESPTobTecStepRKFitterForLoopers" ),
  Smoother = cms.string( "hltESPTobTecStepRKSmootherForLoopers" ),
  EstimateCut = cms.double( 30.0 ),
  MaxFractionOutliers = cms.double( 0.3 ),
  MaxNumberOfOutliers = cms.int32( 3 ),
  MinDof = cms.int32( 2 ),
  NoOutliersBeginEnd = cms.bool( False ),
  MinNumberOfHits = cms.int32( 7 ),
  MinNumberOfHitsHighEta = cms.int32( 5 ),
  HighEtaSwitch = cms.double( 5.0 ),
  RejectTracks = cms.bool( True ),
  BreakTrajWith2ConsecutiveMissing = cms.bool( False ),
  NoInvalidHitsBeginEnd = cms.bool( False ),
  LogPixelProbabilityCut = cms.double( -16.0 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPTobTecStepFittingSmoother = cms.ESProducer( "KFFittingSmootherESProducer",
  ComponentName = cms.string( "hltESPTobTecStepFitterSmoother" ),
  Fitter = cms.string( "hltESPTobTecStepRKFitter" ),
  Smoother = cms.string( "hltESPTobTecStepRKSmoother" ),
  EstimateCut = cms.double( 30.0 ),
  MaxFractionOutliers = cms.double( 0.3 ),
  MaxNumberOfOutliers = cms.int32( 3 ),
  MinDof = cms.int32( 2 ),
  NoOutliersBeginEnd = cms.bool( False ),
  MinNumberOfHits = cms.int32( 7 ),
  MinNumberOfHitsHighEta = cms.int32( 5 ),
  HighEtaSwitch = cms.double( 5.0 ),
  RejectTracks = cms.bool( True ),
  BreakTrajWith2ConsecutiveMissing = cms.bool( False ),
  NoInvalidHitsBeginEnd = cms.bool( False ),
  LogPixelProbabilityCut = cms.double( -16.0 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPTobTecStepClusterShapeHitFilter = cms.ESProducer( "ClusterShapeHitFilterESProducer",
  PixelShapeFile = cms.string( "RecoPixelVertexing/PixelLowPtUtilities/data/pixelShapePhase1_noL1.par" ),
  PixelShapeFileL1 = cms.string( "RecoPixelVertexing/PixelLowPtUtilities/data/pixelShapePhase1_loose.par" ),
  ComponentName = cms.string( "hltESPTobTecStepClusterShapeHitFilter" ),
  isPhase2 = cms.bool( False ),
  doPixelShapeCut = cms.bool( True ),
  doStripShapeCut = cms.bool( True ),
  clusterChargeCut = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutTight" ) ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPTobTecStepChi2ChargeMeasurementEstimator16 = cms.ESProducer( "Chi2ChargeMeasurementEstimatorESProducer",
  MaxChi2 = cms.double( 16.0 ),
  nSigma = cms.double( 3.0 ),
  MaxDisplacement = cms.double( 0.5 ),
  MaxSagitta = cms.double( 2.0 ),
  MinimalTolerance = cms.double( 0.5 ),
  MinPtForHitRecoveryInGluedDet = cms.double( 1000000.0 ),
  ComponentName = cms.string( "hltESPTobTecStepChi2ChargeMeasurementEstimator16" ),
  pTChargeCutThreshold = cms.double( -1.0 ),
  clusterChargeCut = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutTight" ) ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPTTRHBuilderWithoutAngle4PixelTriplets = cms.ESProducer( "TkTransientTrackingRecHitBuilderESProducer",
  ComponentName = cms.string( "hltESPTTRHBuilderWithoutAngle4PixelTriplets" ),
  ComputeCoarseLocalPositionFromDisk = cms.bool( False ),
  StripCPE = cms.string( "Fake" ),
  PixelCPE = cms.string( "hltESPPixelCPEGeneric" ),
  Matcher = cms.string( "StandardMatcher" ),
  Phase2StripCPE = cms.string( "" ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPTTRHBuilderPixelOnly = cms.ESProducer( "TkTransientTrackingRecHitBuilderESProducer",
  ComponentName = cms.string( "hltESPTTRHBuilderPixelOnly" ),
  ComputeCoarseLocalPositionFromDisk = cms.bool( False ),
  StripCPE = cms.string( "Fake" ),
  PixelCPE = cms.string( "hltESPPixelCPEGeneric" ),
  Matcher = cms.string( "StandardMatcher" ),
  Phase2StripCPE = cms.string( "" ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPTTRHBuilderAngleAndTemplate = cms.ESProducer( "TkTransientTrackingRecHitBuilderESProducer",
  ComponentName = cms.string( "hltESPTTRHBuilderAngleAndTemplate" ),
  ComputeCoarseLocalPositionFromDisk = cms.bool( False ),
  StripCPE = cms.string( "hltESPStripCPEfromTrackAngle" ),
  PixelCPE = cms.string( "hltESPPixelCPETemplateReco" ),
  Matcher = cms.string( "StandardMatcher" ),
  Phase2StripCPE = cms.string( "" ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPTTRHBWithTrackAngle = cms.ESProducer( "TkTransientTrackingRecHitBuilderESProducer",
  ComponentName = cms.string( "hltESPTTRHBWithTrackAngle" ),
  ComputeCoarseLocalPositionFromDisk = cms.bool( False ),
  StripCPE = cms.string( "hltESPStripCPEfromTrackAngle" ),
  PixelCPE = cms.string( "hltESPPixelCPEGeneric" ),
  Matcher = cms.string( "StandardMatcher" ),
  Phase2StripCPE = cms.string( "" ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPStripCPEfromTrackAngle = cms.ESProducer( "StripCPEESProducer",
  ComponentName = cms.string( "hltESPStripCPEfromTrackAngle" ),
  ComponentType = cms.string( "StripCPEfromTrackAngle" ),
  parameters = cms.PSet( 
    mTIB_P1 = cms.double( 0.202 ),
    maxChgOneMIP = cms.double( 6000.0 ),
    mTEC_P0 = cms.double( -1.885 ),
    mTOB_P1 = cms.double( 0.253 ),
    mTEC_P1 = cms.double( 0.471 ),
    mLC_P2 = cms.double( 0.3 ),
    mLC_P1 = cms.double( 0.618 ),
    mTOB_P0 = cms.double( -1.026 ),
    mLC_P0 = cms.double( -0.326 ),
    useLegacyError = cms.bool( False ),
    mTIB_P0 = cms.double( -0.742 ),
    mTID_P1 = cms.double( 0.433 ),
    mTID_P0 = cms.double( -1.427 )
  )
)
process.hltESPSteppingHelixPropagatorOpposite = cms.ESProducer( "SteppingHelixPropagatorESProducer",
  ComponentName = cms.string( "hltESPSteppingHelixPropagatorOpposite" ),
  NoErrorPropagation = cms.bool( False ),
  PropagationDirection = cms.string( "oppositeToMomentum" ),
  useTuningForL2Speed = cms.bool( False ),
  useIsYokeFlag = cms.bool( True ),
  endcapShiftInZNeg = cms.double( 0.0 ),
  SetVBFPointer = cms.bool( False ),
  AssumeNoMaterial = cms.bool( False ),
  endcapShiftInZPos = cms.double( 0.0 ),
  useInTeslaFromMagField = cms.bool( False ),
  VBFName = cms.string( "VolumeBasedMagneticField" ),
  useEndcapShiftsInZ = cms.bool( False ),
  sendLogWarning = cms.bool( False ),
  useMatVolumes = cms.bool( True ),
  debug = cms.bool( False ),
  ApplyRadX0Correction = cms.bool( True ),
  useMagVolumes = cms.bool( True ),
  returnTangentPlane = cms.bool( True )
)
process.hltESPSteppingHelixPropagatorAlong = cms.ESProducer( "SteppingHelixPropagatorESProducer",
  ComponentName = cms.string( "hltESPSteppingHelixPropagatorAlong" ),
  NoErrorPropagation = cms.bool( False ),
  PropagationDirection = cms.string( "alongMomentum" ),
  useTuningForL2Speed = cms.bool( False ),
  useIsYokeFlag = cms.bool( True ),
  endcapShiftInZNeg = cms.double( 0.0 ),
  SetVBFPointer = cms.bool( False ),
  AssumeNoMaterial = cms.bool( False ),
  endcapShiftInZPos = cms.double( 0.0 ),
  useInTeslaFromMagField = cms.bool( False ),
  VBFName = cms.string( "VolumeBasedMagneticField" ),
  useEndcapShiftsInZ = cms.bool( False ),
  sendLogWarning = cms.bool( False ),
  useMatVolumes = cms.bool( True ),
  debug = cms.bool( False ),
  ApplyRadX0Correction = cms.bool( True ),
  useMagVolumes = cms.bool( True ),
  returnTangentPlane = cms.bool( True )
)
process.hltESPSoftLeptonByDistance = cms.ESProducer( "LeptonTaggerByDistanceESProducer",
  distance = cms.double( 0.5 )
)
process.hltESPSmartPropagatorAnyOpposite = cms.ESProducer( "SmartPropagatorESProducer",
  ComponentName = cms.string( "hltESPSmartPropagatorAnyOpposite" ),
  TrackerPropagator = cms.string( "PropagatorWithMaterialOpposite" ),
  MuonPropagator = cms.string( "SteppingHelixPropagatorAny" ),
  PropagationDirection = cms.string( "oppositeToMomentum" ),
  Epsilon = cms.double( 5.0 )
)
process.hltESPSmartPropagatorAny = cms.ESProducer( "SmartPropagatorESProducer",
  ComponentName = cms.string( "hltESPSmartPropagatorAny" ),
  TrackerPropagator = cms.string( "PropagatorWithMaterial" ),
  MuonPropagator = cms.string( "SteppingHelixPropagatorAny" ),
  PropagationDirection = cms.string( "alongMomentum" ),
  Epsilon = cms.double( 5.0 )
)
process.hltESPSmartPropagator = cms.ESProducer( "SmartPropagatorESProducer",
  ComponentName = cms.string( "hltESPSmartPropagator" ),
  TrackerPropagator = cms.string( "PropagatorWithMaterial" ),
  MuonPropagator = cms.string( "hltESPSteppingHelixPropagatorAlong" ),
  PropagationDirection = cms.string( "alongMomentum" ),
  Epsilon = cms.double( 5.0 )
)
process.hltESPRungeKuttaTrackerPropagator = cms.ESProducer( "PropagatorWithMaterialESProducer",
  SimpleMagneticField = cms.string( "" ),
  MaxDPhi = cms.double( 1.6 ),
  ComponentName = cms.string( "hltESPRungeKuttaTrackerPropagator" ),
  Mass = cms.double( 0.105 ),
  PropagationDirection = cms.string( "alongMomentum" ),
  useRungeKutta = cms.bool( True ),
  ptMin = cms.double( -1.0 )
)
process.hltESPRKTrajectorySmoother = cms.ESProducer( "KFTrajectorySmootherESProducer",
  ComponentName = cms.string( "hltESPRKTrajectorySmoother" ),
  Propagator = cms.string( "hltESPRungeKuttaTrackerPropagator" ),
  Updator = cms.string( "hltESPKFUpdator" ),
  Estimator = cms.string( "hltESPChi2MeasurementEstimator30" ),
  RecoGeometry = cms.string( "hltESPGlobalDetLayerGeometry" ),
  errorRescaling = cms.double( 100.0 ),
  minHits = cms.int32( 3 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPRKTrajectoryFitter = cms.ESProducer( "KFTrajectoryFitterESProducer",
  ComponentName = cms.string( "hltESPRKTrajectoryFitter" ),
  Propagator = cms.string( "hltESPRungeKuttaTrackerPropagator" ),
  Updator = cms.string( "hltESPKFUpdator" ),
  Estimator = cms.string( "hltESPChi2MeasurementEstimator30" ),
  RecoGeometry = cms.string( "hltESPGlobalDetLayerGeometry" ),
  minHits = cms.int32( 3 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPPixelPairTrajectoryCleanerBySharedHits = cms.ESProducer( "TrajectoryCleanerESProducer",
  ComponentName = cms.string( "hltESPPixelPairTrajectoryCleanerBySharedHits" ),
  ComponentType = cms.string( "TrajectoryCleanerBySharedHits" ),
  fractionShared = cms.double( 0.19 ),
  ValidHitBonus = cms.double( 5.0 ),
  MissingHitPenalty = cms.double( 20.0 ),
  allowSharedFirstHit = cms.bool( True )
)
process.hltESPPixelPairStepChi2MeasurementEstimator25 = cms.ESProducer( "Chi2MeasurementEstimatorESProducer",
  MaxChi2 = cms.double( 25.0 ),
  nSigma = cms.double( 3.0 ),
  MaxDisplacement = cms.double( 100.0 ),
  MaxSagitta = cms.double( -1.0 ),
  MinimalTolerance = cms.double( 10.0 ),
  MinPtForHitRecoveryInGluedDet = cms.double( 1000000.0 ),
  ComponentName = cms.string( "hltESPPixelPairStepChi2MeasurementEstimator25" ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPPixelPairStepChi2ChargeMeasurementEstimator9 = cms.ESProducer( "Chi2ChargeMeasurementEstimatorESProducer",
  MaxChi2 = cms.double( 9.0 ),
  nSigma = cms.double( 3.0 ),
  MaxDisplacement = cms.double( 0.5 ),
  MaxSagitta = cms.double( 2.0 ),
  MinimalTolerance = cms.double( 0.5 ),
  MinPtForHitRecoveryInGluedDet = cms.double( 1.0E12 ),
  ComponentName = cms.string( "hltESPPixelPairStepChi2ChargeMeasurementEstimator9" ),
  pTChargeCutThreshold = cms.double( 15.0 ),
  clusterChargeCut = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutLoose" ) ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPPixelLessStepTrajectoryCleanerBySharedHits = cms.ESProducer( "TrajectoryCleanerESProducer",
  ComponentName = cms.string( "hltESPPixelLessStepTrajectoryCleanerBySharedHits" ),
  ComponentType = cms.string( "TrajectoryCleanerBySharedHits" ),
  fractionShared = cms.double( 0.11 ),
  ValidHitBonus = cms.double( 5.0 ),
  MissingHitPenalty = cms.double( 20.0 ),
  allowSharedFirstHit = cms.bool( True )
)
process.hltESPPixelLessStepClusterShapeHitFilter = cms.ESProducer( "ClusterShapeHitFilterESProducer",
  PixelShapeFile = cms.string( "RecoPixelVertexing/PixelLowPtUtilities/data/pixelShapePhase1_noL1.par" ),
  PixelShapeFileL1 = cms.string( "RecoPixelVertexing/PixelLowPtUtilities/data/pixelShapePhase1_loose.par" ),
  ComponentName = cms.string( "hltESPPixelLessStepClusterShapeHitFilter" ),
  isPhase2 = cms.bool( False ),
  doPixelShapeCut = cms.bool( True ),
  doStripShapeCut = cms.bool( True ),
  clusterChargeCut = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutTight" ) ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPPixelLessStepChi2ChargeMeasurementEstimator16 = cms.ESProducer( "Chi2ChargeMeasurementEstimatorESProducer",
  MaxChi2 = cms.double( 16.0 ),
  nSigma = cms.double( 3.0 ),
  MaxDisplacement = cms.double( 0.5 ),
  MaxSagitta = cms.double( 2.0 ),
  MinimalTolerance = cms.double( 0.5 ),
  MinPtForHitRecoveryInGluedDet = cms.double( 1000000.0 ),
  ComponentName = cms.string( "hltESPPixelLessStepChi2ChargeMeasurementEstimator16" ),
  pTChargeCutThreshold = cms.double( -1.0 ),
  clusterChargeCut = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutTight" ) ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPPixelCPETemplateReco = cms.ESProducer( "PixelCPETemplateRecoESProducer",
  LoadTemplatesFromDB = cms.bool( True ),
  Alpha2Order = cms.bool( True ),
  ClusterProbComputationFlag = cms.int32( 0 ),
  useLAWidthFromDB = cms.bool( True ),
  lAOffset = cms.double( 0.0 ),
  lAWidthBPix = cms.double( 0.0 ),
  lAWidthFPix = cms.double( 0.0 ),
  doLorentzFromAlignment = cms.bool( False ),
  useLAFromDB = cms.bool( True ),
  barrelTemplateID = cms.int32( 0 ),
  forwardTemplateID = cms.int32( 0 ),
  directoryWithTemplates = cms.int32( 0 ),
  speed = cms.int32( -2 ),
  UseClusterSplitter = cms.bool( False ),
  ComponentName = cms.string( "hltESPPixelCPETemplateReco" ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPPixelCPEGeneric = cms.ESProducer( "PixelCPEGenericESProducer",
  LoadTemplatesFromDB = cms.bool( True ),
  Alpha2Order = cms.bool( True ),
  ClusterProbComputationFlag = cms.int32( 0 ),
  useLAWidthFromDB = cms.bool( False ),
  lAOffset = cms.double( 0.0 ),
  lAWidthBPix = cms.double( 0.0 ),
  lAWidthFPix = cms.double( 0.0 ),
  doLorentzFromAlignment = cms.bool( False ),
  useLAFromDB = cms.bool( True ),
  xerr_barrel_l1 = cms.vdouble( 0.00115, 0.0012, 8.8E-4 ),
  yerr_barrel_l1 = cms.vdouble( 0.00375, 0.0023, 0.0025, 0.0025, 0.0023, 0.0023, 0.0021, 0.0021, 0.0024 ),
  xerr_barrel_ln = cms.vdouble( 0.00115, 0.0012, 8.8E-4 ),
  yerr_barrel_ln = cms.vdouble( 0.00375, 0.0023, 0.0025, 0.0025, 0.0023, 0.0023, 0.0021, 0.0021, 0.0024 ),
  xerr_endcap = cms.vdouble( 0.002, 0.002 ),
  yerr_endcap = cms.vdouble( 0.0021 ),
  xerr_barrel_l1_def = cms.double( 0.0103 ),
  yerr_barrel_l1_def = cms.double( 0.0021 ),
  xerr_barrel_ln_def = cms.double( 0.0103 ),
  yerr_barrel_ln_def = cms.double( 0.0021 ),
  xerr_endcap_def = cms.double( 0.002 ),
  yerr_endcap_def = cms.double( 7.5E-4 ),
  eff_charge_cut_highX = cms.double( 1.0 ),
  eff_charge_cut_highY = cms.double( 1.0 ),
  eff_charge_cut_lowX = cms.double( 0.0 ),
  eff_charge_cut_lowY = cms.double( 0.0 ),
  size_cutX = cms.double( 3.0 ),
  size_cutY = cms.double( 3.0 ),
  EdgeClusterErrorX = cms.double( 50.0 ),
  EdgeClusterErrorY = cms.double( 85.0 ),
  inflate_errors = cms.bool( False ),
  inflate_all_errors_no_trk_angle = cms.bool( False ),
  NoTemplateErrorsWhenNoTrkAngles = cms.bool( False ),
  UseErrorsFromTemplates = cms.bool( True ),
  TruncatePixelCharge = cms.bool( True ),
  IrradiationBiasCorrection = cms.bool( True ),
  DoCosmics = cms.bool( False ),
  Upgrade = cms.bool( False ),
  SmallPitch = cms.bool( False ),
  ComponentName = cms.string( "hltESPPixelCPEGeneric" ),
  MagneticFieldRecord = cms.ESInputTag( "","" ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPMuonTransientTrackingRecHitBuilder = cms.ESProducer( "MuonTransientTrackingRecHitBuilderESProducer",
  ComponentName = cms.string( "hltESPMuonTransientTrackingRecHitBuilder" )
)
process.hltESPMuonDetLayerGeometryESProducer = cms.ESProducer( "MuonDetLayerGeometryESProducer" )
process.hltESPMixedTripletStepTrajectoryCleanerBySharedHits = cms.ESProducer( "TrajectoryCleanerESProducer",
  ComponentName = cms.string( "hltESPMixedTripletStepTrajectoryCleanerBySharedHits" ),
  ComponentType = cms.string( "TrajectoryCleanerBySharedHits" ),
  fractionShared = cms.double( 0.11 ),
  ValidHitBonus = cms.double( 5.0 ),
  MissingHitPenalty = cms.double( 20.0 ),
  allowSharedFirstHit = cms.bool( True )
)
process.hltESPMixedTripletStepChi2ChargeMeasurementEstimator16 = cms.ESProducer( "Chi2ChargeMeasurementEstimatorESProducer",
  MaxChi2 = cms.double( 16.0 ),
  nSigma = cms.double( 3.0 ),
  MaxDisplacement = cms.double( 0.5 ),
  MaxSagitta = cms.double( 2.0 ),
  MinimalTolerance = cms.double( 0.5 ),
  MinPtForHitRecoveryInGluedDet = cms.double( 1000000.0 ),
  ComponentName = cms.string( "hltESPMixedTripletStepChi2ChargeMeasurementEstimator16" ),
  pTChargeCutThreshold = cms.double( -1.0 ),
  clusterChargeCut = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutTight" ) ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPMixedStepTrajectoryCleanerBySharedHits = cms.ESProducer( "TrajectoryCleanerESProducer",
  ComponentName = cms.string( "hltESPMixedStepTrajectoryCleanerBySharedHits" ),
  ComponentType = cms.string( "TrajectoryCleanerBySharedHits" ),
  fractionShared = cms.double( 0.11 ),
  ValidHitBonus = cms.double( 5.0 ),
  MissingHitPenalty = cms.double( 20.0 ),
  allowSharedFirstHit = cms.bool( True )
)
process.hltESPMixedStepClusterShapeHitFilter = cms.ESProducer( "ClusterShapeHitFilterESProducer",
  PixelShapeFile = cms.string( "RecoPixelVertexing/PixelLowPtUtilities/data/pixelShapePhase1_noL1.par" ),
  PixelShapeFileL1 = cms.string( "RecoPixelVertexing/PixelLowPtUtilities/data/pixelShapePhase1_loose.par" ),
  ComponentName = cms.string( "hltESPMixedStepClusterShapeHitFilter" ),
  isPhase2 = cms.bool( False ),
  doPixelShapeCut = cms.bool( True ),
  doStripShapeCut = cms.bool( True ),
  clusterChargeCut = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutTight" ) ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPMeasurementTracker = cms.ESProducer( "MeasurementTrackerESProducer",
  ComponentName = cms.string( "hltESPMeasurementTracker" ),
  PixelCPE = cms.string( "hltESPPixelCPEGeneric" ),
  StripCPE = cms.string( "hltESPStripCPEfromTrackAngle" ),
  HitMatcher = cms.string( "StandardMatcher" ),
  Phase2StripCPE = cms.string( "" ),
  SiStripQualityLabel = cms.string( "" ),
  UseStripModuleQualityDB = cms.bool( True ),
  DebugStripModuleQualityDB = cms.untracked.bool( False ),
  UseStripAPVFiberQualityDB = cms.bool( True ),
  DebugStripAPVFiberQualityDB = cms.untracked.bool( False ),
  MaskBadAPVFibers = cms.bool( True ),
  UseStripStripQualityDB = cms.bool( True ),
  DebugStripStripQualityDB = cms.untracked.bool( False ),
  badStripCuts = cms.PSet( 
    TOB = cms.PSet( 
      maxBad = cms.uint32( 4 ),
      maxConsecutiveBad = cms.uint32( 2 )
    ),
    TIB = cms.PSet( 
      maxBad = cms.uint32( 4 ),
      maxConsecutiveBad = cms.uint32( 2 )
    ),
    TID = cms.PSet( 
      maxBad = cms.uint32( 4 ),
      maxConsecutiveBad = cms.uint32( 2 )
    ),
    TEC = cms.PSet( 
      maxBad = cms.uint32( 4 ),
      maxConsecutiveBad = cms.uint32( 2 )
    )
  ),
  UsePixelModuleQualityDB = cms.bool( True ),
  DebugPixelModuleQualityDB = cms.untracked.bool( False ),
  UsePixelROCQualityDB = cms.bool( True ),
  DebugPixelROCQualityDB = cms.untracked.bool( False ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPLowPtTripletStepTrajectoryCleanerBySharedHits = cms.ESProducer( "TrajectoryCleanerESProducer",
  ComponentName = cms.string( "hltESPLowPtTripletStepTrajectoryCleanerBySharedHits" ),
  ComponentType = cms.string( "TrajectoryCleanerBySharedHits" ),
  fractionShared = cms.double( 0.16 ),
  ValidHitBonus = cms.double( 5.0 ),
  MissingHitPenalty = cms.double( 20.0 ),
  allowSharedFirstHit = cms.bool( True )
)
process.hltESPLowPtTripletStepChi2ChargeMeasurementEstimator9 = cms.ESProducer( "Chi2ChargeMeasurementEstimatorESProducer",
  MaxChi2 = cms.double( 9.0 ),
  nSigma = cms.double( 3.0 ),
  MaxDisplacement = cms.double( 0.5 ),
  MaxSagitta = cms.double( 2.0 ),
  MinimalTolerance = cms.double( 0.5 ),
  MinPtForHitRecoveryInGluedDet = cms.double( 1000000.0 ),
  ComponentName = cms.string( "hltESPLowPtTripletStepChi2ChargeMeasurementEstimator9" ),
  pTChargeCutThreshold = cms.double( -1.0 ),
  clusterChargeCut = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutTight" ) ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPLowPtStepTrajectoryCleanerBySharedHits = cms.ESProducer( "TrajectoryCleanerESProducer",
  ComponentName = cms.string( "hltESPLowPtStepTrajectoryCleanerBySharedHits" ),
  ComponentType = cms.string( "TrajectoryCleanerBySharedHits" ),
  fractionShared = cms.double( 0.16 ),
  ValidHitBonus = cms.double( 5.0 ),
  MissingHitPenalty = cms.double( 20.0 ),
  allowSharedFirstHit = cms.bool( True )
)
process.hltESPLowPtQuadStepTrajectoryCleanerBySharedHits = cms.ESProducer( "TrajectoryCleanerESProducer",
  ComponentName = cms.string( "hltESPLowPtQuadStepTrajectoryCleanerBySharedHits" ),
  ComponentType = cms.string( "TrajectoryCleanerBySharedHits" ),
  fractionShared = cms.double( 0.16 ),
  ValidHitBonus = cms.double( 5.0 ),
  MissingHitPenalty = cms.double( 20.0 ),
  allowSharedFirstHit = cms.bool( True )
)
process.hltESPLowPtQuadStepChi2ChargeMeasurementEstimator9 = cms.ESProducer( "Chi2ChargeMeasurementEstimatorESProducer",
  MaxChi2 = cms.double( 9.0 ),
  nSigma = cms.double( 3.0 ),
  MaxDisplacement = cms.double( 0.5 ),
  MaxSagitta = cms.double( 2.0 ),
  MinimalTolerance = cms.double( 0.5 ),
  MinPtForHitRecoveryInGluedDet = cms.double( 1000000.0 ),
  ComponentName = cms.string( "hltESPLowPtQuadStepChi2ChargeMeasurementEstimator9" ),
  pTChargeCutThreshold = cms.double( -1.0 ),
  clusterChargeCut = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutTight" ) ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPL3MuKFTrajectoryFitter = cms.ESProducer( "KFTrajectoryFitterESProducer",
  ComponentName = cms.string( "hltESPL3MuKFTrajectoryFitter" ),
  Propagator = cms.string( "hltESPSmartPropagatorAny" ),
  Updator = cms.string( "hltESPKFUpdator" ),
  Estimator = cms.string( "hltESPChi2MeasurementEstimator30" ),
  RecoGeometry = cms.string( "hltESPDummyDetLayerGeometry" ),
  minHits = cms.int32( 3 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPKullbackLeiblerDistance5D = cms.ESProducer( "DistanceBetweenComponentsESProducer5D",
  ComponentName = cms.string( "hltESPKullbackLeiblerDistance5D" ),
  DistanceMeasure = cms.string( "KullbackLeibler" )
)
process.hltESPKFUpdator = cms.ESProducer( "KFUpdatorESProducer",
  ComponentName = cms.string( "hltESPKFUpdator" )
)
process.hltESPKFTrajectorySmootherForMuonTrackLoader = cms.ESProducer( "KFTrajectorySmootherESProducer",
  ComponentName = cms.string( "hltESPKFTrajectorySmootherForMuonTrackLoader" ),
  Propagator = cms.string( "hltESPSmartPropagatorAnyOpposite" ),
  Updator = cms.string( "hltESPKFUpdator" ),
  Estimator = cms.string( "hltESPChi2MeasurementEstimator30" ),
  RecoGeometry = cms.string( "hltESPDummyDetLayerGeometry" ),
  errorRescaling = cms.double( 10.0 ),
  minHits = cms.int32( 3 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPKFTrajectorySmootherForLoopers = cms.ESProducer( "KFTrajectorySmootherESProducer",
  ComponentName = cms.string( "hltESPKFTrajectorySmootherForLoopers" ),
  Propagator = cms.string( "PropagatorWithMaterialForLoopers" ),
  Updator = cms.string( "hltESPKFUpdator" ),
  Estimator = cms.string( "hltESPChi2MeasurementEstimator30" ),
  RecoGeometry = cms.string( "hltESPGlobalDetLayerGeometry" ),
  errorRescaling = cms.double( 10.0 ),
  minHits = cms.int32( 3 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPKFTrajectorySmootherForL2Muon = cms.ESProducer( "KFTrajectorySmootherESProducer",
  ComponentName = cms.string( "hltESPKFTrajectorySmootherForL2Muon" ),
  Propagator = cms.string( "hltESPFastSteppingHelixPropagatorOpposite" ),
  Updator = cms.string( "hltESPKFUpdator" ),
  Estimator = cms.string( "hltESPChi2MeasurementEstimator30" ),
  RecoGeometry = cms.string( "hltESPDummyDetLayerGeometry" ),
  errorRescaling = cms.double( 100.0 ),
  minHits = cms.int32( 3 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPKFTrajectorySmoother = cms.ESProducer( "KFTrajectorySmootherESProducer",
  ComponentName = cms.string( "hltESPKFTrajectorySmoother" ),
  Propagator = cms.string( "PropagatorWithMaterialParabolicMf" ),
  Updator = cms.string( "hltESPKFUpdator" ),
  Estimator = cms.string( "hltESPChi2MeasurementEstimator30" ),
  RecoGeometry = cms.string( "hltESPDummyDetLayerGeometry" ),
  errorRescaling = cms.double( 100.0 ),
  minHits = cms.int32( 3 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPKFTrajectoryFitterForLoopers = cms.ESProducer( "KFTrajectoryFitterESProducer",
  ComponentName = cms.string( "hltESPKFTrajectoryFitterForLoopers" ),
  Propagator = cms.string( "PropagatorWithMaterialForLoopers" ),
  Updator = cms.string( "hltESPKFUpdator" ),
  Estimator = cms.string( "hltESPChi2MeasurementEstimator30" ),
  RecoGeometry = cms.string( "hltESPGlobalDetLayerGeometry" ),
  minHits = cms.int32( 3 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPKFTrajectoryFitterForL2Muon = cms.ESProducer( "KFTrajectoryFitterESProducer",
  ComponentName = cms.string( "hltESPKFTrajectoryFitterForL2Muon" ),
  Propagator = cms.string( "hltESPFastSteppingHelixPropagatorAny" ),
  Updator = cms.string( "hltESPKFUpdator" ),
  Estimator = cms.string( "hltESPChi2MeasurementEstimator30" ),
  RecoGeometry = cms.string( "hltESPDummyDetLayerGeometry" ),
  minHits = cms.int32( 3 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPKFTrajectoryFitter = cms.ESProducer( "KFTrajectoryFitterESProducer",
  ComponentName = cms.string( "hltESPKFTrajectoryFitter" ),
  Propagator = cms.string( "PropagatorWithMaterialParabolicMf" ),
  Updator = cms.string( "hltESPKFUpdator" ),
  Estimator = cms.string( "hltESPChi2MeasurementEstimator30" ),
  RecoGeometry = cms.string( "hltESPDummyDetLayerGeometry" ),
  minHits = cms.int32( 3 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPKFFittingSmootherWithOutliersRejectionAndRK = cms.ESProducer( "KFFittingSmootherESProducer",
  ComponentName = cms.string( "hltESPKFFittingSmootherWithOutliersRejectionAndRK" ),
  Fitter = cms.string( "hltESPRKTrajectoryFitter" ),
  Smoother = cms.string( "hltESPRKTrajectorySmoother" ),
  EstimateCut = cms.double( 20.0 ),
  MaxFractionOutliers = cms.double( 0.3 ),
  MaxNumberOfOutliers = cms.int32( 3 ),
  MinDof = cms.int32( 2 ),
  NoOutliersBeginEnd = cms.bool( False ),
  MinNumberOfHits = cms.int32( 3 ),
  MinNumberOfHitsHighEta = cms.int32( 5 ),
  HighEtaSwitch = cms.double( 5.0 ),
  RejectTracks = cms.bool( True ),
  BreakTrajWith2ConsecutiveMissing = cms.bool( True ),
  NoInvalidHitsBeginEnd = cms.bool( True ),
  LogPixelProbabilityCut = cms.double( -14.0 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPKFFittingSmootherForLoopers = cms.ESProducer( "KFFittingSmootherESProducer",
  ComponentName = cms.string( "hltESPKFFittingSmootherForLoopers" ),
  Fitter = cms.string( "hltESPKFTrajectoryFitterForLoopers" ),
  Smoother = cms.string( "hltESPKFTrajectorySmootherForLoopers" ),
  EstimateCut = cms.double( 20.0 ),
  MaxFractionOutliers = cms.double( 0.3 ),
  MaxNumberOfOutliers = cms.int32( 3 ),
  MinDof = cms.int32( 2 ),
  NoOutliersBeginEnd = cms.bool( False ),
  MinNumberOfHits = cms.int32( 3 ),
  MinNumberOfHitsHighEta = cms.int32( 5 ),
  HighEtaSwitch = cms.double( 5.0 ),
  RejectTracks = cms.bool( True ),
  BreakTrajWith2ConsecutiveMissing = cms.bool( True ),
  NoInvalidHitsBeginEnd = cms.bool( True ),
  LogPixelProbabilityCut = cms.double( -14.0 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPKFFittingSmootherForL2Muon = cms.ESProducer( "KFFittingSmootherESProducer",
  ComponentName = cms.string( "hltESPKFFittingSmootherForL2Muon" ),
  Fitter = cms.string( "hltESPKFTrajectoryFitterForL2Muon" ),
  Smoother = cms.string( "hltESPKFTrajectorySmootherForL2Muon" ),
  EstimateCut = cms.double( -1.0 ),
  MaxFractionOutliers = cms.double( 0.3 ),
  MaxNumberOfOutliers = cms.int32( 3 ),
  MinDof = cms.int32( 2 ),
  NoOutliersBeginEnd = cms.bool( False ),
  MinNumberOfHits = cms.int32( 5 ),
  MinNumberOfHitsHighEta = cms.int32( 5 ),
  HighEtaSwitch = cms.double( 5.0 ),
  RejectTracks = cms.bool( True ),
  BreakTrajWith2ConsecutiveMissing = cms.bool( False ),
  NoInvalidHitsBeginEnd = cms.bool( False ),
  LogPixelProbabilityCut = cms.double( -16.0 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPKFFittingSmoother = cms.ESProducer( "KFFittingSmootherESProducer",
  ComponentName = cms.string( "hltESPKFFittingSmoother" ),
  Fitter = cms.string( "hltESPKFTrajectoryFitter" ),
  Smoother = cms.string( "hltESPKFTrajectorySmoother" ),
  EstimateCut = cms.double( -1.0 ),
  MaxFractionOutliers = cms.double( 0.3 ),
  MaxNumberOfOutliers = cms.int32( 3 ),
  MinDof = cms.int32( 2 ),
  NoOutliersBeginEnd = cms.bool( False ),
  MinNumberOfHits = cms.int32( 5 ),
  MinNumberOfHitsHighEta = cms.int32( 5 ),
  HighEtaSwitch = cms.double( 5.0 ),
  RejectTracks = cms.bool( True ),
  BreakTrajWith2ConsecutiveMissing = cms.bool( False ),
  NoInvalidHitsBeginEnd = cms.bool( False ),
  LogPixelProbabilityCut = cms.double( -16.0 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPInitialStepChi2MeasurementEstimator36 = cms.ESProducer( "Chi2MeasurementEstimatorESProducer",
  MaxChi2 = cms.double( 36.0 ),
  nSigma = cms.double( 3.0 ),
  MaxDisplacement = cms.double( 100.0 ),
  MaxSagitta = cms.double( -1.0 ),
  MinimalTolerance = cms.double( 10.0 ),
  MinPtForHitRecoveryInGluedDet = cms.double( 1000000.0 ),
  ComponentName = cms.string( "hltESPInitialStepChi2MeasurementEstimator36" ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPInitialStepChi2ChargeMeasurementEstimator30 = cms.ESProducer( "Chi2ChargeMeasurementEstimatorESProducer",
  MaxChi2 = cms.double( 30.0 ),
  nSigma = cms.double( 3.0 ),
  MaxDisplacement = cms.double( 0.5 ),
  MaxSagitta = cms.double( 2.0 ),
  MinimalTolerance = cms.double( 0.5 ),
  MinPtForHitRecoveryInGluedDet = cms.double( 1000000.0 ),
  ComponentName = cms.string( "hltESPInitialStepChi2ChargeMeasurementEstimator30" ),
  pTChargeCutThreshold = cms.double( 15.0 ),
  clusterChargeCut = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutLoose" ) ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPHighPtTripletStepChi2ChargeMeasurementEstimator30 = cms.ESProducer( "Chi2ChargeMeasurementEstimatorESProducer",
  MaxChi2 = cms.double( 30.0 ),
  nSigma = cms.double( 3.0 ),
  MaxDisplacement = cms.double( 0.5 ),
  MaxSagitta = cms.double( 2.0 ),
  MinimalTolerance = cms.double( 0.5 ),
  MinPtForHitRecoveryInGluedDet = cms.double( 1000000.0 ),
  ComponentName = cms.string( "hltESPHighPtTripletStepChi2ChargeMeasurementEstimator30" ),
  pTChargeCutThreshold = cms.double( 15.0 ),
  clusterChargeCut = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutLoose" ) ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPGsfTrajectorySmoother = cms.ESProducer( "GsfTrajectorySmootherESProducer",
  Merger = cms.string( "hltESPCloseComponentsMerger5D" ),
  ComponentName = cms.string( "hltESPGsfTrajectorySmoother" ),
  MaterialEffectsUpdator = cms.string( "hltESPElectronMaterialEffects" ),
  ErrorRescaling = cms.double( 100.0 ),
  GeometricalPropagator = cms.string( "hltESPBwdAnalyticalPropagator" ),
  RecoGeometry = cms.string( "hltESPGlobalDetLayerGeometry" )
)
process.hltESPGsfTrajectoryFitter = cms.ESProducer( "GsfTrajectoryFitterESProducer",
  Merger = cms.string( "hltESPCloseComponentsMerger5D" ),
  ComponentName = cms.string( "hltESPGsfTrajectoryFitter" ),
  MaterialEffectsUpdator = cms.string( "hltESPElectronMaterialEffects" ),
  GeometricalPropagator = cms.string( "hltESPAnalyticalPropagator" ),
  RecoGeometry = cms.string( "hltESPGlobalDetLayerGeometry" )
)
process.hltESPGsfElectronFittingSmoother = cms.ESProducer( "KFFittingSmootherESProducer",
  ComponentName = cms.string( "hltESPGsfElectronFittingSmoother" ),
  Fitter = cms.string( "hltESPGsfTrajectoryFitter" ),
  Smoother = cms.string( "hltESPGsfTrajectorySmoother" ),
  EstimateCut = cms.double( -1.0 ),
  MaxFractionOutliers = cms.double( 0.3 ),
  MaxNumberOfOutliers = cms.int32( 3 ),
  MinDof = cms.int32( 2 ),
  NoOutliersBeginEnd = cms.bool( False ),
  MinNumberOfHits = cms.int32( 5 ),
  MinNumberOfHitsHighEta = cms.int32( 5 ),
  HighEtaSwitch = cms.double( 5.0 ),
  RejectTracks = cms.bool( True ),
  BreakTrajWith2ConsecutiveMissing = cms.bool( True ),
  NoInvalidHitsBeginEnd = cms.bool( True ),
  LogPixelProbabilityCut = cms.double( -16.0 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPGlobalTrackingGeometryESProducer = cms.ESProducer( "GlobalTrackingGeometryESProducer" )
process.hltESPGlobalDetLayerGeometry = cms.ESProducer( "GlobalDetLayerGeometryESProducer",
  ComponentName = cms.string( "hltESPGlobalDetLayerGeometry" )
)
process.hltESPFwdElectronPropagator = cms.ESProducer( "PropagatorWithMaterialESProducer",
  SimpleMagneticField = cms.string( "" ),
  MaxDPhi = cms.double( 1.6 ),
  ComponentName = cms.string( "hltESPFwdElectronPropagator" ),
  Mass = cms.double( 5.11E-4 ),
  PropagationDirection = cms.string( "alongMomentum" ),
  useRungeKutta = cms.bool( False ),
  ptMin = cms.double( -1.0 )
)
process.hltESPFlexibleKFFittingSmoother = cms.ESProducer( "FlexibleKFFittingSmootherESProducer",
  ComponentName = cms.string( "hltESPFlexibleKFFittingSmoother" ),
  standardFitter = cms.string( "hltESPKFFittingSmootherWithOutliersRejectionAndRK" ),
  looperFitter = cms.string( "hltESPKFFittingSmootherForLoopers" ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPFittingSmootherRK = cms.ESProducer( "KFFittingSmootherESProducer",
  ComponentName = cms.string( "hltESPFittingSmootherRK" ),
  Fitter = cms.string( "hltESPTrajectoryFitterRK" ),
  Smoother = cms.string( "hltESPTrajectorySmootherRK" ),
  EstimateCut = cms.double( -1.0 ),
  MaxFractionOutliers = cms.double( 0.3 ),
  MaxNumberOfOutliers = cms.int32( 3 ),
  MinDof = cms.int32( 2 ),
  NoOutliersBeginEnd = cms.bool( False ),
  MinNumberOfHits = cms.int32( 5 ),
  MinNumberOfHitsHighEta = cms.int32( 5 ),
  HighEtaSwitch = cms.double( 5.0 ),
  RejectTracks = cms.bool( True ),
  BreakTrajWith2ConsecutiveMissing = cms.bool( False ),
  NoInvalidHitsBeginEnd = cms.bool( False ),
  LogPixelProbabilityCut = cms.double( -16.0 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPFittingSmootherIT = cms.ESProducer( "KFFittingSmootherESProducer",
  ComponentName = cms.string( "hltESPFittingSmootherIT" ),
  Fitter = cms.string( "hltESPTrajectoryFitterRK" ),
  Smoother = cms.string( "hltESPTrajectorySmootherRK" ),
  EstimateCut = cms.double( -1.0 ),
  MaxFractionOutliers = cms.double( 0.3 ),
  MaxNumberOfOutliers = cms.int32( 3 ),
  MinDof = cms.int32( 2 ),
  NoOutliersBeginEnd = cms.bool( False ),
  MinNumberOfHits = cms.int32( 3 ),
  MinNumberOfHitsHighEta = cms.int32( 5 ),
  HighEtaSwitch = cms.double( 5.0 ),
  RejectTracks = cms.bool( True ),
  BreakTrajWith2ConsecutiveMissing = cms.bool( True ),
  NoInvalidHitsBeginEnd = cms.bool( True ),
  LogPixelProbabilityCut = cms.double( -16.0 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPFastSteppingHelixPropagatorOpposite = cms.ESProducer( "SteppingHelixPropagatorESProducer",
  ComponentName = cms.string( "hltESPFastSteppingHelixPropagatorOpposite" ),
  NoErrorPropagation = cms.bool( False ),
  PropagationDirection = cms.string( "oppositeToMomentum" ),
  useTuningForL2Speed = cms.bool( True ),
  useIsYokeFlag = cms.bool( True ),
  endcapShiftInZNeg = cms.double( 0.0 ),
  SetVBFPointer = cms.bool( False ),
  AssumeNoMaterial = cms.bool( False ),
  endcapShiftInZPos = cms.double( 0.0 ),
  useInTeslaFromMagField = cms.bool( False ),
  VBFName = cms.string( "VolumeBasedMagneticField" ),
  useEndcapShiftsInZ = cms.bool( False ),
  sendLogWarning = cms.bool( False ),
  useMatVolumes = cms.bool( True ),
  debug = cms.bool( False ),
  ApplyRadX0Correction = cms.bool( True ),
  useMagVolumes = cms.bool( True ),
  returnTangentPlane = cms.bool( True )
)
process.hltESPFastSteppingHelixPropagatorAny = cms.ESProducer( "SteppingHelixPropagatorESProducer",
  ComponentName = cms.string( "hltESPFastSteppingHelixPropagatorAny" ),
  NoErrorPropagation = cms.bool( False ),
  PropagationDirection = cms.string( "anyDirection" ),
  useTuningForL2Speed = cms.bool( True ),
  useIsYokeFlag = cms.bool( True ),
  endcapShiftInZNeg = cms.double( 0.0 ),
  SetVBFPointer = cms.bool( False ),
  AssumeNoMaterial = cms.bool( False ),
  endcapShiftInZPos = cms.double( 0.0 ),
  useInTeslaFromMagField = cms.bool( False ),
  VBFName = cms.string( "VolumeBasedMagneticField" ),
  useEndcapShiftsInZ = cms.bool( False ),
  sendLogWarning = cms.bool( False ),
  useMatVolumes = cms.bool( True ),
  debug = cms.bool( False ),
  ApplyRadX0Correction = cms.bool( True ),
  useMagVolumes = cms.bool( True ),
  returnTangentPlane = cms.bool( True )
)
process.hltESPElectronMaterialEffects = cms.ESProducer( "GsfMaterialEffectsESProducer",
  BetheHeitlerParametrization = cms.string( "BetheHeitler_cdfmom_nC6_O5.par" ),
  EnergyLossUpdator = cms.string( "GsfBetheHeitlerUpdator" ),
  ComponentName = cms.string( "hltESPElectronMaterialEffects" ),
  MultipleScatteringUpdator = cms.string( "MultipleScatteringUpdator" ),
  Mass = cms.double( 5.11E-4 ),
  BetheHeitlerCorrection = cms.int32( 2 )
)
process.hltESPEcalTrigTowerConstituentsMapBuilder = cms.ESProducer( "EcalTrigTowerConstituentsMapBuilder",
  MapFile = cms.untracked.string( "Geometry/EcalMapping/data/EndCap_TTMap.txt" )
)
process.hltESPDummyDetLayerGeometry = cms.ESProducer( "DetLayerGeometryESProducer",
  ComponentName = cms.string( "hltESPDummyDetLayerGeometry" )
)
process.hltESPDisplacedDijethltTrackCounting2D2ndLong = cms.ESProducer( "TrackCountingESProducer",
  a_dR = cms.double( -0.001053 ),
  b_dR = cms.double( 0.6263 ),
  a_pT = cms.double( 0.005263 ),
  b_pT = cms.double( 0.3684 ),
  min_pT = cms.double( 120.0 ),
  max_pT = cms.double( 500.0 ),
  min_pT_dRcut = cms.double( 0.5 ),
  max_pT_dRcut = cms.double( 0.1 ),
  max_pT_trackPTcut = cms.double( 3.0 ),
  minimumImpactParameter = cms.double( 0.2 ),
  useSignedImpactParameterSig = cms.bool( True ),
  impactParameterType = cms.int32( 1 ),
  maximumDistanceToJetAxis = cms.double( 9999999.0 ),
  deltaR = cms.double( -1.0 ),
  maximumDecayLength = cms.double( 999999.0 ),
  nthTrack = cms.int32( 2 ),
  trackQualityClass = cms.string( "any" ),
  useVariableJTA = cms.bool( False )
)
process.hltESPDisplacedDijethltTrackCounting2D1stLoose = cms.ESProducer( "TrackCountingESProducer",
  a_dR = cms.double( -0.001053 ),
  b_dR = cms.double( 0.6263 ),
  a_pT = cms.double( 0.005263 ),
  b_pT = cms.double( 0.3684 ),
  min_pT = cms.double( 120.0 ),
  max_pT = cms.double( 500.0 ),
  min_pT_dRcut = cms.double( 0.5 ),
  max_pT_dRcut = cms.double( 0.1 ),
  max_pT_trackPTcut = cms.double( 3.0 ),
  minimumImpactParameter = cms.double( 0.03 ),
  useSignedImpactParameterSig = cms.bool( False ),
  impactParameterType = cms.int32( 1 ),
  maximumDistanceToJetAxis = cms.double( 9999999.0 ),
  deltaR = cms.double( -1.0 ),
  maximumDecayLength = cms.double( 999999.0 ),
  nthTrack = cms.int32( 1 ),
  trackQualityClass = cms.string( "any" ),
  useVariableJTA = cms.bool( False )
)
process.hltESPDisplacedDijethltTrackCounting2D1st = cms.ESProducer( "TrackCountingESProducer",
  a_dR = cms.double( -0.001053 ),
  b_dR = cms.double( 0.6263 ),
  a_pT = cms.double( 0.005263 ),
  b_pT = cms.double( 0.3684 ),
  min_pT = cms.double( 120.0 ),
  max_pT = cms.double( 500.0 ),
  min_pT_dRcut = cms.double( 0.5 ),
  max_pT_dRcut = cms.double( 0.1 ),
  max_pT_trackPTcut = cms.double( 3.0 ),
  minimumImpactParameter = cms.double( 0.05 ),
  useSignedImpactParameterSig = cms.bool( False ),
  impactParameterType = cms.int32( 1 ),
  maximumDistanceToJetAxis = cms.double( 9999999.0 ),
  deltaR = cms.double( -1.0 ),
  maximumDecayLength = cms.double( 999999.0 ),
  nthTrack = cms.int32( 1 ),
  trackQualityClass = cms.string( "any" ),
  useVariableJTA = cms.bool( False )
)
process.hltESPDisplacedDijethltPromptTrackCountingESProducerShortSig5 = cms.ESProducer( "PromptTrackCountingESProducer",
  impactParameterType = cms.int32( 1 ),
  minimumImpactParameter = cms.double( -1.0 ),
  useSignedImpactParameterSig = cms.bool( False ),
  maximumDistanceToJetAxis = cms.double( 999999.0 ),
  deltaR = cms.double( -1.0 ),
  deltaRmin = cms.double( 0.0 ),
  maximumDecayLength = cms.double( 999999.0 ),
  maxImpactParameter = cms.double( 0.05 ),
  maxImpactParameterSig = cms.double( 5.0 ),
  trackQualityClass = cms.string( "any" ),
  nthTrack = cms.int32( -1 )
)
process.hltESPDisplacedDijethltPromptTrackCountingESProducerLong = cms.ESProducer( "PromptTrackCountingESProducer",
  impactParameterType = cms.int32( 1 ),
  minimumImpactParameter = cms.double( -1.0 ),
  useSignedImpactParameterSig = cms.bool( True ),
  maximumDistanceToJetAxis = cms.double( 999999.0 ),
  deltaR = cms.double( -1.0 ),
  deltaRmin = cms.double( 0.0 ),
  maximumDecayLength = cms.double( 999999.0 ),
  maxImpactParameter = cms.double( 0.2 ),
  maxImpactParameterSig = cms.double( 999999.0 ),
  trackQualityClass = cms.string( "any" ),
  nthTrack = cms.int32( -1 )
)
process.hltESPDisplacedDijethltPromptTrackCountingESProducer = cms.ESProducer( "PromptTrackCountingESProducer",
  impactParameterType = cms.int32( 1 ),
  minimumImpactParameter = cms.double( -1.0 ),
  useSignedImpactParameterSig = cms.bool( True ),
  maximumDistanceToJetAxis = cms.double( 999999.0 ),
  deltaR = cms.double( -1.0 ),
  deltaRmin = cms.double( 0.0 ),
  maximumDecayLength = cms.double( 999999.0 ),
  maxImpactParameter = cms.double( 0.1 ),
  maxImpactParameterSig = cms.double( 999999.0 ),
  trackQualityClass = cms.string( "any" ),
  nthTrack = cms.int32( -1 )
)
process.hltESPDetachedTripletStepTrajectoryCleanerBySharedHits = cms.ESProducer( "TrajectoryCleanerESProducer",
  ComponentName = cms.string( "hltESPDetachedTripletStepTrajectoryCleanerBySharedHits" ),
  ComponentType = cms.string( "TrajectoryCleanerBySharedHits" ),
  fractionShared = cms.double( 0.13 ),
  ValidHitBonus = cms.double( 5.0 ),
  MissingHitPenalty = cms.double( 20.0 ),
  allowSharedFirstHit = cms.bool( True )
)
process.hltESPDetachedTripletStepChi2ChargeMeasurementEstimator9 = cms.ESProducer( "Chi2ChargeMeasurementEstimatorESProducer",
  MaxChi2 = cms.double( 9.0 ),
  nSigma = cms.double( 3.0 ),
  MaxDisplacement = cms.double( 0.5 ),
  MaxSagitta = cms.double( 2.0 ),
  MinimalTolerance = cms.double( 0.5 ),
  MinPtForHitRecoveryInGluedDet = cms.double( 1000000.0 ),
  ComponentName = cms.string( "hltESPDetachedTripletStepChi2ChargeMeasurementEstimator9" ),
  pTChargeCutThreshold = cms.double( -1.0 ),
  clusterChargeCut = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutTight" ) ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPDetachedStepTrajectoryCleanerBySharedHits = cms.ESProducer( "TrajectoryCleanerESProducer",
  ComponentName = cms.string( "hltESPDetachedStepTrajectoryCleanerBySharedHits" ),
  ComponentType = cms.string( "TrajectoryCleanerBySharedHits" ),
  fractionShared = cms.double( 0.13 ),
  ValidHitBonus = cms.double( 5.0 ),
  MissingHitPenalty = cms.double( 20.0 ),
  allowSharedFirstHit = cms.bool( True )
)
process.hltESPDetachedQuadStepTrajectoryCleanerBySharedHits = cms.ESProducer( "TrajectoryCleanerESProducer",
  ComponentName = cms.string( "hltESPDetachedQuadStepTrajectoryCleanerBySharedHits" ),
  ComponentType = cms.string( "TrajectoryCleanerBySharedHits" ),
  fractionShared = cms.double( 0.13 ),
  ValidHitBonus = cms.double( 5.0 ),
  MissingHitPenalty = cms.double( 20.0 ),
  allowSharedFirstHit = cms.bool( True )
)
process.hltESPDetachedQuadStepChi2ChargeMeasurementEstimator9 = cms.ESProducer( "Chi2ChargeMeasurementEstimatorESProducer",
  MaxChi2 = cms.double( 9.0 ),
  nSigma = cms.double( 3.0 ),
  MaxDisplacement = cms.double( 0.5 ),
  MaxSagitta = cms.double( 2.0 ),
  MinimalTolerance = cms.double( 0.5 ),
  MinPtForHitRecoveryInGluedDet = cms.double( 1000000.0 ),
  ComponentName = cms.string( "hltESPDetachedQuadStepChi2ChargeMeasurementEstimator9" ),
  pTChargeCutThreshold = cms.double( -1.0 ),
  clusterChargeCut = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutTight" ) ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPCloseComponentsMerger5D = cms.ESProducer( "CloseComponentsMergerESProducer5D",
  ComponentName = cms.string( "hltESPCloseComponentsMerger5D" ),
  MaxComponents = cms.int32( 12 ),
  DistanceMeasure = cms.string( "hltESPKullbackLeiblerDistance5D" )
)
process.hltESPChi2MeasurementEstimator9 = cms.ESProducer( "Chi2MeasurementEstimatorESProducer",
  MaxChi2 = cms.double( 9.0 ),
  nSigma = cms.double( 3.0 ),
  MaxDisplacement = cms.double( 100.0 ),
  MaxSagitta = cms.double( -1.0 ),
  MinimalTolerance = cms.double( 10.0 ),
  MinPtForHitRecoveryInGluedDet = cms.double( 1000000.0 ),
  ComponentName = cms.string( "hltESPChi2MeasurementEstimator9" ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPChi2MeasurementEstimator30 = cms.ESProducer( "Chi2MeasurementEstimatorESProducer",
  MaxChi2 = cms.double( 30.0 ),
  nSigma = cms.double( 3.0 ),
  MaxDisplacement = cms.double( 100.0 ),
  MaxSagitta = cms.double( -1.0 ),
  MinimalTolerance = cms.double( 10.0 ),
  MinPtForHitRecoveryInGluedDet = cms.double( 1000000.0 ),
  ComponentName = cms.string( "hltESPChi2MeasurementEstimator30" ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPChi2MeasurementEstimator16 = cms.ESProducer( "Chi2MeasurementEstimatorESProducer",
  MaxChi2 = cms.double( 16.0 ),
  nSigma = cms.double( 3.0 ),
  MaxDisplacement = cms.double( 100.0 ),
  MaxSagitta = cms.double( -1.0 ),
  MinimalTolerance = cms.double( 10.0 ),
  MinPtForHitRecoveryInGluedDet = cms.double( 1000000.0 ),
  ComponentName = cms.string( "hltESPChi2MeasurementEstimator16" ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPChi2MeasurementEstimator100 = cms.ESProducer( "Chi2MeasurementEstimatorESProducer",
  MaxChi2 = cms.double( 40.0 ),
  nSigma = cms.double( 4.0 ),
  MaxDisplacement = cms.double( 0.5 ),
  MaxSagitta = cms.double( 2.0 ),
  MinimalTolerance = cms.double( 0.5 ),
  MinPtForHitRecoveryInGluedDet = cms.double( 1.0E12 ),
  ComponentName = cms.string( "hltESPChi2MeasurementEstimator100" ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPChi2ChargeTightMeasurementEstimator16 = cms.ESProducer( "Chi2ChargeMeasurementEstimatorESProducer",
  MaxChi2 = cms.double( 16.0 ),
  nSigma = cms.double( 3.0 ),
  MaxDisplacement = cms.double( 0.5 ),
  MaxSagitta = cms.double( 2.0 ),
  MinimalTolerance = cms.double( 0.5 ),
  MinPtForHitRecoveryInGluedDet = cms.double( 1000000.0 ),
  ComponentName = cms.string( "hltESPChi2ChargeTightMeasurementEstimator16" ),
  pTChargeCutThreshold = cms.double( -1.0 ),
  clusterChargeCut = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutTight" ) ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPChi2ChargeMeasurementEstimator9ForHI = cms.ESProducer( "Chi2ChargeMeasurementEstimatorESProducer",
  MaxChi2 = cms.double( 9.0 ),
  nSigma = cms.double( 3.0 ),
  MaxDisplacement = cms.double( 100.0 ),
  MaxSagitta = cms.double( -1.0 ),
  MinimalTolerance = cms.double( 10.0 ),
  MinPtForHitRecoveryInGluedDet = cms.double( 1000000.0 ),
  ComponentName = cms.string( "hltESPChi2ChargeMeasurementEstimator9ForHI" ),
  pTChargeCutThreshold = cms.double( 15.0 ),
  clusterChargeCut = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutForHI" ) ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPChi2ChargeMeasurementEstimator9 = cms.ESProducer( "Chi2ChargeMeasurementEstimatorESProducer",
  MaxChi2 = cms.double( 9.0 ),
  nSigma = cms.double( 3.0 ),
  MaxDisplacement = cms.double( 0.5 ),
  MaxSagitta = cms.double( 2.0 ),
  MinimalTolerance = cms.double( 0.5 ),
  MinPtForHitRecoveryInGluedDet = cms.double( 1000000.0 ),
  ComponentName = cms.string( "hltESPChi2ChargeMeasurementEstimator9" ),
  pTChargeCutThreshold = cms.double( 15.0 ),
  clusterChargeCut = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutLoose" ) ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPChi2ChargeMeasurementEstimator30 = cms.ESProducer( "Chi2ChargeMeasurementEstimatorESProducer",
  MaxChi2 = cms.double( 30.0 ),
  nSigma = cms.double( 3.0 ),
  MaxDisplacement = cms.double( 100.0 ),
  MaxSagitta = cms.double( -1.0 ),
  MinimalTolerance = cms.double( 10.0 ),
  MinPtForHitRecoveryInGluedDet = cms.double( 1000000.0 ),
  ComponentName = cms.string( "hltESPChi2ChargeMeasurementEstimator30" ),
  pTChargeCutThreshold = cms.double( -1.0 ),
  clusterChargeCut = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPChi2ChargeMeasurementEstimator2000 = cms.ESProducer( "Chi2ChargeMeasurementEstimatorESProducer",
  MaxChi2 = cms.double( 2000.0 ),
  nSigma = cms.double( 3.0 ),
  MaxDisplacement = cms.double( 100.0 ),
  MaxSagitta = cms.double( -1.0 ),
  MinimalTolerance = cms.double( 10.0 ),
  MinPtForHitRecoveryInGluedDet = cms.double( 1000000.0 ),
  ComponentName = cms.string( "hltESPChi2ChargeMeasurementEstimator2000" ),
  pTChargeCutThreshold = cms.double( -1.0 ),
  clusterChargeCut = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPChi2ChargeMeasurementEstimator16 = cms.ESProducer( "Chi2ChargeMeasurementEstimatorESProducer",
  MaxChi2 = cms.double( 16.0 ),
  nSigma = cms.double( 3.0 ),
  MaxDisplacement = cms.double( 0.5 ),
  MaxSagitta = cms.double( 2.0 ),
  MinimalTolerance = cms.double( 0.5 ),
  MinPtForHitRecoveryInGluedDet = cms.double( 1000000.0 ),
  ComponentName = cms.string( "hltESPChi2ChargeMeasurementEstimator16" ),
  pTChargeCutThreshold = cms.double( -1.0 ),
  clusterChargeCut = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutLoose" ) ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPChi2ChargeLooseMeasurementEstimator16 = cms.ESProducer( "Chi2ChargeMeasurementEstimatorESProducer",
  MaxChi2 = cms.double( 16.0 ),
  nSigma = cms.double( 3.0 ),
  MaxDisplacement = cms.double( 0.5 ),
  MaxSagitta = cms.double( 2.0 ),
  MinimalTolerance = cms.double( 0.5 ),
  MinPtForHitRecoveryInGluedDet = cms.double( 1000000.0 ),
  ComponentName = cms.string( "hltESPChi2ChargeLooseMeasurementEstimator16" ),
  pTChargeCutThreshold = cms.double( -1.0 ),
  clusterChargeCut = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutLoose" ) ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPBwdElectronPropagator = cms.ESProducer( "PropagatorWithMaterialESProducer",
  SimpleMagneticField = cms.string( "" ),
  MaxDPhi = cms.double( 1.6 ),
  ComponentName = cms.string( "hltESPBwdElectronPropagator" ),
  Mass = cms.double( 5.11E-4 ),
  PropagationDirection = cms.string( "oppositeToMomentum" ),
  useRungeKutta = cms.bool( False ),
  ptMin = cms.double( -1.0 )
)
process.hltESPBwdAnalyticalPropagator = cms.ESProducer( "AnalyticalPropagatorESProducer",
  MaxDPhi = cms.double( 1.6 ),
  ComponentName = cms.string( "hltESPBwdAnalyticalPropagator" ),
  PropagationDirection = cms.string( "oppositeToMomentum" )
)
process.hltESPAnalyticalPropagator = cms.ESProducer( "AnalyticalPropagatorESProducer",
  MaxDPhi = cms.double( 1.6 ),
  ComponentName = cms.string( "hltESPAnalyticalPropagator" ),
  PropagationDirection = cms.string( "alongMomentum" )
)
process.hltDisplacedDijethltESPTrackCounting2D1st = cms.ESProducer( "TrackCountingESProducer",
  a_dR = cms.double( -0.001053 ),
  b_dR = cms.double( 0.6263 ),
  a_pT = cms.double( 0.005263 ),
  b_pT = cms.double( 0.3684 ),
  min_pT = cms.double( 120.0 ),
  max_pT = cms.double( 500.0 ),
  min_pT_dRcut = cms.double( 0.5 ),
  max_pT_dRcut = cms.double( 0.1 ),
  max_pT_trackPTcut = cms.double( 3.0 ),
  minimumImpactParameter = cms.double( 0.05 ),
  useSignedImpactParameterSig = cms.bool( False ),
  impactParameterType = cms.int32( 1 ),
  maximumDistanceToJetAxis = cms.double( 9999999.0 ),
  deltaR = cms.double( -1.0 ),
  maximumDecayLength = cms.double( 999999.0 ),
  nthTrack = cms.int32( 1 ),
  trackQualityClass = cms.string( "any" ),
  useVariableJTA = cms.bool( False )
)
process.hltDisplacedDijethltESPPromptTrackCountingESProducer = cms.ESProducer( "PromptTrackCountingESProducer",
  impactParameterType = cms.int32( 1 ),
  minimumImpactParameter = cms.double( -1.0 ),
  useSignedImpactParameterSig = cms.bool( True ),
  maximumDistanceToJetAxis = cms.double( 999999.0 ),
  deltaR = cms.double( -1.0 ),
  deltaRmin = cms.double( 0.0 ),
  maximumDecayLength = cms.double( 999999.0 ),
  maxImpactParameter = cms.double( 0.1 ),
  maxImpactParameterSig = cms.double( 999999.0 ),
  trackQualityClass = cms.string( "any" ),
  nthTrack = cms.int32( -1 )
)
process.hltCombinedSecondaryVertexV2 = cms.ESProducer( "CombinedSecondaryVertexESProducer",
  trackPseudoSelection = cms.PSet( 
    max_pT_dRcut = cms.double( 0.1 ),
    b_dR = cms.double( 0.6263 ),
    min_pT = cms.double( 120.0 ),
    b_pT = cms.double( 0.3684 ),
    ptMin = cms.double( 0.0 ),
    max_pT_trackPTcut = cms.double( 3.0 ),
    max_pT = cms.double( 500.0 ),
    useVariableJTA = cms.bool( False ),
    maxDecayLen = cms.double( 5.0 ),
    qualityClass = cms.string( "any" ),
    normChi2Max = cms.double( 99999.9 ),
    sip2dValMin = cms.double( -99999.9 ),
    sip3dValMin = cms.double( -99999.9 ),
    a_dR = cms.double( -0.001053 ),
    maxDistToAxis = cms.double( 0.07 ),
    totalHitsMin = cms.uint32( 0 ),
    a_pT = cms.double( 0.005263 ),
    sip2dSigMax = cms.double( 99999.9 ),
    sip2dValMax = cms.double( 99999.9 ),
    sip3dSigMax = cms.double( 99999.9 ),
    sip3dValMax = cms.double( 99999.9 ),
    min_pT_dRcut = cms.double( 0.5 ),
    jetDeltaRMax = cms.double( 0.3 ),
    pixelHitsMin = cms.uint32( 0 ),
    sip3dSigMin = cms.double( -99999.9 ),
    sip2dSigMin = cms.double( 2.0 )
  ),
  trackSelection = cms.PSet( 
    max_pT_dRcut = cms.double( 0.1 ),
    b_dR = cms.double( 0.6263 ),
    min_pT = cms.double( 120.0 ),
    b_pT = cms.double( 0.3684 ),
    ptMin = cms.double( 0.0 ),
    max_pT_trackPTcut = cms.double( 3.0 ),
    max_pT = cms.double( 500.0 ),
    useVariableJTA = cms.bool( False ),
    maxDecayLen = cms.double( 5.0 ),
    qualityClass = cms.string( "any" ),
    normChi2Max = cms.double( 99999.9 ),
    sip2dValMin = cms.double( -99999.9 ),
    sip3dValMin = cms.double( -99999.9 ),
    a_dR = cms.double( -0.001053 ),
    maxDistToAxis = cms.double( 0.07 ),
    totalHitsMin = cms.uint32( 0 ),
    a_pT = cms.double( 0.005263 ),
    sip2dSigMax = cms.double( 99999.9 ),
    sip2dValMax = cms.double( 99999.9 ),
    sip3dSigMax = cms.double( 99999.9 ),
    sip3dValMax = cms.double( 99999.9 ),
    min_pT_dRcut = cms.double( 0.5 ),
    jetDeltaRMax = cms.double( 0.3 ),
    pixelHitsMin = cms.uint32( 0 ),
    sip3dSigMin = cms.double( -99999.9 ),
    sip2dSigMin = cms.double( -99999.9 )
  ),
  trackFlip = cms.bool( False ),
  vertexFlip = cms.bool( False ),
  SoftLeptonFlip = cms.bool( False ),
  useTrackWeights = cms.bool( True ),
  pseudoMultiplicityMin = cms.uint32( 2 ),
  correctVertexMass = cms.bool( True ),
  trackPairV0Filter = cms.PSet(  k0sMassWindow = cms.double( 0.03 ) ),
  charmCut = cms.double( 1.5 ),
  minimumTrackWeight = cms.double( 0.5 ),
  pseudoVertexV0Filter = cms.PSet(  k0sMassWindow = cms.double( 0.05 ) ),
  trackMultiplicityMin = cms.uint32( 3 ),
  trackSort = cms.string( "sip2dSig" ),
  useCategories = cms.bool( True ),
  calibrationRecords = cms.vstring( 'CombinedSVIVFV2RecoVertex',
    'CombinedSVIVFV2PseudoVertex',
    'CombinedSVIVFV2NoVertex' ),
  recordLabel = cms.string( "HLT" ),
  categoryVariableName = cms.string( "vertexCategory" )
)
process.hltCombinedSecondaryVertex = cms.ESProducer( "CombinedSecondaryVertexESProducer",
  trackPseudoSelection = cms.PSet( 
    maxDistToAxis = cms.double( 0.07 ),
    totalHitsMin = cms.uint32( 0 ),
    ptMin = cms.double( 0.0 ),
    sip2dSigMax = cms.double( 99999.9 ),
    sip2dValMax = cms.double( 99999.9 ),
    sip3dSigMax = cms.double( 99999.9 ),
    sip3dValMax = cms.double( 99999.9 ),
    maxDecayLen = cms.double( 5.0 ),
    qualityClass = cms.string( "any" ),
    jetDeltaRMax = cms.double( 0.3 ),
    normChi2Max = cms.double( 99999.9 ),
    pixelHitsMin = cms.uint32( 0 ),
    sip2dSigMin = cms.double( 2.0 ),
    sip2dValMin = cms.double( -99999.9 ),
    sip3dSigMin = cms.double( -99999.9 ),
    sip3dValMin = cms.double( -99999.9 )
  ),
  trackSelection = cms.PSet( 
    maxDistToAxis = cms.double( 0.07 ),
    totalHitsMin = cms.uint32( 0 ),
    ptMin = cms.double( 0.0 ),
    sip2dSigMax = cms.double( 99999.9 ),
    sip2dValMax = cms.double( 99999.9 ),
    sip3dSigMax = cms.double( 99999.9 ),
    sip3dValMax = cms.double( 99999.9 ),
    maxDecayLen = cms.double( 5.0 ),
    qualityClass = cms.string( "any" ),
    jetDeltaRMax = cms.double( 0.3 ),
    normChi2Max = cms.double( 99999.9 ),
    pixelHitsMin = cms.uint32( 0 ),
    sip2dSigMin = cms.double( -99999.9 ),
    sip2dValMin = cms.double( -99999.9 ),
    sip3dSigMin = cms.double( -99999.9 ),
    sip3dValMin = cms.double( -99999.9 )
  ),
  trackFlip = cms.bool( False ),
  vertexFlip = cms.bool( False ),
  SoftLeptonFlip = cms.bool( False ),
  useTrackWeights = cms.bool( True ),
  pseudoMultiplicityMin = cms.uint32( 2 ),
  correctVertexMass = cms.bool( True ),
  trackPairV0Filter = cms.PSet(  k0sMassWindow = cms.double( 0.03 ) ),
  charmCut = cms.double( 1.5 ),
  minimumTrackWeight = cms.double( 0.5 ),
  pseudoVertexV0Filter = cms.PSet(  k0sMassWindow = cms.double( 0.05 ) ),
  trackMultiplicityMin = cms.uint32( 3 ),
  trackSort = cms.string( "sip2dSig" ),
  useCategories = cms.bool( True ),
  calibrationRecords = cms.vstring( 'CombinedSVRecoVertex',
    'CombinedSVPseudoVertex',
    'CombinedSVNoVertex' ),
  recordLabel = cms.string( "HLT" ),
  categoryVariableName = cms.string( "vertexCategory" )
)
process.hltBoostedDoubleSecondaryVertexAK8Computer = cms.ESProducer( "CandidateBoostedDoubleSecondaryVertexESProducer",
  useCondDB = cms.bool( False ),
  weightFile = cms.FileInPath( "RecoBTag/SecondaryVertex/data/BoostedDoubleSV_AK8_BDT_v4.weights.xml.gz" ),
  useGBRForest = cms.bool( True ),
  useAdaBoost = cms.bool( False )
)
process.hcal_db_producer = cms.ESProducer( "HcalDbProducer" )
process.hcalRecAlgos = cms.ESProducer( "HcalRecAlgoESProducer",
  phase = cms.uint32( 1 ),
  RecoveredRecHitBits = cms.vstring(  ),
  SeverityLevels = cms.VPSet( 
    cms.PSet(  ChannelStatus = cms.vstring(  ),
      RecHitFlags = cms.vstring( 'TimingFromTDC' ),
      Level = cms.int32( 0 )
    ),
    cms.PSet(  ChannelStatus = cms.vstring( 'HcalCellCaloTowerProb' ),
      RecHitFlags = cms.vstring(  ),
      Level = cms.int32( 1 )
    ),
    cms.PSet(  ChannelStatus = cms.vstring( 'HcalCellExcludeFromHBHENoiseSummary' ),
      RecHitFlags = cms.vstring(  ),
      Level = cms.int32( 5 )
    ),
    cms.PSet(  ChannelStatus = cms.vstring(  ),
      RecHitFlags = cms.vstring( 'HBHEHpdHitMultiplicity',
        'HBHEIsolatedNoise',
        'HBHEFlatNoise',
        'HBHESpikeNoise',
        'HBHETS4TS5Noise',
        'HBHENegativeNoise',
        'HBHEPulseFitBit',
        'HBHEOOTPU' ),
      Level = cms.int32( 8 )
    ),
    cms.PSet(  ChannelStatus = cms.vstring(  ),
      RecHitFlags = cms.vstring( 'HFLongShort',
        'HFS8S1Ratio',
        'HFPET',
        'HFSignalAsymmetry' ),
      Level = cms.int32( 11 )
    ),
    cms.PSet(  ChannelStatus = cms.vstring( 'HcalCellHot' ),
      RecHitFlags = cms.vstring(  ),
      Level = cms.int32( 15 )
    ),
    cms.PSet(  ChannelStatus = cms.vstring( 'HcalCellOff',
  'HcalCellDead' ),
      RecHitFlags = cms.vstring(  ),
      Level = cms.int32( 20 )
    )
  ),
  DropChannelStatusBits = cms.vstring( 'HcalCellMask',
    'HcalCellOff',
    'HcalCellDead' ),
  appendToDataLabel = cms.string( "" )
)
process.hcalDetIdAssociator = cms.ESProducer( "DetIdAssociatorESProducer",
  ComponentName = cms.string( "HcalDetIdAssociator" ),
  etaBinSize = cms.double( 0.087 ),
  nEta = cms.int32( 70 ),
  nPhi = cms.int32( 72 ),
  hcalRegion = cms.int32( 2 ),
  includeBadChambers = cms.bool( False ),
  includeGEM = cms.bool( False ),
  includeME0 = cms.bool( False )
)
process.hcalDDDSimConstants = cms.ESProducer( "HcalDDDSimConstantsESModule",
  appendToDataLabel = cms.string( "" )
)
process.hcalDDDRecConstants = cms.ESProducer( "HcalDDDRecConstantsESModule",
  appendToDataLabel = cms.string( "" )
)
process.hcalChannelPropertiesESProd = cms.ESProducer( "HcalChannelPropertiesEP" )
process.ecalSeverityLevel = cms.ESProducer( "EcalSeverityLevelESProducer",
  flagMask = cms.PSet( 
    kBad = cms.vstring( 'kFaultyHardware',
      'kDead',
      'kKilled' ),
    kGood = cms.vstring( 'kGood' ),
    kRecovered = cms.vstring( 'kLeadingEdgeRecovered',
      'kTowerRecovered' ),
    kProblematic = cms.vstring( 'kPoorReco',
      'kPoorCalib',
      'kNoisy',
      'kSaturated' ),
    kWeird = cms.vstring( 'kWeird',
      'kDiWeird' ),
    kTime = cms.vstring( 'kOutOfTime' )
  ),
  dbstatusMask = cms.PSet( 
    kBad = cms.vstring( 'kNonRespondingIsolated',
      'kDeadVFE',
      'kDeadFE',
      'kNoDataNoTP' ),
    kGood = cms.vstring( 'kOk' ),
    kRecovered = cms.vstring(  ),
    kProblematic = cms.vstring( 'kDAC',
      'kNoLaser',
      'kNoisy',
      'kNNoisy',
      'kNNNoisy',
      'kNNNNoisy',
      'kNNNNNoisy',
      'kFixedG6',
      'kFixedG1',
      'kFixedG0' ),
    kWeird = cms.vstring(  ),
    kTime = cms.vstring(  )
  ),
  timeThresh = cms.double( 2.0 )
)
process.ecalDetIdAssociator = cms.ESProducer( "DetIdAssociatorESProducer",
  ComponentName = cms.string( "EcalDetIdAssociator" ),
  etaBinSize = cms.double( 0.02 ),
  nEta = cms.int32( 300 ),
  nPhi = cms.int32( 360 ),
  hcalRegion = cms.int32( 2 ),
  includeBadChambers = cms.bool( False ),
  includeGEM = cms.bool( False ),
  includeME0 = cms.bool( False )
)
process.ctppsInterpolatedOpticalFunctionsESSource = cms.ESProducer( "CTPPSInterpolatedOpticalFunctionsESSource",
  lhcInfoLabel = cms.string( "" ),
  opticsLabel = cms.string( "" ),
  appendToDataLabel = cms.string( "" )
)
process.ctppsGeometryESModule = cms.ESProducer( "CTPPSGeometryESModule",
  verbosity = cms.untracked.uint32( 1 ),
  buildMisalignedGeometry = cms.bool( False ),
  isRun2 = cms.bool( False ),
  dbTag = cms.string( "" ),
  compactViewTag = cms.string( "" ),
  fromPreprocessedDB = cms.untracked.bool( True ),
  fromDD4hep = cms.untracked.bool( False ),
  appendToDataLabel = cms.string( "" )
)
process.cosmicsNavigationSchoolESProducer = cms.ESProducer( "NavigationSchoolESProducer",
  ComponentName = cms.string( "CosmicNavigationSchool" ),
  SimpleMagneticField = cms.string( "" )
)
process.caloDetIdAssociator = cms.ESProducer( "DetIdAssociatorESProducer",
  ComponentName = cms.string( "CaloDetIdAssociator" ),
  etaBinSize = cms.double( 0.087 ),
  nEta = cms.int32( 70 ),
  nPhi = cms.int32( 72 ),
  hcalRegion = cms.int32( 2 ),
  includeBadChambers = cms.bool( False ),
  includeGEM = cms.bool( False ),
  includeME0 = cms.bool( False )
)
process.ZdcGeometryFromDBEP = cms.ESProducer( "ZdcGeometryFromDBEP",
  applyAlignment = cms.bool( False )
)
process.VolumeBasedMagneticFieldESProducer = cms.ESProducer( "VolumeBasedMagneticFieldESProducerFromDB",
  label = cms.untracked.string( "" ),
  debugBuilder = cms.untracked.bool( False ),
  valueOverride = cms.int32( -1 )
)
process.TransientTrackBuilderESProducer = cms.ESProducer( "TransientTrackBuilderESProducer",
  ComponentName = cms.string( "TransientTrackBuilder" )
)
process.TrackerGeometricDetESModule = cms.ESProducer( "TrackerGeometricDetESModule",
  fromDDD = cms.bool( False ),
  fromDD4hep = cms.bool( False ),
  appendToDataLabel = cms.string( "" )
)
process.TrackerDigiGeometryESModule = cms.ESProducer( "TrackerDigiGeometryESModule",
  appendToDataLabel = cms.string( "" ),
  fromDDD = cms.bool( False ),
  applyAlignment = cms.bool( True ),
  alignmentsLabel = cms.string( "" )
)
process.TrackerAdditionalParametersPerDetESModule = cms.ESProducer( "TrackerAdditionalParametersPerDetESModule",
  appendToDataLabel = cms.string( "" )
)
process.SteppingHelixPropagatorAny = cms.ESProducer( "SteppingHelixPropagatorESProducer",
  ComponentName = cms.string( "SteppingHelixPropagatorAny" ),
  NoErrorPropagation = cms.bool( False ),
  PropagationDirection = cms.string( "anyDirection" ),
  useTuningForL2Speed = cms.bool( False ),
  useIsYokeFlag = cms.bool( True ),
  endcapShiftInZNeg = cms.double( 0.0 ),
  SetVBFPointer = cms.bool( False ),
  AssumeNoMaterial = cms.bool( False ),
  endcapShiftInZPos = cms.double( 0.0 ),
  useInTeslaFromMagField = cms.bool( False ),
  VBFName = cms.string( "VolumeBasedMagneticField" ),
  useEndcapShiftsInZ = cms.bool( False ),
  sendLogWarning = cms.bool( False ),
  useMatVolumes = cms.bool( True ),
  debug = cms.bool( False ),
  ApplyRadX0Correction = cms.bool( True ),
  useMagVolumes = cms.bool( True ),
  returnTangentPlane = cms.bool( True )
)
process.SimpleSecondaryVertex3TrkComputer = cms.ESProducer( "SimpleSecondaryVertexESProducer",
  use3d = cms.bool( True ),
  unBoost = cms.bool( False ),
  useSignificance = cms.bool( True ),
  minTracks = cms.uint32( 3 ),
  minVertices = cms.uint32( 1 )
)
process.SiStripRegionConnectivity = cms.ESProducer( "SiStripRegionConnectivity",
  EtaDivisions = cms.untracked.uint32( 20 ),
  PhiDivisions = cms.untracked.uint32( 20 ),
  EtaMax = cms.untracked.double( 2.5 )
)
process.SiStripRecHitMatcherESProducer = cms.ESProducer( "SiStripRecHitMatcherESProducer",
  ComponentName = cms.string( "StandardMatcher" ),
  NSigmaInside = cms.double( 3.0 ),
  PreFilter = cms.bool( False )
)
process.SiStripQualityESProducer = cms.ESProducer( "SiStripQualityESProducer",
  appendToDataLabel = cms.string( "" ),
  ListOfRecordToMerge = cms.VPSet( 
    cms.PSet(  record = cms.string( "SiStripDetVOffRcd" ),
      tag = cms.string( "" )
    ),
    cms.PSet(  record = cms.string( "SiStripDetCablingRcd" ),
      tag = cms.string( "" )
    ),
    cms.PSet(  record = cms.string( "SiStripBadChannelRcd" ),
      tag = cms.string( "" )
    ),
    cms.PSet(  record = cms.string( "SiStripBadFiberRcd" ),
      tag = cms.string( "" )
    ),
    cms.PSet(  record = cms.string( "SiStripBadModuleRcd" ),
      tag = cms.string( "" )
    )
  ),
  ReduceGranularity = cms.bool( False ),
  ThresholdForReducedGranularity = cms.double( 0.3 ),
  PrintDebugOutput = cms.bool( False ),
  UseEmptyRunInfo = cms.bool( False )
)
process.SiStripGainESProducer = cms.ESProducer( "SiStripGainESProducer",
  appendToDataLabel = cms.string( "" ),
  printDebug = cms.untracked.bool( False ),
  AutomaticNormalization = cms.bool( False ),
  APVGain = cms.VPSet( 
    cms.PSet(  NormalizationFactor = cms.untracked.double( 1.0 ),
      Label = cms.untracked.string( "" ),
      Record = cms.string( "SiStripApvGainRcd" )
    ),
    cms.PSet(  NormalizationFactor = cms.untracked.double( 1.0 ),
      Label = cms.untracked.string( "" ),
      Record = cms.string( "SiStripApvGain2Rcd" )
    )
  )
)
process.SiStripClusterizerConditionsESProducer = cms.ESProducer( "SiStripClusterizerConditionsESProducer",
  QualityLabel = cms.string( "" ),
  Label = cms.string( "" ),
  appendToDataLabel = cms.string( "" )
)
process.RPCGeometryESModule = cms.ESProducer( "RPCGeometryESModule",
  fromDDD = cms.untracked.bool( False ),
  fromDD4hep = cms.untracked.bool( False ),
  appendToDataLabel = cms.string( "" )
)
process.PropagatorWithMaterialForMixedStep = cms.ESProducer( "PropagatorWithMaterialESProducer",
  SimpleMagneticField = cms.string( "ParabolicMf" ),
  MaxDPhi = cms.double( 1.6 ),
  ComponentName = cms.string( "PropagatorWithMaterialForMixedStep" ),
  Mass = cms.double( 0.105 ),
  PropagationDirection = cms.string( "alongMomentum" ),
  useRungeKutta = cms.bool( False ),
  ptMin = cms.double( 0.1 )
)
process.PropagatorWithMaterialForLoopers = cms.ESProducer( "PropagatorWithMaterialESProducer",
  SimpleMagneticField = cms.string( "ParabolicMf" ),
  MaxDPhi = cms.double( 4.0 ),
  ComponentName = cms.string( "PropagatorWithMaterialForLoopers" ),
  Mass = cms.double( 0.1396 ),
  PropagationDirection = cms.string( "alongMomentum" ),
  useRungeKutta = cms.bool( False ),
  ptMin = cms.double( -1.0 )
)
process.PixelCPEFastESProducer = cms.ESProducer( "PixelCPEFastESProducer",
  LoadTemplatesFromDB = cms.bool( True ),
  Alpha2Order = cms.bool( True ),
  ClusterProbComputationFlag = cms.int32( 0 ),
  useLAWidthFromDB = cms.bool( True ),
  lAOffset = cms.double( 0.0 ),
  lAWidthBPix = cms.double( 0.0 ),
  lAWidthFPix = cms.double( 0.0 ),
  doLorentzFromAlignment = cms.bool( False ),
  useLAFromDB = cms.bool( True ),
  xerr_barrel_l1 = cms.vdouble( 0.00115, 0.0012, 8.8E-4 ),
  yerr_barrel_l1 = cms.vdouble( 0.00375, 0.0023, 0.0025, 0.0025, 0.0023, 0.0023, 0.0021, 0.0021, 0.0024 ),
  xerr_barrel_ln = cms.vdouble( 0.00115, 0.0012, 8.8E-4 ),
  yerr_barrel_ln = cms.vdouble( 0.00375, 0.0023, 0.0025, 0.0025, 0.0023, 0.0023, 0.0021, 0.0021, 0.0024 ),
  xerr_endcap = cms.vdouble( 0.002, 0.002 ),
  yerr_endcap = cms.vdouble( 0.0021 ),
  xerr_barrel_l1_def = cms.double( 0.0103 ),
  yerr_barrel_l1_def = cms.double( 0.0021 ),
  xerr_barrel_ln_def = cms.double( 0.0103 ),
  yerr_barrel_ln_def = cms.double( 0.0021 ),
  xerr_endcap_def = cms.double( 0.002 ),
  yerr_endcap_def = cms.double( 7.5E-4 ),
  isPhase2 = cms.bool( False ),
  EdgeClusterErrorX = cms.double( 50.0 ),
  EdgeClusterErrorY = cms.double( 85.0 ),
  UseErrorsFromTemplates = cms.bool( True ),
  TruncatePixelCharge = cms.bool( True ),
  ComponentName = cms.string( "PixelCPEFast" ),
  MagneticFieldRecord = cms.ESInputTag( "","" ),
  appendToDataLabel = cms.string( "" )
)
process.ParametrizedMagneticFieldProducer = cms.ESProducer( "AutoParametrizedMagneticFieldProducer",
  version = cms.string( "Parabolic" ),
  label = cms.untracked.string( "ParabolicMf" ),
  valueOverride = cms.int32( -1 )
)
process.OppositePropagatorWithMaterialForMixedStep = cms.ESProducer( "PropagatorWithMaterialESProducer",
  SimpleMagneticField = cms.string( "ParabolicMf" ),
  MaxDPhi = cms.double( 1.6 ),
  ComponentName = cms.string( "PropagatorWithMaterialForMixedStepOpposite" ),
  Mass = cms.double( 0.105 ),
  PropagationDirection = cms.string( "oppositeToMomentum" ),
  useRungeKutta = cms.bool( False ),
  ptMin = cms.double( 0.1 )
)
process.OppositeMaterialPropagatorParabolicMF = cms.ESProducer( "PropagatorWithMaterialESProducer",
  SimpleMagneticField = cms.string( "ParabolicMf" ),
  MaxDPhi = cms.double( 1.6 ),
  ComponentName = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  Mass = cms.double( 0.105 ),
  PropagationDirection = cms.string( "oppositeToMomentum" ),
  useRungeKutta = cms.bool( False ),
  ptMin = cms.double( -1.0 )
)
process.OppositeMaterialPropagatorForHI = cms.ESProducer( "PropagatorWithMaterialESProducer",
  SimpleMagneticField = cms.string( "ParabolicMf" ),
  MaxDPhi = cms.double( 1.6 ),
  ComponentName = cms.string( "PropagatorWithMaterialOppositeForHI" ),
  Mass = cms.double( 0.139 ),
  PropagationDirection = cms.string( "oppositeToMomentum" ),
  useRungeKutta = cms.bool( False ),
  ptMin = cms.double( -1.0 )
)
process.OppositeMaterialPropagator = cms.ESProducer( "PropagatorWithMaterialESProducer",
  SimpleMagneticField = cms.string( "" ),
  MaxDPhi = cms.double( 1.6 ),
  ComponentName = cms.string( "PropagatorWithMaterialOpposite" ),
  Mass = cms.double( 0.105 ),
  PropagationDirection = cms.string( "oppositeToMomentum" ),
  useRungeKutta = cms.bool( False ),
  ptMin = cms.double( -1.0 )
)
process.MaterialPropagatorParabolicMF = cms.ESProducer( "PropagatorWithMaterialESProducer",
  SimpleMagneticField = cms.string( "ParabolicMf" ),
  MaxDPhi = cms.double( 1.6 ),
  ComponentName = cms.string( "PropagatorWithMaterialParabolicMf" ),
  Mass = cms.double( 0.105 ),
  PropagationDirection = cms.string( "alongMomentum" ),
  useRungeKutta = cms.bool( False ),
  ptMin = cms.double( -1.0 )
)
process.MaterialPropagatorForHI = cms.ESProducer( "PropagatorWithMaterialESProducer",
  SimpleMagneticField = cms.string( "ParabolicMf" ),
  MaxDPhi = cms.double( 1.6 ),
  ComponentName = cms.string( "PropagatorWithMaterialForHI" ),
  Mass = cms.double( 0.139 ),
  PropagationDirection = cms.string( "alongMomentum" ),
  useRungeKutta = cms.bool( False ),
  ptMin = cms.double( -1.0 )
)
process.MaterialPropagator = cms.ESProducer( "PropagatorWithMaterialESProducer",
  SimpleMagneticField = cms.string( "" ),
  MaxDPhi = cms.double( 1.6 ),
  ComponentName = cms.string( "PropagatorWithMaterial" ),
  Mass = cms.double( 0.105 ),
  PropagationDirection = cms.string( "alongMomentum" ),
  useRungeKutta = cms.bool( False ),
  ptMin = cms.double( -1.0 )
)
process.HcalTopologyIdealEP = cms.ESProducer( "HcalTopologyIdealEP",
  Exclude = cms.untracked.string( "" ),
  MergePosition = cms.untracked.bool( True ),
  appendToDataLabel = cms.string( "" )
)
process.HcalGeometryFromDBEP = cms.ESProducer( "HcalGeometryFromDBEP",
  applyAlignment = cms.bool( False )
)
process.GlobalParameters = cms.ESProducer( "StableParametersTrivialProducer",
  TotalBxInEvent = cms.int32( 5 ),
  NumberPhysTriggers = cms.uint32( 512 ),
  NumberL1Muon = cms.uint32( 8 ),
  NumberL1EGamma = cms.uint32( 12 ),
  NumberL1Jet = cms.uint32( 12 ),
  NumberL1Tau = cms.uint32( 12 ),
  NumberChips = cms.uint32( 1 ),
  PinsOnChip = cms.uint32( 512 ),
  OrderOfChip = cms.vint32( 1 ),
  NumberL1IsoEG = cms.uint32( 4 ),
  NumberL1JetCounts = cms.uint32( 12 ),
  UnitLength = cms.int32( 8 ),
  NumberL1ForJet = cms.uint32( 4 ),
  IfCaloEtaNumberBits = cms.uint32( 4 ),
  IfMuEtaNumberBits = cms.uint32( 6 ),
  NumberL1TauJet = cms.uint32( 4 ),
  NumberL1Mu = cms.uint32( 4 ),
  NumberConditionChips = cms.uint32( 1 ),
  NumberPsbBoards = cms.int32( 7 ),
  NumberL1CenJet = cms.uint32( 4 ),
  PinsOnConditionChip = cms.uint32( 512 ),
  NumberL1NoIsoEG = cms.uint32( 4 ),
  NumberTechnicalTriggers = cms.uint32( 64 ),
  NumberPhysTriggersExtended = cms.uint32( 64 ),
  WordLength = cms.int32( 64 ),
  OrderConditionChip = cms.vint32( 1 ),
  appendToDataLabel = cms.string( "" )
)
process.GEMGeometryESModule = cms.ESProducer( "GEMGeometryESModule",
  fromDDD = cms.bool( False ),
  fromDD4hep = cms.bool( False ),
  applyAlignment = cms.bool( False ),
  alignmentsLabel = cms.string( "" ),
  appendToDataLabel = cms.string( "" )
)
process.EcalPreshowerGeometryFromDBEP = cms.ESProducer( "EcalPreshowerGeometryFromDBEP",
  applyAlignment = cms.bool( True )
)
process.EcalLaserCorrectionService = cms.ESProducer( "EcalLaserCorrectionService",
  maxExtrapolationTimeInSec = cms.uint32( 0 ),
  appendToDataLabel = cms.string( "" )
)
process.EcalEndcapGeometryFromDBEP = cms.ESProducer( "EcalEndcapGeometryFromDBEP",
  applyAlignment = cms.bool( True )
)
process.EcalElectronicsMappingBuilder = cms.ESProducer( "EcalElectronicsMappingBuilder" )
process.EcalBarrelGeometryFromDBEP = cms.ESProducer( "EcalBarrelGeometryFromDBEP",
  applyAlignment = cms.bool( True )
)
process.DTObjectMapESProducer = cms.ESProducer( "DTObjectMapESProducer",
  appendToDataLabel = cms.string( "" )
)
process.DTGeometryESModule = cms.ESProducer( "DTGeometryESModule",
  fromDDD = cms.bool( False ),
  fromDD4hep = cms.bool( False ),
  DDDetector = cms.ESInputTag( "","" ),
  alignmentsLabel = cms.string( "" ),
  appendToDataLabel = cms.string( "" ),
  attribute = cms.string( "MuStructure" ),
  value = cms.string( "MuonBarrelDT" ),
  applyAlignment = cms.bool( True )
)
process.ClusterShapeHitFilterESProducer = cms.ESProducer( "ClusterShapeHitFilterESProducer",
  PixelShapeFile = cms.string( "RecoPixelVertexing/PixelLowPtUtilities/data/pixelShapePhase1_noL1.par" ),
  PixelShapeFileL1 = cms.string( "RecoPixelVertexing/PixelLowPtUtilities/data/pixelShapePhase1_loose.par" ),
  ComponentName = cms.string( "ClusterShapeHitFilter" ),
  isPhase2 = cms.bool( False ),
  doPixelShapeCut = cms.bool( True ),
  doStripShapeCut = cms.bool( True ),
  clusterChargeCut = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
  appendToDataLabel = cms.string( "" )
)
process.CastorDbProducer = cms.ESProducer( "CastorDbProducer",
  appendToDataLabel = cms.string( "" )
)
process.CaloTowerTopologyEP = cms.ESProducer( "CaloTowerTopologyEP",
  appendToDataLabel = cms.string( "" )
)
process.CaloTowerGeometryFromDBEP = cms.ESProducer( "CaloTowerGeometryFromDBEP",
  applyAlignment = cms.bool( False )
)
process.CaloTowerConstituentsMapBuilder = cms.ESProducer( "CaloTowerConstituentsMapBuilder",
  MapFile = cms.untracked.string( "Geometry/CaloTopology/data/CaloTowerEEGeometric.map.gz" ),
  MapAuto = cms.untracked.bool( False ),
  SkipHE = cms.untracked.bool( False ),
  appendToDataLabel = cms.string( "" )
)
process.CaloTopologyBuilder = cms.ESProducer( "CaloTopologyBuilder" )
process.CaloGeometryBuilder = cms.ESProducer( "CaloGeometryBuilder",
  SelectedCalos = cms.vstring( 'HCAL',
    'ZDC',
    'EcalBarrel',
    'EcalEndcap',
    'EcalPreshower',
    'TOWER' )
)
process.CSCObjectMapESProducer = cms.ESProducer( "CSCObjectMapESProducer",
  appendToDataLabel = cms.string( "" )
)
process.CSCIndexerESProducer = cms.ESProducer( "CSCIndexerESProducer",
  AlgoName = cms.string( "CSCIndexerPostls1" )
)
process.CSCGeometryESModule = cms.ESProducer( "CSCGeometryESModule",
  fromDDD = cms.bool( False ),
  fromDD4hep = cms.bool( False ),
  alignmentsLabel = cms.string( "" ),
  appendToDataLabel = cms.string( "" ),
  useRealWireGeometry = cms.bool( True ),
  useOnlyWiresInME1a = cms.bool( False ),
  useGangedStripsInME1a = cms.bool( False ),
  useCentreTIOffsets = cms.bool( False ),
  applyAlignment = cms.bool( True ),
  debugV = cms.untracked.bool( False )
)
process.CSCChannelMapperESProducer = cms.ESProducer( "CSCChannelMapperESProducer",
  AlgoName = cms.string( "CSCChannelMapperPostls1" )
)
process.AnyDirectionAnalyticalPropagator = cms.ESProducer( "AnalyticalPropagatorESProducer",
  MaxDPhi = cms.double( 1.6 ),
  ComponentName = cms.string( "AnyDirectionAnalyticalPropagator" ),
  PropagationDirection = cms.string( "anyDirection" )
)

process.ThroughputService = cms.Service( "ThroughputService",
    eventRange = cms.untracked.uint32( 10000 ),
    eventResolution = cms.untracked.uint32( 1 ),
    printEventSummary = cms.untracked.bool( False ),
    enableDQM = cms.untracked.bool( True ),
    dqmPathByProcesses = cms.untracked.bool( False ),
    dqmPath = cms.untracked.string( "HLT/Throughput" ),
    timeRange = cms.untracked.double( 60000.0 ),
    timeResolution = cms.untracked.double( 5.828 )
)
process.MessageLogger = cms.Service( "MessageLogger",
    suppressWarning = cms.untracked.vstring( 'hltOnlineBeamSpot',
      'hltCtf3HitL1SeededWithMaterialTracks',
      'hltL3MuonsOIState',
      'hltPixelTracksForHighMult',
      'hltHITPixelTracksHE',
      'hltHITPixelTracksHB',
      'hltCtfL1SeededWithMaterialTracks',
      'hltRegionalTracksForL3MuonIsolation',
      'hltSiPixelClusters',
      'hltActivityStartUpElectronPixelSeeds',
      'hltLightPFTracks',
      'hltPixelVertices3DbbPhi',
      'hltL3MuonsIOHit',
      'hltPixelTracks',
      'hltSiPixelDigis',
      'hltL3MuonsOIHit',
      'hltL1SeededElectronGsfTracks',
      'hltL1SeededStartUpElectronPixelSeeds',
      'hltBLifetimeRegionalCtfWithMaterialTracksbbPhiL1FastJetFastPV',
      'hltCtfActivityWithMaterialTracks' ),
    suppressFwkInfo = cms.untracked.vstring(  ),
    suppressInfo = cms.untracked.vstring(  ),
    suppressDebug = cms.untracked.vstring(  ),
    debugModules = cms.untracked.vstring(  ),
    cerr = cms.untracked.PSet( 
      INFO = cms.untracked.PSet(  limit = cms.untracked.int32( 0 ) ),
      noTimeStamps = cms.untracked.bool( False ),
      FwkReport = cms.untracked.PSet( 
        reportEvery = cms.untracked.int32( 1 ),
        limit = cms.untracked.int32( 0 )
      ),
      default = cms.untracked.PSet(  limit = cms.untracked.int32( 10000000 ) ),
      Root_NoDictionary = cms.untracked.PSet(  limit = cms.untracked.int32( 0 ) ),
      FwkJob = cms.untracked.PSet(  limit = cms.untracked.int32( 0 ) ),
      FwkSummary = cms.untracked.PSet( 
        reportEvery = cms.untracked.int32( 1 ),
        limit = cms.untracked.int32( 10000000 )
      ),
      threshold = cms.untracked.string( "INFO" ),
    ),
    suppressError = cms.untracked.vstring( 'hltOnlineBeamSpot',
      'hltL3MuonCandidates',
      'hltL3TkTracksFromL2OIState',
      'hltPFJetCtfWithMaterialTracks',
      'hltL3TkTracksFromL2IOHit',
      'hltL3TkTracksFromL2OIHit' )
)
process.FastTimerService = cms.Service( "FastTimerService",
    printEventSummary = cms.untracked.bool( False ),
    printRunSummary = cms.untracked.bool( True ),
    printJobSummary = cms.untracked.bool( True ),
    writeJSONSummary = cms.untracked.bool( False ),
    jsonFileName = cms.untracked.string( "resources.json" ),
    enableDQM = cms.untracked.bool( True ),
    enableDQMbyModule = cms.untracked.bool( False ),
    enableDQMbyPath = cms.untracked.bool( False ),
    enableDQMbyLumiSection = cms.untracked.bool( True ),
    enableDQMbyProcesses = cms.untracked.bool( True ),
    enableDQMTransitions = cms.untracked.bool( False ),
    dqmTimeRange = cms.untracked.double( 2000.0 ),
    dqmTimeResolution = cms.untracked.double( 5.0 ),
    dqmMemoryRange = cms.untracked.double( 1000000.0 ),
    dqmMemoryResolution = cms.untracked.double( 5000.0 ),
    dqmPathTimeRange = cms.untracked.double( 100.0 ),
    dqmPathTimeResolution = cms.untracked.double( 0.5 ),
    dqmPathMemoryRange = cms.untracked.double( 1000000.0 ),
    dqmPathMemoryResolution = cms.untracked.double( 5000.0 ),
    dqmModuleTimeRange = cms.untracked.double( 40.0 ),
    dqmModuleTimeResolution = cms.untracked.double( 0.2 ),
    dqmModuleMemoryRange = cms.untracked.double( 100000.0 ),
    dqmModuleMemoryResolution = cms.untracked.double( 500.0 ),
    dqmLumiSectionsRange = cms.untracked.uint32( 2500 ),
    dqmPath = cms.untracked.string( "HLT/TimerService" ),
)
process.CUDAService = cms.Service( "CUDAService",
    enabled = cms.untracked.bool( True ),
    verbose = cms.untracked.bool( False ),
    limits = cms.untracked.PSet( 
      cudaLimitDevRuntimePendingLaunchCount = cms.untracked.int32( -1 ),
      cudaLimitDevRuntimeSyncDepth = cms.untracked.int32( -1 ),
      cudaLimitStackSize = cms.untracked.int32( -1 ),
      cudaLimitPrintfFifoSize = cms.untracked.int32( -1 ),
      cudaLimitMallocHeapSize = cms.untracked.int32( -1 )
    ),
    allocator = cms.untracked.PSet( 
      hostPreallocate = cms.untracked.vuint32(  ),
      devicePreallocate = cms.untracked.vuint32(  )
    )
)

process.hltGetConditions = cms.EDAnalyzer( "EventSetupRecordDataGetter",
    verbose = cms.untracked.bool( False ),
    toGet = cms.VPSet( 
    )
)
process.hltGetRaw = cms.EDAnalyzer( "HLTGetRaw",
    RawDataCollection = cms.InputTag( "rawDataCollector" )
)
process.hltPSetMap = cms.EDProducer( "ParameterSetBlobProducer" )
process.hltBoolFalse = cms.EDFilter( "HLTBool",
    result = cms.bool( False )
)
process.hltTriggerType = cms.EDFilter( "HLTTriggerTypeFilter",
    SelectedTriggerType = cms.int32( 1 )
)
process.hltGtStage2Digis = cms.EDProducer( "L1TRawToDigi",
    FedIds = cms.vint32( 1404 ),
    Setup = cms.string( "stage2::GTSetup" ),
    FWId = cms.uint32( 0 ),
    DmxFWId = cms.uint32( 0 ),
    FWOverride = cms.bool( False ),
    TMTCheck = cms.bool( True ),
    CTP7 = cms.untracked.bool( False ),
    MTF7 = cms.untracked.bool( False ),
    InputLabel = cms.InputTag( "rawDataCollector" ),
    lenSlinkHeader = cms.untracked.int32( 8 ),
    lenSlinkTrailer = cms.untracked.int32( 8 ),
    lenAMCHeader = cms.untracked.int32( 8 ),
    lenAMCTrailer = cms.untracked.int32( 0 ),
    lenAMC13Header = cms.untracked.int32( 8 ),
    lenAMC13Trailer = cms.untracked.int32( 8 ),
    debug = cms.untracked.bool( False ),
    MinFeds = cms.uint32( 0 )
)
process.hltGtStage2ObjectMap = cms.EDProducer( "L1TGlobalProducer",
    MuonInputTag = cms.InputTag( 'hltGtStage2Digis','Muon' ),
    MuonShowerInputTag = cms.InputTag( 'hltGtStage2Digis','MuonShower' ),
    EGammaInputTag = cms.InputTag( 'hltGtStage2Digis','EGamma' ),
    TauInputTag = cms.InputTag( 'hltGtStage2Digis','Tau' ),
    JetInputTag = cms.InputTag( 'hltGtStage2Digis','Jet' ),
    EtSumInputTag = cms.InputTag( 'hltGtStage2Digis','EtSum' ),
    ExtInputTag = cms.InputTag( "hltGtStage2Digis" ),
    AlgoBlkInputTag = cms.InputTag( "hltGtStage2Digis" ),
    GetPrescaleColumnFromData = cms.bool( False ),
    AlgorithmTriggersUnprescaled = cms.bool( True ),
    RequireMenuToMatchAlgoBlkInput = cms.bool( True ),
    AlgorithmTriggersUnmasked = cms.bool( True ),
    useMuonShowers = cms.bool( True ),
    resetPSCountersEachLumiSec = cms.bool( True ),
    ProduceL1GtDaqRecord = cms.bool( True ),
    ProduceL1GtObjectMapRecord = cms.bool( True ),
    EmulateBxInEvent = cms.int32( 1 ),
    L1DataBxInEvent = cms.int32( 5 ),
    AlternativeNrBxBoardDaq = cms.uint32( 0 ),
    BstLengthBytes = cms.int32( -1 ),
    PrescaleSet = cms.uint32( 1 ),
    Verbosity = cms.untracked.int32( 0 ),
    PrintL1Menu = cms.untracked.bool( False ),
    TriggerMenuLuminosity = cms.string( "startup" )
)
process.hltScalersRawToDigi = cms.EDProducer( "ScalersRawToDigi",
    scalersInputTag = cms.InputTag( "rawDataCollector" )
)
process.hltOnlineMetaDataDigis = cms.EDProducer( "OnlineMetaDataRawToDigi",
    onlineMetaDataInputLabel = cms.InputTag( "rawDataCollector" )
)
process.hltOnlineBeamSpot = cms.EDProducer( "BeamSpotOnlineProducer",
    changeToCMSCoordinates = cms.bool( False ),
    maxZ = cms.double( 40.0 ),
    setSigmaZ = cms.double( 0.0 ),
    beamMode = cms.untracked.uint32( 11 ),
    src = cms.InputTag( "hltScalersRawToDigi" ),
    gtEvmLabel = cms.InputTag( "" ),
    maxRadius = cms.double( 2.0 ),
    useTransientRecord = cms.bool( True )
)
process.hltL1sDSTRun3PFScoutingPixelTracking = cms.EDFilter( "HLTL1TSeed",
    saveTags = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_DoubleMu_12_5 OR L1_DoubleMu_15_7 OR L1_HTT200er OR L1_HTT255er OR L1_HTT280er OR L1_HTT320er OR L1_HTT360er OR L1_ETT2000 OR L1_HTT400er OR L1_HTT450er OR L1_SingleJet180 OR L1_SingleJet200 OR L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5 OR L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5 OR L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5 OR L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7 OR L1_DoubleMu4_SQ_OS_dR_Max1p2 OR L1_SingleEG36er2p5 OR L1_SingleLooseIsoEG28er2p1" ),
    L1ObjectMapInputTag = cms.InputTag( "hltGtStage2ObjectMap" ),
    L1GlobalInputTag = cms.InputTag( "hltGtStage2Digis" ),
    L1MuonInputTag = cms.InputTag( 'hltGtStage2Digis','Muon' ),
    L1MuonShowerInputTag = cms.InputTag( 'hltGtStage2Digis','MuonShower' ),
    L1EGammaInputTag = cms.InputTag( 'hltGtStage2Digis','EGamma' ),
    L1JetInputTag = cms.InputTag( 'hltGtStage2Digis','Jet' ),
    L1TauInputTag = cms.InputTag( 'hltGtStage2Digis','Tau' ),
    L1EtSumInputTag = cms.InputTag( 'hltGtStage2Digis','EtSum' )
)
process.hltPreDSTRun3PFScoutingFullTracking = cms.EDFilter( "HLTPrescaler",
    offset = cms.uint32( 0 ),
    L1GtReadoutRecordTag = cms.InputTag( "hltGtStage2Digis" )
)
process.hltEcalDigisLegacy = cms.EDProducer( "EcalRawToDigi",
    tccUnpacking = cms.bool( True ),
    FedLabel = cms.InputTag( "listfeds" ),
    srpUnpacking = cms.bool( True ),
    syncCheck = cms.bool( True ),
    feIdCheck = cms.bool( True ),
    silentMode = cms.untracked.bool( True ),
    InputLabel = cms.InputTag( "rawDataCollector" ),
    orderedFedList = cms.vint32( 601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 612, 613, 614, 615, 616, 617, 618, 619, 620, 621, 622, 623, 624, 625, 626, 627, 628, 629, 630, 631, 632, 633, 634, 635, 636, 637, 638, 639, 640, 641, 642, 643, 644, 645, 646, 647, 648, 649, 650, 651, 652, 653, 654 ),
    eventPut = cms.bool( True ),
    numbTriggerTSamples = cms.int32( 1 ),
    numbXtalTSamples = cms.int32( 10 ),
    orderedDCCIdList = cms.vint32( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54 ),
    FEDs = cms.vint32( 601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 612, 613, 614, 615, 616, 617, 618, 619, 620, 621, 622, 623, 624, 625, 626, 627, 628, 629, 630, 631, 632, 633, 634, 635, 636, 637, 638, 639, 640, 641, 642, 643, 644, 645, 646, 647, 648, 649, 650, 651, 652, 653, 654 ),
    DoRegional = cms.bool( False ),
    feUnpacking = cms.bool( True ),
    forceToKeepFRData = cms.bool( False ),
    headerUnpacking = cms.bool( True ),
    memUnpacking = cms.bool( True )
)
process.hltEcalDetIdToBeRecovered = cms.EDProducer( "EcalDetIdToBeRecoveredProducer",
    ebIntegrityChIdErrors = cms.InputTag( 'hltEcalDigis','EcalIntegrityChIdErrors' ),
    ebDetIdToBeRecovered = cms.string( "ebDetId" ),
    integrityTTIdErrors = cms.InputTag( 'hltEcalDigis','EcalIntegrityTTIdErrors' ),
    eeIntegrityGainErrors = cms.InputTag( 'hltEcalDigis','EcalIntegrityGainErrors' ),
    ebFEToBeRecovered = cms.string( "ebFE" ),
    ebIntegrityGainErrors = cms.InputTag( 'hltEcalDigis','EcalIntegrityGainErrors' ),
    eeDetIdToBeRecovered = cms.string( "eeDetId" ),
    eeIntegrityGainSwitchErrors = cms.InputTag( 'hltEcalDigis','EcalIntegrityGainSwitchErrors' ),
    eeIntegrityChIdErrors = cms.InputTag( 'hltEcalDigis','EcalIntegrityChIdErrors' ),
    ebIntegrityGainSwitchErrors = cms.InputTag( 'hltEcalDigis','EcalIntegrityGainSwitchErrors' ),
    ebSrFlagCollection = cms.InputTag( "hltEcalDigis" ),
    eeFEToBeRecovered = cms.string( "eeFE" ),
    integrityBlockSizeErrors = cms.InputTag( 'hltEcalDigis','EcalIntegrityBlockSizeErrors' ),
    eeSrFlagCollection = cms.InputTag( "hltEcalDigis" )
)
process.hltHcalDigis = cms.EDProducer( "HcalRawToDigi",
    HcalFirstFED = cms.untracked.int32( 700 ),
    firstSample = cms.int32( 0 ),
    lastSample = cms.int32( 9 ),
    FilterDataQuality = cms.bool( True ),
    FEDs = cms.untracked.vint32(  ),
    UnpackZDC = cms.untracked.bool( True ),
    UnpackCalib = cms.untracked.bool( True ),
    UnpackUMNio = cms.untracked.bool( True ),
    UnpackTTP = cms.untracked.bool( False ),
    silent = cms.untracked.bool( True ),
    saveQIE10DataNSamples = cms.untracked.vint32(  ),
    saveQIE10DataTags = cms.untracked.vstring(  ),
    saveQIE11DataNSamples = cms.untracked.vint32(  ),
    saveQIE11DataTags = cms.untracked.vstring(  ),
    ComplainEmptyData = cms.untracked.bool( False ),
    UnpackerMode = cms.untracked.int32( 0 ),
    ExpectedOrbitMessageTime = cms.untracked.int32( -1 ),
    InputLabel = cms.InputTag( "rawDataCollector" ),
    ElectronicsMap = cms.string( "" )
)
process.hltHfprereco = cms.EDProducer( "HFPreReconstructor",
    digiLabel = cms.InputTag( "hltHcalDigis" ),
    dropZSmarkedPassed = cms.bool( True ),
    tsFromDB = cms.bool( False ),
    sumAllTimeSlices = cms.bool( False ),
    forceSOI = cms.int32( -1 ),
    soiShift = cms.int32( 0 )
)
process.hltHfreco = cms.EDProducer( "HFPhase1Reconstructor",
    inputLabel = cms.InputTag( "hltHfprereco" ),
    useChannelQualityFromDB = cms.bool( False ),
    checkChannelQualityForDepth3and4 = cms.bool( False ),
    algorithm = cms.PSet( 
      tfallIfNoTDC = cms.double( -101.0 ),
      triseIfNoTDC = cms.double( -100.0 ),
      rejectAllFailures = cms.bool( True ),
      energyWeights = cms.vdouble( 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 2.0, 0.0, 2.0, 0.0, 2.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 2.0, 0.0, 2.0, 0.0, 2.0, 0.0, 1.0 ),
      soiPhase = cms.uint32( 1 ),
      timeShift = cms.double( 0.0 ),
      tlimits = cms.vdouble( -1000.0, 1000.0, -1000.0, 1000.0 ),
      Class = cms.string( "HFFlexibleTimeCheck" )
    ),
    algoConfigClass = cms.string( "HFPhase1PMTParams" ),
    setNoiseFlags = cms.bool( True ),
    runHFStripFilter = cms.bool( False ),
    S9S1stat = cms.PSet( 
      shortEnergyParams = cms.vdouble( 35.1773, 35.37, 35.7933, 36.4472, 37.3317, 38.4468, 39.7925, 41.3688, 43.1757, 45.2132, 47.4813, 49.98, 52.7093 ),
      shortETParams = cms.vdouble( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ),
      long_optimumSlope = cms.vdouble( -99999.0, 0.0164905, 0.0238698, 0.0321383, 0.041296, 0.0513428, 0.0622789, 0.0741041, 0.0868186, 0.100422, 0.135313, 0.136289, 0.0589927 ),
      isS8S1 = cms.bool( False ),
      longETParams = cms.vdouble( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ),
      longEnergyParams = cms.vdouble( 43.5, 45.7, 48.32, 51.36, 54.82, 58.7, 63.0, 67.72, 72.86, 78.42, 84.4, 90.8, 97.62 ),
      short_optimumSlope = cms.vdouble( -99999.0, 0.0164905, 0.0238698, 0.0321383, 0.041296, 0.0513428, 0.0622789, 0.0741041, 0.0868186, 0.100422, 0.135313, 0.136289, 0.0589927 ),
      HcalAcceptSeverityLevel = cms.int32( 9 )
    ),
    S8S1stat = cms.PSet( 
      shortEnergyParams = cms.vdouble( 40.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0 ),
      shortETParams = cms.vdouble( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ),
      long_optimumSlope = cms.vdouble( 0.3, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 ),
      isS8S1 = cms.bool( True ),
      longETParams = cms.vdouble( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ),
      longEnergyParams = cms.vdouble( 40.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0 ),
      short_optimumSlope = cms.vdouble( 0.3, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 ),
      HcalAcceptSeverityLevel = cms.int32( 9 )
    ),
    PETstat = cms.PSet( 
      shortEnergyParams = cms.vdouble( 35.1773, 35.37, 35.7933, 36.4472, 37.3317, 38.4468, 39.7925, 41.3688, 43.1757, 45.2132, 47.4813, 49.98, 52.7093 ),
      shortETParams = cms.vdouble( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ),
      long_R_29 = cms.vdouble( 0.8 ),
      longETParams = cms.vdouble( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ),
      longEnergyParams = cms.vdouble( 43.5, 45.7, 48.32, 51.36, 54.82, 58.7, 63.0, 67.72, 72.86, 78.42, 84.4, 90.8, 97.62 ),
      short_R_29 = cms.vdouble( 0.8 ),
      long_R = cms.vdouble( 0.98 ),
      short_R = cms.vdouble( 0.8 ),
      HcalAcceptSeverityLevel = cms.int32( 9 )
    ),
    HFStripFilter = cms.PSet( 
      timeMax = cms.double( 6.0 ),
      seedHitIetaMax = cms.int32( 35 ),
      gap = cms.int32( 2 ),
      verboseLevel = cms.untracked.int32( 10 ),
      wedgeCut = cms.double( 0.05 ),
      stripThreshold = cms.double( 40.0 ),
      maxStripTime = cms.double( 10.0 ),
      maxThreshold = cms.double( 100.0 ),
      lstrips = cms.int32( 2 )
    )
)
process.hltHoreco = cms.EDProducer( "HcalHitReconstructor",
    correctForPhaseContainment = cms.bool( True ),
    correctionPhaseNS = cms.double( 13.0 ),
    digiLabel = cms.InputTag( "hltHcalDigis" ),
    Subdetector = cms.string( "HO" ),
    correctForTimeslew = cms.bool( True ),
    dropZSmarkedPassed = cms.bool( True ),
    firstSample = cms.int32( 4 ),
    samplesToAdd = cms.int32( 4 ),
    tsFromDB = cms.bool( True ),
    recoParamsFromDB = cms.bool( True ),
    useLeakCorrection = cms.bool( False ),
    dataOOTCorrectionName = cms.string( "" ),
    dataOOTCorrectionCategory = cms.string( "Data" ),
    mcOOTCorrectionName = cms.string( "" ),
    mcOOTCorrectionCategory = cms.string( "MC" ),
    correctTiming = cms.bool( False ),
    firstAuxTS = cms.int32( 4 ),
    setNoiseFlags = cms.bool( False ),
    digiTimeFromDB = cms.bool( True ),
    setHSCPFlags = cms.bool( False ),
    setSaturationFlags = cms.bool( False ),
    setTimingTrustFlags = cms.bool( False ),
    setPulseShapeFlags = cms.bool( False ),
    setNegativeFlags = cms.bool( False ),
    digistat = cms.PSet(  ),
    HFInWindowStat = cms.PSet(  ),
    S9S1stat = cms.PSet(  ),
    S8S1stat = cms.PSet(  ),
    PETstat = cms.PSet(  ),
    saturationParameters = cms.PSet(  maxADCvalue = cms.int32( 127 ) ),
    hfTimingTrustParameters = cms.PSet(  )
)
process.hltTowerMakerForAll = cms.EDProducer( "CaloTowersCreator",
    EBSumThreshold = cms.double( 0.2 ),
    HF2Weight = cms.double( 1.0 ),
    EBWeight = cms.double( 1.0 ),
    hfInput = cms.InputTag( "hltHfreco" ),
    EESumThreshold = cms.double( 0.45 ),
    HOThreshold0 = cms.double( 3.5 ),
    HOThresholdPlus1 = cms.double( 3.5 ),
    HOThresholdMinus1 = cms.double( 3.5 ),
    HOThresholdPlus2 = cms.double( 3.5 ),
    HOThresholdMinus2 = cms.double( 3.5 ),
    HBGrid = cms.vdouble(  ),
    HBThreshold1 = cms.double( 0.1 ),
    HBThreshold2 = cms.double( 0.2 ),
    HBThreshold = cms.double( 0.3 ),
    EEWeights = cms.vdouble(  ),
    HF1Threshold = cms.double( 0.5 ),
    HF2Weights = cms.vdouble(  ),
    HOWeights = cms.vdouble(  ),
    EEGrid = cms.vdouble(  ),
    HEDWeight = cms.double( 1.0 ),
    EEWeight = cms.double( 1.0 ),
    UseHO = cms.bool( False ),
    HBWeights = cms.vdouble(  ),
    HESWeight = cms.double( 1.0 ),
    HF1Weight = cms.double( 1.0 ),
    HF2Grid = cms.vdouble(  ),
    HEDWeights = cms.vdouble(  ),
    HF1Grid = cms.vdouble(  ),
    EBWeights = cms.vdouble(  ),
    HOWeight = cms.double( 1.0E-99 ),
    EBThreshold = cms.double( 0.07 ),
    EEThreshold = cms.double( 0.3 ),
    UseEtEBTreshold = cms.bool( False ),
    UseSymEBTreshold = cms.bool( False ),
    UseEtEETreshold = cms.bool( False ),
    UseSymEETreshold = cms.bool( False ),
    hbheInput = cms.InputTag( "hltHbhereco" ),
    HcalThreshold = cms.double( -1000.0 ),
    HF2Threshold = cms.double( 0.85 ),
    HESThreshold1 = cms.double( 0.1 ),
    HESThreshold = cms.double( 0.2 ),
    HF1Weights = cms.vdouble(  ),
    hoInput = cms.InputTag( "hltHoreco" ),
    HESGrid = cms.vdouble(  ),
    HESWeights = cms.vdouble(  ),
    HEDThreshold1 = cms.double( 0.1 ),
    HEDThreshold = cms.double( 0.2 ),
    EcutTower = cms.double( -1000.0 ),
    HEDGrid = cms.vdouble(  ),
    ecalInputs = cms.VInputTag( 'hltEcalRecHit:EcalRecHitsEB','hltEcalRecHit:EcalRecHitsEE' ),
    HBWeight = cms.double( 1.0 ),
    HOGrid = cms.vdouble(  ),
    EBGrid = cms.vdouble(  ),
    MomConstrMethod = cms.int32( 1 ),
    MomHBDepth = cms.double( 0.2 ),
    MomHEDepth = cms.double( 0.4 ),
    MomEBDepth = cms.double( 0.3 ),
    MomEEDepth = cms.double( 0.0 ),
    HcalAcceptSeverityLevel = cms.uint32( 9 ),
    EcalRecHitSeveritiesToBeExcluded = cms.vstring( 'kTime',
      'kWeird',
      'kBad' ),
    UseHcalRecoveredHits = cms.bool( False ),
    UseEcalRecoveredHits = cms.bool( False ),
    UseRejectedHitsOnly = cms.bool( False ),
    HcalAcceptSeverityLevelForRejectedHit = cms.uint32( 9999 ),
    EcalSeveritiesToBeUsedInBadTowers = cms.vstring(  ),
    UseRejectedRecoveredHcalHits = cms.bool( False ),
    UseRejectedRecoveredEcalHits = cms.bool( False ),
    missingHcalRescaleFactorForEcal = cms.double( 0.0 ),
    AllowMissingInputs = cms.bool( False ),
    HcalPhase = cms.int32( 1 )
)
process.hltAK4CaloJets = cms.EDProducer( "FastjetJetProducer",
    useMassDropTagger = cms.bool( False ),
    useFiltering = cms.bool( False ),
    useDynamicFiltering = cms.bool( False ),
    useTrimming = cms.bool( False ),
    usePruning = cms.bool( False ),
    useCMSBoostedTauSeedingAlgorithm = cms.bool( False ),
    useKtPruning = cms.bool( False ),
    useConstituentSubtraction = cms.bool( False ),
    useSoftDrop = cms.bool( False ),
    correctShape = cms.bool( False ),
    UseOnlyVertexTracks = cms.bool( False ),
    UseOnlyOnePV = cms.bool( False ),
    muCut = cms.double( -1.0 ),
    yCut = cms.double( -1.0 ),
    rFilt = cms.double( -1.0 ),
    rFiltFactor = cms.double( -1.0 ),
    trimPtFracMin = cms.double( -1.0 ),
    zcut = cms.double( -1.0 ),
    rcut_factor = cms.double( -1.0 ),
    csRho_EtaMax = cms.double( -1.0 ),
    csRParam = cms.double( -1.0 ),
    beta = cms.double( -1.0 ),
    R0 = cms.double( -1.0 ),
    gridMaxRapidity = cms.double( -1.0 ),
    gridSpacing = cms.double( -1.0 ),
    DzTrVtxMax = cms.double( 0.0 ),
    DxyTrVtxMax = cms.double( 0.0 ),
    MaxVtxZ = cms.double( 15.0 ),
    subjetPtMin = cms.double( -1.0 ),
    muMin = cms.double( -1.0 ),
    muMax = cms.double( -1.0 ),
    yMin = cms.double( -1.0 ),
    yMax = cms.double( -1.0 ),
    dRMin = cms.double( -1.0 ),
    dRMax = cms.double( -1.0 ),
    maxDepth = cms.int32( -1 ),
    nFilt = cms.int32( -1 ),
    MinVtxNdof = cms.int32( 5 ),
    src = cms.InputTag( "hltTowerMakerForAll" ),
    srcPVs = cms.InputTag( "NotUsed" ),
    jetType = cms.string( "CaloJet" ),
    jetAlgorithm = cms.string( "AntiKt" ),
    rParam = cms.double( 0.4 ),
    inputEtMin = cms.double( 0.3 ),
    inputEMin = cms.double( 0.0 ),
    jetPtMin = cms.double( 1.0 ),
    doPVCorrection = cms.bool( False ),
    doAreaFastjet = cms.bool( False ),
    doRhoFastjet = cms.bool( False ),
    doPUOffsetCorr = cms.bool( False ),
    puPtMin = cms.double( 10.0 ),
    nSigmaPU = cms.double( 1.0 ),
    radiusPU = cms.double( 0.4 ),
    subtractorName = cms.string( "" ),
    useExplicitGhosts = cms.bool( False ),
    doAreaDiskApprox = cms.bool( True ),
    voronoiRfact = cms.double( 0.9 ),
    Rho_EtaMax = cms.double( 4.4 ),
    Ghost_EtaMax = cms.double( 6.0 ),
    Active_Area_Repeats = cms.int32( 5 ),
    GhostArea = cms.double( 0.01 ),
    restrictInputs = cms.bool( False ),
    maxInputs = cms.uint32( 1 ),
    writeCompound = cms.bool( False ),
    writeJetsWithConst = cms.bool( False ),
    doFastJetNonUniform = cms.bool( False ),
    useDeterministicSeed = cms.bool( True ),
    minSeed = cms.uint32( 14327 ),
    verbosity = cms.int32( 0 ),
    puWidth = cms.double( 0.0 ),
    nExclude = cms.uint32( 0 ),
    maxBadEcalCells = cms.uint32( 9999999 ),
    maxBadHcalCells = cms.uint32( 9999999 ),
    maxProblematicEcalCells = cms.uint32( 9999999 ),
    maxProblematicHcalCells = cms.uint32( 9999999 ),
    maxRecoveredEcalCells = cms.uint32( 9999999 ),
    maxRecoveredHcalCells = cms.uint32( 9999999 ),
    puCenters = cms.vdouble(  ),
    applyWeight = cms.bool( False ),
    srcWeights = cms.InputTag( "" ),
    minimumTowersFraction = cms.double( 0.0 ),
    jetCollInstanceName = cms.string( "" ),
    sumRecHits = cms.bool( False )
)
process.hltAK4CaloJetsIDPassed = cms.EDProducer( "HLTCaloJetIDProducer",
    min_N90 = cms.int32( -2 ),
    min_N90hits = cms.int32( 2 ),
    min_EMF = cms.double( 1.0E-6 ),
    max_EMF = cms.double( 999.0 ),
    jetsInput = cms.InputTag( "hltAK4CaloJets" ),
    JetIDParams = cms.PSet( 
      hfRecHitsColl = cms.InputTag( "hltHfreco" ),
      hoRecHitsColl = cms.InputTag( "hltHoreco" ),
      ebRecHitsColl = cms.InputTag( 'hltEcalRecHit','EcalRecHitsEB' ),
      hbheRecHitsColl = cms.InputTag( "hltHbhereco" ),
      useRecHits = cms.bool( True ),
      eeRecHitsColl = cms.InputTag( 'hltEcalRecHit','EcalRecHitsEE' )
    )
)
process.hltFixedGridRhoFastjetAllCalo = cms.EDProducer( "FixedGridRhoProducerFastjet",
    pfCandidatesTag = cms.InputTag( "hltTowerMakerForAll" ),
    maxRapidity = cms.double( 5.0 ),
    gridSpacing = cms.double( 0.55 )
)
process.hltAK4CaloFastJetCorrector = cms.EDProducer( "L1FastjetCorrectorProducer",
    level = cms.string( "L1FastJet" ),
    algorithm = cms.string( "AK4CaloHLT" ),
    srcRho = cms.InputTag( "hltFixedGridRhoFastjetAllCalo" )
)
process.hltAK4CaloRelativeCorrector = cms.EDProducer( "LXXXCorrectorProducer",
    level = cms.string( "L2Relative" ),
    algorithm = cms.string( "AK4CaloHLT" )
)
process.hltAK4CaloAbsoluteCorrector = cms.EDProducer( "LXXXCorrectorProducer",
    level = cms.string( "L3Absolute" ),
    algorithm = cms.string( "AK4CaloHLT" )
)
process.hltAK4CaloResidualCorrector = cms.EDProducer( "LXXXCorrectorProducer",
    level = cms.string( "L2L3Residual" ),
    algorithm = cms.string( "AK4CaloHLT" )
)
process.hltAK4CaloCorrector = cms.EDProducer( "ChainedJetCorrectorProducer",
    correctors = cms.VInputTag( 'hltAK4CaloFastJetCorrector','hltAK4CaloRelativeCorrector','hltAK4CaloAbsoluteCorrector','hltAK4CaloResidualCorrector' )
)
process.hltAK4CaloJetsCorrected = cms.EDProducer( "CorrectedCaloJetProducer",
    src = cms.InputTag( "hltAK4CaloJets" ),
    correctors = cms.VInputTag( 'hltAK4CaloCorrector' )
)
process.hltAK4CaloJetsCorrectedIDPassed = cms.EDProducer( "CorrectedCaloJetProducer",
    src = cms.InputTag( "hltAK4CaloJetsIDPassed" ),
    correctors = cms.VInputTag( 'hltAK4CaloCorrector' )
)
process.hltAK4CaloJetsPF = cms.EDProducer( "FastjetJetProducer",
    useMassDropTagger = cms.bool( False ),
    useFiltering = cms.bool( False ),
    useDynamicFiltering = cms.bool( False ),
    useTrimming = cms.bool( False ),
    usePruning = cms.bool( False ),
    useCMSBoostedTauSeedingAlgorithm = cms.bool( False ),
    useKtPruning = cms.bool( False ),
    useConstituentSubtraction = cms.bool( False ),
    useSoftDrop = cms.bool( False ),
    correctShape = cms.bool( False ),
    UseOnlyVertexTracks = cms.bool( False ),
    UseOnlyOnePV = cms.bool( False ),
    muCut = cms.double( -1.0 ),
    yCut = cms.double( -1.0 ),
    rFilt = cms.double( -1.0 ),
    rFiltFactor = cms.double( -1.0 ),
    trimPtFracMin = cms.double( -1.0 ),
    zcut = cms.double( -1.0 ),
    rcut_factor = cms.double( -1.0 ),
    csRho_EtaMax = cms.double( -1.0 ),
    csRParam = cms.double( -1.0 ),
    beta = cms.double( -1.0 ),
    R0 = cms.double( -1.0 ),
    gridMaxRapidity = cms.double( -1.0 ),
    gridSpacing = cms.double( -1.0 ),
    DzTrVtxMax = cms.double( 0.0 ),
    DxyTrVtxMax = cms.double( 0.0 ),
    MaxVtxZ = cms.double( 15.0 ),
    subjetPtMin = cms.double( -1.0 ),
    muMin = cms.double( -1.0 ),
    muMax = cms.double( -1.0 ),
    yMin = cms.double( -1.0 ),
    yMax = cms.double( -1.0 ),
    dRMin = cms.double( -1.0 ),
    dRMax = cms.double( -1.0 ),
    maxDepth = cms.int32( -1 ),
    nFilt = cms.int32( -1 ),
    MinVtxNdof = cms.int32( 5 ),
    src = cms.InputTag( "hltTowerMakerForAll" ),
    srcPVs = cms.InputTag( "NotUsed" ),
    jetType = cms.string( "CaloJet" ),
    jetAlgorithm = cms.string( "AntiKt" ),
    rParam = cms.double( 0.4 ),
    inputEtMin = cms.double( 0.3 ),
    inputEMin = cms.double( 0.0 ),
    jetPtMin = cms.double( 1.0 ),
    doPVCorrection = cms.bool( False ),
    doAreaFastjet = cms.bool( False ),
    doRhoFastjet = cms.bool( False ),
    doPUOffsetCorr = cms.bool( False ),
    puPtMin = cms.double( 10.0 ),
    nSigmaPU = cms.double( 1.0 ),
    radiusPU = cms.double( 0.4 ),
    subtractorName = cms.string( "" ),
    useExplicitGhosts = cms.bool( False ),
    doAreaDiskApprox = cms.bool( False ),
    voronoiRfact = cms.double( -9.0 ),
    Rho_EtaMax = cms.double( 4.4 ),
    Ghost_EtaMax = cms.double( 6.0 ),
    Active_Area_Repeats = cms.int32( 5 ),
    GhostArea = cms.double( 0.01 ),
    restrictInputs = cms.bool( False ),
    maxInputs = cms.uint32( 1 ),
    writeCompound = cms.bool( False ),
    writeJetsWithConst = cms.bool( False ),
    doFastJetNonUniform = cms.bool( False ),
    useDeterministicSeed = cms.bool( True ),
    minSeed = cms.uint32( 0 ),
    verbosity = cms.int32( 0 ),
    puWidth = cms.double( 0.0 ),
    nExclude = cms.uint32( 0 ),
    maxBadEcalCells = cms.uint32( 9999999 ),
    maxBadHcalCells = cms.uint32( 9999999 ),
    maxProblematicEcalCells = cms.uint32( 9999999 ),
    maxProblematicHcalCells = cms.uint32( 9999999 ),
    maxRecoveredEcalCells = cms.uint32( 9999999 ),
    maxRecoveredHcalCells = cms.uint32( 9999999 ),
    puCenters = cms.vdouble(  ),
    applyWeight = cms.bool( False ),
    srcWeights = cms.InputTag( "" ),
    minimumTowersFraction = cms.double( 0.0 ),
    jetCollInstanceName = cms.string( "" ),
    sumRecHits = cms.bool( False )
)
process.hltAK4CaloJetsPFEt5 = cms.EDFilter( "EtMinCaloJetSelector",
    src = cms.InputTag( "hltAK4CaloJetsPF" ),
    filter = cms.bool( False ),
    etMin = cms.double( 5.0 )
)
process.hltMuonDTDigis = cms.EDProducer( "DTuROSRawToDigi",
    inputLabel = cms.InputTag( "rawDataCollector" ),
    debug = cms.untracked.bool( False )
)
process.hltDt1DRecHits = cms.EDProducer( "DTRecHitProducer",
    recAlgoConfig = cms.PSet( 
      maxTime = cms.double( 420.0 ),
      debug = cms.untracked.bool( False ),
      stepTwoFromDigi = cms.bool( False ),
      tTrigModeConfig = cms.PSet( 
        debug = cms.untracked.bool( False ),
        tofCorrType = cms.int32( 0 ),
        tTrigLabel = cms.string( "" ),
        wirePropCorrType = cms.int32( 0 ),
        doTOFCorrection = cms.bool( True ),
        vPropWire = cms.double( 24.4 ),
        doT0Correction = cms.bool( True ),
        doWirePropCorrection = cms.bool( True ),
        t0Label = cms.string( "" )
      ),
      useUncertDB = cms.bool( True ),
      doVdriftCorr = cms.bool( True ),
      minTime = cms.double( -3.0 ),
      tTrigMode = cms.string( "DTTTrigSyncFromDB" ),
      readLegacyTTrigDB = cms.bool( True ),
      readLegacyVDriftDB = cms.bool( True )
    ),
    recAlgo = cms.string( "DTLinearDriftFromDBAlgo" ),
    debug = cms.untracked.bool( False ),
    dtDigiLabel = cms.InputTag( "hltMuonDTDigis" )
)
process.hltDt4DSegments = cms.EDProducer( "DTRecSegment4DProducer",
    Reco4DAlgoName = cms.string( "DTCombinatorialPatternReco4D" ),
    Reco4DAlgoConfig = cms.PSet( 
      Reco2DAlgoConfig = cms.PSet( 
        AlphaMaxPhi = cms.double( 1.0 ),
        debug = cms.untracked.bool( False ),
        segmCleanerMode = cms.int32( 2 ),
        AlphaMaxTheta = cms.double( 0.9 ),
        hit_afterT0_resolution = cms.double( 0.03 ),
        performT0_vdriftSegCorrection = cms.bool( False ),
        recAlgo = cms.string( "DTLinearDriftFromDBAlgo" ),
        recAlgoConfig = cms.PSet( 
          maxTime = cms.double( 420.0 ),
          debug = cms.untracked.bool( False ),
          stepTwoFromDigi = cms.bool( False ),
          tTrigModeConfig = cms.PSet( 
            debug = cms.untracked.bool( False ),
            tofCorrType = cms.int32( 0 ),
            tTrigLabel = cms.string( "" ),
            wirePropCorrType = cms.int32( 0 ),
            doTOFCorrection = cms.bool( True ),
            vPropWire = cms.double( 24.4 ),
            doT0Correction = cms.bool( True ),
            doWirePropCorrection = cms.bool( True ),
            t0Label = cms.string( "" )
          ),
          useUncertDB = cms.bool( True ),
          doVdriftCorr = cms.bool( True ),
          minTime = cms.double( -3.0 ),
          tTrigMode = cms.string( "DTTTrigSyncFromDB" ),
          readLegacyTTrigDB = cms.bool( True ),
          readLegacyVDriftDB = cms.bool( True )
        ),
        MaxAllowedHits = cms.uint32( 50 ),
        nUnSharedHitsMin = cms.int32( 2 ),
        nSharedHitsMax = cms.int32( 2 ),
        performT0SegCorrection = cms.bool( False ),
        perform_delta_rejecting = cms.bool( False )
      ),
      Reco2DAlgoName = cms.string( "DTCombinatorialPatternReco" ),
      debug = cms.untracked.bool( False ),
      segmCleanerMode = cms.int32( 2 ),
      AllDTRecHits = cms.bool( True ),
      hit_afterT0_resolution = cms.double( 0.03 ),
      performT0_vdriftSegCorrection = cms.bool( False ),
      recAlgo = cms.string( "DTLinearDriftFromDBAlgo" ),
      recAlgoConfig = cms.PSet( 
        maxTime = cms.double( 420.0 ),
        debug = cms.untracked.bool( False ),
        stepTwoFromDigi = cms.bool( False ),
        tTrigModeConfig = cms.PSet( 
          debug = cms.untracked.bool( False ),
          tofCorrType = cms.int32( 0 ),
          tTrigLabel = cms.string( "" ),
          wirePropCorrType = cms.int32( 0 ),
          doTOFCorrection = cms.bool( True ),
          vPropWire = cms.double( 24.4 ),
          doT0Correction = cms.bool( True ),
          doWirePropCorrection = cms.bool( True ),
          t0Label = cms.string( "" )
        ),
        useUncertDB = cms.bool( True ),
        doVdriftCorr = cms.bool( True ),
        minTime = cms.double( -3.0 ),
        tTrigMode = cms.string( "DTTTrigSyncFromDB" ),
        readLegacyTTrigDB = cms.bool( True ),
        readLegacyVDriftDB = cms.bool( True )
      ),
      nUnSharedHitsMin = cms.int32( 2 ),
      nSharedHitsMax = cms.int32( 2 ),
      performT0SegCorrection = cms.bool( False ),
      perform_delta_rejecting = cms.bool( False )
    ),
    debug = cms.untracked.bool( False ),
    recHits1DLabel = cms.InputTag( "hltDt1DRecHits" ),
    recHits2DLabel = cms.InputTag( "dt2DSegments" )
)
process.hltMuonCSCDigis = cms.EDProducer( "CSCDCCUnpacker",
    InputObjects = cms.InputTag( "rawDataCollector" ),
    UseExaminer = cms.bool( True ),
    ExaminerMask = cms.uint32( 535558134 ),
    UseSelectiveUnpacking = cms.bool( True ),
    ErrorMask = cms.uint32( 0 ),
    UnpackStatusDigis = cms.bool( False ),
    UseFormatStatus = cms.bool( True ),
    useRPCs = cms.bool( False ),
    useGEMs = cms.bool( False ),
    useCSCShowers = cms.bool( False ),
    Debug = cms.untracked.bool( False ),
    PrintEventNumber = cms.untracked.bool( False ),
    runDQM = cms.untracked.bool( False ),
    VisualFEDInspect = cms.untracked.bool( False ),
    VisualFEDShort = cms.untracked.bool( False ),
    FormatedEventDump = cms.untracked.bool( False ),
    SuppressZeroLCT = cms.untracked.bool( True ),
    DisableMappingCheck = cms.untracked.bool( False ),
    B904Setup = cms.untracked.bool( False )
)
process.hltCsc2DRecHits = cms.EDProducer( "CSCRecHitDProducer",
    CSCStripPeakThreshold = cms.double( 10.0 ),
    CSCStripClusterChargeCut = cms.double( 25.0 ),
    CSCStripxtalksOffset = cms.double( 0.03 ),
    UseAverageTime = cms.bool( False ),
    UseParabolaFit = cms.bool( False ),
    UseFivePoleFit = cms.bool( True ),
    CSCWireClusterDeltaT = cms.int32( 1 ),
    CSCUseCalibrations = cms.bool( True ),
    CSCUseStaticPedestals = cms.bool( False ),
    CSCNoOfTimeBinsForDynamicPedestal = cms.int32( 2 ),
    wireDigiTag = cms.InputTag( 'hltMuonCSCDigis','MuonCSCWireDigi' ),
    stripDigiTag = cms.InputTag( 'hltMuonCSCDigis','MuonCSCStripDigi' ),
    readBadChannels = cms.bool( False ),
    readBadChambers = cms.bool( True ),
    CSCUseTimingCorrections = cms.bool( True ),
    CSCUseGasGainCorrections = cms.bool( False ),
    CSCDebug = cms.untracked.bool( False ),
    CSCstripWireDeltaTime = cms.int32( 8 ),
    XTasymmetry_ME1a = cms.double( 0.0 ),
    XTasymmetry_ME1b = cms.double( 0.0 ),
    XTasymmetry_ME12 = cms.double( 0.0 ),
    XTasymmetry_ME13 = cms.double( 0.0 ),
    XTasymmetry_ME21 = cms.double( 0.0 ),
    XTasymmetry_ME22 = cms.double( 0.0 ),
    XTasymmetry_ME31 = cms.double( 0.0 ),
    XTasymmetry_ME32 = cms.double( 0.0 ),
    XTasymmetry_ME41 = cms.double( 0.0 ),
    ConstSyst_ME1a = cms.double( 0.022 ),
    ConstSyst_ME1b = cms.double( 0.007 ),
    ConstSyst_ME12 = cms.double( 0.0 ),
    ConstSyst_ME13 = cms.double( 0.0 ),
    ConstSyst_ME21 = cms.double( 0.0 ),
    ConstSyst_ME22 = cms.double( 0.0 ),
    ConstSyst_ME31 = cms.double( 0.0 ),
    ConstSyst_ME32 = cms.double( 0.0 ),
    ConstSyst_ME41 = cms.double( 0.0 ),
    NoiseLevel_ME1a = cms.double( 7.0 ),
    NoiseLevel_ME1b = cms.double( 8.0 ),
    NoiseLevel_ME12 = cms.double( 9.0 ),
    NoiseLevel_ME13 = cms.double( 8.0 ),
    NoiseLevel_ME21 = cms.double( 9.0 ),
    NoiseLevel_ME22 = cms.double( 9.0 ),
    NoiseLevel_ME31 = cms.double( 9.0 ),
    NoiseLevel_ME32 = cms.double( 9.0 ),
    NoiseLevel_ME41 = cms.double( 9.0 ),
    CSCUseReducedWireTimeWindow = cms.bool( False ),
    CSCWireTimeWindowLow = cms.int32( 0 ),
    CSCWireTimeWindowHigh = cms.int32( 15 )
)
process.hltCscSegments = cms.EDProducer( "CSCSegmentProducer",
    inputObjects = cms.InputTag( "hltCsc2DRecHits" ),
    algo_type = cms.int32( 1 ),
    algo_psets = cms.VPSet( 
      cms.PSet(  parameters_per_chamber_type = cms.vint32( 1, 2, 3, 4, 5, 6, 5, 6, 5, 6 ),
        algo_psets = cms.VPSet( 
          cms.PSet(  wideSeg = cms.double( 3.0 ),
            chi2Norm_2D_ = cms.double( 35.0 ),
            dRIntMax = cms.double( 2.0 ),
            doCollisions = cms.bool( True ),
            dPhiMax = cms.double( 0.006 ),
            dRMax = cms.double( 1.5 ),
            dPhiIntMax = cms.double( 0.005 ),
            minLayersApart = cms.int32( 1 ),
            chi2Max = cms.double( 100.0 ),
            chi2_str = cms.double( 50.0 ),
            enlarge = cms.bool( False )
          ),
          cms.PSet(  wideSeg = cms.double( 3.0 ),
            chi2Norm_2D_ = cms.double( 35.0 ),
            dRIntMax = cms.double( 2.0 ),
            doCollisions = cms.bool( True ),
            dPhiMax = cms.double( 0.005 ),
            dRMax = cms.double( 1.5 ),
            dPhiIntMax = cms.double( 0.004 ),
            minLayersApart = cms.int32( 1 ),
            chi2Max = cms.double( 100.0 ),
            chi2_str = cms.double( 50.0 ),
            enlarge = cms.bool( False )
          ),
          cms.PSet(  wideSeg = cms.double( 3.0 ),
            chi2Norm_2D_ = cms.double( 35.0 ),
            dRIntMax = cms.double( 2.0 ),
            doCollisions = cms.bool( True ),
            dPhiMax = cms.double( 0.004 ),
            dRMax = cms.double( 1.5 ),
            dPhiIntMax = cms.double( 0.003 ),
            minLayersApart = cms.int32( 1 ),
            chi2Max = cms.double( 100.0 ),
            chi2_str = cms.double( 50.0 ),
            enlarge = cms.bool( False )
          ),
          cms.PSet(  wideSeg = cms.double( 3.0 ),
            chi2Norm_2D_ = cms.double( 20.0 ),
            dRIntMax = cms.double( 2.0 ),
            doCollisions = cms.bool( True ),
            dPhiMax = cms.double( 0.003 ),
            dRMax = cms.double( 1.5 ),
            dPhiIntMax = cms.double( 0.002 ),
            minLayersApart = cms.int32( 1 ),
            chi2Max = cms.double( 60.0 ),
            chi2_str = cms.double( 30.0 ),
            enlarge = cms.bool( False )
          ),
          cms.PSet(  wideSeg = cms.double( 3.0 ),
            chi2Norm_2D_ = cms.double( 60.0 ),
            dRIntMax = cms.double( 2.0 ),
            doCollisions = cms.bool( True ),
            dPhiMax = cms.double( 0.007 ),
            dRMax = cms.double( 1.5 ),
            dPhiIntMax = cms.double( 0.005 ),
            minLayersApart = cms.int32( 1 ),
            chi2Max = cms.double( 180.0 ),
            chi2_str = cms.double( 80.0 ),
            enlarge = cms.bool( False )
          ),
          cms.PSet(  wideSeg = cms.double( 3.0 ),
            chi2Norm_2D_ = cms.double( 35.0 ),
            dRIntMax = cms.double( 2.0 ),
            doCollisions = cms.bool( True ),
            dPhiMax = cms.double( 0.006 ),
            dRMax = cms.double( 1.5 ),
            dPhiIntMax = cms.double( 0.004 ),
            minLayersApart = cms.int32( 1 ),
            chi2Max = cms.double( 100.0 ),
            chi2_str = cms.double( 50.0 ),
            enlarge = cms.bool( False )
          )
        ),
        algo_name = cms.string( "CSCSegAlgoRU" ),
        chamber_types = cms.vstring( 'ME1/a',
          'ME1/b',
          'ME1/2',
          'ME1/3',
          'ME2/1',
          'ME2/2',
          'ME3/1',
          'ME3/2',
          'ME4/1',
          'ME4/2' )
      )
    )
)
process.hltMuonRPCDigis = cms.EDProducer( "RPCUnpackingModule",
    InputLabel = cms.InputTag( "rawDataCollector" ),
    doSynchro = cms.bool( False )
)
process.hltRpcRecHits = cms.EDProducer( "RPCRecHitProducer",
    recAlgoConfig = cms.PSet(  ),
    recAlgo = cms.string( "RPCRecHitStandardAlgo" ),
    rpcDigiLabel = cms.InputTag( "hltMuonRPCDigis" ),
    maskSource = cms.string( "File" ),
    maskvecfile = cms.FileInPath( "RecoLocalMuon/RPCRecHit/data/RPCMaskVec.dat" ),
    deadSource = cms.string( "File" ),
    deadvecfile = cms.FileInPath( "RecoLocalMuon/RPCRecHit/data/RPCDeadVec.dat" )
)
process.hltMuonGEMDigis = cms.EDProducer( "GEMRawToDigiModule",
    InputLabel = cms.InputTag( "rawDataCollector" ),
    useDBEMap = cms.bool( True ),
    keepDAQStatus = cms.bool( False ),
    readMultiBX = cms.bool( False ),
    ge21Off = cms.bool( True ),
    fedIdStart = cms.uint32( 1467 ),
    fedIdEnd = cms.uint32( 1478 )
)
process.hltGemRecHits = cms.EDProducer( "GEMRecHitProducer",
    recAlgoConfig = cms.PSet(  ),
    recAlgo = cms.string( "GEMRecHitStandardAlgo" ),
    gemDigiLabel = cms.InputTag( "hltMuonGEMDigis" ),
    applyMasking = cms.bool( False )
)
process.hltGemSegments = cms.EDProducer( "GEMSegmentProducer",
    gemRecHitLabel = cms.InputTag( "hltGemRecHits" ),
    ge0_name = cms.string( "GE0SegAlgoRU" ),
    algo_name = cms.string( "GEMSegmentAlgorithm" ),
    ge0_pset = cms.PSet( 
      maxChi2GoodSeg = cms.double( 50.0 ),
      maxChi2Prune = cms.double( 50.0 ),
      maxNumberOfHitsPerLayer = cms.uint32( 100 ),
      maxETASeeds = cms.double( 0.1 ),
      maxPhiAdditional = cms.double( 0.001096605744 ),
      minNumberOfHits = cms.uint32( 4 ),
      doCollisions = cms.bool( True ),
      maxPhiSeeds = cms.double( 0.001096605744 ),
      requireCentralBX = cms.bool( True ),
      maxChi2Additional = cms.double( 100.0 ),
      allowWideSegments = cms.bool( True ),
      maxNumberOfHits = cms.uint32( 300 ),
      maxTOFDiff = cms.double( 25.0 )
    ),
    algo_pset = cms.PSet( 
      dYclusBoxMax = cms.double( 5.0 ),
      dXclusBoxMax = cms.double( 1.0 ),
      maxRecHitsInCluster = cms.int32( 4 ),
      preClustering = cms.bool( True ),
      preClusteringUseChaining = cms.bool( True ),
      dEtaChainBoxMax = cms.double( 0.05 ),
      clusterOnlySameBXRecHits = cms.bool( True ),
      minHitsPerSegment = cms.uint32( 2 ),
      dPhiChainBoxMax = cms.double( 0.02 )
    )
)
process.hltL2OfflineMuonSeeds = cms.EDProducer( "MuonSeedGenerator",
    beamSpotTag = cms.InputTag( "hltOnlineBeamSpot" ),
    scaleDT = cms.bool( True ),
    CSCRecSegmentLabel = cms.InputTag( "hltCscSegments" ),
    DTRecSegmentLabel = cms.InputTag( "hltDt4DSegments" ),
    ME0RecSegmentLabel = cms.InputTag( "me0Segments" ),
    EnableDTMeasurement = cms.bool( True ),
    EnableCSCMeasurement = cms.bool( True ),
    EnableME0Measurement = cms.bool( False ),
    crackEtas = cms.vdouble( 0.2, 1.6, 1.7 ),
    crackWindow = cms.double( 0.04 ),
    deltaPhiSearchWindow = cms.double( 0.25 ),
    deltaEtaSearchWindow = cms.double( 0.2 ),
    deltaEtaCrackSearchWindow = cms.double( 0.25 ),
    CSC_01 = cms.vdouble( 0.166, 0.0, 0.0, 0.031, 0.0, 0.0 ),
    CSC_12 = cms.vdouble( -0.161, 0.254, -0.047, 0.042, -0.007, 0.0 ),
    CSC_02 = cms.vdouble( 0.612, -0.207, 0.0, 0.067, -0.001, 0.0 ),
    CSC_13 = cms.vdouble( 0.901, -1.302, 0.533, 0.045, 0.005, 0.0 ),
    CSC_03 = cms.vdouble( 0.787, -0.338, 0.029, 0.101, -0.008, 0.0 ),
    CSC_14 = cms.vdouble( 0.606, -0.181, -0.002, 0.111, -0.003, 0.0 ),
    CSC_23 = cms.vdouble( -0.081, 0.113, -0.029, 0.015, 0.008, 0.0 ),
    CSC_24 = cms.vdouble( 0.004, 0.021, -0.002, 0.053, 0.0, 0.0 ),
    CSC_34 = cms.vdouble( 0.062, -0.067, 0.019, 0.021, 0.003, 0.0 ),
    DT_12 = cms.vdouble( 0.183, 0.054, -0.087, 0.028, 0.002, 0.0 ),
    DT_13 = cms.vdouble( 0.315, 0.068, -0.127, 0.051, -0.002, 0.0 ),
    DT_14 = cms.vdouble( 0.359, 0.052, -0.107, 0.072, -0.004, 0.0 ),
    DT_23 = cms.vdouble( 0.13, 0.023, -0.057, 0.028, 0.004, 0.0 ),
    DT_24 = cms.vdouble( 0.176, 0.014, -0.051, 0.051, 0.003, 0.0 ),
    DT_34 = cms.vdouble( 0.044, 0.004, -0.013, 0.029, 0.003, 0.0 ),
    OL_1213 = cms.vdouble( 0.96, -0.737, 0.0, 0.052, 0.0, 0.0 ),
    OL_1222 = cms.vdouble( 0.848, -0.591, 0.0, 0.062, 0.0, 0.0 ),
    OL_1232 = cms.vdouble( 0.184, 0.0, 0.0, 0.066, 0.0, 0.0 ),
    OL_2213 = cms.vdouble( 0.117, 0.0, 0.0, 0.044, 0.0, 0.0 ),
    OL_2222 = cms.vdouble( 0.107, 0.0, 0.0, 0.04, 0.0, 0.0 ),
    SME_11 = cms.vdouble( 3.295, -1.527, 0.112, 0.378, 0.02, 0.0 ),
    SME_12 = cms.vdouble( 0.102, 0.599, 0.0, 0.38, 0.0, 0.0 ),
    SME_13 = cms.vdouble( -1.286, 1.711, 0.0, 0.356, 0.0, 0.0 ),
    SME_21 = cms.vdouble( -0.529, 1.194, -0.358, 0.472, 0.086, 0.0 ),
    SME_22 = cms.vdouble( -1.207, 1.491, -0.251, 0.189, 0.243, 0.0 ),
    SME_31 = cms.vdouble( -1.594, 1.482, -0.317, 0.487, 0.097, 0.0 ),
    SME_32 = cms.vdouble( -0.901, 1.333, -0.47, 0.41, 0.073, 0.0 ),
    SME_41 = cms.vdouble( -0.003, 0.005, 0.005, 0.608, 0.076, 0.0 ),
    SME_42 = cms.vdouble( -0.003, 0.005, 0.005, 0.608, 0.076, 0.0 ),
    SMB_10 = cms.vdouble( 1.387, -0.038, 0.0, 0.19, 0.0, 0.0 ),
    SMB_11 = cms.vdouble( 1.247, 0.72, -0.802, 0.229, -0.075, 0.0 ),
    SMB_12 = cms.vdouble( 2.128, -0.956, 0.0, 0.199, 0.0, 0.0 ),
    SMB_20 = cms.vdouble( 1.011, -0.052, 0.0, 0.188, 0.0, 0.0 ),
    SMB_21 = cms.vdouble( 1.043, -0.124, 0.0, 0.183, 0.0, 0.0 ),
    SMB_22 = cms.vdouble( 1.474, -0.758, 0.0, 0.185, 0.0, 0.0 ),
    SMB_30 = cms.vdouble( 0.505, -0.022, 0.0, 0.215, 0.0, 0.0 ),
    SMB_31 = cms.vdouble( 0.549, -0.145, 0.0, 0.207, 0.0, 0.0 ),
    SMB_32 = cms.vdouble( 0.67, -0.327, 0.0, 0.22, 0.0, 0.0 ),
    CSC_01_1_scale = cms.vdouble( -1.915329, 0.0 ),
    CSC_12_1_scale = cms.vdouble( -6.434242, 0.0 ),
    CSC_12_2_scale = cms.vdouble( -1.63622, 0.0 ),
    CSC_12_3_scale = cms.vdouble( -1.63622, 0.0 ),
    CSC_13_2_scale = cms.vdouble( -6.077936, 0.0 ),
    CSC_13_3_scale = cms.vdouble( -1.701268, 0.0 ),
    CSC_14_3_scale = cms.vdouble( -1.969563, 0.0 ),
    CSC_23_1_scale = cms.vdouble( -19.084285, 0.0 ),
    CSC_23_2_scale = cms.vdouble( -6.079917, 0.0 ),
    CSC_24_1_scale = cms.vdouble( -6.055701, 0.0 ),
    CSC_34_1_scale = cms.vdouble( -11.520507, 0.0 ),
    OL_1213_0_scale = cms.vdouble( -4.488158, 0.0 ),
    OL_1222_0_scale = cms.vdouble( -5.810449, 0.0 ),
    OL_1232_0_scale = cms.vdouble( -5.964634, 0.0 ),
    OL_2213_0_scale = cms.vdouble( -7.239789, 0.0 ),
    OL_2222_0_scale = cms.vdouble( -7.667231, 0.0 ),
    DT_12_1_scale = cms.vdouble( -3.692398, 0.0 ),
    DT_12_2_scale = cms.vdouble( -3.518165, 0.0 ),
    DT_13_1_scale = cms.vdouble( -4.520923, 0.0 ),
    DT_13_2_scale = cms.vdouble( -4.257687, 0.0 ),
    DT_14_1_scale = cms.vdouble( -5.644816, 0.0 ),
    DT_14_2_scale = cms.vdouble( -4.808546, 0.0 ),
    DT_23_1_scale = cms.vdouble( -5.320346, 0.0 ),
    DT_23_2_scale = cms.vdouble( -5.117625, 0.0 ),
    DT_24_1_scale = cms.vdouble( -7.490909, 0.0 ),
    DT_24_2_scale = cms.vdouble( -6.63094, 0.0 ),
    DT_34_1_scale = cms.vdouble( -13.783765, 0.0 ),
    DT_34_2_scale = cms.vdouble( -11.901897, 0.0 ),
    SMB_10_0_scale = cms.vdouble( 2.448566, 0.0 ),
    SMB_11_0_scale = cms.vdouble( 2.56363, 0.0 ),
    SMB_12_0_scale = cms.vdouble( 2.283221, 0.0 ),
    SMB_20_0_scale = cms.vdouble( 1.486168, 0.0 ),
    SMB_21_0_scale = cms.vdouble( 1.58384, 0.0 ),
    SMB_22_0_scale = cms.vdouble( 1.346681, 0.0 ),
    SMB_30_0_scale = cms.vdouble( -3.629838, 0.0 ),
    SMB_31_0_scale = cms.vdouble( -3.323768, 0.0 ),
    SMB_32_0_scale = cms.vdouble( -3.054156, 0.0 ),
    SME_11_0_scale = cms.vdouble( 1.325085, 0.0 ),
    SME_12_0_scale = cms.vdouble( 2.279181, 0.0 ),
    SME_13_0_scale = cms.vdouble( 0.104905, 0.0 ),
    SME_21_0_scale = cms.vdouble( -0.040862, 0.0 ),
    SME_22_0_scale = cms.vdouble( -3.457901, 0.0 )
)
process.hltL2MuonSeeds = cms.EDProducer( "L2MuonSeedGeneratorFromL1T",
    GMTReadoutCollection = cms.InputTag( "" ),
    InputObjects = cms.InputTag( 'hltGtStage2Digis','Muon' ),
    Propagator = cms.string( "SteppingHelixPropagatorAny" ),
    L1MinPt = cms.double( 0.0 ),
    L1MaxEta = cms.double( 2.5 ),
    L1MinQuality = cms.uint32( 7 ),
    SetMinPtBarrelTo = cms.double( 3.5 ),
    SetMinPtEndcapTo = cms.double( 1.0 ),
    UseOfflineSeed = cms.untracked.bool( True ),
    UseUnassociatedL1 = cms.bool( False ),
    MatchDR = cms.vdouble( 0.3 ),
    EtaMatchingBins = cms.vdouble( 0.0, 2.5 ),
    CentralBxOnly = cms.bool( True ),
    MatchType = cms.uint32( 0 ),
    SortType = cms.uint32( 0 ),
    OfflineSeedLabel = cms.untracked.InputTag( "hltL2OfflineMuonSeeds" ),
    ServiceParameters = cms.PSet( 
      RPCLayers = cms.bool( True ),
      UseMuonNavigation = cms.untracked.bool( True ),
      Propagators = cms.untracked.vstring( 'SteppingHelixPropagatorAny' )
    )
)
process.hltL2Muons = cms.EDProducer( "L2MuonProducer",
    ServiceParameters = cms.PSet( 
      RPCLayers = cms.bool( True ),
      UseMuonNavigation = cms.untracked.bool( True ),
      Propagators = cms.untracked.vstring( 'hltESPFastSteppingHelixPropagatorAny',
        'hltESPFastSteppingHelixPropagatorOpposite' )
    ),
    InputObjects = cms.InputTag( "hltL2MuonSeeds" ),
    SeedTransformerParameters = cms.PSet( 
      Fitter = cms.string( "hltESPKFFittingSmootherForL2Muon" ),
      NMinRecHits = cms.uint32( 2 ),
      RescaleError = cms.double( 100.0 ),
      Propagator = cms.string( "hltESPFastSteppingHelixPropagatorAny" ),
      UseSubRecHits = cms.bool( False ),
      MuonRecHitBuilder = cms.string( "hltESPMuonTransientTrackingRecHitBuilder" )
    ),
    L2TrajBuilderParameters = cms.PSet( 
      BWFilterParameters = cms.PSet( 
        DTRecSegmentLabel = cms.InputTag( "hltDt4DSegments" ),
        CSCRecSegmentLabel = cms.InputTag( "hltCscSegments" ),
        BWSeedType = cms.string( "fromGenerator" ),
        GEMRecSegmentLabel = cms.InputTag( "hltGemRecHits" ),
        RPCRecSegmentLabel = cms.InputTag( "hltRpcRecHits" ),
        EnableGEMMeasurement = cms.bool( True ),
        EnableRPCMeasurement = cms.bool( True ),
        MuonTrajectoryUpdatorParameters = cms.PSet( 
          ExcludeRPCFromFit = cms.bool( False ),
          Granularity = cms.int32( 0 ),
          MaxChi2 = cms.double( 25.0 ),
          RescaleError = cms.bool( False ),
          RescaleErrorFactor = cms.double( 100.0 ),
          UseInvalidHits = cms.bool( True )
        ),
        EnableCSCMeasurement = cms.bool( True ),
        MaxChi2 = cms.double( 100.0 ),
        FitDirection = cms.string( "outsideIn" ),
        Propagator = cms.string( "hltESPFastSteppingHelixPropagatorAny" ),
        NumberOfSigma = cms.double( 3.0 ),
        EnableDTMeasurement = cms.bool( True )
      ),
      DoSeedRefit = cms.bool( False ),
      FilterParameters = cms.PSet( 
        DTRecSegmentLabel = cms.InputTag( "hltDt4DSegments" ),
        CSCRecSegmentLabel = cms.InputTag( "hltCscSegments" ),
        GEMRecSegmentLabel = cms.InputTag( "hltGemRecHits" ),
        RPCRecSegmentLabel = cms.InputTag( "hltRpcRecHits" ),
        EnableGEMMeasurement = cms.bool( True ),
        EnableRPCMeasurement = cms.bool( True ),
        MuonTrajectoryUpdatorParameters = cms.PSet( 
          ExcludeRPCFromFit = cms.bool( False ),
          Granularity = cms.int32( 0 ),
          MaxChi2 = cms.double( 25.0 ),
          RescaleError = cms.bool( False ),
          RescaleErrorFactor = cms.double( 100.0 ),
          UseInvalidHits = cms.bool( True )
        ),
        EnableCSCMeasurement = cms.bool( True ),
        MaxChi2 = cms.double( 1000.0 ),
        FitDirection = cms.string( "insideOut" ),
        Propagator = cms.string( "hltESPFastSteppingHelixPropagatorAny" ),
        NumberOfSigma = cms.double( 3.0 ),
        EnableDTMeasurement = cms.bool( True )
      ),
      SeedPosition = cms.string( "in" ),
      DoBackwardFilter = cms.bool( True ),
      DoRefit = cms.bool( False ),
      NavigationType = cms.string( "Standard" ),
      SeedTransformerParameters = cms.PSet( 
        Fitter = cms.string( "hltESPKFFittingSmootherForL2Muon" ),
        NMinRecHits = cms.uint32( 2 ),
        RescaleError = cms.double( 100.0 ),
        Propagator = cms.string( "hltESPFastSteppingHelixPropagatorAny" ),
        UseSubRecHits = cms.bool( False ),
        MuonRecHitBuilder = cms.string( "hltESPMuonTransientTrackingRecHitBuilder" )
      ),
      SeedPropagator = cms.string( "hltESPFastSteppingHelixPropagatorAny" )
    ),
    DoSeedRefit = cms.bool( False ),
    TrackLoaderParameters = cms.PSet( 
      TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
      beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
      DoSmoothing = cms.bool( False ),
      VertexConstraint = cms.bool( True ),
      MuonUpdatorAtVertexParameters = cms.PSet( 
        MaxChi2 = cms.double( 1000000.0 ),
        BeamSpotPositionErrors = cms.vdouble( 0.1, 0.1, 5.3 ),
        BeamSpotPosition = cms.vdouble( 0.0, 0.0, 0.0 ),
        Propagator = cms.string( "hltESPFastSteppingHelixPropagatorOpposite" )
      ),
      Smoother = cms.string( "hltESPKFTrajectorySmootherForMuonTrackLoader" )
    ),
    MuonTrajectoryBuilder = cms.string( "Exhaustive" )
)
process.hltL2MuonCandidates = cms.EDProducer( "L2MuonCandidateProducer",
    InputObjects = cms.InputTag( 'hltL2Muons','UpdatedAtVtx' )
)
process.hltSiPixelDigisLegacy = cms.EDProducer( "SiPixelRawToDigi",
    IncludeErrors = cms.bool( True ),
    UseQualityInfo = cms.bool( False ),
    ErrorList = cms.vint32( 29 ),
    UserErrorList = cms.vint32(  ),
    InputLabel = cms.InputTag( "rawDataCollector" ),
    Regions = cms.PSet(  ),
    UsePilotBlade = cms.bool( False ),
    UsePhase1 = cms.bool( True ),
    CablingMapLabel = cms.string( "" ),
    SiPixelQualityLabel = cms.string( "" )
)
process.hltSiPixelClustersLegacy = cms.EDProducer( "SiPixelClusterProducer",
    src = cms.InputTag( "hltSiPixelDigisLegacy" ),
    ClusterMode = cms.string( "PixelThresholdClusterizer" ),
    maxNumberOfClusters = cms.int32( 40000 ),
    payloadType = cms.string( "HLT" ),
    ChannelThreshold = cms.int32( 10 ),
    MissCalibrate = cms.bool( True ),
    SplitClusters = cms.bool( False ),
    VCaltoElectronGain = cms.int32( 1 ),
    VCaltoElectronGain_L1 = cms.int32( 1 ),
    VCaltoElectronOffset = cms.int32( 0 ),
    VCaltoElectronOffset_L1 = cms.int32( 0 ),
    SeedThreshold = cms.int32( 1000 ),
    ClusterThreshold_L1 = cms.int32( 4000 ),
    ClusterThreshold = cms.int32( 4000 ),
    ElectronPerADCGain = cms.double( 135.0 ),
    Phase2Calibration = cms.bool( False ),
    Phase2ReadoutMode = cms.int32( -1 ),
    Phase2DigiBaseline = cms.double( 1200.0 ),
    Phase2KinkADC = cms.int32( 8 )
)
process.hltSiPixelClustersCache = cms.EDProducer( "SiPixelClusterShapeCacheProducer",
    src = cms.InputTag( "hltSiPixelClusters" ),
    onDemand = cms.bool( False )
)
process.hltSiPixelRecHitSoA = cms.EDProducer( "SiPixelRecHitSoAFromLegacy",
    beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
    src = cms.InputTag( "hltSiPixelClusters" ),
    CPE = cms.string( "PixelCPEFast" ),
    convertToLegacy = cms.bool( True ),
    isPhase2 = cms.bool( False )
)
process.hltSiStripExcludedFEDListProducer = cms.EDProducer( "SiStripExcludedFEDListProducer",
    ProductLabel = cms.InputTag( "rawDataCollector" )
)
process.hltSiStripRawToClustersFacility = cms.EDProducer( "SiStripClusterizerFromRaw",
    onDemand = cms.bool( True ),
    Clusterizer = cms.PSet( 
      ConditionsLabel = cms.string( "" ),
      ClusterThreshold = cms.double( 5.0 ),
      SeedThreshold = cms.double( 3.0 ),
      Algorithm = cms.string( "ThreeThresholdAlgorithm" ),
      ChannelThreshold = cms.double( 2.0 ),
      MaxAdjacentBad = cms.uint32( 0 ),
      setDetId = cms.bool( True ),
      MaxSequentialHoles = cms.uint32( 0 ),
      RemoveApvShots = cms.bool( True ),
      clusterChargeCut = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
      MaxSequentialBad = cms.uint32( 1 )
    ),
    Algorithms = cms.PSet( 
      Use10bitsTruncation = cms.bool( False ),
      CommonModeNoiseSubtractionMode = cms.string( "Median" ),
      useCMMeanMap = cms.bool( False ),
      TruncateInSuppressor = cms.bool( True ),
      doAPVRestore = cms.bool( False ),
      SiStripFedZeroSuppressionMode = cms.uint32( 4 ),
      PedestalSubtractionFedMode = cms.bool( True )
    ),
    DoAPVEmulatorCheck = cms.bool( False ),
    HybridZeroSuppressed = cms.bool( False ),
    ProductLabel = cms.InputTag( "rawDataCollector" )
)
process.hltSiStripClusters = cms.EDProducer( "MeasurementTrackerEventProducer",
    measurementTracker = cms.string( "hltESPMeasurementTracker" ),
    skipClusters = cms.InputTag( "" ),
    pixelClusterProducer = cms.string( "hltSiPixelClusters" ),
    stripClusterProducer = cms.string( "hltSiStripRawToClustersFacility" ),
    Phase2TrackerCluster1DProducer = cms.string( "" ),
    vectorHits = cms.InputTag( "" ),
    vectorHitsRej = cms.InputTag( "" ),
    inactivePixelDetectorLabels = cms.VInputTag( 'hltSiPixelDigis' ),
    badPixelFEDChannelCollectionLabels = cms.VInputTag( 'hltSiPixelDigis' ),
    pixelCablingMapLabel = cms.string( "" ),
    inactiveStripDetectorLabels = cms.VInputTag( 'hltSiStripExcludedFEDListProducer' ),
    switchOffPixelsIfEmpty = cms.bool( True )
)
process.hltIterL3OISeedsFromL2Muons = cms.EDProducer( "TSGForOIDNN",
    src = cms.InputTag( 'hltL2Muons','UpdatedAtVtx' ),
    layersToTry = cms.int32( 2 ),
    fixedErrorRescaleFactorForHitless = cms.double( 2.0 ),
    hitsToTry = cms.int32( 1 ),
    MeasurementTrackerEvent = cms.InputTag( "hltSiStripClusters" ),
    estimator = cms.string( "hltESPChi2MeasurementEstimator100" ),
    maxEtaForTOB = cms.double( 1.8 ),
    minEtaForTEC = cms.double( 0.7 ),
    debug = cms.untracked.bool( False ),
    maxSeeds = cms.uint32( 20 ),
    maxHitlessSeeds = cms.uint32( 5 ),
    maxHitSeeds = cms.uint32( 1 ),
    propagatorName = cms.string( "PropagatorWithMaterialParabolicMf" ),
    maxHitlessSeedsIP = cms.uint32( 5 ),
    maxHitlessSeedsMuS = cms.uint32( 0 ),
    maxHitDoubletSeeds = cms.uint32( 0 ),
    getStrategyFromDNN = cms.bool( True ),
    useRegressor = cms.bool( False ),
    dnnMetadataPath = cms.string( "RecoMuon/TrackerSeedGenerator/data/OIseeding/DNNclassifier_Run3_metadata.json" )
)
process.hltIterL3OITrackCandidates = cms.EDProducer( "CkfTrackCandidateMaker",
    cleanTrajectoryAfterInOut = cms.bool( False ),
    doSeedingRegionRebuilding = cms.bool( False ),
    onlyPixelHitsForSeedCleaner = cms.bool( False ),
    reverseTrajectories = cms.bool( True ),
    useHitsSplitting = cms.bool( False ),
    MeasurementTrackerEvent = cms.InputTag( "hltSiStripClusters" ),
    src = cms.InputTag( "hltIterL3OISeedsFromL2Muons" ),
    clustersToSkip = cms.InputTag( "" ),
    phase2clustersToSkip = cms.InputTag( "" ),
    TrajectoryBuilderPSet = cms.PSet(  refToPSet_ = cms.string( "HLTPSetMuonCkfTrajectoryBuilder" ) ),
    TransientInitialStateEstimatorParameters = cms.PSet( 
      propagatorAlongTISE = cms.string( "PropagatorWithMaterial" ),
      numberMeasurementsForFit = cms.int32( 4 ),
      propagatorOppositeTISE = cms.string( "PropagatorWithMaterialOpposite" )
    ),
    numHitsForSeedCleaner = cms.int32( 4 ),
    NavigationSchool = cms.string( "SimpleNavigationSchool" ),
    RedundantSeedCleaner = cms.string( "CachingSeedCleanerBySharedInput" ),
    TrajectoryCleaner = cms.string( "muonSeededTrajectoryCleanerBySharedHits" ),
    maxNSeeds = cms.uint32( 500000 ),
    maxSeedsBeforeCleaning = cms.uint32( 5000 )
)
process.hltIterL3OIMuCtfWithMaterialTracks = cms.EDProducer( "TrackProducer",
    useSimpleMF = cms.bool( False ),
    SimpleMagneticField = cms.string( "" ),
    src = cms.InputTag( "hltIterL3OITrackCandidates" ),
    clusterRemovalInfo = cms.InputTag( "" ),
    beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
    Fitter = cms.string( "hltESPKFFittingSmootherWithOutliersRejectionAndRK" ),
    useHitsSplitting = cms.bool( False ),
    alias = cms.untracked.string( "ctfWithMaterialTracks" ),
    TrajectoryInEvent = cms.bool( False ),
    TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
    AlgorithmName = cms.string( "iter10" ),
    Propagator = cms.string( "PropagatorWithMaterial" ),
    GeometricInnerState = cms.bool( True ),
    NavigationSchool = cms.string( "SimpleNavigationSchool" ),
    MeasurementTracker = cms.string( "hltESPMeasurementTracker" ),
    MeasurementTrackerEvent = cms.InputTag( "hltSiStripClusters" )
)
process.hltIterL3OIMuonTrackCutClassifier = cms.EDProducer( "TrackCutClassifier",
    src = cms.InputTag( "hltIterL3OIMuCtfWithMaterialTracks" ),
    beamspot = cms.InputTag( "hltOnlineBeamSpot" ),
    vertices = cms.InputTag( "Notused" ),
    ignoreVertices = cms.bool( True ),
    qualityCuts = cms.vdouble( -0.7, 0.1, 0.7 ),
    mva = cms.PSet( 
      minPixelHits = cms.vint32( 0, 0, 1 ),
      maxDzWrtBS = cms.vdouble( 3.40282346639E38, 24.0, 100.0 ),
      dr_par = cms.PSet( 
        d0err = cms.vdouble( 0.003, 0.003, 3.40282346639E38 ),
        dr_par2 = cms.vdouble( 0.3, 0.3, 3.40282346639E38 ),
        dr_par1 = cms.vdouble( 0.4, 0.4, 3.40282346639E38 ),
        dr_exp = cms.vint32( 4, 4, 2147483647 ),
        d0err_par = cms.vdouble( 0.001, 0.001, 3.40282346639E38 )
      ),
      maxLostLayers = cms.vint32( 4, 3, 2 ),
      min3DLayers = cms.vint32( 1, 2, 1 ),
      dz_par = cms.PSet( 
        dz_par1 = cms.vdouble( 0.4, 0.4, 3.40282346639E38 ),
        dz_par2 = cms.vdouble( 0.35, 0.35, 3.40282346639E38 ),
        dz_exp = cms.vint32( 4, 4, 2147483647 )
      ),
      minNVtxTrk = cms.int32( 3 ),
      maxDz = cms.vdouble( 0.5, 0.2, 3.40282346639E38 ),
      minNdof = cms.vdouble( 1.0E-5, 1.0E-5, 1.0E-5 ),
      maxChi2 = cms.vdouble( 3.40282346639E38, 3.40282346639E38, 3.40282346639E38 ),
      maxChi2n = cms.vdouble( 10.0, 1.0, 0.4 ),
      maxDr = cms.vdouble( 0.5, 0.03, 3.40282346639E38 ),
      minLayers = cms.vint32( 3, 5, 5 )
    )
)
process.hltIterL3OIMuonTrackSelectionHighPurity = cms.EDProducer( "TrackCollectionFilterCloner",
    originalSource = cms.InputTag( "hltIterL3OIMuCtfWithMaterialTracks" ),
    originalMVAVals = cms.InputTag( 'hltIterL3OIMuonTrackCutClassifier','MVAValues' ),
    originalQualVals = cms.InputTag( 'hltIterL3OIMuonTrackCutClassifier','QualityMasks' ),
    minQuality = cms.string( "highPurity" ),
    copyExtras = cms.untracked.bool( True ),
    copyTrajectories = cms.untracked.bool( False )
)
process.hltL3MuonsIterL3OI = cms.EDProducer( "L3MuonProducer",
    ServiceParameters = cms.PSet( 
      RPCLayers = cms.bool( True ),
      UseMuonNavigation = cms.untracked.bool( True ),
      Propagators = cms.untracked.vstring( 'hltESPSmartPropagatorAny',
        'SteppingHelixPropagatorAny',
        'hltESPSmartPropagator',
        'hltESPSteppingHelixPropagatorOpposite' )
    ),
    MuonCollectionLabel = cms.InputTag( 'hltL2Muons','UpdatedAtVtx' ),
    TrackLoaderParameters = cms.PSet( 
      MuonSeededTracksInstance = cms.untracked.string( "L2Seeded" ),
      TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
      beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
      DoSmoothing = cms.bool( True ),
      SmoothTkTrack = cms.untracked.bool( False ),
      VertexConstraint = cms.bool( False ),
      MuonUpdatorAtVertexParameters = cms.PSet( 
        MaxChi2 = cms.double( 1000000.0 ),
        BeamSpotPositionErrors = cms.vdouble( 0.1, 0.1, 5.3 ),
        Propagator = cms.string( "hltESPSteppingHelixPropagatorOpposite" )
      ),
      PutTkTrackIntoEvent = cms.untracked.bool( False ),
      Smoother = cms.string( "hltESPKFTrajectorySmootherForMuonTrackLoader" )
    ),
    L3TrajBuilderParameters = cms.PSet( 
      PtCut = cms.double( 1.0 ),
      TrackerPropagator = cms.string( "SteppingHelixPropagatorAny" ),
      GlobalMuonTrackMatcher = cms.PSet( 
        Chi2Cut_3 = cms.double( 200.0 ),
        DeltaDCut_2 = cms.double( 10.0 ),
        Eta_threshold = cms.double( 1.2 ),
        Quality_2 = cms.double( 15.0 ),
        DeltaDCut_1 = cms.double( 40.0 ),
        Quality_3 = cms.double( 7.0 ),
        DeltaDCut_3 = cms.double( 15.0 ),
        Quality_1 = cms.double( 20.0 ),
        Pt_threshold1 = cms.double( 0.0 ),
        DeltaRCut_2 = cms.double( 0.2 ),
        DeltaRCut_1 = cms.double( 0.1 ),
        Pt_threshold2 = cms.double( 9.99999999E8 ),
        Chi2Cut_1 = cms.double( 50.0 ),
        Chi2Cut_2 = cms.double( 50.0 ),
        DeltaRCut_3 = cms.double( 1.0 ),
        LocChi2Cut = cms.double( 0.001 ),
        Propagator = cms.string( "hltESPSmartPropagator" ),
        MinPt = cms.double( 1.0 ),
        MinP = cms.double( 2.5 )
      ),
      ScaleTECxFactor = cms.double( -1.0 ),
      tkTrajUseVertex = cms.bool( False ),
      MuonTrackingRegionBuilder = cms.PSet( 
        Rescale_Dz = cms.double( 4.0 ),
        Pt_fixed = cms.bool( False ),
        Eta_fixed = cms.bool( True ),
        Eta_min = cms.double( 0.1 ),
        DeltaZ = cms.double( 24.2 ),
        maxRegions = cms.int32( 2 ),
        EtaR_UpperLimit_Par1 = cms.double( 0.25 ),
        UseVertex = cms.bool( False ),
        Z_fixed = cms.bool( False ),
        PhiR_UpperLimit_Par1 = cms.double( 0.6 ),
        PhiR_UpperLimit_Par2 = cms.double( 0.2 ),
        Rescale_phi = cms.double( 3.0 ),
        DeltaEta = cms.double( 0.2 ),
        precise = cms.bool( True ),
        OnDemand = cms.int32( -1 ),
        EtaR_UpperLimit_Par2 = cms.double( 0.15 ),
        MeasurementTrackerName = cms.InputTag( "hltESPMeasurementTracker" ),
        vertexCollection = cms.InputTag( "pixelVertices" ),
        Pt_min = cms.double( 3.0 ),
        beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
        Phi_fixed = cms.bool( True ),
        DeltaR = cms.double( 0.025 ),
        input = cms.InputTag( 'hltL2Muons','UpdatedAtVtx' ),
        DeltaPhi = cms.double( 0.15 ),
        Phi_min = cms.double( 0.1 ),
        Rescale_eta = cms.double( 3.0 )
      ),
      TrackTransformer = cms.PSet( 
        Fitter = cms.string( "hltESPL3MuKFTrajectoryFitter" ),
        RefitDirection = cms.string( "insideOut" ),
        RefitRPCHits = cms.bool( True ),
        Propagator = cms.string( "hltESPSmartPropagatorAny" ),
        DoPredictionsOnly = cms.bool( False ),
        TrackerRecHitBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
        MuonRecHitBuilder = cms.string( "hltESPMuonTransientTrackingRecHitBuilder" ),
        Smoother = cms.string( "hltESPKFTrajectorySmootherForMuonTrackLoader" )
      ),
      tkTrajBeamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
      RefitRPCHits = cms.bool( True ),
      tkTrajVertex = cms.InputTag( "Notused" ),
      GlbRefitterParameters = cms.PSet( 
        Fitter = cms.string( "hltESPL3MuKFTrajectoryFitter" ),
        DTRecSegmentLabel = cms.InputTag( "hltDt4DSegments" ),
        RefitFlag = cms.bool( True ),
        SkipStation = cms.int32( -1 ),
        Chi2CutRPC = cms.double( 1.0 ),
        PropDirForCosmics = cms.bool( False ),
        CSCRecSegmentLabel = cms.InputTag( "hltCscSegments" ),
        GEMRecHitLabel = cms.InputTag( "hltGemRecHits" ),
        HitThreshold = cms.int32( 1 ),
        Chi2CutGEM = cms.double( 1.0 ),
        DYTthrs = cms.vint32( 30, 15 ),
        TrackerSkipSystem = cms.int32( -1 ),
        RefitDirection = cms.string( "insideOut" ),
        Chi2CutCSC = cms.double( 150.0 ),
        Chi2CutDT = cms.double( 10.0 ),
        RefitRPCHits = cms.bool( True ),
        TrackerSkipSection = cms.int32( -1 ),
        Propagator = cms.string( "hltESPSmartPropagatorAny" ),
        DoPredictionsOnly = cms.bool( False ),
        TrackerRecHitBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
        MuonHitsOption = cms.int32( 1 ),
        MuonRecHitBuilder = cms.string( "hltESPMuonTransientTrackingRecHitBuilder" )
      ),
      PCut = cms.double( 2.5 ),
      tkTrajMaxDXYBeamSpot = cms.double( 9999.0 ),
      TrackerRecHitBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
      tkTrajMaxChi2 = cms.double( 9999.0 ),
      MuonRecHitBuilder = cms.string( "hltESPMuonTransientTrackingRecHitBuilder" ),
      ScaleTECyFactor = cms.double( -1.0 ),
      tkTrajLabel = cms.InputTag( "hltIterL3OIMuonTrackSelectionHighPurity" )
    )
)
process.hltIterL3OIL3MuonsLinksCombination = cms.EDProducer( "L3TrackLinksCombiner",
    labels = cms.VInputTag( 'hltL3MuonsIterL3OI' )
)
process.hltIterL3OIL3Muons = cms.EDProducer( "L3TrackCombiner",
    labels = cms.VInputTag( 'hltL3MuonsIterL3OI' )
)
process.hltIterL3OIL3MuonCandidates = cms.EDProducer( "L3MuonCandidateProducer",
    InputObjects = cms.InputTag( "hltIterL3OIL3Muons" ),
    InputLinksObjects = cms.InputTag( "hltIterL3OIL3MuonsLinksCombination" ),
    MuonPtOption = cms.string( "Tracker" )
)
process.hltL2SelectorForL3IO = cms.EDProducer( "HLTMuonL2SelectorForL3IO",
    l2Src = cms.InputTag( 'hltL2Muons','UpdatedAtVtx' ),
    l3OISrc = cms.InputTag( "hltIterL3OIL3MuonCandidates" ),
    InputLinks = cms.InputTag( "hltIterL3OIL3MuonsLinksCombination" ),
    applyL3Filters = cms.bool( False ),
    MinNhits = cms.int32( 1 ),
    MaxNormalizedChi2 = cms.double( 20.0 ),
    MinNmuonHits = cms.int32( 1 ),
    MaxPtDifference = cms.double( 0.3 )
)
process.hltPixelTracksFitter = cms.EDProducer( "PixelFitterByHelixProjectionsProducer",
    scaleErrorsForBPix1 = cms.bool( False ),
    scaleFactor = cms.double( 0.65 )
)
process.hltPixelTracksFilter = cms.EDProducer( "PixelTrackFilterByKinematicsProducer",
    ptMin = cms.double( 0.1 ),
    nSigmaInvPtTolerance = cms.double( 0.0 ),
    tipMax = cms.double( 1.0 ),
    nSigmaTipMaxTolerance = cms.double( 0.0 ),
    chi2 = cms.double( 1000.0 )
)
process.hltPixelTracksTrackingRegions = cms.EDProducer( "GlobalTrackingRegionFromBeamSpotEDProducer",
    RegionPSet = cms.PSet( 
      nSigmaZ = cms.double( 4.0 ),
      beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
      ptMin = cms.double( 0.8 ),
      originRadius = cms.double( 0.02 ),
      precise = cms.bool( True )
    )
)
process.hltPixelTracks = cms.EDProducer( "PixelTrackProducerFromSoA",
    beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
    trackSrc = cms.InputTag( "hltPixelTracksSoA" ),
    pixelRecHitLegacySrc = cms.InputTag( "hltSiPixelRecHits" ),
    minNumberOfHits = cms.int32( 0 ),
    minQuality = cms.string( "loose" )
)
process.hltPixelVertices = cms.EDProducer( "PixelVertexProducerFromSoA",
    TrackCollection = cms.InputTag( "hltPixelTracks" ),
    beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
    src = cms.InputTag( "hltPixelVerticesSoA" )
)
process.hltTrimmedPixelVertices = cms.EDProducer( "PixelVertexCollectionTrimmer",
    src = cms.InputTag( "hltPixelVertices" ),
    maxVtx = cms.uint32( 100 ),
    fractionSumPt2 = cms.double( 0.3 ),
    minSumPt2 = cms.double( 0.0 ),
    PVcomparer = cms.PSet(  refToPSet_ = cms.string( "HLTPSetPvClusterComparerForIT" ) )
)
process.hltIterL3MuonPixelTracksTrackingRegions = cms.EDProducer( "MuonTrackingRegionByPtEDProducer",
    DeltaR = cms.double( 0.025 ),
    beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
    OnDemand = cms.int32( -1 ),
    vertexCollection = cms.InputTag( "notUsed" ),
    MeasurementTrackerName = cms.InputTag( "" ),
    UseVertex = cms.bool( False ),
    Rescale_Dz = cms.double( 4.0 ),
    Pt_fixed = cms.bool( True ),
    Z_fixed = cms.bool( True ),
    Pt_min = cms.double( 0.0 ),
    DeltaZ = cms.double( 24.2 ),
    ptRanges = cms.vdouble( 0.0, 15.0, 20.0, 1.0E64 ),
    deltaEtas = cms.vdouble( 0.2, 0.2, 0.2 ),
    deltaPhis = cms.vdouble( 0.75, 0.45, 0.225 ),
    maxRegions = cms.int32( 5 ),
    precise = cms.bool( True ),
    input = cms.InputTag( "hltL2SelectorForL3IO" )
)
process.hltPixelTracksInRegionL2 = cms.EDProducer( "TrackSelectorByRegion",
    tracks = cms.InputTag( "hltPixelTracks" ),
    regions = cms.InputTag( "hltIterL3MuonPixelTracksTrackingRegions" ),
    produceTrackCollection = cms.bool( True ),
    produceMask = cms.bool( False )
)
process.hltIter0IterL3MuonPixelSeedsFromPixelTracks = cms.EDProducer( "SeedGeneratorFromProtoTracksEDProducer",
    InputCollection = cms.InputTag( "hltPixelTracksInRegionL2" ),
    InputVertexCollection = cms.InputTag( "" ),
    originHalfLength = cms.double( 0.3 ),
    originRadius = cms.double( 0.1 ),
    useProtoTrackKinematics = cms.bool( False ),
    useEventsWithNoVertex = cms.bool( True ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    usePV = cms.bool( False ),
    includeFourthHit = cms.bool( True ),
    SeedCreatorPSet = cms.PSet(  refToPSet_ = cms.string( "HLTSeedFromProtoTracks" ) )
)
process.hltIter0IterL3MuonPixelSeedsFromPixelTracksFiltered = cms.EDProducer( "MuonHLTSeedMVAClassifier",
    src = cms.InputTag( "hltIter0IterL3MuonPixelSeedsFromPixelTracks" ),
    L1Muon = cms.InputTag( 'hltGtStage2Digis','Muon' ),
    L2Muon = cms.InputTag( "hltL2MuonCandidates" ),
    rejectAll = cms.bool( False ),
    isFromL1 = cms.bool( False ),
    mvaFileBL1 = cms.FileInPath( "RecoMuon/TrackerSeedGenerator/data/xgb_Run3_Iter0FromL1_PatatrackSeeds_barrel_v2.xml" ),
    mvaFileEL1 = cms.FileInPath( "RecoMuon/TrackerSeedGenerator/data/xgb_Run3_Iter0FromL1_PatatrackSeeds_endcap_v2.xml" ),
    mvaFileBL2 = cms.FileInPath( "RecoMuon/TrackerSeedGenerator/data/xgb_Run3_Iter0_PatatrackSeeds_barrel_v2.xml" ),
    mvaFileEL2 = cms.FileInPath( "RecoMuon/TrackerSeedGenerator/data/xgb_Run3_Iter0_PatatrackSeeds_endcap_v2.xml" ),
    mvaScaleMeanBL1 = cms.vdouble( 5.000589710660383E-4, 3.90864688207247E-6, 5.566857527819724E-6, 2.919765272506742E-5, 0.0020554125627069656, 3.7638302502636416E-4, 0.17319245569742575, -0.001928435311705264 ),
    mvaScaleStdBL1 = cms.vdouble( 0.0024104956748896007, 2.65014770862918E-4, 7.474433025576554E-4, 0.07076843921543456, 0.8267930143339328, 0.5908013434667966, 0.16945213299167364, 0.19343332776069666 ),
    mvaScaleMeanEL1 = cms.vdouble( 4.633287331326438E-4, 4.336517990150388E-6, 1.1353478376025763E-5, -9.703999978722195E-4, -0.016733916560431213, 0.0010338859165970008, 0.1664362346287898, -0.002093062022327288 ),
    mvaScaleStdEL1 = cms.vdouble( 0.0017591716932616446, 4.3554653252872314E-4, 0.0018509069108135767, 0.3154887917554714, 1.0554267175028256, 0.46607446165391897, 0.1599804481244901, 0.1802872646042384 ),
    mvaScaleMeanBL2 = cms.vdouble( 5.535405438531338E-4, 5.803137756667417E-6, 9.844857031022169E-6, 8.572205780682784E-6, 0.0015762679722632618, 4.395397563024067E-4, 0.1374567309015325, -0.005115438205054733, 0.11749212456078427, 6.10630569023574E-5 ),
    mvaScaleStdBL2 = cms.vdouble( 0.002626537945220411, 0.004217812355595011, 0.01003037371073428, 0.08030848788317736, 0.8112437912866772, 0.682723351939253, 0.14357759358396366, 0.1714659927555431, 0.12492124275832932, 0.14997195939803123 ),
    mvaScaleMeanEL2 = cms.vdouble( 4.7677144995220935E-4, 3.6745265753688003E-6, 8.865597384997202E-6, 3.035250572254308E-4, -0.017277139191991336, 7.614700645923214E-4, 0.13098441610679598, -0.004062335797232639, 0.10313449419117496, 8.091070691633311E-5 ),
    mvaScaleStdEL2 = cms.vdouble( 0.0017644245219973625, 1.4737438911735824E-4, 5.815417377127688E-4, 0.32486307205392545, 0.9554477166210787, 0.6368300706332602, 0.18228633778650377, 0.17472670187127687, 0.1069305400952516, 0.12396350430658183 ),
    doSort = cms.bool( False ),
    nSeedsMaxB = cms.int32( 99999 ),
    nSeedsMaxE = cms.int32( 99999 ),
    etaEdge = cms.double( 1.2 ),
    mvaCutB = cms.double( 0.01 ),
    mvaCutE = cms.double( 0.01 ),
    minL1Qual = cms.int32( 7 ),
    baseScore = cms.double( 0.5 )
)
process.hltIter0IterL3MuonCkfTrackCandidates = cms.EDProducer( "CkfTrackCandidateMaker",
    cleanTrajectoryAfterInOut = cms.bool( False ),
    doSeedingRegionRebuilding = cms.bool( True ),
    onlyPixelHitsForSeedCleaner = cms.bool( False ),
    reverseTrajectories = cms.bool( False ),
    useHitsSplitting = cms.bool( True ),
    MeasurementTrackerEvent = cms.InputTag( "hltSiStripClusters" ),
    src = cms.InputTag( "hltIter0IterL3MuonPixelSeedsFromPixelTracksFiltered" ),
    clustersToSkip = cms.InputTag( "" ),
    phase2clustersToSkip = cms.InputTag( "" ),
    TrajectoryBuilderPSet = cms.PSet(  refToPSet_ = cms.string( "HLTIter0IterL3MuonPSetGroupedCkfTrajectoryBuilderIT" ) ),
    TransientInitialStateEstimatorParameters = cms.PSet( 
      propagatorAlongTISE = cms.string( "PropagatorWithMaterialParabolicMf" ),
      numberMeasurementsForFit = cms.int32( 4 ),
      propagatorOppositeTISE = cms.string( "PropagatorWithMaterialParabolicMfOpposite" )
    ),
    numHitsForSeedCleaner = cms.int32( 4 ),
    NavigationSchool = cms.string( "SimpleNavigationSchool" ),
    RedundantSeedCleaner = cms.string( "none" ),
    TrajectoryCleaner = cms.string( "hltESPTrajectoryCleanerBySharedHits" ),
    maxNSeeds = cms.uint32( 100000 ),
    maxSeedsBeforeCleaning = cms.uint32( 1000 )
)
process.hltIter0IterL3MuonCtfWithMaterialTracks = cms.EDProducer( "TrackProducer",
    useSimpleMF = cms.bool( True ),
    SimpleMagneticField = cms.string( "ParabolicMf" ),
    src = cms.InputTag( "hltIter0IterL3MuonCkfTrackCandidates" ),
    clusterRemovalInfo = cms.InputTag( "" ),
    beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
    Fitter = cms.string( "hltESPFittingSmootherIT" ),
    useHitsSplitting = cms.bool( False ),
    alias = cms.untracked.string( "ctfWithMaterialTracks" ),
    TrajectoryInEvent = cms.bool( False ),
    TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
    AlgorithmName = cms.string( "hltIter0" ),
    Propagator = cms.string( "hltESPRungeKuttaTrackerPropagator" ),
    GeometricInnerState = cms.bool( True ),
    NavigationSchool = cms.string( "" ),
    MeasurementTracker = cms.string( "" ),
    MeasurementTrackerEvent = cms.InputTag( "hltSiStripClusters" )
)
process.hltIter0IterL3MuonTrackCutClassifier = cms.EDProducer( "TrackCutClassifier",
    src = cms.InputTag( "hltIter0IterL3MuonCtfWithMaterialTracks" ),
    beamspot = cms.InputTag( "hltOnlineBeamSpot" ),
    vertices = cms.InputTag( "hltTrimmedPixelVertices" ),
    ignoreVertices = cms.bool( False ),
    qualityCuts = cms.vdouble( -0.7, 0.1, 0.7 ),
    mva = cms.PSet( 
      minPixelHits = cms.vint32( 0, 0, 0 ),
      maxDzWrtBS = cms.vdouble( 3.40282346639E38, 24.0, 100.0 ),
      dr_par = cms.PSet( 
        d0err = cms.vdouble( 0.003, 0.003, 3.40282346639E38 ),
        dr_par2 = cms.vdouble( 0.3, 0.3, 3.40282346639E38 ),
        dr_par1 = cms.vdouble( 0.4, 0.4, 3.40282346639E38 ),
        dr_exp = cms.vint32( 4, 4, 2147483647 ),
        d0err_par = cms.vdouble( 0.001, 0.001, 3.40282346639E38 )
      ),
      maxLostLayers = cms.vint32( 1, 1, 1 ),
      min3DLayers = cms.vint32( 0, 0, 0 ),
      dz_par = cms.PSet( 
        dz_par1 = cms.vdouble( 0.4, 0.4, 3.40282346639E38 ),
        dz_par2 = cms.vdouble( 0.35, 0.35, 3.40282346639E38 ),
        dz_exp = cms.vint32( 4, 4, 2147483647 )
      ),
      minNVtxTrk = cms.int32( 3 ),
      maxDz = cms.vdouble( 0.5, 0.2, 3.40282346639E38 ),
      minNdof = cms.vdouble( 1.0E-5, 1.0E-5, 1.0E-5 ),
      maxChi2 = cms.vdouble( 3.40282346639E38, 3.40282346639E38, 3.40282346639E38 ),
      maxChi2n = cms.vdouble( 1.2, 1.0, 0.7 ),
      maxDr = cms.vdouble( 0.5, 0.03, 3.40282346639E38 ),
      minLayers = cms.vint32( 3, 3, 4 )
    )
)
process.hltIter0IterL3MuonTrackSelectionHighPurity = cms.EDProducer( "TrackCollectionFilterCloner",
    originalSource = cms.InputTag( "hltIter0IterL3MuonCtfWithMaterialTracks" ),
    originalMVAVals = cms.InputTag( 'hltIter0IterL3MuonTrackCutClassifier','MVAValues' ),
    originalQualVals = cms.InputTag( 'hltIter0IterL3MuonTrackCutClassifier','QualityMasks' ),
    minQuality = cms.string( "highPurity" ),
    copyExtras = cms.untracked.bool( True ),
    copyTrajectories = cms.untracked.bool( False )
)
process.hltL3MuonsIterL3IO = cms.EDProducer( "L3MuonProducer",
    ServiceParameters = cms.PSet( 
      RPCLayers = cms.bool( True ),
      UseMuonNavigation = cms.untracked.bool( True ),
      Propagators = cms.untracked.vstring( 'hltESPSmartPropagatorAny',
        'SteppingHelixPropagatorAny',
        'hltESPSmartPropagator',
        'hltESPSteppingHelixPropagatorOpposite' )
    ),
    MuonCollectionLabel = cms.InputTag( 'hltL2Muons','UpdatedAtVtx' ),
    TrackLoaderParameters = cms.PSet( 
      MuonSeededTracksInstance = cms.untracked.string( "L2Seeded" ),
      beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
      DoSmoothing = cms.bool( False ),
      SmoothTkTrack = cms.untracked.bool( False ),
      VertexConstraint = cms.bool( False ),
      MuonUpdatorAtVertexParameters = cms.PSet( 
        MaxChi2 = cms.double( 1000000.0 ),
        BeamSpotPositionErrors = cms.vdouble( 0.1, 0.1, 5.3 ),
        Propagator = cms.string( "hltESPSteppingHelixPropagatorOpposite" )
      ),
      PutTkTrackIntoEvent = cms.untracked.bool( False ),
      Smoother = cms.string( "hltESPKFTrajectorySmootherForMuonTrackLoader" )
    ),
    L3TrajBuilderParameters = cms.PSet( 
      PtCut = cms.double( 1.0 ),
      TrackerPropagator = cms.string( "SteppingHelixPropagatorAny" ),
      GlobalMuonTrackMatcher = cms.PSet( 
        Chi2Cut_3 = cms.double( 200.0 ),
        DeltaDCut_2 = cms.double( 10.0 ),
        Eta_threshold = cms.double( 1.2 ),
        Quality_2 = cms.double( 15.0 ),
        DeltaDCut_1 = cms.double( 40.0 ),
        Quality_3 = cms.double( 7.0 ),
        DeltaDCut_3 = cms.double( 15.0 ),
        Quality_1 = cms.double( 20.0 ),
        Pt_threshold1 = cms.double( 0.0 ),
        DeltaRCut_2 = cms.double( 0.2 ),
        DeltaRCut_1 = cms.double( 0.1 ),
        Pt_threshold2 = cms.double( 9.99999999E8 ),
        Chi2Cut_1 = cms.double( 50.0 ),
        Chi2Cut_2 = cms.double( 50.0 ),
        DeltaRCut_3 = cms.double( 1.0 ),
        LocChi2Cut = cms.double( 0.001 ),
        Propagator = cms.string( "hltESPSmartPropagator" ),
        MinPt = cms.double( 1.0 ),
        MinP = cms.double( 2.5 )
      ),
      ScaleTECxFactor = cms.double( -1.0 ),
      tkTrajUseVertex = cms.bool( False ),
      MuonTrackingRegionBuilder = cms.PSet( 
        Rescale_Dz = cms.double( 4.0 ),
        Pt_fixed = cms.bool( True ),
        Eta_fixed = cms.bool( True ),
        Eta_min = cms.double( 0.1 ),
        DeltaZ = cms.double( 24.2 ),
        maxRegions = cms.int32( 2 ),
        EtaR_UpperLimit_Par1 = cms.double( 0.25 ),
        UseVertex = cms.bool( False ),
        Z_fixed = cms.bool( True ),
        PhiR_UpperLimit_Par1 = cms.double( 0.6 ),
        PhiR_UpperLimit_Par2 = cms.double( 0.2 ),
        Rescale_phi = cms.double( 3.0 ),
        DeltaEta = cms.double( 0.04 ),
        precise = cms.bool( True ),
        OnDemand = cms.int32( -1 ),
        EtaR_UpperLimit_Par2 = cms.double( 0.15 ),
        MeasurementTrackerName = cms.InputTag( "hltESPMeasurementTracker" ),
        vertexCollection = cms.InputTag( "pixelVertices" ),
        Pt_min = cms.double( 3.0 ),
        beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
        Phi_fixed = cms.bool( True ),
        DeltaR = cms.double( 0.025 ),
        input = cms.InputTag( "hltL2SelectorForL3IO" ),
        DeltaPhi = cms.double( 0.15 ),
        Phi_min = cms.double( 0.1 ),
        Rescale_eta = cms.double( 3.0 )
      ),
      TrackTransformer = cms.PSet( 
        Fitter = cms.string( "hltESPL3MuKFTrajectoryFitter" ),
        RefitDirection = cms.string( "insideOut" ),
        RefitRPCHits = cms.bool( True ),
        Propagator = cms.string( "hltESPSmartPropagatorAny" ),
        DoPredictionsOnly = cms.bool( False ),
        TrackerRecHitBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
        MuonRecHitBuilder = cms.string( "hltESPMuonTransientTrackingRecHitBuilder" ),
        Smoother = cms.string( "hltESPKFTrajectorySmootherForMuonTrackLoader" )
      ),
      tkTrajBeamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
      RefitRPCHits = cms.bool( True ),
      tkTrajVertex = cms.InputTag( "hltTrimmedPixelVertices" ),
      GlbRefitterParameters = cms.PSet( 
        Fitter = cms.string( "hltESPL3MuKFTrajectoryFitter" ),
        DTRecSegmentLabel = cms.InputTag( "hltDt4DSegments" ),
        RefitFlag = cms.bool( True ),
        SkipStation = cms.int32( -1 ),
        Chi2CutRPC = cms.double( 1.0 ),
        PropDirForCosmics = cms.bool( False ),
        CSCRecSegmentLabel = cms.InputTag( "hltCscSegments" ),
        GEMRecHitLabel = cms.InputTag( "hltGemRecHits" ),
        HitThreshold = cms.int32( 1 ),
        Chi2CutGEM = cms.double( 1.0 ),
        DYTthrs = cms.vint32( 30, 15 ),
        TrackerSkipSystem = cms.int32( -1 ),
        RefitDirection = cms.string( "insideOut" ),
        Chi2CutCSC = cms.double( 150.0 ),
        Chi2CutDT = cms.double( 10.0 ),
        RefitRPCHits = cms.bool( True ),
        TrackerSkipSection = cms.int32( -1 ),
        Propagator = cms.string( "hltESPSmartPropagatorAny" ),
        DoPredictionsOnly = cms.bool( False ),
        TrackerRecHitBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
        MuonHitsOption = cms.int32( 1 ),
        MuonRecHitBuilder = cms.string( "hltESPMuonTransientTrackingRecHitBuilder" )
      ),
      PCut = cms.double( 2.5 ),
      tkTrajMaxDXYBeamSpot = cms.double( 9999.0 ),
      TrackerRecHitBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
      matchToSeeds = cms.bool( True ),
      tkTrajMaxChi2 = cms.double( 9999.0 ),
      MuonRecHitBuilder = cms.string( "hltESPMuonTransientTrackingRecHitBuilder" ),
      ScaleTECyFactor = cms.double( -1.0 ),
      tkTrajLabel = cms.InputTag( "hltIter0IterL3MuonTrackSelectionHighPurity" )
    )
)
process.hltIterL3MuonsFromL2LinksCombination = cms.EDProducer( "L3TrackLinksCombiner",
    labels = cms.VInputTag( 'hltL3MuonsIterL3OI','hltL3MuonsIterL3IO' )
)
process.hltL1MuonsPt0 = cms.EDProducer( "HLTL1TMuonSelector",
    InputObjects = cms.InputTag( 'hltGtStage2Digis','Muon' ),
    L1MinPt = cms.double( -1.0 ),
    L1MaxEta = cms.double( 5.0 ),
    L1MinQuality = cms.uint32( 7 ),
    CentralBxOnly = cms.bool( True )
)
process.hltIterL3FromL1MuonPixelTracksTrackingRegions = cms.EDProducer( "L1MuonSeededTrackingRegionsEDProducer",
    Propagator = cms.string( "SteppingHelixPropagatorAny" ),
    L1MinPt = cms.double( 0.0 ),
    L1MaxEta = cms.double( 2.5 ),
    L1MinQuality = cms.uint32( 7 ),
    SetMinPtBarrelTo = cms.double( 3.5 ),
    SetMinPtEndcapTo = cms.double( 1.0 ),
    CentralBxOnly = cms.bool( True ),
    RegionPSet = cms.PSet( 
      vertexCollection = cms.InputTag( "notUsed" ),
      deltaEtas = cms.vdouble( 0.35, 0.35, 0.35, 0.35 ),
      zErrorVetex = cms.double( 0.2 ),
      beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
      zErrorBeamSpot = cms.double( 24.2 ),
      maxNVertices = cms.int32( 1 ),
      maxNRegions = cms.int32( 5 ),
      nSigmaZVertex = cms.double( 3.0 ),
      nSigmaZBeamSpot = cms.double( 4.0 ),
      ptMin = cms.double( 0.0 ),
      mode = cms.string( "BeamSpotSigma" ),
      input = cms.InputTag( "hltL1MuonsPt0" ),
      ptRanges = cms.vdouble( 0.0, 10.0, 15.0, 20.0, 1.0E64 ),
      searchOpt = cms.bool( False ),
      deltaPhis = cms.vdouble( 1.0, 0.8, 0.6, 0.3 ),
      whereToUseMeasurementTracker = cms.string( "Never" ),
      originRadius = cms.double( 0.2 ),
      measurementTrackerName = cms.InputTag( "" ),
      precise = cms.bool( True )
    ),
    ServiceParameters = cms.PSet( 
      RPCLayers = cms.bool( True ),
      UseMuonNavigation = cms.untracked.bool( True ),
      Propagators = cms.untracked.vstring( 'SteppingHelixPropagatorAny' )
    )
)
process.hltPixelTracksInRegionL1 = cms.EDProducer( "TrackSelectorByRegion",
    tracks = cms.InputTag( "hltPixelTracks" ),
    regions = cms.InputTag( "hltIterL3FromL1MuonPixelTracksTrackingRegions" ),
    produceTrackCollection = cms.bool( True ),
    produceMask = cms.bool( False )
)
process.hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks = cms.EDProducer( "SeedGeneratorFromProtoTracksEDProducer",
    InputCollection = cms.InputTag( "hltPixelTracksInRegionL1" ),
    InputVertexCollection = cms.InputTag( "" ),
    originHalfLength = cms.double( 0.3 ),
    originRadius = cms.double( 0.1 ),
    useProtoTrackKinematics = cms.bool( False ),
    useEventsWithNoVertex = cms.bool( True ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    usePV = cms.bool( False ),
    includeFourthHit = cms.bool( True ),
    SeedCreatorPSet = cms.PSet(  refToPSet_ = cms.string( "HLTSeedFromProtoTracks" ) )
)
process.hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracksFiltered = cms.EDProducer( "MuonHLTSeedMVAClassifier",
    src = cms.InputTag( "hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks" ),
    L1Muon = cms.InputTag( 'hltGtStage2Digis','Muon' ),
    L2Muon = cms.InputTag( "hltL2MuonCandidates" ),
    rejectAll = cms.bool( False ),
    isFromL1 = cms.bool( True ),
    mvaFileBL1 = cms.FileInPath( "RecoMuon/TrackerSeedGenerator/data/xgb_Run3_Iter0FromL1_PatatrackSeeds_barrel_v2.xml" ),
    mvaFileEL1 = cms.FileInPath( "RecoMuon/TrackerSeedGenerator/data/xgb_Run3_Iter0FromL1_PatatrackSeeds_endcap_v2.xml" ),
    mvaFileBL2 = cms.FileInPath( "RecoMuon/TrackerSeedGenerator/data/xgb_Run3_Iter0_PatatrackSeeds_barrel_v2.xml" ),
    mvaFileEL2 = cms.FileInPath( "RecoMuon/TrackerSeedGenerator/data/xgb_Run3_Iter0_PatatrackSeeds_endcap_v2.xml" ),
    mvaScaleMeanBL1 = cms.vdouble( 5.000589710660383E-4, 3.90864688207247E-6, 5.566857527819724E-6, 2.919765272506742E-5, 0.0020554125627069656, 3.7638302502636416E-4, 0.17319245569742575, -0.001928435311705264 ),
    mvaScaleStdBL1 = cms.vdouble( 0.0024104956748896007, 2.65014770862918E-4, 7.474433025576554E-4, 0.07076843921543456, 0.8267930143339328, 0.5908013434667966, 0.16945213299167364, 0.19343332776069666 ),
    mvaScaleMeanEL1 = cms.vdouble( 4.633287331326438E-4, 4.336517990150388E-6, 1.1353478376025763E-5, -9.703999978722195E-4, -0.016733916560431213, 0.0010338859165970008, 0.1664362346287898, -0.002093062022327288 ),
    mvaScaleStdEL1 = cms.vdouble( 0.0017591716932616446, 4.3554653252872314E-4, 0.0018509069108135767, 0.3154887917554714, 1.0554267175028256, 0.46607446165391897, 0.1599804481244901, 0.1802872646042384 ),
    mvaScaleMeanBL2 = cms.vdouble( 5.535405438531338E-4, 5.803137756667417E-6, 9.844857031022169E-6, 8.572205780682784E-6, 0.0015762679722632618, 4.395397563024067E-4, 0.1374567309015325, -0.005115438205054733, 0.11749212456078427, 6.10630569023574E-5 ),
    mvaScaleStdBL2 = cms.vdouble( 0.002626537945220411, 0.004217812355595011, 0.01003037371073428, 0.08030848788317736, 0.8112437912866772, 0.682723351939253, 0.14357759358396366, 0.1714659927555431, 0.12492124275832932, 0.14997195939803123 ),
    mvaScaleMeanEL2 = cms.vdouble( 4.7677144995220935E-4, 3.6745265753688003E-6, 8.865597384997202E-6, 3.035250572254308E-4, -0.017277139191991336, 7.614700645923214E-4, 0.13098441610679598, -0.004062335797232639, 0.10313449419117496, 8.091070691633311E-5 ),
    mvaScaleStdEL2 = cms.vdouble( 0.0017644245219973625, 1.4737438911735824E-4, 5.815417377127688E-4, 0.32486307205392545, 0.9554477166210787, 0.6368300706332602, 0.18228633778650377, 0.17472670187127687, 0.1069305400952516, 0.12396350430658183 ),
    doSort = cms.bool( False ),
    nSeedsMaxB = cms.int32( 99999 ),
    nSeedsMaxE = cms.int32( 99999 ),
    etaEdge = cms.double( 1.2 ),
    mvaCutB = cms.double( 0.01 ),
    mvaCutE = cms.double( 0.01 ),
    minL1Qual = cms.int32( 7 ),
    baseScore = cms.double( 0.5 )
)
process.hltIter0IterL3FromL1MuonCkfTrackCandidates = cms.EDProducer( "CkfTrackCandidateMaker",
    cleanTrajectoryAfterInOut = cms.bool( False ),
    doSeedingRegionRebuilding = cms.bool( True ),
    onlyPixelHitsForSeedCleaner = cms.bool( False ),
    reverseTrajectories = cms.bool( False ),
    useHitsSplitting = cms.bool( True ),
    MeasurementTrackerEvent = cms.InputTag( "hltSiStripClusters" ),
    src = cms.InputTag( "hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracksFiltered" ),
    clustersToSkip = cms.InputTag( "" ),
    phase2clustersToSkip = cms.InputTag( "" ),
    TrajectoryBuilderPSet = cms.PSet(  refToPSet_ = cms.string( "HLTIter0IterL3FromL1MuonPSetGroupedCkfTrajectoryBuilderIT" ) ),
    TransientInitialStateEstimatorParameters = cms.PSet( 
      propagatorAlongTISE = cms.string( "PropagatorWithMaterialParabolicMf" ),
      numberMeasurementsForFit = cms.int32( 4 ),
      propagatorOppositeTISE = cms.string( "PropagatorWithMaterialParabolicMfOpposite" )
    ),
    numHitsForSeedCleaner = cms.int32( 4 ),
    NavigationSchool = cms.string( "SimpleNavigationSchool" ),
    RedundantSeedCleaner = cms.string( "none" ),
    TrajectoryCleaner = cms.string( "hltESPTrajectoryCleanerBySharedHits" ),
    maxNSeeds = cms.uint32( 100000 ),
    maxSeedsBeforeCleaning = cms.uint32( 1000 )
)
process.hltIter0IterL3FromL1MuonCtfWithMaterialTracks = cms.EDProducer( "TrackProducer",
    useSimpleMF = cms.bool( True ),
    SimpleMagneticField = cms.string( "ParabolicMf" ),
    src = cms.InputTag( "hltIter0IterL3FromL1MuonCkfTrackCandidates" ),
    clusterRemovalInfo = cms.InputTag( "" ),
    beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
    Fitter = cms.string( "hltESPFittingSmootherIT" ),
    useHitsSplitting = cms.bool( False ),
    alias = cms.untracked.string( "ctfWithMaterialTracks" ),
    TrajectoryInEvent = cms.bool( False ),
    TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
    AlgorithmName = cms.string( "hltIter0" ),
    Propagator = cms.string( "hltESPRungeKuttaTrackerPropagator" ),
    GeometricInnerState = cms.bool( True ),
    NavigationSchool = cms.string( "" ),
    MeasurementTracker = cms.string( "" ),
    MeasurementTrackerEvent = cms.InputTag( "hltSiStripClusters" )
)
process.hltIter0IterL3FromL1MuonTrackCutClassifier = cms.EDProducer( "TrackCutClassifier",
    src = cms.InputTag( "hltIter0IterL3FromL1MuonCtfWithMaterialTracks" ),
    beamspot = cms.InputTag( "hltOnlineBeamSpot" ),
    vertices = cms.InputTag( "hltTrimmedPixelVertices" ),
    ignoreVertices = cms.bool( False ),
    qualityCuts = cms.vdouble( -0.7, 0.1, 0.7 ),
    mva = cms.PSet( 
      minPixelHits = cms.vint32( 0, 0, 0 ),
      maxDzWrtBS = cms.vdouble( 3.40282346639E38, 24.0, 100.0 ),
      dr_par = cms.PSet( 
        d0err = cms.vdouble( 0.003, 0.003, 3.40282346639E38 ),
        dr_par2 = cms.vdouble( 0.3, 0.3, 3.40282346639E38 ),
        dr_par1 = cms.vdouble( 0.4, 0.4, 3.40282346639E38 ),
        dr_exp = cms.vint32( 4, 4, 2147483647 ),
        d0err_par = cms.vdouble( 0.001, 0.001, 3.40282346639E38 )
      ),
      maxLostLayers = cms.vint32( 1, 1, 1 ),
      min3DLayers = cms.vint32( 0, 0, 0 ),
      dz_par = cms.PSet( 
        dz_par1 = cms.vdouble( 0.4, 0.4, 3.40282346639E38 ),
        dz_par2 = cms.vdouble( 0.35, 0.35, 3.40282346639E38 ),
        dz_exp = cms.vint32( 4, 4, 2147483647 )
      ),
      minNVtxTrk = cms.int32( 3 ),
      maxDz = cms.vdouble( 0.5, 0.2, 3.40282346639E38 ),
      minNdof = cms.vdouble( 1.0E-5, 1.0E-5, 1.0E-5 ),
      maxChi2 = cms.vdouble( 3.40282346639E38, 3.40282346639E38, 3.40282346639E38 ),
      maxChi2n = cms.vdouble( 1.2, 1.0, 0.7 ),
      maxDr = cms.vdouble( 0.5, 0.03, 3.40282346639E38 ),
      minLayers = cms.vint32( 3, 3, 4 )
    )
)
process.hltIter0IterL3FromL1MuonTrackSelectionHighPurity = cms.EDProducer( "TrackCollectionFilterCloner",
    originalSource = cms.InputTag( "hltIter0IterL3FromL1MuonCtfWithMaterialTracks" ),
    originalMVAVals = cms.InputTag( 'hltIter0IterL3FromL1MuonTrackCutClassifier','MVAValues' ),
    originalQualVals = cms.InputTag( 'hltIter0IterL3FromL1MuonTrackCutClassifier','QualityMasks' ),
    minQuality = cms.string( "highPurity" ),
    copyExtras = cms.untracked.bool( True ),
    copyTrajectories = cms.untracked.bool( False )
)
process.hltIterL3MuonMerged = cms.EDProducer( "TrackListMerger",
    ShareFrac = cms.double( 0.19 ),
    FoundHitBonus = cms.double( 5.0 ),
    LostHitPenalty = cms.double( 20.0 ),
    MinPT = cms.double( 0.05 ),
    Epsilon = cms.double( -0.001 ),
    MaxNormalizedChisq = cms.double( 1000.0 ),
    MinFound = cms.int32( 3 ),
    TrackProducers = cms.VInputTag( 'hltIterL3OIMuonTrackSelectionHighPurity','hltIter0IterL3MuonTrackSelectionHighPurity' ),
    hasSelector = cms.vint32( 0, 0 ),
    indivShareFrac = cms.vdouble( 1.0, 1.0 ),
    selectedTrackQuals = cms.VInputTag( 'hltIterL3OIMuonTrackSelectionHighPurity','hltIter0IterL3MuonTrackSelectionHighPurity' ),
    setsToMerge = cms.VPSet( 
      cms.PSet(  pQual = cms.bool( False ),
        tLists = cms.vint32( 0, 1 )
      )
    ),
    trackAlgoPriorityOrder = cms.string( "hltESPTrackAlgoPriorityOrder" ),
    allowFirstHitShare = cms.bool( True ),
    newQuality = cms.string( "confirmed" ),
    copyExtras = cms.untracked.bool( True ),
    writeOnlyTrkQuals = cms.bool( False ),
    copyMVA = cms.bool( False )
)
process.hltIterL3MuonAndMuonFromL1Merged = cms.EDProducer( "TrackListMerger",
    ShareFrac = cms.double( 0.19 ),
    FoundHitBonus = cms.double( 5.0 ),
    LostHitPenalty = cms.double( 20.0 ),
    MinPT = cms.double( 0.05 ),
    Epsilon = cms.double( -0.001 ),
    MaxNormalizedChisq = cms.double( 1000.0 ),
    MinFound = cms.int32( 3 ),
    TrackProducers = cms.VInputTag( 'hltIterL3MuonMerged','hltIter0IterL3FromL1MuonTrackSelectionHighPurity' ),
    hasSelector = cms.vint32( 0, 0 ),
    indivShareFrac = cms.vdouble( 1.0, 1.0 ),
    selectedTrackQuals = cms.VInputTag( 'hltIterL3MuonMerged','hltIter0IterL3FromL1MuonTrackSelectionHighPurity' ),
    setsToMerge = cms.VPSet( 
      cms.PSet(  pQual = cms.bool( False ),
        tLists = cms.vint32( 0, 1 )
      )
    ),
    trackAlgoPriorityOrder = cms.string( "hltESPTrackAlgoPriorityOrder" ),
    allowFirstHitShare = cms.bool( True ),
    newQuality = cms.string( "confirmed" ),
    copyExtras = cms.untracked.bool( True ),
    writeOnlyTrkQuals = cms.bool( False ),
    copyMVA = cms.bool( False )
)
process.hltIterL3GlbMuon = cms.EDProducer( "L3MuonProducer",
    ServiceParameters = cms.PSet( 
      RPCLayers = cms.bool( True ),
      UseMuonNavigation = cms.untracked.bool( True ),
      Propagators = cms.untracked.vstring( 'hltESPSmartPropagatorAny',
        'SteppingHelixPropagatorAny',
        'hltESPSmartPropagator',
        'hltESPSteppingHelixPropagatorOpposite' )
    ),
    MuonCollectionLabel = cms.InputTag( 'hltL2Muons','UpdatedAtVtx' ),
    TrackLoaderParameters = cms.PSet( 
      MuonSeededTracksInstance = cms.untracked.string( "L2Seeded" ),
      TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
      beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
      DoSmoothing = cms.bool( True ),
      SmoothTkTrack = cms.untracked.bool( False ),
      VertexConstraint = cms.bool( False ),
      MuonUpdatorAtVertexParameters = cms.PSet( 
        MaxChi2 = cms.double( 1000000.0 ),
        BeamSpotPositionErrors = cms.vdouble( 0.1, 0.1, 5.3 ),
        Propagator = cms.string( "hltESPSteppingHelixPropagatorOpposite" )
      ),
      PutTkTrackIntoEvent = cms.untracked.bool( False ),
      Smoother = cms.string( "hltESPKFTrajectorySmootherForMuonTrackLoader" )
    ),
    L3TrajBuilderParameters = cms.PSet( 
      PtCut = cms.double( 1.0 ),
      TrackerPropagator = cms.string( "SteppingHelixPropagatorAny" ),
      GlobalMuonTrackMatcher = cms.PSet( 
        Chi2Cut_3 = cms.double( 200.0 ),
        DeltaDCut_2 = cms.double( 10.0 ),
        Eta_threshold = cms.double( 1.2 ),
        Quality_2 = cms.double( 15.0 ),
        DeltaDCut_1 = cms.double( 40.0 ),
        Quality_3 = cms.double( 7.0 ),
        DeltaDCut_3 = cms.double( 15.0 ),
        Quality_1 = cms.double( 20.0 ),
        Pt_threshold1 = cms.double( 0.0 ),
        DeltaRCut_2 = cms.double( 0.2 ),
        DeltaRCut_1 = cms.double( 0.1 ),
        Pt_threshold2 = cms.double( 9.99999999E8 ),
        Chi2Cut_1 = cms.double( 50.0 ),
        Chi2Cut_2 = cms.double( 50.0 ),
        DeltaRCut_3 = cms.double( 1.0 ),
        LocChi2Cut = cms.double( 0.001 ),
        Propagator = cms.string( "hltESPSmartPropagator" ),
        MinPt = cms.double( 1.0 ),
        MinP = cms.double( 2.5 )
      ),
      ScaleTECxFactor = cms.double( -1.0 ),
      tkTrajUseVertex = cms.bool( False ),
      MuonTrackingRegionBuilder = cms.PSet( 
        Rescale_Dz = cms.double( 4.0 ),
        Pt_fixed = cms.bool( False ),
        Eta_fixed = cms.bool( True ),
        Eta_min = cms.double( 0.1 ),
        DeltaZ = cms.double( 24.2 ),
        maxRegions = cms.int32( 2 ),
        EtaR_UpperLimit_Par1 = cms.double( 0.25 ),
        UseVertex = cms.bool( False ),
        Z_fixed = cms.bool( False ),
        PhiR_UpperLimit_Par1 = cms.double( 0.6 ),
        PhiR_UpperLimit_Par2 = cms.double( 0.2 ),
        Rescale_phi = cms.double( 3.0 ),
        DeltaEta = cms.double( 0.2 ),
        precise = cms.bool( True ),
        OnDemand = cms.int32( -1 ),
        EtaR_UpperLimit_Par2 = cms.double( 0.15 ),
        MeasurementTrackerName = cms.InputTag( "hltESPMeasurementTracker" ),
        vertexCollection = cms.InputTag( "pixelVertices" ),
        Pt_min = cms.double( 3.0 ),
        beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
        Phi_fixed = cms.bool( True ),
        DeltaR = cms.double( 0.025 ),
        input = cms.InputTag( 'hltL2Muons','UpdatedAtVtx' ),
        DeltaPhi = cms.double( 0.15 ),
        Phi_min = cms.double( 0.1 ),
        Rescale_eta = cms.double( 3.0 )
      ),
      TrackTransformer = cms.PSet( 
        Fitter = cms.string( "hltESPL3MuKFTrajectoryFitter" ),
        RefitDirection = cms.string( "insideOut" ),
        RefitRPCHits = cms.bool( True ),
        Propagator = cms.string( "hltESPSmartPropagatorAny" ),
        DoPredictionsOnly = cms.bool( False ),
        TrackerRecHitBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
        MuonRecHitBuilder = cms.string( "hltESPMuonTransientTrackingRecHitBuilder" ),
        Smoother = cms.string( "hltESPKFTrajectorySmootherForMuonTrackLoader" )
      ),
      tkTrajBeamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
      RefitRPCHits = cms.bool( True ),
      tkTrajVertex = cms.InputTag( "Notused" ),
      GlbRefitterParameters = cms.PSet( 
        Fitter = cms.string( "hltESPL3MuKFTrajectoryFitter" ),
        DTRecSegmentLabel = cms.InputTag( "hltDt4DSegments" ),
        RefitFlag = cms.bool( True ),
        SkipStation = cms.int32( -1 ),
        Chi2CutRPC = cms.double( 1.0 ),
        PropDirForCosmics = cms.bool( False ),
        CSCRecSegmentLabel = cms.InputTag( "hltCscSegments" ),
        GEMRecHitLabel = cms.InputTag( "hltGemRecHits" ),
        HitThreshold = cms.int32( 1 ),
        Chi2CutGEM = cms.double( 1.0 ),
        DYTthrs = cms.vint32( 30, 15 ),
        TrackerSkipSystem = cms.int32( -1 ),
        RefitDirection = cms.string( "insideOut" ),
        Chi2CutCSC = cms.double( 150.0 ),
        Chi2CutDT = cms.double( 10.0 ),
        RefitRPCHits = cms.bool( True ),
        TrackerSkipSection = cms.int32( -1 ),
        Propagator = cms.string( "hltESPSmartPropagatorAny" ),
        DoPredictionsOnly = cms.bool( False ),
        TrackerRecHitBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
        MuonHitsOption = cms.int32( 1 ),
        MuonRecHitBuilder = cms.string( "hltESPMuonTransientTrackingRecHitBuilder" )
      ),
      PCut = cms.double( 2.5 ),
      tkTrajMaxDXYBeamSpot = cms.double( 9999.0 ),
      TrackerRecHitBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
      tkTrajMaxChi2 = cms.double( 9999.0 ),
      MuonRecHitBuilder = cms.string( "hltESPMuonTransientTrackingRecHitBuilder" ),
      ScaleTECyFactor = cms.double( -1.0 ),
      tkTrajLabel = cms.InputTag( "hltIterL3MuonAndMuonFromL1Merged" )
    )
)
process.hltIterL3MuonsNoID = cms.EDProducer( "MuonIdProducer",
    MuonCaloCompatibility = cms.PSet( 
      delta_eta = cms.double( 0.02 ),
      delta_phi = cms.double( 0.02 ),
      allSiPMHO = cms.bool( False ),
      MuonTemplateFileName = cms.FileInPath( "RecoMuon/MuonIdentification/data/MuID_templates_muons_lowPt_3_1_norm.root" ),
      PionTemplateFileName = cms.FileInPath( "RecoMuon/MuonIdentification/data/MuID_templates_pions_lowPt_3_1_norm.root" )
    ),
    TrackAssociatorParameters = cms.PSet( 
      useMuon = cms.bool( True ),
      truthMatch = cms.bool( False ),
      usePreshower = cms.bool( False ),
      dRPreshowerPreselection = cms.double( 0.2 ),
      muonMaxDistanceSigmaY = cms.double( 0.0 ),
      useEcal = cms.bool( False ),
      muonMaxDistanceSigmaX = cms.double( 0.0 ),
      dRMuon = cms.double( 9999.0 ),
      dREcal = cms.double( 9999.0 ),
      CSCSegmentCollectionLabel = cms.InputTag( "hltCscSegments" ),
      DTRecSegment4DCollectionLabel = cms.InputTag( "hltDt4DSegments" ),
      EBRecHitCollectionLabel = cms.InputTag( "Notused" ),
      useGEM = cms.bool( True ),
      GEMSegmentCollectionLabel = cms.InputTag( "hltGemSegments" ),
      CaloTowerCollectionLabel = cms.InputTag( "Notused" ),
      propagateAllDirections = cms.bool( True ),
      muonMaxDistanceY = cms.double( 5.0 ),
      useHO = cms.bool( False ),
      muonMaxDistanceX = cms.double( 5.0 ),
      trajectoryUncertaintyTolerance = cms.double( -1.0 ),
      useHcal = cms.bool( False ),
      HBHERecHitCollectionLabel = cms.InputTag( "Notused" ),
      accountForTrajectoryChangeCalo = cms.bool( False ),
      dREcalPreselection = cms.double( 0.05 ),
      useCalo = cms.bool( False ),
      dRMuonPreselection = cms.double( 0.2 ),
      EERecHitCollectionLabel = cms.InputTag( "Notused" ),
      dRHcal = cms.double( 9999.0 ),
      dRHcalPreselection = cms.double( 0.2 ),
      HORecHitCollectionLabel = cms.InputTag( "Notused" )
    ),
    CaloExtractorPSet = cms.PSet( 
      DR_Veto_H = cms.double( 0.1 ),
      CenterConeOnCalIntersection = cms.bool( False ),
      NoiseTow_EE = cms.double( 0.15 ),
      Noise_EB = cms.double( 0.025 ),
      Noise_HE = cms.double( 0.2 ),
      DR_Veto_E = cms.double( 0.07 ),
      NoiseTow_EB = cms.double( 0.04 ),
      Noise_EE = cms.double( 0.1 ),
      UseRecHitsFlag = cms.bool( False ),
      DR_Max = cms.double( 1.0 ),
      DepositLabel = cms.untracked.string( "Cal" ),
      Noise_HO = cms.double( 0.2 ),
      DR_Veto_HO = cms.double( 0.1 ),
      Threshold_H = cms.double( 0.5 ),
      PrintTimeReport = cms.untracked.bool( False ),
      Threshold_E = cms.double( 0.2 ),
      PropagatorName = cms.string( "hltESPFastSteppingHelixPropagatorAny" ),
      ComponentName = cms.string( "CaloExtractorByAssociator" ),
      Threshold_HO = cms.double( 0.5 ),
      DepositInstanceLabels = cms.vstring( 'ecal',
        'hcal',
        'ho' ),
      ServiceParameters = cms.PSet( 
        RPCLayers = cms.bool( False ),
        UseMuonNavigation = cms.untracked.bool( False ),
        Propagators = cms.untracked.vstring( 'hltESPFastSteppingHelixPropagatorAny' )
      ),
      TrackAssociatorParameters = cms.PSet( 
        useMuon = cms.bool( False ),
        truthMatch = cms.bool( False ),
        usePreshower = cms.bool( False ),
        dRPreshowerPreselection = cms.double( 0.2 ),
        muonMaxDistanceSigmaY = cms.double( 0.0 ),
        useEcal = cms.bool( False ),
        muonMaxDistanceSigmaX = cms.double( 0.0 ),
        dRMuon = cms.double( 9999.0 ),
        dREcal = cms.double( 1.0 ),
        CSCSegmentCollectionLabel = cms.InputTag( "hltCscSegments" ),
        DTRecSegment4DCollectionLabel = cms.InputTag( "hltDt4DSegments" ),
        EBRecHitCollectionLabel = cms.InputTag( "Notused" ),
        CaloTowerCollectionLabel = cms.InputTag( "Notused" ),
        propagateAllDirections = cms.bool( True ),
        muonMaxDistanceY = cms.double( 5.0 ),
        useHO = cms.bool( False ),
        muonMaxDistanceX = cms.double( 5.0 ),
        trajectoryUncertaintyTolerance = cms.double( -1.0 ),
        useHcal = cms.bool( False ),
        HBHERecHitCollectionLabel = cms.InputTag( "Notused" ),
        accountForTrajectoryChangeCalo = cms.bool( False ),
        dREcalPreselection = cms.double( 1.0 ),
        useCalo = cms.bool( True ),
        dRMuonPreselection = cms.double( 0.2 ),
        EERecHitCollectionLabel = cms.InputTag( "Notused" ),
        dRHcal = cms.double( 1.0 ),
        dRHcalPreselection = cms.double( 1.0 ),
        HORecHitCollectionLabel = cms.InputTag( "Notused" )
      ),
      Noise_HB = cms.double( 0.2 )
    ),
    TrackExtractorPSet = cms.PSet( 
      Diff_z = cms.double( 0.2 ),
      inputTrackCollection = cms.InputTag( "hltIter0IterL3FromL1MuonTrackSelectionHighPurity" ),
      Chi2Ndof_Max = cms.double( 1.0E64 ),
      BeamSpotLabel = cms.InputTag( "hltOnlineBeamSpot" ),
      DR_Veto = cms.double( 0.01 ),
      Pt_Min = cms.double( -1.0 ),
      DR_Max = cms.double( 1.0 ),
      NHits_Min = cms.uint32( 0 ),
      Chi2Prob_Min = cms.double( -1.0 ),
      Diff_r = cms.double( 0.1 ),
      BeamlineOption = cms.string( "BeamSpotFromEvent" ),
      ComponentName = cms.string( "TrackExtractor" )
    ),
    JetExtractorPSet = cms.PSet( 
      JetCollectionLabel = cms.InputTag( "Notused" ),
      DR_Veto = cms.double( 0.1 ),
      DR_Max = cms.double( 1.0 ),
      ExcludeMuonVeto = cms.bool( True ),
      PrintTimeReport = cms.untracked.bool( False ),
      PropagatorName = cms.string( "hltESPFastSteppingHelixPropagatorAny" ),
      ComponentName = cms.string( "JetExtractor" ),
      ServiceParameters = cms.PSet( 
        RPCLayers = cms.bool( False ),
        UseMuonNavigation = cms.untracked.bool( False ),
        Propagators = cms.untracked.vstring( 'hltESPFastSteppingHelixPropagatorAny' )
      ),
      TrackAssociatorParameters = cms.PSet( 
        useMuon = cms.bool( False ),
        truthMatch = cms.bool( False ),
        usePreshower = cms.bool( False ),
        dRPreshowerPreselection = cms.double( 0.2 ),
        muonMaxDistanceSigmaY = cms.double( 0.0 ),
        useEcal = cms.bool( False ),
        muonMaxDistanceSigmaX = cms.double( 0.0 ),
        dRMuon = cms.double( 9999.0 ),
        dREcal = cms.double( 0.5 ),
        CSCSegmentCollectionLabel = cms.InputTag( "hltCscSegments" ),
        DTRecSegment4DCollectionLabel = cms.InputTag( "hltDt4DSegments" ),
        EBRecHitCollectionLabel = cms.InputTag( "Notused" ),
        CaloTowerCollectionLabel = cms.InputTag( "Notused" ),
        propagateAllDirections = cms.bool( True ),
        muonMaxDistanceY = cms.double( 5.0 ),
        useHO = cms.bool( False ),
        muonMaxDistanceX = cms.double( 5.0 ),
        trajectoryUncertaintyTolerance = cms.double( -1.0 ),
        useHcal = cms.bool( False ),
        HBHERecHitCollectionLabel = cms.InputTag( "Notused" ),
        accountForTrajectoryChangeCalo = cms.bool( False ),
        dREcalPreselection = cms.double( 0.5 ),
        useCalo = cms.bool( True ),
        dRMuonPreselection = cms.double( 0.2 ),
        EERecHitCollectionLabel = cms.InputTag( "Notused" ),
        dRHcal = cms.double( 0.5 ),
        dRHcalPreselection = cms.double( 0.5 ),
        HORecHitCollectionLabel = cms.InputTag( "Notused" )
      ),
      Threshold = cms.double( 5.0 )
    ),
    trackDepositName = cms.string( "tracker" ),
    ecalDepositName = cms.string( "ecal" ),
    hcalDepositName = cms.string( "hcal" ),
    hoDepositName = cms.string( "ho" ),
    jetDepositName = cms.string( "jets" ),
    TimingFillerParameters = cms.PSet( 
      DTTimingParameters = cms.PSet( 
        HitError = cms.double( 6.0 ),
        MatchParameters = cms.PSet( 
          TightMatchDT = cms.bool( False ),
          DTradius = cms.double( 0.01 ),
          TightMatchCSC = cms.bool( True ),
          CSCsegments = cms.InputTag( "hltCscSegments" ),
          DTsegments = cms.InputTag( "hltDt4DSegments" )
        ),
        debug = cms.bool( False ),
        DoWireCorr = cms.bool( False ),
        RequireBothProjections = cms.bool( False ),
        DTTimeOffset = cms.double( 2.7 ),
        PruneCut = cms.double( 10000.0 ),
        DTsegments = cms.InputTag( "hltDt4DSegments" ),
        UseSegmentT0 = cms.bool( False ),
        HitsMin = cms.int32( 5 ),
        DropTheta = cms.bool( True ),
        ServiceParameters = cms.PSet( 
          RPCLayers = cms.bool( True ),
          Propagators = cms.untracked.vstring( 'hltESPFastSteppingHelixPropagatorAny' )
        )
      ),
      UseCSC = cms.bool( True ),
      CSCTimingParameters = cms.PSet( 
        MatchParameters = cms.PSet( 
          TightMatchDT = cms.bool( False ),
          DTradius = cms.double( 0.01 ),
          TightMatchCSC = cms.bool( True ),
          CSCsegments = cms.InputTag( "hltCscSegments" ),
          DTsegments = cms.InputTag( "hltDt4DSegments" )
        ),
        debug = cms.bool( False ),
        CSCWireTimeOffset = cms.double( 0.0 ),
        CSCStripError = cms.double( 7.0 ),
        CSCTimeOffset = cms.double( 0.0 ),
        CSCWireError = cms.double( 8.6 ),
        PruneCut = cms.double( 100.0 ),
        CSCsegments = cms.InputTag( "hltCscSegments" ),
        UseStripTime = cms.bool( True ),
        CSCStripTimeOffset = cms.double( 0.0 ),
        UseWireTime = cms.bool( True ),
        ServiceParameters = cms.PSet( 
          RPCLayers = cms.bool( True ),
          Propagators = cms.untracked.vstring( 'hltESPFastSteppingHelixPropagatorAny' )
        )
      ),
      ErrorDT = cms.double( 6.0 ),
      EcalEnergyCut = cms.double( 0.4 ),
      UseECAL = cms.bool( True ),
      ErrorEB = cms.double( 2.085 ),
      UseDT = cms.bool( True ),
      ErrorEE = cms.double( 6.95 ),
      ErrorCSC = cms.double( 7.4 )
    ),
    ShowerDigiFillerParameters = cms.PSet( 
      cscDigiCollectionLabel = cms.InputTag( 'muonCSCDigis','MuonCSCStripDigi' ),
      digiMaxDistanceX = cms.double( 25.0 ),
      dtDigiCollectionLabel = cms.InputTag( "muonDTDigis" )
    ),
    TrackerKinkFinderParameters = cms.PSet( 
      usePosition = cms.bool( False ),
      diagonalOnly = cms.bool( False )
    ),
    fillEnergy = cms.bool( False ),
    storeCrossedHcalRecHits = cms.bool( False ),
    maxAbsPullX = cms.double( 4.0 ),
    maxAbsEta = cms.double( 3.0 ),
    minPt = cms.double( 2.0 ),
    inputCollectionTypes = cms.vstring( 'inner tracks',
      'links',
      'outer tracks' ),
    addExtraSoftMuons = cms.bool( False ),
    fillGlobalTrackRefits = cms.bool( False ),
    debugWithTruthMatching = cms.bool( False ),
    inputCollectionLabels = cms.VInputTag( 'hltIterL3MuonAndMuonFromL1Merged','hltIterL3GlbMuon','hltL2Muons:UpdatedAtVtx' ),
    fillCaloCompatibility = cms.bool( False ),
    maxAbsPullY = cms.double( 9999.0 ),
    maxAbsDy = cms.double( 9999.0 ),
    minP = cms.double( 0.0 ),
    minPCaloMuon = cms.double( 1.0E9 ),
    maxAbsDx = cms.double( 3.0 ),
    fillIsolation = cms.bool( False ),
    writeIsoDeposits = cms.bool( False ),
    minNumberOfMatches = cms.int32( 1 ),
    fillMatching = cms.bool( True ),
    fillShowerDigis = cms.bool( False ),
    ptThresholdToFillCandidateP4WithGlobalFit = cms.double( 200.0 ),
    sigmaThresholdToFillCandidateP4WithGlobalFit = cms.double( 2.0 ),
    fillGlobalTrackQuality = cms.bool( False ),
    globalTrackQualityInputTag = cms.InputTag( "glbTrackQual" ),
    selectHighPurity = cms.bool( False ),
    pvInputTag = cms.InputTag( "offlinePrimaryVertices" ),
    fillTrackerKink = cms.bool( False ),
    minCaloCompatibility = cms.double( 0.6 ),
    runArbitrationCleaner = cms.bool( False ),
    arbitrationCleanerOptions = cms.PSet( 
      OverlapDTheta = cms.double( 0.02 ),
      Overlap = cms.bool( True ),
      Clustering = cms.bool( True ),
      ME1a = cms.bool( True ),
      ClusterDTheta = cms.double( 0.02 ),
      ClusterDPhi = cms.double( 0.6 ),
      OverlapDPhi = cms.double( 0.0786 )
    ),
    arbitrateTrackerMuons = cms.bool( True )
)
process.hltIterL3Muons = cms.EDProducer( "MuonIDFilterProducerForHLT",
    inputMuonCollection = cms.InputTag( "hltIterL3MuonsNoID" ),
    applyTriggerIdLoose = cms.bool( True ),
    typeMuon = cms.uint32( 0 ),
    allowedTypeMask = cms.uint32( 0 ),
    requiredTypeMask = cms.uint32( 0 ),
    minNMuonHits = cms.int32( 0 ),
    minNMuonStations = cms.int32( 0 ),
    minNTrkLayers = cms.int32( 0 ),
    minTrkHits = cms.int32( 0 ),
    minPixLayer = cms.int32( 0 ),
    minPixHits = cms.int32( 0 ),
    minPt = cms.double( 0.0 ),
    maxNormalizedChi2 = cms.double( 9999.0 )
)
process.hltL3MuonsIterL3Links = cms.EDProducer( "MuonLinksProducer",
    inputCollection = cms.InputTag( "hltIterL3Muons" )
)
process.hltIterL3MuonTracks = cms.EDProducer( "HLTMuonTrackSelector",
    track = cms.InputTag( "hltIterL3MuonAndMuonFromL1Merged" ),
    muon = cms.InputTag( "hltIterL3Muons" ),
    originalMVAVals = cms.InputTag( "none" ),
    copyMVA = cms.bool( False ),
    copyExtras = cms.untracked.bool( True ),
    copyTrajectories = cms.untracked.bool( False )
)
process.hltIterL3MuonCandidates = cms.EDProducer( "L3MuonCandidateProducerFromMuons",
    InputObjects = cms.InputTag( "hltIterL3Muons" )
)
process.hltPixelTracksZetaClean = cms.EDProducer( "TrackWithVertexSelector",
    src = cms.InputTag( "hltPixelTracks" ),
    etaMin = cms.double( 0.0 ),
    etaMax = cms.double( 5.0 ),
    ptMin = cms.double( 0.3 ),
    ptMax = cms.double( 500.0 ),
    d0Max = cms.double( 999.0 ),
    dzMax = cms.double( 999.0 ),
    normalizedChi2 = cms.double( 999999.0 ),
    numberOfValidHits = cms.uint32( 0 ),
    numberOfLostHits = cms.uint32( 999 ),
    numberOfValidPixelHits = cms.uint32( 3 ),
    ptErrorCut = cms.double( 5.0 ),
    quality = cms.string( "highPurity" ),
    useVtx = cms.bool( True ),
    vertexTag = cms.InputTag( "hltPixelVertices" ),
    timesTag = cms.InputTag( "" ),
    timeResosTag = cms.InputTag( "" ),
    nVertices = cms.uint32( 2 ),
    vtxFallback = cms.bool( True ),
    zetaVtx = cms.double( 0.3 ),
    rhoVtx = cms.double( 0.2 ),
    nSigmaDtVertex = cms.double( 0.0 ),
    copyExtras = cms.untracked.bool( False ),
    copyTrajectories = cms.untracked.bool( False )
)
process.hltPixelOnlyPFMuonMerging = cms.EDProducer( "TrackListMerger",
    ShareFrac = cms.double( 0.19 ),
    FoundHitBonus = cms.double( 5.0 ),
    LostHitPenalty = cms.double( 20.0 ),
    MinPT = cms.double( 0.05 ),
    Epsilon = cms.double( -0.001 ),
    MaxNormalizedChisq = cms.double( 1000.0 ),
    MinFound = cms.int32( 3 ),
    TrackProducers = cms.VInputTag( 'hltIterL3MuonTracks','hltPixelTracksZetaClean' ),
    hasSelector = cms.vint32( 0, 0 ),
    indivShareFrac = cms.vdouble( 1.0, 1.0 ),
    selectedTrackQuals = cms.VInputTag( 'hltIterL3MuonTracks','hltPixelTracksZetaClean' ),
    setsToMerge = cms.VPSet( 
      cms.PSet(  pQual = cms.bool( False ),
        tLists = cms.vint32( 0, 1 )
      )
    ),
    trackAlgoPriorityOrder = cms.string( "hltESPTrackAlgoPriorityOrder" ),
    allowFirstHitShare = cms.bool( True ),
    newQuality = cms.string( "confirmed" ),
    copyExtras = cms.untracked.bool( True ),
    writeOnlyTrkQuals = cms.bool( False ),
    copyMVA = cms.bool( False )
)
process.hltPixelOnlyMuonLinks = cms.EDProducer( "MuonLinksProducerForHLT",
    InclusiveTrackerTrackCollection = cms.InputTag( "hltPixelOnlyPFMuonMerging" ),
    LinkCollection = cms.InputTag( "hltL3MuonsIterL3Links" ),
    ptMin = cms.double( 2.5 ),
    pMin = cms.double( 2.5 ),
    shareHitFraction = cms.double( 0.8 )
)
process.hltPixelOnlyMuons = cms.EDProducer( "MuonIdProducer",
    MuonCaloCompatibility = cms.PSet( 
      delta_eta = cms.double( 0.02 ),
      delta_phi = cms.double( 0.02 ),
      allSiPMHO = cms.bool( False ),
      MuonTemplateFileName = cms.FileInPath( "RecoMuon/MuonIdentification/data/MuID_templates_muons_lowPt_3_1_norm.root" ),
      PionTemplateFileName = cms.FileInPath( "RecoMuon/MuonIdentification/data/MuID_templates_pions_lowPt_3_1_norm.root" )
    ),
    TrackAssociatorParameters = cms.PSet( 
      useMuon = cms.bool( True ),
      truthMatch = cms.bool( False ),
      usePreshower = cms.bool( False ),
      dRPreshowerPreselection = cms.double( 0.2 ),
      muonMaxDistanceSigmaY = cms.double( 0.0 ),
      useEcal = cms.bool( True ),
      muonMaxDistanceSigmaX = cms.double( 0.0 ),
      dRMuon = cms.double( 9999.0 ),
      dREcal = cms.double( 9999.0 ),
      CSCSegmentCollectionLabel = cms.InputTag( "hltCscSegments" ),
      DTRecSegment4DCollectionLabel = cms.InputTag( "hltDt4DSegments" ),
      EBRecHitCollectionLabel = cms.InputTag( 'hltEcalRecHit','EcalRecHitsEB' ),
      CaloTowerCollectionLabel = cms.InputTag( "hltTowerMakerForAll" ),
      propagateAllDirections = cms.bool( True ),
      muonMaxDistanceY = cms.double( 5.0 ),
      useHO = cms.bool( True ),
      muonMaxDistanceX = cms.double( 5.0 ),
      trajectoryUncertaintyTolerance = cms.double( -1.0 ),
      useHcal = cms.bool( True ),
      HBHERecHitCollectionLabel = cms.InputTag( "hltHbhereco" ),
      accountForTrajectoryChangeCalo = cms.bool( False ),
      dREcalPreselection = cms.double( 0.05 ),
      useCalo = cms.bool( False ),
      dRMuonPreselection = cms.double( 0.2 ),
      EERecHitCollectionLabel = cms.InputTag( 'hltEcalRecHit','EcalRecHitsEE' ),
      dRHcal = cms.double( 9999.0 ),
      dRHcalPreselection = cms.double( 0.2 ),
      HORecHitCollectionLabel = cms.InputTag( "hltHoreco" )
    ),
    CaloExtractorPSet = cms.PSet( 
      DR_Veto_H = cms.double( 0.1 ),
      CenterConeOnCalIntersection = cms.bool( False ),
      NoiseTow_EE = cms.double( 0.15 ),
      Noise_EB = cms.double( 0.025 ),
      Noise_HE = cms.double( 0.2 ),
      DR_Veto_E = cms.double( 0.07 ),
      NoiseTow_EB = cms.double( 0.04 ),
      Noise_EE = cms.double( 0.1 ),
      UseRecHitsFlag = cms.bool( False ),
      DR_Max = cms.double( 1.0 ),
      DepositLabel = cms.untracked.string( "Cal" ),
      Noise_HO = cms.double( 0.2 ),
      DR_Veto_HO = cms.double( 0.1 ),
      Threshold_H = cms.double( 0.5 ),
      PrintTimeReport = cms.untracked.bool( False ),
      Threshold_E = cms.double( 0.2 ),
      PropagatorName = cms.string( "hltESPFastSteppingHelixPropagatorAny" ),
      ComponentName = cms.string( "CaloExtractorByAssociator" ),
      Threshold_HO = cms.double( 0.5 ),
      DepositInstanceLabels = cms.vstring( 'ecal',
        'hcal',
        'ho' ),
      ServiceParameters = cms.PSet( 
        RPCLayers = cms.bool( False ),
        UseMuonNavigation = cms.untracked.bool( False ),
        Propagators = cms.untracked.vstring( 'hltESPFastSteppingHelixPropagatorAny' )
      ),
      TrackAssociatorParameters = cms.PSet( 
        useMuon = cms.bool( False ),
        truthMatch = cms.bool( False ),
        usePreshower = cms.bool( False ),
        dRPreshowerPreselection = cms.double( 0.2 ),
        muonMaxDistanceSigmaY = cms.double( 0.0 ),
        useEcal = cms.bool( False ),
        muonMaxDistanceSigmaX = cms.double( 0.0 ),
        dRMuon = cms.double( 9999.0 ),
        dREcal = cms.double( 1.0 ),
        CSCSegmentCollectionLabel = cms.InputTag( "hltCscSegments" ),
        DTRecSegment4DCollectionLabel = cms.InputTag( "hltDt4DSegments" ),
        EBRecHitCollectionLabel = cms.InputTag( 'hltEcalRecHit','EcalRecHitsEB' ),
        CaloTowerCollectionLabel = cms.InputTag( "hltTowerMakerForAll" ),
        propagateAllDirections = cms.bool( True ),
        muonMaxDistanceY = cms.double( 5.0 ),
        useHO = cms.bool( False ),
        muonMaxDistanceX = cms.double( 5.0 ),
        trajectoryUncertaintyTolerance = cms.double( -1.0 ),
        useHcal = cms.bool( False ),
        HBHERecHitCollectionLabel = cms.InputTag( "hltHbhereco" ),
        accountForTrajectoryChangeCalo = cms.bool( False ),
        dREcalPreselection = cms.double( 1.0 ),
        useCalo = cms.bool( True ),
        dRMuonPreselection = cms.double( 0.2 ),
        EERecHitCollectionLabel = cms.InputTag( 'hltEcalRecHit','EcalRecHitsEE' ),
        dRHcal = cms.double( 1.0 ),
        dRHcalPreselection = cms.double( 1.0 ),
        HORecHitCollectionLabel = cms.InputTag( "hltHoreco" )
      ),
      Noise_HB = cms.double( 0.2 )
    ),
    TrackExtractorPSet = cms.PSet( 
      Diff_z = cms.double( 0.2 ),
      inputTrackCollection = cms.InputTag( "hltPixelOnlyPFMuonMerging" ),
      Chi2Ndof_Max = cms.double( 1.0E64 ),
      BeamSpotLabel = cms.InputTag( "hltOnlineBeamSpot" ),
      DR_Veto = cms.double( 0.01 ),
      Pt_Min = cms.double( -1.0 ),
      DR_Max = cms.double( 1.0 ),
      DepositLabel = cms.untracked.string( "" ),
      NHits_Min = cms.uint32( 0 ),
      Chi2Prob_Min = cms.double( -1.0 ),
      Diff_r = cms.double( 0.1 ),
      BeamlineOption = cms.string( "BeamSpotFromEvent" ),
      ComponentName = cms.string( "TrackExtractor" )
    ),
    JetExtractorPSet = cms.PSet( 
      JetCollectionLabel = cms.InputTag( "hltAK4CaloJetsPFEt5" ),
      DR_Veto = cms.double( 0.1 ),
      DR_Max = cms.double( 1.0 ),
      ExcludeMuonVeto = cms.bool( True ),
      PrintTimeReport = cms.untracked.bool( False ),
      PropagatorName = cms.string( "hltESPFastSteppingHelixPropagatorAny" ),
      ComponentName = cms.string( "JetExtractor" ),
      ServiceParameters = cms.PSet( 
        RPCLayers = cms.bool( False ),
        UseMuonNavigation = cms.untracked.bool( False ),
        Propagators = cms.untracked.vstring( 'hltESPFastSteppingHelixPropagatorAny' )
      ),
      TrackAssociatorParameters = cms.PSet( 
        useMuon = cms.bool( False ),
        truthMatch = cms.bool( False ),
        usePreshower = cms.bool( False ),
        dRPreshowerPreselection = cms.double( 0.2 ),
        muonMaxDistanceSigmaY = cms.double( 0.0 ),
        useEcal = cms.bool( False ),
        muonMaxDistanceSigmaX = cms.double( 0.0 ),
        dRMuon = cms.double( 9999.0 ),
        dREcal = cms.double( 0.5 ),
        CSCSegmentCollectionLabel = cms.InputTag( "hltCscSegments" ),
        DTRecSegment4DCollectionLabel = cms.InputTag( "hltDt4DSegments" ),
        EBRecHitCollectionLabel = cms.InputTag( 'hltEcalRecHit','EcalRecHitsEB' ),
        CaloTowerCollectionLabel = cms.InputTag( "hltTowerMakerForAll" ),
        propagateAllDirections = cms.bool( True ),
        muonMaxDistanceY = cms.double( 5.0 ),
        useHO = cms.bool( False ),
        muonMaxDistanceX = cms.double( 5.0 ),
        trajectoryUncertaintyTolerance = cms.double( -1.0 ),
        useHcal = cms.bool( False ),
        HBHERecHitCollectionLabel = cms.InputTag( "hltHbhereco" ),
        accountForTrajectoryChangeCalo = cms.bool( False ),
        dREcalPreselection = cms.double( 0.5 ),
        useCalo = cms.bool( True ),
        dRMuonPreselection = cms.double( 0.2 ),
        EERecHitCollectionLabel = cms.InputTag( 'hltEcalRecHit','EcalRecHitsEE' ),
        dRHcal = cms.double( 0.5 ),
        dRHcalPreselection = cms.double( 0.5 ),
        HORecHitCollectionLabel = cms.InputTag( "hltHoreco" )
      ),
      Threshold = cms.double( 5.0 )
    ),
    trackDepositName = cms.string( "tracker" ),
    ecalDepositName = cms.string( "ecal" ),
    hcalDepositName = cms.string( "hcal" ),
    hoDepositName = cms.string( "ho" ),
    jetDepositName = cms.string( "jets" ),
    TimingFillerParameters = cms.PSet( 
      DTTimingParameters = cms.PSet( 
        HitError = cms.double( 6.0 ),
        MatchParameters = cms.PSet( 
          TightMatchDT = cms.bool( False ),
          DTradius = cms.double( 0.01 ),
          TightMatchCSC = cms.bool( True ),
          CSCsegments = cms.InputTag( "hltCscSegments" ),
          DTsegments = cms.InputTag( "hltDt4DSegments" )
        ),
        debug = cms.bool( False ),
        DoWireCorr = cms.bool( False ),
        RequireBothProjections = cms.bool( False ),
        DTTimeOffset = cms.double( 2.7 ),
        PruneCut = cms.double( 10000.0 ),
        DTsegments = cms.InputTag( "hltDt4DSegments" ),
        UseSegmentT0 = cms.bool( False ),
        HitsMin = cms.int32( 5 ),
        DropTheta = cms.bool( True ),
        ServiceParameters = cms.PSet( 
          RPCLayers = cms.bool( True ),
          Propagators = cms.untracked.vstring( 'hltESPFastSteppingHelixPropagatorAny' )
        )
      ),
      UseCSC = cms.bool( True ),
      CSCTimingParameters = cms.PSet( 
        MatchParameters = cms.PSet( 
          TightMatchDT = cms.bool( False ),
          DTradius = cms.double( 0.01 ),
          TightMatchCSC = cms.bool( True ),
          CSCsegments = cms.InputTag( "hltCscSegments" ),
          DTsegments = cms.InputTag( "hltDt4DSegments" )
        ),
        debug = cms.bool( False ),
        CSCWireTimeOffset = cms.double( 0.0 ),
        CSCStripError = cms.double( 7.0 ),
        CSCTimeOffset = cms.double( 0.0 ),
        CSCWireError = cms.double( 8.6 ),
        PruneCut = cms.double( 100.0 ),
        CSCsegments = cms.InputTag( "hltCscSegments" ),
        UseStripTime = cms.bool( True ),
        CSCStripTimeOffset = cms.double( 0.0 ),
        UseWireTime = cms.bool( True ),
        ServiceParameters = cms.PSet( 
          RPCLayers = cms.bool( True ),
          Propagators = cms.untracked.vstring( 'hltESPFastSteppingHelixPropagatorAny' )
        )
      ),
      ErrorDT = cms.double( 6.0 ),
      EcalEnergyCut = cms.double( 0.4 ),
      UseECAL = cms.bool( True ),
      ErrorEB = cms.double( 2.085 ),
      UseDT = cms.bool( True ),
      ErrorEE = cms.double( 6.95 ),
      ErrorCSC = cms.double( 7.4 )
    ),
    ShowerDigiFillerParameters = cms.PSet( 
      cscDigiCollectionLabel = cms.InputTag( 'muonCSCDigis','MuonCSCStripDigi' ),
      dtDigiCollectionLabel = cms.InputTag( "muonDTDigis" ),
      digiMaxDistanceX = cms.double( 25.0 )
    ),
    TrackerKinkFinderParameters = cms.PSet( 
      usePosition = cms.bool( False ),
      diagonalOnly = cms.bool( False )
    ),
    fillEnergy = cms.bool( True ),
    storeCrossedHcalRecHits = cms.bool( False ),
    maxAbsPullX = cms.double( 4.0 ),
    maxAbsEta = cms.double( 3.0 ),
    minPt = cms.double( 10.0 ),
    inputCollectionTypes = cms.vstring( 'inner tracks',
      'links',
      'outer tracks' ),
    addExtraSoftMuons = cms.bool( False ),
    fillGlobalTrackRefits = cms.bool( False ),
    debugWithTruthMatching = cms.bool( False ),
    inputCollectionLabels = cms.VInputTag( 'hltPixelOnlyPFMuonMerging','hltPixelOnlyMuonLinks','hltL2Muons' ),
    fillCaloCompatibility = cms.bool( True ),
    maxAbsPullY = cms.double( 9999.0 ),
    maxAbsDy = cms.double( 9999.0 ),
    minP = cms.double( 10.0 ),
    minPCaloMuon = cms.double( 1.0E9 ),
    maxAbsDx = cms.double( 3.0 ),
    fillIsolation = cms.bool( True ),
    writeIsoDeposits = cms.bool( False ),
    minNumberOfMatches = cms.int32( 1 ),
    fillMatching = cms.bool( True ),
    fillShowerDigis = cms.bool( False ),
    ptThresholdToFillCandidateP4WithGlobalFit = cms.double( 200.0 ),
    sigmaThresholdToFillCandidateP4WithGlobalFit = cms.double( 2.0 ),
    fillGlobalTrackQuality = cms.bool( False ),
    globalTrackQualityInputTag = cms.InputTag( "glbTrackQual" ),
    selectHighPurity = cms.bool( False ),
    pvInputTag = cms.InputTag( "offlinePrimaryVertices" ),
    fillTrackerKink = cms.bool( False ),
    minCaloCompatibility = cms.double( 0.6 ),
    runArbitrationCleaner = cms.bool( False ),
    arbitrationCleanerOptions = cms.PSet( 
      OverlapDTheta = cms.double( 0.02 ),
      Overlap = cms.bool( True ),
      Clustering = cms.bool( True ),
      ME1a = cms.bool( True ),
      ClusterDTheta = cms.double( 0.02 ),
      ClusterDPhi = cms.double( 0.6 ),
      OverlapDPhi = cms.double( 0.0786 )
    ),
    arbitrateTrackerMuons = cms.bool( False )
)
process.hltEcalPreshowerDigis = cms.EDProducer( "ESRawToDigi",
    sourceTag = cms.InputTag( "rawDataCollector" ),
    debugMode = cms.untracked.bool( False ),
    InstanceES = cms.string( "" ),
    LookupTable = cms.FileInPath( "EventFilter/ESDigiToRaw/data/ES_lookup_table.dat" ),
    ESdigiCollection = cms.string( "" )
)
process.hltEcalPreshowerRecHit = cms.EDProducer( "ESRecHitProducer",
    ESrechitCollection = cms.string( "EcalRecHitsES" ),
    ESdigiCollection = cms.InputTag( "hltEcalPreshowerDigis" ),
    algo = cms.string( "ESRecHitWorker" ),
    ESRecoAlgo = cms.int32( 0 )
)
process.hltParticleFlowRecHitECALUnseeded = cms.EDProducer( "PFRecHitProducer",
    navigator = cms.PSet( 
      barrel = cms.PSet(  ),
      endcap = cms.PSet(  ),
      name = cms.string( "PFRecHitECALNavigator" )
    ),
    producers = cms.VPSet( 
      cms.PSet(  src = cms.InputTag( 'hltEcalRecHit','EcalRecHitsEB' ),
        srFlags = cms.InputTag( "" ),
        name = cms.string( "PFEBRecHitCreator" ),
        qualityTests = cms.VPSet( 
          cms.PSet(  name = cms.string( "PFRecHitQTestDBThreshold" ),
            applySelectionsToAllCrystals = cms.bool( True )
          ),
          cms.PSet(  topologicalCleaning = cms.bool( True ),
            skipTTRecoveredHits = cms.bool( True ),
            cleaningThreshold = cms.double( 2.0 ),
            name = cms.string( "PFRecHitQTestECAL" ),
            timingCleaning = cms.bool( True )
          )
        )
      ),
      cms.PSet(  src = cms.InputTag( 'hltEcalRecHit','EcalRecHitsEE' ),
        srFlags = cms.InputTag( "" ),
        name = cms.string( "PFEERecHitCreator" ),
        qualityTests = cms.VPSet( 
          cms.PSet(  name = cms.string( "PFRecHitQTestDBThreshold" ),
            applySelectionsToAllCrystals = cms.bool( True )
          ),
          cms.PSet(  topologicalCleaning = cms.bool( True ),
            skipTTRecoveredHits = cms.bool( True ),
            cleaningThreshold = cms.double( 2.0 ),
            name = cms.string( "PFRecHitQTestECAL" ),
            timingCleaning = cms.bool( True )
          )
        )
      )
    )
)
process.hltParticleFlowRecHitHBHE = cms.EDProducer( "PFRecHitProducer",
    navigator = cms.PSet( 
      name = cms.string( "PFRecHitHCALDenseIdNavigator" ),
      hcalEnums = cms.vint32( 1, 2 )
    ),
    producers = cms.VPSet( 
      cms.PSet(  src = cms.InputTag( "hltHbhereco" ),
        name = cms.string( "PFHBHERecHitCreator" ),
        qualityTests = cms.VPSet( 
          cms.PSet(  threshold = cms.double( 0.8 ),
            name = cms.string( "PFRecHitQTestThreshold" ),
            cuts = cms.VPSet( 
              cms.PSet(  depth = cms.vint32( 1, 2, 3, 4 ),
                threshold = cms.vdouble( 0.1, 0.2, 0.3, 0.3 ),
                detectorEnum = cms.int32( 1 )
              ),
              cms.PSet(  depth = cms.vint32( 1, 2, 3, 4, 5, 6, 7 ),
                threshold = cms.vdouble( 0.1, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2 ),
                detectorEnum = cms.int32( 2 )
              )
            )
          ),
          cms.PSet(  flags = cms.vstring( 'Standard' ),
            cleaningThresholds = cms.vdouble( 0.0 ),
            name = cms.string( "PFRecHitQTestHCALChannel" ),
            maxSeverities = cms.vint32( 11 )
          )
        )
      )
    )
)
process.hltParticleFlowRecHitHF = cms.EDProducer( "PFRecHitProducer",
    navigator = cms.PSet( 
      name = cms.string( "PFRecHitHCALDenseIdNavigator" ),
      hcalEnums = cms.vint32( 4 )
    ),
    producers = cms.VPSet( 
      cms.PSet(  thresh_HF = cms.double( 0.4 ),
        LongFibre_Fraction = cms.double( 0.1 ),
        src = cms.InputTag( "hltHfreco" ),
        EMDepthCorrection = cms.double( 22.0 ),
        ShortFibre_Fraction = cms.double( 0.01 ),
        HADDepthCorrection = cms.double( 25.0 ),
        HFCalib29 = cms.double( 1.07 ),
        LongFibre_Cut = cms.double( 120.0 ),
        name = cms.string( "PFHFRecHitCreator" ),
        qualityTests = cms.VPSet( 
          cms.PSet(  flags = cms.vstring( 'Standard',
  'HFLong',
  'HFShort' ),
            cleaningThresholds = cms.vdouble( 0.0, 120.0, 60.0 ),
            name = cms.string( "PFRecHitQTestHCALChannel" ),
            maxSeverities = cms.vint32( 11, 9, 9 )
          ),
          cms.PSet(  name = cms.string( "PFRecHitQTestHCALThresholdVsDepth" ),
            cuts = cms.VPSet( 
              cms.PSet(  depth = cms.vint32( 1, 2 ),
                threshold = cms.vdouble( 1.2, 1.8 ),
                detectorEnum = cms.int32( 4 )
              )
            )
          )
        ),
        ShortFibre_Cut = cms.double( 60.0 )
      )
    )
)
process.hltParticleFlowRecHitPSUnseeded = cms.EDProducer( "PFRecHitProducer",
    navigator = cms.PSet(  name = cms.string( "PFRecHitPreshowerNavigator" ) ),
    producers = cms.VPSet( 
      cms.PSet(  src = cms.InputTag( 'hltEcalPreshowerRecHit','EcalRecHitsES' ),
        name = cms.string( "PFPSRecHitCreator" ),
        qualityTests = cms.VPSet( 
          cms.PSet(  threshold = cms.double( 7.0E-6 ),
            name = cms.string( "PFRecHitQTestThreshold" )
          )
        )
      )
    )
)
process.hltParticleFlowClusterECALUncorrectedUnseeded = cms.EDProducer( "PFClusterProducer",
    recHitsSource = cms.InputTag( "hltParticleFlowRecHitECALUnseeded" ),
    recHitCleaners = cms.VPSet( 
    ),
    seedCleaners = cms.VPSet( 
    ),
    seedFinder = cms.PSet( 
      thresholdsByDetector = cms.VPSet( 
        cms.PSet(  seedingThresholdPt = cms.double( 0.15 ),
          seedingThreshold = cms.double( 0.6 ),
          detector = cms.string( "ECAL_ENDCAP" )
        ),
        cms.PSet(  seedingThresholdPt = cms.double( 0.0 ),
          seedingThreshold = cms.double( 0.23 ),
          detector = cms.string( "ECAL_BARREL" )
        )
      ),
      algoName = cms.string( "LocalMaximumSeedFinder" ),
      nNeighbours = cms.int32( 8 )
    ),
    initialClusteringStep = cms.PSet( 
      thresholdsByDetector = cms.VPSet( 
        cms.PSet(  gatheringThreshold = cms.double( 0.08 ),
          gatheringThresholdPt = cms.double( 0.0 ),
          detector = cms.string( "ECAL_BARREL" )
        ),
        cms.PSet(  gatheringThreshold = cms.double( 0.3 ),
          gatheringThresholdPt = cms.double( 0.0 ),
          detector = cms.string( "ECAL_ENDCAP" )
        )
      ),
      algoName = cms.string( "Basic2DGenericTopoClusterizer" ),
      useCornerCells = cms.bool( True )
    ),
    pfClusterBuilder = cms.PSet( 
      minFracTot = cms.double( 1.0E-20 ),
      stoppingTolerance = cms.double( 1.0E-8 ),
      positionCalc = cms.PSet( 
        minAllowedNormalization = cms.double( 1.0E-9 ),
        posCalcNCrystals = cms.int32( 9 ),
        algoName = cms.string( "Basic2DGenericPFlowPositionCalc" ),
        logWeightDenominator = cms.double( 0.08 ),
        minFractionInCalc = cms.double( 1.0E-9 ),
        timeResolutionCalcBarrel = cms.PSet( 
          corrTermLowE = cms.double( 0.0510871 ),
          threshLowE = cms.double( 0.5 ),
          noiseTerm = cms.double( 1.10889 ),
          constantTermLowE = cms.double( 0.0 ),
          noiseTermLowE = cms.double( 1.31883 ),
          threshHighE = cms.double( 5.0 ),
          constantTerm = cms.double( 0.428192 )
        ),
        timeResolutionCalcEndcap = cms.PSet( 
          corrTermLowE = cms.double( 0.0 ),
          threshLowE = cms.double( 1.0 ),
          noiseTerm = cms.double( 5.72489999999 ),
          constantTermLowE = cms.double( 0.0 ),
          noiseTermLowE = cms.double( 6.92683000001 ),
          threshHighE = cms.double( 10.0 ),
          constantTerm = cms.double( 0.0 )
        )
      ),
      maxIterations = cms.uint32( 50 ),
      positionCalcForConvergence = cms.PSet( 
        minAllowedNormalization = cms.double( 0.0 ),
        T0_ES = cms.double( 1.2 ),
        algoName = cms.string( "ECAL2DPositionCalcWithDepthCorr" ),
        T0_EE = cms.double( 3.1 ),
        T0_EB = cms.double( 7.4 ),
        X0 = cms.double( 0.89 ),
        minFractionInCalc = cms.double( 0.0 ),
        W0 = cms.double( 4.2 )
      ),
      allCellsPositionCalc = cms.PSet( 
        minAllowedNormalization = cms.double( 1.0E-9 ),
        posCalcNCrystals = cms.int32( -1 ),
        algoName = cms.string( "Basic2DGenericPFlowPositionCalc" ),
        logWeightDenominator = cms.double( 0.08 ),
        minFractionInCalc = cms.double( 1.0E-9 ),
        timeResolutionCalcBarrel = cms.PSet( 
          corrTermLowE = cms.double( 0.0510871 ),
          threshLowE = cms.double( 0.5 ),
          noiseTerm = cms.double( 1.10889 ),
          constantTermLowE = cms.double( 0.0 ),
          noiseTermLowE = cms.double( 1.31883 ),
          threshHighE = cms.double( 5.0 ),
          constantTerm = cms.double( 0.428192 )
        ),
        timeResolutionCalcEndcap = cms.PSet( 
          corrTermLowE = cms.double( 0.0 ),
          threshLowE = cms.double( 1.0 ),
          noiseTerm = cms.double( 5.72489999999 ),
          constantTermLowE = cms.double( 0.0 ),
          noiseTermLowE = cms.double( 6.92683000001 ),
          threshHighE = cms.double( 10.0 ),
          constantTerm = cms.double( 0.0 )
        )
      ),
      algoName = cms.string( "Basic2DGenericPFlowClusterizer" ),
      recHitEnergyNorms = cms.VPSet( 
        cms.PSet(  recHitEnergyNorm = cms.double( 0.08 ),
          detector = cms.string( "ECAL_BARREL" )
        ),
        cms.PSet(  recHitEnergyNorm = cms.double( 0.3 ),
          detector = cms.string( "ECAL_ENDCAP" )
        )
      ),
      showerSigma = cms.double( 1.5 ),
      minFractionToKeep = cms.double( 1.0E-7 ),
      excludeOtherSeeds = cms.bool( True )
    ),
    positionReCalc = cms.PSet( 
      minAllowedNormalization = cms.double( 0.0 ),
      T0_ES = cms.double( 1.2 ),
      algoName = cms.string( "ECAL2DPositionCalcWithDepthCorr" ),
      T0_EE = cms.double( 3.1 ),
      T0_EB = cms.double( 7.4 ),
      X0 = cms.double( 0.89 ),
      minFractionInCalc = cms.double( 0.0 ),
      W0 = cms.double( 4.2 )
    ),
    energyCorrector = cms.PSet(  )
)
process.hltParticleFlowClusterPSUnseeded = cms.EDProducer( "PFClusterProducer",
    recHitsSource = cms.InputTag( "hltParticleFlowRecHitPSUnseeded" ),
    recHitCleaners = cms.VPSet( 
    ),
    seedCleaners = cms.VPSet( 
    ),
    seedFinder = cms.PSet( 
      thresholdsByDetector = cms.VPSet( 
        cms.PSet(  seedingThresholdPt = cms.double( 0.0 ),
          seedingThreshold = cms.double( 1.2E-4 ),
          detector = cms.string( "PS1" )
        ),
        cms.PSet(  seedingThresholdPt = cms.double( 0.0 ),
          seedingThreshold = cms.double( 1.2E-4 ),
          detector = cms.string( "PS2" )
        )
      ),
      algoName = cms.string( "LocalMaximumSeedFinder" ),
      nNeighbours = cms.int32( 4 )
    ),
    initialClusteringStep = cms.PSet( 
      thresholdsByDetector = cms.VPSet( 
        cms.PSet(  gatheringThreshold = cms.double( 6.0E-5 ),
          gatheringThresholdPt = cms.double( 0.0 ),
          detector = cms.string( "PS1" )
        ),
        cms.PSet(  gatheringThreshold = cms.double( 6.0E-5 ),
          gatheringThresholdPt = cms.double( 0.0 ),
          detector = cms.string( "PS2" )
        )
      ),
      algoName = cms.string( "Basic2DGenericTopoClusterizer" ),
      useCornerCells = cms.bool( False )
    ),
    pfClusterBuilder = cms.PSet( 
      minFracTot = cms.double( 1.0E-20 ),
      stoppingTolerance = cms.double( 1.0E-8 ),
      positionCalc = cms.PSet( 
        minAllowedNormalization = cms.double( 1.0E-9 ),
        posCalcNCrystals = cms.int32( -1 ),
        algoName = cms.string( "Basic2DGenericPFlowPositionCalc" ),
        logWeightDenominator = cms.double( 6.0E-5 ),
        minFractionInCalc = cms.double( 1.0E-9 )
      ),
      maxIterations = cms.uint32( 50 ),
      algoName = cms.string( "Basic2DGenericPFlowClusterizer" ),
      recHitEnergyNorms = cms.VPSet( 
        cms.PSet(  recHitEnergyNorm = cms.double( 6.0E-5 ),
          detector = cms.string( "PS1" )
        ),
        cms.PSet(  recHitEnergyNorm = cms.double( 6.0E-5 ),
          detector = cms.string( "PS2" )
        )
      ),
      showerSigma = cms.double( 0.3 ),
      minFractionToKeep = cms.double( 1.0E-7 ),
      excludeOtherSeeds = cms.bool( True )
    ),
    positionReCalc = cms.PSet(  ),
    energyCorrector = cms.PSet(  )
)
process.hltParticleFlowClusterECALUnseeded = cms.EDProducer( "CorrectedECALPFClusterProducer",
    minimumPSEnergy = cms.double( 0.0 ),
    skipPS = cms.bool( False ),
    inputPS = cms.InputTag( "hltParticleFlowClusterPSUnseeded" ),
    energyCorrector = cms.PSet( 
      applyCrackCorrections = cms.bool( False ),
      srfAwareCorrection = cms.bool( True ),
      applyMVACorrections = cms.bool( True ),
      maxPtForMVAEvaluation = cms.double( 300.0 ),
      recHitsEBLabel = cms.InputTag( 'hltEcalRecHit','EcalRecHitsEB' ),
      recHitsEELabel = cms.InputTag( 'hltEcalRecHit','EcalRecHitsEE' ),
      ebSrFlagLabel = cms.InputTag( "hltEcalDigis" ),
      eeSrFlagLabel = cms.InputTag( "hltEcalDigis" )
    ),
    inputECAL = cms.InputTag( "hltParticleFlowClusterECALUncorrectedUnseeded" )
)
process.hltParticleFlowClusterHBHE = cms.EDProducer( "PFClusterProducer",
    recHitsSource = cms.InputTag( "hltParticleFlowRecHitHBHE" ),
    recHitCleaners = cms.VPSet( 
    ),
    seedCleaners = cms.VPSet( 
    ),
    seedFinder = cms.PSet( 
      thresholdsByDetector = cms.VPSet( 
        cms.PSet(  detector = cms.string( "HCAL_BARREL1" ),
          depths = cms.vint32( 1, 2, 3, 4 ),
          seedingThreshold = cms.vdouble( 0.125, 0.25, 0.35, 0.35 ),
          seedingThresholdPt = cms.vdouble( 0.0, 0.0, 0.0, 0.0 )
        ),
        cms.PSet(  detector = cms.string( "HCAL_ENDCAP" ),
          depths = cms.vint32( 1, 2, 3, 4, 5, 6, 7 ),
          seedingThreshold = cms.vdouble( 0.1375, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275 ),
          seedingThresholdPt = cms.vdouble( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 )
        )
      ),
      algoName = cms.string( "LocalMaximumSeedFinder" ),
      nNeighbours = cms.int32( 4 )
    ),
    initialClusteringStep = cms.PSet( 
      thresholdsByDetector = cms.VPSet( 
        cms.PSet(  detector = cms.string( "HCAL_BARREL1" ),
          depths = cms.vint32( 1, 2, 3, 4 ),
          gatheringThreshold = cms.vdouble( 0.1, 0.2, 0.3, 0.3 ),
          gatheringThresholdPt = cms.vdouble( 0.0, 0.0, 0.0, 0.0 )
        ),
        cms.PSet(  detector = cms.string( "HCAL_ENDCAP" ),
          depths = cms.vint32( 1, 2, 3, 4, 5, 6, 7 ),
          gatheringThreshold = cms.vdouble( 0.1, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2 ),
          gatheringThresholdPt = cms.vdouble( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 )
        )
      ),
      algoName = cms.string( "Basic2DGenericTopoClusterizer" ),
      useCornerCells = cms.bool( True )
    ),
    pfClusterBuilder = cms.PSet( 
      minFracTot = cms.double( 1.0E-20 ),
      stoppingTolerance = cms.double( 1.0E-8 ),
      positionCalc = cms.PSet( 
        minAllowedNormalization = cms.double( 1.0E-9 ),
        posCalcNCrystals = cms.int32( 5 ),
        algoName = cms.string( "Basic2DGenericPFlowPositionCalc" ),
        minFractionInCalc = cms.double( 1.0E-9 ),
        logWeightDenominatorByDetector = cms.VPSet( 
          cms.PSet(  depths = cms.vint32( 1, 2, 3, 4 ),
            detector = cms.string( "HCAL_BARREL1" ),
            logWeightDenominator = cms.vdouble( 0.1, 0.2, 0.3, 0.3 )
          ),
          cms.PSet(  depths = cms.vint32( 1, 2, 3, 4, 5, 6, 7 ),
            detector = cms.string( "HCAL_ENDCAP" ),
            logWeightDenominator = cms.vdouble( 0.1, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2 )
          )
        )
      ),
      maxIterations = cms.uint32( 50 ),
      minChi2Prob = cms.double( 0.0 ),
      allCellsPositionCalc = cms.PSet( 
        minAllowedNormalization = cms.double( 1.0E-9 ),
        posCalcNCrystals = cms.int32( -1 ),
        algoName = cms.string( "Basic2DGenericPFlowPositionCalc" ),
        minFractionInCalc = cms.double( 1.0E-9 ),
        logWeightDenominatorByDetector = cms.VPSet( 
          cms.PSet(  depths = cms.vint32( 1, 2, 3, 4 ),
            detector = cms.string( "HCAL_BARREL1" ),
            logWeightDenominator = cms.vdouble( 0.1, 0.2, 0.3, 0.3 )
          ),
          cms.PSet(  depths = cms.vint32( 1, 2, 3, 4, 5, 6, 7 ),
            detector = cms.string( "HCAL_ENDCAP" ),
            logWeightDenominator = cms.vdouble( 0.1, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2 )
          )
        )
      ),
      algoName = cms.string( "Basic2DGenericPFlowClusterizer" ),
      recHitEnergyNorms = cms.VPSet( 
        cms.PSet(  detector = cms.string( "HCAL_BARREL1" ),
          depths = cms.vint32( 1, 2, 3, 4 ),
          recHitEnergyNorm = cms.vdouble( 0.1, 0.2, 0.3, 0.3 )
        ),
        cms.PSet(  detector = cms.string( "HCAL_ENDCAP" ),
          depths = cms.vint32( 1, 2, 3, 4, 5, 6, 7 ),
          recHitEnergyNorm = cms.vdouble( 0.1, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2 )
        )
      ),
      maxNSigmaTime = cms.double( 10.0 ),
      showerSigma = cms.double( 10.0 ),
      timeSigmaEE = cms.double( 10.0 ),
      clusterTimeResFromSeed = cms.bool( False ),
      minFractionToKeep = cms.double( 1.0E-7 ),
      excludeOtherSeeds = cms.bool( True ),
      timeResolutionCalcBarrel = cms.PSet( 
        corrTermLowE = cms.double( 0.0 ),
        threshLowE = cms.double( 6.0 ),
        noiseTerm = cms.double( 21.86 ),
        constantTermLowE = cms.double( 4.24 ),
        noiseTermLowE = cms.double( 8.0 ),
        threshHighE = cms.double( 15.0 ),
        constantTerm = cms.double( 2.82 )
      ),
      timeResolutionCalcEndcap = cms.PSet( 
        corrTermLowE = cms.double( 0.0 ),
        threshLowE = cms.double( 6.0 ),
        noiseTerm = cms.double( 21.86 ),
        constantTermLowE = cms.double( 4.24 ),
        noiseTermLowE = cms.double( 8.0 ),
        threshHighE = cms.double( 15.0 ),
        constantTerm = cms.double( 2.82 )
      ),
      timeSigmaEB = cms.double( 10.0 )
    ),
    positionReCalc = cms.PSet(  ),
    energyCorrector = cms.PSet(  )
)
process.hltParticleFlowClusterHCAL = cms.EDProducer( "PFMultiDepthClusterProducer",
    clustersSource = cms.InputTag( "hltParticleFlowClusterHBHE" ),
    pfClusterBuilder = cms.PSet( 
      allCellsPositionCalc = cms.PSet( 
        minAllowedNormalization = cms.double( 1.0E-9 ),
        posCalcNCrystals = cms.int32( -1 ),
        algoName = cms.string( "Basic2DGenericPFlowPositionCalc" ),
        minFractionInCalc = cms.double( 1.0E-9 ),
        logWeightDenominatorByDetector = cms.VPSet( 
          cms.PSet(  depths = cms.vint32( 1, 2, 3, 4 ),
            detector = cms.string( "HCAL_BARREL1" ),
            logWeightDenominator = cms.vdouble( 0.1, 0.2, 0.3, 0.3 )
          ),
          cms.PSet(  depths = cms.vint32( 1, 2, 3, 4, 5, 6, 7 ),
            detector = cms.string( "HCAL_ENDCAP" ),
            logWeightDenominator = cms.vdouble( 0.1, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2 )
          )
        )
      ),
      algoName = cms.string( "PFMultiDepthClusterizer" ),
      nSigmaPhi = cms.double( 2.0 ),
      minFractionToKeep = cms.double( 1.0E-7 ),
      nSigmaEta = cms.double( 2.0 )
    ),
    positionReCalc = cms.PSet(  ),
    energyCorrector = cms.PSet(  )
)
process.hltParticleFlowClusterHF = cms.EDProducer( "PFClusterProducer",
    recHitsSource = cms.InputTag( "hltParticleFlowRecHitHF" ),
    recHitCleaners = cms.VPSet( 
    ),
    seedCleaners = cms.VPSet( 
    ),
    seedFinder = cms.PSet( 
      thresholdsByDetector = cms.VPSet( 
        cms.PSet(  seedingThresholdPt = cms.double( 0.0 ),
          seedingThreshold = cms.double( 1.4 ),
          detector = cms.string( "HF_EM" )
        ),
        cms.PSet(  seedingThresholdPt = cms.double( 0.0 ),
          seedingThreshold = cms.double( 1.4 ),
          detector = cms.string( "HF_HAD" )
        )
      ),
      algoName = cms.string( "LocalMaximumSeedFinder" ),
      nNeighbours = cms.int32( 0 )
    ),
    initialClusteringStep = cms.PSet( 
      thresholdsByDetector = cms.VPSet( 
      ),
      algoName = cms.string( "Basic2DClusterForEachSeed" )
    ),
    pfClusterBuilder = cms.PSet(  ),
    positionReCalc = cms.PSet(  ),
    energyCorrector = cms.PSet(  )
)
process.hltLightPixelOnlyPFTracks = cms.EDProducer( "LightPFTrackProducer",
    TrackQuality = cms.string( "none" ),
    UseQuality = cms.bool( False ),
    TkColList = cms.VInputTag( 'hltPixelOnlyPFMuonMerging' )
)
process.hltPixelOnlyParticleFlowBlock = cms.EDProducer( "PFBlockProducer",
    verbose = cms.untracked.bool( False ),
    debug = cms.untracked.bool( False ),
    elementImporters = cms.VPSet( 
      cms.PSet(  muonSrc = cms.InputTag( "hltPixelOnlyMuons" ),
        source = cms.InputTag( "hltLightPixelOnlyPFTracks" ),
        NHitCuts_byTrackAlgo = cms.vuint32( 3, 3, 3, 3, 3, 3 ),
        useIterativeTracking = cms.bool( False ),
        importerName = cms.string( "GeneralTracksImporter" ),
        DPtOverPtCuts_byTrackAlgo = cms.vdouble( 5.0, 5.0, 5.0, 5.0, 5.0, 5.0 ),
        muonMaxDPtOPt = cms.double( 1.0 ),
        trackQuality = cms.string( "highPurity" ),
        cleanBadConvertedBrems = cms.bool( False ),
        vetoEndcap = cms.bool( False )
      ),
      cms.PSet(  source = cms.InputTag( "hltParticleFlowClusterECALUnseeded" ),
        importerName = cms.string( "ECALClusterImporter" ),
        BCtoPFCMap = cms.InputTag( "" )
      ),
      cms.PSet(  source = cms.InputTag( "hltParticleFlowClusterHCAL" ),
        importerName = cms.string( "GenericClusterImporter" )
      ),
      cms.PSet(  source = cms.InputTag( "hltParticleFlowClusterHF" ),
        importerName = cms.string( "GenericClusterImporter" )
      )
    ),
    linkDefinitions = cms.VPSet( 
      cms.PSet(  linkType = cms.string( "TRACK:ECAL" ),
        useKDTree = cms.bool( True ),
        linkerName = cms.string( "TrackAndECALLinker" )
      ),
      cms.PSet(  linkType = cms.string( "TRACK:HCAL" ),
        useKDTree = cms.bool( True ),
        linkerName = cms.string( "TrackAndHCALLinker" ),
        trajectoryLayerEntrance = cms.string( "HCALEntrance" ),
        trajectoryLayerExit = cms.string( "HCALExit" ),
        nMaxHcalLinksPerTrack = cms.int32( 1 )
      ),
      cms.PSet(  linkType = cms.string( "ECAL:HCAL" ),
        useKDTree = cms.bool( False ),
        linkerName = cms.string( "ECALAndHCALLinker" ),
        minAbsEtaEcal = cms.double( 2.5 )
      ),
      cms.PSet(  linkType = cms.string( "HFEM:HFHAD" ),
        useKDTree = cms.bool( False ),
        linkerName = cms.string( "HFEMAndHFHADLinker" )
      )
    )
)
process.hltPixelOnlyParticleFlow = cms.EDProducer( "PFProducer",
    verbose = cms.untracked.bool( False ),
    debug = cms.untracked.bool( False ),
    blocks = cms.InputTag( "hltPixelOnlyParticleFlowBlock" ),
    muons = cms.InputTag( "hltPixelOnlyMuons" ),
    postMuonCleaning = cms.bool( True ),
    vetoEndcap = cms.bool( False ),
    vertexCollection = cms.InputTag( "hltPixelVertices" ),
    useVerticesForNeutral = cms.bool( True ),
    useHO = cms.bool( False ),
    PFEGammaCandidates = cms.InputTag( "particleFlowEGamma" ),
    GedElectronValueMap = cms.InputTag( "gedGsfElectronsTmp" ),
    GedPhotonValueMap = cms.InputTag( 'tmpGedPhotons','valMapPFEgammaCandToPhoton' ),
    useEGammaElectrons = cms.bool( False ),
    egammaElectrons = cms.InputTag( "" ),
    useEGammaFilters = cms.bool( False ),
    useProtectionsForJetMET = cms.bool( True ),
    PFEGammaFiltersParameters = cms.PSet( 
      electron_protectionsForJetMET = cms.PSet( 
        maxE = cms.double( 50.0 ),
        maxTrackPOverEele = cms.double( 1.0 ),
        maxEcalEOverP_2 = cms.double( 0.2 ),
        maxHcalEOverEcalE = cms.double( 0.1 ),
        maxEcalEOverP_1 = cms.double( 0.5 ),
        maxHcalEOverP = cms.double( 1.0 ),
        maxEcalEOverPRes = cms.double( 0.2 ),
        maxHcalE = cms.double( 10.0 ),
        maxEeleOverPout = cms.double( 0.2 ),
        maxNtracks = cms.double( 3.0 ),
        maxEleHcalEOverEcalE = cms.double( 0.1 ),
        maxDPhiIN = cms.double( 0.1 ),
        maxEeleOverPoutRes = cms.double( 0.5 )
      ),
      electron_maxElePtForOnlyMVAPresel = cms.double( 50.0 ),
      photon_SigmaiEtaiEta_endcap = cms.double( 0.034 ),
      electron_iso_combIso_endcap = cms.double( 10.0 ),
      photon_protectionsForBadHcal = cms.PSet( 
        solidConeTrkIsoSlope = cms.double( 0.3 ),
        enableProtections = cms.bool( False ),
        solidConeTrkIsoOffset = cms.double( 10.0 )
      ),
      electron_missinghits = cms.uint32( 1 ),
      photon_MinEt = cms.double( 10.0 ),
      electron_iso_pt = cms.double( 10.0 ),
      electron_ecalDrivenHademPreselCut = cms.double( 0.15 ),
      electron_iso_mva_endcap = cms.double( -0.1075 ),
      electron_iso_combIso_barrel = cms.double( 10.0 ),
      photon_protectionsForJetMET = cms.PSet( 
        sumPtTrackIsoSlope = cms.double( 0.001 ),
        sumPtTrackIso = cms.double( 4.0 )
      ),
      electron_protectionsForBadHcal = cms.PSet( 
        dEta = cms.vdouble( 0.0064, 0.01264 ),
        dPhi = cms.vdouble( 0.0547, 0.0394 ),
        enableProtections = cms.bool( False ),
        eInvPInv = cms.vdouble( 0.184, 0.0721 ),
        full5x5_sigmaIetaIeta = cms.vdouble( 0.0106, 0.0387 )
      ),
      electron_noniso_mvaCut = cms.double( -0.1 ),
      electron_iso_mva_barrel = cms.double( -0.1875 ),
      photon_SigmaiEtaiEta_barrel = cms.double( 0.0125 ),
      photon_combIso = cms.double( 10.0 ),
      photon_HoE = cms.double( 0.05 )
    ),
    muon_HCAL = cms.vdouble( 3.0, 3.0 ),
    muon_ECAL = cms.vdouble( 0.5, 0.5 ),
    muon_HO = cms.vdouble( 0.9, 0.9 ),
    PFMuonAlgoParameters = cms.PSet(  ),
    rejectTracks_Bad = cms.bool( False ),
    rejectTracks_Step45 = cms.bool( False ),
    usePFNuclearInteractions = cms.bool( False ),
    usePFConversions = cms.bool( False ),
    usePFDecays = cms.bool( False ),
    dptRel_DispVtx = cms.double( 10.0 ),
    iCfgCandConnector = cms.PSet( 
      nuclCalibFactors = cms.vdouble( 0.8, 0.15, 0.5, 0.5, 0.05 ),
      bCorrect = cms.bool( False ),
      bCalibPrimary = cms.bool( False )
    ),
    nsigma_TRACK = cms.double( 1.0 ),
    pt_Error = cms.double( 1.0 ),
    factors_45 = cms.vdouble( 10.0, 100.0 ),
    goodTrackDeadHcal_ptErrRel = cms.double( 0.2 ),
    goodTrackDeadHcal_chi2n = cms.double( 5.0 ),
    goodTrackDeadHcal_layers = cms.uint32( 4 ),
    goodTrackDeadHcal_validFr = cms.double( 0.5 ),
    goodTrackDeadHcal_dxy = cms.double( 0.5 ),
    goodPixelTrackDeadHcal_minEta = cms.double( 2.3 ),
    goodPixelTrackDeadHcal_maxPt = cms.double( 50.0 ),
    goodPixelTrackDeadHcal_ptErrRel = cms.double( 1.0 ),
    goodPixelTrackDeadHcal_chi2n = cms.double( 2.0 ),
    goodPixelTrackDeadHcal_maxLost3Hit = cms.int32( 0 ),
    goodPixelTrackDeadHcal_maxLost4Hit = cms.int32( 1 ),
    goodPixelTrackDeadHcal_dxy = cms.double( 0.02 ),
    goodPixelTrackDeadHcal_dz = cms.double( 0.05 ),
    pf_nsigma_ECAL = cms.double( 0.0 ),
    pf_nsigma_HCAL = cms.double( 1.0 ),
    pf_nsigma_HFEM = cms.double( 1.0 ),
    pf_nsigma_HFHAD = cms.double( 1.0 ),
    useCalibrationsFromDB = cms.bool( True ),
    calibrationsLabel = cms.string( "HLT" ),
    postHFCleaning = cms.bool( False ),
    PFHFCleaningParameters = cms.PSet( 
      minSignificance = cms.double( 2.5 ),
      maxSignificance = cms.double( 2.5 ),
      minDeltaMet = cms.double( 0.4 ),
      maxDeltaPhiPt = cms.double( 7.0 ),
      minHFCleaningPt = cms.double( 5.0 ),
      minSignificanceReduction = cms.double( 1.4 )
    ),
    cleanedHF = cms.VInputTag( 'hltParticleFlowRecHitHF:Cleaned','hltParticleFlowClusterHF:Cleaned' ),
    calibHF_use = cms.bool( False ),
    calibHF_eta_step = cms.vdouble( 0.0, 2.9, 3.0, 3.2, 4.2, 4.4, 4.6, 4.8, 5.2, 5.4 ),
    calibHF_a_EMonly = cms.vdouble( 0.96945, 0.96701, 0.76309, 0.82268, 0.87583, 0.89718, 0.98674, 1.4681, 1.458, 1.458 ),
    calibHF_a_EMHAD = cms.vdouble( 1.42215, 1.00496, 0.68961, 0.81656, 0.98504, 0.98504, 1.00802, 1.0593, 1.4576, 1.4576 ),
    calibHF_b_HADonly = cms.vdouble( 1.27541, 0.85361, 0.86333, 0.89091, 0.94348, 0.94348, 0.9437, 1.0034, 1.0444, 1.0444 ),
    calibHF_b_EMHAD = cms.vdouble( 1.27541, 0.85361, 0.86333, 0.89091, 0.94348, 0.94348, 0.9437, 1.0034, 1.0444, 1.0444 ),
    resolHF_square = cms.vdouble( 7.834401, 0.012996, 0.0 )
)
process.hltAK4PixelOnlyPFJets = cms.EDProducer( "FastjetJetProducer",
    useMassDropTagger = cms.bool( False ),
    useFiltering = cms.bool( False ),
    useDynamicFiltering = cms.bool( False ),
    useTrimming = cms.bool( False ),
    usePruning = cms.bool( False ),
    useCMSBoostedTauSeedingAlgorithm = cms.bool( False ),
    useKtPruning = cms.bool( False ),
    useConstituentSubtraction = cms.bool( False ),
    useSoftDrop = cms.bool( False ),
    correctShape = cms.bool( False ),
    UseOnlyVertexTracks = cms.bool( False ),
    UseOnlyOnePV = cms.bool( False ),
    muCut = cms.double( -1.0 ),
    yCut = cms.double( -1.0 ),
    rFilt = cms.double( -1.0 ),
    rFiltFactor = cms.double( -1.0 ),
    trimPtFracMin = cms.double( -1.0 ),
    zcut = cms.double( -1.0 ),
    rcut_factor = cms.double( -1.0 ),
    csRho_EtaMax = cms.double( -1.0 ),
    csRParam = cms.double( -1.0 ),
    beta = cms.double( -1.0 ),
    R0 = cms.double( -1.0 ),
    gridMaxRapidity = cms.double( -1.0 ),
    gridSpacing = cms.double( -1.0 ),
    DzTrVtxMax = cms.double( 0.0 ),
    DxyTrVtxMax = cms.double( 0.0 ),
    MaxVtxZ = cms.double( 15.0 ),
    subjetPtMin = cms.double( -1.0 ),
    muMin = cms.double( -1.0 ),
    muMax = cms.double( -1.0 ),
    yMin = cms.double( -1.0 ),
    yMax = cms.double( -1.0 ),
    dRMin = cms.double( -1.0 ),
    dRMax = cms.double( -1.0 ),
    maxDepth = cms.int32( -1 ),
    nFilt = cms.int32( -1 ),
    MinVtxNdof = cms.int32( 0 ),
    src = cms.InputTag( "hltPixelOnlyParticleFlow" ),
    srcPVs = cms.InputTag( "hltPixelVertices" ),
    jetType = cms.string( "PFJet" ),
    jetAlgorithm = cms.string( "AntiKt" ),
    rParam = cms.double( 0.4 ),
    inputEtMin = cms.double( 0.0 ),
    inputEMin = cms.double( 0.0 ),
    jetPtMin = cms.double( 0.0 ),
    doPVCorrection = cms.bool( False ),
    doAreaFastjet = cms.bool( False ),
    doRhoFastjet = cms.bool( False ),
    doPUOffsetCorr = cms.bool( False ),
    puPtMin = cms.double( 10.0 ),
    nSigmaPU = cms.double( 1.0 ),
    radiusPU = cms.double( 0.4 ),
    subtractorName = cms.string( "" ),
    useExplicitGhosts = cms.bool( False ),
    doAreaDiskApprox = cms.bool( True ),
    voronoiRfact = cms.double( -9.0 ),
    Rho_EtaMax = cms.double( 4.4 ),
    Ghost_EtaMax = cms.double( 6.0 ),
    Active_Area_Repeats = cms.int32( 5 ),
    GhostArea = cms.double( 0.01 ),
    restrictInputs = cms.bool( False ),
    maxInputs = cms.uint32( 1 ),
    writeCompound = cms.bool( False ),
    writeJetsWithConst = cms.bool( False ),
    doFastJetNonUniform = cms.bool( False ),
    useDeterministicSeed = cms.bool( True ),
    minSeed = cms.uint32( 0 ),
    verbosity = cms.int32( 0 ),
    puWidth = cms.double( 0.0 ),
    nExclude = cms.uint32( 0 ),
    maxBadEcalCells = cms.uint32( 9999999 ),
    maxBadHcalCells = cms.uint32( 9999999 ),
    maxProblematicEcalCells = cms.uint32( 9999999 ),
    maxProblematicHcalCells = cms.uint32( 9999999 ),
    maxRecoveredEcalCells = cms.uint32( 9999999 ),
    maxRecoveredHcalCells = cms.uint32( 9999999 ),
    puCenters = cms.vdouble(  ),
    applyWeight = cms.bool( False ),
    srcWeights = cms.InputTag( "" ),
    minimumTowersFraction = cms.double( 0.0 ),
    jetCollInstanceName = cms.string( "" ),
    sumRecHits = cms.bool( False )
)
process.hltAK4PixelOnlyPFJetsLooseID = cms.EDProducer( "HLTPFJetIDProducer",
    minPt = cms.double( 20.0 ),
    maxEta = cms.double( 1.0E99 ),
    CHF = cms.double( 0.0 ),
    NHF = cms.double( 0.99 ),
    CEF = cms.double( 0.99 ),
    NEF = cms.double( 0.99 ),
    maxCF = cms.double( 99.0 ),
    NCH = cms.int32( 0 ),
    NTOT = cms.int32( 1 ),
    jetsInput = cms.InputTag( "hltAK4PixelOnlyPFJets" )
)
process.hltAK4PixelOnlyPFJetsTightID = cms.EDProducer( "HLTPFJetIDProducer",
    minPt = cms.double( 20.0 ),
    maxEta = cms.double( 1.0E99 ),
    CHF = cms.double( 0.0 ),
    NHF = cms.double( 0.9 ),
    CEF = cms.double( 0.99 ),
    NEF = cms.double( 0.99 ),
    maxCF = cms.double( 99.0 ),
    NCH = cms.int32( 0 ),
    NTOT = cms.int32( 1 ),
    jetsInput = cms.InputTag( "hltAK4PixelOnlyPFJets" )
)
process.hltFixedGridRhoFastjetPixelOnlyAll = cms.EDProducer( "FixedGridRhoProducerFastjet",
    pfCandidatesTag = cms.InputTag( "hltPixelOnlyParticleFlow" ),
    maxRapidity = cms.double( 5.0 ),
    gridSpacing = cms.double( 0.55 )
)
process.hltAK4PixelOnlyPFFastJetCorrector = cms.EDProducer( "L1FastjetCorrectorProducer",
    level = cms.string( "L1FastJet" ),
    algorithm = cms.string( "AK4PFHLT" ),
    srcRho = cms.InputTag( "hltFixedGridRhoFastjetPixelOnlyAll" )
)
process.hltAK4PixelOnlyPFRelativeCorrector = cms.EDProducer( "LXXXCorrectorProducer",
    level = cms.string( "L2Relative" ),
    algorithm = cms.string( "AK4PFHLT" )
)
process.hltAK4PixelOnlyPFAbsoluteCorrector = cms.EDProducer( "LXXXCorrectorProducer",
    level = cms.string( "L3Absolute" ),
    algorithm = cms.string( "AK4PFHLT" )
)
process.hltAK4PixelOnlyPFResidualCorrector = cms.EDProducer( "LXXXCorrectorProducer",
    level = cms.string( "L2L3Residual" ),
    algorithm = cms.string( "AK4PFHLT" )
)
process.hltAK4PixelOnlyPFCorrector = cms.EDProducer( "ChainedJetCorrectorProducer",
    correctors = cms.VInputTag( 'hltAK4PixelOnlyPFFastJetCorrector','hltAK4PixelOnlyPFRelativeCorrector','hltAK4PixelOnlyPFAbsoluteCorrector','hltAK4PixelOnlyPFResidualCorrector' )
)
process.hltAK4PixelOnlyPFJetsCorrected = cms.EDProducer( "CorrectedPFJetProducer",
    src = cms.InputTag( "hltAK4PixelOnlyPFJets" ),
    correctors = cms.VInputTag( 'hltAK4PixelOnlyPFCorrector' )
)
process.hltAK4PixelOnlyPFJetsLooseIDCorrected = cms.EDProducer( "CorrectedPFJetProducer",
    src = cms.InputTag( "hltAK4PixelOnlyPFJetsLooseID" ),
    correctors = cms.VInputTag( 'hltAK4PixelOnlyPFCorrector' )
)
process.hltAK4PixelOnlyPFJetsTightIDCorrected = cms.EDProducer( "CorrectedPFJetProducer",
    src = cms.InputTag( "hltAK4PixelOnlyPFJetsTightID" ),
    correctors = cms.VInputTag( 'hltAK4PixelOnlyPFCorrector' )
)
process.hltPixelOnlyPFMETProducer = cms.EDProducer( "PFMETProducer",
    src = cms.InputTag( "hltPixelOnlyParticleFlow" ),
    globalThreshold = cms.double( 0.0 ),
    alias = cms.string( "hltPFMet" ),
    calculateSignificance = cms.bool( False ),
    parameters = cms.PSet(  ),
    applyWeight = cms.bool( False ),
    srcWeights = cms.InputTag( "" )
)
process.hltFixedGridRhoFastjetECALMFForMuons = cms.EDProducer( "FixedGridRhoProducerFastjetFromRecHit",
    hbheRecHitsTag = cms.InputTag( "hltHbhereco" ),
    ebRecHitsTag = cms.InputTag( 'hltEcalRecHit','EcalRecHitsEB' ),
    eeRecHitsTag = cms.InputTag( 'hltEcalRecHit','EcalRecHitsEE' ),
    skipHCAL = cms.bool( True ),
    skipECAL = cms.bool( False ),
    eThresHB = cms.vdouble( 0.1, 0.2, 0.3, 0.3 ),
    eThresHE = cms.vdouble( 0.1, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2 ),
    maxRapidity = cms.double( 2.5 ),
    gridSpacing = cms.double( 0.55 )
)
process.hltFixedGridRhoFastjetHCAL = cms.EDProducer( "FixedGridRhoProducerFastjetFromRecHit",
    hbheRecHitsTag = cms.InputTag( "hltHbhereco" ),
    ebRecHitsTag = cms.InputTag( 'hltEcalRecHit','EcalRecHitsEB' ),
    eeRecHitsTag = cms.InputTag( 'hltEcalRecHit','EcalRecHitsEE' ),
    skipHCAL = cms.bool( False ),
    skipECAL = cms.bool( True ),
    eThresHB = cms.vdouble( 0.1, 0.2, 0.3, 0.3 ),
    eThresHE = cms.vdouble( 0.1, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2 ),
    maxRapidity = cms.double( 2.5 ),
    gridSpacing = cms.double( 0.55 )
)
process.hltRecHitInRegionForMuonsMF = cms.EDProducer( "MuonHLTRechitInRegionsProducer",
    productLabels = cms.vstring( 'EcalRegionalRecHitsEB',
      'EcalRegionalRecHitsEE' ),
    ecalhitLabels = cms.VInputTag( 'hltEcalRecHit:EcalRecHitsEB','hltEcalRecHit:EcalRecHitsEE' ),
    l1TagIsolated = cms.InputTag( "hltIterL3MuonCandidates" ),
    l1TagNonIsolated = cms.InputTag( "NotUsed" ),
    useUncalib = cms.bool( False ),
    doIsolated = cms.bool( True ),
    l1LowerThr = cms.double( 0.0 ),
    l1UpperThr = cms.double( 999.0 ),
    l1LowerThrIgnoreIsolation = cms.double( 100.0 ),
    regionEtaMargin = cms.double( 0.4 ),
    regionPhiMargin = cms.double( 0.4 )
)
process.hltRecHitInRegionForMuonsES = cms.EDProducer( "MuonHLTRechitInRegionsProducer",
    productLabels = cms.vstring( 'EcalRegionalRecHitsES' ),
    ecalhitLabels = cms.VInputTag( 'hltEcalPreshowerRecHit:EcalRecHitsES' ),
    l1TagIsolated = cms.InputTag( "hltIterL3MuonCandidates" ),
    l1TagNonIsolated = cms.InputTag( "NotUsed" ),
    useUncalib = cms.bool( False ),
    doIsolated = cms.bool( True ),
    l1LowerThr = cms.double( 0.0 ),
    l1UpperThr = cms.double( 999.0 ),
    l1LowerThrIgnoreIsolation = cms.double( 100.0 ),
    regionEtaMargin = cms.double( 0.4 ),
    regionPhiMargin = cms.double( 0.4 )
)
process.hltParticleFlowRecHitECALForMuonsMF = cms.EDProducer( "PFRecHitProducer",
    navigator = cms.PSet( 
      barrel = cms.PSet(  ),
      endcap = cms.PSet(  ),
      name = cms.string( "PFRecHitECALNavigator" )
    ),
    producers = cms.VPSet( 
      cms.PSet(  src = cms.InputTag( 'hltRecHitInRegionForMuonsMF','EcalRegionalRecHitsEB' ),
        srFlags = cms.InputTag( "" ),
        name = cms.string( "PFEBRecHitCreator" ),
        qualityTests = cms.VPSet( 
          cms.PSet(  name = cms.string( "PFRecHitQTestDBThreshold" ),
            applySelectionsToAllCrystals = cms.bool( True )
          ),
          cms.PSet(  topologicalCleaning = cms.bool( True ),
            skipTTRecoveredHits = cms.bool( True ),
            cleaningThreshold = cms.double( 2.0 ),
            name = cms.string( "PFRecHitQTestECAL" ),
            timingCleaning = cms.bool( True )
          )
        )
      ),
      cms.PSet(  src = cms.InputTag( 'hltRecHitInRegionForMuonsMF','EcalRegionalRecHitsEE' ),
        srFlags = cms.InputTag( "" ),
        name = cms.string( "PFEERecHitCreator" ),
        qualityTests = cms.VPSet( 
          cms.PSet(  name = cms.string( "PFRecHitQTestDBThreshold" ),
            applySelectionsToAllCrystals = cms.bool( True )
          ),
          cms.PSet(  topologicalCleaning = cms.bool( True ),
            skipTTRecoveredHits = cms.bool( True ),
            cleaningThreshold = cms.double( 2.0 ),
            name = cms.string( "PFRecHitQTestECAL" ),
            timingCleaning = cms.bool( True )
          )
        )
      )
    )
)
process.hltParticleFlowRecHitPSForMuons = cms.EDProducer( "PFRecHitProducer",
    navigator = cms.PSet(  name = cms.string( "PFRecHitPreshowerNavigator" ) ),
    producers = cms.VPSet( 
      cms.PSet(  src = cms.InputTag( 'hltRecHitInRegionForMuonsES','EcalRegionalRecHitsES' ),
        name = cms.string( "PFPSRecHitCreator" ),
        qualityTests = cms.VPSet( 
          cms.PSet(  threshold = cms.double( 7.0E-6 ),
            name = cms.string( "PFRecHitQTestThreshold" )
          )
        )
      )
    )
)
process.hltParticleFlowClusterECALUncorrectedForMuonsMF = cms.EDProducer( "PFClusterProducer",
    recHitsSource = cms.InputTag( "hltParticleFlowRecHitECALForMuonsMF" ),
    recHitCleaners = cms.VPSet( 
    ),
    seedCleaners = cms.VPSet( 
    ),
    seedFinder = cms.PSet( 
      thresholdsByDetector = cms.VPSet( 
        cms.PSet(  seedingThresholdPt = cms.double( 0.15 ),
          seedingThreshold = cms.double( 0.6 ),
          detector = cms.string( "ECAL_ENDCAP" )
        ),
        cms.PSet(  seedingThresholdPt = cms.double( 0.0 ),
          seedingThreshold = cms.double( 0.23 ),
          detector = cms.string( "ECAL_BARREL" )
        )
      ),
      algoName = cms.string( "LocalMaximumSeedFinder" ),
      nNeighbours = cms.int32( 8 )
    ),
    initialClusteringStep = cms.PSet( 
      thresholdsByDetector = cms.VPSet( 
        cms.PSet(  gatheringThreshold = cms.double( 0.08 ),
          gatheringThresholdPt = cms.double( 0.0 ),
          detector = cms.string( "ECAL_BARREL" )
        ),
        cms.PSet(  gatheringThreshold = cms.double( 0.3 ),
          gatheringThresholdPt = cms.double( 0.0 ),
          detector = cms.string( "ECAL_ENDCAP" )
        )
      ),
      algoName = cms.string( "Basic2DGenericTopoClusterizer" ),
      useCornerCells = cms.bool( True )
    ),
    pfClusterBuilder = cms.PSet( 
      minFracTot = cms.double( 1.0E-20 ),
      stoppingTolerance = cms.double( 1.0E-8 ),
      positionCalc = cms.PSet( 
        minAllowedNormalization = cms.double( 1.0E-9 ),
        posCalcNCrystals = cms.int32( 9 ),
        algoName = cms.string( "Basic2DGenericPFlowPositionCalc" ),
        logWeightDenominator = cms.double( 0.08 ),
        minFractionInCalc = cms.double( 1.0E-9 ),
        timeResolutionCalcBarrel = cms.PSet( 
          corrTermLowE = cms.double( 0.0510871 ),
          threshLowE = cms.double( 0.5 ),
          noiseTerm = cms.double( 1.10889 ),
          constantTermLowE = cms.double( 0.0 ),
          noiseTermLowE = cms.double( 1.31883 ),
          threshHighE = cms.double( 5.0 ),
          constantTerm = cms.double( 0.428192 )
        ),
        timeResolutionCalcEndcap = cms.PSet( 
          corrTermLowE = cms.double( 0.0 ),
          threshLowE = cms.double( 1.0 ),
          noiseTerm = cms.double( 5.72489999999 ),
          constantTermLowE = cms.double( 0.0 ),
          noiseTermLowE = cms.double( 6.92683000001 ),
          threshHighE = cms.double( 10.0 ),
          constantTerm = cms.double( 0.0 )
        )
      ),
      maxIterations = cms.uint32( 50 ),
      positionCalcForConvergence = cms.PSet( 
        minAllowedNormalization = cms.double( 0.0 ),
        T0_ES = cms.double( 1.2 ),
        algoName = cms.string( "ECAL2DPositionCalcWithDepthCorr" ),
        T0_EE = cms.double( 3.1 ),
        T0_EB = cms.double( 7.4 ),
        X0 = cms.double( 0.89 ),
        minFractionInCalc = cms.double( 0.0 ),
        W0 = cms.double( 4.2 )
      ),
      allCellsPositionCalc = cms.PSet( 
        minAllowedNormalization = cms.double( 1.0E-9 ),
        posCalcNCrystals = cms.int32( -1 ),
        algoName = cms.string( "Basic2DGenericPFlowPositionCalc" ),
        logWeightDenominator = cms.double( 0.08 ),
        minFractionInCalc = cms.double( 1.0E-9 ),
        timeResolutionCalcBarrel = cms.PSet( 
          corrTermLowE = cms.double( 0.0510871 ),
          threshLowE = cms.double( 0.5 ),
          noiseTerm = cms.double( 1.10889 ),
          constantTermLowE = cms.double( 0.0 ),
          noiseTermLowE = cms.double( 1.31883 ),
          threshHighE = cms.double( 5.0 ),
          constantTerm = cms.double( 0.428192 )
        ),
        timeResolutionCalcEndcap = cms.PSet( 
          corrTermLowE = cms.double( 0.0 ),
          threshLowE = cms.double( 1.0 ),
          noiseTerm = cms.double( 5.72489999999 ),
          constantTermLowE = cms.double( 0.0 ),
          noiseTermLowE = cms.double( 6.92683000001 ),
          threshHighE = cms.double( 10.0 ),
          constantTerm = cms.double( 0.0 )
        )
      ),
      algoName = cms.string( "Basic2DGenericPFlowClusterizer" ),
      recHitEnergyNorms = cms.VPSet( 
        cms.PSet(  recHitEnergyNorm = cms.double( 0.08 ),
          detector = cms.string( "ECAL_BARREL" )
        ),
        cms.PSet(  recHitEnergyNorm = cms.double( 0.3 ),
          detector = cms.string( "ECAL_ENDCAP" )
        )
      ),
      showerSigma = cms.double( 1.5 ),
      minFractionToKeep = cms.double( 1.0E-7 ),
      excludeOtherSeeds = cms.bool( True )
    ),
    positionReCalc = cms.PSet( 
      minAllowedNormalization = cms.double( 0.0 ),
      T0_ES = cms.double( 1.2 ),
      algoName = cms.string( "ECAL2DPositionCalcWithDepthCorr" ),
      T0_EE = cms.double( 3.1 ),
      T0_EB = cms.double( 7.4 ),
      X0 = cms.double( 0.89 ),
      minFractionInCalc = cms.double( 0.0 ),
      W0 = cms.double( 4.2 )
    ),
    energyCorrector = cms.PSet(  )
)
process.hltParticleFlowClusterPSForMuons = cms.EDProducer( "PFClusterProducer",
    recHitsSource = cms.InputTag( "hltParticleFlowRecHitPSForMuons" ),
    recHitCleaners = cms.VPSet( 
    ),
    seedCleaners = cms.VPSet( 
    ),
    seedFinder = cms.PSet( 
      thresholdsByDetector = cms.VPSet( 
        cms.PSet(  seedingThresholdPt = cms.double( 0.0 ),
          seedingThreshold = cms.double( 1.2E-4 ),
          detector = cms.string( "PS1" )
        ),
        cms.PSet(  seedingThresholdPt = cms.double( 0.0 ),
          seedingThreshold = cms.double( 1.2E-4 ),
          detector = cms.string( "PS2" )
        )
      ),
      algoName = cms.string( "LocalMaximumSeedFinder" ),
      nNeighbours = cms.int32( 4 )
    ),
    initialClusteringStep = cms.PSet( 
      thresholdsByDetector = cms.VPSet( 
        cms.PSet(  gatheringThreshold = cms.double( 6.0E-5 ),
          gatheringThresholdPt = cms.double( 0.0 ),
          detector = cms.string( "PS1" )
        ),
        cms.PSet(  gatheringThreshold = cms.double( 6.0E-5 ),
          gatheringThresholdPt = cms.double( 0.0 ),
          detector = cms.string( "PS2" )
        )
      ),
      algoName = cms.string( "Basic2DGenericTopoClusterizer" ),
      useCornerCells = cms.bool( False )
    ),
    pfClusterBuilder = cms.PSet( 
      minFracTot = cms.double( 1.0E-20 ),
      stoppingTolerance = cms.double( 1.0E-8 ),
      positionCalc = cms.PSet( 
        minAllowedNormalization = cms.double( 1.0E-9 ),
        posCalcNCrystals = cms.int32( -1 ),
        algoName = cms.string( "Basic2DGenericPFlowPositionCalc" ),
        logWeightDenominator = cms.double( 6.0E-5 ),
        minFractionInCalc = cms.double( 1.0E-9 )
      ),
      maxIterations = cms.uint32( 50 ),
      algoName = cms.string( "Basic2DGenericPFlowClusterizer" ),
      recHitEnergyNorms = cms.VPSet( 
        cms.PSet(  recHitEnergyNorm = cms.double( 6.0E-5 ),
          detector = cms.string( "PS1" )
        ),
        cms.PSet(  recHitEnergyNorm = cms.double( 6.0E-5 ),
          detector = cms.string( "PS2" )
        )
      ),
      showerSigma = cms.double( 0.3 ),
      minFractionToKeep = cms.double( 1.0E-7 ),
      excludeOtherSeeds = cms.bool( True )
    ),
    positionReCalc = cms.PSet(  ),
    energyCorrector = cms.PSet(  )
)
process.hltParticleFlowClusterECALForMuonsMF = cms.EDProducer( "CorrectedECALPFClusterProducer",
    minimumPSEnergy = cms.double( 0.0 ),
    skipPS = cms.bool( False ),
    inputPS = cms.InputTag( "hltParticleFlowClusterPSForMuons" ),
    energyCorrector = cms.PSet( 
      applyCrackCorrections = cms.bool( False ),
      srfAwareCorrection = cms.bool( True ),
      applyMVACorrections = cms.bool( True ),
      maxPtForMVAEvaluation = cms.double( 300.0 ),
      recHitsEBLabel = cms.InputTag( 'hltEcalRecHit','EcalRecHitsEB' ),
      recHitsEELabel = cms.InputTag( 'hltEcalRecHit','EcalRecHitsEE' ),
      ebSrFlagLabel = cms.InputTag( "hltEcalDigis" ),
      eeSrFlagLabel = cms.InputTag( "hltEcalDigis" )
    ),
    inputECAL = cms.InputTag( "hltParticleFlowClusterECALUncorrectedForMuonsMF" )
)
process.hltMuonEcalMFPFClusterIsoForMuons = cms.EDProducer( "MuonHLTEcalPFClusterIsolationProducer",
    recoCandidateProducer = cms.InputTag( "hltIterL3MuonCandidates" ),
    pfClusterProducer = cms.InputTag( "hltParticleFlowClusterECALForMuonsMF" ),
    rhoProducer = cms.InputTag( "hltFixedGridRhoFastjetECALMFForMuons" ),
    doRhoCorrection = cms.bool( True ),
    rhoMax = cms.double( 9.9999999E7 ),
    rhoScale = cms.double( 1.0 ),
    drMax = cms.double( 0.3 ),
    drVetoBarrel = cms.double( 0.05 ),
    drVetoEndcap = cms.double( 0.05 ),
    etaStripBarrel = cms.double( 0.0 ),
    etaStripEndcap = cms.double( 0.0 ),
    energyBarrel = cms.double( 0.0 ),
    energyEndcap = cms.double( 0.0 ),
    effectiveAreas = cms.vdouble( 0.35, 0.193 ),
    absEtaLowEdges = cms.vdouble( 0.0, 1.479 )
)
process.hltMuonHcalRegPFClusterIsoForMuons = cms.EDProducer( "MuonHLTHcalPFClusterIsolationProducer",
    recoCandidateProducer = cms.InputTag( "hltIterL3MuonCandidates" ),
    pfClusterProducerHCAL = cms.InputTag( "hltParticleFlowClusterHCAL" ),
    useHF = cms.bool( False ),
    pfClusterProducerHFEM = cms.InputTag( "" ),
    pfClusterProducerHFHAD = cms.InputTag( "" ),
    rhoProducer = cms.InputTag( "hltFixedGridRhoFastjetHCAL" ),
    doRhoCorrection = cms.bool( True ),
    rhoMax = cms.double( 9.9999999E7 ),
    rhoScale = cms.double( 1.0 ),
    drMax = cms.double( 0.3 ),
    drVetoBarrel = cms.double( 0.1 ),
    drVetoEndcap = cms.double( 0.1 ),
    etaStripBarrel = cms.double( 0.0 ),
    etaStripEndcap = cms.double( 0.0 ),
    energyBarrel = cms.double( 0.0 ),
    energyEndcap = cms.double( 0.0 ),
    useEt = cms.bool( True ),
    effectiveAreas = cms.vdouble( 0.227, 0.372 ),
    absEtaLowEdges = cms.vdouble( 0.0, 1.479 )
)
process.hltPixelTracksTrackingRegionsForSeedsL3Muon = cms.EDProducer( "CandidateSeededTrackingRegionsEDProducer",
    RegionPSet = cms.PSet( 
      vertexCollection = cms.InputTag( "hltPixelVertices" ),
      zErrorVetex = cms.double( 0.2 ),
      beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
      zErrorBeamSpot = cms.double( 24.2 ),
      maxNVertices = cms.int32( 4 ),
      maxNRegions = cms.int32( 40 ),
      nSigmaZVertex = cms.double( 3.0 ),
      nSigmaZBeamSpot = cms.double( 4.0 ),
      ptMin = cms.double( 0.3 ),
      mode = cms.string( "VerticesFixed" ),
      input = cms.InputTag( "hltIterL3MuonCandidates" ),
      searchOpt = cms.bool( False ),
      whereToUseMeasurementTracker = cms.string( "Never" ),
      originRadius = cms.double( 0.1 ),
      measurementTrackerName = cms.InputTag( "" ),
      precise = cms.bool( True ),
      deltaEta = cms.double( 0.3 ),
      deltaPhi = cms.double( 0.3 )
    )
)
process.hltPixelTracksInRegionIter0L3Muon = cms.EDProducer( "TrackSelectorByRegion",
    tracks = cms.InputTag( "hltPixelTracks" ),
    regions = cms.InputTag( "hltPixelTracksTrackingRegionsForSeedsL3Muon" ),
    produceTrackCollection = cms.bool( True ),
    produceMask = cms.bool( False )
)
process.hltIter0L3MuonPixelSeedsFromPixelTracks = cms.EDProducer( "SeedGeneratorFromProtoTracksEDProducer",
    InputCollection = cms.InputTag( "hltPixelTracksInRegionIter0L3Muon" ),
    InputVertexCollection = cms.InputTag( "hltTrimmedPixelVertices" ),
    originHalfLength = cms.double( 0.3 ),
    originRadius = cms.double( 0.1 ),
    useProtoTrackKinematics = cms.bool( False ),
    useEventsWithNoVertex = cms.bool( True ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    usePV = cms.bool( False ),
    includeFourthHit = cms.bool( True ),
    SeedCreatorPSet = cms.PSet(  refToPSet_ = cms.string( "HLTSeedFromProtoTracks" ) )
)
process.hltIter0L3MuonCkfTrackCandidates = cms.EDProducer( "CkfTrackCandidateMaker",
    cleanTrajectoryAfterInOut = cms.bool( False ),
    doSeedingRegionRebuilding = cms.bool( False ),
    onlyPixelHitsForSeedCleaner = cms.bool( False ),
    reverseTrajectories = cms.bool( False ),
    useHitsSplitting = cms.bool( False ),
    MeasurementTrackerEvent = cms.InputTag( "hltSiStripClusters" ),
    src = cms.InputTag( "hltIter0L3MuonPixelSeedsFromPixelTracks" ),
    clustersToSkip = cms.InputTag( "" ),
    phase2clustersToSkip = cms.InputTag( "" ),
    TrajectoryBuilderPSet = cms.PSet(  refToPSet_ = cms.string( "HLTIter0GroupedCkfTrajectoryBuilderIT" ) ),
    TransientInitialStateEstimatorParameters = cms.PSet( 
      propagatorAlongTISE = cms.string( "PropagatorWithMaterialParabolicMf" ),
      numberMeasurementsForFit = cms.int32( 4 ),
      propagatorOppositeTISE = cms.string( "PropagatorWithMaterialParabolicMfOpposite" )
    ),
    numHitsForSeedCleaner = cms.int32( 4 ),
    NavigationSchool = cms.string( "SimpleNavigationSchool" ),
    RedundantSeedCleaner = cms.string( "CachingSeedCleanerBySharedInput" ),
    TrajectoryCleaner = cms.string( "hltESPTrajectoryCleanerBySharedHits" ),
    maxNSeeds = cms.uint32( 100000 ),
    maxSeedsBeforeCleaning = cms.uint32( 1000 )
)
process.hltIter0L3MuonCtfWithMaterialTracks = cms.EDProducer( "TrackProducer",
    useSimpleMF = cms.bool( True ),
    SimpleMagneticField = cms.string( "ParabolicMf" ),
    src = cms.InputTag( "hltIter0L3MuonCkfTrackCandidates" ),
    clusterRemovalInfo = cms.InputTag( "" ),
    beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
    Fitter = cms.string( "hltESPFittingSmootherIT" ),
    useHitsSplitting = cms.bool( False ),
    alias = cms.untracked.string( "ctfWithMaterialTracks" ),
    TrajectoryInEvent = cms.bool( False ),
    TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
    AlgorithmName = cms.string( "hltIter0" ),
    Propagator = cms.string( "hltESPRungeKuttaTrackerPropagator" ),
    GeometricInnerState = cms.bool( True ),
    NavigationSchool = cms.string( "" ),
    MeasurementTracker = cms.string( "" ),
    MeasurementTrackerEvent = cms.InputTag( "hltSiStripClusters" )
)
process.hltIter0L3MuonTrackCutClassifier = cms.EDProducer( "TrackCutClassifier",
    src = cms.InputTag( "hltIter0L3MuonCtfWithMaterialTracks" ),
    beamspot = cms.InputTag( "hltOnlineBeamSpot" ),
    vertices = cms.InputTag( "hltTrimmedPixelVertices" ),
    ignoreVertices = cms.bool( False ),
    qualityCuts = cms.vdouble( -0.7, 0.1, 0.7 ),
    mva = cms.PSet( 
      minPixelHits = cms.vint32( 0, 0, 0 ),
      maxDzWrtBS = cms.vdouble( 3.40282346639E38, 24.0, 15.0 ),
      dr_par = cms.PSet( 
        d0err = cms.vdouble( 0.003, 0.003, 0.003 ),
        dr_par2 = cms.vdouble( 3.40282346639E38, 0.6, 0.6 ),
        dr_par1 = cms.vdouble( 3.40282346639E38, 0.8, 0.8 ),
        dr_exp = cms.vint32( 4, 4, 4 ),
        d0err_par = cms.vdouble( 0.001, 0.001, 0.001 )
      ),
      maxLostLayers = cms.vint32( 1, 1, 1 ),
      min3DLayers = cms.vint32( 0, 0, 0 ),
      dz_par = cms.PSet( 
        dz_par1 = cms.vdouble( 3.40282346639E38, 0.75, 0.75 ),
        dz_par2 = cms.vdouble( 3.40282346639E38, 0.5, 0.5 ),
        dz_exp = cms.vint32( 4, 4, 4 )
      ),
      minNVtxTrk = cms.int32( 3 ),
      maxDz = cms.vdouble( 0.5, 0.2, 3.40282346639E38 ),
      minNdof = cms.vdouble( 1.0E-5, 1.0E-5, 1.0E-5 ),
      maxChi2 = cms.vdouble( 9999.0, 25.0, 16.0 ),
      maxChi2n = cms.vdouble( 1.2, 1.0, 0.7 ),
      maxDr = cms.vdouble( 0.5, 0.03, 3.40282346639E38 ),
      minLayers = cms.vint32( 3, 3, 3 )
    )
)
process.hltIter0L3MuonTrackSelectionHighPurity = cms.EDProducer( "TrackCollectionFilterCloner",
    originalSource = cms.InputTag( "hltIter0L3MuonCtfWithMaterialTracks" ),
    originalMVAVals = cms.InputTag( 'hltIter0L3MuonTrackCutClassifier','MVAValues' ),
    originalQualVals = cms.InputTag( 'hltIter0L3MuonTrackCutClassifier','QualityMasks' ),
    minQuality = cms.string( "highPurity" ),
    copyExtras = cms.untracked.bool( True ),
    copyTrajectories = cms.untracked.bool( False )
)
process.hltMuonTkRelIsolationCut0p08Map = cms.EDProducer( "L3MuonCombinedRelativeIsolationProducer",
    UseRhoCorrectedCaloDeposits = cms.bool( False ),
    UseCaloIso = cms.bool( False ),
    CaloDepositsLabel = cms.InputTag( "notUsed" ),
    inputMuonCollection = cms.InputTag( "hltIterL3MuonCandidates" ),
    OutputMuIsoDeposits = cms.bool( True ),
    TrackPt_Min = cms.double( -1.0 ),
    printDebug = cms.bool( False ),
    CutsPSet = cms.PSet( 
      applyCutsORmaxNTracks = cms.bool( False ),
      maxNTracks = cms.int32( -1 ),
      Thresholds = cms.vdouble( 0.08 ),
      EtaBounds = cms.vdouble( 2.411 ),
      ComponentName = cms.string( "SimpleCuts" ),
      ConeSizes = cms.vdouble( 0.3 )
    ),
    TrkExtractorPSet = cms.PSet( 
      Diff_z = cms.double( 0.2 ),
      inputTrackCollection = cms.InputTag( "hltIter0L3MuonTrackSelectionHighPurity" ),
      Chi2Ndof_Max = cms.double( 1.0E64 ),
      BeamSpotLabel = cms.InputTag( "hltOnlineBeamSpot" ),
      DR_Veto = cms.double( 0.01 ),
      Pt_Min = cms.double( -1.0 ),
      VetoLeadingTrack = cms.bool( True ),
      DR_Max = cms.double( 0.3 ),
      DepositLabel = cms.untracked.string( "PXLS" ),
      PtVeto_Min = cms.double( 2.0 ),
      NHits_Min = cms.uint32( 0 ),
      PropagateTracksToRadius = cms.bool( True ),
      ReferenceRadius = cms.double( 6.0 ),
      Chi2Prob_Min = cms.double( -1.0 ),
      Diff_r = cms.double( 0.1 ),
      BeamlineOption = cms.string( "BeamSpotFromEvent" ),
      ComponentName = cms.string( "PixelTrackExtractor" ),
      DR_VetoPt = cms.double( 0.025 )
    ),
    CaloExtractorPSet = cms.PSet( 
      DR_Veto_H = cms.double( 0.1 ),
      Vertex_Constraint_Z = cms.bool( False ),
      DR_Veto_E = cms.double( 0.07 ),
      Weight_H = cms.double( 1.0 ),
      CaloTowerCollectionLabel = cms.InputTag( "hltTowerMakerForAll" ),
      DR_Max = cms.double( 0.3 ),
      DepositLabel = cms.untracked.string( "EcalPlusHcal" ),
      Vertex_Constraint_XY = cms.bool( False ),
      Threshold_H = cms.double( 0.5 ),
      Threshold_E = cms.double( 0.2 ),
      ComponentName = cms.string( "CaloExtractor" ),
      Weight_E = cms.double( 1.0 )
    )
)
process.hltRechitInRegionsECAL = cms.EDProducer( "HLTEcalRecHitInAllL1RegionsProducer",
    productLabels = cms.vstring( 'EcalRecHitsEB',
      'EcalRecHitsEE' ),
    recHitLabels = cms.VInputTag( 'hltEcalRecHit:EcalRecHitsEB','hltEcalRecHit:EcalRecHitsEE' ),
    l1InputRegions = cms.VPSet( 
      cms.PSet(  inputColl = cms.InputTag( 'hltGtStage2Digis','EGamma' ),
        regionEtaMargin = cms.double( 0.9 ),
        type = cms.string( "EGamma" ),
        minEt = cms.double( 5.0 ),
        regionPhiMargin = cms.double( 1.2 ),
        maxEt = cms.double( 999999.0 )
      ),
      cms.PSet(  inputColl = cms.InputTag( 'hltGtStage2Digis','Jet' ),
        regionEtaMargin = cms.double( 0.9 ),
        type = cms.string( "Jet" ),
        minEt = cms.double( 170.0 ),
        regionPhiMargin = cms.double( 1.2 ),
        maxEt = cms.double( 999999.0 )
      ),
      cms.PSet(  inputColl = cms.InputTag( 'hltGtStage2Digis','Tau' ),
        regionEtaMargin = cms.double( 0.9 ),
        type = cms.string( "Tau" ),
        minEt = cms.double( 100.0 ),
        regionPhiMargin = cms.double( 1.2 ),
        maxEt = cms.double( 999999.0 )
      )
    )
)
process.hltRechitInRegionsES = cms.EDProducer( "HLTEcalRecHitInAllL1RegionsProducer",
    productLabels = cms.vstring( 'EcalRecHitsES' ),
    recHitLabels = cms.VInputTag( 'hltEcalPreshowerRecHit:EcalRecHitsES' ),
    l1InputRegions = cms.VPSet( 
      cms.PSet(  inputColl = cms.InputTag( 'hltGtStage2Digis','EGamma' ),
        regionEtaMargin = cms.double( 0.9 ),
        type = cms.string( "EGamma" ),
        minEt = cms.double( 5.0 ),
        regionPhiMargin = cms.double( 1.2 ),
        maxEt = cms.double( 999999.0 )
      ),
      cms.PSet(  inputColl = cms.InputTag( 'hltGtStage2Digis','Jet' ),
        regionEtaMargin = cms.double( 0.9 ),
        type = cms.string( "Jet" ),
        minEt = cms.double( 170.0 ),
        regionPhiMargin = cms.double( 1.2 ),
        maxEt = cms.double( 999999.0 )
      ),
      cms.PSet(  inputColl = cms.InputTag( 'hltGtStage2Digis','Tau' ),
        regionEtaMargin = cms.double( 0.9 ),
        type = cms.string( "Tau" ),
        minEt = cms.double( 100.0 ),
        regionPhiMargin = cms.double( 1.2 ),
        maxEt = cms.double( 999999.0 )
      )
    )
)
process.hltParticleFlowRecHitECALL1Seeded = cms.EDProducer( "PFRecHitProducer",
    navigator = cms.PSet( 
      barrel = cms.PSet(  ),
      endcap = cms.PSet(  ),
      name = cms.string( "PFRecHitECALNavigator" )
    ),
    producers = cms.VPSet( 
      cms.PSet(  src = cms.InputTag( 'hltRechitInRegionsECAL','EcalRecHitsEB' ),
        srFlags = cms.InputTag( "" ),
        name = cms.string( "PFEBRecHitCreator" ),
        qualityTests = cms.VPSet( 
          cms.PSet(  name = cms.string( "PFRecHitQTestDBThreshold" ),
            applySelectionsToAllCrystals = cms.bool( True )
          ),
          cms.PSet(  topologicalCleaning = cms.bool( True ),
            skipTTRecoveredHits = cms.bool( True ),
            cleaningThreshold = cms.double( 2.0 ),
            name = cms.string( "PFRecHitQTestECAL" ),
            timingCleaning = cms.bool( True )
          )
        )
      ),
      cms.PSet(  src = cms.InputTag( 'hltRechitInRegionsECAL','EcalRecHitsEE' ),
        srFlags = cms.InputTag( "" ),
        name = cms.string( "PFEERecHitCreator" ),
        qualityTests = cms.VPSet( 
          cms.PSet(  name = cms.string( "PFRecHitQTestDBThreshold" ),
            applySelectionsToAllCrystals = cms.bool( True )
          ),
          cms.PSet(  topologicalCleaning = cms.bool( True ),
            skipTTRecoveredHits = cms.bool( True ),
            cleaningThreshold = cms.double( 2.0 ),
            name = cms.string( "PFRecHitQTestECAL" ),
            timingCleaning = cms.bool( True )
          )
        )
      )
    )
)
process.hltParticleFlowRecHitPSL1Seeded = cms.EDProducer( "PFRecHitProducer",
    navigator = cms.PSet(  name = cms.string( "PFRecHitPreshowerNavigator" ) ),
    producers = cms.VPSet( 
      cms.PSet(  src = cms.InputTag( 'hltRechitInRegionsES','EcalRecHitsES' ),
        name = cms.string( "PFPSRecHitCreator" ),
        qualityTests = cms.VPSet( 
          cms.PSet(  threshold = cms.double( 7.0E-6 ),
            name = cms.string( "PFRecHitQTestThreshold" )
          )
        )
      )
    )
)
process.hltParticleFlowClusterPSL1Seeded = cms.EDProducer( "PFClusterProducer",
    recHitsSource = cms.InputTag( "hltParticleFlowRecHitPSL1Seeded" ),
    recHitCleaners = cms.VPSet( 
    ),
    seedCleaners = cms.VPSet( 
    ),
    seedFinder = cms.PSet( 
      thresholdsByDetector = cms.VPSet( 
        cms.PSet(  seedingThresholdPt = cms.double( 0.0 ),
          seedingThreshold = cms.double( 1.2E-4 ),
          detector = cms.string( "PS1" )
        ),
        cms.PSet(  seedingThresholdPt = cms.double( 0.0 ),
          seedingThreshold = cms.double( 1.2E-4 ),
          detector = cms.string( "PS2" )
        )
      ),
      algoName = cms.string( "LocalMaximumSeedFinder" ),
      nNeighbours = cms.int32( 4 )
    ),
    initialClusteringStep = cms.PSet( 
      thresholdsByDetector = cms.VPSet( 
        cms.PSet(  gatheringThreshold = cms.double( 6.0E-5 ),
          gatheringThresholdPt = cms.double( 0.0 ),
          detector = cms.string( "PS1" )
        ),
        cms.PSet(  gatheringThreshold = cms.double( 6.0E-5 ),
          gatheringThresholdPt = cms.double( 0.0 ),
          detector = cms.string( "PS2" )
        )
      ),
      algoName = cms.string( "Basic2DGenericTopoClusterizer" ),
      useCornerCells = cms.bool( False )
    ),
    pfClusterBuilder = cms.PSet( 
      minFracTot = cms.double( 1.0E-20 ),
      stoppingTolerance = cms.double( 1.0E-8 ),
      positionCalc = cms.PSet( 
        minAllowedNormalization = cms.double( 1.0E-9 ),
        posCalcNCrystals = cms.int32( -1 ),
        algoName = cms.string( "Basic2DGenericPFlowPositionCalc" ),
        logWeightDenominator = cms.double( 6.0E-5 ),
        minFractionInCalc = cms.double( 1.0E-9 )
      ),
      maxIterations = cms.uint32( 50 ),
      algoName = cms.string( "Basic2DGenericPFlowClusterizer" ),
      recHitEnergyNorms = cms.VPSet( 
        cms.PSet(  recHitEnergyNorm = cms.double( 6.0E-5 ),
          detector = cms.string( "PS1" )
        ),
        cms.PSet(  recHitEnergyNorm = cms.double( 6.0E-5 ),
          detector = cms.string( "PS2" )
        )
      ),
      showerSigma = cms.double( 0.3 ),
      minFractionToKeep = cms.double( 1.0E-7 ),
      excludeOtherSeeds = cms.bool( True )
    ),
    positionReCalc = cms.PSet(  ),
    energyCorrector = cms.PSet(  )
)
process.hltParticleFlowClusterECALUncorrectedL1Seeded = cms.EDProducer( "PFClusterProducer",
    recHitsSource = cms.InputTag( "hltParticleFlowRecHitECALL1Seeded" ),
    recHitCleaners = cms.VPSet( 
    ),
    seedCleaners = cms.VPSet( 
    ),
    seedFinder = cms.PSet( 
      thresholdsByDetector = cms.VPSet( 
        cms.PSet(  seedingThresholdPt = cms.double( 0.15 ),
          seedingThreshold = cms.double( 0.6 ),
          detector = cms.string( "ECAL_ENDCAP" )
        ),
        cms.PSet(  seedingThresholdPt = cms.double( 0.0 ),
          seedingThreshold = cms.double( 0.23 ),
          detector = cms.string( "ECAL_BARREL" )
        )
      ),
      algoName = cms.string( "LocalMaximumSeedFinder" ),
      nNeighbours = cms.int32( 8 )
    ),
    initialClusteringStep = cms.PSet( 
      thresholdsByDetector = cms.VPSet( 
        cms.PSet(  gatheringThreshold = cms.double( 0.08 ),
          gatheringThresholdPt = cms.double( 0.0 ),
          detector = cms.string( "ECAL_BARREL" )
        ),
        cms.PSet(  gatheringThreshold = cms.double( 0.3 ),
          gatheringThresholdPt = cms.double( 0.0 ),
          detector = cms.string( "ECAL_ENDCAP" )
        )
      ),
      algoName = cms.string( "Basic2DGenericTopoClusterizer" ),
      useCornerCells = cms.bool( True )
    ),
    pfClusterBuilder = cms.PSet( 
      minFracTot = cms.double( 1.0E-20 ),
      stoppingTolerance = cms.double( 1.0E-8 ),
      positionCalc = cms.PSet( 
        minAllowedNormalization = cms.double( 1.0E-9 ),
        posCalcNCrystals = cms.int32( 9 ),
        algoName = cms.string( "Basic2DGenericPFlowPositionCalc" ),
        logWeightDenominator = cms.double( 0.08 ),
        minFractionInCalc = cms.double( 1.0E-9 ),
        timeResolutionCalcBarrel = cms.PSet( 
          corrTermLowE = cms.double( 0.0510871 ),
          threshLowE = cms.double( 0.5 ),
          noiseTerm = cms.double( 1.10889 ),
          constantTermLowE = cms.double( 0.0 ),
          noiseTermLowE = cms.double( 1.31883 ),
          threshHighE = cms.double( 5.0 ),
          constantTerm = cms.double( 0.428192 )
        ),
        timeResolutionCalcEndcap = cms.PSet( 
          corrTermLowE = cms.double( 0.0 ),
          threshLowE = cms.double( 1.0 ),
          noiseTerm = cms.double( 5.72489999999 ),
          constantTermLowE = cms.double( 0.0 ),
          noiseTermLowE = cms.double( 6.92683000001 ),
          threshHighE = cms.double( 10.0 ),
          constantTerm = cms.double( 0.0 )
        )
      ),
      maxIterations = cms.uint32( 50 ),
      positionCalcForConvergence = cms.PSet( 
        minAllowedNormalization = cms.double( 0.0 ),
        T0_ES = cms.double( 1.2 ),
        algoName = cms.string( "ECAL2DPositionCalcWithDepthCorr" ),
        T0_EE = cms.double( 3.1 ),
        T0_EB = cms.double( 7.4 ),
        X0 = cms.double( 0.89 ),
        minFractionInCalc = cms.double( 0.0 ),
        W0 = cms.double( 4.2 )
      ),
      allCellsPositionCalc = cms.PSet( 
        minAllowedNormalization = cms.double( 1.0E-9 ),
        posCalcNCrystals = cms.int32( -1 ),
        algoName = cms.string( "Basic2DGenericPFlowPositionCalc" ),
        logWeightDenominator = cms.double( 0.08 ),
        minFractionInCalc = cms.double( 1.0E-9 ),
        timeResolutionCalcBarrel = cms.PSet( 
          corrTermLowE = cms.double( 0.0510871 ),
          threshLowE = cms.double( 0.5 ),
          noiseTerm = cms.double( 1.10889 ),
          constantTermLowE = cms.double( 0.0 ),
          noiseTermLowE = cms.double( 1.31883 ),
          threshHighE = cms.double( 5.0 ),
          constantTerm = cms.double( 0.428192 )
        ),
        timeResolutionCalcEndcap = cms.PSet( 
          corrTermLowE = cms.double( 0.0 ),
          threshLowE = cms.double( 1.0 ),
          noiseTerm = cms.double( 5.72489999999 ),
          constantTermLowE = cms.double( 0.0 ),
          noiseTermLowE = cms.double( 6.92683000001 ),
          threshHighE = cms.double( 10.0 ),
          constantTerm = cms.double( 0.0 )
        )
      ),
      algoName = cms.string( "Basic2DGenericPFlowClusterizer" ),
      recHitEnergyNorms = cms.VPSet( 
        cms.PSet(  recHitEnergyNorm = cms.double( 0.08 ),
          detector = cms.string( "ECAL_BARREL" )
        ),
        cms.PSet(  recHitEnergyNorm = cms.double( 0.3 ),
          detector = cms.string( "ECAL_ENDCAP" )
        )
      ),
      showerSigma = cms.double( 1.5 ),
      minFractionToKeep = cms.double( 1.0E-7 ),
      excludeOtherSeeds = cms.bool( True )
    ),
    positionReCalc = cms.PSet( 
      minAllowedNormalization = cms.double( 0.0 ),
      T0_ES = cms.double( 1.2 ),
      algoName = cms.string( "ECAL2DPositionCalcWithDepthCorr" ),
      T0_EE = cms.double( 3.1 ),
      T0_EB = cms.double( 7.4 ),
      X0 = cms.double( 0.89 ),
      minFractionInCalc = cms.double( 0.0 ),
      W0 = cms.double( 4.2 )
    ),
    energyCorrector = cms.PSet(  )
)
process.hltParticleFlowClusterECALL1Seeded = cms.EDProducer( "CorrectedECALPFClusterProducer",
    minimumPSEnergy = cms.double( 0.0 ),
    skipPS = cms.bool( False ),
    inputPS = cms.InputTag( "hltParticleFlowClusterPSL1Seeded" ),
    energyCorrector = cms.PSet( 
      applyCrackCorrections = cms.bool( False ),
      srfAwareCorrection = cms.bool( True ),
      applyMVACorrections = cms.bool( True ),
      maxPtForMVAEvaluation = cms.double( 300.0 ),
      recHitsEBLabel = cms.InputTag( 'hltEcalRecHit','EcalRecHitsEB' ),
      recHitsEELabel = cms.InputTag( 'hltEcalRecHit','EcalRecHitsEE' ),
      ebSrFlagLabel = cms.InputTag( "hltEcalDigis" ),
      eeSrFlagLabel = cms.InputTag( "hltEcalDigis" )
    ),
    inputECAL = cms.InputTag( "hltParticleFlowClusterECALUncorrectedL1Seeded" )
)
process.hltParticleFlowSuperClusterECALL1Seeded = cms.EDProducer( "PFECALSuperClusterProducer",
    PFSuperClusterCollectionEndcap = cms.string( "hltParticleFlowSuperClusterECALEndcap" ),
    doSatelliteClusterMerge = cms.bool( False ),
    thresh_PFClusterBarrel = cms.double( 0.5 ),
    PFBasicClusterCollectionBarrel = cms.string( "hltParticleFlowBasicClusterECALBarrel" ),
    useRegression = cms.bool( True ),
    satelliteMajorityFraction = cms.double( 0.5 ),
    thresh_PFClusterEndcap = cms.double( 0.5 ),
    ESAssociation = cms.InputTag( "hltParticleFlowClusterECALL1Seeded" ),
    PFBasicClusterCollectionPreshower = cms.string( "hltParticleFlowBasicClusterECALPreshower" ),
    use_preshower = cms.bool( True ),
    verbose = cms.untracked.bool( False ),
    thresh_SCEt = cms.double( 4.0 ),
    etawidth_SuperClusterEndcap = cms.double( 0.04 ),
    phiwidth_SuperClusterEndcap = cms.double( 0.6 ),
    useDynamicDPhiWindow = cms.bool( True ),
    PFSuperClusterCollectionBarrel = cms.string( "hltParticleFlowSuperClusterECALBarrel" ),
    regressionConfig = cms.PSet( 
      uncertaintyKeyEB = cms.string( "pfscecal_EBUncertainty_online" ),
      ecalRecHitsEE = cms.InputTag( 'hltEcalRecHit','EcalRecHitsEE' ),
      ecalRecHitsEB = cms.InputTag( 'hltEcalRecHit','EcalRecHitsEB' ),
      regressionKeyEE = cms.string( "pfscecal_EECorrection_online" ),
      regressionKeyEB = cms.string( "pfscecal_EBCorrection_online" ),
      uncertaintyKeyEE = cms.string( "pfscecal_EEUncertainty_online" ),
      isHLT = cms.bool( True )
    ),
    applyCrackCorrections = cms.bool( False ),
    satelliteClusterSeedThreshold = cms.double( 50.0 ),
    etawidth_SuperClusterBarrel = cms.double( 0.04 ),
    PFBasicClusterCollectionEndcap = cms.string( "hltParticleFlowBasicClusterECALEndcap" ),
    PFClusters = cms.InputTag( "hltParticleFlowClusterECALL1Seeded" ),
    thresh_PFClusterSeedBarrel = cms.double( 1.0 ),
    ClusteringType = cms.string( "Mustache" ),
    EnergyWeight = cms.string( "Raw" ),
    BeamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
    thresh_PFClusterSeedEndcap = cms.double( 1.0 ),
    phiwidth_SuperClusterBarrel = cms.double( 0.6 ),
    thresh_PFClusterES = cms.double( 0.5 ),
    seedThresholdIsET = cms.bool( True ),
    isOOTCollection = cms.bool( False ),
    barrelRecHits = cms.InputTag( 'ecalRecHit','EcalRecHitsEE' ),
    endcapRecHits = cms.InputTag( 'ecalRecHit','EcalRecHitsEB' ),
    PFSuperClusterCollectionEndcapWithPreshower = cms.string( "hltParticleFlowSuperClusterECALEndcapWithPreshower" ),
    dropUnseedable = cms.bool( False )
)
process.hltEgammaCandidates = cms.EDProducer( "EgammaHLTRecoEcalCandidateProducers",
    scHybridBarrelProducer = cms.InputTag( 'hltParticleFlowSuperClusterECALL1Seeded','hltParticleFlowSuperClusterECALBarrel' ),
    scIslandEndcapProducer = cms.InputTag( 'hltParticleFlowSuperClusterECALL1Seeded','hltParticleFlowSuperClusterECALEndcapWithPreshower' ),
    recoEcalCandidateCollection = cms.string( "" )
)
process.hltEgammaClusterShape = cms.EDProducer( "EgammaHLTClusterShapeProducer",
    recoEcalCandidateProducer = cms.InputTag( "hltEgammaCandidates" ),
    ecalRechitEB = cms.InputTag( 'hltRechitInRegionsECAL','EcalRecHitsEB' ),
    ecalRechitEE = cms.InputTag( 'hltRechitInRegionsECAL','EcalRecHitsEE' ),
    isIeta = cms.bool( True ),
    multThresEB = cms.double( 1.0 ),
    multThresEE = cms.double( 1.25 )
)
process.hltFixedGridRhoFastjetAllCaloForMuons = cms.EDProducer( "FixedGridRhoProducerFastjet",
    pfCandidatesTag = cms.InputTag( "hltTowerMakerForAll" ),
    maxRapidity = cms.double( 2.5 ),
    gridSpacing = cms.double( 0.55 )
)
process.hltEgammaHoverE = cms.EDProducer( "EgammaHLTBcHcalIsolationProducersRegional",
    recoEcalCandidateProducer = cms.InputTag( "hltEgammaCandidates" ),
    caloTowerProducer = cms.InputTag( "hltTowerMakerForAll" ),
    rhoProducer = cms.InputTag( "hltFixedGridRhoFastjetAllCaloForMuons" ),
    doRhoCorrection = cms.bool( False ),
    rhoMax = cms.double( 9.9999999E7 ),
    rhoScale = cms.double( 1.0 ),
    etMin = cms.double( 0.0 ),
    innerCone = cms.double( 0.0 ),
    outerCone = cms.double( 0.14 ),
    depth = cms.int32( -1 ),
    doEtSum = cms.bool( False ),
    useSingleTower = cms.bool( False ),
    effectiveAreas = cms.vdouble( 0.105, 0.17 ),
    absEtaLowEdges = cms.vdouble( 0.0, 1.479 )
)
process.hltEgammaEcalPFClusterIso = cms.EDProducer( "EgammaHLTEcalPFClusterIsolationProducer",
    recoEcalCandidateProducer = cms.InputTag( "hltEgammaCandidates" ),
    pfClusterProducer = cms.InputTag( "hltParticleFlowClusterECALL1Seeded" ),
    rhoProducer = cms.InputTag( "hltFixedGridRhoFastjetAllCaloForMuons" ),
    doRhoCorrection = cms.bool( False ),
    rhoMax = cms.double( 9.9999999E7 ),
    rhoScale = cms.double( 1.0 ),
    drMax = cms.double( 0.3 ),
    drVetoBarrel = cms.double( 0.0 ),
    drVetoEndcap = cms.double( 0.0 ),
    etaStripBarrel = cms.double( 0.0 ),
    etaStripEndcap = cms.double( 0.0 ),
    energyBarrel = cms.double( 0.0 ),
    energyEndcap = cms.double( 0.0 ),
    effectiveAreas = cms.vdouble( 0.29, 0.21 ),
    absEtaLowEdges = cms.vdouble( 0.0, 1.479 )
)
process.hltEgammaHcalPFClusterIso = cms.EDProducer( "EgammaHLTHcalPFClusterIsolationProducer",
    recoEcalCandidateProducer = cms.InputTag( "hltEgammaCandidates" ),
    pfClusterProducerHCAL = cms.InputTag( "hltParticleFlowClusterHCAL" ),
    useHF = cms.bool( False ),
    pfClusterProducerHFEM = cms.InputTag( "" ),
    pfClusterProducerHFHAD = cms.InputTag( "" ),
    rhoProducer = cms.InputTag( "hltFixedGridRhoFastjetAllCaloForMuons" ),
    doRhoCorrection = cms.bool( False ),
    rhoMax = cms.double( 9.9999999E7 ),
    rhoScale = cms.double( 1.0 ),
    drMax = cms.double( 0.3 ),
    drVetoBarrel = cms.double( 0.0 ),
    drVetoEndcap = cms.double( 0.0 ),
    etaStripBarrel = cms.double( 0.0 ),
    etaStripEndcap = cms.double( 0.0 ),
    energyBarrel = cms.double( 0.0 ),
    energyEndcap = cms.double( 0.0 ),
    useEt = cms.bool( True ),
    effectiveAreas = cms.vdouble( 0.2, 0.25 ),
    absEtaLowEdges = cms.vdouble( 0.0, 1.479 )
)
process.hltPixelLayerPairs = cms.EDProducer( "SeedingLayersEDProducer",
    layerList = cms.vstring( 'BPix1+BPix2',
      'BPix1+BPix3',
      'BPix1+BPix4',
      'BPix2+BPix3',
      'BPix2+BPix4',
      'BPix3+BPix4',
      'FPix1_pos+FPix2_pos',
      'FPix1_pos+FPix3_pos',
      'FPix2_pos+FPix3_pos',
      'BPix1+FPix1_pos',
      'BPix1+FPix2_pos',
      'BPix1+FPix3_pos',
      'BPix2+FPix1_pos',
      'BPix2+FPix2_pos',
      'BPix2+FPix3_pos',
      'BPix3+FPix1_pos',
      'BPix3+FPix2_pos',
      'BPix3+FPix3_pos',
      'BPix4+FPix1_pos',
      'BPix4+FPix2_pos',
      'BPix4+FPix3_pos',
      'FPix1_neg+FPix2_neg',
      'FPix1_neg+FPix3_neg',
      'FPix2_neg+FPix3_neg',
      'BPix1+FPix1_neg',
      'BPix1+FPix2_neg',
      'BPix1+FPix3_neg',
      'BPix2+FPix1_neg',
      'BPix2+FPix2_neg',
      'BPix2+FPix3_neg',
      'BPix3+FPix1_neg',
      'BPix3+FPix2_neg',
      'BPix3+FPix3_neg',
      'BPix4+FPix1_neg',
      'BPix4+FPix2_neg',
      'BPix4+FPix3_neg' ),
    BPix = cms.PSet( 
      hitErrorRPhi = cms.double( 0.0027 ),
      TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
      useErrorsFromParam = cms.bool( True ),
      hitErrorRZ = cms.double( 0.006 ),
      HitProducer = cms.string( "hltSiPixelRecHits" )
    ),
    FPix = cms.PSet( 
      hitErrorRPhi = cms.double( 0.0051 ),
      TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
      useErrorsFromParam = cms.bool( True ),
      hitErrorRZ = cms.double( 0.0036 ),
      HitProducer = cms.string( "hltSiPixelRecHits" )
    ),
    TIB = cms.PSet(  ),
    TID = cms.PSet(  ),
    TOB = cms.PSet(  ),
    TEC = cms.PSet(  ),
    MTIB = cms.PSet(  ),
    MTID = cms.PSet(  ),
    MTOB = cms.PSet(  ),
    MTEC = cms.PSet(  )
)
process.hltPixelLayerTriplets = cms.EDProducer( "SeedingLayersEDProducer",
    layerList = cms.vstring( 'BPix1+BPix2+BPix3',
      'BPix2+BPix3+BPix4',
      'BPix1+BPix3+BPix4',
      'BPix1+BPix2+BPix4',
      'BPix2+BPix3+FPix1_pos',
      'BPix2+BPix3+FPix1_neg',
      'BPix1+BPix2+FPix1_pos',
      'BPix1+BPix2+FPix1_neg',
      'BPix2+FPix1_pos+FPix2_pos',
      'BPix2+FPix1_neg+FPix2_neg',
      'BPix1+FPix1_pos+FPix2_pos',
      'BPix1+FPix1_neg+FPix2_neg',
      'FPix1_pos+FPix2_pos+FPix3_pos',
      'FPix1_neg+FPix2_neg+FPix3_neg',
      'BPix1+BPix3+FPix1_pos',
      'BPix1+BPix2+FPix2_pos',
      'BPix1+BPix3+FPix1_neg',
      'BPix1+BPix2+FPix2_neg',
      'BPix1+FPix2_neg+FPix3_neg',
      'BPix1+FPix1_neg+FPix3_neg',
      'BPix1+FPix2_pos+FPix3_pos',
      'BPix1+FPix1_pos+FPix3_pos' ),
    BPix = cms.PSet( 
      hitErrorRPhi = cms.double( 0.0027 ),
      TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
      useErrorsFromParam = cms.bool( True ),
      hitErrorRZ = cms.double( 0.006 ),
      HitProducer = cms.string( "hltSiPixelRecHits" )
    ),
    FPix = cms.PSet( 
      hitErrorRPhi = cms.double( 0.0051 ),
      TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
      useErrorsFromParam = cms.bool( True ),
      hitErrorRZ = cms.double( 0.0036 ),
      HitProducer = cms.string( "hltSiPixelRecHits" )
    ),
    TIB = cms.PSet(  ),
    TID = cms.PSet(  ),
    TOB = cms.PSet(  ),
    TEC = cms.PSet(  ),
    MTIB = cms.PSet(  ),
    MTID = cms.PSet(  ),
    MTOB = cms.PSet(  ),
    MTEC = cms.PSet(  )
)
process.hltEgammaSuperClustersToPixelMatch = cms.EDProducer( "EgammaHLTFilteredSuperClusterProducer",
    cands = cms.InputTag( "hltEgammaCandidates" ),
    minEtCutEB = cms.double( 0.0 ),
    minEtCutEE = cms.double( 0.0 ),
    cuts = cms.VPSet( 
      cms.PSet(  endcapCut = cms.PSet( 
  useEt = cms.bool( False ),
  cutOverE = cms.double( 0.2 )
),
        var = cms.InputTag( "hltEgammaHoverE" ),
        barrelCut = cms.PSet( 
          useEt = cms.bool( False ),
          cutOverE = cms.double( 0.2 )
        )
      )
    )
)
process.hltEleSeedsTrackingRegions = cms.EDProducer( "TrackingRegionsFromSuperClustersEDProducer",
    RegionPSet = cms.PSet( 
      minBSDeltaZ = cms.double( 0.0 ),
      useZInVertex = cms.bool( False ),
      vertices = cms.InputTag( "" ),
      beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
      useZInBeamspot = cms.bool( False ),
      ptMin = cms.double( 1.5 ),
      deltaEtaRegion = cms.double( 0.1 ),
      nrSigmaForBSDeltaZ = cms.double( 4.0 ),
      originHalfLength = cms.double( 12.5 ),
      measurementTrackerEvent = cms.InputTag( "" ),
      originRadius = cms.double( 0.2 ),
      precise = cms.bool( True ),
      superClusters = cms.VInputTag( 'hltEgammaSuperClustersToPixelMatch' ),
      whereToUseMeasTracker = cms.string( "kNever" ),
      deltaPhiRegion = cms.double( 0.4 ),
      defaultZ = cms.double( 0.0 )
    )
)
process.hltElePixelHitDoublets = cms.EDProducer( "HitPairEDProducer",
    seedingLayers = cms.InputTag( "hltPixelLayerPairs" ),
    trackingRegions = cms.InputTag( "hltEleSeedsTrackingRegions" ),
    trackingRegionsSeedingLayers = cms.InputTag( "" ),
    clusterCheck = cms.InputTag( "" ),
    produceSeedingHitSets = cms.bool( True ),
    produceIntermediateHitDoublets = cms.bool( True ),
    maxElement = cms.uint32( 0 ),
    maxElementTotal = cms.uint32( 50000000 ),
    layerPairs = cms.vuint32( 0 )
)
process.hltElePixelHitDoubletsForTriplets = cms.EDProducer( "HitPairEDProducer",
    seedingLayers = cms.InputTag( "hltPixelLayerTriplets" ),
    trackingRegions = cms.InputTag( "hltEleSeedsTrackingRegions" ),
    trackingRegionsSeedingLayers = cms.InputTag( "" ),
    clusterCheck = cms.InputTag( "" ),
    produceSeedingHitSets = cms.bool( True ),
    produceIntermediateHitDoublets = cms.bool( True ),
    maxElement = cms.uint32( 0 ),
    maxElementTotal = cms.uint32( 50000000 ),
    layerPairs = cms.vuint32( 0, 1 )
)
process.hltElePixelHitTriplets = cms.EDProducer( "CAHitTripletEDProducer",
    doublets = cms.InputTag( "hltElePixelHitDoubletsForTriplets" ),
    extraHitRPhitolerance = cms.double( 0.032 ),
    useBendingCorrection = cms.bool( True ),
    CAThetaCut = cms.double( 0.004 ),
    CAPhiCut = cms.double( 0.1 ),
    CAThetaCut_byTriplets = cms.VPSet( 
      cms.PSet(  seedingLayers = cms.string( "" ),
        cut = cms.double( -1.0 )
      )
    ),
    CAPhiCut_byTriplets = cms.VPSet( 
      cms.PSet(  seedingLayers = cms.string( "" ),
        cut = cms.double( -1.0 )
      )
    ),
    CAHardPtCut = cms.double( 0.3 ),
    maxChi2 = cms.PSet( 
      value2 = cms.double( 6.0 ),
      value1 = cms.double( 100.0 ),
      pt1 = cms.double( 0.8 ),
      enabled = cms.bool( True ),
      pt2 = cms.double( 8.0 )
    ),
    SeedComparitorPSet = cms.PSet(  ComponentName = cms.string( "none" ) )
)
process.hltElePixelSeedsDoublets = cms.EDProducer( "SeedCreatorFromRegionConsecutiveHitsEDProducer",
    seedingHitSets = cms.InputTag( "hltElePixelHitDoublets" ),
    propagator = cms.string( "PropagatorWithMaterialParabolicMf" ),
    SeedMomentumForBOFF = cms.double( 5.0 ),
    OriginTransverseErrorMultiplier = cms.double( 1.0 ),
    MinOneOverPtError = cms.double( 1.0 ),
    TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
    magneticField = cms.string( "ParabolicMf" ),
    forceKinematicWithRegionDirection = cms.bool( False ),
    SeedComparitorPSet = cms.PSet(  ComponentName = cms.string( "none" ) )
)
process.hltElePixelSeedsTriplets = cms.EDProducer( "SeedCreatorFromRegionConsecutiveHitsEDProducer",
    seedingHitSets = cms.InputTag( "hltElePixelHitTriplets" ),
    propagator = cms.string( "PropagatorWithMaterialParabolicMf" ),
    SeedMomentumForBOFF = cms.double( 5.0 ),
    OriginTransverseErrorMultiplier = cms.double( 1.0 ),
    MinOneOverPtError = cms.double( 1.0 ),
    TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
    magneticField = cms.string( "ParabolicMf" ),
    forceKinematicWithRegionDirection = cms.bool( False ),
    SeedComparitorPSet = cms.PSet(  ComponentName = cms.string( "none" ) )
)
process.hltElePixelSeedsCombined = cms.EDProducer( "SeedCombiner",
    seedCollections = cms.VInputTag( 'hltElePixelSeedsDoublets','hltElePixelSeedsTriplets' )
)
process.hltEgammaElectronPixelSeeds = cms.EDProducer( "ElectronNHitSeedProducer",
    initialSeeds = cms.InputTag( "hltElePixelSeedsCombined" ),
    vertices = cms.InputTag( "" ),
    beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
    measTkEvt = cms.InputTag( "hltSiStripClusters" ),
    superClusters = cms.VInputTag( 'hltEgammaSuperClustersToPixelMatch' ),
    matcherConfig = cms.PSet( 
      useRecoVertex = cms.bool( False ),
      minNrHits = cms.vuint32( 2, 3 ),
      matchingCuts = cms.VPSet( 
        cms.PSet(  dPhiMaxHighEt = cms.vdouble( 0.05 ),
          version = cms.int32( 2 ),
          dRZMaxHighEt = cms.vdouble( 9999.0 ),
          dRZMaxLowEtGrad = cms.vdouble( 0.0 ),
          dPhiMaxLowEtGrad = cms.vdouble( -0.002 ),
          dPhiMaxHighEtThres = cms.vdouble( 20.0 ),
          dRZMaxHighEtThres = cms.vdouble( 0.0 )
        ),
        cms.PSet(  etaBins = cms.vdouble(  ),
          dPhiMaxHighEt = cms.vdouble( 0.003 ),
          version = cms.int32( 2 ),
          dRZMaxHighEt = cms.vdouble( 0.05 ),
          dRZMaxLowEtGrad = cms.vdouble( -0.002 ),
          dPhiMaxLowEtGrad = cms.vdouble( 0.0 ),
          dPhiMaxHighEtThres = cms.vdouble( 0.0 ),
          dRZMaxHighEtThres = cms.vdouble( 30.0 )
        ),
        cms.PSet(  etaBins = cms.vdouble(  ),
          dPhiMaxHighEt = cms.vdouble( 0.003 ),
          version = cms.int32( 2 ),
          dRZMaxHighEt = cms.vdouble( 0.05 ),
          dRZMaxLowEtGrad = cms.vdouble( -0.002 ),
          dPhiMaxLowEtGrad = cms.vdouble( 0.0 ),
          dPhiMaxHighEtThres = cms.vdouble( 0.0 ),
          dRZMaxHighEtThres = cms.vdouble( 30.0 )
        )
      ),
      minNrHitsValidLayerBins = cms.vint32( 4 ),
      detLayerGeom = cms.ESInputTag( "","hltESPGlobalDetLayerGeometry" ),
      navSchool = cms.ESInputTag( "","SimpleNavigationSchool" ),
      paramMagField = cms.ESInputTag( "","ParabolicMf" )
    )
)
process.hltEgammaPixelMatchVars = cms.EDProducer( "EgammaHLTPixelMatchVarProducer",
    recoEcalCandidateProducer = cms.InputTag( "hltEgammaCandidates" ),
    pixelSeedsProducer = cms.InputTag( "hltEgammaElectronPixelSeeds" ),
    dPhi1SParams = cms.PSet(  bins = cms.VPSet( 
  cms.PSet(  yMin = cms.int32( 1 ),
    binType = cms.string( "AbsEtaClus" ),
    funcParams = cms.vdouble( 0.00112, 7.52E-4, -0.00122, 0.00109 ),
    xMin = cms.double( 0.0 ),
    yMax = cms.int32( 1 ),
    xMax = cms.double( 1.5 ),
    funcType = cms.string( "TF1:=pol3" )
  ),
  cms.PSet(  yMin = cms.int32( 2 ),
    binType = cms.string( "AbsEtaClus" ),
    funcParams = cms.vdouble( 0.00222, 1.96E-4, -2.03E-4, 4.47E-4 ),
    xMin = cms.double( 0.0 ),
    yMax = cms.int32( 2 ),
    xMax = cms.double( 1.5 ),
    funcType = cms.string( "TF1:=pol3" )
  ),
  cms.PSet(  yMin = cms.int32( 3 ),
    binType = cms.string( "AbsEtaClus" ),
    funcParams = cms.vdouble( 0.00236, 6.91E-4, 1.99E-4, 4.16E-4 ),
    xMin = cms.double( 0.0 ),
    yMax = cms.int32( 99999 ),
    xMax = cms.double( 1.5 ),
    funcType = cms.string( "TF1:=pol3" )
  ),
  cms.PSet(  yMin = cms.int32( 1 ),
    binType = cms.string( "AbsEtaClus" ),
    funcParams = cms.vdouble( 0.00823, -0.0029 ),
    xMin = cms.double( 1.5 ),
    yMax = cms.int32( 1 ),
    xMax = cms.double( 2.0 ),
    funcType = cms.string( "TF1:=pol1" )
  ),
  cms.PSet(  yMin = cms.int32( 1 ),
    binType = cms.string( "AbsEtaClus" ),
    funcParams = cms.vdouble( 0.00282 ),
    xMin = cms.double( 2.0 ),
    yMax = cms.int32( 1 ),
    xMax = cms.double( 3.0 ),
    funcType = cms.string( "TF1:=pol0" )
  ),
  cms.PSet(  yMin = cms.int32( 2 ),
    binType = cms.string( "AbsEtaClus" ),
    funcParams = cms.vdouble( 0.010838, -0.00345 ),
    xMin = cms.double( 1.5 ),
    yMax = cms.int32( 2 ),
    xMax = cms.double( 2.0 ),
    funcType = cms.string( "TF1:=pol1" )
  ),
  cms.PSet(  yMin = cms.int32( 2 ),
    binType = cms.string( "AbsEtaClus" ),
    funcParams = cms.vdouble( 0.0043 ),
    xMin = cms.double( 2.0 ),
    yMax = cms.int32( 2 ),
    xMax = cms.double( 3.0 ),
    funcType = cms.string( "TF1:=pol0" )
  ),
  cms.PSet(  yMin = cms.int32( 3 ),
    binType = cms.string( "AbsEtaClus" ),
    funcParams = cms.vdouble( 0.0208, -0.0125, 0.00231 ),
    xMin = cms.double( 1.5 ),
    yMax = cms.int32( 99999 ),
    xMax = cms.double( 3.0 ),
    funcType = cms.string( "TF1:=pol2" )
  )
) ),
    dPhi2SParams = cms.PSet(  bins = cms.VPSet( 
  cms.PSet(  yMin = cms.int32( 1 ),
    binType = cms.string( "AbsEtaClus" ),
    funcParams = cms.vdouble( 1.3E-4 ),
    xMin = cms.double( 0.0 ),
    yMax = cms.int32( 99999 ),
    xMax = cms.double( 1.6 ),
    funcType = cms.string( "TF1:=pol0" )
  ),
  cms.PSet(  yMin = cms.int32( 1 ),
    binType = cms.string( "AbsEtaClus" ),
    funcParams = cms.vdouble( 4.5E-4, -1.99E-4 ),
    xMin = cms.double( 1.5 ),
    yMax = cms.int32( 99999 ),
    xMax = cms.double( 1.9 ),
    funcType = cms.string( "TF1:=pol1" )
  ),
  cms.PSet(  yMin = cms.int32( 1 ),
    binType = cms.string( "AbsEtaClus" ),
    funcParams = cms.vdouble( 7.94E-5 ),
    xMin = cms.double( 1.9 ),
    yMax = cms.int32( 99999 ),
    xMax = cms.double( 3.0 ),
    funcType = cms.string( "TF1:=pol0" )
  )
) ),
    dRZ2SParams = cms.PSet(  bins = cms.VPSet( 
  cms.PSet(  yMin = cms.int32( 1 ),
    binType = cms.string( "AbsEtaClus" ),
    funcParams = cms.vdouble( 0.00299, 2.99E-4, -4.13E-6, 0.00191 ),
    xMin = cms.double( 0.0 ),
    yMax = cms.int32( 99999 ),
    xMax = cms.double( 1.5 ),
    funcType = cms.string( "TF1:=pol3" )
  ),
  cms.PSet(  yMin = cms.int32( 1 ),
    binType = cms.string( "AbsEtaClus" ),
    funcParams = cms.vdouble( 0.248, -0.329, 0.148, -0.0222 ),
    xMin = cms.double( 1.5 ),
    yMax = cms.int32( 99999 ),
    xMax = cms.double( 3.0 ),
    funcType = cms.string( "TF1:=pol3" )
  )
) ),
    productsToWrite = cms.int32( 0 )
)
process.hltEgammaCkfTrackCandidatesForGSF = cms.EDProducer( "CkfTrackCandidateMaker",
    cleanTrajectoryAfterInOut = cms.bool( True ),
    doSeedingRegionRebuilding = cms.bool( True ),
    onlyPixelHitsForSeedCleaner = cms.bool( False ),
    reverseTrajectories = cms.bool( False ),
    useHitsSplitting = cms.bool( True ),
    MeasurementTrackerEvent = cms.InputTag( "hltSiStripClusters" ),
    src = cms.InputTag( "hltEgammaElectronPixelSeeds" ),
    clustersToSkip = cms.InputTag( "" ),
    phase2clustersToSkip = cms.InputTag( "" ),
    TrajectoryBuilderPSet = cms.PSet(  refToPSet_ = cms.string( "HLTPSetTrajectoryBuilderForGsfElectrons" ) ),
    TransientInitialStateEstimatorParameters = cms.PSet( 
      propagatorAlongTISE = cms.string( "PropagatorWithMaterial" ),
      numberMeasurementsForFit = cms.int32( 4 ),
      propagatorOppositeTISE = cms.string( "PropagatorWithMaterialOpposite" )
    ),
    numHitsForSeedCleaner = cms.int32( 4 ),
    NavigationSchool = cms.string( "SimpleNavigationSchool" ),
    RedundantSeedCleaner = cms.string( "CachingSeedCleanerBySharedInput" ),
    TrajectoryCleaner = cms.string( "hltESPTrajectoryCleanerBySharedHits" ),
    maxNSeeds = cms.uint32( 1000000 ),
    maxSeedsBeforeCleaning = cms.uint32( 1000 )
)
process.hltEgammaGsfTracks = cms.EDProducer( "GsfTrackProducer",
    src = cms.InputTag( "hltEgammaCkfTrackCandidatesForGSF" ),
    beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
    producer = cms.string( "" ),
    Fitter = cms.string( "hltESPGsfElectronFittingSmoother" ),
    useHitsSplitting = cms.bool( False ),
    TrajectoryInEvent = cms.bool( True ),
    TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
    Propagator = cms.string( "hltESPFwdElectronPropagator" ),
    NavigationSchool = cms.string( "SimpleNavigationSchool" ),
    MeasurementTracker = cms.string( "hltESPMeasurementTracker" ),
    MeasurementTrackerEvent = cms.InputTag( "hltSiStripClusters" ),
    GeometricInnerState = cms.bool( True ),
    AlgorithmName = cms.string( "gsf" )
)
process.hltEgammaGsfElectrons = cms.EDProducer( "EgammaHLTPixelMatchElectronProducers",
    TrackProducer = cms.InputTag( "" ),
    GsfTrackProducer = cms.InputTag( "hltEgammaGsfTracks" ),
    UseGsfTracks = cms.bool( True ),
    BSProducer = cms.InputTag( "hltOnlineBeamSpot" )
)
process.hltEgammaGsfTrackVars = cms.EDProducer( "EgammaHLTGsfTrackVarProducer",
    recoEcalCandidateProducer = cms.InputTag( "hltEgammaCandidates" ),
    inputCollection = cms.InputTag( "hltEgammaGsfTracks" ),
    beamSpotProducer = cms.InputTag( "hltOnlineBeamSpot" ),
    upperTrackNrToRemoveCut = cms.int32( 9999 ),
    lowerTrackNrToRemoveCut = cms.int32( -1 ),
    useDefaultValuesForBarrel = cms.bool( False ),
    useDefaultValuesForEndcap = cms.bool( False )
)
process.hltEgammaEleGsfTrackIsoPixelOnly = cms.EDProducer( "EgammaHLTElectronTrackIsolationProducers",
    electronProducer = cms.InputTag( "hltEgammaGsfElectrons" ),
    trackProducer = cms.InputTag( "hltPixelTracksZetaClean" ),
    recoEcalCandidateProducer = cms.InputTag( "hltEgammaCandidates" ),
    beamSpotProducer = cms.InputTag( "hltOnlineBeamSpot" ),
    egTrkIsoPtMin = cms.double( 1.0 ),
    egTrkIsoConeSize = cms.double( 0.2 ),
    egTrkIsoZSpan = cms.double( 0.15 ),
    egTrkIsoRSpan = cms.double( 999999.0 ),
    egTrkIsoVetoConeSizeBarrel = cms.double( 0.03 ),
    egTrkIsoVetoConeSizeEndcap = cms.double( 0.03 ),
    egTrkIsoStripBarrel = cms.double( 0.01 ),
    egTrkIsoStripEndcap = cms.double( 0.01 ),
    useGsfTrack = cms.bool( True ),
    useSCRefs = cms.bool( True )
)
process.hltIter0PFLowPixelSeedsFromPixelTracks = cms.EDProducer( "SeedGeneratorFromProtoTracksEDProducer",
    InputCollection = cms.InputTag( "hltPixelTracks" ),
    InputVertexCollection = cms.InputTag( "hltTrimmedPixelVertices" ),
    originHalfLength = cms.double( 0.3 ),
    originRadius = cms.double( 0.1 ),
    useProtoTrackKinematics = cms.bool( False ),
    useEventsWithNoVertex = cms.bool( True ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    usePV = cms.bool( False ),
    includeFourthHit = cms.bool( True ),
    SeedCreatorPSet = cms.PSet(  refToPSet_ = cms.string( "HLTSeedFromProtoTracks" ) )
)
process.hltIter0PFlowCkfTrackCandidates = cms.EDProducer( "CkfTrackCandidateMaker",
    cleanTrajectoryAfterInOut = cms.bool( False ),
    doSeedingRegionRebuilding = cms.bool( False ),
    onlyPixelHitsForSeedCleaner = cms.bool( False ),
    reverseTrajectories = cms.bool( False ),
    useHitsSplitting = cms.bool( False ),
    MeasurementTrackerEvent = cms.InputTag( "hltSiStripClusters" ),
    src = cms.InputTag( "hltIter0PFLowPixelSeedsFromPixelTracks" ),
    clustersToSkip = cms.InputTag( "" ),
    phase2clustersToSkip = cms.InputTag( "" ),
    TrajectoryBuilderPSet = cms.PSet(  refToPSet_ = cms.string( "HLTIter0GroupedCkfTrajectoryBuilderIT" ) ),
    TransientInitialStateEstimatorParameters = cms.PSet( 
      propagatorAlongTISE = cms.string( "PropagatorWithMaterialParabolicMf" ),
      numberMeasurementsForFit = cms.int32( 4 ),
      propagatorOppositeTISE = cms.string( "PropagatorWithMaterialParabolicMfOpposite" )
    ),
    numHitsForSeedCleaner = cms.int32( 4 ),
    NavigationSchool = cms.string( "SimpleNavigationSchool" ),
    RedundantSeedCleaner = cms.string( "CachingSeedCleanerBySharedInput" ),
    TrajectoryCleaner = cms.string( "hltESPTrajectoryCleanerBySharedHits" ),
    maxNSeeds = cms.uint32( 100000 ),
    maxSeedsBeforeCleaning = cms.uint32( 1000 )
)
process.hltIter0PFlowCtfWithMaterialTracks = cms.EDProducer( "TrackProducer",
    useSimpleMF = cms.bool( True ),
    SimpleMagneticField = cms.string( "ParabolicMf" ),
    src = cms.InputTag( "hltIter0PFlowCkfTrackCandidates" ),
    clusterRemovalInfo = cms.InputTag( "" ),
    beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
    Fitter = cms.string( "hltESPFittingSmootherIT" ),
    useHitsSplitting = cms.bool( False ),
    alias = cms.untracked.string( "ctfWithMaterialTracks" ),
    TrajectoryInEvent = cms.bool( False ),
    TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
    AlgorithmName = cms.string( "hltIter0" ),
    Propagator = cms.string( "hltESPRungeKuttaTrackerPropagator" ),
    GeometricInnerState = cms.bool( True ),
    NavigationSchool = cms.string( "" ),
    MeasurementTracker = cms.string( "" ),
    MeasurementTrackerEvent = cms.InputTag( "hltSiStripClusters" )
)
process.hltIter0PFlowTrackCutClassifier = cms.EDProducer( "TrackCutClassifier",
    src = cms.InputTag( "hltIter0PFlowCtfWithMaterialTracks" ),
    beamspot = cms.InputTag( "hltOnlineBeamSpot" ),
    vertices = cms.InputTag( "hltTrimmedPixelVertices" ),
    ignoreVertices = cms.bool( False ),
    qualityCuts = cms.vdouble( -0.7, 0.1, 0.7 ),
    mva = cms.PSet( 
      minPixelHits = cms.vint32( 0, 0, 0 ),
      maxDzWrtBS = cms.vdouble( 3.40282346639E38, 24.0, 15.0 ),
      dr_par = cms.PSet( 
        d0err = cms.vdouble( 0.003, 0.003, 0.003 ),
        dr_par2 = cms.vdouble( 3.40282346639E38, 0.6, 0.6 ),
        dr_par1 = cms.vdouble( 3.40282346639E38, 0.8, 0.8 ),
        dr_exp = cms.vint32( 4, 4, 4 ),
        d0err_par = cms.vdouble( 0.001, 0.001, 0.001 )
      ),
      maxLostLayers = cms.vint32( 1, 1, 1 ),
      min3DLayers = cms.vint32( 0, 0, 0 ),
      dz_par = cms.PSet( 
        dz_par1 = cms.vdouble( 3.40282346639E38, 0.75, 0.75 ),
        dz_par2 = cms.vdouble( 3.40282346639E38, 0.5, 0.5 ),
        dz_exp = cms.vint32( 4, 4, 4 )
      ),
      minNVtxTrk = cms.int32( 3 ),
      maxDz = cms.vdouble( 0.5, 0.2, 3.40282346639E38 ),
      minNdof = cms.vdouble( 1.0E-5, 1.0E-5, 1.0E-5 ),
      maxChi2 = cms.vdouble( 9999.0, 25.0, 16.0 ),
      maxChi2n = cms.vdouble( 1.2, 1.0, 0.7 ),
      maxDr = cms.vdouble( 0.5, 0.03, 3.40282346639E38 ),
      minLayers = cms.vint32( 3, 3, 3 )
    )
)
process.hltMergedTracks = cms.EDProducer( "TrackCollectionFilterCloner",
    originalSource = cms.InputTag( "hltIter0PFlowCtfWithMaterialTracks" ),
    originalMVAVals = cms.InputTag( 'hltIter0PFlowTrackCutClassifier','MVAValues' ),
    originalQualVals = cms.InputTag( 'hltIter0PFlowTrackCutClassifier','QualityMasks' ),
    minQuality = cms.string( "highPurity" ),
    copyExtras = cms.untracked.bool( True ),
    copyTrajectories = cms.untracked.bool( False )
)
process.hltPFMuonMerging = cms.EDProducer( "TrackListMerger",
    ShareFrac = cms.double( 0.19 ),
    FoundHitBonus = cms.double( 5.0 ),
    LostHitPenalty = cms.double( 20.0 ),
    MinPT = cms.double( 0.05 ),
    Epsilon = cms.double( -0.001 ),
    MaxNormalizedChisq = cms.double( 1000.0 ),
    MinFound = cms.int32( 3 ),
    TrackProducers = cms.VInputTag( 'hltIterL3MuonTracks','hltMergedTracks' ),
    hasSelector = cms.vint32( 0, 0 ),
    indivShareFrac = cms.vdouble( 1.0, 1.0 ),
    selectedTrackQuals = cms.VInputTag( 'hltIterL3MuonTracks','hltMergedTracks' ),
    setsToMerge = cms.VPSet( 
      cms.PSet(  pQual = cms.bool( False ),
        tLists = cms.vint32( 0, 1 )
      )
    ),
    trackAlgoPriorityOrder = cms.string( "hltESPTrackAlgoPriorityOrder" ),
    allowFirstHitShare = cms.bool( True ),
    newQuality = cms.string( "confirmed" ),
    copyExtras = cms.untracked.bool( True ),
    writeOnlyTrkQuals = cms.bool( False ),
    copyMVA = cms.bool( False )
)
process.hltMuonLinks = cms.EDProducer( "MuonLinksProducerForHLT",
    InclusiveTrackerTrackCollection = cms.InputTag( "hltPFMuonMerging" ),
    LinkCollection = cms.InputTag( "hltL3MuonsIterL3Links" ),
    ptMin = cms.double( 2.5 ),
    pMin = cms.double( 2.5 ),
    shareHitFraction = cms.double( 0.8 )
)
process.hltMuons = cms.EDProducer( "MuonIdProducer",
    MuonCaloCompatibility = cms.PSet( 
      delta_eta = cms.double( 0.02 ),
      delta_phi = cms.double( 0.02 ),
      allSiPMHO = cms.bool( False ),
      MuonTemplateFileName = cms.FileInPath( "RecoMuon/MuonIdentification/data/MuID_templates_muons_lowPt_3_1_norm.root" ),
      PionTemplateFileName = cms.FileInPath( "RecoMuon/MuonIdentification/data/MuID_templates_pions_lowPt_3_1_norm.root" )
    ),
    TrackAssociatorParameters = cms.PSet( 
      useMuon = cms.bool( True ),
      truthMatch = cms.bool( False ),
      usePreshower = cms.bool( False ),
      dRPreshowerPreselection = cms.double( 0.2 ),
      muonMaxDistanceSigmaY = cms.double( 0.0 ),
      useEcal = cms.bool( True ),
      muonMaxDistanceSigmaX = cms.double( 0.0 ),
      dRMuon = cms.double( 9999.0 ),
      dREcal = cms.double( 9999.0 ),
      CSCSegmentCollectionLabel = cms.InputTag( "hltCscSegments" ),
      DTRecSegment4DCollectionLabel = cms.InputTag( "hltDt4DSegments" ),
      EBRecHitCollectionLabel = cms.InputTag( 'hltEcalRecHit','EcalRecHitsEB' ),
      CaloTowerCollectionLabel = cms.InputTag( "hltTowerMakerForAll" ),
      propagateAllDirections = cms.bool( True ),
      muonMaxDistanceY = cms.double( 5.0 ),
      useHO = cms.bool( True ),
      muonMaxDistanceX = cms.double( 5.0 ),
      trajectoryUncertaintyTolerance = cms.double( -1.0 ),
      useHcal = cms.bool( True ),
      HBHERecHitCollectionLabel = cms.InputTag( "hltHbhereco" ),
      accountForTrajectoryChangeCalo = cms.bool( False ),
      dREcalPreselection = cms.double( 0.05 ),
      useCalo = cms.bool( False ),
      dRMuonPreselection = cms.double( 0.2 ),
      EERecHitCollectionLabel = cms.InputTag( 'hltEcalRecHit','EcalRecHitsEE' ),
      dRHcal = cms.double( 9999.0 ),
      dRHcalPreselection = cms.double( 0.2 ),
      HORecHitCollectionLabel = cms.InputTag( "hltHoreco" )
    ),
    CaloExtractorPSet = cms.PSet( 
      DR_Veto_H = cms.double( 0.1 ),
      CenterConeOnCalIntersection = cms.bool( False ),
      NoiseTow_EE = cms.double( 0.15 ),
      Noise_EB = cms.double( 0.025 ),
      Noise_HE = cms.double( 0.2 ),
      DR_Veto_E = cms.double( 0.07 ),
      NoiseTow_EB = cms.double( 0.04 ),
      Noise_EE = cms.double( 0.1 ),
      UseRecHitsFlag = cms.bool( False ),
      DR_Max = cms.double( 1.0 ),
      DepositLabel = cms.untracked.string( "Cal" ),
      Noise_HO = cms.double( 0.2 ),
      DR_Veto_HO = cms.double( 0.1 ),
      Threshold_H = cms.double( 0.5 ),
      PrintTimeReport = cms.untracked.bool( False ),
      Threshold_E = cms.double( 0.2 ),
      PropagatorName = cms.string( "hltESPFastSteppingHelixPropagatorAny" ),
      ComponentName = cms.string( "CaloExtractorByAssociator" ),
      Threshold_HO = cms.double( 0.5 ),
      DepositInstanceLabels = cms.vstring( 'ecal',
        'hcal',
        'ho' ),
      ServiceParameters = cms.PSet( 
        RPCLayers = cms.bool( False ),
        UseMuonNavigation = cms.untracked.bool( False ),
        Propagators = cms.untracked.vstring( 'hltESPFastSteppingHelixPropagatorAny' )
      ),
      TrackAssociatorParameters = cms.PSet( 
        useMuon = cms.bool( False ),
        truthMatch = cms.bool( False ),
        usePreshower = cms.bool( False ),
        dRPreshowerPreselection = cms.double( 0.2 ),
        muonMaxDistanceSigmaY = cms.double( 0.0 ),
        useEcal = cms.bool( False ),
        muonMaxDistanceSigmaX = cms.double( 0.0 ),
        dRMuon = cms.double( 9999.0 ),
        dREcal = cms.double( 1.0 ),
        CSCSegmentCollectionLabel = cms.InputTag( "hltCscSegments" ),
        DTRecSegment4DCollectionLabel = cms.InputTag( "hltDt4DSegments" ),
        EBRecHitCollectionLabel = cms.InputTag( 'hltEcalRecHit','EcalRecHitsEB' ),
        CaloTowerCollectionLabel = cms.InputTag( "hltTowerMakerForAll" ),
        propagateAllDirections = cms.bool( True ),
        muonMaxDistanceY = cms.double( 5.0 ),
        useHO = cms.bool( False ),
        muonMaxDistanceX = cms.double( 5.0 ),
        trajectoryUncertaintyTolerance = cms.double( -1.0 ),
        useHcal = cms.bool( False ),
        HBHERecHitCollectionLabel = cms.InputTag( "hltHbhereco" ),
        accountForTrajectoryChangeCalo = cms.bool( False ),
        dREcalPreselection = cms.double( 1.0 ),
        useCalo = cms.bool( True ),
        dRMuonPreselection = cms.double( 0.2 ),
        EERecHitCollectionLabel = cms.InputTag( 'hltEcalRecHit','EcalRecHitsEE' ),
        dRHcal = cms.double( 1.0 ),
        dRHcalPreselection = cms.double( 1.0 ),
        HORecHitCollectionLabel = cms.InputTag( "hltHoreco" )
      ),
      Noise_HB = cms.double( 0.2 )
    ),
    TrackExtractorPSet = cms.PSet( 
      Diff_z = cms.double( 0.2 ),
      inputTrackCollection = cms.InputTag( "hltPFMuonMerging" ),
      Chi2Ndof_Max = cms.double( 1.0E64 ),
      BeamSpotLabel = cms.InputTag( "hltOnlineBeamSpot" ),
      DR_Veto = cms.double( 0.01 ),
      Pt_Min = cms.double( -1.0 ),
      DR_Max = cms.double( 1.0 ),
      DepositLabel = cms.untracked.string( "" ),
      NHits_Min = cms.uint32( 0 ),
      Chi2Prob_Min = cms.double( -1.0 ),
      Diff_r = cms.double( 0.1 ),
      BeamlineOption = cms.string( "BeamSpotFromEvent" ),
      ComponentName = cms.string( "TrackExtractor" )
    ),
    JetExtractorPSet = cms.PSet( 
      JetCollectionLabel = cms.InputTag( "hltAK4CaloJetsPFEt5" ),
      DR_Veto = cms.double( 0.1 ),
      DR_Max = cms.double( 1.0 ),
      ExcludeMuonVeto = cms.bool( True ),
      PrintTimeReport = cms.untracked.bool( False ),
      PropagatorName = cms.string( "hltESPFastSteppingHelixPropagatorAny" ),
      ComponentName = cms.string( "JetExtractor" ),
      ServiceParameters = cms.PSet( 
        RPCLayers = cms.bool( False ),
        UseMuonNavigation = cms.untracked.bool( False ),
        Propagators = cms.untracked.vstring( 'hltESPFastSteppingHelixPropagatorAny' )
      ),
      TrackAssociatorParameters = cms.PSet( 
        useMuon = cms.bool( False ),
        truthMatch = cms.bool( False ),
        usePreshower = cms.bool( False ),
        dRPreshowerPreselection = cms.double( 0.2 ),
        muonMaxDistanceSigmaY = cms.double( 0.0 ),
        useEcal = cms.bool( False ),
        muonMaxDistanceSigmaX = cms.double( 0.0 ),
        dRMuon = cms.double( 9999.0 ),
        dREcal = cms.double( 0.5 ),
        CSCSegmentCollectionLabel = cms.InputTag( "hltCscSegments" ),
        DTRecSegment4DCollectionLabel = cms.InputTag( "hltDt4DSegments" ),
        EBRecHitCollectionLabel = cms.InputTag( 'hltEcalRecHit','EcalRecHitsEB' ),
        CaloTowerCollectionLabel = cms.InputTag( "hltTowerMakerForAll" ),
        propagateAllDirections = cms.bool( True ),
        muonMaxDistanceY = cms.double( 5.0 ),
        useHO = cms.bool( False ),
        muonMaxDistanceX = cms.double( 5.0 ),
        trajectoryUncertaintyTolerance = cms.double( -1.0 ),
        useHcal = cms.bool( False ),
        HBHERecHitCollectionLabel = cms.InputTag( "hltHbhereco" ),
        accountForTrajectoryChangeCalo = cms.bool( False ),
        dREcalPreselection = cms.double( 0.5 ),
        useCalo = cms.bool( True ),
        dRMuonPreselection = cms.double( 0.2 ),
        EERecHitCollectionLabel = cms.InputTag( 'hltEcalRecHit','EcalRecHitsEE' ),
        dRHcal = cms.double( 0.5 ),
        dRHcalPreselection = cms.double( 0.5 ),
        HORecHitCollectionLabel = cms.InputTag( "hltHoreco" )
      ),
      Threshold = cms.double( 5.0 )
    ),
    trackDepositName = cms.string( "tracker" ),
    ecalDepositName = cms.string( "ecal" ),
    hcalDepositName = cms.string( "hcal" ),
    hoDepositName = cms.string( "ho" ),
    jetDepositName = cms.string( "jets" ),
    TimingFillerParameters = cms.PSet( 
      DTTimingParameters = cms.PSet( 
        HitError = cms.double( 6.0 ),
        MatchParameters = cms.PSet( 
          TightMatchDT = cms.bool( False ),
          DTradius = cms.double( 0.01 ),
          TightMatchCSC = cms.bool( True ),
          CSCsegments = cms.InputTag( "hltCscSegments" ),
          DTsegments = cms.InputTag( "hltDt4DSegments" )
        ),
        debug = cms.bool( False ),
        DoWireCorr = cms.bool( False ),
        RequireBothProjections = cms.bool( False ),
        DTTimeOffset = cms.double( 2.7 ),
        PruneCut = cms.double( 10000.0 ),
        DTsegments = cms.InputTag( "hltDt4DSegments" ),
        UseSegmentT0 = cms.bool( False ),
        HitsMin = cms.int32( 5 ),
        DropTheta = cms.bool( True ),
        ServiceParameters = cms.PSet( 
          RPCLayers = cms.bool( True ),
          Propagators = cms.untracked.vstring( 'hltESPFastSteppingHelixPropagatorAny' )
        )
      ),
      UseCSC = cms.bool( True ),
      CSCTimingParameters = cms.PSet( 
        MatchParameters = cms.PSet( 
          TightMatchDT = cms.bool( False ),
          DTradius = cms.double( 0.01 ),
          TightMatchCSC = cms.bool( True ),
          CSCsegments = cms.InputTag( "hltCscSegments" ),
          DTsegments = cms.InputTag( "hltDt4DSegments" )
        ),
        debug = cms.bool( False ),
        CSCWireTimeOffset = cms.double( 0.0 ),
        CSCStripError = cms.double( 7.0 ),
        CSCTimeOffset = cms.double( 0.0 ),
        CSCWireError = cms.double( 8.6 ),
        PruneCut = cms.double( 100.0 ),
        CSCsegments = cms.InputTag( "hltCscSegments" ),
        UseStripTime = cms.bool( True ),
        CSCStripTimeOffset = cms.double( 0.0 ),
        UseWireTime = cms.bool( True ),
        ServiceParameters = cms.PSet( 
          RPCLayers = cms.bool( True ),
          Propagators = cms.untracked.vstring( 'hltESPFastSteppingHelixPropagatorAny' )
        )
      ),
      ErrorDT = cms.double( 6.0 ),
      EcalEnergyCut = cms.double( 0.4 ),
      UseECAL = cms.bool( True ),
      ErrorEB = cms.double( 2.085 ),
      UseDT = cms.bool( True ),
      ErrorEE = cms.double( 6.95 ),
      ErrorCSC = cms.double( 7.4 )
    ),
    ShowerDigiFillerParameters = cms.PSet( 
      cscDigiCollectionLabel = cms.InputTag( 'muonCSCDigis','MuonCSCStripDigi' ),
      dtDigiCollectionLabel = cms.InputTag( "muonDTDigis" ),
      digiMaxDistanceX = cms.double( 25.0 )
    ),
    TrackerKinkFinderParameters = cms.PSet( 
      usePosition = cms.bool( False ),
      diagonalOnly = cms.bool( False )
    ),
    fillEnergy = cms.bool( True ),
    storeCrossedHcalRecHits = cms.bool( False ),
    maxAbsPullX = cms.double( 4.0 ),
    maxAbsEta = cms.double( 3.0 ),
    minPt = cms.double( 10.0 ),
    inputCollectionTypes = cms.vstring( 'inner tracks',
      'links',
      'outer tracks' ),
    addExtraSoftMuons = cms.bool( False ),
    fillGlobalTrackRefits = cms.bool( False ),
    debugWithTruthMatching = cms.bool( False ),
    inputCollectionLabels = cms.VInputTag( 'hltPFMuonMerging','hltMuonLinks','hltL2Muons' ),
    fillCaloCompatibility = cms.bool( True ),
    maxAbsPullY = cms.double( 9999.0 ),
    maxAbsDy = cms.double( 9999.0 ),
    minP = cms.double( 10.0 ),
    minPCaloMuon = cms.double( 1.0E9 ),
    maxAbsDx = cms.double( 3.0 ),
    fillIsolation = cms.bool( True ),
    writeIsoDeposits = cms.bool( False ),
    minNumberOfMatches = cms.int32( 1 ),
    fillMatching = cms.bool( True ),
    fillShowerDigis = cms.bool( False ),
    ptThresholdToFillCandidateP4WithGlobalFit = cms.double( 200.0 ),
    sigmaThresholdToFillCandidateP4WithGlobalFit = cms.double( 2.0 ),
    fillGlobalTrackQuality = cms.bool( False ),
    globalTrackQualityInputTag = cms.InputTag( "glbTrackQual" ),
    selectHighPurity = cms.bool( False ),
    pvInputTag = cms.InputTag( "offlinePrimaryVertices" ),
    fillTrackerKink = cms.bool( False ),
    minCaloCompatibility = cms.double( 0.6 ),
    runArbitrationCleaner = cms.bool( False ),
    arbitrationCleanerOptions = cms.PSet( 
      OverlapDTheta = cms.double( 0.02 ),
      Overlap = cms.bool( True ),
      Clustering = cms.bool( True ),
      ME1a = cms.bool( True ),
      ClusterDTheta = cms.double( 0.02 ),
      ClusterDPhi = cms.double( 0.6 ),
      OverlapDPhi = cms.double( 0.0786 )
    ),
    arbitrateTrackerMuons = cms.bool( False )
)
process.hltLightPFTracks = cms.EDProducer( "LightPFTrackProducer",
    TrackQuality = cms.string( "none" ),
    UseQuality = cms.bool( False ),
    TkColList = cms.VInputTag( 'hltPFMuonMerging' )
)
process.hltParticleFlowBlock = cms.EDProducer( "PFBlockProducer",
    verbose = cms.untracked.bool( False ),
    debug = cms.untracked.bool( False ),
    elementImporters = cms.VPSet( 
      cms.PSet(  muonSrc = cms.InputTag( "hltMuons" ),
        source = cms.InputTag( "hltLightPFTracks" ),
        NHitCuts_byTrackAlgo = cms.vuint32( 3, 3, 3, 3, 3, 3 ),
        useIterativeTracking = cms.bool( False ),
        importerName = cms.string( "GeneralTracksImporter" ),
        DPtOverPtCuts_byTrackAlgo = cms.vdouble( 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 ),
        muonMaxDPtOPt = cms.double( 1.0 ),
        trackQuality = cms.string( "highPurity" ),
        cleanBadConvertedBrems = cms.bool( False ),
        vetoEndcap = cms.bool( False )
      ),
      cms.PSet(  source = cms.InputTag( "hltParticleFlowClusterECALUnseeded" ),
        importerName = cms.string( "ECALClusterImporter" ),
        BCtoPFCMap = cms.InputTag( "" )
      ),
      cms.PSet(  source = cms.InputTag( "hltParticleFlowClusterHCAL" ),
        importerName = cms.string( "GenericClusterImporter" )
      ),
      cms.PSet(  source = cms.InputTag( "hltParticleFlowClusterHF" ),
        importerName = cms.string( "GenericClusterImporter" )
      )
    ),
    linkDefinitions = cms.VPSet( 
      cms.PSet(  linkType = cms.string( "TRACK:ECAL" ),
        useKDTree = cms.bool( True ),
        linkerName = cms.string( "TrackAndECALLinker" )
      ),
      cms.PSet(  linkType = cms.string( "TRACK:HCAL" ),
        useKDTree = cms.bool( True ),
        linkerName = cms.string( "TrackAndHCALLinker" ),
        trajectoryLayerEntrance = cms.string( "HCALEntrance" ),
        trajectoryLayerExit = cms.string( "HCALExit" ),
        nMaxHcalLinksPerTrack = cms.int32( 1 )
      ),
      cms.PSet(  linkType = cms.string( "ECAL:HCAL" ),
        useKDTree = cms.bool( False ),
        linkerName = cms.string( "ECALAndHCALLinker" ),
        minAbsEtaEcal = cms.double( 2.5 )
      ),
      cms.PSet(  linkType = cms.string( "HFEM:HFHAD" ),
        useKDTree = cms.bool( False ),
        linkerName = cms.string( "HFEMAndHFHADLinker" )
      )
    )
)
process.hltParticleFlow = cms.EDProducer( "PFProducer",
    verbose = cms.untracked.bool( False ),
    debug = cms.untracked.bool( False ),
    blocks = cms.InputTag( "hltParticleFlowBlock" ),
    muons = cms.InputTag( "hltMuons" ),
    postMuonCleaning = cms.bool( True ),
    vetoEndcap = cms.bool( False ),
    vertexCollection = cms.InputTag( "hltPixelVertices" ),
    useVerticesForNeutral = cms.bool( True ),
    useHO = cms.bool( False ),
    PFEGammaCandidates = cms.InputTag( "particleFlowEGamma" ),
    GedElectronValueMap = cms.InputTag( "gedGsfElectronsTmp" ),
    GedPhotonValueMap = cms.InputTag( 'tmpGedPhotons','valMapPFEgammaCandToPhoton' ),
    useEGammaElectrons = cms.bool( False ),
    egammaElectrons = cms.InputTag( "" ),
    useEGammaFilters = cms.bool( False ),
    useProtectionsForJetMET = cms.bool( True ),
    PFEGammaFiltersParameters = cms.PSet( 
      electron_protectionsForJetMET = cms.PSet( 
        maxE = cms.double( 50.0 ),
        maxTrackPOverEele = cms.double( 1.0 ),
        maxEcalEOverP_2 = cms.double( 0.2 ),
        maxHcalEOverEcalE = cms.double( 0.1 ),
        maxEcalEOverP_1 = cms.double( 0.5 ),
        maxHcalEOverP = cms.double( 1.0 ),
        maxEcalEOverPRes = cms.double( 0.2 ),
        maxHcalE = cms.double( 10.0 ),
        maxEeleOverPout = cms.double( 0.2 ),
        maxNtracks = cms.double( 3.0 ),
        maxEleHcalEOverEcalE = cms.double( 0.1 ),
        maxDPhiIN = cms.double( 0.1 ),
        maxEeleOverPoutRes = cms.double( 0.5 )
      ),
      electron_maxElePtForOnlyMVAPresel = cms.double( 50.0 ),
      photon_SigmaiEtaiEta_endcap = cms.double( 0.034 ),
      electron_iso_combIso_endcap = cms.double( 10.0 ),
      photon_protectionsForBadHcal = cms.PSet( 
        solidConeTrkIsoSlope = cms.double( 0.3 ),
        enableProtections = cms.bool( False ),
        solidConeTrkIsoOffset = cms.double( 10.0 )
      ),
      electron_missinghits = cms.uint32( 1 ),
      photon_MinEt = cms.double( 10.0 ),
      electron_iso_pt = cms.double( 10.0 ),
      electron_ecalDrivenHademPreselCut = cms.double( 0.15 ),
      electron_iso_mva_endcap = cms.double( -0.1075 ),
      electron_iso_combIso_barrel = cms.double( 10.0 ),
      photon_protectionsForJetMET = cms.PSet( 
        sumPtTrackIsoSlope = cms.double( 0.001 ),
        sumPtTrackIso = cms.double( 4.0 )
      ),
      electron_protectionsForBadHcal = cms.PSet( 
        dEta = cms.vdouble( 0.0064, 0.01264 ),
        dPhi = cms.vdouble( 0.0547, 0.0394 ),
        enableProtections = cms.bool( False ),
        eInvPInv = cms.vdouble( 0.184, 0.0721 ),
        full5x5_sigmaIetaIeta = cms.vdouble( 0.0106, 0.0387 )
      ),
      electron_noniso_mvaCut = cms.double( -0.1 ),
      electron_iso_mva_barrel = cms.double( -0.1875 ),
      photon_SigmaiEtaiEta_barrel = cms.double( 0.0125 ),
      photon_combIso = cms.double( 10.0 ),
      photon_HoE = cms.double( 0.05 )
    ),
    muon_HCAL = cms.vdouble( 3.0, 3.0 ),
    muon_ECAL = cms.vdouble( 0.5, 0.5 ),
    muon_HO = cms.vdouble( 0.9, 0.9 ),
    PFMuonAlgoParameters = cms.PSet(  ),
    rejectTracks_Bad = cms.bool( False ),
    rejectTracks_Step45 = cms.bool( False ),
    usePFNuclearInteractions = cms.bool( False ),
    usePFConversions = cms.bool( False ),
    usePFDecays = cms.bool( False ),
    dptRel_DispVtx = cms.double( 10.0 ),
    iCfgCandConnector = cms.PSet( 
      nuclCalibFactors = cms.vdouble( 0.8, 0.15, 0.5, 0.5, 0.05 ),
      bCorrect = cms.bool( False ),
      bCalibPrimary = cms.bool( False )
    ),
    nsigma_TRACK = cms.double( 1.0 ),
    pt_Error = cms.double( 1.0 ),
    factors_45 = cms.vdouble( 10.0, 100.0 ),
    goodTrackDeadHcal_ptErrRel = cms.double( 0.2 ),
    goodTrackDeadHcal_chi2n = cms.double( 5.0 ),
    goodTrackDeadHcal_layers = cms.uint32( 4 ),
    goodTrackDeadHcal_validFr = cms.double( 0.5 ),
    goodTrackDeadHcal_dxy = cms.double( 0.5 ),
    goodPixelTrackDeadHcal_minEta = cms.double( 2.3 ),
    goodPixelTrackDeadHcal_maxPt = cms.double( 50.0 ),
    goodPixelTrackDeadHcal_ptErrRel = cms.double( 1.0 ),
    goodPixelTrackDeadHcal_chi2n = cms.double( 2.0 ),
    goodPixelTrackDeadHcal_maxLost3Hit = cms.int32( 0 ),
    goodPixelTrackDeadHcal_maxLost4Hit = cms.int32( 1 ),
    goodPixelTrackDeadHcal_dxy = cms.double( 0.02 ),
    goodPixelTrackDeadHcal_dz = cms.double( 0.05 ),
    pf_nsigma_ECAL = cms.double( 0.0 ),
    pf_nsigma_HCAL = cms.double( 1.0 ),
    pf_nsigma_HFEM = cms.double( 1.0 ),
    pf_nsigma_HFHAD = cms.double( 1.0 ),
    useCalibrationsFromDB = cms.bool( True ),
    calibrationsLabel = cms.string( "HLT" ),
    postHFCleaning = cms.bool( False ),
    PFHFCleaningParameters = cms.PSet( 
      minSignificance = cms.double( 2.5 ),
      maxSignificance = cms.double( 2.5 ),
      minDeltaMet = cms.double( 0.4 ),
      maxDeltaPhiPt = cms.double( 7.0 ),
      minHFCleaningPt = cms.double( 5.0 ),
      minSignificanceReduction = cms.double( 1.4 )
    ),
    cleanedHF = cms.VInputTag( 'hltParticleFlowRecHitHF:Cleaned','hltParticleFlowClusterHF:Cleaned' ),
    calibHF_use = cms.bool( False ),
    calibHF_eta_step = cms.vdouble( 0.0, 2.9, 3.0, 3.2, 4.2, 4.4, 4.6, 4.8, 5.2, 5.4 ),
    calibHF_a_EMonly = cms.vdouble( 0.96945, 0.96701, 0.76309, 0.82268, 0.87583, 0.89718, 0.98674, 1.4681, 1.458, 1.458 ),
    calibHF_a_EMHAD = cms.vdouble( 1.42215, 1.00496, 0.68961, 0.81656, 0.98504, 0.98504, 1.00802, 1.0593, 1.4576, 1.4576 ),
    calibHF_b_HADonly = cms.vdouble( 1.27541, 0.85361, 0.86333, 0.89091, 0.94348, 0.94348, 0.9437, 1.0034, 1.0444, 1.0444 ),
    calibHF_b_EMHAD = cms.vdouble( 1.27541, 0.85361, 0.86333, 0.89091, 0.94348, 0.94348, 0.9437, 1.0034, 1.0444, 1.0444 ),
    resolHF_square = cms.vdouble( 7.834401, 0.012996, 0.0 )
)
process.hltAK4PFJets = cms.EDProducer( "FastjetJetProducer",
    useMassDropTagger = cms.bool( False ),
    useFiltering = cms.bool( False ),
    useDynamicFiltering = cms.bool( False ),
    useTrimming = cms.bool( False ),
    usePruning = cms.bool( False ),
    useCMSBoostedTauSeedingAlgorithm = cms.bool( False ),
    useKtPruning = cms.bool( False ),
    useConstituentSubtraction = cms.bool( False ),
    useSoftDrop = cms.bool( False ),
    correctShape = cms.bool( False ),
    UseOnlyVertexTracks = cms.bool( False ),
    UseOnlyOnePV = cms.bool( False ),
    muCut = cms.double( -1.0 ),
    yCut = cms.double( -1.0 ),
    rFilt = cms.double( -1.0 ),
    rFiltFactor = cms.double( -1.0 ),
    trimPtFracMin = cms.double( -1.0 ),
    zcut = cms.double( -1.0 ),
    rcut_factor = cms.double( -1.0 ),
    csRho_EtaMax = cms.double( -1.0 ),
    csRParam = cms.double( -1.0 ),
    beta = cms.double( -1.0 ),
    R0 = cms.double( -1.0 ),
    gridMaxRapidity = cms.double( -1.0 ),
    gridSpacing = cms.double( -1.0 ),
    DzTrVtxMax = cms.double( 0.0 ),
    DxyTrVtxMax = cms.double( 0.0 ),
    MaxVtxZ = cms.double( 15.0 ),
    subjetPtMin = cms.double( -1.0 ),
    muMin = cms.double( -1.0 ),
    muMax = cms.double( -1.0 ),
    yMin = cms.double( -1.0 ),
    yMax = cms.double( -1.0 ),
    dRMin = cms.double( -1.0 ),
    dRMax = cms.double( -1.0 ),
    maxDepth = cms.int32( -1 ),
    nFilt = cms.int32( -1 ),
    MinVtxNdof = cms.int32( 0 ),
    src = cms.InputTag( "hltParticleFlow" ),
    srcPVs = cms.InputTag( "hltPixelVertices" ),
    jetType = cms.string( "PFJet" ),
    jetAlgorithm = cms.string( "AntiKt" ),
    rParam = cms.double( 0.4 ),
    inputEtMin = cms.double( 0.0 ),
    inputEMin = cms.double( 0.0 ),
    jetPtMin = cms.double( 0.0 ),
    doPVCorrection = cms.bool( False ),
    doAreaFastjet = cms.bool( False ),
    doRhoFastjet = cms.bool( False ),
    doPUOffsetCorr = cms.bool( False ),
    puPtMin = cms.double( 10.0 ),
    nSigmaPU = cms.double( 1.0 ),
    radiusPU = cms.double( 0.4 ),
    subtractorName = cms.string( "" ),
    useExplicitGhosts = cms.bool( False ),
    doAreaDiskApprox = cms.bool( True ),
    voronoiRfact = cms.double( -9.0 ),
    Rho_EtaMax = cms.double( 4.4 ),
    Ghost_EtaMax = cms.double( 6.0 ),
    Active_Area_Repeats = cms.int32( 5 ),
    GhostArea = cms.double( 0.01 ),
    restrictInputs = cms.bool( False ),
    maxInputs = cms.uint32( 1 ),
    writeCompound = cms.bool( False ),
    writeJetsWithConst = cms.bool( False ),
    doFastJetNonUniform = cms.bool( False ),
    useDeterministicSeed = cms.bool( True ),
    minSeed = cms.uint32( 0 ),
    verbosity = cms.int32( 0 ),
    puWidth = cms.double( 0.0 ),
    nExclude = cms.uint32( 0 ),
    maxBadEcalCells = cms.uint32( 9999999 ),
    maxBadHcalCells = cms.uint32( 9999999 ),
    maxProblematicEcalCells = cms.uint32( 9999999 ),
    maxProblematicHcalCells = cms.uint32( 9999999 ),
    maxRecoveredEcalCells = cms.uint32( 9999999 ),
    maxRecoveredHcalCells = cms.uint32( 9999999 ),
    puCenters = cms.vdouble(  ),
    applyWeight = cms.bool( False ),
    srcWeights = cms.InputTag( "" ),
    minimumTowersFraction = cms.double( 0.0 ),
    jetCollInstanceName = cms.string( "" ),
    sumRecHits = cms.bool( False )
)
process.hltAK4PFJetsLooseID = cms.EDProducer( "HLTPFJetIDProducer",
    minPt = cms.double( 20.0 ),
    maxEta = cms.double( 1.0E99 ),
    CHF = cms.double( 0.0 ),
    NHF = cms.double( 0.99 ),
    CEF = cms.double( 0.99 ),
    NEF = cms.double( 0.99 ),
    maxCF = cms.double( 99.0 ),
    NCH = cms.int32( 0 ),
    NTOT = cms.int32( 1 ),
    jetsInput = cms.InputTag( "hltAK4PFJets" )
)
process.hltAK4PFJetsTightID = cms.EDProducer( "HLTPFJetIDProducer",
    minPt = cms.double( 20.0 ),
    maxEta = cms.double( 1.0E99 ),
    CHF = cms.double( 0.0 ),
    NHF = cms.double( 0.9 ),
    CEF = cms.double( 0.99 ),
    NEF = cms.double( 0.99 ),
    maxCF = cms.double( 99.0 ),
    NCH = cms.int32( 0 ),
    NTOT = cms.int32( 1 ),
    jetsInput = cms.InputTag( "hltAK4PFJets" )
)
process.hltFixedGridRhoFastjetAll = cms.EDProducer( "FixedGridRhoProducerFastjet",
    pfCandidatesTag = cms.InputTag( "hltParticleFlow" ),
    maxRapidity = cms.double( 5.0 ),
    gridSpacing = cms.double( 0.55 )
)
process.hltAK4PFFastJetCorrector = cms.EDProducer( "L1FastjetCorrectorProducer",
    level = cms.string( "L1FastJet" ),
    algorithm = cms.string( "AK4PFHLT" ),
    srcRho = cms.InputTag( "hltFixedGridRhoFastjetAll" )
)
process.hltAK4PFRelativeCorrector = cms.EDProducer( "LXXXCorrectorProducer",
    level = cms.string( "L2Relative" ),
    algorithm = cms.string( "AK4PFHLT" )
)
process.hltAK4PFAbsoluteCorrector = cms.EDProducer( "LXXXCorrectorProducer",
    level = cms.string( "L3Absolute" ),
    algorithm = cms.string( "AK4PFHLT" )
)
process.hltAK4PFResidualCorrector = cms.EDProducer( "LXXXCorrectorProducer",
    level = cms.string( "L2L3Residual" ),
    algorithm = cms.string( "AK4PFHLT" )
)
process.hltAK4PFCorrector = cms.EDProducer( "ChainedJetCorrectorProducer",
    correctors = cms.VInputTag( 'hltAK4PFFastJetCorrector','hltAK4PFRelativeCorrector','hltAK4PFAbsoluteCorrector','hltAK4PFResidualCorrector' )
)
process.hltAK4PFJetsCorrected = cms.EDProducer( "CorrectedPFJetProducer",
    src = cms.InputTag( "hltAK4PFJets" ),
    correctors = cms.VInputTag( 'hltAK4PFCorrector' )
)
process.hltAK4PFJetsLooseIDCorrected = cms.EDProducer( "CorrectedPFJetProducer",
    src = cms.InputTag( "hltAK4PFJetsLooseID" ),
    correctors = cms.VInputTag( 'hltAK4PFCorrector' )
)
process.hltAK4PFJetsTightIDCorrected = cms.EDProducer( "CorrectedPFJetProducer",
    src = cms.InputTag( "hltAK4PFJetsTightID" ),
    correctors = cms.VInputTag( 'hltAK4PFCorrector' )
)
process.hltEgammaR9ID = cms.EDProducer( "EgammaHLTR9IDProducer",
    recoEcalCandidateProducer = cms.InputTag( "hltEgammaCandidates" ),
    ecalRechitEB = cms.InputTag( 'hltEcalRecHit','EcalRecHitsEB' ),
    ecalRechitEE = cms.InputTag( 'hltEcalRecHit','EcalRecHitsEE' )
)
process.hltL2MuonCandidatesNoVtx = cms.EDProducer( "L2MuonCandidateProducer",
    InputObjects = cms.InputTag( "hltL2Muons" )
)
process.hltIterL3OISeedsFromL2MuonsNoVtx = cms.EDProducer( "TSGForOIFromL2",
    src = cms.InputTag( "hltL2Muons" ),
    layersToTry = cms.int32( 2 ),
    fixedErrorRescaleFactorForHitless = cms.double( 2.0 ),
    hitsToTry = cms.int32( 1 ),
    adjustErrorsDynamicallyForHits = cms.bool( False ),
    adjustErrorsDynamicallyForHitless = cms.bool( True ),
    MeasurementTrackerEvent = cms.InputTag( "hltSiStripClusters" ),
    UseHitLessSeeds = cms.bool( True ),
    estimator = cms.string( "hltESPChi2MeasurementEstimator100" ),
    maxEtaForTOB = cms.double( 1.8 ),
    minEtaForTEC = cms.double( 0.7 ),
    debug = cms.untracked.bool( False ),
    fixedErrorRescaleFactorForHits = cms.double( 1.0 ),
    maxSeeds = cms.uint32( 20 ),
    maxHitlessSeeds = cms.uint32( 5 ),
    maxHitSeeds = cms.uint32( 1 ),
    numL2ValidHitsCutAllEta = cms.uint32( 20 ),
    numL2ValidHitsCutAllEndcap = cms.uint32( 30 ),
    pT1 = cms.double( 13.0 ),
    pT2 = cms.double( 30.0 ),
    pT3 = cms.double( 70.0 ),
    eta1 = cms.double( 0.2 ),
    eta2 = cms.double( 0.3 ),
    eta3 = cms.double( 1.0 ),
    eta4 = cms.double( 1.2 ),
    eta5 = cms.double( 1.6 ),
    eta6 = cms.double( 1.4 ),
    eta7 = cms.double( 2.1 ),
    SF1 = cms.double( 3.0 ),
    SF2 = cms.double( 4.0 ),
    SF3 = cms.double( 5.0 ),
    SF4 = cms.double( 7.0 ),
    SF5 = cms.double( 10.0 ),
    SF6 = cms.double( 2.0 ),
    tsosDiff1 = cms.double( 0.2 ),
    tsosDiff2 = cms.double( 0.02 ),
    propagatorName = cms.string( "PropagatorWithMaterialParabolicMf" )
)
process.hltIterL3OITrackCandidatesNoVtx = cms.EDProducer( "CkfTrackCandidateMaker",
    cleanTrajectoryAfterInOut = cms.bool( False ),
    doSeedingRegionRebuilding = cms.bool( False ),
    onlyPixelHitsForSeedCleaner = cms.bool( False ),
    reverseTrajectories = cms.bool( True ),
    useHitsSplitting = cms.bool( False ),
    MeasurementTrackerEvent = cms.InputTag( "hltSiStripClusters" ),
    src = cms.InputTag( "hltIterL3OISeedsFromL2MuonsNoVtx" ),
    clustersToSkip = cms.InputTag( "" ),
    phase2clustersToSkip = cms.InputTag( "" ),
    TrajectoryBuilderPSet = cms.PSet(  refToPSet_ = cms.string( "HLTPSetMuonCkfTrajectoryBuilder" ) ),
    TransientInitialStateEstimatorParameters = cms.PSet( 
      propagatorAlongTISE = cms.string( "PropagatorWithMaterial" ),
      numberMeasurementsForFit = cms.int32( 4 ),
      propagatorOppositeTISE = cms.string( "PropagatorWithMaterialOpposite" )
    ),
    numHitsForSeedCleaner = cms.int32( 4 ),
    NavigationSchool = cms.string( "SimpleNavigationSchool" ),
    RedundantSeedCleaner = cms.string( "CachingSeedCleanerBySharedInput" ),
    TrajectoryCleaner = cms.string( "muonSeededTrajectoryCleanerBySharedHits" ),
    maxNSeeds = cms.uint32( 500000 ),
    maxSeedsBeforeCleaning = cms.uint32( 5000 )
)
process.hltIterL3OIMuCtfWithMaterialTracksNoVtx = cms.EDProducer( "TrackProducer",
    useSimpleMF = cms.bool( False ),
    SimpleMagneticField = cms.string( "" ),
    src = cms.InputTag( "hltIterL3OITrackCandidatesNoVtx" ),
    clusterRemovalInfo = cms.InputTag( "" ),
    beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
    Fitter = cms.string( "hltESPKFFittingSmootherWithOutliersRejectionAndRK" ),
    useHitsSplitting = cms.bool( False ),
    alias = cms.untracked.string( "ctfWithMaterialTracks" ),
    TrajectoryInEvent = cms.bool( False ),
    TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
    AlgorithmName = cms.string( "iter10" ),
    Propagator = cms.string( "PropagatorWithMaterial" ),
    GeometricInnerState = cms.bool( True ),
    NavigationSchool = cms.string( "SimpleNavigationSchool" ),
    MeasurementTracker = cms.string( "hltESPMeasurementTracker" ),
    MeasurementTrackerEvent = cms.InputTag( "hltSiStripClusters" )
)
process.hltIterL3OIMuonTrackCutClassifierNoVtx = cms.EDProducer( "TrackCutClassifier",
    src = cms.InputTag( "hltIterL3OIMuCtfWithMaterialTracksNoVtx" ),
    beamspot = cms.InputTag( "hltOnlineBeamSpot" ),
    vertices = cms.InputTag( "Notused" ),
    ignoreVertices = cms.bool( True ),
    qualityCuts = cms.vdouble( -0.7, 0.1, 0.7 ),
    mva = cms.PSet( 
      minPixelHits = cms.vint32( 0, 0, 0 ),
      maxDzWrtBS = cms.vdouble( 3.40282346639E38, 24.0, 100.0 ),
      dr_par = cms.PSet( 
        d0err = cms.vdouble( 0.003, 0.003, 3.40282346639E38 ),
        dr_par2 = cms.vdouble( 0.3, 0.3, 3.40282346639E38 ),
        dr_par1 = cms.vdouble( 0.4, 0.4, 3.40282346639E38 ),
        dr_exp = cms.vint32( 4, 4, 2147483647 ),
        d0err_par = cms.vdouble( 0.001, 0.001, 3.40282346639E38 )
      ),
      maxLostLayers = cms.vint32( 4, 3, 2 ),
      min3DLayers = cms.vint32( 0, 0, 0 ),
      dz_par = cms.PSet( 
        dz_par1 = cms.vdouble( 0.4, 0.4, 3.40282346639E38 ),
        dz_par2 = cms.vdouble( 0.35, 0.35, 3.40282346639E38 ),
        dz_exp = cms.vint32( 4, 4, 2147483647 )
      ),
      minNVtxTrk = cms.int32( 3 ),
      maxDz = cms.vdouble( 0.5, 0.2, 3.40282346639E38 ),
      minNdof = cms.vdouble( 1.0E-5, 1.0E-5, 1.0E-5 ),
      maxChi2 = cms.vdouble( 3.40282346639E38, 3.40282346639E38, 3.40282346639E38 ),
      maxChi2n = cms.vdouble( 10.0, 1.0, 0.4 ),
      maxDr = cms.vdouble( 0.5, 0.03, 3.40282346639E38 ),
      minLayers = cms.vint32( 3, 5, 5 )
    )
)
process.hltIterL3OIMuonTrackSelectionHighPurityNoVtx = cms.EDProducer( "TrackCollectionFilterCloner",
    originalSource = cms.InputTag( "hltIterL3OIMuCtfWithMaterialTracksNoVtx" ),
    originalMVAVals = cms.InputTag( 'hltIterL3OIMuonTrackCutClassifierNoVtx','MVAValues' ),
    originalQualVals = cms.InputTag( 'hltIterL3OIMuonTrackCutClassifierNoVtx','QualityMasks' ),
    minQuality = cms.string( "highPurity" ),
    copyExtras = cms.untracked.bool( True ),
    copyTrajectories = cms.untracked.bool( False )
)
process.hltL3MuonsIterL3OINoVtx = cms.EDProducer( "L3MuonProducer",
    ServiceParameters = cms.PSet( 
      RPCLayers = cms.bool( True ),
      UseMuonNavigation = cms.untracked.bool( True ),
      Propagators = cms.untracked.vstring( 'hltESPSmartPropagatorAny',
        'SteppingHelixPropagatorAny',
        'hltESPSmartPropagator',
        'hltESPSteppingHelixPropagatorOpposite' )
    ),
    MuonCollectionLabel = cms.InputTag( "hltL2Muons" ),
    TrackLoaderParameters = cms.PSet( 
      MuonSeededTracksInstance = cms.untracked.string( "L2Seeded" ),
      TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
      beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
      DoSmoothing = cms.bool( True ),
      SmoothTkTrack = cms.untracked.bool( False ),
      VertexConstraint = cms.bool( False ),
      MuonUpdatorAtVertexParameters = cms.PSet( 
        MaxChi2 = cms.double( 1000000.0 ),
        BeamSpotPositionErrors = cms.vdouble( 0.1, 0.1, 5.3 ),
        Propagator = cms.string( "hltESPSteppingHelixPropagatorOpposite" )
      ),
      PutTkTrackIntoEvent = cms.untracked.bool( False ),
      Smoother = cms.string( "hltESPKFTrajectorySmootherForMuonTrackLoader" )
    ),
    L3TrajBuilderParameters = cms.PSet( 
      PtCut = cms.double( 1.0 ),
      TrackerPropagator = cms.string( "SteppingHelixPropagatorAny" ),
      GlobalMuonTrackMatcher = cms.PSet( 
        Chi2Cut_3 = cms.double( 200.0 ),
        DeltaDCut_2 = cms.double( 10.0 ),
        Eta_threshold = cms.double( 1.2 ),
        Quality_2 = cms.double( 15.0 ),
        DeltaDCut_1 = cms.double( 40.0 ),
        Quality_3 = cms.double( 7.0 ),
        DeltaDCut_3 = cms.double( 15.0 ),
        Quality_1 = cms.double( 20.0 ),
        Pt_threshold1 = cms.double( 0.0 ),
        DeltaRCut_2 = cms.double( 0.2 ),
        DeltaRCut_1 = cms.double( 0.1 ),
        Pt_threshold2 = cms.double( 9.99999999E8 ),
        Chi2Cut_1 = cms.double( 50.0 ),
        Chi2Cut_2 = cms.double( 50.0 ),
        DeltaRCut_3 = cms.double( 1.0 ),
        LocChi2Cut = cms.double( 0.001 ),
        Propagator = cms.string( "hltESPSmartPropagator" ),
        MinPt = cms.double( 1.0 ),
        MinP = cms.double( 2.5 )
      ),
      ScaleTECxFactor = cms.double( -1.0 ),
      tkTrajUseVertex = cms.bool( False ),
      MuonTrackingRegionBuilder = cms.PSet( 
        Rescale_Dz = cms.double( 4.0 ),
        Pt_fixed = cms.bool( False ),
        Eta_fixed = cms.bool( True ),
        Eta_min = cms.double( 0.1 ),
        DeltaZ = cms.double( 24.2 ),
        maxRegions = cms.int32( 2 ),
        EtaR_UpperLimit_Par1 = cms.double( 0.25 ),
        UseVertex = cms.bool( False ),
        Z_fixed = cms.bool( False ),
        PhiR_UpperLimit_Par1 = cms.double( 0.6 ),
        PhiR_UpperLimit_Par2 = cms.double( 0.2 ),
        Rescale_phi = cms.double( 3.0 ),
        DeltaEta = cms.double( 0.2 ),
        precise = cms.bool( True ),
        OnDemand = cms.int32( -1 ),
        EtaR_UpperLimit_Par2 = cms.double( 0.15 ),
        MeasurementTrackerName = cms.InputTag( "hltESPMeasurementTracker" ),
        vertexCollection = cms.InputTag( "pixelVertices" ),
        Pt_min = cms.double( 3.0 ),
        beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
        Phi_fixed = cms.bool( True ),
        DeltaR = cms.double( 0.025 ),
        input = cms.InputTag( "hltL2Muons" ),
        DeltaPhi = cms.double( 0.15 ),
        Phi_min = cms.double( 0.1 ),
        Rescale_eta = cms.double( 3.0 )
      ),
      TrackTransformer = cms.PSet( 
        Fitter = cms.string( "hltESPL3MuKFTrajectoryFitter" ),
        RefitDirection = cms.string( "insideOut" ),
        RefitRPCHits = cms.bool( True ),
        Propagator = cms.string( "hltESPSmartPropagatorAny" ),
        DoPredictionsOnly = cms.bool( False ),
        TrackerRecHitBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
        MuonRecHitBuilder = cms.string( "hltESPMuonTransientTrackingRecHitBuilder" ),
        Smoother = cms.string( "hltESPKFTrajectorySmootherForMuonTrackLoader" )
      ),
      tkTrajBeamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
      RefitRPCHits = cms.bool( True ),
      tkTrajVertex = cms.InputTag( "Notused" ),
      GlbRefitterParameters = cms.PSet( 
        Fitter = cms.string( "hltESPL3MuKFTrajectoryFitter" ),
        DTRecSegmentLabel = cms.InputTag( "hltDt4DSegments" ),
        RefitFlag = cms.bool( True ),
        SkipStation = cms.int32( -1 ),
        Chi2CutRPC = cms.double( 1.0 ),
        PropDirForCosmics = cms.bool( False ),
        CSCRecSegmentLabel = cms.InputTag( "hltCscSegments" ),
        HitThreshold = cms.int32( 1 ),
        DYTthrs = cms.vint32( 30, 15 ),
        TrackerSkipSystem = cms.int32( -1 ),
        RefitDirection = cms.string( "insideOut" ),
        Chi2CutCSC = cms.double( 150.0 ),
        Chi2CutDT = cms.double( 10.0 ),
        RefitRPCHits = cms.bool( True ),
        TrackerSkipSection = cms.int32( -1 ),
        Propagator = cms.string( "hltESPSmartPropagatorAny" ),
        DoPredictionsOnly = cms.bool( False ),
        TrackerRecHitBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
        MuonHitsOption = cms.int32( 1 ),
        MuonRecHitBuilder = cms.string( "hltESPMuonTransientTrackingRecHitBuilder" )
      ),
      PCut = cms.double( 2.5 ),
      tkTrajMaxDXYBeamSpot = cms.double( 9999.0 ),
      TrackerRecHitBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
      tkTrajMaxChi2 = cms.double( 9999.0 ),
      MuonRecHitBuilder = cms.string( "hltESPMuonTransientTrackingRecHitBuilder" ),
      ScaleTECyFactor = cms.double( -1.0 ),
      tkTrajLabel = cms.InputTag( "hltIterL3OIMuonTrackSelectionHighPurityNoVtx" )
    )
)
process.hltIterL3OIL3MuonsLinksCombinationNoVtx = cms.EDProducer( "L3TrackLinksCombiner",
    labels = cms.VInputTag( 'hltL3MuonsIterL3OINoVtx' )
)
process.hltIterL3OIL3MuonsNoVtx = cms.EDProducer( "L3TrackCombiner",
    labels = cms.VInputTag( 'hltL3MuonsIterL3OINoVtx' )
)
process.hltIterL3OIL3MuonCandidatesNoVtx = cms.EDProducer( "L3MuonCandidateProducer",
    InputObjects = cms.InputTag( "hltIterL3OIL3MuonsNoVtx" ),
    InputLinksObjects = cms.InputTag( "hltIterL3OIL3MuonsLinksCombinationNoVtx" ),
    MuonPtOption = cms.string( "Tracker" )
)
process.hltL2SelectorForL3IONoVtx = cms.EDProducer( "HLTMuonL2SelectorForL3IO",
    l2Src = cms.InputTag( "hltL2Muons" ),
    l3OISrc = cms.InputTag( "hltIterL3OIL3MuonCandidatesNoVtx" ),
    InputLinks = cms.InputTag( "hltIterL3OIL3MuonsLinksCombinationNoVtx" ),
    applyL3Filters = cms.bool( False ),
    MinNhits = cms.int32( 1 ),
    MaxNormalizedChi2 = cms.double( 20.0 ),
    MinNmuonHits = cms.int32( 1 ),
    MaxPtDifference = cms.double( 0.3 )
)
process.hltIterL3MuonPixelTracksTrackingRegionsNoVtx = cms.EDProducer( "MuonTrackingRegionEDProducer",
    EtaR_UpperLimit_Par1 = cms.double( 0.25 ),
    DeltaR = cms.double( 0.025 ),
    beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
    OnDemand = cms.int32( -1 ),
    vertexCollection = cms.InputTag( "notUsed" ),
    Rescale_phi = cms.double( 3.0 ),
    Eta_fixed = cms.bool( True ),
    Rescale_eta = cms.double( 3.0 ),
    PhiR_UpperLimit_Par2 = cms.double( 0.2 ),
    Eta_min = cms.double( 0.0 ),
    Phi_fixed = cms.bool( True ),
    Phi_min = cms.double( 0.0 ),
    PhiR_UpperLimit_Par1 = cms.double( 0.6 ),
    EtaR_UpperLimit_Par2 = cms.double( 0.15 ),
    MeasurementTrackerName = cms.InputTag( "" ),
    UseVertex = cms.bool( False ),
    Rescale_Dz = cms.double( 4.0 ),
    Pt_fixed = cms.bool( True ),
    Z_fixed = cms.bool( True ),
    Pt_min = cms.double( 0.0 ),
    DeltaZ = cms.double( 24.2 ),
    DeltaEta = cms.double( 0.2 ),
    DeltaPhi = cms.double( 0.15 ),
    maxRegions = cms.int32( 5 ),
    precise = cms.bool( True ),
    input = cms.InputTag( "hltL2SelectorForL3IONoVtx" )
)
process.hltPixelTracksInRegionL2NoVtx = cms.EDProducer( "TrackSelectorByRegion",
    tracks = cms.InputTag( "hltPixelTracks" ),
    regions = cms.InputTag( "hltIterL3MuonPixelTracksTrackingRegionsNoVtx" ),
    produceTrackCollection = cms.bool( True ),
    produceMask = cms.bool( False )
)
process.hltIter0IterL3MuonPixelSeedsFromPixelTracksNoVtx = cms.EDProducer( "SeedGeneratorFromProtoTracksEDProducer",
    InputCollection = cms.InputTag( "hltPixelTracksInRegionL2NoVtx" ),
    InputVertexCollection = cms.InputTag( "" ),
    originHalfLength = cms.double( 0.3 ),
    originRadius = cms.double( 0.1 ),
    useProtoTrackKinematics = cms.bool( False ),
    useEventsWithNoVertex = cms.bool( True ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    usePV = cms.bool( False ),
    includeFourthHit = cms.bool( True ),
    SeedCreatorPSet = cms.PSet(  refToPSet_ = cms.string( "HLTSeedFromProtoTracks" ) )
)
process.hltIter0IterL3MuonCkfTrackCandidatesNoVtx = cms.EDProducer( "CkfTrackCandidateMaker",
    cleanTrajectoryAfterInOut = cms.bool( False ),
    doSeedingRegionRebuilding = cms.bool( True ),
    onlyPixelHitsForSeedCleaner = cms.bool( False ),
    reverseTrajectories = cms.bool( False ),
    useHitsSplitting = cms.bool( True ),
    MeasurementTrackerEvent = cms.InputTag( "hltSiStripClusters" ),
    src = cms.InputTag( "hltIter0IterL3MuonPixelSeedsFromPixelTracksNoVtx" ),
    clustersToSkip = cms.InputTag( "" ),
    phase2clustersToSkip = cms.InputTag( "" ),
    TrajectoryBuilderPSet = cms.PSet(  refToPSet_ = cms.string( "HLTIter0IterL3MuonPSetGroupedCkfTrajectoryBuilderIT" ) ),
    TransientInitialStateEstimatorParameters = cms.PSet( 
      propagatorAlongTISE = cms.string( "PropagatorWithMaterialParabolicMf" ),
      numberMeasurementsForFit = cms.int32( 4 ),
      propagatorOppositeTISE = cms.string( "PropagatorWithMaterialParabolicMfOpposite" )
    ),
    numHitsForSeedCleaner = cms.int32( 4 ),
    NavigationSchool = cms.string( "SimpleNavigationSchool" ),
    RedundantSeedCleaner = cms.string( "none" ),
    TrajectoryCleaner = cms.string( "hltESPTrajectoryCleanerBySharedHits" ),
    maxNSeeds = cms.uint32( 100000 ),
    maxSeedsBeforeCleaning = cms.uint32( 1000 )
)
process.hltIter0IterL3MuonCtfWithMaterialTracksNoVtx = cms.EDProducer( "TrackProducer",
    useSimpleMF = cms.bool( True ),
    SimpleMagneticField = cms.string( "ParabolicMf" ),
    src = cms.InputTag( "hltIter0IterL3MuonCkfTrackCandidatesNoVtx" ),
    clusterRemovalInfo = cms.InputTag( "" ),
    beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
    Fitter = cms.string( "hltESPFittingSmootherIT" ),
    useHitsSplitting = cms.bool( False ),
    alias = cms.untracked.string( "ctfWithMaterialTracks" ),
    TrajectoryInEvent = cms.bool( False ),
    TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
    AlgorithmName = cms.string( "hltIter0" ),
    Propagator = cms.string( "hltESPRungeKuttaTrackerPropagator" ),
    GeometricInnerState = cms.bool( True ),
    NavigationSchool = cms.string( "" ),
    MeasurementTracker = cms.string( "" ),
    MeasurementTrackerEvent = cms.InputTag( "hltSiStripClusters" )
)
process.hltIter0IterL3MuonTrackCutClassifierNoVtx = cms.EDProducer( "TrackCutClassifier",
    src = cms.InputTag( "hltIter0IterL3MuonCtfWithMaterialTracksNoVtx" ),
    beamspot = cms.InputTag( "hltOnlineBeamSpot" ),
    vertices = cms.InputTag( "hltTrimmedPixelVertices" ),
    ignoreVertices = cms.bool( False ),
    qualityCuts = cms.vdouble( -0.7, 0.1, 0.7 ),
    mva = cms.PSet( 
      minPixelHits = cms.vint32( 0, 0, 0 ),
      maxDzWrtBS = cms.vdouble( 3.40282346639E38, 24.0, 100.0 ),
      dr_par = cms.PSet( 
        d0err = cms.vdouble( 0.003, 0.003, 3.40282346639E38 ),
        dr_par2 = cms.vdouble( 0.3, 0.3, 3.40282346639E38 ),
        dr_par1 = cms.vdouble( 0.4, 0.4, 3.40282346639E38 ),
        dr_exp = cms.vint32( 4, 4, 2147483647 ),
        d0err_par = cms.vdouble( 0.001, 0.001, 3.40282346639E38 )
      ),
      maxLostLayers = cms.vint32( 1, 1, 1 ),
      min3DLayers = cms.vint32( 0, 0, 0 ),
      dz_par = cms.PSet( 
        dz_par1 = cms.vdouble( 0.4, 0.4, 3.40282346639E38 ),
        dz_par2 = cms.vdouble( 0.35, 0.35, 3.40282346639E38 ),
        dz_exp = cms.vint32( 4, 4, 2147483647 )
      ),
      minNVtxTrk = cms.int32( 3 ),
      maxDz = cms.vdouble( 0.5, 0.2, 3.40282346639E38 ),
      minNdof = cms.vdouble( 1.0E-5, 1.0E-5, 1.0E-5 ),
      maxChi2 = cms.vdouble( 3.40282346639E38, 3.40282346639E38, 3.40282346639E38 ),
      maxChi2n = cms.vdouble( 1.2, 1.0, 0.7 ),
      maxDr = cms.vdouble( 0.5, 0.03, 3.40282346639E38 ),
      minLayers = cms.vint32( 3, 3, 4 )
    )
)
process.hltIter0IterL3MuonTrackSelectionHighPurityNoVtx = cms.EDProducer( "TrackCollectionFilterCloner",
    originalSource = cms.InputTag( "hltIter0IterL3MuonCtfWithMaterialTracksNoVtx" ),
    originalMVAVals = cms.InputTag( 'hltIter0IterL3MuonTrackCutClassifierNoVtx','MVAValues' ),
    originalQualVals = cms.InputTag( 'hltIter0IterL3MuonTrackCutClassifierNoVtx','QualityMasks' ),
    minQuality = cms.string( "highPurity" ),
    copyExtras = cms.untracked.bool( True ),
    copyTrajectories = cms.untracked.bool( False )
)
process.hltL3MuonsIterL3IONoVtx = cms.EDProducer( "L3MuonProducer",
    ServiceParameters = cms.PSet( 
      RPCLayers = cms.bool( True ),
      UseMuonNavigation = cms.untracked.bool( True ),
      Propagators = cms.untracked.vstring( 'hltESPSmartPropagatorAny',
        'SteppingHelixPropagatorAny',
        'hltESPSmartPropagator',
        'hltESPSteppingHelixPropagatorOpposite' )
    ),
    MuonCollectionLabel = cms.InputTag( "hltL2Muons" ),
    TrackLoaderParameters = cms.PSet( 
      MuonSeededTracksInstance = cms.untracked.string( "L2Seeded" ),
      beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
      DoSmoothing = cms.bool( False ),
      SmoothTkTrack = cms.untracked.bool( False ),
      VertexConstraint = cms.bool( False ),
      MuonUpdatorAtVertexParameters = cms.PSet( 
        MaxChi2 = cms.double( 1000000.0 ),
        BeamSpotPositionErrors = cms.vdouble( 0.1, 0.1, 5.3 ),
        Propagator = cms.string( "hltESPSteppingHelixPropagatorOpposite" )
      ),
      PutTkTrackIntoEvent = cms.untracked.bool( False ),
      Smoother = cms.string( "hltESPKFTrajectorySmootherForMuonTrackLoader" )
    ),
    L3TrajBuilderParameters = cms.PSet( 
      PtCut = cms.double( 1.0 ),
      TrackerPropagator = cms.string( "SteppingHelixPropagatorAny" ),
      GlobalMuonTrackMatcher = cms.PSet( 
        Chi2Cut_3 = cms.double( 200.0 ),
        DeltaDCut_2 = cms.double( 10.0 ),
        Eta_threshold = cms.double( 1.2 ),
        Quality_2 = cms.double( 15.0 ),
        DeltaDCut_1 = cms.double( 40.0 ),
        Quality_3 = cms.double( 7.0 ),
        DeltaDCut_3 = cms.double( 15.0 ),
        Quality_1 = cms.double( 20.0 ),
        Pt_threshold1 = cms.double( 0.0 ),
        DeltaRCut_2 = cms.double( 0.2 ),
        DeltaRCut_1 = cms.double( 0.1 ),
        Pt_threshold2 = cms.double( 9.99999999E8 ),
        Chi2Cut_1 = cms.double( 50.0 ),
        Chi2Cut_2 = cms.double( 50.0 ),
        DeltaRCut_3 = cms.double( 1.0 ),
        LocChi2Cut = cms.double( 0.001 ),
        Propagator = cms.string( "hltESPSmartPropagator" ),
        MinPt = cms.double( 1.0 ),
        MinP = cms.double( 2.5 )
      ),
      ScaleTECxFactor = cms.double( -1.0 ),
      tkTrajUseVertex = cms.bool( False ),
      MuonTrackingRegionBuilder = cms.PSet( 
        Rescale_Dz = cms.double( 4.0 ),
        Pt_fixed = cms.bool( True ),
        Eta_fixed = cms.bool( True ),
        Eta_min = cms.double( 0.1 ),
        DeltaZ = cms.double( 24.2 ),
        maxRegions = cms.int32( 2 ),
        EtaR_UpperLimit_Par1 = cms.double( 0.25 ),
        UseVertex = cms.bool( False ),
        Z_fixed = cms.bool( True ),
        PhiR_UpperLimit_Par1 = cms.double( 0.6 ),
        PhiR_UpperLimit_Par2 = cms.double( 0.2 ),
        Rescale_phi = cms.double( 3.0 ),
        DeltaEta = cms.double( 0.04 ),
        precise = cms.bool( True ),
        OnDemand = cms.int32( -1 ),
        EtaR_UpperLimit_Par2 = cms.double( 0.15 ),
        MeasurementTrackerName = cms.InputTag( "hltESPMeasurementTracker" ),
        vertexCollection = cms.InputTag( "pixelVertices" ),
        Pt_min = cms.double( 3.0 ),
        beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
        Phi_fixed = cms.bool( True ),
        DeltaR = cms.double( 0.025 ),
        input = cms.InputTag( "hltL2SelectorForL3IONoVtx" ),
        DeltaPhi = cms.double( 0.15 ),
        Phi_min = cms.double( 0.1 ),
        Rescale_eta = cms.double( 3.0 )
      ),
      TrackTransformer = cms.PSet( 
        Fitter = cms.string( "hltESPL3MuKFTrajectoryFitter" ),
        RefitDirection = cms.string( "insideOut" ),
        RefitRPCHits = cms.bool( True ),
        Propagator = cms.string( "hltESPSmartPropagatorAny" ),
        DoPredictionsOnly = cms.bool( False ),
        TrackerRecHitBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
        MuonRecHitBuilder = cms.string( "hltESPMuonTransientTrackingRecHitBuilder" ),
        Smoother = cms.string( "hltESPKFTrajectorySmootherForMuonTrackLoader" )
      ),
      tkTrajBeamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
      RefitRPCHits = cms.bool( True ),
      tkTrajVertex = cms.InputTag( "hltTrimmedPixelVertices" ),
      GlbRefitterParameters = cms.PSet( 
        Fitter = cms.string( "hltESPL3MuKFTrajectoryFitter" ),
        DTRecSegmentLabel = cms.InputTag( "hltDt4DSegments" ),
        RefitFlag = cms.bool( True ),
        SkipStation = cms.int32( -1 ),
        Chi2CutRPC = cms.double( 1.0 ),
        PropDirForCosmics = cms.bool( False ),
        CSCRecSegmentLabel = cms.InputTag( "hltCscSegments" ),
        HitThreshold = cms.int32( 1 ),
        DYTthrs = cms.vint32( 30, 15 ),
        TrackerSkipSystem = cms.int32( -1 ),
        RefitDirection = cms.string( "insideOut" ),
        Chi2CutCSC = cms.double( 150.0 ),
        Chi2CutDT = cms.double( 10.0 ),
        RefitRPCHits = cms.bool( True ),
        TrackerSkipSection = cms.int32( -1 ),
        Propagator = cms.string( "hltESPSmartPropagatorAny" ),
        DoPredictionsOnly = cms.bool( False ),
        TrackerRecHitBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
        MuonHitsOption = cms.int32( 1 ),
        MuonRecHitBuilder = cms.string( "hltESPMuonTransientTrackingRecHitBuilder" )
      ),
      PCut = cms.double( 2.5 ),
      tkTrajMaxDXYBeamSpot = cms.double( 9999.0 ),
      TrackerRecHitBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
      matchToSeeds = cms.bool( True ),
      tkTrajMaxChi2 = cms.double( 9999.0 ),
      MuonRecHitBuilder = cms.string( "hltESPMuonTransientTrackingRecHitBuilder" ),
      ScaleTECyFactor = cms.double( -1.0 ),
      tkTrajLabel = cms.InputTag( "hltIter0IterL3MuonTrackSelectionHighPurityNoVtx" )
    )
)
process.hltIterL3MuonsFromL2LinksCombinationNoVtx = cms.EDProducer( "L3TrackLinksCombiner",
    labels = cms.VInputTag( 'hltL3MuonsIterL3OINoVtx','hltL3MuonsIterL3IONoVtx' )
)
process.hltIterL3MuonsFromL2NoVtx = cms.EDProducer( "L3TrackCombiner",
    labels = cms.VInputTag( 'hltL3MuonsIterL3OINoVtx','hltL3MuonsIterL3IONoVtx' )
)
process.hltIterL3MuonL1MuonNoL2SelectorNoVtx = cms.EDProducer( "HLTL1MuonNoL2Selector",
    InputObjects = cms.InputTag( 'hltGtStage2Digis','Muon' ),
    L2CandTag = cms.InputTag( "hltL2MuonCandidatesNoVtx" ),
    SeedMapTag = cms.InputTag( "hltL2Muons" ),
    L1MinPt = cms.double( -1.0 ),
    L1MaxEta = cms.double( 5.0 ),
    L1MinQuality = cms.uint32( 7 ),
    CentralBxOnly = cms.bool( True )
)
process.hltIterL3FromL1MuonPixelTracksTrackingRegionsNoVtx = cms.EDProducer( "CandidateSeededTrackingRegionsEDProducer",
    RegionPSet = cms.PSet( 
      vertexCollection = cms.InputTag( "notUsed" ),
      zErrorVetex = cms.double( 0.2 ),
      beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
      zErrorBeamSpot = cms.double( 24.2 ),
      maxNVertices = cms.int32( 1 ),
      maxNRegions = cms.int32( 5 ),
      nSigmaZVertex = cms.double( 3.0 ),
      nSigmaZBeamSpot = cms.double( 4.0 ),
      ptMin = cms.double( 0.0 ),
      mode = cms.string( "BeamSpotSigma" ),
      input = cms.InputTag( "hltIterL3MuonL1MuonNoL2SelectorNoVtx" ),
      searchOpt = cms.bool( False ),
      whereToUseMeasurementTracker = cms.string( "Never" ),
      originRadius = cms.double( 0.2 ),
      measurementTrackerName = cms.InputTag( "" ),
      precise = cms.bool( True ),
      deltaEta = cms.double( 0.35 ),
      deltaPhi = cms.double( 0.2 )
    )
)
process.hltPixelTracksInRegionL1NoVtx = cms.EDProducer( "TrackSelectorByRegion",
    tracks = cms.InputTag( "hltPixelTracks" ),
    regions = cms.InputTag( "hltIterL3FromL1MuonPixelTracksTrackingRegionsNoVtx" ),
    produceTrackCollection = cms.bool( True ),
    produceMask = cms.bool( False )
)
process.hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracksNoVtx = cms.EDProducer( "SeedGeneratorFromProtoTracksEDProducer",
    InputCollection = cms.InputTag( "hltPixelTracksInRegionL1NoVtx" ),
    InputVertexCollection = cms.InputTag( "" ),
    originHalfLength = cms.double( 0.3 ),
    originRadius = cms.double( 0.1 ),
    useProtoTrackKinematics = cms.bool( False ),
    useEventsWithNoVertex = cms.bool( True ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    usePV = cms.bool( False ),
    includeFourthHit = cms.bool( True ),
    SeedCreatorPSet = cms.PSet(  refToPSet_ = cms.string( "HLTSeedFromProtoTracks" ) )
)
process.hltIter0IterL3FromL1MuonCkfTrackCandidatesNoVtx = cms.EDProducer( "CkfTrackCandidateMaker",
    cleanTrajectoryAfterInOut = cms.bool( False ),
    doSeedingRegionRebuilding = cms.bool( True ),
    onlyPixelHitsForSeedCleaner = cms.bool( False ),
    reverseTrajectories = cms.bool( False ),
    useHitsSplitting = cms.bool( True ),
    MeasurementTrackerEvent = cms.InputTag( "hltSiStripClusters" ),
    src = cms.InputTag( "hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracksNoVtx" ),
    clustersToSkip = cms.InputTag( "" ),
    phase2clustersToSkip = cms.InputTag( "" ),
    TrajectoryBuilderPSet = cms.PSet(  refToPSet_ = cms.string( "HLTIter0IterL3FromL1MuonPSetGroupedCkfTrajectoryBuilderIT" ) ),
    TransientInitialStateEstimatorParameters = cms.PSet( 
      propagatorAlongTISE = cms.string( "PropagatorWithMaterialParabolicMf" ),
      numberMeasurementsForFit = cms.int32( 4 ),
      propagatorOppositeTISE = cms.string( "PropagatorWithMaterialParabolicMfOpposite" )
    ),
    numHitsForSeedCleaner = cms.int32( 4 ),
    NavigationSchool = cms.string( "SimpleNavigationSchool" ),
    RedundantSeedCleaner = cms.string( "none" ),
    TrajectoryCleaner = cms.string( "hltESPTrajectoryCleanerBySharedHits" ),
    maxNSeeds = cms.uint32( 100000 ),
    maxSeedsBeforeCleaning = cms.uint32( 1000 )
)
process.hltIter0IterL3FromL1MuonCtfWithMaterialTracksNoVtx = cms.EDProducer( "TrackProducer",
    useSimpleMF = cms.bool( True ),
    SimpleMagneticField = cms.string( "ParabolicMf" ),
    src = cms.InputTag( "hltIter0IterL3FromL1MuonCkfTrackCandidatesNoVtx" ),
    clusterRemovalInfo = cms.InputTag( "" ),
    beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
    Fitter = cms.string( "hltESPFittingSmootherIT" ),
    useHitsSplitting = cms.bool( False ),
    alias = cms.untracked.string( "ctfWithMaterialTracks" ),
    TrajectoryInEvent = cms.bool( False ),
    TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
    AlgorithmName = cms.string( "hltIter0" ),
    Propagator = cms.string( "hltESPRungeKuttaTrackerPropagator" ),
    GeometricInnerState = cms.bool( True ),
    NavigationSchool = cms.string( "" ),
    MeasurementTracker = cms.string( "" ),
    MeasurementTrackerEvent = cms.InputTag( "hltSiStripClusters" )
)
process.hltIter0IterL3FromL1MuonTrackCutClassifierNoVtx = cms.EDProducer( "TrackCutClassifier",
    src = cms.InputTag( "hltIter0IterL3FromL1MuonCtfWithMaterialTracksNoVtx" ),
    beamspot = cms.InputTag( "hltOnlineBeamSpot" ),
    vertices = cms.InputTag( "hltTrimmedPixelVertices" ),
    ignoreVertices = cms.bool( False ),
    qualityCuts = cms.vdouble( -0.7, 0.1, 0.7 ),
    mva = cms.PSet( 
      minPixelHits = cms.vint32( 0, 0, 0 ),
      maxDzWrtBS = cms.vdouble( 3.40282346639E38, 24.0, 100.0 ),
      dr_par = cms.PSet( 
        d0err = cms.vdouble( 0.003, 0.003, 3.40282346639E38 ),
        dr_par2 = cms.vdouble( 0.3, 0.3, 3.40282346639E38 ),
        dr_par1 = cms.vdouble( 0.4, 0.4, 3.40282346639E38 ),
        dr_exp = cms.vint32( 4, 4, 2147483647 ),
        d0err_par = cms.vdouble( 0.001, 0.001, 3.40282346639E38 )
      ),
      maxLostLayers = cms.vint32( 1, 1, 1 ),
      min3DLayers = cms.vint32( 0, 0, 0 ),
      dz_par = cms.PSet( 
        dz_par1 = cms.vdouble( 0.4, 0.4, 3.40282346639E38 ),
        dz_par2 = cms.vdouble( 0.35, 0.35, 3.40282346639E38 ),
        dz_exp = cms.vint32( 4, 4, 2147483647 )
      ),
      minNVtxTrk = cms.int32( 3 ),
      maxDz = cms.vdouble( 0.5, 0.2, 3.40282346639E38 ),
      minNdof = cms.vdouble( 1.0E-5, 1.0E-5, 1.0E-5 ),
      maxChi2 = cms.vdouble( 3.40282346639E38, 3.40282346639E38, 3.40282346639E38 ),
      maxChi2n = cms.vdouble( 1.2, 1.0, 0.7 ),
      maxDr = cms.vdouble( 0.5, 0.03, 3.40282346639E38 ),
      minLayers = cms.vint32( 3, 3, 4 )
    )
)
process.hltIter0IterL3FromL1MuonTrackSelectionHighPurityNoVtx = cms.EDProducer( "TrackCollectionFilterCloner",
    originalSource = cms.InputTag( "hltIter0IterL3FromL1MuonCtfWithMaterialTracksNoVtx" ),
    originalMVAVals = cms.InputTag( 'hltIter0IterL3FromL1MuonTrackCutClassifierNoVtx','MVAValues' ),
    originalQualVals = cms.InputTag( 'hltIter0IterL3FromL1MuonTrackCutClassifierNoVtx','QualityMasks' ),
    minQuality = cms.string( "highPurity" ),
    copyExtras = cms.untracked.bool( True ),
    copyTrajectories = cms.untracked.bool( False )
)
process.hltIterL3MuonMergedNoVtx = cms.EDProducer( "TrackListMerger",
    ShareFrac = cms.double( 0.19 ),
    FoundHitBonus = cms.double( 5.0 ),
    LostHitPenalty = cms.double( 20.0 ),
    MinPT = cms.double( 0.05 ),
    Epsilon = cms.double( -0.001 ),
    MaxNormalizedChisq = cms.double( 1000.0 ),
    MinFound = cms.int32( 3 ),
    TrackProducers = cms.VInputTag( 'hltIterL3OIMuonTrackSelectionHighPurityNoVtx','hltIter0IterL3MuonTrackSelectionHighPurityNoVtx' ),
    hasSelector = cms.vint32( 0, 0 ),
    indivShareFrac = cms.vdouble( 1.0, 1.0 ),
    selectedTrackQuals = cms.VInputTag( 'hltIterL3OIMuonTrackSelectionHighPurityNoVtx','hltIter0IterL3MuonTrackSelectionHighPurityNoVtx' ),
    setsToMerge = cms.VPSet( 
      cms.PSet(  pQual = cms.bool( False ),
        tLists = cms.vint32( 0, 1 )
      )
    ),
    trackAlgoPriorityOrder = cms.string( "hltESPTrackAlgoPriorityOrder" ),
    allowFirstHitShare = cms.bool( True ),
    newQuality = cms.string( "confirmed" ),
    copyExtras = cms.untracked.bool( True ),
    writeOnlyTrkQuals = cms.bool( False ),
    copyMVA = cms.bool( False )
)
process.hltIterL3MuonAndMuonFromL1MergedNoVtx = cms.EDProducer( "TrackListMerger",
    ShareFrac = cms.double( 0.19 ),
    FoundHitBonus = cms.double( 5.0 ),
    LostHitPenalty = cms.double( 20.0 ),
    MinPT = cms.double( 0.05 ),
    Epsilon = cms.double( -0.001 ),
    MaxNormalizedChisq = cms.double( 1000.0 ),
    MinFound = cms.int32( 3 ),
    TrackProducers = cms.VInputTag( 'hltIterL3MuonMergedNoVtx','hltIter0IterL3FromL1MuonTrackSelectionHighPurityNoVtx' ),
    hasSelector = cms.vint32( 0, 0 ),
    indivShareFrac = cms.vdouble( 1.0, 1.0 ),
    selectedTrackQuals = cms.VInputTag( 'hltIterL3MuonMergedNoVtx','hltIter0IterL3FromL1MuonTrackSelectionHighPurityNoVtx' ),
    setsToMerge = cms.VPSet( 
      cms.PSet(  pQual = cms.bool( False ),
        tLists = cms.vint32( 0, 1 )
      )
    ),
    trackAlgoPriorityOrder = cms.string( "hltESPTrackAlgoPriorityOrder" ),
    allowFirstHitShare = cms.bool( True ),
    newQuality = cms.string( "confirmed" ),
    copyExtras = cms.untracked.bool( True ),
    writeOnlyTrkQuals = cms.bool( False ),
    copyMVA = cms.bool( False )
)
process.hltL3MuonsIterL3LinksNoVtx = cms.EDProducer( "MuonLinksProducerForHLT",
    InclusiveTrackerTrackCollection = cms.InputTag( "hltIterL3MuonAndMuonFromL1MergedNoVtx" ),
    LinkCollection = cms.InputTag( "hltIterL3MuonsFromL2LinksCombinationNoVtx" ),
    ptMin = cms.double( 2.5 ),
    pMin = cms.double( 2.5 ),
    shareHitFraction = cms.double( 0.19 )
)
process.hltIterL3MuonsNoVtx = cms.EDProducer( "MuonIdProducer",
    MuonCaloCompatibility = cms.PSet( 
      delta_eta = cms.double( 0.02 ),
      delta_phi = cms.double( 0.02 ),
      allSiPMHO = cms.bool( False ),
      MuonTemplateFileName = cms.FileInPath( "RecoMuon/MuonIdentification/data/MuID_templates_muons_lowPt_3_1_norm.root" ),
      PionTemplateFileName = cms.FileInPath( "RecoMuon/MuonIdentification/data/MuID_templates_pions_lowPt_3_1_norm.root" )
    ),
    TrackAssociatorParameters = cms.PSet( 
      useMuon = cms.bool( True ),
      truthMatch = cms.bool( False ),
      usePreshower = cms.bool( False ),
      dRPreshowerPreselection = cms.double( 0.2 ),
      muonMaxDistanceSigmaY = cms.double( 0.0 ),
      useEcal = cms.bool( False ),
      muonMaxDistanceSigmaX = cms.double( 0.0 ),
      dRMuon = cms.double( 9999.0 ),
      dREcal = cms.double( 9999.0 ),
      CSCSegmentCollectionLabel = cms.InputTag( "hltCscSegments" ),
      DTRecSegment4DCollectionLabel = cms.InputTag( "hltDt4DSegments" ),
      EBRecHitCollectionLabel = cms.InputTag( "Notused" ),
      CaloTowerCollectionLabel = cms.InputTag( "Notused" ),
      propagateAllDirections = cms.bool( True ),
      muonMaxDistanceY = cms.double( 5.0 ),
      useHO = cms.bool( False ),
      muonMaxDistanceX = cms.double( 5.0 ),
      trajectoryUncertaintyTolerance = cms.double( -1.0 ),
      useHcal = cms.bool( False ),
      HBHERecHitCollectionLabel = cms.InputTag( "Notused" ),
      accountForTrajectoryChangeCalo = cms.bool( False ),
      dREcalPreselection = cms.double( 0.05 ),
      useCalo = cms.bool( False ),
      dRMuonPreselection = cms.double( 0.2 ),
      EERecHitCollectionLabel = cms.InputTag( "Notused" ),
      dRHcal = cms.double( 9999.0 ),
      dRHcalPreselection = cms.double( 0.2 ),
      HORecHitCollectionLabel = cms.InputTag( "Notused" )
    ),
    CaloExtractorPSet = cms.PSet( 
      DR_Veto_H = cms.double( 0.1 ),
      CenterConeOnCalIntersection = cms.bool( False ),
      NoiseTow_EE = cms.double( 0.15 ),
      Noise_EB = cms.double( 0.025 ),
      Noise_HE = cms.double( 0.2 ),
      DR_Veto_E = cms.double( 0.07 ),
      NoiseTow_EB = cms.double( 0.04 ),
      Noise_EE = cms.double( 0.1 ),
      UseRecHitsFlag = cms.bool( False ),
      DR_Max = cms.double( 1.0 ),
      DepositLabel = cms.untracked.string( "Cal" ),
      Noise_HO = cms.double( 0.2 ),
      DR_Veto_HO = cms.double( 0.1 ),
      Threshold_H = cms.double( 0.5 ),
      PrintTimeReport = cms.untracked.bool( False ),
      Threshold_E = cms.double( 0.2 ),
      PropagatorName = cms.string( "hltESPFastSteppingHelixPropagatorAny" ),
      ComponentName = cms.string( "CaloExtractorByAssociator" ),
      Threshold_HO = cms.double( 0.5 ),
      DepositInstanceLabels = cms.vstring( 'ecal',
        'hcal',
        'ho' ),
      ServiceParameters = cms.PSet( 
        RPCLayers = cms.bool( False ),
        UseMuonNavigation = cms.untracked.bool( False ),
        Propagators = cms.untracked.vstring( 'hltESPFastSteppingHelixPropagatorAny' )
      ),
      TrackAssociatorParameters = cms.PSet( 
        useMuon = cms.bool( False ),
        truthMatch = cms.bool( False ),
        usePreshower = cms.bool( False ),
        dRPreshowerPreselection = cms.double( 0.2 ),
        muonMaxDistanceSigmaY = cms.double( 0.0 ),
        useEcal = cms.bool( False ),
        muonMaxDistanceSigmaX = cms.double( 0.0 ),
        dRMuon = cms.double( 9999.0 ),
        dREcal = cms.double( 1.0 ),
        CSCSegmentCollectionLabel = cms.InputTag( "hltCscSegments" ),
        DTRecSegment4DCollectionLabel = cms.InputTag( "hltDt4DSegments" ),
        EBRecHitCollectionLabel = cms.InputTag( "Notused" ),
        CaloTowerCollectionLabel = cms.InputTag( "Notused" ),
        propagateAllDirections = cms.bool( True ),
        muonMaxDistanceY = cms.double( 5.0 ),
        useHO = cms.bool( False ),
        muonMaxDistanceX = cms.double( 5.0 ),
        trajectoryUncertaintyTolerance = cms.double( -1.0 ),
        useHcal = cms.bool( False ),
        HBHERecHitCollectionLabel = cms.InputTag( "Notused" ),
        accountForTrajectoryChangeCalo = cms.bool( False ),
        dREcalPreselection = cms.double( 1.0 ),
        useCalo = cms.bool( True ),
        dRMuonPreselection = cms.double( 0.2 ),
        EERecHitCollectionLabel = cms.InputTag( "Notused" ),
        dRHcal = cms.double( 1.0 ),
        dRHcalPreselection = cms.double( 1.0 ),
        HORecHitCollectionLabel = cms.InputTag( "Notused" )
      ),
      Noise_HB = cms.double( 0.2 )
    ),
    TrackExtractorPSet = cms.PSet( 
      Diff_z = cms.double( 0.2 ),
      inputTrackCollection = cms.InputTag( "hltIter0IterL3FromL1MuonTrackSelectionHighPurityNoVtx" ),
      Chi2Ndof_Max = cms.double( 1.0E64 ),
      BeamSpotLabel = cms.InputTag( "hltOnlineBeamSpot" ),
      DR_Veto = cms.double( 0.01 ),
      Pt_Min = cms.double( -1.0 ),
      DR_Max = cms.double( 1.0 ),
      NHits_Min = cms.uint32( 0 ),
      Chi2Prob_Min = cms.double( -1.0 ),
      Diff_r = cms.double( 0.1 ),
      BeamlineOption = cms.string( "BeamSpotFromEvent" ),
      ComponentName = cms.string( "TrackExtractor" )
    ),
    JetExtractorPSet = cms.PSet( 
      JetCollectionLabel = cms.InputTag( "Notused" ),
      DR_Veto = cms.double( 0.1 ),
      DR_Max = cms.double( 1.0 ),
      ExcludeMuonVeto = cms.bool( True ),
      PrintTimeReport = cms.untracked.bool( False ),
      PropagatorName = cms.string( "hltESPFastSteppingHelixPropagatorAny" ),
      ComponentName = cms.string( "JetExtractor" ),
      ServiceParameters = cms.PSet( 
        RPCLayers = cms.bool( False ),
        UseMuonNavigation = cms.untracked.bool( False ),
        Propagators = cms.untracked.vstring( 'hltESPFastSteppingHelixPropagatorAny' )
      ),
      TrackAssociatorParameters = cms.PSet( 
        useMuon = cms.bool( False ),
        truthMatch = cms.bool( False ),
        usePreshower = cms.bool( False ),
        dRPreshowerPreselection = cms.double( 0.2 ),
        muonMaxDistanceSigmaY = cms.double( 0.0 ),
        useEcal = cms.bool( False ),
        muonMaxDistanceSigmaX = cms.double( 0.0 ),
        dRMuon = cms.double( 9999.0 ),
        dREcal = cms.double( 0.5 ),
        CSCSegmentCollectionLabel = cms.InputTag( "hltCscSegments" ),
        DTRecSegment4DCollectionLabel = cms.InputTag( "hltDt4DSegments" ),
        EBRecHitCollectionLabel = cms.InputTag( "Notused" ),
        CaloTowerCollectionLabel = cms.InputTag( "Notused" ),
        propagateAllDirections = cms.bool( True ),
        muonMaxDistanceY = cms.double( 5.0 ),
        useHO = cms.bool( False ),
        muonMaxDistanceX = cms.double( 5.0 ),
        trajectoryUncertaintyTolerance = cms.double( -1.0 ),
        useHcal = cms.bool( False ),
        HBHERecHitCollectionLabel = cms.InputTag( "Notused" ),
        accountForTrajectoryChangeCalo = cms.bool( False ),
        dREcalPreselection = cms.double( 0.5 ),
        useCalo = cms.bool( True ),
        dRMuonPreselection = cms.double( 0.2 ),
        EERecHitCollectionLabel = cms.InputTag( "Notused" ),
        dRHcal = cms.double( 0.5 ),
        dRHcalPreselection = cms.double( 0.5 ),
        HORecHitCollectionLabel = cms.InputTag( "Notused" )
      ),
      Threshold = cms.double( 5.0 )
    ),
    trackDepositName = cms.string( "tracker" ),
    ecalDepositName = cms.string( "ecal" ),
    hcalDepositName = cms.string( "hcal" ),
    hoDepositName = cms.string( "ho" ),
    jetDepositName = cms.string( "jets" ),
    TimingFillerParameters = cms.PSet( 
      DTTimingParameters = cms.PSet( 
        HitError = cms.double( 6.0 ),
        MatchParameters = cms.PSet( 
          TightMatchDT = cms.bool( False ),
          DTradius = cms.double( 0.01 ),
          TightMatchCSC = cms.bool( True ),
          CSCsegments = cms.InputTag( "hltCscSegments" ),
          DTsegments = cms.InputTag( "hltDt4DSegments" )
        ),
        debug = cms.bool( False ),
        DoWireCorr = cms.bool( False ),
        RequireBothProjections = cms.bool( False ),
        DTTimeOffset = cms.double( 2.7 ),
        PruneCut = cms.double( 10000.0 ),
        DTsegments = cms.InputTag( "hltDt4DSegments" ),
        UseSegmentT0 = cms.bool( False ),
        HitsMin = cms.int32( 5 ),
        DropTheta = cms.bool( True ),
        ServiceParameters = cms.PSet( 
          RPCLayers = cms.bool( True ),
          Propagators = cms.untracked.vstring( 'hltESPFastSteppingHelixPropagatorAny' )
        )
      ),
      UseCSC = cms.bool( True ),
      CSCTimingParameters = cms.PSet( 
        MatchParameters = cms.PSet( 
          TightMatchDT = cms.bool( False ),
          DTradius = cms.double( 0.01 ),
          TightMatchCSC = cms.bool( True ),
          CSCsegments = cms.InputTag( "hltCscSegments" ),
          DTsegments = cms.InputTag( "hltDt4DSegments" )
        ),
        debug = cms.bool( False ),
        CSCWireTimeOffset = cms.double( 0.0 ),
        CSCStripError = cms.double( 7.0 ),
        CSCTimeOffset = cms.double( 0.0 ),
        CSCWireError = cms.double( 8.6 ),
        PruneCut = cms.double( 100.0 ),
        CSCsegments = cms.InputTag( "hltCscSegments" ),
        UseStripTime = cms.bool( True ),
        CSCStripTimeOffset = cms.double( 0.0 ),
        UseWireTime = cms.bool( True ),
        ServiceParameters = cms.PSet( 
          RPCLayers = cms.bool( True ),
          Propagators = cms.untracked.vstring( 'hltESPFastSteppingHelixPropagatorAny' )
        )
      ),
      ErrorDT = cms.double( 6.0 ),
      EcalEnergyCut = cms.double( 0.4 ),
      UseECAL = cms.bool( True ),
      ErrorEB = cms.double( 2.085 ),
      UseDT = cms.bool( True ),
      ErrorEE = cms.double( 6.95 ),
      ErrorCSC = cms.double( 7.4 )
    ),
    ShowerDigiFillerParameters = cms.PSet( 
      cscDigiCollectionLabel = cms.InputTag( 'muonCSCDigis','MuonCSCStripDigi' ),
      digiMaxDistanceX = cms.double( 25.0 ),
      dtDigiCollectionLabel = cms.InputTag( "muonDTDigis" )
    ),
    TrackerKinkFinderParameters = cms.PSet( 
      usePosition = cms.bool( False ),
      diagonalOnly = cms.bool( False )
    ),
    fillEnergy = cms.bool( False ),
    storeCrossedHcalRecHits = cms.bool( False ),
    maxAbsPullX = cms.double( 4.0 ),
    maxAbsEta = cms.double( 3.0 ),
    minPt = cms.double( 8.0 ),
    inputCollectionTypes = cms.vstring( 'inner tracks',
      'links' ),
    addExtraSoftMuons = cms.bool( False ),
    fillGlobalTrackRefits = cms.bool( False ),
    debugWithTruthMatching = cms.bool( False ),
    inputCollectionLabels = cms.VInputTag( 'hltIterL3MuonAndMuonFromL1MergedNoVtx','hltL3MuonsIterL3LinksNoVtx' ),
    fillCaloCompatibility = cms.bool( False ),
    maxAbsPullY = cms.double( 9999.0 ),
    maxAbsDy = cms.double( 9999.0 ),
    minP = cms.double( 0.0 ),
    minPCaloMuon = cms.double( 1.0E9 ),
    maxAbsDx = cms.double( 3.0 ),
    fillIsolation = cms.bool( False ),
    writeIsoDeposits = cms.bool( False ),
    minNumberOfMatches = cms.int32( 1 ),
    fillMatching = cms.bool( True ),
    fillShowerDigis = cms.bool( False ),
    ptThresholdToFillCandidateP4WithGlobalFit = cms.double( 200.0 ),
    sigmaThresholdToFillCandidateP4WithGlobalFit = cms.double( 2.0 ),
    fillGlobalTrackQuality = cms.bool( False ),
    globalTrackQualityInputTag = cms.InputTag( "glbTrackQual" ),
    selectHighPurity = cms.bool( False ),
    pvInputTag = cms.InputTag( "offlinePrimaryVertices" ),
    fillTrackerKink = cms.bool( False ),
    minCaloCompatibility = cms.double( 0.6 ),
    runArbitrationCleaner = cms.bool( False ),
    arbitrationCleanerOptions = cms.PSet( 
      OverlapDTheta = cms.double( 0.02 ),
      Overlap = cms.bool( True ),
      Clustering = cms.bool( True ),
      ME1a = cms.bool( True ),
      ClusterDTheta = cms.double( 0.02 ),
      ClusterDPhi = cms.double( 0.6 ),
      OverlapDPhi = cms.double( 0.0786 )
    ),
    arbitrateTrackerMuons = cms.bool( True )
)
process.hltIterL3MuonCandidatesNoVtx = cms.EDProducer( "L3MuonCandidateProducerFromMuons",
    InputObjects = cms.InputTag( "hltIterL3MuonsNoVtx" )
)
process.hltDisplacedmumuVtxNoMatchingProducer = cms.EDProducer( "HLTDisplacedmumuVtxProducer",
    Src = cms.InputTag( "hltIterL3MuonCandidatesNoVtx" ),
    PreviousCandTag = cms.InputTag( "hltDoubleMu3L3FilteredNoVtx" ),
    matchToPrevious = cms.bool( False ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 0.0 ),
    MinPtPair = cms.double( 0.0 ),
    MinInvMass = cms.double( 0.0 ),
    MaxInvMass = cms.double( 99999.0 ),
    ChargeOpt = cms.int32( 0 )
)
process.hltRecHitInRegionForMuonsMFnoVtx = cms.EDProducer( "MuonHLTRechitInRegionsProducer",
    productLabels = cms.vstring( 'EcalRegionalRecHitsEB',
      'EcalRegionalRecHitsEE' ),
    ecalhitLabels = cms.VInputTag( 'hltEcalRecHit:EcalRecHitsEB','hltEcalRecHit:EcalRecHitsEE' ),
    l1TagIsolated = cms.InputTag( "hltIterL3MuonCandidatesNoVtx" ),
    l1TagNonIsolated = cms.InputTag( "NotUsed" ),
    useUncalib = cms.bool( False ),
    doIsolated = cms.bool( True ),
    l1LowerThr = cms.double( 0.0 ),
    l1UpperThr = cms.double( 999.0 ),
    l1LowerThrIgnoreIsolation = cms.double( 100.0 ),
    regionEtaMargin = cms.double( 0.4 ),
    regionPhiMargin = cms.double( 0.4 )
)
process.hltRecHitInRegionForMuonsESNoVtx = cms.EDProducer( "MuonHLTRechitInRegionsProducer",
    productLabels = cms.vstring( 'EcalRegionalRecHitsES' ),
    ecalhitLabels = cms.VInputTag( 'hltEcalPreshowerRecHit:EcalRecHitsES' ),
    l1TagIsolated = cms.InputTag( "hltIterL3MuonCandidatesNoVtx" ),
    l1TagNonIsolated = cms.InputTag( "NotUsed" ),
    useUncalib = cms.bool( False ),
    doIsolated = cms.bool( True ),
    l1LowerThr = cms.double( 0.0 ),
    l1UpperThr = cms.double( 999.0 ),
    l1LowerThrIgnoreIsolation = cms.double( 100.0 ),
    regionEtaMargin = cms.double( 0.4 ),
    regionPhiMargin = cms.double( 0.4 )
)
process.hltParticleFlowRecHitECALForMuonsMFNoVtx = cms.EDProducer( "PFRecHitProducer",
    navigator = cms.PSet( 
      barrel = cms.PSet(  ),
      endcap = cms.PSet(  ),
      name = cms.string( "PFRecHitECALNavigator" )
    ),
    producers = cms.VPSet( 
      cms.PSet(  src = cms.InputTag( 'hltRecHitInRegionForMuonsMFnoVtx','EcalRegionalRecHitsEB' ),
        srFlags = cms.InputTag( "" ),
        name = cms.string( "PFEBRecHitCreator" ),
        qualityTests = cms.VPSet( 
          cms.PSet(  name = cms.string( "PFRecHitQTestDBThreshold" ),
            applySelectionsToAllCrystals = cms.bool( True )
          ),
          cms.PSet(  topologicalCleaning = cms.bool( True ),
            skipTTRecoveredHits = cms.bool( True ),
            cleaningThreshold = cms.double( 2.0 ),
            name = cms.string( "PFRecHitQTestECAL" ),
            timingCleaning = cms.bool( True )
          )
        )
      ),
      cms.PSet(  src = cms.InputTag( 'hltRecHitInRegionForMuonsMFnoVtx','EcalRegionalRecHitsEE' ),
        srFlags = cms.InputTag( "" ),
        name = cms.string( "PFEERecHitCreator" ),
        qualityTests = cms.VPSet( 
          cms.PSet(  name = cms.string( "PFRecHitQTestDBThreshold" ),
            applySelectionsToAllCrystals = cms.bool( True )
          ),
          cms.PSet(  topologicalCleaning = cms.bool( True ),
            skipTTRecoveredHits = cms.bool( True ),
            cleaningThreshold = cms.double( 2.0 ),
            name = cms.string( "PFRecHitQTestECAL" ),
            timingCleaning = cms.bool( True )
          )
        )
      )
    )
)
process.hltParticleFlowRecHitPSForMuonsNoVtx = cms.EDProducer( "PFRecHitProducer",
    navigator = cms.PSet(  name = cms.string( "PFRecHitPreshowerNavigator" ) ),
    producers = cms.VPSet( 
      cms.PSet(  src = cms.InputTag( 'hltRecHitInRegionForMuonsESNoVtx','EcalRegionalRecHitsES' ),
        name = cms.string( "PFPSRecHitCreator" ),
        qualityTests = cms.VPSet( 
          cms.PSet(  threshold = cms.double( 7.0E-6 ),
            name = cms.string( "PFRecHitQTestThreshold" )
          )
        )
      )
    )
)
process.hltParticleFlowClusterECALUncorrectedForMuonsMFNoVtx = cms.EDProducer( "PFClusterProducer",
    recHitsSource = cms.InputTag( "hltParticleFlowRecHitECALForMuonsMFNoVtx" ),
    recHitCleaners = cms.VPSet( 
    ),
    seedCleaners = cms.VPSet( 
    ),
    seedFinder = cms.PSet( 
      thresholdsByDetector = cms.VPSet( 
        cms.PSet(  seedingThresholdPt = cms.double( 0.15 ),
          seedingThreshold = cms.double( 0.6 ),
          detector = cms.string( "ECAL_ENDCAP" )
        ),
        cms.PSet(  seedingThresholdPt = cms.double( 0.0 ),
          seedingThreshold = cms.double( 0.23 ),
          detector = cms.string( "ECAL_BARREL" )
        )
      ),
      algoName = cms.string( "LocalMaximumSeedFinder" ),
      nNeighbours = cms.int32( 8 )
    ),
    initialClusteringStep = cms.PSet( 
      thresholdsByDetector = cms.VPSet( 
        cms.PSet(  gatheringThreshold = cms.double( 0.08 ),
          gatheringThresholdPt = cms.double( 0.0 ),
          detector = cms.string( "ECAL_BARREL" )
        ),
        cms.PSet(  gatheringThreshold = cms.double( 0.3 ),
          gatheringThresholdPt = cms.double( 0.0 ),
          detector = cms.string( "ECAL_ENDCAP" )
        )
      ),
      algoName = cms.string( "Basic2DGenericTopoClusterizer" ),
      useCornerCells = cms.bool( True )
    ),
    pfClusterBuilder = cms.PSet( 
      minFracTot = cms.double( 1.0E-20 ),
      stoppingTolerance = cms.double( 1.0E-8 ),
      positionCalc = cms.PSet( 
        minAllowedNormalization = cms.double( 1.0E-9 ),
        posCalcNCrystals = cms.int32( 9 ),
        algoName = cms.string( "Basic2DGenericPFlowPositionCalc" ),
        logWeightDenominator = cms.double( 0.08 ),
        minFractionInCalc = cms.double( 1.0E-9 ),
        timeResolutionCalcBarrel = cms.PSet( 
          corrTermLowE = cms.double( 0.0510871 ),
          threshLowE = cms.double( 0.5 ),
          noiseTerm = cms.double( 1.10889 ),
          constantTermLowE = cms.double( 0.0 ),
          noiseTermLowE = cms.double( 1.31883 ),
          threshHighE = cms.double( 5.0 ),
          constantTerm = cms.double( 0.428192 )
        ),
        timeResolutionCalcEndcap = cms.PSet( 
          corrTermLowE = cms.double( 0.0 ),
          threshLowE = cms.double( 1.0 ),
          noiseTerm = cms.double( 5.72489999999 ),
          constantTermLowE = cms.double( 0.0 ),
          noiseTermLowE = cms.double( 6.92683000001 ),
          threshHighE = cms.double( 10.0 ),
          constantTerm = cms.double( 0.0 )
        )
      ),
      maxIterations = cms.uint32( 50 ),
      positionCalcForConvergence = cms.PSet( 
        minAllowedNormalization = cms.double( 0.0 ),
        T0_ES = cms.double( 1.2 ),
        algoName = cms.string( "ECAL2DPositionCalcWithDepthCorr" ),
        T0_EE = cms.double( 3.1 ),
        T0_EB = cms.double( 7.4 ),
        X0 = cms.double( 0.89 ),
        minFractionInCalc = cms.double( 0.0 ),
        W0 = cms.double( 4.2 )
      ),
      allCellsPositionCalc = cms.PSet( 
        minAllowedNormalization = cms.double( 1.0E-9 ),
        posCalcNCrystals = cms.int32( -1 ),
        algoName = cms.string( "Basic2DGenericPFlowPositionCalc" ),
        logWeightDenominator = cms.double( 0.08 ),
        minFractionInCalc = cms.double( 1.0E-9 ),
        timeResolutionCalcBarrel = cms.PSet( 
          corrTermLowE = cms.double( 0.0510871 ),
          threshLowE = cms.double( 0.5 ),
          noiseTerm = cms.double( 1.10889 ),
          constantTermLowE = cms.double( 0.0 ),
          noiseTermLowE = cms.double( 1.31883 ),
          threshHighE = cms.double( 5.0 ),
          constantTerm = cms.double( 0.428192 )
        ),
        timeResolutionCalcEndcap = cms.PSet( 
          corrTermLowE = cms.double( 0.0 ),
          threshLowE = cms.double( 1.0 ),
          noiseTerm = cms.double( 5.72489999999 ),
          constantTermLowE = cms.double( 0.0 ),
          noiseTermLowE = cms.double( 6.92683000001 ),
          threshHighE = cms.double( 10.0 ),
          constantTerm = cms.double( 0.0 )
        )
      ),
      algoName = cms.string( "Basic2DGenericPFlowClusterizer" ),
      recHitEnergyNorms = cms.VPSet( 
        cms.PSet(  recHitEnergyNorm = cms.double( 0.08 ),
          detector = cms.string( "ECAL_BARREL" )
        ),
        cms.PSet(  recHitEnergyNorm = cms.double( 0.3 ),
          detector = cms.string( "ECAL_ENDCAP" )
        )
      ),
      showerSigma = cms.double( 1.5 ),
      minFractionToKeep = cms.double( 1.0E-7 ),
      excludeOtherSeeds = cms.bool( True )
    ),
    positionReCalc = cms.PSet( 
      minAllowedNormalization = cms.double( 0.0 ),
      T0_ES = cms.double( 1.2 ),
      algoName = cms.string( "ECAL2DPositionCalcWithDepthCorr" ),
      T0_EE = cms.double( 3.1 ),
      T0_EB = cms.double( 7.4 ),
      X0 = cms.double( 0.89 ),
      minFractionInCalc = cms.double( 0.0 ),
      W0 = cms.double( 4.2 )
    ),
    energyCorrector = cms.PSet(  )
)
process.hltParticleFlowClusterPSForMuonsNoVtx = cms.EDProducer( "PFClusterProducer",
    recHitsSource = cms.InputTag( "hltParticleFlowRecHitPSForMuonsNoVtx" ),
    recHitCleaners = cms.VPSet( 
    ),
    seedCleaners = cms.VPSet( 
    ),
    seedFinder = cms.PSet( 
      thresholdsByDetector = cms.VPSet( 
        cms.PSet(  seedingThresholdPt = cms.double( 0.0 ),
          seedingThreshold = cms.double( 1.2E-4 ),
          detector = cms.string( "PS1" )
        ),
        cms.PSet(  seedingThresholdPt = cms.double( 0.0 ),
          seedingThreshold = cms.double( 1.2E-4 ),
          detector = cms.string( "PS2" )
        )
      ),
      algoName = cms.string( "LocalMaximumSeedFinder" ),
      nNeighbours = cms.int32( 4 )
    ),
    initialClusteringStep = cms.PSet( 
      thresholdsByDetector = cms.VPSet( 
        cms.PSet(  gatheringThreshold = cms.double( 6.0E-5 ),
          gatheringThresholdPt = cms.double( 0.0 ),
          detector = cms.string( "PS1" )
        ),
        cms.PSet(  gatheringThreshold = cms.double( 6.0E-5 ),
          gatheringThresholdPt = cms.double( 0.0 ),
          detector = cms.string( "PS2" )
        )
      ),
      algoName = cms.string( "Basic2DGenericTopoClusterizer" ),
      useCornerCells = cms.bool( False )
    ),
    pfClusterBuilder = cms.PSet( 
      minFracTot = cms.double( 1.0E-20 ),
      stoppingTolerance = cms.double( 1.0E-8 ),
      positionCalc = cms.PSet( 
        minAllowedNormalization = cms.double( 1.0E-9 ),
        posCalcNCrystals = cms.int32( -1 ),
        algoName = cms.string( "Basic2DGenericPFlowPositionCalc" ),
        logWeightDenominator = cms.double( 6.0E-5 ),
        minFractionInCalc = cms.double( 1.0E-9 )
      ),
      maxIterations = cms.uint32( 50 ),
      algoName = cms.string( "Basic2DGenericPFlowClusterizer" ),
      recHitEnergyNorms = cms.VPSet( 
        cms.PSet(  recHitEnergyNorm = cms.double( 6.0E-5 ),
          detector = cms.string( "PS1" )
        ),
        cms.PSet(  recHitEnergyNorm = cms.double( 6.0E-5 ),
          detector = cms.string( "PS2" )
        )
      ),
      showerSigma = cms.double( 0.3 ),
      minFractionToKeep = cms.double( 1.0E-7 ),
      excludeOtherSeeds = cms.bool( True )
    ),
    positionReCalc = cms.PSet(  ),
    energyCorrector = cms.PSet(  )
)
process.hltParticleFlowClusterECALForMuonsMFNoVtx = cms.EDProducer( "CorrectedECALPFClusterProducer",
    minimumPSEnergy = cms.double( 0.0 ),
    skipPS = cms.bool( False ),
    inputPS = cms.InputTag( "hltParticleFlowClusterPSForMuonsNoVtx" ),
    energyCorrector = cms.PSet( 
      applyCrackCorrections = cms.bool( False ),
      srfAwareCorrection = cms.bool( True ),
      applyMVACorrections = cms.bool( True ),
      maxPtForMVAEvaluation = cms.double( 300.0 ),
      recHitsEBLabel = cms.InputTag( 'hltEcalRecHit','EcalRecHitsEB' ),
      recHitsEELabel = cms.InputTag( 'hltEcalRecHit','EcalRecHitsEE' ),
      ebSrFlagLabel = cms.InputTag( "hltEcalDigis" ),
      eeSrFlagLabel = cms.InputTag( "hltEcalDigis" )
    ),
    inputECAL = cms.InputTag( "hltParticleFlowClusterECALUncorrectedForMuonsMFNoVtx" )
)
process.hltMuonEcalMFPFClusterIsoForMuonsNoVtx = cms.EDProducer( "MuonHLTEcalPFClusterIsolationProducer",
    recoCandidateProducer = cms.InputTag( "hltIterL3MuonCandidatesNoVtx" ),
    pfClusterProducer = cms.InputTag( "hltParticleFlowClusterECALForMuonsMFNoVtx" ),
    rhoProducer = cms.InputTag( "hltFixedGridRhoFastjetECALMFForMuons" ),
    doRhoCorrection = cms.bool( True ),
    rhoMax = cms.double( 9.9999999E7 ),
    rhoScale = cms.double( 1.0 ),
    drMax = cms.double( 0.3 ),
    drVetoBarrel = cms.double( 0.05 ),
    drVetoEndcap = cms.double( 0.05 ),
    etaStripBarrel = cms.double( 0.0 ),
    etaStripEndcap = cms.double( 0.0 ),
    energyBarrel = cms.double( 0.0 ),
    energyEndcap = cms.double( 0.0 ),
    effectiveAreas = cms.vdouble( 0.35, 0.193 ),
    absEtaLowEdges = cms.vdouble( 0.0, 1.479 )
)
process.hltMuonHcalPFClusterIsoForMuonsNoVtx = cms.EDProducer( "MuonHLTHcalPFClusterIsolationProducer",
    recoCandidateProducer = cms.InputTag( "hltIterL3MuonCandidatesNoVtx" ),
    pfClusterProducerHCAL = cms.InputTag( "hltParticleFlowClusterHCAL" ),
    useHF = cms.bool( False ),
    pfClusterProducerHFEM = cms.InputTag( "" ),
    pfClusterProducerHFHAD = cms.InputTag( "" ),
    rhoProducer = cms.InputTag( "hltFixedGridRhoFastjetAllCalo" ),
    doRhoCorrection = cms.bool( True ),
    rhoMax = cms.double( 9.9999999E7 ),
    rhoScale = cms.double( 1.0 ),
    drMax = cms.double( 0.3 ),
    drVetoBarrel = cms.double( 0.1 ),
    drVetoEndcap = cms.double( 0.1 ),
    etaStripBarrel = cms.double( 0.0 ),
    etaStripEndcap = cms.double( 0.0 ),
    energyBarrel = cms.double( 0.0 ),
    energyEndcap = cms.double( 0.0 ),
    useEt = cms.bool( True ),
    effectiveAreas = cms.vdouble( 0.11, 0.163 ),
    absEtaLowEdges = cms.vdouble( 0.0, 1.479 )
)
process.hltPixelTracksTrackingRegionsForSeedsL3MuonNoVtx = cms.EDProducer( "CandidateSeededTrackingRegionsEDProducer",
    RegionPSet = cms.PSet( 
      vertexCollection = cms.InputTag( "hltPixelVertices" ),
      zErrorVetex = cms.double( 0.2 ),
      beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
      zErrorBeamSpot = cms.double( 24.2 ),
      maxNVertices = cms.int32( 4 ),
      maxNRegions = cms.int32( 40 ),
      nSigmaZVertex = cms.double( 3.0 ),
      nSigmaZBeamSpot = cms.double( 4.0 ),
      ptMin = cms.double( 0.3 ),
      mode = cms.string( "VerticesFixed" ),
      input = cms.InputTag( "hltIterL3MuonCandidatesNoVtx" ),
      searchOpt = cms.bool( False ),
      whereToUseMeasurementTracker = cms.string( "Never" ),
      originRadius = cms.double( 0.1 ),
      measurementTrackerName = cms.InputTag( "" ),
      precise = cms.bool( True ),
      deltaEta = cms.double( 0.3 ),
      deltaPhi = cms.double( 0.3 )
    )
)
process.hltPixelTracksInRegionIter0L3MuonNoVtx = cms.EDProducer( "TrackSelectorByRegion",
    tracks = cms.InputTag( "hltPixelTracks" ),
    regions = cms.InputTag( "hltPixelTracksTrackingRegionsForSeedsL3MuonNoVtx" ),
    produceTrackCollection = cms.bool( True ),
    produceMask = cms.bool( False )
)
process.hltIter0L3MuonPixelSeedsFromPixelTracksNoVtx = cms.EDProducer( "SeedGeneratorFromProtoTracksEDProducer",
    InputCollection = cms.InputTag( "hltPixelTracksInRegionIter0L3MuonNoVtx" ),
    InputVertexCollection = cms.InputTag( "hltTrimmedPixelVertices" ),
    originHalfLength = cms.double( 0.3 ),
    originRadius = cms.double( 0.1 ),
    useProtoTrackKinematics = cms.bool( False ),
    useEventsWithNoVertex = cms.bool( True ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    usePV = cms.bool( False ),
    includeFourthHit = cms.bool( True ),
    SeedCreatorPSet = cms.PSet(  refToPSet_ = cms.string( "HLTSeedFromProtoTracks" ) )
)
process.hltIter0L3MuonCkfTrackCandidatesNoVtx = cms.EDProducer( "CkfTrackCandidateMaker",
    cleanTrajectoryAfterInOut = cms.bool( False ),
    doSeedingRegionRebuilding = cms.bool( False ),
    onlyPixelHitsForSeedCleaner = cms.bool( False ),
    reverseTrajectories = cms.bool( False ),
    useHitsSplitting = cms.bool( False ),
    MeasurementTrackerEvent = cms.InputTag( "hltSiStripClusters" ),
    src = cms.InputTag( "hltIter0L3MuonPixelSeedsFromPixelTracksNoVtx" ),
    clustersToSkip = cms.InputTag( "" ),
    phase2clustersToSkip = cms.InputTag( "" ),
    TrajectoryBuilderPSet = cms.PSet(  refToPSet_ = cms.string( "HLTIter0GroupedCkfTrajectoryBuilderIT" ) ),
    TransientInitialStateEstimatorParameters = cms.PSet( 
      propagatorAlongTISE = cms.string( "PropagatorWithMaterialParabolicMf" ),
      numberMeasurementsForFit = cms.int32( 4 ),
      propagatorOppositeTISE = cms.string( "PropagatorWithMaterialParabolicMfOpposite" )
    ),
    numHitsForSeedCleaner = cms.int32( 4 ),
    NavigationSchool = cms.string( "SimpleNavigationSchool" ),
    RedundantSeedCleaner = cms.string( "CachingSeedCleanerBySharedInput" ),
    TrajectoryCleaner = cms.string( "hltESPTrajectoryCleanerBySharedHits" ),
    maxNSeeds = cms.uint32( 100000 ),
    maxSeedsBeforeCleaning = cms.uint32( 1000 )
)
process.hltIter0L3MuonCtfWithMaterialTracksNoVtx = cms.EDProducer( "TrackProducer",
    useSimpleMF = cms.bool( True ),
    SimpleMagneticField = cms.string( "ParabolicMf" ),
    src = cms.InputTag( "hltIter0L3MuonCkfTrackCandidatesNoVtx" ),
    clusterRemovalInfo = cms.InputTag( "" ),
    beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
    Fitter = cms.string( "hltESPFittingSmootherIT" ),
    useHitsSplitting = cms.bool( False ),
    alias = cms.untracked.string( "ctfWithMaterialTracks" ),
    TrajectoryInEvent = cms.bool( False ),
    TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
    AlgorithmName = cms.string( "hltIter0" ),
    Propagator = cms.string( "hltESPRungeKuttaTrackerPropagator" ),
    GeometricInnerState = cms.bool( True ),
    NavigationSchool = cms.string( "" ),
    MeasurementTracker = cms.string( "" ),
    MeasurementTrackerEvent = cms.InputTag( "hltSiStripClusters" )
)
process.hltIter0L3MuonTrackCutClassifierNoVtx = cms.EDProducer( "TrackCutClassifier",
    src = cms.InputTag( "hltIter0L3MuonCtfWithMaterialTracksNoVtx" ),
    beamspot = cms.InputTag( "hltOnlineBeamSpot" ),
    vertices = cms.InputTag( "hltTrimmedPixelVertices" ),
    ignoreVertices = cms.bool( False ),
    qualityCuts = cms.vdouble( -0.7, 0.1, 0.7 ),
    mva = cms.PSet( 
      minPixelHits = cms.vint32( 0, 0, 0 ),
      maxDzWrtBS = cms.vdouble( 3.40282346639E38, 24.0, 15.0 ),
      dr_par = cms.PSet( 
        d0err = cms.vdouble( 0.003, 0.003, 0.003 ),
        dr_par2 = cms.vdouble( 3.40282346639E38, 0.6, 0.6 ),
        dr_par1 = cms.vdouble( 3.40282346639E38, 0.8, 0.8 ),
        dr_exp = cms.vint32( 4, 4, 4 ),
        d0err_par = cms.vdouble( 0.001, 0.001, 0.001 )
      ),
      maxLostLayers = cms.vint32( 1, 1, 1 ),
      min3DLayers = cms.vint32( 0, 0, 0 ),
      dz_par = cms.PSet( 
        dz_par1 = cms.vdouble( 3.40282346639E38, 0.75, 0.75 ),
        dz_par2 = cms.vdouble( 3.40282346639E38, 0.5, 0.5 ),
        dz_exp = cms.vint32( 4, 4, 4 )
      ),
      minNVtxTrk = cms.int32( 3 ),
      maxDz = cms.vdouble( 0.5, 0.2, 3.40282346639E38 ),
      minNdof = cms.vdouble( 1.0E-5, 1.0E-5, 1.0E-5 ),
      maxChi2 = cms.vdouble( 9999.0, 25.0, 16.0 ),
      maxChi2n = cms.vdouble( 1.2, 1.0, 0.7 ),
      maxDr = cms.vdouble( 0.5, 0.03, 3.40282346639E38 ),
      minLayers = cms.vint32( 3, 3, 3 )
    )
)
process.hltIter0L3MuonTrackSelectionHighPurityNoVtx = cms.EDProducer( "TrackCollectionFilterCloner",
    originalSource = cms.InputTag( "hltIter0L3MuonCtfWithMaterialTracksNoVtx" ),
    originalMVAVals = cms.InputTag( 'hltIter0L3MuonTrackCutClassifierNoVtx','MVAValues' ),
    originalQualVals = cms.InputTag( 'hltIter0L3MuonTrackCutClassifierNoVtx','QualityMasks' ),
    minQuality = cms.string( "highPurity" ),
    copyExtras = cms.untracked.bool( True ),
    copyTrajectories = cms.untracked.bool( False )
)
process.hltMuonTkRelIsolationCut0p09MapNoVtx = cms.EDProducer( "L3MuonCombinedRelativeIsolationProducer",
    UseRhoCorrectedCaloDeposits = cms.bool( False ),
    UseCaloIso = cms.bool( False ),
    CaloDepositsLabel = cms.InputTag( "notUsed" ),
    inputMuonCollection = cms.InputTag( "hltIterL3MuonCandidatesNoVtx" ),
    OutputMuIsoDeposits = cms.bool( True ),
    TrackPt_Min = cms.double( -1.0 ),
    printDebug = cms.bool( False ),
    CutsPSet = cms.PSet( 
      applyCutsORmaxNTracks = cms.bool( False ),
      maxNTracks = cms.int32( -1 ),
      Thresholds = cms.vdouble( 0.09 ),
      EtaBounds = cms.vdouble( 2.411 ),
      ComponentName = cms.string( "SimpleCuts" ),
      ConeSizes = cms.vdouble( 0.3 )
    ),
    TrkExtractorPSet = cms.PSet( 
      Diff_z = cms.double( 0.2 ),
      inputTrackCollection = cms.InputTag( "hltIter0L3MuonTrackSelectionHighPurityNoVtx" ),
      Chi2Ndof_Max = cms.double( 1.0E64 ),
      BeamSpotLabel = cms.InputTag( "hltOnlineBeamSpot" ),
      DR_Veto = cms.double( 0.01 ),
      Pt_Min = cms.double( -1.0 ),
      VetoLeadingTrack = cms.bool( True ),
      DR_Max = cms.double( 0.3 ),
      DepositLabel = cms.untracked.string( "PXLS" ),
      PtVeto_Min = cms.double( 2.0 ),
      NHits_Min = cms.uint32( 0 ),
      PropagateTracksToRadius = cms.bool( True ),
      ReferenceRadius = cms.double( 6.0 ),
      Chi2Prob_Min = cms.double( -1.0 ),
      Diff_r = cms.double( 0.1 ),
      BeamlineOption = cms.string( "BeamSpotFromEvent" ),
      ComponentName = cms.string( "PixelTrackExtractor" ),
      DR_VetoPt = cms.double( 0.025 )
    ),
    CaloExtractorPSet = cms.PSet( 
      DR_Veto_H = cms.double( 0.1 ),
      Vertex_Constraint_Z = cms.bool( False ),
      DR_Veto_E = cms.double( 0.07 ),
      Weight_H = cms.double( 1.0 ),
      CaloTowerCollectionLabel = cms.InputTag( "hltTowerMakerForAll" ),
      DR_Max = cms.double( 0.3 ),
      DepositLabel = cms.untracked.string( "EcalPlusHcal" ),
      Vertex_Constraint_XY = cms.bool( False ),
      Threshold_H = cms.double( 0.5 ),
      Threshold_E = cms.double( 0.2 ),
      ComponentName = cms.string( "CaloExtractor" ),
      Weight_E = cms.double( 1.0 )
    )
)
process.hltFEDSelectorL1 = cms.EDProducer( "EvFFEDSelector",
    inputTag = cms.InputTag( "rawDataCollector" ),
    fedList = cms.vuint32( 1404 )
)

### Full tracking track packer
process.hltScoutingFullTrackPacker = cms.EDProducer( "HLTScoutingTrackProducer",
    OtherTracks = cms.InputTag( "hltPFMuonMerging" ),
    vertexCollection = cms.InputTag( "hltPixelVertices" ),
    mantissaPrecision = cms.int32( 10 ),
    vtxMinDist = cms.double( 0.01 ),
    ptMin = cms.double( 0.3 )
)
process.hltScoutingPrimaryVertexPacker = cms.EDProducer( "HLTScoutingPrimaryVertexProducer",
    vertexCollection = cms.InputTag( 'hltPixelVertices','','@currentProcess' ),
    mantissaPrecision = cms.int32( 10 )
)

### Full tracking PF packer
process.hltScoutingFullTrackPFPacker = cms.EDProducer( "HLTScoutingPFProducer",
    pfJetCollection = cms.InputTag( 'hltAK4PFJets','','@currentProcess' ),
    pfJetTagCollection = cms.InputTag( 'hltCombinedSecondaryVertexBJetTagsPF','','@currentProcess' ),
    pfCandidateCollection = cms.InputTag( 'hltParticleFlow','','@currentProcess' ),
    vertexCollection = cms.InputTag( 'hltPixelVertices','','@currentProcess' ),
    metCollection = cms.InputTag( 'hltPFMETProducer','','@currentProcess' ),
    rho = cms.InputTag( 'hltFixedGridRhoFastjetAll','','@currentProcess' ),
    pfJetPtCut = cms.double( 20.0 ),
    pfJetEtaCut = cms.double( 3.0 ),
    pfCandidatePtCut = cms.double( 0.6 ),
    pfCandidateEtaCut = cms.double( 3.0 ),
    mantissaPrecision = cms.int32( 10 ),
    doJetTags = cms.bool( True ),
    doCandidates = cms.bool( True ),
    doMet = cms.bool( True ),
    doTrackRelVars = cms.bool( True ),
    doCandIndsForJets = cms.bool( False )
)
process.hltScoutingMuonPacker = cms.EDProducer( "HLTScoutingMuonProducer",
    ChargedCandidates = cms.InputTag( 'hltIterL3MuonCandidatesNoVtx','','@currentProcess' ),
    displacedvertexCollection = cms.InputTag( 'hltDisplacedmumuVtxNoMatchingProducer','','@currentProcess' ),
    InputMuons = cms.InputTag( "hltIterL3MuonsNoVtx" ),
    InputLinks = cms.InputTag( 'hltL3MuonsIterL3LinksNoVtx','','@currentProcess' ),
    Tracks = cms.InputTag( "hltIterL3MuonAndMuonFromL1MergedNoVtx" ),
    EcalPFClusterIsoMap = cms.InputTag( 'hltMuonEcalMFPFClusterIsoForMuonsNoVtx','','@currentProcess' ),
    HcalPFClusterIsoMap = cms.InputTag( 'hltMuonHcalPFClusterIsoForMuonsNoVtx','','@currentProcess' ),
    TrackIsoMap = cms.InputTag( 'hltMuonTkRelIsolationCut0p09MapNoVtx','combinedRelativeIsoDeposits','@currentProcess' ),
    muonPtCut = cms.double( 0.0 ),
    muonEtaCut = cms.double( 2.4 ),
    minVtxProbCut = cms.double( 0.001 )
)
process.hltScoutingEgammaPacker = cms.EDProducer( "HLTScoutingEgammaProducer",
    EgammaCandidates = cms.InputTag( 'hltEgammaCandidates','','@currentProcess' ),
    EgammaGsfTracks = cms.InputTag( 'hltEgammaGsfTracks','','@currentProcess' ),
    SigmaIEtaIEtaMap = cms.InputTag( 'hltEgammaClusterShape','sigmaIEtaIEta5x5','@currentProcess' ),
    r9Map = cms.InputTag( 'hltEgammaR9ID','r95x5' ),
    HoverEMap = cms.InputTag( 'hltEgammaHoverE','','@currentProcess' ),
    DetaMap = cms.InputTag( 'hltEgammaGsfTrackVars','DetaSeed','@currentProcess' ),
    DphiMap = cms.InputTag( 'hltEgammaGsfTrackVars','Dphi','@currentProcess' ),
    MissingHitsMap = cms.InputTag( 'hltEgammaGsfTrackVars','MissingHits','@currentProcess' ),
    OneOEMinusOneOPMap = cms.InputTag( 'hltEgammaGsfTrackVars','OneOESuperMinusOneOP','@currentProcess' ),
    EcalPFClusterIsoMap = cms.InputTag( 'hltEgammaEcalPFClusterIso','','@currentProcess' ),
    EleGsfTrackIsoMap = cms.InputTag( 'hltEgammaEleGsfTrackIsoPixelOnly','','@currentProcess' ),
    HcalPFClusterIsoMap = cms.InputTag( 'hltEgammaHcalPFClusterIso','','@currentProcess' ),
    egammaPtCut = cms.double( 2.0 ),
    egammaEtaCut = cms.double( 2.5 ),
    egammaHoverECut = cms.double( 1.0 ),
    saveRecHitTiming = cms.bool( False ),
    mantissaPrecision = cms.int32( 10 ),
    rechitMatrixSize = cms.int32( 15 ),
    ecalRechitEB = cms.InputTag( 'hltEcalRecHit','EcalRecHitsEB' ),
    ecalRechitEE = cms.InputTag( 'hltEcalRecHit','EcalRecHitsEE' )
)
process.hltBoolEnd = cms.EDFilter( "HLTBool",
    result = cms.bool( True )
)
process.hltPreDSTRun3PFScoutingPixelTracking = cms.EDFilter( "HLTPrescaler",
    offset = cms.uint32( 0 ),
    L1GtReadoutRecordTag = cms.InputTag( "hltGtStage2Digis" )
)
process.hltScoutingTrackPacker = cms.EDProducer( "HLTScoutingTrackProducer",
    OtherTracks = cms.InputTag( "hltPixelOnlyPFMuonMerging" ),
    vertexCollection = cms.InputTag( "hltPixelVertices" ),
    mantissaPrecision = cms.int32( 10 ),
    vtxMinDist = cms.double( 0.01 ),
    ptMin = cms.double( 0.3 )
)
process.hltScoutingPFPacker = cms.EDProducer( "HLTScoutingPFProducer",
    pfJetCollection = cms.InputTag( 'hltAK4PixelOnlyPFJets','','@currentProcess' ),
    pfJetTagCollection = cms.InputTag( 'hltCombinedSecondaryVertexBJetTagsPF','','@currentProcess' ),
    pfCandidateCollection = cms.InputTag( 'hltPixelOnlyParticleFlow','','@currentProcess' ),
    vertexCollection = cms.InputTag( 'hltPixelVertices','','@currentProcess' ),
    metCollection = cms.InputTag( 'hltPixelOnlyPFMETProducer','','@currentProcess' ),
    rho = cms.InputTag( 'hltFixedGridRhoFastjetPixelOnlyAll','','@currentProcess' ),
    pfJetPtCut = cms.double( 20.0 ),
    pfJetEtaCut = cms.double( 3.0 ),
    pfCandidatePtCut = cms.double( 0.6 ),
    pfCandidateEtaCut = cms.double( 3.0 ),
    mantissaPrecision = cms.int32( 10 ),
    doJetTags = cms.bool( True ),
    doCandidates = cms.bool( True ),
    doMet = cms.bool( True ),
    doTrackRelVars = cms.bool( True ),
    doCandIndsForJets = cms.bool( False )
)
process.hltFEDSelector = cms.EDProducer( "EvFFEDSelector",
    inputTag = cms.InputTag( "rawDataCollector" ),
    fedList = cms.vuint32( 1023, 1024 )
)
process.hltTriggerSummaryAOD = cms.EDProducer( "TriggerSummaryProducerAOD",
    throw = cms.bool( False ),
    processName = cms.string( "@" ),
    moduleLabelPatternsToMatch = cms.vstring( 'hlt*' ),
    moduleLabelPatternsToSkip = cms.vstring(  )
)
process.hltTriggerSummaryRAW = cms.EDProducer( "TriggerSummaryProducerRAW",
    processName = cms.string( "@" )
)

process.hltEcalDigis = SwitchProducerCUDA(
   cpu = cms.EDAlias(
       hltEcalDigisLegacy = cms.VPSet( 
         cms.PSet(  type = cms.string( "EBDigiCollection" )         ),
         cms.PSet(  type = cms.string( "EEDigiCollection" )         ),
         cms.PSet(  type = cms.string( "EBDetIdedmEDCollection" )         ),
         cms.PSet(  type = cms.string( "EEDetIdedmEDCollection" )         ),
         cms.PSet(  type = cms.string( "EBSrFlagsSorted" )         ),
         cms.PSet(  type = cms.string( "EESrFlagsSorted" )         ),
         cms.PSet(  type = cms.string( "EcalElectronicsIdedmEDCollection" ),
           fromProductInstance = cms.string( "EcalIntegrityBlockSizeErrors" )
         ),
         cms.PSet(  type = cms.string( "EcalElectronicsIdedmEDCollection" ),
           fromProductInstance = cms.string( "EcalIntegrityTTIdErrors" )
         ),
         cms.PSet(  type = cms.string( "EcalElectronicsIdedmEDCollection" ),
           fromProductInstance = cms.string( "EcalIntegrityZSXtalIdErrors" )
         ),
         cms.PSet(  type = cms.string( "EcalPnDiodeDigisSorted" )         ),
         cms.PSet(  type = cms.string( "EcalPseudoStripInputDigisSorted" ),
           fromProductInstance = cms.string( "EcalPseudoStripInputs" )
         ),
         cms.PSet(  type = cms.string( "EcalTriggerPrimitiveDigisSorted" ),
           fromProductInstance = cms.string( "EcalTriggerPrimitives" )
         )
       )
   ),
 )
process.hltEcalUncalibRecHit = SwitchProducerCUDA(
   cpu = cms.EDProducer( "EcalUncalibRecHitProducer",
       EBdigiCollection = cms.InputTag( 'hltEcalDigis','ebDigis' ),
       EEhitCollection = cms.string( "EcalUncalibRecHitsEE" ),
       EEdigiCollection = cms.InputTag( 'hltEcalDigis','eeDigis' ),
       EBhitCollection = cms.string( "EcalUncalibRecHitsEB" ),
       algo = cms.string( "EcalUncalibRecHitWorkerMultiFit" ),
       algoPSet = cms.PSet( 
         EBtimeFitLimits_Upper = cms.double( 1.4 ),
         EEtimeFitLimits_Lower = cms.double( 0.2 ),
         timealgo = cms.string( "RatioMethod" ),
         EBtimeNconst = cms.double( 28.5 ),
         prefitMaxChiSqEE = cms.double( 10.0 ),
         outOfTimeThresholdGain12mEB = cms.double( 1000.0 ),
         outOfTimeThresholdGain12mEE = cms.double( 1000.0 ),
         EEtimeFitParameters = cms.vdouble( -2.390548, 3.553628, -17.62341, 67.67538, -133.213, 140.7432, -75.41106, 16.20277 ),
         prefitMaxChiSqEB = cms.double( 25.0 ),
         simplifiedNoiseModelForGainSwitch = cms.bool( True ),
         EBtimeFitParameters = cms.vdouble( -2.015452, 3.130702, -12.3473, 41.88921, -82.83944, 91.01147, -50.35761, 11.05621 ),
         selectiveBadSampleCriteriaEB = cms.bool( False ),
         dynamicPedestalsEB = cms.bool( False ),
         useLumiInfoRunHeader = cms.bool( False ),
         EBamplitudeFitParameters = cms.vdouble( 1.138, 1.652 ),
         doPrefitEE = cms.bool( False ),
         dynamicPedestalsEE = cms.bool( False ),
         selectiveBadSampleCriteriaEE = cms.bool( False ),
         outOfTimeThresholdGain61pEE = cms.double( 1000.0 ),
         outOfTimeThresholdGain61pEB = cms.double( 1000.0 ),
         activeBXs = cms.vint32( -5, -4, -3, -2, -1, 0, 1, 2, 3, 4 ),
         doPrefitEB = cms.bool( False ),
         addPedestalUncertaintyEE = cms.double( 0.0 ),
         addPedestalUncertaintyEB = cms.double( 0.0 ),
         gainSwitchUseMaxSampleEB = cms.bool( True ),
         EEtimeNconst = cms.double( 31.8 ),
         EEamplitudeFitParameters = cms.vdouble( 1.89, 1.4 ),
         outOfTimeThresholdGain12pEB = cms.double( 1000.0 ),
         gainSwitchUseMaxSampleEE = cms.bool( False ),
         mitigateBadSamplesEB = cms.bool( False ),
         outOfTimeThresholdGain12pEE = cms.double( 1000.0 ),
         ampErrorCalculation = cms.bool( False ),
         mitigateBadSamplesEE = cms.bool( False ),
         amplitudeThresholdEB = cms.double( 10.0 ),
         amplitudeThresholdEE = cms.double( 10.0 ),
         EBtimeFitLimits_Lower = cms.double( 0.2 ),
         EEtimeFitLimits_Upper = cms.double( 1.4 ),
         outOfTimeThresholdGain61mEE = cms.double( 1000.0 ),
         EEtimeConstantTerm = cms.double( 1.0 ),
         EBtimeConstantTerm = cms.double( 0.6 ),
         outOfTimeThresholdGain61mEB = cms.double( 1000.0 )
       )
   ),
 )
process.hltEcalRecHit = SwitchProducerCUDA(
   cpu = cms.EDProducer( "EcalRecHitProducer",
       recoverEEVFE = cms.bool( False ),
       EErechitCollection = cms.string( "EcalRecHitsEE" ),
       recoverEBIsolatedChannels = cms.bool( False ),
       recoverEBVFE = cms.bool( False ),
       laserCorrection = cms.bool( True ),
       EBLaserMIN = cms.double( 0.5 ),
       killDeadChannels = cms.bool( True ),
       dbStatusToBeExcludedEB = cms.vint32( 14, 78, 142 ),
       EEuncalibRecHitCollection = cms.InputTag( 'hltEcalUncalibRecHit','EcalUncalibRecHitsEE' ),
       dbStatusToBeExcludedEE = cms.vint32( 14, 78, 142 ),
       EELaserMIN = cms.double( 0.5 ),
       ebFEToBeRecovered = cms.InputTag( 'hltEcalDetIdToBeRecovered','ebFE' ),
       cleaningConfig = cms.PSet( 
         e6e2thresh = cms.double( 0.04 ),
         tightenCrack_e6e2_double = cms.double( 3.0 ),
         e4e1Threshold_endcap = cms.double( 0.3 ),
         tightenCrack_e4e1_single = cms.double( 3.0 ),
         tightenCrack_e1_double = cms.double( 2.0 ),
         cThreshold_barrel = cms.double( 4.0 ),
         e4e1Threshold_barrel = cms.double( 0.08 ),
         tightenCrack_e1_single = cms.double( 2.0 ),
         e4e1_b_barrel = cms.double( -0.024 ),
         e4e1_a_barrel = cms.double( 0.04 ),
         ignoreOutOfTimeThresh = cms.double( 1.0E9 ),
         cThreshold_endcap = cms.double( 15.0 ),
         e4e1_b_endcap = cms.double( -0.0125 ),
         e4e1_a_endcap = cms.double( 0.02 ),
         cThreshold_double = cms.double( 10.0 )
       ),
       logWarningEtThreshold_EE_FE = cms.double( 50.0 ),
       eeDetIdToBeRecovered = cms.InputTag( 'hltEcalDetIdToBeRecovered','eeDetId' ),
       recoverEBFE = cms.bool( True ),
       eeFEToBeRecovered = cms.InputTag( 'hltEcalDetIdToBeRecovered','eeFE' ),
       ebDetIdToBeRecovered = cms.InputTag( 'hltEcalDetIdToBeRecovered','ebDetId' ),
       singleChannelRecoveryThreshold = cms.double( 8.0 ),
       sum8ChannelRecoveryThreshold = cms.double( 0.0 ),
       bdtWeightFileNoCracks = cms.FileInPath( "RecoLocalCalo/EcalDeadChannelRecoveryAlgos/data/BDTWeights/bdtgAllRH_8GT700MeV_noCracks_ZskimData2017_v1.xml" ),
       bdtWeightFileCracks = cms.FileInPath( "RecoLocalCalo/EcalDeadChannelRecoveryAlgos/data/BDTWeights/bdtgAllRH_8GT700MeV_onlyCracks_ZskimData2017_v1.xml" ),
       ChannelStatusToBeExcluded = cms.vstring(  ),
       EBrechitCollection = cms.string( "EcalRecHitsEB" ),
       triggerPrimitiveDigiCollection = cms.InputTag( 'hltEcalDigisLegacy','EcalTriggerPrimitives' ),
       recoverEEFE = cms.bool( True ),
       singleChannelRecoveryMethod = cms.string( "NeuralNetworks" ),
       EBLaserMAX = cms.double( 3.0 ),
       flagsMapDBReco = cms.PSet( 
         kGood = cms.vstring( 'kOk',
           'kDAC',
           'kNoLaser',
           'kNoisy' ),
         kNeighboursRecovered = cms.vstring( 'kFixedG0',
           'kNonRespondingIsolated',
           'kDeadVFE' ),
         kDead = cms.vstring( 'kNoDataNoTP' ),
         kNoisy = cms.vstring( 'kNNoisy',
           'kFixedG6',
           'kFixedG1' ),
         kTowerRecovered = cms.vstring( 'kDeadFE' )
       ),
       EBuncalibRecHitCollection = cms.InputTag( 'hltEcalUncalibRecHit','EcalUncalibRecHitsEB' ),
       algoRecover = cms.string( "EcalRecHitWorkerRecover" ),
       algo = cms.string( "EcalRecHitWorkerSimple" ),
       EELaserMAX = cms.double( 8.0 ),
       logWarningEtThreshold_EB_FE = cms.double( 50.0 ),
       recoverEEIsolatedChannels = cms.bool( False ),
       skipTimeCalib = cms.bool( False )
   ),
 )
process.hltHbhereco = SwitchProducerCUDA(
   cpu = cms.EDProducer( "HBHEPhase1Reconstructor",
       digiLabelQIE8 = cms.InputTag( "hltHcalDigis" ),
       processQIE8 = cms.bool( False ),
       digiLabelQIE11 = cms.InputTag( "hltHcalDigis" ),
       processQIE11 = cms.bool( True ),
       tsFromDB = cms.bool( False ),
       recoParamsFromDB = cms.bool( True ),
       saveEffectivePedestal = cms.bool( True ),
       dropZSmarkedPassed = cms.bool( True ),
       makeRecHits = cms.bool( True ),
       saveInfos = cms.bool( False ),
       saveDroppedInfos = cms.bool( False ),
       use8ts = cms.bool( True ),
       sipmQTSShift = cms.int32( 0 ),
       sipmQNTStoSum = cms.int32( 3 ),
       algorithm = cms.PSet( 
         ts4Thresh = cms.double( 0.0 ),
         meanTime = cms.double( 0.0 ),
         nnlsThresh = cms.double( 1.0E-11 ),
         nMaxItersMin = cms.int32( 50 ),
         timeSigmaSiPM = cms.double( 2.5 ),
         applyTimeSlew = cms.bool( True ),
         timeSlewParsType = cms.int32( 3 ),
         ts4Max = cms.vdouble( 100.0, 20000.0, 30000.0 ),
         samplesToAdd = cms.int32( 2 ),
         deltaChiSqThresh = cms.double( 0.001 ),
         applyTimeConstraint = cms.bool( False ),
         timeSigmaHPD = cms.double( 5.0 ),
         useMahi = cms.bool( True ),
         correctForPhaseContainment = cms.bool( True ),
         respCorrM3 = cms.double( 1.0 ),
         pulseJitter = cms.double( 1.0 ),
         applyPedConstraint = cms.bool( False ),
         fitTimes = cms.int32( 1 ),
         nMaxItersNNLS = cms.int32( 500 ),
         applyTimeSlewM3 = cms.bool( True ),
         meanPed = cms.double( 0.0 ),
         ts4Min = cms.double( 0.0 ),
         applyPulseJitter = cms.bool( False ),
         useM2 = cms.bool( False ),
         timeMin = cms.double( -12.5 ),
         useM3 = cms.bool( False ),
         chiSqSwitch = cms.double( -1.0 ),
         dynamicPed = cms.bool( False ),
         tdcTimeShift = cms.double( 0.0 ),
         correctionPhaseNS = cms.double( 6.0 ),
         firstSampleShift = cms.int32( 0 ),
         activeBXs = cms.vint32( -3, -2, -1, 0, 1, 2, 3, 4 ),
         ts4chi2 = cms.vdouble( 15.0, 15.0 ),
         timeMax = cms.double( 12.5 ),
         Class = cms.string( "SimpleHBHEPhase1Algo" ),
         calculateArrivalTime = cms.bool( False ),
         applyLegacyHBMCorrection = cms.bool( False )
       ),
       algoConfigClass = cms.string( "" ),
       setNegativeFlagsQIE8 = cms.bool( False ),
       setNegativeFlagsQIE11 = cms.bool( False ),
       setNoiseFlagsQIE8 = cms.bool( False ),
       setNoiseFlagsQIE11 = cms.bool( False ),
       setPulseShapeFlagsQIE8 = cms.bool( False ),
       setPulseShapeFlagsQIE11 = cms.bool( False ),
       setLegacyFlagsQIE8 = cms.bool( False ),
       setLegacyFlagsQIE11 = cms.bool( False ),
       flagParametersQIE8 = cms.PSet( 
         hitEnergyMinimum = cms.double( 1.0 ),
         pulseShapeParameterSets = cms.VPSet( 
           cms.PSet(  pulseShapeParameters = cms.vdouble( 0.0, 100.0, -50.0, 0.0, -15.0, 0.15 )           ),
           cms.PSet(  pulseShapeParameters = cms.vdouble( 100.0, 2000.0, -50.0, 0.0, -5.0, 0.05 )           ),
           cms.PSet(  pulseShapeParameters = cms.vdouble( 2000.0, 1000000.0, -50.0, 0.0, 95.0, 0.0 )           ),
           cms.PSet(  pulseShapeParameters = cms.vdouble( -1000000.0, 1000000.0, 45.0, 0.1, 1000000.0, 0.0 )           )
         ),
         nominalPedestal = cms.double( 3.0 ),
         hitMultiplicityThreshold = cms.int32( 17 )
       ),
       flagParametersQIE11 = cms.PSet(  ),
       pulseShapeParametersQIE8 = cms.PSet( 
         UseDualFit = cms.bool( True ),
         LinearCut = cms.vdouble( -3.0, -0.054, -0.054 ),
         TriangleIgnoreSlow = cms.bool( False ),
         TS4TS5LowerThreshold = cms.vdouble( 100.0, 120.0, 160.0, 200.0, 300.0, 500.0 ),
         LinearThreshold = cms.vdouble( 20.0, 100.0, 100000.0 ),
         RightSlopeSmallCut = cms.vdouble( 1.08, 1.16, 1.16 ),
         TS4TS5UpperThreshold = cms.vdouble( 70.0, 90.0, 100.0, 400.0 ),
         TS3TS4ChargeThreshold = cms.double( 70.0 ),
         R45PlusOneRange = cms.double( 0.2 ),
         TS4TS5LowerCut = cms.vdouble( -1.0, -0.7, -0.5, -0.4, -0.3, 0.1 ),
         RightSlopeThreshold = cms.vdouble( 250.0, 400.0, 100000.0 ),
         TS3TS4UpperChargeThreshold = cms.double( 20.0 ),
         MinimumChargeThreshold = cms.double( 20.0 ),
         RightSlopeCut = cms.vdouble( 5.0, 4.15, 4.15 ),
         RMS8MaxThreshold = cms.vdouble( 20.0, 100.0, 100000.0 ),
         MinimumTS4TS5Threshold = cms.double( 100.0 ),
         LeftSlopeThreshold = cms.vdouble( 250.0, 500.0, 100000.0 ),
         TS5TS6ChargeThreshold = cms.double( 70.0 ),
         TrianglePeakTS = cms.uint32( 10000 ),
         TS5TS6UpperChargeThreshold = cms.double( 20.0 ),
         RightSlopeSmallThreshold = cms.vdouble( 150.0, 200.0, 100000.0 ),
         RMS8MaxCut = cms.vdouble( -13.5, -11.5, -11.5 ),
         TS4TS5ChargeThreshold = cms.double( 70.0 ),
         R45MinusOneRange = cms.double( 0.2 ),
         LeftSlopeCut = cms.vdouble( 5.0, 2.55, 2.55 ),
         TS4TS5UpperCut = cms.vdouble( 1.0, 0.8, 0.75, 0.72 )
       ),
       pulseShapeParametersQIE11 = cms.PSet(  )
   ),
 )
process.hltSiPixelDigis = SwitchProducerCUDA(
   cpu = cms.EDAlias(
       hltSiPixelDigisLegacy = cms.VPSet( 
         cms.PSet(  type = cms.string( "DetIdedmEDCollection" )         ),
         cms.PSet(  type = cms.string( "SiPixelRawDataErroredmDetSetVector" )         ),
         cms.PSet(  type = cms.string( "PixelFEDChanneledmNewDetSetVector" )         )
       )
   ),
 )
process.hltSiPixelClusters = SwitchProducerCUDA(
   cpu = cms.EDAlias(
       hltSiPixelClustersLegacy = cms.VPSet( 
         cms.PSet(  type = cms.string( "SiPixelClusteredmNewDetSetVector" )         )
       )
   ),
 )
process.hltSiPixelRecHits = SwitchProducerCUDA(
   cpu = cms.EDAlias(
       hltSiPixelRecHitSoA = cms.VPSet( 
         cms.PSet(  type = cms.string( "SiPixelRecHitedmNewDetSetVector" )         ),
         cms.PSet(  type = cms.string( "uintAsHostProduct" )         )
       )
   ),
 )
process.hltPixelTracksSoA = SwitchProducerCUDA(
   cpu = cms.EDProducer( "CAHitNtupletCUDA",
       onGPU = cms.bool( False ),
       pixelRecHitSrc = cms.InputTag( "hltSiPixelRecHitSoA" ),
       ptmin = cms.double( 0.899999976158 ),
       CAThetaCutBarrel = cms.double( 0.00200000009499 ),
       CAThetaCutForward = cms.double( 0.00300000002608 ),
       hardCurvCut = cms.double( 0.0328407224959 ),
       dcaCutInnerTriplet = cms.double( 0.15000000596 ),
       dcaCutOuterTriplet = cms.double( 0.25 ),
       earlyFishbone = cms.bool( True ),
       lateFishbone = cms.bool( False ),
       idealConditions = cms.bool( False ),
       fillStatistics = cms.bool( False ),
       minHitsPerNtuplet = cms.uint32( 3 ),
       maxNumberOfDoublets = cms.uint32( 524288 ),
       minHitsForSharingCut = cms.uint32( 10 ),
       includeJumpingForwardDoublets = cms.bool( True ),
       fitNas4 = cms.bool( False ),
       doClusterCut = cms.bool( True ),
       doZ0Cut = cms.bool( True ),
       doPtCut = cms.bool( True ),
       useRiemannFit = cms.bool( False ),
       doSharedHitCut = cms.bool( True ),
       dupPassThrough = cms.bool( False ),
       useSimpleTripletCleaner = cms.bool( True ),
       trackQualityCuts = cms.PSet( 
         chi2MaxPt = cms.double( 10.0 ),
         tripletMaxTip = cms.double( 0.3 ),
         chi2Scale = cms.double( 8.0 ),
         quadrupletMaxTip = cms.double( 0.5 ),
         quadrupletMinPt = cms.double( 0.3 ),
         quadrupletMaxZip = cms.double( 12.0 ),
         tripletMaxZip = cms.double( 12.0 ),
         tripletMinPt = cms.double( 0.5 ),
         chi2Coeff = cms.vdouble( 0.9, 1.8 )
       )
   ),
 )
process.hltPixelVerticesSoA = SwitchProducerCUDA(
   cpu = cms.EDProducer( "PixelVertexProducerCUDA",
       onGPU = cms.bool( False ),
       oneKernel = cms.bool( True ),
       useDensity = cms.bool( True ),
       useDBSCAN = cms.bool( False ),
       useIterative = cms.bool( False ),
       minT = cms.int32( 2 ),
       eps = cms.double( 0.07 ),
       errmax = cms.double( 0.01 ),
       chi2max = cms.double( 9.0 ),
       PtMin = cms.double( 0.5 ),
       PtMax = cms.double( 75.0 ),
       pixelTrackSrc = cms.InputTag( "hltPixelTracksSoA" )
   ),
 )

process.HLTL1UnpackerSequence = cms.Sequence( process.hltGtStage2Digis + process.hltGtStage2ObjectMap )
process.HLTBeamSpot = cms.Sequence( process.hltScalersRawToDigi + process.hltOnlineMetaDataDigis + process.hltOnlineBeamSpot )
process.HLTBeginSequence = cms.Sequence( process.hltTriggerType + process.HLTL1UnpackerSequence + process.HLTBeamSpot )
process.HLTDoFullUnpackingEgammaEcalWithoutPreshowerSequence = cms.Sequence( process.hltEcalDigisLegacy + process.hltEcalDigis + process.hltEcalUncalibRecHit + process.hltEcalDetIdToBeRecovered + process.hltEcalRecHit )
process.HLTDoLocalHcalSequence = cms.Sequence( process.hltHcalDigis + process.hltHbhereco + process.hltHfprereco + process.hltHfreco + process.hltHoreco )
process.HLTDoCaloSequence = cms.Sequence( process.HLTDoFullUnpackingEgammaEcalWithoutPreshowerSequence + process.HLTDoLocalHcalSequence + process.hltTowerMakerForAll )
process.HLTAK4CaloJetsReconstructionSequence = cms.Sequence( process.HLTDoCaloSequence + process.hltAK4CaloJets + process.hltAK4CaloJetsIDPassed )
process.HLTAK4CaloCorrectorProducersSequence = cms.Sequence( process.hltAK4CaloFastJetCorrector + process.hltAK4CaloRelativeCorrector + process.hltAK4CaloAbsoluteCorrector + process.hltAK4CaloResidualCorrector + process.hltAK4CaloCorrector )
process.HLTAK4CaloJetsCorrectionSequence = cms.Sequence( process.hltFixedGridRhoFastjetAllCalo + process.HLTAK4CaloCorrectorProducersSequence + process.hltAK4CaloJetsCorrected + process.hltAK4CaloJetsCorrectedIDPassed )
process.HLTAK4CaloJetsSequence = cms.Sequence( process.HLTAK4CaloJetsReconstructionSequence + process.HLTAK4CaloJetsCorrectionSequence )
process.HLTDoCaloSequencePF = cms.Sequence( process.HLTDoFullUnpackingEgammaEcalWithoutPreshowerSequence + process.HLTDoLocalHcalSequence + process.hltTowerMakerForAll )
process.HLTAK4CaloJetsPrePFRecoSequence = cms.Sequence( process.HLTDoCaloSequencePF + process.hltAK4CaloJetsPF )
process.HLTPreAK4PFJetsRecoSequence = cms.Sequence( process.HLTAK4CaloJetsPrePFRecoSequence + process.hltAK4CaloJetsPFEt5 )
process.HLTMuonLocalRecoSequence = cms.Sequence( process.hltMuonDTDigis + process.hltDt1DRecHits + process.hltDt4DSegments + process.hltMuonCSCDigis + process.hltCsc2DRecHits + process.hltCscSegments + process.hltMuonRPCDigis + process.hltRpcRecHits + process.hltMuonGEMDigis + process.hltGemRecHits + process.hltGemSegments )
process.HLTL2muonrecoNocandSequence = cms.Sequence( process.HLTMuonLocalRecoSequence + process.hltL2OfflineMuonSeeds + process.hltL2MuonSeeds + process.hltL2Muons )
process.HLTL2muonrecoSequence = cms.Sequence( process.HLTL2muonrecoNocandSequence + process.hltL2MuonCandidates )
process.HLTDoLocalPixelSequence = cms.Sequence( process.hltSiPixelDigisLegacy + process.hltSiPixelDigis + process.hltSiPixelClustersLegacy + process.hltSiPixelClusters + process.hltSiPixelClustersCache + process.hltSiPixelRecHitSoA + process.hltSiPixelRecHits )
process.HLTDoLocalStripSequence = cms.Sequence( process.hltSiStripExcludedFEDListProducer + process.hltSiStripRawToClustersFacility + process.hltSiStripClusters )
process.HLTIterL3OImuonTkCandidateSequence = cms.Sequence( process.hltIterL3OISeedsFromL2Muons + process.hltIterL3OITrackCandidates + process.hltIterL3OIMuCtfWithMaterialTracks + process.hltIterL3OIMuonTrackCutClassifier + process.hltIterL3OIMuonTrackSelectionHighPurity + process.hltL3MuonsIterL3OI )
process.HLTRecoPixelTracksSequence = cms.Sequence( process.hltPixelTracksTrackingRegions + process.hltPixelTracksSoA + process.hltPixelTracks )
process.HLTRecopixelvertexingSequence = cms.Sequence( process.hltPixelTracksFitter + process.hltPixelTracksFilter + process.HLTRecoPixelTracksSequence + process.hltPixelVerticesSoA + process.hltPixelVertices + process.hltTrimmedPixelVertices )
process.HLTIterL3MuonRecopixelvertexingSequence = cms.Sequence( process.HLTRecopixelvertexingSequence + process.hltIterL3MuonPixelTracksTrackingRegions + process.hltPixelTracksInRegionL2 )
process.HLTIterativeTrackingIteration0ForIterL3Muon = cms.Sequence( process.hltIter0IterL3MuonPixelSeedsFromPixelTracks + process.hltIter0IterL3MuonPixelSeedsFromPixelTracksFiltered + process.hltIter0IterL3MuonCkfTrackCandidates + process.hltIter0IterL3MuonCtfWithMaterialTracks + process.hltIter0IterL3MuonTrackCutClassifier + process.hltIter0IterL3MuonTrackSelectionHighPurity )
process.HLTIterL3IOmuonTkCandidateSequence = cms.Sequence( process.HLTIterL3MuonRecopixelvertexingSequence + process.HLTIterativeTrackingIteration0ForIterL3Muon + process.hltL3MuonsIterL3IO )
process.HLTIterL3OIAndIOFromL2muonTkCandidateSequence = cms.Sequence( process.HLTIterL3OImuonTkCandidateSequence + process.hltIterL3OIL3MuonsLinksCombination + process.hltIterL3OIL3Muons + process.hltIterL3OIL3MuonCandidates + process.hltL2SelectorForL3IO + process.HLTIterL3IOmuonTkCandidateSequence + process.hltIterL3MuonsFromL2LinksCombination )
process.HLTRecopixelvertexingSequenceForIterL3FromL1Muon = cms.Sequence( process.HLTRecopixelvertexingSequence + process.hltIterL3FromL1MuonPixelTracksTrackingRegions + process.hltPixelTracksInRegionL1 )
process.HLTIterativeTrackingIteration0ForIterL3FromL1Muon = cms.Sequence( process.hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks + process.hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracksFiltered + process.hltIter0IterL3FromL1MuonCkfTrackCandidates + process.hltIter0IterL3FromL1MuonCtfWithMaterialTracks + process.hltIter0IterL3FromL1MuonTrackCutClassifier + process.hltIter0IterL3FromL1MuonTrackSelectionHighPurity )
process.HLTIterL3IOmuonFromL1TkCandidateSequence = cms.Sequence( process.HLTRecopixelvertexingSequenceForIterL3FromL1Muon + process.HLTIterativeTrackingIteration0ForIterL3FromL1Muon )
process.HLTIterL3muonTkCandidateSequence = cms.Sequence( process.HLTDoLocalPixelSequence + process.HLTDoLocalStripSequence + process.HLTIterL3OIAndIOFromL2muonTkCandidateSequence + process.hltL1MuonsPt0 + process.HLTIterL3IOmuonFromL1TkCandidateSequence )
process.HLTL3muonrecoNocandSequence = cms.Sequence( process.HLTIterL3muonTkCandidateSequence + process.hltIterL3MuonMerged + process.hltIterL3MuonAndMuonFromL1Merged + process.hltIterL3GlbMuon + process.hltIterL3MuonsNoID + process.hltIterL3Muons + process.hltL3MuonsIterL3Links + process.hltIterL3MuonTracks )
process.HLTL3muonrecoSequence = cms.Sequence( process.HLTL3muonrecoNocandSequence + process.hltIterL3MuonCandidates )
process.HLTTrackReconstructionForPixelOnlyPF = cms.Sequence( process.HLTDoLocalPixelSequence + process.HLTRecopixelvertexingSequence + process.hltPixelTracksZetaClean + process.hltPixelOnlyPFMuonMerging + process.hltPixelOnlyMuonLinks + process.hltPixelOnlyMuons )
process.HLTPreshowerSequence = cms.Sequence( process.hltEcalPreshowerDigis + process.hltEcalPreshowerRecHit )
process.HLTPixelOnlyParticleFlowSequence = cms.Sequence( process.HLTPreshowerSequence + process.hltParticleFlowRecHitECALUnseeded + process.hltParticleFlowRecHitHBHE + process.hltParticleFlowRecHitHF + process.hltParticleFlowRecHitPSUnseeded + process.hltParticleFlowClusterECALUncorrectedUnseeded + process.hltParticleFlowClusterPSUnseeded + process.hltParticleFlowClusterECALUnseeded + process.hltParticleFlowClusterHBHE + process.hltParticleFlowClusterHCAL + process.hltParticleFlowClusterHF + process.hltLightPixelOnlyPFTracks + process.hltPixelOnlyParticleFlowBlock + process.hltPixelOnlyParticleFlow )
process.HLTAK4PixelOnlyPFJetsReconstructionSequence = cms.Sequence( process.HLTL2muonrecoSequence + process.HLTL3muonrecoSequence + process.HLTTrackReconstructionForPixelOnlyPF + process.HLTPixelOnlyParticleFlowSequence + process.hltAK4PixelOnlyPFJets + process.hltAK4PixelOnlyPFJetsLooseID + process.hltAK4PixelOnlyPFJetsTightID )
process.HLTAK4PixelOnlyPFCorrectorProducersSequence = cms.Sequence( process.hltAK4PixelOnlyPFFastJetCorrector + process.hltAK4PixelOnlyPFRelativeCorrector + process.hltAK4PixelOnlyPFAbsoluteCorrector + process.hltAK4PixelOnlyPFResidualCorrector + process.hltAK4PixelOnlyPFCorrector )
process.HLTAK4PixelOnlyPFJetsCorrectionSequence = cms.Sequence( process.hltFixedGridRhoFastjetPixelOnlyAll + process.HLTAK4PixelOnlyPFCorrectorProducersSequence + process.hltAK4PixelOnlyPFJetsCorrected + process.hltAK4PixelOnlyPFJetsLooseIDCorrected + process.hltAK4PixelOnlyPFJetsTightIDCorrected )
process.HLTAK4PixelOnlyPFJetsSequence = cms.Sequence( process.HLTPreAK4PFJetsRecoSequence + process.HLTAK4PixelOnlyPFJetsReconstructionSequence + process.HLTAK4PixelOnlyPFJetsCorrectionSequence )
process.HLTDoFullUnpackingEgammaEcalMFSequence = cms.Sequence( process.HLTDoFullUnpackingEgammaEcalWithoutPreshowerSequence + process.HLTPreshowerSequence )
process.HLTPFClusteringEcalMFForMuons = cms.Sequence( process.hltRecHitInRegionForMuonsMF + process.hltRecHitInRegionForMuonsES + process.hltParticleFlowRecHitECALForMuonsMF + process.hltParticleFlowRecHitPSForMuons + process.hltParticleFlowClusterECALUncorrectedForMuonsMF + process.hltParticleFlowClusterPSForMuons + process.hltParticleFlowClusterECALForMuonsMF )
process.HLTL3muonEcalPFisorecoSequenceNoBoolsForMuons = cms.Sequence( process.HLTDoFullUnpackingEgammaEcalMFSequence + process.HLTDoLocalHcalSequence + process.hltFixedGridRhoFastjetECALMFForMuons + process.hltFixedGridRhoFastjetHCAL + process.HLTPFClusteringEcalMFForMuons + process.hltMuonEcalMFPFClusterIsoForMuons )
process.HLTPFHcalClustering = cms.Sequence( process.hltParticleFlowRecHitHBHE + process.hltParticleFlowClusterHBHE + process.hltParticleFlowClusterHCAL )
process.HLTL3muonHcalPFisorecoSequenceNoBoolsForMuons = cms.Sequence( process.HLTPFHcalClustering + process.hltMuonHcalRegPFClusterIsoForMuons )
process.HLTIterativeTrackingL3MuonIteration0 = cms.Sequence( process.hltPixelTracksTrackingRegionsForSeedsL3Muon + process.hltPixelTracksInRegionIter0L3Muon + process.hltIter0L3MuonPixelSeedsFromPixelTracks + process.hltIter0L3MuonCkfTrackCandidates + process.hltIter0L3MuonCtfWithMaterialTracks + process.hltIter0L3MuonTrackCutClassifier + process.hltIter0L3MuonTrackSelectionHighPurity )
process.HLTTrackReconstructionForIsoL3MuonIter02 = cms.Sequence( process.HLTDoLocalPixelSequence + process.HLTDoLocalStripSequence + process.HLTIterativeTrackingL3MuonIteration0 )
process.HLTMuIsolationSequence = cms.Sequence( process.HLTL3muonEcalPFisorecoSequenceNoBoolsForMuons + process.HLTL3muonHcalPFisorecoSequenceNoBoolsForMuons + process.HLTTrackReconstructionForIsoL3MuonIter02 + process.hltMuonTkRelIsolationCut0p08Map )
process.HLTDoFullUnpackingEgammaEcalSequence = cms.Sequence( process.HLTDoFullUnpackingEgammaEcalWithoutPreshowerSequence + process.HLTPreshowerSequence )
process.HLTPFClusteringForEgamma = cms.Sequence( process.hltRechitInRegionsECAL + process.hltRechitInRegionsES + process.hltParticleFlowRecHitECALL1Seeded + process.hltParticleFlowRecHitPSL1Seeded + process.hltParticleFlowClusterPSL1Seeded + process.hltParticleFlowClusterECALUncorrectedL1Seeded + process.hltParticleFlowClusterECALL1Seeded + process.hltParticleFlowSuperClusterECALL1Seeded )
process.HLTFastJetForEgamma = cms.Sequence( process.hltTowerMakerForAll + process.hltFixedGridRhoFastjetAllCaloForMuons )
process.HLTElePixelMatchSequence = cms.Sequence( process.HLTDoLocalPixelSequence + process.HLTDoLocalStripSequence + process.hltPixelLayerPairs + process.hltPixelLayerTriplets + process.hltEgammaHoverE + process.hltEgammaSuperClustersToPixelMatch + process.hltEleSeedsTrackingRegions + process.hltElePixelHitDoublets + process.hltElePixelHitDoubletsForTriplets + process.hltElePixelHitTriplets + process.hltElePixelSeedsDoublets + process.hltElePixelSeedsTriplets + process.hltElePixelSeedsCombined + process.hltEgammaElectronPixelSeeds + process.hltEgammaPixelMatchVars )
process.HLTGsfElectronSequence = cms.Sequence( process.hltEgammaCkfTrackCandidatesForGSF + process.hltEgammaGsfTracks + process.hltEgammaGsfElectrons + process.hltEgammaGsfTrackVars )
process.HLTTrackReconstructionForPixelOnlyPFNoMu = cms.Sequence( process.HLTDoLocalPixelSequence + process.HLTRecopixelvertexingSequence )
process.HLTTrackReconstructionForIsoElectronPixelOnly = cms.Sequence( process.HLTPreAK4PFJetsRecoSequence + process.HLTTrackReconstructionForPixelOnlyPFNoMu + process.hltPixelTracksZetaClean )
process.HLTPixelOnlyPFScoutingSequence = cms.Sequence( process.HLTAK4PixelOnlyPFJetsSequence + process.hltPixelOnlyPFMETProducer + process.HLTMuIsolationSequence + process.HLTDoFullUnpackingEgammaEcalSequence + process.HLTPFClusteringForEgamma + process.hltEgammaCandidates + process.hltEgammaClusterShape + process.HLTDoLocalHcalSequence + process.HLTFastJetForEgamma + process.hltEgammaHoverE + process.hltEgammaEcalPFClusterIso + process.HLTPFHcalClustering + process.hltEgammaHcalPFClusterIso + process.HLTElePixelMatchSequence + process.HLTGsfElectronSequence + process.HLTTrackReconstructionForIsoElectronPixelOnly + process.hltEgammaEleGsfTrackIsoPixelOnly )


process.HLTIterativeTrackingIteration0 = cms.Sequence( process.hltIter0PFLowPixelSeedsFromPixelTracks + process.hltIter0PFlowCkfTrackCandidates + process.hltIter0PFlowCtfWithMaterialTracks + process.hltIter0PFlowTrackCutClassifier + process.hltMergedTracks )
process.HLTIterativeTrackingIter02 = cms.Sequence( process.HLTIterativeTrackingIteration0 )
process.HLTTrackReconstructionForPF = cms.Sequence( process.HLTDoLocalPixelSequence + process.HLTRecopixelvertexingSequence + process.HLTDoLocalStripSequence + process.HLTIterativeTrackingIter02 + process.hltPFMuonMerging + process.hltMuonLinks + process.hltMuons )
process.HLTParticleFlowSequence = cms.Sequence( process.HLTPreshowerSequence + process.hltParticleFlowRecHitECALUnseeded + process.hltParticleFlowRecHitHBHE + process.hltParticleFlowRecHitHF + process.hltParticleFlowRecHitPSUnseeded + process.hltParticleFlowClusterECALUncorrectedUnseeded + process.hltParticleFlowClusterPSUnseeded + process.hltParticleFlowClusterECALUnseeded + process.hltParticleFlowClusterHBHE + process.hltParticleFlowClusterHCAL + process.hltParticleFlowClusterHF + process.hltLightPFTracks + process.hltParticleFlowBlock + process.hltParticleFlow )
process.HLTAK4PFJetsReconstructionSequence = cms.Sequence( process.HLTL2muonrecoSequence + process.HLTL3muonrecoSequence + process.HLTTrackReconstructionForPF + process.HLTParticleFlowSequence + process.hltAK4PFJets + process.hltAK4PFJetsLooseID + process.hltAK4PFJetsTightID )
process.HLTAK4PFCorrectorProducersSequence = cms.Sequence( process.hltAK4PFFastJetCorrector + process.hltAK4PFRelativeCorrector + process.hltAK4PFAbsoluteCorrector + process.hltAK4PFResidualCorrector + process.hltAK4PFCorrector )
process.HLTAK4PFJetsCorrectionSequence = cms.Sequence( process.hltFixedGridRhoFastjetAll + process.HLTAK4PFCorrectorProducersSequence + process.hltAK4PFJetsCorrected + process.hltAK4PFJetsLooseIDCorrected + process.hltAK4PFJetsTightIDCorrected )
process.HLTAK4PFJetsSequence = cms.Sequence( process.HLTPreAK4PFJetsRecoSequence + process.HLTAK4PFJetsReconstructionSequence + process.HLTAK4PFJetsCorrectionSequence )
process.HLTL2muonrecoSequenceNoVtx = cms.Sequence( process.HLTL2muonrecoNocandSequence + process.hltL2MuonCandidatesNoVtx )
process.HLTIterL3OImuonTkCandidateSequenceNoVtx = cms.Sequence( process.hltIterL3OISeedsFromL2MuonsNoVtx + process.hltIterL3OITrackCandidatesNoVtx + process.hltIterL3OIMuCtfWithMaterialTracksNoVtx + process.hltIterL3OIMuonTrackCutClassifierNoVtx + process.hltIterL3OIMuonTrackSelectionHighPurityNoVtx + process.hltL3MuonsIterL3OINoVtx )
process.HLTIterL3MuonRecoPixelTracksSequenceNoVtx = cms.Sequence( process.HLTRecopixelvertexingSequence + process.hltIterL3MuonPixelTracksTrackingRegionsNoVtx + process.hltPixelTracksInRegionL2NoVtx )
process.HLTIterL3MuonRecopixelvertexingSequenceNoVtx = cms.Sequence( process.HLTIterL3MuonRecoPixelTracksSequenceNoVtx )
process.HLTIterativeTrackingIteration0ForIterL3MuonNoVtx = cms.Sequence( process.hltIter0IterL3MuonPixelSeedsFromPixelTracksNoVtx + process.hltIter0IterL3MuonCkfTrackCandidatesNoVtx + process.hltIter0IterL3MuonCtfWithMaterialTracksNoVtx + process.hltIter0IterL3MuonTrackCutClassifierNoVtx + process.hltIter0IterL3MuonTrackSelectionHighPurityNoVtx )
process.HLTIterL3IOmuonTkCandidateSequenceNoVtx = cms.Sequence( process.HLTIterL3MuonRecopixelvertexingSequenceNoVtx + process.HLTIterativeTrackingIteration0ForIterL3MuonNoVtx + process.hltL3MuonsIterL3IONoVtx )
process.HLTIterL3OIAndIOFromL2muonTkCandidateSequenceNoVtx = cms.Sequence( process.HLTIterL3OImuonTkCandidateSequenceNoVtx + process.hltIterL3OIL3MuonsLinksCombinationNoVtx + process.hltIterL3OIL3MuonsNoVtx + process.hltIterL3OIL3MuonCandidatesNoVtx + process.hltL2SelectorForL3IONoVtx + process.HLTIterL3IOmuonTkCandidateSequenceNoVtx + process.hltIterL3MuonsFromL2LinksCombinationNoVtx + process.hltIterL3MuonsFromL2NoVtx )
process.HLTRecopixelvertexingSequenceForIterL3FromL1MuonNoVtx = cms.Sequence( process.HLTRecopixelvertexingSequence + process.hltIterL3FromL1MuonPixelTracksTrackingRegionsNoVtx + process.hltPixelTracksInRegionL1NoVtx )
process.HLTIterativeTrackingIteration0ForIterL3FromL1MuonNoVtx = cms.Sequence( process.hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracksNoVtx + process.hltIter0IterL3FromL1MuonCkfTrackCandidatesNoVtx + process.hltIter0IterL3FromL1MuonCtfWithMaterialTracksNoVtx + process.hltIter0IterL3FromL1MuonTrackCutClassifierNoVtx + process.hltIter0IterL3FromL1MuonTrackSelectionHighPurityNoVtx )
process.HLTIterL3IOmuonFromL1TkCandidateSequenceNoVtx = cms.Sequence( process.HLTRecopixelvertexingSequenceForIterL3FromL1MuonNoVtx + process.HLTIterativeTrackingIteration0ForIterL3FromL1MuonNoVtx )
process.HLTIterL3muonTkCandidateSequenceNoVtx = cms.Sequence( process.HLTDoLocalPixelSequence + process.HLTDoLocalStripSequence + process.HLTIterL3OIAndIOFromL2muonTkCandidateSequenceNoVtx + process.hltIterL3MuonL1MuonNoL2SelectorNoVtx + process.HLTIterL3IOmuonFromL1TkCandidateSequenceNoVtx )
process.HLTL3muonrecoNocandSequenceNoVtx = cms.Sequence( process.HLTIterL3muonTkCandidateSequenceNoVtx + process.hltIterL3MuonMergedNoVtx + process.hltIterL3MuonAndMuonFromL1MergedNoVtx + process.hltL3MuonsIterL3LinksNoVtx + process.hltIterL3MuonsNoVtx )
process.HLTL3muonrecoSequenceNoVtx = cms.Sequence( process.HLTL3muonrecoNocandSequenceNoVtx + process.hltIterL3MuonCandidatesNoVtx )
process.HLTPFClusteringEcalMFForMuonsNoVtx = cms.Sequence( process.hltRecHitInRegionForMuonsMFnoVtx + process.hltRecHitInRegionForMuonsESNoVtx + process.hltParticleFlowRecHitECALForMuonsMFNoVtx + process.hltParticleFlowRecHitPSForMuonsNoVtx + process.hltParticleFlowClusterECALUncorrectedForMuonsMFNoVtx + process.hltParticleFlowClusterPSForMuonsNoVtx + process.hltParticleFlowClusterECALForMuonsMFNoVtx )
process.HLTL3muonEcalPFisorecoSequenceNoBoolsForMuonsNoVtx = cms.Sequence( process.HLTDoFullUnpackingEgammaEcalMFSequence + process.HLTDoLocalHcalSequence + process.hltFixedGridRhoFastjetECALMFForMuons + process.hltFixedGridRhoFastjetHCAL + process.HLTPFClusteringEcalMFForMuonsNoVtx + process.hltMuonEcalMFPFClusterIsoForMuonsNoVtx )
process.HLTL3muonHcalPFisorecoSequenceNoBoolsForMuonsNoVtx = cms.Sequence( process.HLTPFHcalClustering + process.hltMuonHcalPFClusterIsoForMuonsNoVtx )
process.HLTIterativeTrackingL3MuonIteration0NoVtx = cms.Sequence( process.hltPixelTracksTrackingRegionsForSeedsL3MuonNoVtx + process.hltPixelTracksInRegionIter0L3MuonNoVtx + process.hltIter0L3MuonPixelSeedsFromPixelTracksNoVtx + process.hltIter0L3MuonCkfTrackCandidatesNoVtx + process.hltIter0L3MuonCtfWithMaterialTracksNoVtx + process.hltIter0L3MuonTrackCutClassifierNoVtx + process.hltIter0L3MuonTrackSelectionHighPurityNoVtx )
process.HLTTrackReconstructionForIsoL3MuonIter02NoVtx = cms.Sequence( process.HLTDoLocalPixelSequence + process.HLTDoLocalStripSequence + process.HLTIterativeTrackingL3MuonIteration0NoVtx )
process.HLTMuIsolationSequenceNoVtx = cms.Sequence( process.HLTL3muonEcalPFisorecoSequenceNoBoolsForMuonsNoVtx + process.HLTL3muonHcalPFisorecoSequenceNoBoolsForMuonsNoVtx + process.HLTTrackReconstructionForIsoL3MuonIter02NoVtx + process.hltMuonTkRelIsolationCut0p09MapNoVtx )

### Full tracking packing sequence
process.HLTFullTrackingPFScoutingPackingSequence = cms.Sequence( process.hltScoutingFullTrackPacker + process.hltScoutingPrimaryVertexPacker + process.hltScoutingFullTrackPFPacker + process.hltScoutingMuonPacker + process.hltScoutingEgammaPacker )
process.HLTEndSequence = cms.Sequence( process.hltBoolEnd )
process.HLTPixelOnlyPFScoutingPackingSequence = cms.Sequence( process.hltScoutingTrackPacker + process.hltScoutingPrimaryVertexPacker + process.hltScoutingPFPacker + process.hltScoutingMuonPacker + process.hltScoutingEgammaPacker )

process.HLTriggerFirstPath = cms.Path( process.hltGetConditions + process.hltGetRaw + process.hltPSetMap + process.hltBoolFalse )

process.hltPFMETProducer = cms.EDProducer( "PFMETProducer",
    src = cms.InputTag( "hltParticleFlow" ),
    globalThreshold = cms.double( 0.0 ),
    alias = cms.string( "hltPFMet" ),
    calculateSignificance = cms.bool( False ),
    parameters = cms.PSet(  ),
    applyWeight = cms.bool( False ),
    srcWeights = cms.InputTag( "" )
)

process.HLTPFScoutingSequence = cms.Sequence( process.HLTAK4PFJetsSequence + process.hltPFMETProducer + process.HLTMuIsolationSequence + process.HLTDoFullUnpackingEgammaEcalSequence + process.HLTPFClusteringForEgamma + process.hltEgammaCandidates + process.hltEgammaClusterShape + process.HLTDoLocalHcalSequence + process.HLTFastJetForEgamma + process.hltEgammaHoverE + process.hltEgammaEcalPFClusterIso + process.HLTPFHcalClustering + process.hltEgammaHcalPFClusterIso + process.HLTElePixelMatchSequence + process.HLTGsfElectronSequence + process.HLTTrackReconstructionForIsoElectronPixelOnly + process.hltEgammaEleGsfTrackIsoPixelOnly )

### Full tracking PF sequence with HLTPFScoutingSequence instead of HLTPixelOnlyPFScoutingSequence
process.DST_Run3_PFScoutingFullTracking_v16 = cms.Path( process.HLTBeginSequence + process.hltL1sDSTRun3PFScoutingPixelTracking + process.hltPreDSTRun3PFScoutingFullTracking + process.HLTAK4CaloJetsSequence + process.HLTPFScoutingSequence + process.hltEgammaR9ID + process.HLTL2muonrecoSequenceNoVtx + process.HLTL3muonrecoSequenceNoVtx + process.hltDisplacedmumuVtxNoMatchingProducer + process.HLTMuIsolationSequenceNoVtx + process.hltFEDSelectorL1 + process.HLTFullTrackingPFScoutingPackingSequence + process.HLTEndSequence )

process.DST_Run3_PFScoutingPixelTracking_v16 = cms.Path( process.HLTBeginSequence + process.hltL1sDSTRun3PFScoutingPixelTracking + process.hltPreDSTRun3PFScoutingPixelTracking + process.HLTAK4CaloJetsSequence + process.HLTPixelOnlyPFScoutingSequence + process.hltEgammaR9ID + process.HLTL2muonrecoSequenceNoVtx + process.HLTL3muonrecoSequenceNoVtx + process.hltDisplacedmumuVtxNoMatchingProducer + process.HLTMuIsolationSequenceNoVtx + process.hltFEDSelectorL1 + process.HLTPixelOnlyPFScoutingPackingSequence + process.HLTEndSequence )

process.HLTriggerFinalPath = cms.Path( process.hltGtStage2Digis + process.hltScalersRawToDigi + process.hltFEDSelector + process.hltTriggerSummaryAOD + process.hltTriggerSummaryRAW + process.hltBoolFalse )


process.schedule = cms.Schedule( *(process.HLTriggerFirstPath, process.DST_Run3_PFScoutingFullTracking_v16, process.DST_Run3_PFScoutingPixelTracking_v16, process.HLTriggerFinalPath, ))


# source module (EDM inputs)
process.source = cms.Source( "PoolSource",
    fileNames = cms.untracked.vstring(
        #'/store/mc/Run3Summer21DRPremix/TT_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW/120X_mcRun3_2021_realistic_v6-v2/2540000/b354245e-d8bc-424d-b527-58815586a6a5.root',
        #'file:/eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/Scouting/Run3/ML_210512/BulkGraviton_hh_GF_HH_14TeV_TuneCP5_pythia8/ML_210512/210610_115224/0000/scouting_1.root',
         'file:/eos/cms/store/relval/CMSSW_12_3_0_pre6/RelValQCD_FlatPt_15_3000HS_14/GEN-SIM-DIGI-RAW/PU_123X_mcRun3_2021_realistic_v11-v1/10000/f9430408-7789-446a-b3f2-39892427afb4.root'
    ),
    inputCommands = cms.untracked.vstring(
        'keep *'
    )
)

# run the Full L1T emulator, then repack the data into a new RAW collection, to be used by the HLT
from HLTrigger.Configuration.CustomConfigs import L1REPACK
process = L1REPACK(process,"uGT")

# add a single "keep *" output
process.hltOutputMinimal = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "output.root" ),
    fastCloning = cms.untracked.bool( False ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string( 'AOD' ),
        filterName = cms.untracked.string( '' )
    ),
    outputCommands = cms.untracked.vstring( 'drop *',
        'keep edmTriggerResults_*_*_*',
        'keep triggerTriggerEvent_*_*_*',
        'keep GlobalAlgBlkBXVector_*_*_*',                  
        'keep GlobalExtBlkBXVector_*_*_*',
        'keep l1tEGammaBXVector_*_EGamma_*',
        'keep l1tEtSumBXVector_*_EtSum_*',
        'keep l1tJetBXVector_*_Jet_*',
        'keep l1tMuonBXVector_*_Muon_*',
        'keep l1tTauBXVector_*_Tau_*',
    )
)
process.MinimalOutput = cms.FinalPath( process.hltOutputMinimal )
process.schedule.append( process.MinimalOutput )

# limit the number of events to be processed
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32( 100 )
)

# enable TrigReport, TimeReport and MultiThreading
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool( True ),
    numberOfThreads = cms.untracked.uint32( 1 ),
    numberOfStreams = cms.untracked.uint32( 0 ),
)

# override the GlobalTag, connection string and pfnPrefix
if 'GlobalTag' in process.__dict__:
    from Configuration.AlCa.GlobalTag import GlobalTag as customiseGlobalTag
    process.GlobalTag = customiseGlobalTag(process.GlobalTag, globaltag = 'auto:phase1_2021_realistic', conditions = 'L1Menu_Collisions2022_v1_0_0_xml,L1TUtmTriggerMenuRcd,,,9999-12-31 23:59:59.000')

# show summaries from trigger analysers used at HLT
if 'MessageLogger' in process.__dict__:
    process.MessageLogger.TriggerSummaryProducerAOD = cms.untracked.PSet()
    process.MessageLogger.L1GtTrigReport = cms.untracked.PSet()
    process.MessageLogger.L1TGlobalSummary = cms.untracked.PSet()
    process.MessageLogger.HLTrigReport = cms.untracked.PSet()
    process.MessageLogger.FastReport = cms.untracked.PSet()
    process.MessageLogger.ThroughputService = cms.untracked.PSet()

# load the DQMStore and DQMRootOutputModule
#process.load( "DQMServices.Core.DQMStore_cfi" )
#
#process.dqmOutput = cms.OutputModule("DQMRootOutputModule",
#    fileName = cms.untracked.string("DQMIO.root")
#)
#
#process.DQMOutput = cms.FinalPath( process.dqmOutput )
#process.schedule.append( process.DQMOutput )

from HLTrigger.Configuration.customizeHLTforCMSSW import customizeHLTforCMSSW
process = customizeHLTforCMSSW(process,"GRun")

# Eras-based customisations
from HLTrigger.Configuration.Eras import modifyHLTforEras
modifyHLTforEras(process)

#User-defined customization functions
from HLTrigger.Configuration.customizeHLTforPatatrack import customizeHLTforPatatrackTriplets
process = customizeHLTforPatatrackTriplets(process)

#scouting event content
#note: the scouting format  of 11_2 looks incompatible with 12_3
#      hence it is dropped (this is bad btw ...)
process.hltOutputMinimal.outputCommands.extend([
     'keep *_hltFEDSelectorL1_*_*',
     'keep *_hltScouting*Packer_*_*',
     'drop *_hltFEDSelectorL1_*_HLT',
     'drop *_hltScoutingEgammaPacker_*_HLT',
     'drop *_hltScoutingMuonPacker_*_HLT',
     'drop *_hltScoutingPFPacker_*_HLT',
     'drop *_hltScoutingPrimaryVertexPacker_*_HLT',
     'drop *_hltScoutingTrackPacker_*_HLT',
    ])


#gen content from RAWMINIAODSIM
process.hltOutputMinimal.outputCommands.extend([
    'keep LHEEventProduct_*_*_*',
    'keep GenFilterInfo_*_*_*',
    'keep GenLumiInfoHeader_generator_*_*',
    'keep GenLumiInfoProduct_*_*_*',
    'keep GenEventInfoProduct_generator_*_*',
    'keep recoGenParticles_genPUProtons_*_*',
    'keep *_slimmedGenJetsFlavourInfos_*_*',
    'keep *_slimmedGenJets__*',
    'keep *_slimmedGenJetsAK8__*',
    'keep *_slimmedGenJetsAK8SoftDropSubJets__*',
    'keep *_genMetTrue_*_*',
    'keep LHERunInfoProduct_*_*_*',
    'keep GenRunInfoProduct_*_*_*',
    'keep *_genParticles_xyz0_*',
    'keep *_genParticles_t0_*',
    'keep PileupSummaryInfos_slimmedAddPileupInfo_*_*',
    'keep L1GtTriggerMenuLite_l1GtTriggerMenuLite__*',
    'keep SimVertexs_g4SimHits_*_*',
    'keep *_packedGenParticles_*_*',
    'keep *_prunedGenParticles_*_*'
])

#to get truth label
process.hltOutputMinimal.outputCommands.extend([
    'keep *_slimmedJets_*_*'
])

#this just allows us to specify an output file
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ('analysis') 
options.parseArguments()
process.hltOutputMinimal.fileName = options.outputFile

#this sets the L1 filter to ignore so the path continues even if it fails
process.DST_Run3_PFScoutingPixelTracking_v16.replace(process.hltL1sDSTRun3PFScoutingPixelTracking,cms.ignore(process.hltL1sDSTRun3PFScoutingPixelTracking))
process.DST_Run3_PFScoutingFullTracking_v16.replace(process.hltL1sDSTRun3PFScoutingPixelTracking,cms.ignore(process.hltL1sDSTRun3PFScoutingPixelTracking))
#this needs to be set to false elsewhere than p5
process.hltOnlineBeamSpot.useTransientRecord = False
