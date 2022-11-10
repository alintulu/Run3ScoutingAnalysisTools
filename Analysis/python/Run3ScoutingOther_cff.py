import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *

def addScouting(process):

   process.photonScoutingTable = cms.EDProducer("SimpleRun3ScoutingPhotonFlatTableProducer",
        src = cms.InputTag("hltScoutingEgammaPacker"),
        cut = cms.string(""),
        name = cms.string("ScoutingPhoton"),
        doc  = cms.string("Photon Scouting Informations"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(False), # this is the main table for the muons
        variables = cms.PSet(
            pt = Var('pt', 'float', precision=14, doc='pt'),
            eta = Var('eta', 'float', precision=14, doc='eta'),
            phi = Var('phi', 'float', precision=14, doc='phi'),
            m = Var('m', 'float', precision=14, doc='m'),
            sigmaIetaIeta = Var('sigmaIetaIeta', 'float', precision=14, doc='sigmaIetaIeta'),
            hOverE = Var('hOverE', 'float', precision=14, doc='hOverE'),
            ecalIso = Var('ecalIso', 'float', precision=14, doc='ecalIso'),
            hcalIso = Var('hcalIso', 'float', precision=14, doc='hcalIso'),
            trkIso = Var('trkIso', 'float', precision=14, doc='trkIso'),
            r9 = Var('r9', 'float', precision=14, doc='r9'),
            sMin = Var('sMin', 'float', precision=14, doc='sMin'),
            sMaj = Var('sMaj', 'float', precision=14, doc='sMaj'),
            seedId = Var('seedId', 'int', doc='seedId'),
        )
   )

   process.electronScoutingTable = cms.EDProducer("SimpleRun3ScoutingElectronFlatTableProducer",
        src = cms.InputTag("hltScoutingEgammaPacker"),
        cut = cms.string(""),
        name = cms.string("ScoutingElectron"),
        doc  = cms.string("Electron Scouting Informations"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(False), # this is the main table for the muons
        variables = cms.PSet(
            pt = Var('pt', 'float', precision=14, doc='pt'),
            eta = Var('eta', 'float', precision=14, doc='eta'),
            phi = Var('phi', 'float', precision=14, doc='phi'),
            m = Var('m', 'float', precision=14, doc='m'),
            d0 = Var('d0', 'float', precision=14, doc='d0'),
            dz = Var('dz', 'float', precision=14, doc='dz'),
            dEtaIn = Var('dEtaIn', 'float', precision=14, doc='dEtaIn'),
            dPhiIn = Var('dPhiIn', 'float', precision=14, doc='dPhiIn'),
            sigmaIetaIeta = Var('sigmaIetaIeta', 'float', precision=14, doc='sigmaIetaIeta'),
            hOverE = Var('hOverE', 'float', precision=14, doc='hOverE'),
            ooEMOOp = Var('ooEMOOp', 'float', precision=14, doc='ooEMOOp'),
            missingHits = Var('missingHits', 'int', doc='missingHits'),
            charge = Var('charge', 'int', doc='charge'),
            ecalIso = Var('ecalIso', 'float', precision=14, doc='ecalIso'),
            hcalIso = Var('hcalIso', 'float', precision=14, doc='hcalIso'),
            trkIso = Var('trkIso', 'float', precision=14, doc='trkIso'),
            r9 = Var('r9', 'float', precision=14, doc='r9'),
            sMin = Var('sMin', 'float', precision=14, doc='sMin'),
            sMaj = Var('sMaj', 'float', precision=14, doc='sMaj'),
            seedId = Var('seedId', 'int', doc='seedId'),
        )
   )

   process.muonScoutingTable = cms.EDProducer("SimpleRun3ScoutingMuonFlatTableProducer",
        src = cms.InputTag("hltScoutingMuonPacker"),
        cut = cms.string(""),
        name = cms.string("ScoutingMuon"),
        doc  = cms.string("Muon Scouting Informations"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(False), # this is the main table for the muons
        variables = cms.PSet(
            pt = Var('pt', 'float', precision=14, doc='pt'),
            eta = Var('eta', 'float', precision=14, doc='eta'),
            phi = Var('phi', 'float', precision=14, doc='phi'),
            m = Var('m', 'float', precision=14, doc='m'),
            type = Var('type', 'int', doc='type'),
            charge = Var('charge', 'int', doc='charge'),
            normchi2 = Var('normalizedChi2', 'float', precision=14, doc='normalizedChi2'),
            ecalIso = Var('ecalIso', 'float', precision=14, doc='ecalIso'),
            hcalIso = Var('hcalIso', 'float', precision=14, doc='hcalIso'),
            nValidStandAloneMuonHits = Var('nValidStandAloneMuonHits', 'int', doc='nValidStandAloneMuonHits'),
            nStandAloneMuonMatchedStations = Var('nStandAloneMuonMatchedStations', 'int', doc='nStandAloneMuonMatchedStations'),
            nValidRecoMuonHits = Var('nValidRecoMuonHits', 'int', doc='nValidRecoMuonHits'),
            nRecoMuonChambers = Var('nRecoMuonChambers', 'int', doc='nRecoMuonChambers'),
            nRecoMuonChambersCSCorDT = Var('nRecoMuonChambersCSCorDT', 'int', doc='nRecoMuonChambersCSCorDT'),
            nRecoMuonMatches = Var('nRecoMuonMatches', 'int', doc='nRecoMuonMatches'),
            nRecoMuonMatchedStations = Var('nRecoMuonMatchedStations', 'int', doc='nRecoMuonMatchedStations'),
            nRecoMuonExpectedMatchedStations = Var('nRecoMuonExpectedMatchedStations', 'int', doc='nRecoMuonExpectedMatchedStations'),
            recoMuonStationMask = Var('recoMuonStationMask', 'int', doc='recoMuonStationMask'),
            nRecoMuonMatchedRPCLayers = Var('nRecoMuonMatchedRPCLayers', 'int', doc='nRecoMuonMatchedRPCLayers'),
            recoMuonRPClayerMask = Var('recoMuonRPClayerMask', 'int', doc='recoMuonRPClayerMask'),
            nValidPixelHits = Var('nValidPixelHits', 'int', doc='nValidPixelHits'),
            nValidStripHits = Var('nValidStripHits', 'int', doc='nValidStripHits'),
            nPixelLayersWithMeasurement = Var('nPixelLayersWithMeasurement', 'int', doc='nPixelLayersWithMeasurement'),
            nTrackerLayersWithMeasurement = Var('nTrackerLayersWithMeasurement', 'int', doc='nTrackerLayersWithMeasurement'),
            trk_chi2 = Var('trk_chi2', 'float', precision=14, doc='trk_chi2'),
            trk_ndof = Var('trk_ndof', 'float', precision=14, doc='trk_ndof'),
            trk_dxy = Var('trk_dxy', 'float', precision=14, doc='trk_dxy'),
            trk_dz = Var('trk_dz', 'float', precision=14, doc='trk_dz'),
            trk_qoverp = Var('trk_qoverp', 'float', precision=14, doc='trk_qoverp'),
            trk_lambda = Var('trk_lambda', 'float', precision=14, doc='trk_lambda'),
            trk_pt = Var('trk_pt', 'float', precision=14, doc='trk_pt'),
            trk_phi = Var('trk_phi', 'float', precision=14, doc='trk_phi'),
            trk_eta = Var('trk_eta', 'float', precision=14, doc='trk_eta'),
            trk_dxyError = Var('trk_dxyError', 'float', precision=14, doc='trk_dxyError'),
            trk_dzError = Var('trk_dzError', 'float', precision=14, doc='trk_dzError'),
            trk_qoverpError = Var('trk_qoverpError', 'float', precision=14, doc='trk_qoverpError'),
            trk_lambdaError = Var('trk_lambdaError', 'float', precision=14, doc='trk_lambdaError'),
            trk_phiError = Var('trk_phiError', 'float', precision=14, doc='trk_phiError'),
            trk_dsz = Var('trk_dsz', 'float', precision=14, doc='trk_dsz'),
            trk_dszError = Var('trk_dszError', 'float', precision=14, doc='trk_dszError'),
            trk_qoverp_lambda_cov = Var('trk_qoverp_lambda_cov', 'float', precision=14, doc='trk_qoverp_lambda_cov'),
            trk_qoverp_phi_cov = Var('trk_qoverp_phi_cov', 'float', precision=14, doc='trk_qoverp_phi_cov'),
            trk_qoverp_dxy_cov = Var('trk_qoverp_dxy_cov', 'float', precision=14, doc='trk_qoverp_dxy_cov'),
            trk_qoverp_dsz_cov = Var('trk_qoverp_dsz_cov', 'float', precision=14, doc='trk_qoverp_dsz_cov'),
            trk_lambda_phi_cov = Var('trk_lambda_phi_cov', 'float', precision=14, doc='trk_lambda_phi_cov'),
            trk_lambda_dxy_cov = Var('trk_lambda_dxy_cov', 'float', precision=14, doc='trk_lambda_dxy_cov'),
            trk_lambda_dsz_cov = Var('trk_lambda_dsz_cov', 'float', precision=14, doc='trk_lambda_dsz_cov'),
            trk_phi_dxy_cov = Var('trk_phi_dxy_cov', 'float', precision=14, doc='trk_phi_dxy_cov'),
            trk_phi_dsz_cov = Var('trk_phi_dsz_cov', 'float', precision=14, doc='trk_phi_dsz_cov'),
            trk_dxy_dsz_cov = Var('trk_dxy_dsz_cov', 'float', precision=14, doc='trk_dxy_dsz_cov'),
            trk_vx = Var('trk_vx', 'float', precision=14, doc='trk_vx'),
            trk_vy = Var('trk_vy', 'float', precision=14, doc='trk_vy'),
            trk_vz = Var('trk_vz', 'float', precision=14, doc='trk_vz'),
        )
   )

   process.trackScoutingTable = cms.EDProducer("SimpleRun3ScoutingTrackFlatTableProducer",
        src = cms.InputTag("hltScoutingTrackPacker"),
        cut = cms.string(""),
        name = cms.string("ScoutingTrack"),
        doc  = cms.string("Track Scouting Informations"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(False), # this is the main table for the muons
        variables = cms.PSet(
            pt = Var('tk_pt', 'float', precision=14, doc='pt'),
            eta = Var('tk_eta', 'float', precision=14, doc='eta'),
            phi = Var('tk_phi', 'float', precision=14, doc='phi'),
            chi2 = Var('tk_chi2', 'float', precision=14, doc='chi2'),
            ndof = Var('tk_ndof', 'float', precision=14, doc='ndof'),
            charge = Var('tk_charge', 'int', doc='charge'),
            dxy = Var('tk_dxy', 'float', precision=14, doc='dxy'),
            dz = Var('tk_dz', 'float', precision=14, doc='dz'),
            nValidPixelHits = Var('tk_nValidPixelHits', 'int', doc='nValidPixelHits'),
            nValidStripHits = Var('tk_nValidStripHits', 'int', doc='nValidStripHits'),
            nPixelLayersWithMeasurement = Var('nPixelLayersWithMeasurement', 'int', doc='nPixelLayersWithMeasurement'),
            nTrackerLayersWithMeasurement = Var('nTrackerLayersWithMeasurement', 'int', doc='nTrackerLayersWithMeasurement'),
            qoverp = Var('tk_qoverp', 'float', precision=14, doc='qoverp'),
            lambda = Var('tk_lambda', 'float', precision=14, doc='lambda'),
            dxyError = Var('tk_dxy_Error', 'float', precision=14, doc='dxyError'),
            dzError = Var('tk_dz_Error', 'float', precision=14, doc='dzError'),
            qoverpError = Var('tk_qoverp_Error', 'float', precision=14, doc='qoverpError'),
            lambdaError = Var('tk_lambda_Error', 'float', precision=14, doc='lambdaError'),
            phiError = Var('tk_phi_Error', 'float', precision=14, doc='phiError'),
            dsz = Var('tk_dsz', 'float', precision=14, doc='dsz'),
            dszError = Var('tk_dsz_Error', 'float', precision=14, doc='dszError'),
            qoverp_lambda_cov = Var('tk_qoverp_lambda_cov', 'float', precision=14, doc='qoverp_lambda_cov'),
            qoverp_phi_cov = Var('tk_qoverp_phi_cov', 'float', precision=14, doc='qoverp_phi_cov'),
            qoverp_dxy_cov = Var('tk_qoverp_dxy_cov', 'float', precision=14, doc='qoverp_dxy_cov'),
            qoverp_dsz_cov = Var('tk_qoverp_dsz_cov', 'float', precision=14, doc='qoverp_dsz_cov'),
            lambda_phi_cov = Var('tk_lambda_phi_cov', 'float', precision=14, doc='lambda_phi_cov'),
            lambda_dxy_cov = Var('tk_lambda_dxy_cov', 'float', precision=14, doc='lambda_dxy_cov'),
            lambda_dsz_cov = Var('tk_lambda_dsz_cov', 'float', precision=14, doc='lambda_dsz_cov'),
            phi_dxy_cov = Var('tk_phi_dxy_cov', 'float', precision=14, doc='phi_dxy_cov'),
            phi_dsz_cov = Var('tk_phi_dsz_cov', 'float', precision=14, doc='phi_dsz_cov'),
            dxy_dsz_cov = Var('tk_dxy_dsz_cov', 'float', precision=14, doc='dxy_dsz_cov'),
            vtxInd = Var('tk_vtxInd', 'int', doc='vtxInd'),
            vx = Var('tk_vx', 'float', precision=14, doc='vx'),
            vy = Var('tk_vy', 'float', precision=14, doc='vy'),
            vz = Var('tk_vz', 'float', precision=14, doc='vz'),
        )
   )

   process.primaryvertexScoutingTable = cms.EDProducer("SimpleRun3ScoutingVertexFlatTableProducer",
        src = cms.Input("hltScoutingPrimaryVertexPacker", "primaryVtx"),
        cut = cms.string(""),
        name = cms.string("ScoutingPrimaryVertex"),
        doc  = cms.string("PrimaryVertex Scouting Informations"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(False), # this is the main table for the muons
        variables = cms.PSet(
            x = Var('x', 'float', precision=14, doc='x'),
            y = Var('y', 'float', precision=14, doc='y'),
            z = Var('z', 'float', precision=14, doc='z'),
            xError = Var('xError', 'float', precision=14, doc='xError'),
            yError = Var('yError', 'float', precision=14, doc='yError'),
            zError = Var('zError', 'float', precision=14, doc='zError'),
            tracksSize = Var('tracksSize', 'int', doc='tracksSize'),
            chi2 = Var('chi2', 'float', precision=14, doc='chi2'),
            ndof = Var('ndof', 'int', doc='ndof'),
            isValidVtx = Var('isValidVtx', 'bool', doc='isValidVtx'),
        )
   )

   process.displacedvertexScoutingTable = cms.EDProducer("SimpleRun3ScoutingVertexFlatTableProducer",
        src = cms.Input("hltScoutingMuonPacker","displacedVtx"),
        cut = cms.string(""),
        name = cms.string("ScoutingDisplacedVertex"),
        doc  = cms.string("DisplacedVertex Scouting Informations"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(False), # this is the main table for the muons
        variables = cms.PSet(
            x = Var('x', 'float', precision=14, doc='x'),
            y = Var('y', 'float', precision=14, doc='y'),
            z = Var('z', 'float', precision=14, doc='z'),
            xError = Var('xError', 'float', precision=14, doc='xError'),
            yError = Var('yError', 'float', precision=14, doc='yError'),
            zError = Var('zError', 'float', precision=14, doc='zError'),
            tracksSize = Var('tracksSize', 'int', doc='tracksSize'),
            chi2 = Var('chi2', 'float', precision=14, doc='chi2'),
            ndof = Var('ndof', 'int', doc='ndof'),
            isValidVtx = Var('isValidVtx', 'bool', doc='isValidVtx'),
        )
   )

   process.run3ScoutingTask = cms.Task(
       process.photonScoutingTable,
       process.electronScoutingTable,
       process.muonScoutingTable,
   )
   process.schedule.associate(process.run3ScoutingTask)

#def addScouting(process):
#
#   process.run3ScoutingTable = cms.EDProducer("Run3ScoutingTableProducer",
#       primaryvertex = cms.InputTag("hltScoutingPrimaryVertexPacker", "primaryVtx"),
#       displacedvertex = cms.InputTag("hltScoutingMuonPacker","displacedVtx"),
#       photon = cms.InputTag("hltScoutingEgammaPacker"),
#       muon = cms.InputTag("hltScoutingMuonPacker"),
#       electron = cms.InputTag("hltScoutingEgammaPacker"),
#       track = cms.InputTag("hltScoutingTrackPacker"),
#       metpt = cms.InputTag("hltScoutingPFPacker","pfMetPt"),
#       metphi = cms.InputTag("hltScoutingPFPacker","pfMetPhi"),
#       rho = cms.InputTag("hltScoutingPFPacker", "rho"),
#   )
#
#   process.run3ScoutingTask = cms.Task(process.run3ScoutingTable)
#   process.schedule.associate(process.run3ScoutingTask)
