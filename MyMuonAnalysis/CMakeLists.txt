# From Tutorial PWGMM cmakelists
# Changed to my muon tasks
o2physics_add_dpl_workflow(my-centrality-check
                  SOURCES centrality.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(my-collision-associate-central-fwd
                  SOURCES collision_associate_central_fwd.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(my-collision-associate
                  SOURCES collision_associate.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(my-muon-fake-comparater
                  SOURCES fake_muon_comparater.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(my-mc-fwd-dca
                  SOURCES fwd_dca.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(my-mc-fwd-allgen
                  SOURCES fwd_particles_allgen.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(my-mc-allgen
                  SOURCES my_mc_allgen.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(my-muon-task
                  SOURCES my_muonPt.cxx MuonTrackCut.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(my-muon-eta-phi-task-types
                  SOURCES my_muon_etaPtPhibytrack.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(my-muon-mcgen
                  SOURCES my_muon_mcgen.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)
#o2physics_add_dpl_workflow(my-muon-cut
#                  SOURCES my_muon_o2cut.cxx
#                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::PWGDQCore #O2::Framework O2::DCAFitter O2::GlobalTracking O2Physics::AnalysisCore KFParticle::KFParticle O2Physics::MLCore
#                  COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(my-mc-reco
                  SOURCES my_muon_recoMc.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(my-muon-tracking-efficiency
                  SOURCES my_muon_tracking_efficiency.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(my-mc-res-track
                  SOURCES my_muon_resto_track.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore 
                  COMPONENT_NAME Analysis)

#o2physics_add_dpl_workflow(my-mch-mft-residuals
#                  SOURCES MchMftResiduals.cxx
#                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore 
#                  COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(my-mch-mft-residuals
                  SOURCES MchMftResiduals_copy.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework
                                        O2Physics::AnalysisCore 
                                        O2Physics::PWGDQCore 
                                        O2::ReconstructionDataFormats
                                        O2::DetectorsBase
                                        O2Physics::AnalysisCCDB
                  COMPONENT_NAME Analysis)
# O2Physics/MyMuonAnalysis/fwdtrackextension.cxx
o2physics_add_dpl_workflow(my-mch-mft-fwdtrackextension
                    SOURCES fwdtrackextension.cxx
                    PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                                          O2::ReconstructionDataFormats
                                          O2::DetectorsBase
                                          O2::DetectorsCommonDataFormats
                    COMPONENT_NAME Analysis)
