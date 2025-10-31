# From Tutorial PWGMM cmakelists
# Changed to my muon tasks

o2physics_add_dpl_workflow(my-mc-allgen
                  SOURCES my_mc_allgen.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(my-muon-task
                  SOURCES my_muonPt.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(my-muon-eta-task-types
                  SOURCES my_muon_etaPtbytrack.cxx
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
o2physics_add_dpl_workflow(my-mc-reco
                  SOURCES my_muon_recoMc.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)
                  