# From Tutorial PWGMM cmakelists
# Changed to my muon tasks


o2physics_add_dpl_workflow(my-muon-task
                  SOURCES my_muonPt.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(my-muon-task-types
                  SOURCES my_muon_etaPtbytrack.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                  COMPONENT_NAME Analysis)
                  