o2physics_add_dpl_workflow(my-mch-mid-resolution
                           SOURCES MCHMIDResolution.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework
                                                 O2Physics::AnalysisCore
                                                 O2Physics::PWGDQCore
                                                 O2::DetectorsBase
                                                 O2::CCDB
                           COMPONENT_NAME Analysis)