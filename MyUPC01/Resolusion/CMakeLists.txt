o2physics_add_dpl_workflow(my-upc-resolution
                           SOURCES MCHMIDResolution.cxx
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)