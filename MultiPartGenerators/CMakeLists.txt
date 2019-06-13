# Add the subfolder containing the deconvolution tools
art_make( BASENAME_ONLY
          MODULE_LIBRARIES larcorealg_Geometry
            larcore_Geometry_Geometry_service
            lardata_Utilities
            larevt_Filters
            lardataobj_RawData
            larevt_CalibrationDBI_IOVData
            larevt_CalibrationDBI_Providers
            lardataobj_RecoBase
            lardata_ArtDataHelper
            larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
            nusimdata_SimulationBase
            lardata_Utilities
            lardataobj_MCBase
            larcoreobj_SummaryData
            nusimdata_SimulationBase
            nutools_RandomUtils_NuRandomService_service
            lardataobj_MCBase
	    #SignalShapingServiceDUNE_service
            ${ART_FRAMEWORK_CORE}
            ${ART_FRAMEWORK_PRINCIPAL}
            ${ART_FRAMEWORK_SERVICES_REGISTRY}
            ${ART_FRAMEWORK_SERVICES_OPTIONAL}
            ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
            art_Persistency_Common
            art_Persistency_Provenance
            art_Utilities
            canvas
            ${MF_MESSAGELOGGER}
            ${MF_UTILITIES}
            ${FHICLCPP}
            ${CETLIB}
            ${ROOT_FFTW}
            ${ROOT_BASIC_LIB_LIST}
	    ${ROOT_EG}
	    #lardata_Utilities_LArFFT_service
            #larcorealg_Geometry
            #larevt_CalibrationDBI_IOVData
            #larevt_CalibrationDBI_Providers
            #larreco_HitFinder
            #lardataobj_RecoBase
            #larcore_Geometry_Geometry_service
            #lardata_Utilities
            #nutools_RandomUtils_NuRandomService_service
            ${ART_FRAMEWORK_SERVICES_BASIC}
            ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
            ${CLHEP}
        )

install_headers()
install_fhicl()
install_source()

