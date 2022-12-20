# If larcv2 is found via UPS product, LARCV2_* will be defined.
# If larcv2 is built locally, more likely that LARCV_* will be defined.
# Prioritize the local build over UPS product.

if(DEFINED ENV{LARCV_LIBDIR})
  set( LARCV_LIB -L$ENV{LARCV_LIBDIR} -llarcv )
else()
  set( LARCV_LIB -L$ENV{LARCV2_LIBDIR} -llarcv )
endif()

set( ICARUS_LIB -L$ENV{ICARUS_LIB} )

if(DEFINED ENV{LARCV_INCDIR})
  include_directories ( $ENV{LARCV_INCDIR} )
else()
  include_directories ( $ENV{LARCV2_INCDIR} )
endif()

#add_subdirectory(job)

art_make_library(
    LIBRARY_NAME
        icarus_LArCVImageMaker
            LIBRARIES
		sbnobj::Common_CRT
		nurandom::RandomUtils_NuRandomService_service
		larcore::Geometry_Geometry_service
		lardataalg::DetectorInfo
		larcorealg::Geometry
		larcoreobj::SummaryData
		lardataobj::Simulation
		lardataobj::MCBase
		lardataobj::RawData
		lardataobj::RecoBase
		lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
		larevt::SpaceCharge
		nusimdata::SimulationBase
		art::Framework_Services_Registry
		art::Persistency_Common
		art::Persistency_Provenance
		art::Utilities
		canvas::canvas
		cetlib::cetlib
		cetlib_except::cetlib_except
		ROOT::Core ROOT::Tree
		CLHEP::CLHEP
		art::Framework_Core
		art::Framework_Principal
		Boost::system
		${LARCV_LIB}
)

cet_build_plugin(LArSoftSuperaDriver art::module BASENAME_ONLY
            LIBRARIES
		icarus_LArCVImageMaker
		nurandom::RandomUtils_NuRandomService_service
		larcore::Geometry_Geometry_service
		larcorealg::Geometry
		lardata::DetectorInfoServices_DetectorPropertiesServiceStandard_service
		larevt::SpaceCharge
		lardataalg::DetectorInfo
		larcoreobj::SummaryData
		lardataobj::Simulation
		lardataobj::MCBase
		lardataobj::RawData
		lardataobj::RecoBase
		nusimdata::SimulationBase
		art::Framework_Core
		art::Framework_Principal
		art::Framework_Services_Registry
		art::Framework_Services_Optional_RandomNumberGenerator_service
		art::Persistency_Common
		art::Persistency_Provenance
		art::Utilities
		canvas::canvas
		messagefacility::MF_MessageLogger
		fhiclcpp::fhiclcpp
		cetlib::cetlib
		cetlib_except::cetlib_except
		ROOT::Core ROOT::Tree
		Boost::system
		CLHEP::CLHEP
		${LARCV_LIB}
		)

cet_build_plugin(RSEFilter art::module BASENAME_ONLY
            LIBRARIES
		icarus_LArCVImageMaker
		nurandom::RandomUtils_NuRandomService_service
		larcore::Geometry_Geometry_service
		larcorealg::Geometry
		lardata::DetectorInfoServices_DetectorPropertiesServiceStandard_service
		larevt::SpaceCharge
		lardataalg::DetectorInfo
		larcoreobj::SummaryData
		lardataobj::Simulation
		lardataobj::MCBase
		lardataobj::RawData
		lardataobj::RecoBase
		nusimdata::SimulationBase
		art::Framework_Core
		art::Framework_Principal
		art::Framework_Services_Registry
		art::Framework_Services_Optional_RandomNumberGenerator_service
		art::Persistency_Common
		art::Persistency_Provenance
		art::Utilities
		canvas::canvas
		messagefacility::MF_MessageLogger
		fhiclcpp::fhiclcpp
		cetlib::cetlib
		cetlib_except::cetlib_except
		ROOT::Core ROOT::Tree
		Boost::system
		CLHEP::CLHEP
		${LARCV_LIB}
		)

install_headers()
install_source()
add_subdirectory(job)
#add_subdirectory(MultiPartGenerators)
