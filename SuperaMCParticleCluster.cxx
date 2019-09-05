#ifndef __SUPERAMCPARTICLECLUSTER_CXX__
#define __SUPERAMCPARTICLECLUSTER_CXX__

#include "SuperaMCParticleCluster.h"
#include "SuperaMCParticleClusterData.h"
#include "larcv/core/DataFormat/EventParticle.h"
#include "larcv/core/DataFormat/EventVoxel3D.h"

namespace larcv {
  
  static SuperaMCParticleClusterProcessFactory __global_SuperaMCParticleClusterProcessFactory__;
  
  SuperaMCParticleCluster::SuperaMCParticleCluster(const std::string name)
    : SuperaBase(name)
  {}
  
  void SuperaMCParticleCluster::configure(const PSet& cfg)
  {
    SuperaBase::configure(cfg);
    _output_label = cfg.get<std::string>("OutputLabel");
    _ref_meta3d_cluster3d = cfg.get<std::string>("Meta3DFromCluster3D","mcst");
    _ref_meta3d_tensor3d = cfg.get<std::string>("Meta3DFromTensor3D","");
    _delta_size = cfg.get<size_t>("DeltaSize",5);
    _eioni_size = cfg.get<size_t>("IonizationSize",5);
    _compton_size = cfg.get<size_t>("ComptonSize",10);
    _compton_energy = cfg.get<double>("ComptonEnergy",4);
    _edep_threshold = cfg.get<double>("EnergyDepositThreshold",0.05);
    _projection_id = cfg.get<int>("ProjectionID",0);
    _use_sed = cfg.get<bool>("UseSimEnergyDeposit");
    _use_true_pos = cfg.get<bool>("UseTruePosition",true);
    
    auto tpc_v = cfg.get<std::vector<unsigned short> >("TPCList");
    larcv::Point3D min_pt(1.e9,1.e9,1.e9);
    larcv::Point3D max_pt(-1.e9,-1.e9,-1.e9);
    for(auto const& tpc_id : tpc_v) {
      auto geop = lar::providerFrom<geo::Geometry>();
      for(size_t c=0; c<geop->Ncryostats(); ++c) {
	auto const& cryostat = geop->Cryostat(c);
	if(!cryostat.HasTPC(tpc_id)) continue;
	auto const& tpcabox = cryostat.TPC(tpc_id).ActiveBoundingBox();
	if(min_pt.x > tpcabox.MinX()) min_pt.x = tpcabox.MinX();
	if(min_pt.y > tpcabox.MinY()) min_pt.y = tpcabox.MinY();
	if(min_pt.z > tpcabox.MinZ()) min_pt.z = tpcabox.MinZ();
	if(max_pt.x < tpcabox.MaxX()) max_pt.x = tpcabox.MaxX();
	if(max_pt.y < tpcabox.MaxY()) max_pt.y = tpcabox.MaxY();
	if(max_pt.z < tpcabox.MaxZ()) max_pt.z = tpcabox.MaxZ();
	break;
      }
    }
    _world_bounds.update(min_pt,max_pt);
  }
  
  
  void SuperaMCParticleCluster::initialize()
  {
    SuperaBase::initialize();
  }
  
  void SuperaMCParticleCluster::FillParticleGroups(const std::vector<simb::MCParticle>& larmcp_v,
						   std::vector<supera::ParticleGroup>& shower_grp_v,
						   std::vector<supera::ParticleGroup>& track_grp_v)
    
  {
    const larcv::Particle invalid_part;
    auto const& parent_pdg_v = _mcpl.ParentPdgCode();
    auto const& trackid2index = _mcpl.TrackIdToIndex();
    for(size_t index=0; index<larmcp_v.size(); ++index) {

      auto const& mcpart = larmcp_v[index];
      int pdg_code       = abs(mcpart.PdgCode());
      int mother_index   = -1;
      int track_id       = mcpart.TrackId();
      if(mcpart.Mother() < ((int)(trackid2index.size())))
	mother_index = trackid2index[mcpart.Mother()];

      //if(pdg_code != -11 && pdg_code != 11 && pdg_code != 22) continue;
      if(pdg_code > 1000000) continue;

      supera::ParticleGroup grp;
      grp.part = this->MakeParticle(mcpart);
      if(mother_index >= 0)
	grp.part.parent_pdg_code(parent_pdg_v[index]);
      grp.valid=true;

      if(pdg_code == 22 || pdg_code == 11) {
	if(pdg_code == 22) {
	  // photon ... reset first, last, and end position
	  grp.type = supera::kPhoton;
	  grp.part.first_step(invalid_part.first_step());
	  grp.part.last_step(invalid_part.last_step());
	  grp.part.end_position(invalid_part.end_position());
	}
	else if(pdg_code == 11) {
	  
	  std::string prc = mcpart.Process();
	  if( prc == "muIoni" || prc == "hIoni")
	    grp.type = supera::kDelta;
	  else if( prc == "muMinusCaptureAtRest" || prc == "muPlusCaptureAtRest" || prc == "Decay" )
	    grp.type = supera::kDecay;
	  else if( prc == "compt" )
	    grp.type = supera::kCompton;
	  else if( prc == "phot")
	    grp.type = supera::kPhotoElectron;
	  else if( prc == "eIoni" )
	    grp.type = supera::kIonization;
	  else if( prc == "conv" )
	    grp.type = supera::kConversion;
	  else if( prc == "primary")
	    grp.type = supera::kPrimary;
	}
	shower_grp_v[track_id] = grp;
      }
      else { 
	grp.type = supera::kTrack;
	if(grp.part.pdg_code() == 2112)
	  grp.type = supera::kNeutron;
	track_grp_v[track_id] = grp;
      }
    }
  }

  void SuperaMCParticleCluster::AnalyzeSimEnergyDeposit(const larcv::Voxel3DMeta& meta,
							std::vector<supera::ParticleGroup>& shower_grp_v,
							std::vector<supera::ParticleGroup>& track_grp_v)
  {
    auto const& larmcp_v = LArData<supera::LArMCParticle_t>();
    auto const& sedep_v = LArData<supera::LArSimEnergyDeposit_t>();
    auto const& trackid2index = _mcpl.TrackIdToIndex();
    LARCV_INFO() << "Processing SimEnergyDeposit array: " << sedep_v.size() << std::endl;
    size_t bad_sedep_counter = 0;
    std::set<size_t> missing_trackid;
    for(size_t sedep_idx=0; sedep_idx<sedep_v.size(); ++sedep_idx) {
      auto const& sedep = sedep_v.at(sedep_idx);

      VoxelID_t vox_id = meta.id(sedep.X(), sedep.Y(), sedep.Z());
      if(vox_id == larcv::kINVALID_VOXELID || !_world_bounds.contains(sedep.X(),sedep.Y(),sedep.Z())) {
	LARCV_DEBUG() << "Skipping sedep from track id " << sedep.TrackID() 
		      << " E=" << sedep.Energy()
		      << " pos=(" << sedep.X() << "," << sedep.Y() << "," << sedep.Z() << ")" << std::endl;
	continue;
      }

      LARCV_DEBUG() << "Recording sedep from track id " << sedep.TrackID() 
		    << " E=" << sedep.Energy() << std::endl;

      supera::EDep pt;
      pt.x = sedep.X();
      pt.y = sedep.Y();
      pt.z = sedep.Z();
      pt.t = sedep.T();
      pt.e = sedep.Energy();

      int track_id = abs(sedep.TrackID());
      if(track_id >= ((int)(trackid2index.size()))) {
	bad_sedep_counter++;
	missing_trackid.insert(track_id);
	continue;
      }
      
      auto const& mcpart = larmcp_v[trackid2index[track_id]];
      int pdg_code = abs(mcpart.PdgCode());
      //if(pdg_code == 2112) std::cout << "Neutron: " << sedep.Energy() << std::endl;

      if(pdg_code == 22 || pdg_code == 11) {
	auto& grp = shower_grp_v[track_id];
	if(!grp.valid) continue;
	grp.vs.emplace(vox_id,pt.e,true);
	//grp.vs.emplace(vox_id,sedep.Energy(),true);
	grp.AddEDep(pt);
      }else{
	auto& grp = track_grp_v[track_id];
	if(!grp.valid) continue;
	grp.vs.emplace(vox_id,pt.e,true);
	grp.AddEDep(pt);
      }
    }
    if(bad_sedep_counter) {
      LARCV_WARNING() << bad_sedep_counter << " SimEnergyDeposit "
		      << "(" << missing_trackid.size() << " track IDs) did not find corresponding MCParticle!"
		      << std::endl;
    }
  }


  void SuperaMCParticleCluster::AnalyzeSimChannel(const larcv::Voxel3DMeta& meta,
						  std::vector<supera::ParticleGroup>& shower_grp_v,
						  std::vector<supera::ParticleGroup>& track_grp_v)
  {
    auto const& sch_v = LArData<supera::LArSimCh_t>();
    auto const& larmcp_v = LArData<supera::LArMCParticle_t>();
    auto const& trackid2index = _mcpl.TrackIdToIndex();
    LARCV_INFO() << "Processing SimChannel array: " << sch_v.size() << std::endl;
    size_t bad_sch_counter = 0;
    std::set<size_t> missing_trackid;

    for(auto const& sch : sch_v) {

      auto ch = sch.Channel();
      // Only use specified projection id, if specified
      if(_projection_id != ::supera::ChannelToProjectionID(ch)) {
	//std::cout << "Skipping channel " << ch << " projection " << ::supera::ChannelToProjectionID(ch) << std::endl;
        continue;
      }
      else{
	//std::cout << "Analyzing channel " << ch << " projection " << ::supera::ChannelToProjectionID(ch) << std::endl;
      }
      for (auto const tick_ides : sch.TDCIDEMap()) {
	
	//x_pos = (supera::TPCTDC2Tick(tick_ides.first) + time_offset) * supera::TPCTickPeriod()  * supera::DriftVelocity();
	double x_pos = supera::TPCTDC2Tick(tick_ides.first) * supera::TPCTickPeriod()  * supera::DriftVelocity();
        //std::cout << tick_ides.first << " : " << supera::TPCTDC2Tick(tick_ides.first) << " : " << time_offset << " : " << x_pos << std::flush;
        //std::cout << (supera::TPCTDC2Tick(tick_ides.first) + time_offset) << std::endl
        //<< (supera::TPCTDC2Tick(tick_ides.first) + time_offset) * supera::TPCTickPeriod() << std::endl
        //<< (supera::TPCTDC2Tick(tick_ides.first) + time_offset) * supera::TPCTickPeriod() * supera::DriftVelocity() << std::endl
        //<< std::endl;

        for (auto const& edep : tick_ides.second) {

	  if(_use_true_pos)
	    x_pos = edep.x;

	  supera::EDep pt;
	  pt.x = x_pos;
	  pt.y = edep.y;
	  pt.z = edep.z;
	  pt.e = edep.energy;

          //supera::ApplySCE(x,y,z);
          //std::cout << " ... " << x << std::endl;

	  VoxelID_t vox_id = meta.id(pt.x, pt.y, pt.z);
	  if (vox_id == larcv::kINVALID_VOXELID) {
	    LARCV_DEBUG() << "Skipping IDE from track id " << edep.trackID
			  << " E=" << edep.energy << " Q=" << edep.numElectrons
			  << " pos=(" << pt.x << "," << pt.y << "," << pt.z << ")" << std::endl;
	    continue;
	  }
	  
	  int track_id = abs(edep.trackID);
	  if(track_id >= ((int)(trackid2index.size()))) {
	    bad_sch_counter++;
	    missing_trackid.insert(track_id);
	    continue;
	  }

	  LARCV_DEBUG() << "Recording IDE from track id " << edep.trackID
			<< " E=" << edep.energy << std::endl;
	  auto const& mcpart = larmcp_v[trackid2index[track_id]];
	  int pdg_code = abs(mcpart.PdgCode());
	  //if(pdg_code == 2112) std::cout << "Neutron: " << edep.energy << std::endl;
	  
	  if(pdg_code == 22 || pdg_code == 11) {
	    auto& grp = shower_grp_v[track_id];
	    if(!grp.valid) continue;
	    grp.vs.emplace(vox_id,pt.e,true);
	    //grp.vs.emplace(vox_id,edep.energy,true);
	    grp.AddEDep(pt);
	  }else{
	    auto& grp = track_grp_v[track_id];
	    if(!grp.valid) continue;
	    grp.vs.emplace(vox_id,pt.e,true);
	    grp.AddEDep(pt);
	  }
	}
      }
    }
    if(bad_sch_counter) {
      LARCV_WARNING() << bad_sch_counter << " SimChannel "
		      << "(" << missing_trackid.size() << " track IDs) did not find corresponding MCParticle!"
		      << std::endl;
    }
    size_t ctr=0;
    for(auto const& grp : shower_grp_v) ctr += grp.vs.size();
    for(auto const& grp : track_grp_v ) ctr += grp.vs.size();
    std::cout<<"voxel count " << ctr << std::endl;
  }


  void SuperaMCParticleCluster::MergeShowerIonizations(std::vector<supera::ParticleGroup>& shower_grp_v)
  {

    int merge_ctr = 0;
    int invalid_ctr = 0;
    do {
      merge_ctr = 0;
      for(auto& grp : shower_grp_v) {
	if(!grp.valid) continue;
	if(grp.type != supera::kIonization) continue;
	// merge to a valid "mother"
	bool parent_found = false;
	int parent_index = grp.part.parent_track_id();
	int parent_index_before = grp.part.track_id();
	while(1) {
	  //std::cout<< "Inspecting: " << grp.part.track_id() << " => " << parent_index << std::endl;
	  if(parent_index <0) {
	    LARCV_DEBUG() << "Invalid parent track id " << parent_index
			  << " Could not find a parent for " << grp.part.track_id() << " PDG " << grp.part.pdg_code()
			  << " " << grp.part.creation_process() << " E = " << grp.part.energy_init() 
			  << " (" << grp.part.energy_deposit() << ") MeV" << std::endl;
	    auto const& parent = shower_grp_v[parent_index_before].part;
	    LARCV_DEBUG() << "Previous parent: " << parent.track_id() << " PDG " << parent.pdg_code() 
			  << " " << parent.creation_process()
			  << std::endl;
	    parent_found=false;
	    invalid_ctr++;
	    break;
	    //throw std::exception();
	  }
	  auto const& parent = shower_grp_v[parent_index];
	  parent_found = parent.valid;
	  if(parent_found) break;
	  else{
	    int ancestor_index = parent.part.parent_track_id();
	    if(ancestor_index == parent_index) {
	      LARCV_INFO() << "Particle " << parent_index << " is root and invalid particle..." << std::endl
			   << "PDG " << parent.part.pdg_code() << " " << parent.part.creation_process() << std::endl;
	      break;
	    }
	    parent_index_before = parent_index;
	    parent_index = ancestor_index;
	  }
	}
	// if parent is found, merge
	if(parent_found) {
	  auto& parent = shower_grp_v[parent_index];
	  for(auto const& vox : grp.vs.as_vector())
	    parent.vs.emplace(vox.id(),vox.value(),true);
	  for(auto const& pt : grp.start.pts)
	    parent.AddEDep(pt);
	  parent.AddEDep(grp.last_pt);
	  parent.trackid_v.push_back(grp.part.track_id());
	  for(auto const& trackid : grp.trackid_v)
	    parent.trackid_v.push_back(trackid);
	  grp.vs.clear_data();
	  grp.valid=false;
	  merge_ctr++;
	  /*
	  std::cout<<"Track " << grp.part.track_id() << " PDG " << grp.part.pdg_code()
		  <<" ... parent found " << parent.part.track_id() << " PDG " << parent.part.pdg_code() << std::endl;
	  */
	}
      }
      LARCV_INFO() << "Ionization merge counter: " << merge_ctr << " invalid counter: " << invalid_ctr << std::endl;
    }while(merge_ctr>0);
  }

  void SuperaMCParticleCluster::ApplyEnergyThreshold(std::vector<supera::ParticleGroup>& shower_grp_v,
						     std::vector<supera::ParticleGroup>& track_grp_v)
  {
    // Loop again and eliminate voxels that has energy below threshold
    for(auto& grp : shower_grp_v) {
      larcv::VoxelSet vs;
      vs.reserve(grp.vs.size());
      for(auto const& vox : grp.vs.as_vector()) {
	if(vox.value() < _edep_threshold) continue;
	vs.emplace(vox.id(),vox.value(),true);
      }
      grp.vs = vs;
      // If compton, here decide whether it should be supera::kComptonHE (high energy)
      if(grp.type == supera::kCompton && grp.vs.size() > _compton_size && grp.vs.sum() > _compton_energy)
        grp.type = supera::kComptonHE;
    }

    for(auto& grp : track_grp_v) {
      larcv::VoxelSet vs;
      vs.reserve(grp.vs.size());
      for(auto const& vox : grp.vs.as_vector()) {
	if(vox.value() < _edep_threshold) continue;
	vs.emplace(vox.id(),vox.value(),true);
      }
      grp.vs = vs;
    }
  }


  void SuperaMCParticleCluster::MergeShowerConversion(std::vector<supera::ParticleGroup>& shower_grp_v)
  {
    int merge_ctr = 0;
    int invalid_ctr = 0;
    do {
      merge_ctr = 0;
      for(auto& grp : shower_grp_v) {
	if(!grp.valid) continue;
	//if(grp.type != supera::kIonization && grp.type != supera::kConversion && grp.type != supera::kComptonHE) continue;
	if(grp.type != supera::kIonization && grp.type != supera::kConversion) continue;
	// merge to a valid "mother"
	bool parent_found = false;
	int parent_index = grp.part.parent_track_id();
	int parent_index_before = grp.part.track_id();
	while(1) {
	  //std::cout<< "Inspecting: " << grp.part.track_id() << " => " << parent_index << std::endl;
	  if(parent_index <0) {
	    LARCV_DEBUG() << "Invalid parent track id " << parent_index
			  << " Could not find a parent for " << grp.part.track_id() << " PDG " << grp.part.pdg_code()
			  << " " << grp.part.creation_process() << " E = " << grp.part.energy_init() 
			  << " (" << grp.part.energy_deposit() << ") MeV" << std::endl;
	    auto const& parent = shower_grp_v[parent_index_before].part;
	    LARCV_DEBUG() << "Previous parent: " << parent.track_id() << " PDG " << parent.pdg_code() 
			  << " " << parent.creation_process()
			  << std::endl;
	    parent_found=false;
	    invalid_ctr++;
	    break;
	    //throw std::exception();
	  }
	  auto const& parent = shower_grp_v[parent_index];
	  parent_found = parent.valid;
	  if(parent_found) break;
	  else{
	    int ancestor_index = parent.part.parent_track_id();
	    if(ancestor_index == parent_index) {
	      LARCV_INFO() << "Particle " << parent_index << " is root and invalid particle..." << std::endl
			   << "PDG " << parent.part.pdg_code() << " " << parent.part.creation_process() << std::endl;
	      break;
	    }
	    parent_index_before = parent_index;
	    parent_index = ancestor_index;
	  }
	}
	// if parent is found, merge
	if(parent_found) {
	  auto& parent = shower_grp_v[parent_index];
	  for(auto const& vox : grp.vs.as_vector())
	    parent.vs.emplace(vox.id(),vox.value(),true);
	  for(auto const& pt : grp.start.pts)
	    parent.AddEDep(pt);
	  parent.AddEDep(grp.last_pt);
	  parent.trackid_v.push_back(grp.part.track_id());
	  for(auto const& trackid : grp.trackid_v)
	    parent.trackid_v.push_back(trackid);
	  grp.vs.clear_data();
	  grp.valid=false;
	  merge_ctr++;
	  /*
	  std::cout<<"Track " << grp.part.track_id() << " PDG " << grp.part.pdg_code() << " " << grp.part.creation_process()
		   <<" ... parent found " << parent.part.track_id() << " PDG " << parent.part.pdg_code() << " " << parent.part.creation_process()
		   << std::endl;
	  */
	}
      }
      LARCV_INFO() << "Merge counter: " << merge_ctr << " invalid counter: " << invalid_ctr << std::endl;
    }while(merge_ctr>0);
  }


  void SuperaMCParticleCluster::MergeShowerTouching(const larcv::Voxel3DMeta& meta,
						    std::vector<supera::ParticleGroup>& shower_grp_v)
  {
    int merge_ctr = 0;
    int invalid_ctr = 0;
    do {
      merge_ctr = 0;
      for(auto& grp : shower_grp_v) {
	if(!grp.valid) continue;
	//if(grp.type != supera::kIonization && grp.type != supera::kConversion && grp.type != supera::kComptonHE) continue;
	if(grp.type == supera::kTrack || grp.type == supera::kDelta) continue;
	// search for a possible parent
	int parent_trackid = -1;
	// a direct parent ?
	if(shower_grp_v[grp.part.parent_track_id()].valid)
	  parent_trackid = grp.part.parent_track_id();
	else {
	  for(size_t shower_trackid = 0; shower_trackid<shower_grp_v.size(); ++shower_trackid) {
	    auto const& candidate_grp = shower_grp_v[shower_trackid];
	    if(shower_trackid == grp.part.parent_track_id() || !candidate_grp.valid) continue;
	    for(auto const& trackid : candidate_grp.trackid_v) {
	      if(trackid != grp.part.parent_track_id()) continue;
	      parent_trackid = shower_trackid;
	      break;
	    }
	    if(parent_trackid >= 0) break;
	  }
	}
	if(parent_trackid >= 0 && parent_trackid != ((int)(grp.part.track_id()))) {
	  auto& parent = shower_grp_v[parent_trackid];
	  if(this->IsTouching(meta,grp.vs,parent.vs)) {
	    // merge
	    for(auto const& vox : grp.vs.as_vector())
	      parent.vs.emplace(vox.id(),vox.value(),true);
	    for(auto const& pt : grp.start.pts)
	      parent.AddEDep(pt);
	    parent.AddEDep(grp.last_pt);
	    parent.trackid_v.push_back(grp.part.track_id());
	    for(auto const& trackid : grp.trackid_v)
	      parent.trackid_v.push_back(trackid);
	    grp.vs.clear_data();
	    grp.valid=false;
	    merge_ctr++;
	    /*
	    std::cout<<"Track " << grp.part.track_id() << " PDG " << grp.part.pdg_code() << " " << grp.part.creation_process()
		     <<" ... parent found " << parent.part.track_id() << " PDG " << parent.part.pdg_code() << " " << parent.part.creation_process()
		     << std::endl;
	    */
	  }
	}
      }
      LARCV_INFO() << "Merge counter: " << merge_ctr << " invalid counter: " << invalid_ctr << std::endl;
    }while(merge_ctr>0);
  }

  void SuperaMCParticleCluster::MergeShowerDeltas(std::vector<supera::ParticleGroup>& shower_grp_v,
						  std::vector<supera::ParticleGroup>& track_grp_v)
  {
    for(auto& grp : shower_grp_v) {
      if(grp.type != supera::kDelta) continue;
      int parent_trackid = grp.part.parent_track_id();
      auto& parent = track_grp_v[parent_trackid];
      if(!parent.valid) continue;
      // if voxel count is smaller than delta ray requirement, simply merge
      if(grp.vs.size() < _delta_size) {
	for(auto const& vox : grp.vs.as_vector())
	  parent.vs.emplace(vox.id(),vox.value(),true);
	for(auto const& pt : grp.start.pts)
	  parent.AddEDep(pt);
	parent.AddEDep(grp.last_pt);
	parent.trackid_v.push_back(grp.part.track_id());
	for(auto const& trackid : grp.trackid_v)
	  parent.trackid_v.push_back(trackid);
	grp.vs.clear_data();
	grp.valid=false;
      }else{
	// check unique number of voxels
	size_t unique_voxel_count = 0;
	for(auto const& vox : grp.vs.as_vector()) {
	  if(parent.vs.find(vox.id()).id() != larcv::kINVALID_VOXELID)
	    continue;
	  ++unique_voxel_count;
	}
	if(unique_voxel_count < _delta_size) {
	  for(auto const& vox : grp.vs.as_vector())
	    parent.vs.emplace(vox.id(),vox.value(),true);
	  for(auto const& pt : grp.start.pts)
	    parent.AddEDep(pt);
	  parent.AddEDep(grp.last_pt);
	  parent.trackid_v.push_back(grp.part.track_id());
	  for(auto const& trackid : grp.trackid_v)
	    parent.trackid_v.push_back(trackid);
	  grp.vs.clear_data();
	  grp.valid=false;
	}
      }
    }
  }

  bool SuperaMCParticleCluster::process(IOManager& mgr)
  {

    SuperaBase::process(mgr);
    larcv::Voxel3DMeta meta;

    if(!_ref_meta3d_cluster3d.empty()) {
      auto const& ev_cluster3d = mgr.get_data<larcv::EventClusterVoxel3D>(_ref_meta3d_cluster3d);
      meta = ev_cluster3d.meta();
    }
    else if(!_ref_meta3d_tensor3d.empty()) {
      auto const& ev_tensor3d = mgr.get_data<larcv::EventSparseTensor3D>(_ref_meta3d_tensor3d);
      meta = ev_tensor3d.meta();
    }
    // Build MCParticle List
    auto const& larmcp_v = LArData<supera::LArMCParticle_t>();
    auto const *ev = GetEvent();
    _mcpl.Update(larmcp_v,ev->id().run(),ev->id().event());

    auto const& trackid2index = _mcpl.TrackIdToIndex();
    // Create ParticleGroup
    std::vector<supera::ParticleGroup> shower_grp_v(trackid2index.size());
    std::vector<supera::ParticleGroup> track_grp_v(trackid2index.size());
    this->FillParticleGroups(larmcp_v, shower_grp_v, track_grp_v);

    // Fill Voxel Information
    if(_use_sed)
      this->AnalyzeSimEnergyDeposit(meta,shower_grp_v, track_grp_v);
    else
      this->AnalyzeSimChannel(meta,shower_grp_v, track_grp_v);

    // Merge fragments of showers
    this->MergeShowerIonizations(shower_grp_v);

    // Apply energy threshold
    this->ApplyEnergyThreshold(shower_grp_v,track_grp_v);

    // Merge fragments of showers
    this->MergeShowerConversion(shower_grp_v);

    // Merge touching shower fragments
    this->MergeShowerTouching(meta,shower_grp_v);

    // merge too small deltas into tracks
    this->MergeShowerDeltas(shower_grp_v,track_grp_v);

    // Assign output IDs
    std::vector<int> trackid2output(trackid2index.size(),-1);
    std::vector<int> output2trackid;
    std::vector<size_t> output_track_trackid, output_shower_trackid;
    output2trackid.reserve(trackid2index.size());
    output_track_trackid.reserve(trackid2index.size());
    output_shower_trackid.reserve(trackid2index.size());
    for(auto& grp : track_grp_v) {
      if(!grp.valid) continue;
      if(grp.vs.size()<1) continue;
      grp.part.energy_deposit(grp.vs.sum());
      if(grp.type == supera::kNeutron) continue;
      size_t output_counter = output2trackid.size();
      grp.part.id(output_counter);
      trackid2output[grp.part.track_id()] = output_counter;
      for(auto const& child : grp.trackid_v)
	trackid2output[child] = output_counter;
      output2trackid.push_back(grp.part.track_id());
      output_track_trackid.push_back(grp.part.track_id());
      ++output_counter;
    }
    for(auto& grp : shower_grp_v) {
      if(!grp.valid) continue;
      if(grp.vs.size()<1) continue;
      grp.part.energy_deposit(grp.vs.sum());
      if(grp.type == supera::kCompton || grp.type == supera::kPhotoElectron) continue;
      size_t output_counter = output2trackid.size();
      grp.part.id(output_counter);
      trackid2output[grp.part.track_id()] = output_counter;
      for(auto const& child : grp.trackid_v)
	trackid2output[child] = output_counter;
      output2trackid.push_back(grp.part.track_id());
      output_shower_trackid.push_back(grp.part.track_id());
    }

    // Assign relationships
    for(auto const& trackid : output_shower_trackid) {
      auto& grp = shower_grp_v[trackid];
      int parent_trackid = grp.part.parent_track_id();
      if(trackid2output[parent_trackid] >= 0) {
      /*
      if(trackid2output[parent_trackid] < 0)
	grp.part.parent_id(grp.part.id());
      else {
      */
	grp.part.parent_id(trackid2output[parent_trackid]);
	int parent_output_id = trackid2output[parent_trackid];
	int parent_id = output2trackid[parent_output_id]; 
	if(shower_grp_v[parent_id].valid)
	  shower_grp_v[parent_id].part.children_id(grp.part.id());
	else if(track_grp_v[parent_id].valid)
	  track_grp_v[parent_id].part.children_id(grp.part.id());
      }
    }

    // At this point, count total number of voxels (will be used for x-check later)
    size_t total_vs_size = 0;
    for(auto const& grp : shower_grp_v) {
      if(!grp.valid) continue;
      total_vs_size += grp.vs.size();
    }
    for(auto const& grp : track_grp_v) {
      if(!grp.valid) continue;
      total_vs_size += grp.vs.size();
    }

    // Combine shower and track groups
    std::vector<supera::ParticleGroup> combined_grp_v(trackid2index.size());
    for(size_t index=0; index < combined_grp_v.size(); ++index) {
      if(track_grp_v[index].valid) std::swap(track_grp_v[index],combined_grp_v[index]);
      else if(shower_grp_v[index].valid) std::swap(shower_grp_v[index],combined_grp_v[index]);
      if(!combined_grp_v[index].valid) continue;
      auto& grp = combined_grp_v[index];
      if(grp.vs.as_vector().empty()) continue;
      auto& part = grp.part;
      auto const& first_pt = grp.start.FirstPoint();
      auto const& last_pt  = grp.last_pt;
      //std::cout<<first_pt.x<< " " << first_pt.y << " " << first_pt.z << std::endl;
      part.first_step(first_pt.x,first_pt.y,first_pt.z,first_pt.t);
      part.last_step(last_pt.x,last_pt.y,last_pt.z,last_pt.t);
    }

    // Unique group id for output particles
    int group_counter = -1;
    
    // loop over MCShower to assign parent/ancestor information
    auto const& mcs_v = LArData<supera::LArMCShower_t>();
    LARCV_INFO() << "Processing MCShower array: " << mcs_v.size() << std::endl;
    for(auto const& mcs : mcs_v) {
      int track_id = mcs.TrackID();
      if(track_id >= ((int)(trackid2output.size()))) {
	LARCV_INFO() << "MCShower " << track_id << " PDG " << mcs.PdgCode() 
		     << " not found in output group..." << std::endl;
	continue;
      }
      
      int output_id = trackid2output[track_id];
      int group_id  = -1;
      if(output_id >= 0) {
	auto& grp = combined_grp_v[track_id];
	if(grp.part.group_id() == larcv::kINVALID_INSTANCEID) {
	  if(group_id < 0) { 
	    group_id = group_counter + 1;
	    ++group_counter;
	  }
	  grp.part.group_id(group_id);
	}
	grp.part.parent_position(mcs.MotherStart().X(),
				 mcs.MotherStart().Y(),
				 mcs.MotherStart().Z(),
				 mcs.MotherStart().T());
	grp.part.parent_creation_process(mcs.MotherProcess());
	grp.part.ancestor_position(mcs.AncestorStart().X(),
				   mcs.AncestorStart().Y(),
				   mcs.AncestorStart().Z(),
				   mcs.AncestorStart().T());
	grp.part.ancestor_track_id(mcs.AncestorTrackID());
	grp.part.ancestor_pdg_code(mcs.AncestorPdgCode());
	grp.part.ancestor_creation_process(mcs.AncestorProcess());
      }

      for(auto const& child : mcs.DaughterTrackID()) {
	//if(child < trackid2output.size() && trackid2output[child] < 0)
	if(child < trackid2output.size() && trackid2output[child] >= 0) {
	  trackid2output[child] = output_id;
	  auto& grp = combined_grp_v[child];
	  if(grp.part.group_id() == larcv::kINVALID_INSTANCEID) {
	    if(group_id < 0) {
	      group_id = group_counter + 1;
	      ++group_counter;
	    }
	    grp.part.group_id(group_id);
	  }
	  grp.part.ancestor_position(mcs.AncestorStart().X(),
				     mcs.AncestorStart().Y(),
				     mcs.AncestorStart().Z(),
				     mcs.AncestorStart().T());
	  grp.part.ancestor_track_id(mcs.AncestorTrackID());
	  grp.part.ancestor_pdg_code(mcs.AncestorPdgCode());
	  grp.part.ancestor_creation_process(mcs.AncestorProcess());
	}
      }
    }

    // loop over MCTrack to assign parent/ancestor information
    auto const& mct_v = LArData<supera::LArMCTrack_t>();
    LARCV_INFO() << "Processing MCTrack array: " << mct_v.size() << std::endl;
    for(auto const& mct : mct_v) {
      int track_id = mct.TrackID();
      int output_id = trackid2output[track_id];
      int group_id = -1;
      if(output_id >= 0) {
	auto& grp = combined_grp_v[track_id];
	if(group_id < 0) {
	  group_id = group_counter + 1;
	  ++group_counter;
	}
	grp.part.group_id(group_id);
	grp.part.parent_position(mct.MotherStart().X(),
				 mct.MotherStart().Y(),
				 mct.MotherStart().Z(),
				 mct.MotherStart().T());
	grp.part.parent_creation_process(mct.MotherProcess());
	grp.part.ancestor_position(mct.AncestorStart().X(),
				   mct.AncestorStart().Y(),
				   mct.AncestorStart().Z(),
				   mct.AncestorStart().T());
	grp.part.ancestor_track_id(mct.AncestorTrackID());
	grp.part.ancestor_pdg_code(mct.AncestorPdgCode());
	grp.part.ancestor_creation_process(mct.AncestorProcess());
      }
      for(size_t output_index=0; output_index<output2trackid.size(); ++output_index) {
	int output_trackid = output2trackid[output_index];
	auto& grp = combined_grp_v[output_trackid];
	if((int)(grp.part.parent_track_id()) != track_id) continue;
	/*
	if(group_id < 0) {
	  group_id = group_counter + 1;
	  ++group_counter;
	}
	grp.part.group_id(group_id);
	*/
	grp.part.ancestor_position(mct.AncestorStart().X(),
				   mct.AncestorStart().Y(),
				   mct.AncestorStart().Z(),
				   mct.AncestorStart().T());
	grp.part.ancestor_track_id(mct.AncestorTrackID());
	grp.part.ancestor_pdg_code(mct.AncestorPdgCode());
	grp.part.ancestor_creation_process(mct.AncestorProcess());
      }
    }

    // for shower particles with invalid parent ID, attempt a search
    for(size_t out_index=0; out_index<output2trackid.size(); ++out_index) {
      int trackid = output2trackid[out_index];
      auto& grp = combined_grp_v[trackid];
      if(grp.part.parent_id() != larcv::kINVALID_INSTANCEID) continue;
      if(grp.type != supera::kComptonHE && grp.type != supera::kPhoton && grp.type != supera::kComptonHE && grp.type != supera::kConversion)
	continue;
      int own_partid = grp.part.id();
      // initiate a search of parent in the valid output particle
      int parent_trackid = grp.part.parent_track_id();
      int parent_partid  = -1;
      while(1) {
	if(parent_trackid >= ((int)(trackid2index.size())) || trackid2index[parent_trackid] <0)
	  break;
	if(parent_trackid < ((int)(trackid2output.size())) && 
	   trackid2output[parent_trackid] >= 0 &&
	   combined_grp_v[parent_trackid].valid ) {
	  //parent_partid = trackid2output[parent_trackid];
	  parent_partid = combined_grp_v[parent_trackid].part.id();
	  break;
	}
	parent_trackid = larmcp_v[trackid2index[parent_trackid]].Mother();
      }
      if(parent_partid >=0) {
	grp.part.parent_id(parent_partid);
	combined_grp_v[parent_trackid].part.children_id(own_partid);
	LARCV_INFO() << "PartID " << own_partid << " (output index " << out_index << ") assigning parent " << parent_partid << std::endl;
      }else{
	grp.part.parent_id(grp.part.id());
	LARCV_INFO() << "PartID " << own_partid << " (output index " << out_index << ") assigning itself as a parent..." << std::endl;
      }
    }


    // Now loop over otuput particle list and check if any remaining group id needs to be assigned
    // Use ancestor to group...
    std::map<size_t,size_t> group_id_by_ancestor;
    for(size_t output_index=0; output_index<output2trackid.size(); ++output_index) {
      auto& grp = combined_grp_v[output2trackid[output_index]];
      if(grp.part.group_id() != larcv::kINVALID_INSTANCEID) continue;
      // If delta, its own grouping
      if(grp.type == supera::kDelta) {
	grp.part.group_id(group_counter + 1);
	for(auto const& child_index : grp.part.children_id())
	  combined_grp_v[output2trackid[child_index]].part.group_id(group_counter+1);
	++group_counter;
      }else{
	if(group_id_by_ancestor.find(grp.part.ancestor_track_id()) == group_id_by_ancestor.end()) {
	  grp.part.group_id(group_counter + 1);
	  group_id_by_ancestor[grp.part.ancestor_track_id()] = group_counter + 1;
	  ++group_counter;
	}else
	  grp.part.group_id(group_id_by_ancestor[grp.part.ancestor_track_id()]);
      }
    }

    // now loop over to find any particle for which first_step is not defined
    for(size_t output_index=0; output_index<output2trackid.size(); ++output_index) {
      auto& grp = combined_grp_v[output2trackid[output_index]];
      auto const& fs = grp.part.first_step();
      if(fs.x()!=0. || fs.y()!=0. || fs.z()!=0.) continue;
      auto const vtx = grp.part.position().as_point3d();
      double min_dist = std::fabs(larcv::kINVALID_DOUBLE);
      larcv::Point3D min_pt;
      for(auto const& vox : grp.vs.as_vector()) {
	auto const pt = meta.position(vox.id());
	double dist = pt.squared_distance(vtx);
	if(dist > min_dist) continue;
	min_dist = dist;
	min_pt = pt;
      }
      if(min_dist > (sqrt(3.) + 1.e-3)) grp.part.first_step(larcv::kINVALID_DOUBLE,larcv::kINVALID_DOUBLE,larcv::kINVALID_DOUBLE,larcv::kINVALID_DOUBLE);
      else grp.part.first_step(min_pt.x, min_pt.y, min_pt.z, grp.part.position().t());

    }

    // now loop over to create VoxelSet for compton/photoelectron
    std::vector<larcv::Particle> part_v; part_v.resize(output2trackid.size());
    auto event_cluster_he = (EventClusterVoxel3D*)(mgr.get_data("cluster3d",_output_label));
    auto event_cluster_le = (EventClusterVoxel3D*)(mgr.get_data("cluster3d",_output_label + "_lowE"));
    auto event_leftover   = (EventSparseTensor3D*)(mgr.get_data("sparse3d",_output_label + "_leftover"));
    event_cluster_he->resize(output2trackid.size());
    event_cluster_le->resize(output2trackid.size());
    for(size_t index=0; index<output2trackid.size(); ++index) {
      int trackid = output2trackid[index];
      auto& grp   = combined_grp_v[trackid];
      std::swap(grp.part, part_v[index]);
      std::swap(grp.vs, event_cluster_he->writeable_voxel_set(index));
      grp.valid=false;
    }

    for(auto& grp : combined_grp_v) {
      if(grp.type != supera::kCompton && grp.type != supera::kPhotoElectron && grp.type != supera::kNeutron) continue;
      if(!grp.valid) continue;
      if(grp.vs.size()<1) continue;
      int trackid = grp.part.parent_track_id();
      int output_index = -1;
      if(trackid < ((int)(trackid2output.size()))) output_index = trackid2output[trackid];
      if(output_index<0) {
	// search the first direct, valid parent
	while(1) {
	  if(trackid >= ((int)(trackid2index.size())) || trackid2index[trackid] < 0) break;
	  trackid = larmcp_v[trackid2index[trackid]].Mother();
	  if(trackid < ((int)(trackid2output.size()))) {
	    output_index = trackid2output[trackid];
	    break;
	  }
	}
      }
      if(output_index<0) continue;

      auto& vs = event_cluster_le->writeable_voxel_set(output_index);
      for(auto const& vox : grp.vs.as_vector()) vs.emplace(vox.id(),vox.value(),true);
      grp.vs.clear_data();
    }

    // create particle ID vs ... overlapped voxel gets higher id number
    auto const& main_vs = event_cluster_he->as_vector();
    auto const& lowe_vs = event_cluster_le->as_vector();
    auto event_cindex = (EventSparseTensor3D*)(mgr.get_data("sparse3d",_output_label + "_index"));
    larcv::VoxelSet cid_vs;
    for(size_t index=0; index<main_vs.size(); ++index) {
      auto const& vs = main_vs[index];
      for(auto const& vox : vs.as_vector()) cid_vs.emplace(vox.id(),(float)(index),false);
    }
    event_cindex->emplace(std::move(cid_vs),meta);

    // Count output voxel count and x-check
    size_t output_vs_size = 0;
    for(auto const& vs : main_vs) output_vs_size += vs.size();
    for(auto const& vs : lowe_vs) output_vs_size += vs.size();
    LARCV_INFO() << "Voxel count x-check: output = " << output_vs_size << " ... total = " << total_vs_size << std::endl;

    LARCV_INFO()<<"Combined reminders..."<<std::endl;
    larcv::VoxelSet leftover_vs; leftover_vs.reserve(total_vs_size - output_vs_size);
    if(total_vs_size > output_vs_size) {
      int ctr= 0;
      for(auto& grp : combined_grp_v) {
	if(grp.vs.as_vector().empty()) continue;
	for(auto const& vox : grp.vs.as_vector()) leftover_vs.emplace(vox.id(),vox.value(),true);
	ctr++;
	auto const& part = grp.part;
	LARCV_INFO() << "Particle ID " << part.id() << " Type " << grp.type << " Valid " << grp.valid << " Track ID " << part.track_id() << " PDG " << part.pdg_code() 
		     << " " << part.creation_process() << " ... " << part.energy_init() << " MeV => " << part.energy_deposit() << " MeV "
		     << grp.trackid_v.size() << " children " << grp.vs.size() << " voxels " << grp.vs.sum() << " MeV" << std::endl;
	LARCV_INFO() << "  Parent " << part.parent_track_id() << " PDG " << part.parent_pdg_code() << " " << part.parent_creation_process() 
		     << " Ancestor " << part.ancestor_track_id() << " PDG " << part.ancestor_pdg_code() << " " << part.ancestor_creation_process() << std::endl;
	LARCV_INFO() << "  Group ID: " << part.group_id() << std::endl;
	std::stringstream ss1,ss2;
	
	ss1 << "  Children particle IDs: " << std::flush;
	for(auto const& child : part.children_id()) ss1 << child << " " << std::flush;
	ss1 << std::endl;
	LARCV_INFO() << ss1.str();

	ss2 << "  Children track IDs: " << std::flush;
	for(auto const& child : grp.trackid_v) ss2 << child << " " << std::flush;
	ss2 << std::endl;
	LARCV_INFO() << ss2.str();

	if(grp.type != supera::kIonization && grp.type != supera::kConversion && grp.type != supera::kComptonHE) continue;
	LARCV_INFO() << "Above was supposed to be merged..." << std::endl;
	
      }
      LARCV_INFO() <<"Shower reminders..."<<std::endl;
      for(auto& grp : shower_grp_v) {
	if(grp.vs.as_vector().empty()) continue;
	for(auto const& vox : grp.vs.as_vector()) leftover_vs.emplace(vox.id(),vox.value(),true);
	ctr++;
	auto const& part = grp.part;
	LARCV_INFO() << "Particle ID " << part.id() << " Track ID " << part.track_id() << " PDG " << part.pdg_code() 
		     << " " << part.creation_process() << " ... " << part.energy_init() << " MeV => " << part.energy_deposit() << " MeV "
		     << grp.trackid_v.size() << " children " << grp.vs.size() << " voxels " << grp.vs.sum() << " MeV" << std::endl;
	LARCV_INFO() << "  Parent " << part.parent_track_id() << " PDG " << part.parent_pdg_code() << " " << part.parent_creation_process() 
		     << " Ancestor " << part.ancestor_track_id() << " PDG " << part.ancestor_pdg_code() << " " << part.ancestor_creation_process() << std::endl;
	LARCV_INFO() << "  Group ID: " << part.group_id() << std::endl;
	std::stringstream ss1,ss2;

	ss1 << "  Children particle IDs: " << std::flush;
	for(auto const& child : part.children_id()) ss1 << child << " " << std::flush;
	ss1 << std::endl;
	LARCV_INFO() << ss1.str();

	ss2 << "  Children track IDs: " << std::flush;
	for(auto const& child : grp.trackid_v) ss2 << child << " " << std::flush;
	ss2 << std::endl;
	LARCV_INFO() << ss2.str();

	if(grp.type != supera::kIonization && grp.type != supera::kConversion && grp.type != supera::kComptonHE) continue;
	LARCV_INFO() << "Above was supposed to be merged..." << std::endl;
	
      }
      LARCV_INFO()<<"Track reminders..."<<std::endl;
      for(auto& grp : track_grp_v) {
	if(grp.vs.as_vector().empty()) continue;
	for(auto const& vox : grp.vs.as_vector()) leftover_vs.emplace(vox.id(),vox.value(),true);
	ctr++;
	auto const& part = grp.part;
	LARCV_INFO() << "Particle ID " << part.id() << " Track ID " << part.track_id() << " PDG " << part.pdg_code() 
		     << " " << part.creation_process() << " ... " << part.energy_init() << " MeV => " << part.energy_deposit() << " MeV "
		     << grp.trackid_v.size() << " children " << grp.vs.size() << " voxels " << grp.vs.sum() << " MeV" << std::endl;
	LARCV_INFO() << "  Parent " << part.parent_track_id() << " PDG " << part.parent_pdg_code() << " " << part.parent_creation_process() 
		     << " Ancestor " << part.ancestor_track_id() << " PDG " << part.ancestor_pdg_code() << " " << part.ancestor_creation_process() << std::endl;
	LARCV_INFO() << "  Group ID: " << part.group_id() << std::endl;
	std::stringstream ss1,ss2;
	
	ss1 << "  Children particle IDs: " << std::flush;
	for(auto const& child : part.children_id()) ss1 << child << " " << std::flush;
	ss1 << std::endl;
	LARCV_INFO() << ss1.str();

	ss2 << "  Children track IDs: " << std::flush;
	for(auto const& child : grp.trackid_v) ss2 << child << " " << std::flush;
	ss2 << std::endl;
	LARCV_INFO() << ss2.str();

	if(grp.type != supera::kIonization && grp.type != supera::kConversion && grp.type != supera::kComptonHE) continue;
	LARCV_INFO() << "Above was supposed to be merged..." << std::endl;
	
      }
      LARCV_INFO() << "... " << ctr << " particles" << std::endl;
    }
    event_leftover->emplace(std::move(leftover_vs),meta);

    // Loop over to find any "stil valid" supera::kIonization supera::kConversion supera::kComptonHE
    LARCV_INFO() << "Particle list" << std::endl;
    for(size_t index = 0; index < part_v.size(); ++index) {
      int trackid = output2trackid[index];
      auto const& grp  = combined_grp_v[trackid];
      auto const& part = part_v[index];
      auto const& vs0  = main_vs[index];
      auto const& vs1  = lowe_vs[index];
      LARCV_INFO() << "Particle ID " << part.id() << " Track ID " << part.track_id() << " PDG " << part.pdg_code() 
		   << " " << part.creation_process() << " ... " << part.energy_init() << " MeV => " << part.energy_deposit() << " MeV "
		   << grp.trackid_v.size() << " children " << vs0.size() << " (" << vs1.size() << ") voxels" << std::endl;
      LARCV_INFO() << "  Parent TrackID " << part.parent_track_id() << " PartID " << part.parent_id() 
		   << " PDG " << part.parent_pdg_code() << " " << part.parent_creation_process() 
		   << " Ancestor TrackID " << part.ancestor_track_id() 
		   << " PDG " << part.ancestor_pdg_code() << " " << part.ancestor_creation_process() << std::endl;
      LARCV_INFO() << "  Group ID: " << part.group_id() << std::endl;
      std::stringstream ss1,ss2;

      ss1 << "  Children particle IDs: " << std::flush;
      for(auto const& child : part.children_id()) ss1 << child << " " << std::flush;
      ss1 << std::endl;
      LARCV_INFO() << ss1.str();

      ss2 << "  Children track IDs: " << std::flush;
      for(auto const& child : grp.trackid_v) ss2 << child << " " << std::flush;
      ss2 << std::endl;
      LARCV_INFO() << ss2.str();

      LARCV_INFO() << "  Start: " << part.first_step().x() << " " << part.first_step().y() << " " << part.first_step().z() << std::endl;
    }
    LARCV_INFO() << "... " << part_v.size() << " particles" << std::endl;

    // Store output
    auto event_mcp = (EventParticle*)(mgr.get_data("particle",_output_label));
    event_mcp->emplace(std::move(part_v));

    // Construct the segmentation
    enum SemanticType_t {
      kSemanticTrack,
      kSemanticShower,
      kSemanticDelta,
      kSemanticMichel,
      kSemanticCompton
    };
    auto event_segment = (EventSparseTensor3D*)(mgr.get_data("sparse3d",_output_label + "_semantics"));
    larcv::VoxelSet semantic_vs; semantic_vs.reserve(total_vs_size);
    for(auto const& vs : event_cluster_le->as_vector()) {
      for(auto const& vox : vs.as_vector()) semantic_vs.emplace(vox.id(),(float)(kSemanticCompton),false);
    }
    for(auto const& vox : event_leftover->as_vector()) semantic_vs.emplace(vox.id(),(float)(kSemanticCompton),false);

    for(size_t index=0; index<output2trackid.size(); ++index) {
      int trackid = output2trackid[index];
      auto const& grp = combined_grp_v[trackid];
      auto const& vs  = event_cluster_he->as_vector()[index];
      float semantic = -1;
      switch(grp.type) {
      case supera::kTrack:
	semantic=(float)kSemanticTrack; break;
      case supera::kPrimary:
      case supera::kPhoton:
      case supera::kComptonHE:
      case supera::kConversion:
	semantic=(float)kSemanticShower; break;
      case supera::kDelta:
	semantic=(float)kSemanticDelta; break;
      case supera::kDecay:
	semantic=(float)kSemanticMichel; break;
	//case supera::kNeutron:
	//semantic=(float)kSemanticCompton; break;
      default:
	LARCV_CRITICAL() << "Unexpected type while assigning semantic class: " << grp.type << std::endl;
	break;
      }
      if(semantic<0) throw std::exception();
      
      for(auto const& vox : vs.as_vector()) {
	auto const& prev = semantic_vs.find(vox.id());
	if(prev.id() == larcv::kINVALID_VOXELID || prev.value() > semantic) 
	  semantic_vs.emplace(vox.id(),semantic,false);
      }
    }
    event_segment->emplace(std::move(semantic_vs),meta);

    return true;
  }

  larcv::Particle 
  SuperaMCParticleCluster::MakeParticle(const supera::LArMCParticle_t& larmcp)
  {
    
    larcv::Particle res;
    res.pdg_code(larmcp.PdgCode());
    res.shape(larcv::kShapeShower);
    res.track_id(larmcp.TrackId());
    res.pdg_code(larmcp.PdgCode());
    res.momentum(larmcp.Px()*1000,larmcp.Py()*1000,larmcp.Pz()*1000);
    res.position(larmcp.Vx(),larmcp.Vy(), larmcp.Vz(), larmcp.T());
    res.first_step(larmcp.Vx(),larmcp.Vy(), larmcp.Vz(), larmcp.T());
    res.end_position(larmcp.EndX(),larmcp.EndY(),larmcp.EndZ(),larmcp.EndT());
    res.last_step(larmcp.EndX(),larmcp.EndY(),larmcp.EndZ(),larmcp.EndT());
    res.energy_init(larmcp.E()*1000);
    res.creation_process(larmcp.Process());

    if(larmcp.Mother() > 0)
      res.parent_track_id(larmcp.Mother());
    else
      res.parent_track_id(res.track_id());
    // parent pdg code
    // parent position
    // ancestor track id
    // ancestor position
    return res;
  }

  bool SuperaMCParticleCluster::IsTouching(const Voxel3DMeta& meta, const VoxelSet& vs1, const VoxelSet& vs2) const
  {
    
    bool touching = false;
    size_t ix1,iy1,iz1;
    size_t ix2,iy2,iz2;
    size_t diffx, diffy, diffz;
    for(auto const& vox1 : vs1.as_vector()) {
      meta.id_to_xyz_index(vox1.id(), ix1, iy1, iz1);
      for(auto const& vox2 : vs2.as_vector()) {
	meta.id_to_xyz_index(vox2.id(), ix2, iy2, iz2);
	if(ix1>ix2) diffx = ix1-ix2; else diffx = ix2-ix1;
	if(iy1>iy2) diffy = iy1-iy2; else diffy = iy2-iy1;
	if(iz1>iz2) diffz = iz1-iz2; else diffz = iz2-iz1;
	touching = diffx<=1 && diffy<=1 && diffz <=1;
	if(touching) break;
      }
      if(touching) break;
    }
    return touching;
  }
  
  void SuperaMCParticleCluster::finalize()
  {}
}

#endif
