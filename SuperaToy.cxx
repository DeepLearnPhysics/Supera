#ifndef __SUPERATOY_CXX__
#define __SUPERATOY_CXX__

#include <unordered_set>
#include <fstream>
#include "SuperaToy.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "larcv/core/DataFormat/EventVoxel3D.h"

namespace larcv {

  static SuperaToyProcessFactory __global_SuperaToyProcessFactory__;

  SuperaToy::SuperaToy(const std::string name)
    : SuperaBase(name)
  {}

  void SuperaToy::configure(const PSet& cfg)
  {
    SuperaBase::configure(cfg);
    _debug = cfg.get<bool> ("DebugMode",false);
    _hit_producer = cfg.get<std::string>("LArHitProducer","gaushit");
    _sps_producer = cfg.get<std::string>("LArSpacePointProducer","cluster3d");
    _use_true_pos = cfg.get<bool>("UseTruePosition",true);
    _ref_meta3d_cluster3d = cfg.get<std::string>("Meta3DFromCluster3D","pcluster");
    _ref_meta3d_tensor3d = cfg.get<std::string>("Meta3DFromTensor3D","");
    _nsigma_match_time = cfg.get<double>("NSigmaMatchTime",1.0);
  }


  void SuperaToy::initialize()
  {
    SuperaBase::initialize();
  }

  larcv::Voxel3DMeta SuperaToy::get_meta3d(IOManager& mgr) const {

    larcv::Voxel3DMeta meta3d;
    if(!_ref_meta3d_cluster3d.empty()) {
      auto const& ev_cluster3d = mgr.get_data<larcv::EventClusterVoxel3D>(_ref_meta3d_cluster3d);
      meta3d = ev_cluster3d.meta();
    }
    else if(!_ref_meta3d_tensor3d.empty()) {
      auto const& ev_tensor3d = mgr.get_data<larcv::EventSparseTensor3D>(_ref_meta3d_tensor3d);
      meta3d = ev_tensor3d.meta();
    }
    else
      LARCV_CRITICAL() << "ref_meta3d_cluster3d nor ref_meta3d_tensor3d set!" << std::endl;
    return meta3d;
  }

  bool SuperaToy::process(IOManager& mgr)
  {
    SuperaBase::process(mgr);
    //
    // Step 0. get 3D meta
    // Step 1. create a list of true hits
    //   - store an array of hits per channel, each hit contains time and 3D voxel ID)
    // Step 2. Create a map of reco 3D voxel <=> true 3D voxel (per plane, as a representation of matched hits)
    // Step 3. Identify valid reco3D voxels
    //   - Find an overlapping true 3D voxel ID across planes per 3D reco voxel ID
    //

    //
    // Step 0. ... get meta3d
    //
		auto event_id = mgr.event_id().event();
    LARCV_INFO() << "Retrieving 3D meta..." << std::endl;
    auto meta3d = get_meta3d(mgr);

		std::unordered_set<unsigned int> true_vox_id_v;
    //
    // Step 1. ... create a list of true hits
    //
    struct TrueHit_t {
      double time;
      std::set<unsigned int> vox_id_s;
    };
    // Get geometry info handler
    auto geop = lar::providerFrom<geo::Geometry>();
    // Create a hit list container (outer = channel number, inner = hit list per channel)
    std::vector<std::vector<TrueHit_t> > true_hit_vv;
    true_hit_vv.resize(geop->Nchannels());
    // Create anohter list w/ same shape to monitor whether true hit is matched to a reco hit
    std::vector<std::vector<bool> > true_hit_monitor_vv;
    true_hit_monitor_vv.resize(geop->Nchannels());
    // Fill a hit list container
    for(auto const& sch : LArData<supera::LArSimCh_t>()){
      // Get a unique readout channel number
      auto ch = sch.Channel();
      // Get a corresponding container
      auto& true_hit_v = true_hit_vv[ch];
      auto& true_hit_monitor_v = true_hit_monitor_vv[ch];
      // Loop over hits and store
      for (auto const tick_ides : sch.TDCIDEMap()) {
	double x_pos = (-supera::TPCTDC2Tick(tick_ides.first) + supera::TriggerOffsetTPC() - supera::TriggerTime()) * supera::TPCTickPeriod()  * supera::DriftVelocity();
	TrueHit_t hit;
	hit.time = supera::TPCTDC2Tick(tick_ides.first);
	for (auto const& edep : tick_ides.second) {
	//std::cout << supera::TPCTDC2Tick(tick_ides.first) << " " << supera::TPCTickPeriod() << " " << supera::DriftVelocity() << " " << edep.x << " " << supera::TriggerOffsetTPC() << " " << supera::PlaneTickOffset(0, 1) << " " << supera::PlaneTickOffset(1, 2) << std::endl;
	//std::cout << edep.x << " " << (-supera::TPCTDC2Tick(tick_ides.first)+supera::TriggerOffsetTPC() - supera::TriggerTime()) * supera::TPCTickPeriod()  * supera::DriftVelocity() << std::endl;
	  if(_use_true_pos) x_pos = edep.x;
	  auto vox_id = meta3d.id(x_pos, edep.y, edep.z);
	  if(vox_id == larcv::kINVALID_VOXELID) continue;
	  hit.vox_id_s.insert(vox_id);
		true_vox_id_v.insert(vox_id);
	}
	if(hit.vox_id_s.size()) {
	  true_hit_v.push_back(hit);
	  true_hit_monitor_v.push_back(false);
	}
      }
    }

    LARCV_INFO() << "Created a list of true hits: " << true_hit_vv.size() << " channels" << std::endl;
    if(_debug) {
      size_t num_hits = 0;
      for(auto const& true_hit_v : true_hit_vv) num_hits += true_hit_v.size();
      LARCV_INFO() << " ... corresponds to " << num_hits << " true hits!" << std::endl;
    }

    //
    // Step 2. ... create a map of reco 3D voxel ID <=> true 3D voxel ID per plane
    //             through reco hit <=> true hit correlation
    //
    // Key = reco 3D voxel, value = 3 set (per plane) of true 3D voxel ID
    std::map<unsigned int, std::vector<std::set<unsigned int> > > reco2true_m;
		// Key = reco 3D voxel, value = whether to keep it or not
		std::map<unsigned int, bool> reco2ghost_m;
    // create data handler (art::Handle) to retrieve reconstructed hits
    //art::Handle<std::vector<recob::Hit> > hit_h;
    //GetEvent()->getByLabel(_hit_producer, hit_h);
    // Get the list of associated pointers
    //art::FindManyP<recob::SpacePoint> ass_v(hit_h, *(GetEvent()), _sps_producer);
		art::Handle<std::vector<recob::SpacePoint> > sps_h;
		GetEvent()-> getByLabel(_sps_producer, sps_h);
		art::FindManyP<recob::Hit> ass_v(sps_h, *(GetEvent()), _sps_producer);

		std::ofstream outfile;
		outfile.open("out_"+std::to_string(event_id)+".csv");
		outfile << "id,keep,x,y,z\n";
		std::ofstream outfile_overlap;
		outfile_overlap.open("overlap_"+std::to_string(event_id)+".csv");
		outfile_overlap << "true,reco\n";

		// Loop over 3D reco space points
		for (size_t idx = 0; idx < sps_h->size(); ++idx) {
			art::Ptr<recob::SpacePoint> sps_ptr(sps_h, idx);
			// Make 3-planes vector
			std::vector<std::set<unsigned int> > reco_hits(3);
			// Get associated reco hit information
			const std::vector<art::Ptr<recob::Hit> >& ptr_coll = ass_v.at(idx);
			assert(ptr_coll.size() > 0);
			assert(ptr_coll.size() < 4);
			// Looping over planes (1 hit/plane)
			for (auto const& hit_ptr : ptr_coll) {
				double tstart = hit_ptr->PeakTime() - hit_ptr->RMS() * _nsigma_match_time;
				double tend   = hit_ptr->PeakTime() + hit_ptr->RMS() * _nsigma_match_time;
				// Get all true hits associated with this reco hit
				auto const& true_hit_v   = true_hit_vv[hit_ptr->Channel()];
				auto& true_hit_monitor_v = true_hit_monitor_vv[hit_ptr->Channel()];
				std::vector<size_t> true_hit_index_v;
				for (size_t true_hit_index = 0; true_hit_index < true_hit_v.size(); ++true_hit_index) {
					auto const& true_hit = true_hit_v[true_hit_index];
					if (true_hit.time < tstart || true_hit.time > tend) continue;
					// Now this true hit meets our time requirements.
					true_hit_monitor_v[true_hit_index] = true;
					true_hit_index_v.push_back(true_hit_index);
					for (auto const& vox_id : true_hit.vox_id_s) {
						reco_hits[hit_ptr->WireID().Plane].insert(vox_id);
					}
				}
			}
			// Find if same true vox ID appears at least across 2 different planes
			auto * xyz = sps_ptr->XYZ();
			VoxelID_t reco_vox_id = meta3d.id(xyz[0], xyz[1], xyz[2]);
			if (reco_vox_id == larcv::kINVALID_VOXELID) {
				LARCV_WARNING() << "Invalid voxel ID at " << xyz[0] << " " << xyz[1] << " " << xyz[2] << std::endl;
				continue;
			}
			if (reco2ghost_m.find(reco_vox_id) == reco2ghost_m.end() || !reco2ghost_m[reco_vox_id]) {
				reco2ghost_m[reco_vox_id] = this->findOverlap<unsigned int>(reco_hits);
			}

			// Save in CSV
			larcv::Point3D p = meta3d.position(reco_vox_id);
			outfile << reco_vox_id << "," << reco2ghost_m[reco_vox_id] << ",";
			outfile << (p.x - meta3d.min_x())/meta3d.size_voxel_x() << ",";
			outfile << (p.y - meta3d.min_y())/meta3d.size_voxel_y() << ",";
			outfile << (p.z - meta3d.min_z())/meta3d.size_voxel_z() << "\n";

			// Save overlaps
			std::vector<unsigned int> edges_true;
			std::vector<unsigned int> edges_reco;
			for (auto vox_id: reco_hits[0]) {
				if (reco_hits[1].find(vox_id) != reco_hits[1].end() && reco_hits[2].find(vox_id) != reco_hits[2].end()) {
					edges_true.push_back(vox_id);
					edges_reco.push_back(reco_vox_id);
				}
			}
			for (size_t i_edge = 0; i_edge < edges_true.size(); ++i_edge) {
				outfile_overlap << edges_true[i_edge] << "," << edges_reco[i_edge] << "\n";
			}
		}
		outfile.close();
		outfile_overlap.close();

    //LARCV_INFO() << "Created " << reco2true_m.size() << " recob::SpacePoint <=> true 3D voxel ID mapping" << std::endl;

    if(_debug) {
      size_t matched_true_hit = 0;
      size_t total_true_hit   = 0;
      for(auto const& true_hit_monitor_v : true_hit_monitor_vv) {
	for(auto const& used : true_hit_monitor_v) {
	  if(used) matched_true_hit++;
	}
	total_true_hit += true_hit_monitor_v.size();
      }
      LARCV_INFO() << "... where " << matched_true_hit << " / " << total_true_hit << " true hits are matched w/ reco hit" << std::endl;
    }

		// Save to CSV true voxels
		std::ofstream outfile2;
		outfile2.open("true_"+std::to_string(event_id)+".csv");
		outfile2 << "id,x,y,z\n";
		//auto& tensor = mgr.get_data<larcv::EventClusterVoxel3D>(_sps_producer);
		for (auto vox : true_vox_id_v) {
			larcv::Point3D p = meta3d.position(vox);
			outfile2 << vox << ",";
			outfile2 << (p.x - meta3d.min_x())/meta3d.size_voxel_x() << ",";
			outfile2 << (p.y - meta3d.min_y())/meta3d.size_voxel_y() << ",";
			outfile2 << (p.z - meta3d.min_z())/meta3d.size_voxel_z() << "\n";
		}
		outfile2.close();
    return true;
  }

	template <class T> bool SuperaToy::findOverlap(std::vector<std::set<T> > planes) {
		if (planes.size() != 3) {
			LARCV_CRITICAL() << "Need 3 planes" << std::endl;
			throw larbys();
		}

		std::unordered_set<T> plane0;
		std::unordered_set<T> plane1;
		for (auto trueVoxId : planes[0]) {
			plane0.insert(trueVoxId);
		}
		for (auto trueVoxId : planes[1]) {
			if (plane0.find(trueVoxId) == plane0.end()) {
				plane1.insert(trueVoxId);
			}
			else {
				return true;
			}
		}
		for (auto trueVoxId : planes[2]) {
			if (plane0.find(trueVoxId) != plane0.end()) {
				return true;
			}
			if (plane1.find(trueVoxId) != plane1.end()) {
				return true;
			}
		}
		return false;
	}

	//template <VoxelID_t> bool SuperaToy::findOverlap(std::vector<std::set<VoxelID_t> > planes);

  void SuperaToy::finalize()
  {}
}

#endif
