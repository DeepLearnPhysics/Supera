#ifndef SUPERA_MCPARTICLELIST_CXX
#define SUPERA_MCPARTICLELIST_CXX

#include "MCParticleList.h"
namespace supera{ 

  void MCParticleList::Update(const std::vector<simb::MCParticle>& larmcp_v, int run, int event)
  {
    if(run == _run && event == _event) return;
    
    _run = run;
    _event = event;
    
    _trackid_v.resize(larmcp_v.size());
    _pdgcode_v.resize(larmcp_v.size());
    _parent_index_v.resize(larmcp_v.size());
    _parent_trackid_v.resize(larmcp_v.size());
    _parent_pdg_v.resize(larmcp_v.size());

    _trackid2index.resize(std::max(_trackid2index.size(),larmcp_v.size()));
    for(auto& v : _trackid2index) v = -1;
    
    for(size_t index=0; index<larmcp_v.size(); ++index) {
      auto const& mcpart = larmcp_v[index];
      _trackid_v[index] = abs(mcpart.TrackId());
      _pdgcode_v[index] = abs(mcpart.PdgCode());
      _parent_trackid_v[index] = mcpart.Mother();
      if(mcpart.TrackId() >= ((int)(_trackid2index.size()))) _trackid2index.resize(mcpart.TrackId()+1,-1);
      _trackid2index[mcpart.TrackId()] = index;
    }

    for(size_t index=0; index<larmcp_v.size(); ++index) {
      auto const& mcpart = larmcp_v[index];
      int mother_id      = mcpart.Mother();
      int mother_index   = -1;
      if(mother_id == 0) mother_id = mcpart.TrackId();
      if(mother_id < ((int)(_trackid2index.size()))) {
	mother_index = _trackid2index[mother_id];
	if(mother_index >= 0) {
	  _parent_pdg_v[index] = larmcp_v[mother_index].PdgCode();
	  _parent_index_v[index] = mother_index;
	}
      }
    }
  }
}
#endif
