#ifndef SUPERA_MCPARTICLELIST_CXX
#define SUPERA_MCPARTICLELIST_CXX

#include "MCParticleList.h"
namespace supera {

  void MCParticleList::Update(const std::vector<simb::MCParticle>& larmcp_v, int run, int event)
  {
    if (run == _run && event == _event) return;

    _run = run;
    _event = event;

    _trackid_v.resize(larmcp_v.size());
    _pdgcode_v.resize(larmcp_v.size());
    _parent_index_v.resize(larmcp_v.size());
    _parent_trackid_v.resize(larmcp_v.size());
    _parent_pdg_v.resize(larmcp_v.size());
    _ancestor_trackid_v.resize(larmcp_v.size());
    _ancestor_index_v.resize(larmcp_v.size());

    for (size_t i = 0; i < larmcp_v.size(); ++i)
      _trackid_v[i] = _pdgcode_v[i] = _parent_index_v[i] = _parent_trackid_v[i] = _parent_pdg_v[i] =
        _ancestor_trackid_v[i] = _ancestor_index_v[i] = -1;

    _trackid2index.resize(std::max(_trackid2index.size(), larmcp_v.size()));
    for (auto& v : _trackid2index)
      v = -1;

    for (size_t index = 0; index < larmcp_v.size(); ++index) {
      auto const& mcpart = larmcp_v[index];
      _trackid_v[index] = abs(mcpart.TrackId());
      _pdgcode_v[index] = abs(mcpart.PdgCode());
      _parent_trackid_v[index] = mcpart.Mother();
      if (mcpart.TrackId() >= ((int)(_trackid2index.size())))
        _trackid2index.resize(mcpart.TrackId() + 1, -1);
      _trackid2index[mcpart.TrackId()] = index;
    }
    // Parent/Ancestor info
    for (size_t index = 0; index < larmcp_v.size(); ++index) {
      auto const& mcpart = larmcp_v[index];
      int mother_id = mcpart.Mother();
      int mother_index = -1;
      if (mother_id == 0) mother_id = abs(mcpart.TrackId());
      if (mother_id < ((int)(_trackid2index.size()))) {
        mother_index = _trackid2index[mother_id];
        if (mother_index >= 0) {
          _parent_pdg_v[index] = larmcp_v[mother_index].PdgCode();
          _parent_index_v[index] = mother_index;
        }
      }

      int subject_track_id = abs(mcpart.TrackId());
      int parent_track_id = abs(mcpart.Mother());
      int ancestor_index = -1;
      int ancestor_track_id = -1;
      while (1) {
        if ((size_t)(parent_track_id) >= _trackid2index.size()) break;
        if (parent_track_id == 0 || parent_track_id == subject_track_id) {
          ancestor_index = _trackid2index[subject_track_id];
          ancestor_track_id = subject_track_id;
          break;
        }
        auto const& parent_index = _trackid2index[parent_track_id];
        if (parent_index < 0) break;
        auto const& parent = larmcp_v[parent_index];
        subject_track_id = abs(parent.TrackId());
        parent_track_id = abs(parent.Mother());
      }
      _ancestor_index_v[index] = ancestor_index;
      _ancestor_trackid_v[index] = ancestor_track_id;
    }
  }
}
#endif
