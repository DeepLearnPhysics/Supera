#ifndef __SUPERAPARMSParticle_CXX__
#define __SUPERAPARMSParticle_CXX__

#include "ParamsParticle.h"

namespace supera {

	void ParamsParticle::configure(const supera::Config_t& cfg)
	{
		_out_roi_producer = cfg.get<std::string>("OutParticleLabel", "");
	}
}
#endif
