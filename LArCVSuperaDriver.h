/**
 * \file LArCVSuperaDriver.h
 *
 * \ingroup MeatSlicer
 *
 * \brief Class def header for a class LArCVSuperaDriver
 *
 * @author kazuhiro
 */

/** \addtogroup MeatSlicer

    @{*/
#ifndef LARCVSUPERADRIVER_H
#define LARCVSUPERADRIVER_H
//#ifndef __CINT__
//#ifndef __CLING__
#include "FMWKInterface.h"
#include "SuperaBase.h"
#include "SuperaTypes.h"
#include "larcv/core/Base/larcv_base.h"
#include "larcv/core/Processor/ProcessDriver.h"
#include <iostream>
//#include "SuperaChStatus.h"

namespace larcv {
  /**
     \class LArCVSuperaDriver
     User defined class LArCVSuperaDriver ... these comments are used to generate
     doxygen documentation!
  */
  class LArCVSuperaDriver : public larcv_base {

  public:
    /// Default constructor
    LArCVSuperaDriver();

    /// Default destructor
    ~LArCVSuperaDriver() {}

    void configure(const std::string cfg_file);

    void configure(const larcv::PSet& cfg);

    void override_input_file(const std::vector<std::string>& flist);

    void override_output_file(const std::string fname);

    const std::set<std::string>& DataLabels(supera::LArDataType_t type) const;

    template <class T>
    void SetDataPointer(const std::vector<T>& data, const std::string label)
    {
      LARCV_INFO() << "Received LArDataType_t " << (int)(supera::LArDataType<T>()) << " by label "
                   << label << std::endl;
      for (auto& idx : _supera_idx_v) {
        auto supera_ptr = (SuperaBase*)(_driver.process_ptr(idx));
        if (supera_ptr->LArDataLabel(supera::LArDataType<T>()) != label) {
          LARCV_DEBUG() << supera_ptr->name() << " does not need LArDataType_t "
                        << (int)(supera::LArDataType<T>()) << " label " << label << std::endl;
          continue;
        }
        LARCV_INFO() << "Data pointer provided to: " << supera_ptr->name() << std::endl;
        supera_ptr->LArData(data);
      }
    }

    //SuperaChStatus* SuperaChStatusPointer();

    void initialize();
    bool process(size_t run, size_t subrun, size_t event);
    void finalize();

    void SetCSV(std::string proc_name, std::string fname);

    inline std::vector<std::string> ProcessNames() const { return _driver.process_names(); }

    // FIXME(kvtsang) Temporary solution to access associations
    void SetEvent(const art::Event* ev)
    {
      LARCV_INFO() << "Passing art::Event to SuperaBase classes...";
      for (auto& idx : _supera_idx_v) {
        auto supera_ptr = (SuperaBase*)(_driver.process_ptr(idx));
        supera_ptr->SetEvent(ev);
      }
    }

  private:
    ProcessDriver _driver;
    std::vector<size_t> _supera_idx_v;

    std::map<supera::LArDataType_t, std::set<std::string>> _data_request_m;
    //SuperaChStatus* _supera_chstatus_ptr;
  };

}
#endif
//#endif
//#endif
/** @} */ // end of doxygen group
