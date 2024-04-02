#ifndef __SUPERACSVREADER_H__
#define __SUPERACSVREADER_H__
#include "SuperaTypes.h"
#include <array>
#include <fstream>
#include <iostream>
#include <map>
namespace supera {

  namespace csvreader {

    inline bool line_to_point(const std::string& line_data,
                              supera::RSEID& id,
                              std::array<double, 3>& point)
    {
      size_t run, subrun, event, x, y, z;
      run = subrun = event = x = y = z = 0;
      try {
        run = line_data.find(",");
        subrun = line_data.find(",", run + 1);
        event = line_data.find(",", subrun + 1);
        x = line_data.find(",", event + 1);
        y = line_data.find(",", x + 1);
        z = line_data.find(" ", y + 1);

        id.run = std::atoi(line_data.substr(0, run).c_str());
        id.subrun = std::atoi(line_data.substr(run + 1, subrun - (run + 1)).c_str());
        id.event = std::atoi(line_data.substr(subrun + 1, event - (subrun + 1)).c_str());
        point[0] = std::atof(line_data.substr(event + 1, x - (event + 1)).c_str());
        point[1] = std::atof(line_data.substr(x + 1, y - (x + 1)).c_str());
        point[2] = std::atof(line_data.substr(y + 1, z - (y + 1)).c_str());
      }
      catch (const std::exception& err) {
        std::cerr << "Failed to convert: \"" << line_data << "\"" << std::endl;
        return false;
      }
      return true;
    }

    inline void read_constraint_file(std::string fname,
                                     std::map<supera::RSEID, std::array<double, 3>>& data)
    {

      std::ifstream fstrm(fname.c_str());
      std::string line_data;
      supera::RSEID id;
      std::array<double, 3> position;

      std::getline(fstrm, line_data); // ignore 1st line
      while (std::getline(fstrm, line_data)) {

        if (line_data.empty()) continue;

        if (supera::csvreader::line_to_point(line_data, id, position)) {

          auto iter = data.find(id);
          if (iter != data.end()) {
            std::cerr << "Run=" << id.run << " Subrun=" << id.subrun << " Event=" << id.event
                      << " is duplicated!" << std::endl;
            throw std::exception();
          }
          data[id] = position;
        }
      }
    }

  }
}
#endif
