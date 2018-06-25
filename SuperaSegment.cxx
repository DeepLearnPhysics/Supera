#ifndef __SUPERASEGMENT_CXX__
#define __SUPERASEGMENT_CXX__


#include "SuperaSegment.h"
#include "Instance2Image.h"
#include "ImageMetaMakerFactory.h"
#include "PulledPork3DSlicer.h"
#include "larcv/core/DataFormat/EventImage2D.h"
   //#include "larcv/core/DataFormat/DataFormatUtil.h"

namespace larcv {

   static SuperaSegmentProcessFactory __global_SuperaSegmentProcessFactory__;

   SuperaSegment::SuperaSegment(const std::string name): SuperaBase(name) {}

   void SuperaSegment::configure(const PSet & cfg) {
      SuperaBase::configure(cfg);
      supera::ParamsImage2D::configure(cfg);
      supera::ImageMetaMaker::configure(cfg);
      _origin = cfg.get <unsigned short>  ("Origin", 0);
      m_ancestor_label = cfg.get <std::string> ("AncestorImageLabel");
      m_instance_label = cfg.get <std::string> ("InstanceImageLabel");
      m_wire_label = cfg.get <std::string> ("WireImageLabel");
   }

   void SuperaSegment::initialize() {}

   bool SuperaSegment::process(IOManager & mgr) {

      std::cout << "Hello World!";
      //Part One, Get the Images we Need
      SuperaBase::process(mgr);

      if (supera::PulledPork3DSlicer::Is(supera::ImageMetaMaker::MetaMakerPtr())) {
         auto ptr = (supera::PulledPork3DSlicer * )(supera::ImageMetaMaker::MetaMakerPtr());
         ptr->ClearEventData();
         ptr->AddConstraint(LArData<supera::LArMCTruth_t> ());
         ptr->GenerateMeta(LArData<supera::LArSimCh_t> (), TimeOffset());
      }

      // get meta for image we initially fill
      auto
      const & filler_meta_v = Meta();
      // get meta of output image (the shared meta)

      if (filler_meta_v.empty()) {
         LARCV_CRITICAL() << "Filler meta not created!" << std::endl;
         throw larbys();
      }

      auto
      const ev_image = (EventImage2D * )(mgr.get_data("image2d", OutImageLabel()));
      if (!ev_image) {
         LARCV_CRITICAL() << "Output image could not be created!" << std::endl;
         throw larbys();

      }

      if (!(ev_image->image2d_array().empty())) {
         LARCV_CRITICAL() << "Output image array not empty!" << std::endl;
         throw larbys();
      }

      // change to import instance and wire
      const EventImage2D * ev_instance = (EventImage2D * )(mgr.get_data("image2d", m_instance_label));
      if (!ev_instance) {
         LARCV_CRITICAL() << "ev_instance image could not be created!" << std::endl;
         throw larbys();
      }
      const std::vector<larcv::Image2D> & instance_v = ev_instance-> image2d_array();
      //Not sure this should be empty, commented out warning --J
      // if(!(ev_ancestor->image2d_array().empty())) {
      //   LARCV_CRITICAL() << "Ancestor image array not empty!" << std::endl;
      //   throw larbys();
      // }

      //add another get data block to import wire image
      const EventImage2D * ev_wire = (EventImage2D * )(mgr.get_data("image2d", m_wire_label));
      if (!ev_wire) {
         LARCV_CRITICAL() << "ev_wire image could not be created!" << std::endl;
         throw larbys();
      }
      const std::vector<larcv::Image2D> & wire_v = ev_wire->image2d_array();

      //Part 2 Make Blank Images for Labeling
      //This loop is meant to go through and make a blank image available for labels
      // for(auto const& i:instance_v)
      //   {
      //     larcv::Image2D x(i.meta());
      //     x.paint(0.0);
      //     ev_image->emplace(std::move(x));
      //   }

      //auto image_v = ev_image->image2d_array();   //   const std::vector<larcv::Image2D>&  ?

      //Part 3&4 Fill image with correct labels

      int nbackground = 0;
      int ntrack = 0;
      int nshower = 0;

      std::vector <larcv::Image2D> image_v = {
         instance_v[0].meta(),
         instance_v[1].meta(),
         instance_v[2].meta()
      };

      int rows2 = image_v[0].meta().rows();
      int cols2 = image_v[0].meta().cols();

      // Let's make a map! Pix/Track ID to label (Track, shower, background)
      std::map <int, int> trackid2label;
      for (auto
         const & mctrack: LArData <supera::LArMCTrack_t> ()) {
         int trackID = static_cast <int> (mctrack.TrackID());
         trackid2label[trackID] = 1;
      }

      for (int plane = 0; plane < 3; plane++) {
         // larcv::Image2D image(instance_v[plane].meta());
         image_v[plane].paint(0.0);
         // size_t Instance_rows = instance_v[plane].meta().rows();
         // size_t Instance_cols = instance_v[plane].meta().cols();

         // std::cout << "Instance rows = " << Instance_rows  << std::endl;
         // std::cout << "Instance cols = " << Instance_cols<< std::endl;

         // size_t Wire_rows = wire_v[plane].meta().rows();
         // size_t Wire_cols = wire_v[plane].meta().cols();

         // std::cout << "Wire rows = " << Wire_rows  << std::endl;
         // std::cout << "Wire cols = " << Wire_cols<< std::endl;

         size_t rows = image_v[plane].meta().rows();
         size_t cols = image_v[plane].meta().cols();

         // std::cout << "rows = " << rows  << std::endl;
         // std::cout << "cols = " << cols;

         //Loop through the rows

         for (unsigned int r = 0; r < rows; r++) {

            //loop through Columns
            for (unsigned int c = 0; c < cols; c++) {

               if (wire_v[plane].pixel(r, c) > 5) //Threshold checker
               {

                  //Execute Labeling
                  int pix_id = static_cast<int> (instance_v[plane].pixel(r, c));

                  bool is_track = 0;

                  if (trackid2label.find(pix_id) != trackid2label.end()) {
                     is_track = 1;
                     image_v[plane].set_pixel(r, c, 1); //Track
                     ntrack += 1;
                  }

                  if (is_track == 0 && pix_id >= 0) {
                     image_v[plane].set_pixel(r, c, 2); // Shower

                  }

               } else {
                  image_v[plane].set_pixel(r, c, 0); // Background

               }

               //counter pixels in each category
               if (image_v[plane].pixel(r, c) == 0) {
                  nbackground += 1;
               } else if (image_v[plane].pixel(r, c) != 1) {
                  nshower += 1;
               }

            } //end column loop
         } //end row loop

      } //end plane loop
      //
      // std::cout << "RIGHT AFTER LABELING: Num Tracks = " << ntrack << "     Num Shower = " << nshower << "        Num Background = " << nbackground << std::endl;
      //
      //     ntrack = 0;
      //     nshower = 0;
      //     nbackground = 0;
      //
      //     int pix_label5;
      //
      //     for (int p5 = 0; p5<3 ; p5++)
      //     {
      //       size_t rows5 = image_v[p5].meta().rows();
      //       size_t cols5 = image_v[p5].meta().cols();
      //       for (unsigned int r5 =0; r5<rows5;r5++)
      //         {
      //         for (unsigned int c5 =0; c5<cols5;c5++)
      //           {
      //           pix_label5 = static_cast<int>(image_v[p5].pixel(r5,c5));
      //           if (pix_label5 ==0) nbackground +=1;
      //           if (pix_label5 ==1) ntrack +=1;
      //           if (pix_label5 ==2) nshower +=1;
      //
      //
      //           }
      //         }
      //     }
      //std::cout << "AFTER FIRST RESET & RUN THROUGH Num Tracks = " << ntrack << "     Num Shower = " << nshower << "        Num Background = " << nbackground << std::endl;
      //Now Go Through tracks and label track ends
      int ntracks = 0;

      for (auto const & mctrack2: LArData <supera::LArMCTrack_t> ()) {
         ntracks += 1;

         //std::cout << ntracks << std::endl;

         if ((mctrack2.empty() != 1)) //(mctrack2.Origin() == 2) &&
         {
            //std::cout << ntracks << "  Origin: " << mctrack2.Origin()<< std::endl;
            //std::cout << std::endl << "I'm labeling track ends! at points:    " ;
            //istart and iend are 4 vectors of ints, (r, ucol, vcol, ycol)
            std::vector<int> istart = {};
            std::vector<int> iend = {};

            if (doesTrackCrossImageBoundary(mctrack2, image_v[0].meta(), 3200) != -2) {

               istart = getFirstStepPosInsideImage(mctrack2, image_v[0].meta(), 4050, true, 0.15, 2.0); //ev_trigger->TriggerTime()
               if (istart.empty() !=1) {std::cout << istart[0] << "   " << istart[1] << "   " << istart[2] << "   " << istart[3] << std::endl;}
               else {std::cout << "Start was empty" <<std::endl;}
               iend = getFirstStepPosInsideImage(mctrack2, image_v[0].meta(), 4050, false, 0.15, 2.0); //ev_trigger->TriggerTime()
               if (iend.empty() != 1) {std::cout << iend[0] << "   " << iend[1] << "   " << iend[2] << "   " << iend[3] << std::endl;}
               else {std::cout << "End was empty" <<std::endl;}
            }
            //Get track ends if is wholly in detector
            // else if (doesTrackCrossImageBoundary(mctrack2, image_v[0].meta(), 3200) == -1) {
            //   std::cout << "A Track that doesn't cross: " << mctrack2.TrackID() << std::endl;
              //
              // const double cm_per_tick = supera::DriftVelocity() * 0.5;
              // int length = mctrack2.size();
              // double xstart = mctrack2.at(0).X();
              // double ystart = mctrack2.at(0).Y();
              // double zstart = mctrack2.at(0).Z();
              // double xend = mctrack2.at(length-1).X();
              // double yend = mctrack2.at(length-1).Y();
              // double zend = mctrack2.at(length-1).Z();
              //
              // std::vector<double> pos_start(3);
              // pos_start[0] = xstart;
              // pos_start[1] = ystart;
              // pos_start[2] = zstart;
              // supera::ApplySCE(pos_start[0], pos_start[1], pos_start[2]);
              // std::vector<double> pos_end(3);
              // pos_end[0] = xend;
              // pos_end[1] = yend;
              // pos_end[2] = zend;
              // supera::ApplySCE(pos_end[0], pos_end[1], pos_end[2]);
              // std::vector<double> pos4vstart(4, 0);
              // pos4vstart = {mctrack2.at(0).T(), pos_start[0], pos_start[1], pos_start[2]};
              // pos4vend = {mctrack2.at(length-1).T(), pos_end[0], pos_end[1], pos_end[2]};
              //
              // std::vector<double> pos4vend(4, 0);
              // double tickstart = getTick( , 4050.0);
              // double tickend = getTick(pos4v, 4050.0);
              // // if (tick < meta.min_y() + 20.0 || tick > meta.max_y() - 20.0)
              // //    continue;
              //
              // pos_start[0] = (tickstart - 3200.0) * cm_per_tick; // we have to give the apparent-x (relative to the trigger) because we need to test the position in the image
              // pos_end[0] = (tickend - 3200.0) * cm_per_tick
              // std::vector<int> imgcoords;
              // istart = getProjectedImagePixel(pos_start, image_v[0].meta(), 3, 0);
              // iend = getProjectedImagePixel(pos_end, image_v[].meta(), 3, 0);



            // }

            if ((istart.empty() != 1) && (iend.empty() != 1)) {

               // std::cout << istart[0] << "   " << istart[1] << "   " << istart[2] << "   " << istart[3] << std::endl;
               // std::cout << iend[0] << "   " << iend[1] << "   " << iend[2] << "   " << iend[3] << std::endl;

               for (int k = 0; k < 3; k++) {
                  //Actual Point
                  image_v[k].set_pixel(istart[0], istart[k + 1], 3); //Track End
                  image_v[k].set_pixel(iend[0], iend[k + 1], 3); //Track End

                  // Make some more pixels around it labeled that
                  int box_dim = 10;
                  for (int a = 0 - (box_dim / 2); a < 0 + (box_dim / 2) + 1; a++) {
                     for (int b = 0 - (box_dim / 2); b < 0 + (box_dim / 2) + 1; b++) {
                        if ((istart[0] + a >= 0) && (istart[0] + a < rows2) && (istart[k + 1] + b >= 0) && (istart[k + 1] + b < cols2)) { //Check to make sure pixel in image bounds
                          if (wire_v[k].pixel(istart[0] + a, istart[k + 1] + b) >5)
                           {image_v[k].set_pixel(istart[0] + a, istart[k + 1] + b, 3);}
                        }
                        if ((iend[0] + a >= 0) && (iend[0] + a < rows2) && (iend[k + 1] + b >= 0) && (iend[k + 1] + b < cols2)) { //Check to make sure pixel in image bounds
                          if (wire_v[k].pixel(iend[0] + a, iend[k + 1] + b) > 5)
                           {image_v[k].set_pixel(iend[0] + a, iend[k + 1] + b, 3);}
                        }
                     }
                  }

                  // image_v[k].set_pixel(istart[0]+1,istart[k+1],3);
                  // image_v[k].set_pixel(istart[0]-1,istart[k+1],3);
                  // image_v[k].set_pixel(istart[0],istart[k+1]+1,3);
                  // image_v[k].set_pixel(istart[0],istart[k+1]-1,3);
                  // image_v[k].set_pixel(iend[0]+1,iend[k+1],3);
                  // image_v[k].set_pixel(iend[0]-1,iend[k+1],3);
                  // image_v[k].set_pixel(iend[0],iend[k+1]+1,3);
                  // image_v[k].set_pixel(iend[0],iend[k+1]-1,3);

               }//End of Loop through planes
            }//End of IF ISTART and IEND not empty

         }//End of if track not empty

      }//End of loop through all tracks
      std::cout << "1st. Num Tracks = " << ntrack << "     Num Shower = " << nshower << "        Num Background = " << nbackground << std::endl;
      ntrack = 0;
      nshower = 0;
      nbackground = 0;
      int nends = 0;
      int pix_label;

      for (int p = 0; p < 3; p++) {

         for (int r2 = 0; r2 < rows2; r2++) {
            for (int c2 = 0; c2 < cols2; c2++) {
               pix_label = static_cast<int> (image_v[p].pixel(r2, c2));
               if (pix_label == 0) nbackground += 1;
               if (pix_label == 1) ntrack += 1;
               if (pix_label == 2) nshower += 1;
               if (pix_label == 3) nends += 1;

            }
         }
      }

      ev_image->emplace(std::move(image_v));

      std::cout << "2nd. Num Tracks = " << ntrack << "     Num Shower = " << nshower << "        Num Background = " << nbackground << "      Num End = " << nends << std::endl;
      std::cout << ntrack + nshower + nbackground + nends << std::endl;

      return true;

   }

   void SuperaSegment::finalize() {}

   double SuperaSegment::getTick(const supera::LArMCStep_t & step,
      const double trig_time) {
      std::vector<double> pos(4, 0);
      pos[0] = step.T();
      pos[1] = step.X();
      pos[2] = step.Y();
      pos[3] = step.Z();
      return getTick(pos, trig_time);
   }

   double SuperaSegment::getTick(const std::vector<double> & step,
      const double trig_time) {
      // Function returns the tick time of a MC step point
      // if SCE pointer is null, we do not correct for the space charge

      std::vector<double> scepos(3, 0);
      scepos[0] = step[1];
      scepos[1] = step[2];
      scepos[2] = step[3];
      supera::ApplySCE(scepos[0], scepos[1], scepos[2]);

      const double cm_per_tick = supera::DriftVelocity() * 0.5;
      double tick = (step[0] * 1.0e-3 - (trig_time - 4050.0)) / 0.5 + scepos[0] / cm_per_tick + 3200.0;

      return tick;
   }

   std::vector<int> SuperaSegment::getProjectedImagePixel(const std::vector<double> & pos3d,
      const larcv::ImageMeta & meta,
         const int nplanes,
            const double fracpixborder) {
      /* -----------------------------------------------------------
       * returns (row,colp1,colp2,..,colppN) corresponding to
       *  projection of pos3d onto images
       *
       * inputs
       * ------
       * pos3d: expect length 3, position in 3D space
       * meta: ImageMeta, defining image size and  coordinates
       * nplanes: number of planes
       * fracpixborder: 3D positions that fall outside of image but within
       *  fracpixboarder*pixel width or height return col or row of nearest
       *  valid pixel. otherwise return -1 value for outside the image
       *
       * output: 1+numplanes vector providing row and col of images for the planes
       *  if outside the image, values of -1 are returned
       *
       * weaknesses: uses hard-coded values for position to tick conversion, should use service
       *
       * ----------------------------------------------------------*/

      std::vector<int> img_coords(nplanes + 1, -1);
      double row_border = fabs(fracpixborder) * meta.pixel_height();
      double col_border = fabs(fracpixborder) * meta.pixel_width();

      // tick/row
      double tick = pos3d[0] / (0.5 * supera::DriftVelocity()) + 3200.0;
      if (tick < meta.min_y()) {
         if (tick > meta.min_y() - row_border)
         // below min_y-border, out of image
            img_coords[0] = meta.rows() - 1; // note that tick axis and row indicies are in inverse order
         else
         // outside of image and border
            img_coords[0] = -1;
      } else if (tick > meta.max_y()) {
         if (tick < meta.max_y() + row_border)
         // within upper border
            img_coords[0] = 0;
         else
         // outside of image and border
            img_coords[0] = -1;
      } else {
         // within the image
         img_coords[0] = meta.row(tick);
      }

      // Columns
      Double_t xyz[3] = {
         pos3d[0],
         pos3d[1],
         pos3d[2]
      };
      // there is a corner where the V plane wire number causes an error
      if ((pos3d[1] > -117.0 && pos3d[1] < -116.0) && pos3d[2] < 2.0) {
         std::cout << __PRETTY_FUNCTION__ << ": v-plane corner hack (" << xyz[0] << "," << xyz[1] << "," << xyz[2] << ")" << std::endl;
         xyz[1] = -116.0;
      }

      for (int p = 0; p < nplanes; p++) {
         int wire = supera::NearestWire(xyz, p);
         // round wire
         //wire = std::roundf(wire);

         // get image coordinates
         if (wire < meta.min_x()) {
            if (wire > meta.min_x() - col_border) {
               std::cout << __PRETTY_FUNCTION__ << " plane=" << p << " wire=" << wire << "<" << meta.min_x() - col_border << std::endl;
               // within lower border
               img_coords[p + 1] = 0;
            } else
               img_coords[p + 1] = -1;
         } else if (wire >= meta.max_x()) {
            if (wire < meta.max_x() + col_border) {
               std::cout << __PRETTY_FUNCTION__ << " plane=" << p << " wire=" << wire << ">" << meta.max_x() + col_border << std::endl;
               // within border
               img_coords[p + 1] = meta.cols() - 1;
            } else
            // outside border
               img_coords[p + 1] = -1;
         } else
         // inside image
            img_coords[p + 1] = meta.col(wire);
      } //end of plane loop

      // there is a corner where the V plane wire number causes an error
      if (pos3d[1] < -116.3 && pos3d[2] < 2.0 && img_coords[1 + 1] == -1) {
         img_coords[1 + 1] = 0;
      }

      return img_coords;
   }

   std::vector < int > SuperaSegment::getFirstStepPosInsideImage(const supera::LArMCTrack_t & track, const larcv::ImageMeta & meta,
            const double trig_time, const bool startAtstart, const double max_step_size, const double fv_border) {
      // This function returns the (SCE-corrected) position where a MC track first is inside the image bounds
      std::cout << "Track ID: " << track.TrackID() << std::endl;
      double dwall_min = 1.0e9;
      std::vector < double > pos_min(3);
      const double cm_per_tick = supera::DriftVelocity() * 0.5;
      int npts = (int) track.size();
      //Begin Loop through points of track
      for (int ipt = 1; ipt < npts; ipt++) {

         int thispt = ipt;
         int lastpt = thispt - 1;

         if (!startAtstart) {
            thispt = npts - 1 - ipt;
            lastpt = thispt + 1;
         }

         const auto & this_step = track.at(thispt);
         const auto & last_step = track.at(lastpt);

         double dir[3] = {
            double(this_step.X() - last_step.X()),
            double(this_step.Y() - last_step.Y()),
            double(this_step.Z() - last_step.Z())
         };
         double dirnorm = 0;

         for (int i = 0; i < 3; i++) {
            dirnorm += dir[i] * dir[i];
         }

         dirnorm = sqrt(dirnorm);
         if (dirnorm < 1.0e-3)
            continue;

         for (int i = 0; i < 3; i++) {
            dir[i] /= dirnorm;
         }

         int nsteps = dirnorm / max_step_size + 1;
         if (nsteps <= 0)
            nsteps = 1;
         double stepsize = dirnorm / double(nsteps);
         for (int istep = 0; istep < nsteps; istep++) {
            std::vector<double> pos(4, 0);
            std::vector<double> pos4v(4, 0);
            pos4v[0] = last_step.T();
            pos[0] = pos4v[1] = last_step.X() + stepsize * double(istep) * dir[0];
            pos[1] = pos4v[2] = last_step.Y() + stepsize * double(istep) * dir[1];
            pos[2] = pos4v[3] = last_step.Z() + stepsize * double(istep) * dir[2];

            // get sce-corrected 3d pos
            std::vector<double> pos_sce(3);
            pos_sce[0] = pos[0];
            pos_sce[1] = pos[1];
            pos_sce[2] = pos[2];
            supera::ApplySCE(pos_sce[0], pos_sce[1], pos_sce[2]);

            int boundary_type = -1;

            double fdwall = dwall(pos_sce, boundary_type); // use apparent distance...
            //indetector unused at present
            //bool indetector = ((pos_sce[0] <258.0) && (pos_sce[0] > 0.0) && (pos_sce[1] > -117.0) && (pos_sce[1] < 117.0) && (pos_sce[2] < 1036.0) && (pos_sce[2] >0.0));
            if (fdwall < fv_border) {
               if (fdwall < dwall_min){
                  dwall_min = fdwall;
                  pos_min = pos_sce;
                  //std::cout << "Boundary " << boundary_type << " was within the distance cutoff of " << fv_border << " Away by " << dwall_min << " With coordinates of:   " << pos_min[0] << "   " << pos_min[1] << "   " << pos_min[2] << std::endl;

               }
               continue;
            }
            // convert to image coordinates
            double tick = getTick(pos4v, trig_time);
            if (tick < meta.min_y() + 20.0 || tick > meta.max_y() - 20.0)
               continue;

            pos_sce[0] = (tick - 3200.0) * cm_per_tick; // we have to give the apparent-x (relative to the trigger) because we need to test the position in the image

            std::vector<int> imgcoords;
            imgcoords = getProjectedImagePixel(pos_sce, meta, 3, 0);

            // try {
            //   imgcoords = getProjectedImagePixel( pos_sce, meta, 3 );
            //   //std::cout << " imgcoords=(row=" << imgcoords[0] << "," << imgcoords[1] << "," << imgcoords[2] << "," << imgcoords[3] << ")" << std::endl;
            // }
            // catch (...) {
            //   //std::cout << std::endl;
            //   continue;
            // }

            // return pos_sce;
            return imgcoords;
         }

      } //end of track point loop

      // didn't find the crossing boundary
      std::cout << "XXXXXXXXX DIDNT FIND BOUNDARY XXXXXXXXXXXXXXX" << std::endl;
      std::vector <int> empty;
      return empty;
   }

   int SuperaSegment::doesTrackCrossImageBoundary(const supera::LArMCTrack_t & track,
      const larcv::ImageMeta & meta,
         const double trig_time) {
      //Returns -1 if all inside
      //Returns 0 if crosses boundary 1 time
      //Return +1 if crosses boundary 2 times
      double tick_start = getTick(track.front(), trig_time);
      double tick_end = getTick(track.back(), trig_time);
      if (tick_start > meta.min_y() && tick_start < meta.max_y() && tick_end > meta.min_y() && tick_end < meta.max_y())
         return -1;

      if (tick_start < meta.min_y() && tick_end > meta.min_y())
         return 0; // start out -> end in
      else if (tick_start > meta.min_y() && tick_end < meta.min_y())
         return 1; // start in -> end out
      else if (tick_start < meta.max_y() && tick_end > meta.max_y())
         return 1; // start in -> end out;
      else if (tick_start > meta.max_y() && tick_end < meta.max_y())
         return 0; // start out -> end in

      return -1;
   }

   double SuperaSegment::dwall(const std::vector<double> & pos, int & boundary_type) {
      std::vector < float > fpos(3);
      for (int i = 0; i < 3; i++)
         fpos[i] = (float)pos[i];
      return dwall(fpos, boundary_type);
   }

   float SuperaSegment::dwall(const std::vector<float> & pos, int & boundary_type) {

      float dx1 = pos[0];
      float dx2 = 258 - pos[0];
      float dy1 = 117.0 - pos[1];
      float dy2 = pos[1] + 117.0;
      float dz1 = pos[2];
      float dz2 = 1036.0 - pos[2];

      float dwall = 1.0e9;

      if (dy1 < dwall) {
         dwall = dy1;
         boundary_type = 0; // top
      }
      if (dy2 < dwall) {
         dwall = dy2;
         boundary_type = 1; // bottom
      }
      if (dz1 < dwall) {
         dwall = dz1;
         boundary_type = 2; // upstream
      }
      if (dz2 < dwall) {
         dwall = dz2;
         boundary_type = 3; // downstream
      }
      if (dx1 < dwall) {
         dwall = dx1;
         boundary_type = 4; // anode
      }
      if (dx2 < dwall) {
         dwall = dx2;
         boundary_type = 5; // cathode
      }

      return dwall;

   }

} //End of namespace
#endif
