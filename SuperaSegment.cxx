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

  SuperaSegment::SuperaSegment(const std::string name)
    : SuperaBase(name)
  {}

  void SuperaSegment::configure(const PSet& cfg)
  {
    SuperaBase::configure(cfg);
    supera::ParamsImage2D::configure(cfg);
    supera::ImageMetaMaker::configure(cfg);
    _origin = cfg.get<unsigned short>("Origin",0);
    m_ancestor_label = cfg.get<std::string>("AncestorImageLabel");
    m_instance_label = cfg.get<std::string>("InstanceImageLabel");
    m_wire_label = cfg.get<std::string>("WireImageLabel");
  }

  void SuperaSegment::initialize()
  {}

  bool SuperaSegment::process(IOManager& mgr)
  {

    std::cout << "Hello World!";
//Part One, Get the Images we Need
    SuperaBase::process(mgr);


    if(supera::PulledPork3DSlicer::Is(supera::ImageMetaMaker::MetaMakerPtr())) {
      auto ptr = (supera::PulledPork3DSlicer*)(supera::ImageMetaMaker::MetaMakerPtr());
      ptr->ClearEventData();
      ptr->AddConstraint(LArData<supera::LArMCTruth_t>());
      ptr->GenerateMeta(LArData<supera::LArSimCh_t>(),TimeOffset());
    }

    // get meta for image we initially fill
    auto const& filler_meta_v = Meta();
    // get meta of output image (the shared meta)


    if(filler_meta_v.empty()) {
      LARCV_CRITICAL() << "Filler meta not created!" << std::endl;
      throw larbys();
    }

    auto const ev_image = (EventImage2D*)(mgr.get_data("image2d",OutImageLabel()));
    if(!ev_image) {
      LARCV_CRITICAL() << "Output image could not be created!" << std::endl;
      throw larbys();

    }

    if(!(ev_image->image2d_array().empty())) {
      LARCV_CRITICAL() << "Output image array not empty!" << std::endl;
      throw larbys();
    }



// change to import instance and wire
    const EventImage2D* ev_instance = (EventImage2D*)(mgr.get_data("image2d",m_instance_label));
    if(!ev_instance) {
      LARCV_CRITICAL() << "ev_instance image could not be created!" << std::endl;
      throw larbys();
    }
    const std::vector<larcv::Image2D>& instance_v = ev_instance->image2d_array();
    //Not sure this should be empty, commented out warning --J
    // if(!(ev_ancestor->image2d_array().empty())) {
    //   LARCV_CRITICAL() << "Ancestor image array not empty!" << std::endl;
    //   throw larbys();
    // }

    //add another get data block to import wire image
    const EventImage2D* ev_wire = (EventImage2D*)(mgr.get_data("image2d",m_wire_label));
    if(!ev_wire) {
      LARCV_CRITICAL() << "ev_wire image could not be created!" << std::endl;
      throw larbys();
    }
    //const std::vector<larcv::Image2D>& wire_v = ev_wire->image2d_array();

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

    int nbackground =0;
    int ntrack =0;
    int nshower =0;

    std::vector<larcv::Image2D> image_v = {instance_v[0].meta() , instance_v[1].meta() , instance_v[2].meta()};

    int rows2 = image_v[0].meta().rows();
    int cols2 = image_v[0].meta().cols();

    // Let's make a map! Pix/Track ID to label (Track, shower, background)
    std::map<int,int> trackid2label;
    for(auto const& mctrack : LArData<supera::LArMCTrack_t>())
      {
        int trackID = static_cast<int>(mctrack.TrackID());
      trackid2label[trackID] = 1;
      }

    for(int plane =0; plane<3; plane++)
      {
      // larcv::Image2D image(instance_v[plane].meta());
      image_v[plane].paint(0.0);

      }//end plane loop




//
//Now Go Through tracks and label track ends
int ntracks = 0;

for(auto const& mctrack2 : LArData<supera::LArMCTrack_t>())
  {
    ntracks +=1;

    //std::cout << ntracks << std::endl;

    if ((mctrack2.Origin() == 2) && (mctrack2.empty() != 1)) //mctrack2.Origin() == 1
      {
      //std::cout << ntracks << "  Origin: " << mctrack2.Origin()<< std::endl;
      //std::cout << std::endl << "I'm labeling track ends! at points:    " ;
      //istart and iend are 4 vectors of ints, (r, ucol, vcol, ycol)
      std::vector<int> istart ={};
      std::vector<int> iend ={};


      if (doesTrackCrossImageBoundary(mctrack2, image_v[0].meta(), 3200) != -1)
        {

        istart  = getFirstStepPosInsideImage( mctrack2, image_v[0].meta(), 3200, true,  0.15, 2.0);     //ev_trigger->TriggerTime()
        iend    = getFirstStepPosInsideImage( mctrack2, image_v[0].meta(), 3200, false, 0.15, 2.0);    //ev_trigger->TriggerTime()
        }

      ntracks = istart.size();
      ntracks = iend.size();


      if ((istart.empty() != 1) && (iend.empty() != 1))
        {
        int t_offset = -200;
        std::cout << istart[0] << "   " << istart[1] <<  "   " << istart[2] <<  "   " << istart[3] <<std::endl;
        std::cout << iend[0] << "   " << iend[1] <<  "   " << iend[2] <<  "   "  << iend[3] <<std::endl;

        for (int k = 0; k<3; k++)
          {
          //Actual Point
          if ((istart[0]+ t_offset >= 0) && (istart[0]+ t_offset <rows2) && (istart[k+1] >=0) && (istart[k+1] <cols2))
          image_v[k].set_pixel(istart[0]+ t_offset,istart[k+1],3); //Track End
          if ((iend[0]+ t_offset >= 0) && (iend[0]+ t_offset <rows2) && (iend[k+1] >=0) && (iend[k+1] <cols2))
          image_v[k].set_pixel(iend[0]+ t_offset,iend[k+1],3); //Track End

          // Make some more pixels around it labeled that
          int box_dim = 50;

          for (int a = 0-(box_dim/2); a<0+(box_dim/2)+1 ; a++)
            {
            for (int b = 0-(box_dim/2);b<0+(box_dim/2)+1; b++)
              {
              if ((istart[0]+a+ t_offset >= 0) && (istart[0]+a+ t_offset <rows2) && (istart[k+1]+b >=0) && (istart[k+1]+b <cols2))
                {image_v[k].set_pixel(istart[0]+a+ t_offset,istart[k+1]+b,3);}
              if ((iend[0]+a+ t_offset >= 0) && (iend[0]+a+ t_offset <rows2) && (iend[k+1]+b >=0) && (iend[k+1]+b <cols2))
                {image_v[k].set_pixel(iend[0]+a+ t_offset,iend[k+1]+b,3);}
              }
            }

          }
          int numsteps = (int)mctrack2.size();
          for (int n =0; n<numsteps; n++)
              {
              auto const& point = mctrack2.at(n);
              std::vector<double> xyzcoord(3,0);
              xyzcoord = {point.X(), point.Y(),point.Z()};
              std::vector<double> pos4v2 = {point.T(), point.X(), point.Y(),point.Z()};
              supera::ApplySCE(xyzcoord[0],xyzcoord[1],xyzcoord[2]);

              double tick2 = getTick( pos4v2, 3200);
              xyzcoord[0] = (tick2-3200.0)*supera::DriftVelocity()*0.5;

              std::vector<int> wirecoord(3,0);
              wirecoord = getProjectedImagePixel(xyzcoord, image_v[0].meta(), 3, 0);
              for (int k = 0; k<3; k++)
                {
                //Actual Point
                if ((wirecoord[0]+ t_offset >= 0) && (wirecoord[0]+ t_offset <rows2) && (wirecoord[k+1] >=0) && (wirecoord[k+1] <cols2))
                  {image_v[k].set_pixel(wirecoord[0]+ t_offset,wirecoord[k+1],3);}
                int box_dim = 25;
                for (int a = 0-(box_dim/2); a<0+(box_dim/2)+1 ; a++)
                  {
                  for (int b = 0-(box_dim/2);b<0+(box_dim/2)+1; b++)
                    {
                    if ((wirecoord[0]+a+ t_offset >= 0) && (wirecoord[0]+a+ t_offset <rows2) && (wirecoord[k+1]+b >=0) && (wirecoord[k+1]+b <cols2))
                      {
                      image_v[k].set_pixel(wirecoord[0]+a+ t_offset,wirecoord[k+1]+b,3);

                      }

                    }
                  }

                }
              }


        }

      }

  }
  std::cout << "1st. Num Tracks = " << ntrack << "     Num Shower = " << nshower << "        Num Background = " << nbackground << std::endl;
  ntrack = 0;
  nshower = 0;
  nbackground = 0;
  int nends =0;
  int pix_label;

  for (int p = 0; p<3 ; p++)
  {

    for (int r2 =0; r2<rows2;r2++)
      {
      for (int c2 =0; c2<cols2;c2++)
        {
        pix_label = static_cast<int>(image_v[p].pixel(r2,c2));
        if (pix_label ==0) nbackground +=1;
        if (pix_label ==1) ntrack +=1;
        if (pix_label ==2) nshower +=1;
        if (pix_label ==3) nends +=1;


        }
      }
  }

    ev_image->emplace(std::move(image_v));

    std::cout << "2nd. Num Tracks = " << ntrack << "     Num Shower = " << nshower << "        Num Background = " << nbackground << "      Num End = " << nends << std::endl;
    std::cout << ntrack+nshower+nbackground+nends << std::endl;

    return true;

  }

  void SuperaSegment::finalize()
  {}


  double SuperaSegment::getTick( const supera::LArMCStep_t& step, const double trig_time) {
    std::vector<double> pos(4,0);
    pos[0] = step.T();
    pos[1] = step.X();
    pos[2] = step.Y();
    pos[3] = step.Z();
    return getTick( pos, trig_time);
  }


  double SuperaSegment::getTick( const std::vector<double>& step, const double trig_time) {
    // Function returns the tick time of a MC step point
    // if SCE pointer is null, we do not correct for the space charge

    std::vector<double> scepos(3,0);
      scepos[0] = step[1];
      scepos[1] = step[2];
      scepos[2] = step[3];
      supera::ApplySCE(scepos[0],scepos[1],scepos[2]);

    const double cm_per_tick =supera::DriftVelocity()*0.5;
    double tick = ( step[0]*1.0e-3 - (trig_time-4050.0) )/0.5 + scepos[0]/cm_per_tick + 3200.0;

    return tick;
  }

  std::vector<int> SuperaSegment::getProjectedImagePixel( const std::vector<double>& pos3d, const larcv::ImageMeta& meta, const int nplanes, const double fracpixborder ) {
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

    std::vector<int> img_coords( nplanes+1, -1 );
    double row_border = fabs(fracpixborder)*meta.pixel_height();
    double col_border = fabs(fracpixborder)*meta.pixel_width();

    // tick/row
    double tick = pos3d[0]/(0.5*supera::DriftVelocity()) + 3200.0;
    if ( tick<meta.min_y() ) {
      if ( tick>meta.min_y()-row_border )
	// below min_y-border, out of image
	img_coords[0] = meta.rows()-1; // note that tick axis and row indicies are in inverse order
      else
	// outside of image and border
	img_coords[0] = -1;
    }
    else if ( tick>meta.max_y() ) {
      if ( tick<meta.max_y()+row_border )
	// within upper border
	img_coords[0] = 0;
      else
	// outside of image and border
	img_coords[0] = -1;
    }
    else {
      // within the image
      img_coords[0] = meta.row( tick );
    }

    // Columns
    Double_t xyz[3] = { pos3d[0], pos3d[1], pos3d[2] };
    // there is a corner where the V plane wire number causes an error
    if ( (pos3d[1]>-117.0 && pos3d[1]<-116.0) && pos3d[2]<2.0 ) {
      std::cout << __PRETTY_FUNCTION__ << ": v-plane corner hack (" << xyz[0] << "," << xyz[1] << "," << xyz[2] << ")" << std::endl;
      xyz[1] = -116.0;
    }

    for (int p=0; p<nplanes; p++) {
      int wire = supera::NearestWire( xyz, p );
      // round wire
      //wire = std::roundf(wire);

      // get image coordinates
      if ( wire<meta.min_x() ) {
	if ( wire>meta.min_x()-col_border ) {
	  std::cout << __PRETTY_FUNCTION__ << " plane=" << p << " wire=" << wire << "<" << meta.min_x()-col_border << std::endl;
	  // within lower border
	  img_coords[p+1] = 0;
	}
	else
	  img_coords[p+1] = -1;
      }
      else if ( wire>=meta.max_x() ) {
	if ( wire<meta.max_x()+col_border ) {
	  std::cout << __PRETTY_FUNCTION__ << " plane=" << p << " wire=" << wire << ">" << meta.max_x()+col_border << std::endl;
	  // within border
	  img_coords[p+1] = meta.cols()-1;
	}
	else
	  // outside border
	  img_coords[p+1] = -1;
      }
      else
	// inside image
	img_coords[p+1] = meta.col( wire );
    }//end of plane loop

    // there is a corner where the V plane wire number causes an error
    if ( pos3d[1]<-116.3 && pos3d[2]<2.0 && img_coords[1+1]==-1 ) {
      img_coords[1+1] = 0;
    }

    return img_coords;
  }

  std::vector<int> SuperaSegment::getFirstStepPosInsideImage( const supera::LArMCTrack_t& track, const larcv::ImageMeta& meta, const double trig_time,
  		const bool startAtstart, const double max_step_size, const double fv_border) { //const double fv_border Another argument to function?
      // This function returns the (SCE-corrected) position where a MC track first is inside the image bounds
      std::cout << "Track ID: " << track.TrackID() << std::endl;
      double dwall_min = 1.0e9;
      std::vector<double> pos_min(3);
      const double cm_per_tick = supera::DriftVelocity()*0.5;
      int npts = (int)track.size();
      //Begin Loop through points of track
      for ( int ipt=1; ipt<npts; ipt++ )
              {

              int thispt = ipt;
              int lastpt = thispt-1;

              if ( !startAtstart )
                   {
          	       thispt = npts-1-ipt;
          	       lastpt = thispt+1;
                   }

              const auto& this_step = track.at( thispt );
              const auto& last_step = track.at( lastpt );

              double dir[3] = { double(this_step.X()-last_step.X()), double(this_step.Y()-last_step.Y()), double(this_step.Z()-last_step.Z()) };
              double dirnorm = 0;

              for (int i=0; i<3; i++)
                    {
            	      dirnorm += dir[i]*dir[i];
                    }

              dirnorm = sqrt(dirnorm);
              if ( dirnorm<1.0e-3 )
        	      continue;

              for (int i=0; i<3; i++)
                    {
        	          dir[i] /= dirnorm;
                    }

              int nsteps=dirnorm/max_step_size+1;
              if ( nsteps<= 0 )
        	           nsteps = 1;
              double stepsize = dirnorm/double(nsteps);
              for (int istep=0; istep<nsteps; istep++)
                     {
            	       std::vector<double> pos(4,0);
            	       std::vector<double> pos4v(4,0);
            	       pos4v[0] = last_step.T();
            	       pos[0] = pos4v[1] = last_step.X() + stepsize*double(istep)*dir[0];
            	       pos[1] = pos4v[2] = last_step.Y() + stepsize*double(istep)*dir[1];
            	       pos[2] = pos4v[3] = last_step.Z() + stepsize*double(istep)*dir[2];



                     // get sce-corrected 3d pos
            	       std::vector<double> pos_sce(3);
            	       pos_sce[0] = pos[0];
            	       pos_sce[1] = pos[1];
            	       pos_sce[2] = pos[2];
                     supera::ApplySCE(pos_sce[0],pos_sce[1],pos_sce[2]);

                     int boundary_type = -1;

          	         double fdwall = dwall(pos_sce, boundary_type); // use apparent distance...
          	         if ( fdwall<fv_border )
                          {
                          if (fdwall < dwall_min)
                            {
                            dwall_min = fdwall;
                            pos_min = pos_sce;
                            std::cout << "Boundary " << boundary_type << " was within the distance cutoff of " << fv_border <<" Away by " << dwall_min << " With coordinates of:   " << pos_min[0]  << "   "<< pos_min[1]  << "   " << pos_min[2] << std::endl;

                            }
                          continue;
                          }
                     // convert to image coordinates
            	       double tick = getTick( pos4v, trig_time);
            	       if ( tick<meta.min_y()+20.0 || tick>meta.max_y()-20.0 )
            	        continue;



                     pos_sce[0] = (tick-3200.0)*cm_per_tick; // we have to give the apparent-x (relative to the trigger) because we need to test the position in the image


                     std::vector<int> imgcoords;
                     imgcoords = getProjectedImagePixel( pos_sce, meta, 3 );

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

            }//end of track point loop

      // didn't find the crossing boundary
      std::cout << "XXXXXXXXX DIDNT FIND BOUNDARY XXXXXXXXXXXXXXX";
      std::vector<int> empty;
      return empty;
    }

    int SuperaSegment::doesTrackCrossImageBoundary( const  supera::LArMCTrack_t& track, const larcv::ImageMeta& meta, const double trig_time ) {
     double tick_start = getTick( track.front(), trig_time );
     double tick_end   = getTick( track.back(), trig_time );
     if ( tick_start>meta.min_y() && tick_start<meta.max_y() && tick_end>meta.min_y() && tick_end<meta.max_y() )
       return -1;

     if ( tick_start < meta.min_y() && tick_end > meta.min_y() )
       return 0; // start out -> end in
     else if ( tick_start > meta.min_y() && tick_end < meta.min_y())
       return 1; // start in -> end out
     else if ( tick_start < meta.max_y() && tick_end > meta.max_y())
       return 1; // start in -> end out;
     else if ( tick_start > meta.max_y() && tick_end < meta.max_y() )
       return 0; // start out -> end in

     return -1;
   }

   double SuperaSegment::dwall( const std::vector<double>& pos, int& boundary_type ) {
  std::vector<float> fpos(3);
  for (int i=0; i<3; i++)
    fpos[i] = (float)pos[i];
  return dwall( fpos, boundary_type);
}

float SuperaSegment::dwall( const std::vector<float>& pos, int& boundary_type ) {

  float dx1 = pos[0];
  float dx2 = 258-pos[0];
  float dy1 = 117.0-pos[1];
  float dy2 = pos[1]+117.0;
  float dz1 = pos[2];
  float dz2 = 1036.0-pos[2];

  float dwall = 1.0e9;

  if ( dy1<dwall ) {
    dwall = dy1;
    boundary_type = 0; // top
  }
  if ( dy2<dwall ) {
    dwall = dy2;
    boundary_type = 1; // bottom
  }
  if ( dz1<dwall ) {
    dwall = dz1;
    boundary_type = 2; // upstream
  }
  if ( dz2<dwall ) {
    dwall = dz2;
    boundary_type = 3; // downstream
  }
  if ( dx1<dwall ) {
    dwall = dx1;
    boundary_type = 4; // anode
  }
  if ( dx2<dwall ) {
    dwall = dx2;
    boundary_type = 5; // cathode
  }

  return dwall;

}


} //End of namespace
#endif
