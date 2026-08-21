
//#include "Core/UpgradeTags.h"
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DD4hep/detail/DetectorInterna.h"
#include "TClass.h"
#include "XML/Utilities.h"
#include "DD4hep/OpticalSurfaces.h"
#include <cmath>


using namespace std;

#define HELPER_CLASS_DEFAULTS( x )                                                                                     \
  x()                      = default;                                                                                  \
  x( x&& )                 = default;                                                                                  \
  x( const x& )            = default;                                                                                  \
  x& operator=( const x& ) = default;                                                                                  \
  ~x()                     = default

namespace {

  /// Helper class to build the SciFi tracker geometry
  struct SciFiBuild : public dd4hep::xml::tools::VolumeBuilder {
    // Debug flags: To be enables/disbaled by the <debug/> section in the detector.xml file.
    bool        m_build_layers      = false;
    bool        m_build_stations    = false;


    // General constants used by several member functions to build the SciFi geometry

    double fiber_radius        = 0.0;
    double core_frac = 0.0;
    double clad1_frac = 0.0;
    double clad2_frac = 0.0;
    double mat_width        = 0.0;
    double fiber_length        = 0.0;
    int    n_layers        = 0;
    double barrel_radius   = 0.0;

    double strip_width         = 0.0;
    double strip_height        = 0.0;
    double strip_thickness     = 0.0;
    double pixel_thickness     = 0.0;
    double air_gap_sensitive   = 0.0;
    int npixel_strip_width  = 0;
    int npixel_strip_height = 0;
    double pixel_sensitive_ratio = 0.0;



    struct Module {
      double      sign = 1.0, leftHoleSizeY = 0, rightHoleSizeY = 0;
      std::string left, right;
      Module( double sgn, double hl, double hr, const std::string& l, const std::string& r )
          : sign( sgn ), leftHoleSizeY( hl ), rightHoleSizeY( hr ), left( l ), right( r ) {}
      HELPER_CLASS_DEFAULTS( Module );
    };
    std::map<std::string, Module>     m_modules;
    // Unused: nothing calls registerDetElement() anymore 
    // Kept commented in case a future cloning-style build
    // (e.g. repeated stations) needs it.
    // std::map<std::string, DetElement> m_detElements;
    std::map<std::string, double>     HoleType;

    // /// Utility function to register detector elements for cloning. Will all be deleted!
    // void registerDetElement( const std::string& nam, DetElement de );
    // /// Utility function to access detector elements for cloning. Will all be deleted!
    // dd4hep::DetElement detElement( const std::string& nam ) const;

    /// Initializing constructor
    SciFiBuild( dd4hep::Detector& description, xml_elt_t e, dd4hep::SensitiveDetector sens );
    SciFiBuild()                              = delete;
    SciFiBuild( SciFiBuild&& )                 = delete;
    SciFiBuild( const SciFiBuild& )            = delete;
    SciFiBuild& operator=( const SciFiBuild& ) = delete;

    /// Default destructor
    virtual ~SciFiBuild();
    void build_fiber( dd4hep::Volume motherVol, dd4hep::Position pos );
    dd4hep::Volume build_sipm();
    void build_mat_slab2( dd4hep::Volume motherVol );
    void build_layers();
    void build_stations();
    void build_detector();
  };

  // void SciFiBuild::registerDetElement( const std::string& nam, dd4hep::DetElement de ) { m_detElements[nam] = de; }
  // dd4hep::DetElement SciFiBuild::detElement( const std::string& nam ) const {
  //   auto i = m_detElements.find( nam );
  //   if ( i == m_detElements.end() ) {
  //     dd4hep::printout( dd4hep::ERROR, "SciFi-geo", "Attempt to access non-existing detector element %s", nam.c_str() );
  //   }
  //   return ( *i ).second;
  // }

  /// Initializing constructor
  SciFiBuild::SciFiBuild( dd4hep::Detector& dsc, xml_elt_t e, dd4hep::SensitiveDetector sens )
      : dd4hep::xml::tools::VolumeBuilder( dsc, e, sens ) {
    // Process debug flags
    xml_comp_t x_dbg = x_det.child( _U( debug ), false );
    if ( x_dbg ) {
      for ( xml_coll_t i( x_dbg, _U( item ) ); i; ++i ) {
        xml_comp_t  c( i );
        std::string n = c.nameStr();
        if ( n == "build_layers" ) m_build_layers = c.attr<bool>( _U( value ) );
        if ( n == "build_stations" ) m_build_stations = c.attr<bool>( _U( value ) );
      }
    }
    
    fiber_radius               = dd4hep::_toDouble( "SciFi:fibreRadius" );
    core_frac       = dd4hep::_toDouble( "SciFi:coreFraction" );
    clad1_frac        = dd4hep::_toDouble( "SciFi:clad1Fraction" );
    clad2_frac        = dd4hep::_toDouble( "SciFi:clad2Fraction" );
    mat_width       = dd4hep::_toDouble( "SciFi:matWidth" );
    fiber_length        = dd4hep::_toDouble( "SciFi:fibreLength" );
    n_layers          = dd4hep::_toInt( "SciFi:nLayers" );
    barrel_radius     = dd4hep::_toDouble( "SciFi:barrelRadius" );
    strip_width         = dd4hep::_toDouble( "SciFi:stripWidth" );
    strip_height        = dd4hep::_toDouble( "SciFi:stripHeight" );
    strip_thickness     = dd4hep::_toDouble( "SciFi:stripThickness" );
    pixel_thickness     = dd4hep::_toDouble( "SciFi:pixelThickness" );
    air_gap_sensitive   = dd4hep::_toDouble( "SciFi:airGapSensitive" );
    npixel_strip_width  = dd4hep::_toInt( "SciFi:npixelPerStripWidth" );
    npixel_strip_height = dd4hep::_toInt( "SciFi:npixelPerStripHeight" );
    pixel_sensitive_ratio = dd4hep::_toDouble( "SciFi:pixelSensitiveRatio" );


  }

  SciFiBuild::~SciFiBuild() {
    // for ( auto& d : m_detElements ) dd4hep::detail::destroyHandle( d.second );
  }

  void SciFiBuild::build_fiber( dd4hep::Volume motherVol, dd4hep::Position pos ) {
    // Places one fiber (core/clad1/clad2) into motherVol at pos.

    static dd4hep::Volume         coreVol, clad1Vol, clad2Vol;
    static dd4hep::OpticalSurface coreCladSurf, cladCladSurf;

    if ( !coreVol.isValid() ) {  //this block is executed only once 
      dd4hep::Material fibre_core  = description.material( "ScintCoreMaterial" );
      dd4hep::Material fibre_clad1 = description.material( "InnerCladdingMaterial" );
      dd4hep::Material fibre_clad2 = description.material( "OuterCladdingMaterial" );

      double half_len = fiber_length / 2.0;
      double core_r   = fiber_radius * core_frac;
      double clad1_r  = fiber_radius * clad1_frac;
      double clad2_r  = fiber_radius * clad2_frac;

      dd4hep::Tube coreSolid( 0.0, core_r, half_len );    //Equivalent of G4Solid and G4Tube
      dd4hep::Tube clad1Solid( core_r, clad1_r, half_len );
      dd4hep::Tube clad2Solid( clad1_r, clad2_r, half_len );

      coreVol  = dd4hep::Volume( "lvSciFiFiberCore", coreSolid, fibre_core ); //Equivalent of G4LogicalVolume
      clad1Vol = dd4hep::Volume( "lvSciFiFiberClad1", clad1Solid, fibre_clad1 );
      clad2Vol = dd4hep::Volume( "lvSciFiFiberClad2", clad2Solid, fibre_clad2 );

      coreVol.setVisAttributes( description, "SciFi:FiberCoreVis" );
      clad1Vol.setVisAttributes( description, "SciFi:FiberClad1Vis" );
      clad2Vol.setVisAttributes( description, "SciFi:FiberClad2Vis" );

      dd4hep::OpticalSurfaceManager surfMgr   = description.surfaceManager();
      coreCladSurf                            = surfMgr.opticalSurface( "Core_InnerClad_Surface" );
      cladCladSurf                            = surfMgr.opticalSurface( "InnerClad_OuterClad_Surface" );
      dd4hep::OpticalSurface outerSurf        = surfMgr.opticalSurface( "OuterClad_Surface" );


      dd4hep::SkinSurface( description, detector, "clad2_outer_surface", outerSurf, clad2Vol );

      dd4hep::printout( dd4hep::DEBUG, "SciFi-geo",
                        "+++ Built (once) fiber cross-section: core r=%.4f clad1 r=%.4f clad2 r=%.4f (mm)", core_r,
                        clad1_r, clad2_r );
    }

    dd4hep::PlacedVolume corePV  = motherVol.placeVolume( coreVol, pos );  //Equivalent of G4PVPlacement
    dd4hep::PlacedVolume clad1PV = motherVol.placeVolume( clad1Vol, pos );
    dd4hep::PlacedVolume clad2PV = motherVol.placeVolume( clad2Vol, pos );

    // Unique per call, since every fiber placement needs its own DetElement
    static int fiberIdx = 0;
    ++fiberIdx;
    std::string fName = "SciFiFiber_" + std::to_string( fiberIdx );

    dd4hep::DetElement coreDE( detector, fName + "CoreDE", 2 * fiberIdx );
    dd4hep::DetElement clad1DE( detector, fName + "Clad1DE", 2 * fiberIdx + 1 );
    coreDE.setPlacement( corePV );
    clad1DE.setPlacement( clad1PV );

    // BorderSurface on the two internal interfaces to account for roughness
    dd4hep::BorderSurface( description, coreDE, fName + "_core_to_clad1", coreCladSurf, corePV, clad1PV );
    dd4hep::BorderSurface( description, coreDE, fName + "_clad1_to_core", coreCladSurf, clad1PV, corePV );
    dd4hep::BorderSurface( description, clad1DE, fName + "_clad1_to_clad2", cladCladSurf, clad1PV, clad2PV );
    dd4hep::BorderSurface( description, clad1DE, fName + "_clad2_to_clad1", cladCladSurf, clad2PV, clad1PV );
  }

  dd4hep::Volume SciFiBuild::build_sipm() {
    static dd4hep::Volume sipmVol, stripVol, pixelVol, sensitiveVol;
    static int n_strips_x = 0;

    if ( !sipmVol.isValid()) {  //this block is executed only once

      n_strips_x = (int)std::floor( mat_width / strip_width );
      double sipm_width  = strip_width * n_strips_x;
      double sipm_half_z = ( pixel_thickness + strip_thickness + air_gap_sensitive ) / 2.0;

      // Outer container volume for the SiPMs
      dd4hep::Box    sipmSolid( sipm_width/2.0, (strip_height*1.2)/2.0, sipm_half_z );
      sipmVol = dd4hep::Volume( "lvSciFiSiPMArray", sipmSolid, description.air() );
      sipmVol.setVisAttributes( description, "SciFi:SiPMVis" );

      // Solid + volume for one strip
      dd4hep::Box    stripSolid( strip_width/2.0, strip_height/2.0, strip_thickness/2.0 );
      stripVol = dd4hep::Volume( "lvSciFiSiPMStrip", stripSolid, description.material("Epoxy") );
      stripVol.setVisAttributes( description.invisible() );

      // Solid + volume for one pixel
      double pixel_dim_x = strip_width  / npixel_strip_width;
      double pixel_dim_y = strip_height / npixel_strip_height;
      dd4hep::Box pixelBase( pixel_dim_x/2.0, pixel_dim_y/2.0, pixel_thickness/2.0 );
      dd4hep::Box sensitiveBase( pixel_dim_x*std::sqrt(pixel_sensitive_ratio)/2.0,
                                pixel_dim_y*std::sqrt(pixel_sensitive_ratio)/2.0, pixel_thickness/2.0 );
      dd4hep::SubtractionSolid pixelSolid( pixelBase, sensitiveBase );  // "dead" pixel shell

      pixelVol = dd4hep::Volume( "lvSciFiPixel", pixelSolid, description.material("Epoxy") );
      pixelVol.setVisAttributes( description.invisible() );

      // Sensitive pixels to produce hits

      sensitiveVol = dd4hep::Volume( "lvSciFiPixelSensitive", sensitiveBase, description.material("Epoxy") );
      sensitiveVol.setVisAttributes( description.invisible() );
      sensitiveVol.setSensitiveDetector( sensitive );
      
      // Vertical alignment of the strip/pixel row against the fiber stack

      double y_dist_sipm = ( (0.210/0.250) * 2.0 * fiber_radius ) / 2.0;  //the local fiber-layer pitch used only inside the SiPM-placement 
      double y_strip    = y_dist_sipm * ( n_layers / 2.0 ) - fiber_radius / 2.0 - y_dist_sipm * (double)( n_layers / 2 ); //vertical alignment of the strip row
      double y_pix_base = y_dist_sipm * ( n_layers / 2.0 ) - fiber_radius + pixel_dim_y - y_dist_sipm * (double)( n_layers / 2 ); //vertical alignment of the pixel grid

      double pz = -sipm_half_z + strip_thickness + pixel_thickness/2.0 + air_gap_sensitive;

      dd4hep::Box    stripContainerSolid( strip_width/2.0, (strip_height*1.2)/2.0, pixel_thickness/2.0 );
      dd4hep::Volume stripContainer( "lvSciFiStripContainer", stripContainerSolid, description.air() );
      stripContainer.setVisAttributes( description.invisible() );

      for ( int i = 0; i < npixel_strip_width; ++i ) {
        for ( int j = 0; j < npixel_strip_height; ++j ) {
          double px = -strip_width/2.0 + pixel_dim_x/2.0 + i*pixel_dim_x;   // local to the strip now
          double py = y_pix_base - strip_height/2.0 + pixel_dim_y/2.0 + j*pixel_dim_y;
          dd4hep::Position pixPos( px, py, 0.0 );
          dd4hep::PlacedVolume pixelPV = stripContainer.placeVolume( pixelVol, pixPos );
          int sipm_copy_no = i + npixel_strip_width*j;
          dd4hep::PlacedVolume pixelsensitiveareaPV = stripContainer.placeVolume( sensitiveVol, sipm_copy_no, pixPos );
          pixelsensitiveareaPV.addPhysVolID( "pixelX", i );
          pixelsensitiveareaPV.addPhysVolID( "pixelY", j );
        }
      }

      // Strip placement

      double x_offset = strip_width * (n_strips_x - 1) / 2.0;
      for ( int k = 0; k < n_strips_x; ++k ) {
        double x_strip = strip_width*k - x_offset;
        dd4hep::Position stripPos( x_strip, y_strip, strip_thickness/2.0 + air_gap_sensitive - sipm_half_z );
        dd4hep::PlacedVolume stripPV = sipmVol.placeVolume( stripVol, stripPos );

        dd4hep::Position stripContainerPos( x_strip, 0.0, pz );
        dd4hep::PlacedVolume stripContainerPV = sipmVol.placeVolume( stripContainer, stripContainerPos );
        stripContainerPV.addPhysVolID( "strip", k );
      }

    }
    return sipmVol;
  }

  void SciFiBuild::build_mat_slab2( dd4hep::Volume motherVol ) {

    // Mat placement
    double               phi0 = 0.0;
    dd4hep::RotationZ    rot( phi0 + M_PI/2.0 );
    dd4hep::Position      pos( barrel_radius*cos(phi0), barrel_radius*sin(phi0), 0.0 );
    dd4hep::Transform3D  tr( rot, pos );

    // Mat volume

    double mat_halfw = mat_width / 2.0;
    // Distance between fibers centers (horizontal and vertical)
    double x_dist       = (0.275/0.250) * 2.0 * fiber_radius;
    double y_dist       = (0.210/0.250) * 2.0 * fiber_radius;
    double mat_half_h   = 0.5*( n_layers*y_dist + 2.0*fiber_radius ) + 0.2*dd4hep::mm;
    double mat_half_len = fiber_length / 2.0;

    dd4hep::Material matMat  = description.material( "EpoxyTi02" );
    dd4hep::Box      matSolid( mat_halfw, mat_half_h, mat_half_len );
    dd4hep::Volume   matVol( "lvSciFiMatSlab2", matSolid, matMat );
    matVol.setVisAttributes( description, "SciFi:MatVis" );

    // Fiber placement in the mat volume

    int    n_fibers_x   = (int)std::floor( 2.0*(mat_halfw - fiber_radius/2.0) / x_dist );
    int    n_fibers_y   = n_layers;
    double center_idx_x = (n_fibers_x - 1.0) / 2.0;
    double center_idx_y = (n_fibers_y - 1.0) / 2.0;


    for (int j = 0; j < n_fibers_y; ++j) {
      double x_displ = (j % 2 == 0) ? -x_dist/4.0 : x_dist/4.0;
      double y_pos   = (center_idx_y - j) * y_dist;
      for (int k = 0; k < n_fibers_x; ++k) {
        double x_pos = (k - center_idx_x) * x_dist + x_displ;
        build_fiber( matVol, dd4hep::Position(x_pos, y_pos, 0.0) );
      }
    }

    // SiPM arrays at the two ends of the mat

    double sipm_half_z = ( pixel_thickness + strip_thickness + air_gap_sensitive ) / 2.0;
    dd4hep::Volume sipmVol = build_sipm();

    dd4hep::Transform3D tr_ep1 = tr * dd4hep::Transform3D( dd4hep::RotationY(M_PI),
                                  dd4hep::Position(0, 0, -(mat_half_len + sipm_half_z)) );  // 180deg-flipped end
    dd4hep::Transform3D tr_ep2 = tr * dd4hep::Transform3D( dd4hep::RotationY(0.0),
                                  dd4hep::Position(0, 0,  (mat_half_len + sipm_half_z)) );  // unflipped end

    dd4hep::PlacedVolume sipm1PV = motherVol.placeVolume( sipmVol, tr_ep1 );
    sipm1PV.addPhysVolID( "side", 0 ); //back readout
    dd4hep::PlacedVolume sipm2PV = motherVol.placeVolume( sipmVol, tr_ep2 );
    sipm2PV.addPhysVolID( "side", 1 ); //front readout


    // Place the mat volume in the mother volume
    motherVol.placeVolume( matVol, tr );

  }

  void SciFiBuild::build_detector() {
    dd4hep::PlacedVolume pv;
    dd4hep::Assembly     scifi_vol( "lvSciFi" );

    scifi_vol.setVisAttributes( description.invisible() );
    sensitive.setType( "tracker" ); //the G4 sensitive detector class used will be Geant4SimpleTrackerAction
    load( x_det, "include" );

    // Detailed per-fiber model
    build_mat_slab2( scifi_vol );

    //xml_h x_tr = x_det.child( _U( transformation ) );
    pv = placeDetector( scifi_vol);
    pv.addPhysVolID( "system", id ); 
    //"system" is the name of the field in the xml file, "id" is the bitfield value of the field
    //In this case the ID of the SciFiBarrel is set to 3 in the xml file of ALLEGRO
  }
} // namespace

static dd4hep::Ref_t create_element( dd4hep::Detector& description, xml_h e, dd4hep::SensitiveDetector sens_det ) {
  SciFiBuild builder( description, e, sens_det );
  builder.build_detector();
  return builder.detector;
}
DECLARE_DETELEMENT( LHCb_SciFi_geo_v1_0, create_element )

