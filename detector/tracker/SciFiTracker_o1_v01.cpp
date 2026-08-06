#include "Core/UpgradeTags.h"
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DD4hep/detail/DetectorInterna.h"
#include "TClass.h"
#include "XML/Utilities.h"
#include "DD4hep/OpticalSurfaces.h"
#include <cmath>

using namespace dd4hep;

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
    double station_gap    = 0.0;
    int    n_layers        = 0;
    double barrel_radius   = 0.0;


    struct Module {
      double      sign = 1.0, leftHoleSizeY = 0, rightHoleSizeY = 0;
      std::string left, right;
      Module( double sgn, double hl, double hr, const std::string& l, const std::string& r )
          : sign( sgn ), leftHoleSizeY( hl ), rightHoleSizeY( hr ), left( l ), right( r ) {}
      HELPER_CLASS_DEFAULTS( Module );
    };
    std::map<std::string, Module>     m_modules;
    std::map<std::string, DetElement> m_detElements;
    std::map<std::string, double>     HoleType;

    /// Utility function to register detector elements for cloning. Will all be deleted!
    void registerDetElement( const std::string& nam, DetElement de );
    /// Utility function to access detector elements for cloning. Will all be deleted!
    dd4hep::DetElement detElement( const std::string& nam ) const;

    /// Initializing constructor
    SciFiBuild( dd4hep::Detector& description, xml_elt_t e, dd4hep::SensitiveDetector sens );
    SciFiBuild()                              = delete;
    SciFiBuild( SciFiBuild&& )                 = delete;
    SciFiBuild( const SciFiBuild& )            = delete;
    SciFiBuild& operator=( const SciFiBuild& ) = delete;

    /// Default destructor
    virtual ~SciFiBuild();
    void build_fiber( dd4hep::Volume motherVol, dd4hep::Position pos );
    void build_mat_slab2( dd4hep::Volume motherVol );
    void build_layers();
    void build_stations();
    void build_detector();
  };

  void SciFiBuild::registerDetElement( const std::string& nam, dd4hep::DetElement de ) { m_detElements[nam] = de; }
  dd4hep::DetElement SciFiBuild::detElement( const std::string& nam ) const {
    auto i = m_detElements.find( nam );
    if ( i == m_detElements.end() ) {
      dd4hep::printout( dd4hep::ERROR, "SciFi-geo", "Attempt to access non-existing detector element %s", nam.c_str() );
    }
    return ( *i ).second;
  }

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
    station_gap       = dd4hep::_toDouble( "SciFi:StationGap" );
    n_layers          = dd4hep::_toInt( "SciFi:nLayers" );
    barrel_radius     = dd4hep::_toDouble( "SciFi:barrelRadius" );

  }

  SciFiBuild::~SciFiBuild() {
    for ( auto& d : m_detElements ) dd4hep::detail::destroyHandle( d.second );
  }

  void SciFiBuild::build_fiber( dd4hep::Volume motherVol, dd4hep::Position pos ) {
    // Places one fiber (core/clad1/clad2) into motherVol at pos.

    static dd4hep::Volume         coreVol, clad1Vol, clad2Vol;
    static dd4hep::OpticalSurface coreCladSurf, cladCladSurf;

    if ( !coreVol.isValid() ) {
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

    // Place the mat volume in the mother volume
    motherVol.placeVolume( matVol, tr );

  }

  void SciFiBuild::build_detector() {
    dd4hep::PlacedVolume pv;
    dd4hep::Assembly     scifi_vol( "lvSciFi" );

    scifi_vol.setVisAttributes( description.invisible() );
    sensitive.setType( "tracker" );
    load( x_det, "include" );

    // Detailed per-fiber model, replacing the old homogenized build_mats().
    build_mat_slab2( scifi_vol );

    xml_h x_tr = x_det.child( _U( transformation ) );
    pv         = placeDetector( scifi_vol, x_tr );
    pv.addPhysVolID( "system", id );
  }
} // namespace

static dd4hep::Ref_t create_element( dd4hep::Detector& description, xml_h e, dd4hep::SensitiveDetector sens_det ) {
  SciFiBuild builder( description, e, sens_det );
  builder.build_detector();
  return builder.detector;
}
DECLARE_DETELEMENT( LHCb_SciFi_geo_v1_0, create_element )

