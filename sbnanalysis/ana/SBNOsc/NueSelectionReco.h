#ifndef __sbnanalysis_ana_SBNOsc_NueSelectionReco__
#define __sbnanalysis_ana_SBNOsc_NueSelectionReco__

/**
 * \file NueSelectionReco.h
 *
 * SBN nue selection.
 *
 * Author:
 */

//C++ Includes 
#include <iostream>
#include <vector> 
#include <ctime> 
#include <map>

//Framework Includes 
#include "canvas/Utilities/InputTag.h"

//SBN Includes 
#include "core/SelectionBase.hh"
#include "core/Event.hh"
#include "core/PostProcessorBase.hh"
#include "core/ProviderManager.hh"

//Larsoft Includes 
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataalg/DetectorInfo//DetectorPropertiesStandard.h"

// take the geobox stuff from uboonecode                          
#include "ubcore/LLBasicTool/GeoAlgo/GeoAABox.h"
#include "ubcore/LLBasicTool/GeoAlgo/GeoAlgo.h"

//Root Includes
#include "TRandom3.h"
#include "THStack.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Shower.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOne.h"


class TH2D;

namespace ana {
  namespace SBNOsc {

/**
 * \class NueSelectionReco
 * \brief Electron neutrino event selection
 */
class NueSelectionReco : public core::SelectionBase {
public:
  /** Constructor. */
  NueSelectionReco();

  /**
   * Initialization.
   *
   * \param config A configuration, as a FHiCL ParameterSet object
   */
  void Initialize(fhicl::ParameterSet* config=NULL);

  /** Finalize and write objects to the output file. */
  void Finalize();

  /**
   * Process one event.
   *
   * \param ev A single event, as a gallery::Event
   * \param Reconstructed interactions
   * \return True to keep event
   */
  bool ProcessEvent(const gallery::Event& ev, const std::vector<event::Interaction> &truth, std::vector<event::RecoInteraction>& reco);

protected:

  unsigned EventCounter;  //!< Count processed events
  unsigned NuCount;  //!< Count selected events

  /** Configuration parameters */
  struct Config {
    //    art::InputTag fTruthTag;  //!< art tag for MCTruth information
    std::vector<geoalgo::AABox> fiducial_volumes; //!< List of FV containers -- set by "fiducial_volumes"
    std::vector<geoalgo::AABox> active_volumes; //!< List of active volumes

    bool ApplyFVCut;
    bool ApplyOneShowerEnergyCut;
    bool ApplyShowerResidualCut;
    bool ApplyConversionGapCut;
    bool ApplydEdxCut;
    bool ApplyLengthCut;
    bool ApplyOpeningAngleCut;
    bool ApplyShowerEnergyCut;
    bool ApplyShowerDensityGradientCut;
    bool ApplyShowerDensityPowerCut;
    bool ApplyShowerTrackLengthCut;
    bool ApplyShowerTrackWidthCut;
    bool ApplyMaxTrackLengthCut;
    bool ApplyMaxTrackPIDACut;
    bool ApplyMaxTrackLengthPIDACut;
    bool ApplyMVACut;
    bool ApplyNeutrinoPdgCodeCut;
    bool ApplyNumneutrinosCut;
     
    bool IncludeCosmics;
    bool IncludeDirt;
    bool CosmicsOnly;
    bool DirtOnly;
    bool Verbose;
    bool IntrinsicOnly;
    bool NuMuOnly;
    bool OscOnly; 
    bool NoDirt;
    bool FillHistograms;
    bool FillModeHistograms;
    bool FillIntTypeHistograms;
    bool UseAllCosmics;
    bool DontUseSimChannels;
    bool UseGenieHists;
    
    float SecondaryShowerEnergyCut;
    float ShowerResidualEnergyCut;
    float ResiudalCut;
    float NumShowersCut;
    float dEdxCut;
    float ConversionGapCut;
    float LengthCut;
    float OpeningAngleCut;
    float ShowerEnergyCut;
    float ShowerDensityGradientCut;
    float ShowerDensityPowerCut;
    float ShowerTrackLengthCut;
    float ShowerTrackWidthCut;
    float MaxTrackLengthCut;
    float MaxTrackPIDACut;
    float MaxTrackLengthPIDACut;
    float MVACut;
    float NeutrinoPdgCodeCut;
    float NumneutrinosCut;
    float XOffset;

    float GlobalWeight;
    float POTWeight;
    float ElectronWeight;
    float PhotonWeight;

    float GlobalTimeOffset;
    float SpillTime;
    float ReadoutWindowSize;

    int NSegments;
    std::vector<float> EnergyConversion;

    float TrackEnergyCorrection;
    float ShowerEnergyCorrection;
    float NutrinoEnergyCorrection;

    bool    SetupMVA;
    TString Weightfile;
    TString MVAMethod;
    std::vector<std::string> MVAValuesNames;

    std::string GeneratorTag;
    std::string PandoraTag;
    std::string HitTag;
    std::string ShowerTag;
    std::string FluxTag;
    std::string TrackTag;
    std::string TrackPIDTag;
    std::string CalorimetryTag;
  };

  //Additional information used by the selection per neutrino interaction/ 
  //C++ is smart and doesn't need the constructor.
  struct NueInteraction {
    double hadronic_energy = 0;
    double shower_energy   = 0;
    double leptonic_energy = 0;
    double other_energy    = 0; 
    double weight;
    int initnu;
    int fnd;
    double True_energy = 0;
    const art::Ptr< simb::MCTruth > mctruth;
    bool booldirtevent = false;
    double GetNueEnergy(){return hadronic_energy + shower_energy + leptonic_energy + other_energy;}
    double GetRecoEnergy(){return hadronic_energy + shower_energy;}

    
  };
  
  //Output Root Histogram map for background selections
  struct RootHistograms{

    double emin = 0.0, emax = 5.0;
    int ebins = 200;

    std::map<std::string,TH1D*> VisibleEnergy_MCTruth_Hist;
    std::map<std::string,TH1D*> TrueNumber_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_BeforeSel_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_AfterSel_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_AfterSelOne_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_AfterSelExtra_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_FVRemoved_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_FVPassed_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_GtrOneShowerRemoved_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_GtrOneShowerPassed_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_OneShowerECutRemoved_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_OneShowerECutPassed_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_OneShowerResidualRemoved_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_OneShowerResidualPassed_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_ConversionGapRemoved_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_ConversionGapPassed_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_dEdxRemoved_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_dEdxPassed_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_LengthRemoved_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_LengthPassed_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_OpeningAngleRemoved_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_OpeningAnglePassed_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_ShowerEnergyRemoved_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_ShowerEnergyPassed_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_ShowerDensityGradientRemoved_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_ShowerDensityGradientPassed_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_ShowerDensityPowerRemoved_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_ShowerDensityPowerPassed_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_ShowerTrackLengthRemoved_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_ShowerTrackLengthPassed_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_ShowerTrackWidthRemoved_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_ShowerTrackWidthPassed_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_MaxTrackLengthRemoved_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_MaxTrackLengthPassed_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_MaxTrackPIDARemoved_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_MaxTrackPIDAPassed_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_MaxTrackLengthPIDARemoved_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_MaxTrackLengthPIDAPassed_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_MVARemoved_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_MVAPassed_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_NeutrinoPdgCodeRemoved_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_NeutrinoPdgCodePassed_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_NumneutrinosRemoved_Hist;
    std::map<std::string,TH1D*> VisibleEnergy_NumneutrinosPassed_Hist;

    TH1D* XDiff;
    TH1D* YDiff; 
    TH1D* ZDiff;
    std::vector<std::string> HistTypes = {"NuMu","InNuE","OscNuE","NCInNuE","NCOscNuE","NCNuMu","DirtNuMu","DirtInNuE",
                                          "DirtOscNuE","DirtNCInNuE","DirtNCOscNuE","DirtNCNuMu","AllSignal","AllBackground"};
  };

  Config fConfig; //!< The config
  RootHistograms fRootHists;


  bool Select(const gallery::Event& ev,
	      const art::Ptr<recob::PFParticle>& neutrino,
	      std::map<int, art::Ptr<recob::PFParticle> >& pfp_map,
	      NueSelectionReco::NueInteraction& intInfo);

  void CalculateEnergy(const gallery::Event& ev,
		       const art::Ptr<recob::PFParticle>& neutrino,
		       std::map<int, art::Ptr<recob::PFParticle> >& pfp_map,
		       NueSelectionReco::NueInteraction& intInfo);


  //Selection Functions 
  bool passFV(const TVector3 &v) {if(fConfig.ApplyFVCut){return containedInFV(v);} else return true;}
  bool passOneShowerEnergyCut(std::vector<art::Ptr<recob::Shower> > neutrino_showers){
    if(fConfig.ApplyOneShowerEnergyCut){return OneShowerEnergyCut(neutrino_showers);}
    return true; 
  }
  bool passShowerResiudalCut(std::vector<art::Ptr<recob::Shower> > const& neutrino_showers,double& numshowers, const TVector3& vertex_position){
    return OneShowerResiudalCut(neutrino_showers,numshowers,vertex_position);
  }
  bool passNeutrinoPdgCodeCut(const double& NeutrinoPdgCode) const {
    if(fConfig.ApplyNeutrinoPdgCodeCut){return NeutrinoPdgCodeCut(NeutrinoPdgCode);}
    return true;
  }
  bool passNumneutrinosCut(const double& Numneutrinos) const {
    if(fConfig.ApplyNumneutrinosCut){return NumneutrinosCut(Numneutrinos);}
    return true;
  }

  bool passConversionGapCut(const TVector3& ShowerStart, double& conversion_gap,const TVector3& vertex_position){
    if(fConfig.ApplyConversionGapCut){return ConversionGapCut(ShowerStart,conversion_gap,vertex_position);}
    ConversionGapCut(ShowerStart,conversion_gap,vertex_position);
    return true;
  }
  bool passdEdxCut(const double& dEdx) const {
    if(fConfig.ApplydEdxCut){return dEdxCut(dEdx);}
    return true;
  }
  bool passLengthCut(const double& Length) const{
    if(fConfig.ApplyLengthCut){return LengthCut(Length);}
    return true;
  }
  bool passOpeningAngleCut(const double& OpeningAngle) const{
    if(fConfig.ApplyOpeningAngleCut){return OpeningAngleCut(OpeningAngle);}
    return true;
  }
  bool passShowerEnergyCut(const double& ShowerEnergy) const{
    if(fConfig.ApplyShowerEnergyCut){return ShowerEnergyCut(ShowerEnergy);}
    return true;
  }
  bool passShowerDensityGradientCut(const double& ShowerDensityGradient) const{
    if(fConfig.ApplyShowerDensityGradientCut){return ShowerDensityGradientCut(ShowerDensityGradient);}
    return true;
  }
  bool passShowerDensityPowerCut(const double& ShowerDensityPower) const{
    if(fConfig.ApplyShowerDensityPowerCut){return ShowerDensityPowerCut(ShowerDensityPower);}
    return true;
  }
  bool passShowerTrackLengthCut(const double& ShowerTrackLength) const{
    if(fConfig.ApplyShowerTrackLengthCut){return ShowerTrackLengthCut(ShowerTrackLength);}
    return true;
  }
  bool passShowerTrackWidthCut(const double& ShowerTrackWidth) const{
    if(fConfig.ApplyShowerTrackWidthCut){return ShowerTrackWidthCut(ShowerTrackWidth);}
    return true;
  }
  bool passMaxTrackLengthCut(const double& MaxTrackLength) const{
    if(fConfig.ApplyMaxTrackLengthCut){return MaxTrackLengthCut(MaxTrackLength);}
    return true;
  }
  bool passMaxTrackPIDACut(const double& MaxTrackPIDA) const{
    if(fConfig.ApplyMaxTrackPIDACut){return MaxTrackPIDACut(MaxTrackPIDA);}
    return true;
  }
  bool passMaxTrackLengthPIDACut(const double& MaxTrackLength, const double& MaxTrackPIDA) const{
    if(fConfig.ApplyMaxTrackLengthPIDACut){return MaxTrackLengthPIDACut(MaxTrackLength,MaxTrackPIDA);}
    return true;
  }
  bool passMVACut() const{
    if(fConfig.ApplyMVACut){return MVACut();}
    return true;
  }


  bool OneShowerEnergyCut(std::vector<art::Ptr<recob::Shower> > const& neutrino_showers) const;

  bool OneShowerResiudalCut(std::vector<art::Ptr<recob::Shower> > const& neutrino_showers,double& numshowers, const TVector3& vertex_position);

  bool ConversionGapCut(const TVector3 ShowerStart, double& conversion_gap,const TVector3& vertex_position);

  bool NeutrinoPdgCodeCut(const double& NeutrinoPdgCode) const;
  bool NumneutrinosCut(const double& Numneutrinos) const;
  bool dEdxCut(const double& dEdx) const;
  bool LengthCut(const double& Length) const;
  bool OpeningAngleCut(const double& OpeningAngle) const;
  bool ShowerEnergyCut(const double& ShowerEnergy) const;
  bool ShowerDensityGradientCut(const double& ShowerDensityGradient) const;
  bool ShowerDensityPowerCut(const double& ShowerDensityPower) const;
  bool ShowerTrackLengthCut(const double& ShowerTrackLength) const;
  bool ShowerTrackWidthCut(const double& ShowerTrackWidth) const;
  bool MaxTrackLengthCut(const double& MaxTrackLength) const;
  bool MaxTrackPIDACut(const double& MaxTrackPIDA) const;
  bool MaxTrackLengthPIDACut(const double& MaxTrackLength, const double& MaxTrackPIDA) const;
  bool MVACut() const;

      
  //Cut Functions
  bool containedInFV(const TVector3 &v);
  bool containedInAV(const TVector3 &v);

  
  void SetupMVA();
     
  void InitialiseHistograms();
  void FillHistograms(std::map<std::string,TH1D*>& HistMap, NueSelectionReco::NueInteraction& intInfo);
  void FillHistograms(std::map<std::string,TH1D*>& HistMap, const simb::MCNeutrino& nu,
                      NueSelectionReco::NueInteraction& intInfo, double Energy, bool &booldirtevent);


  void GetSliceHits(const art::Ptr<recob::PFParticle>& pfp, 
		    std::vector< art::Ptr< recob::Hit> >& pfpHits,
		    std::map<int, art::Ptr<recob::PFParticle> >& pfpMap,
		    art::FindManyP<recob::Cluster>& fmpfc,
		    art::FindManyP<recob::Hit>& fmch);

  double ShowerDensityGradient(std::vector<art::Ptr<recob::SpacePoint> >& sps, 
			       const TVector3& ShowerStartPosition, 
			       const TVector3& ShowerDirection,
			       const double ShowerLength, 
			       const double& OpenAngle, 
			       art::FindManyP<recob::Hit> const& fmh, 
			       const double ShowerEnergy, 
			       double& pw);
  
  double TotalEnergy(const std::vector<art::Ptr<recob::SpacePoint> >&sps, 
		     art::FindManyP<recob::Hit> const& fmh);
  
  int SpacePointPlane(art::Ptr<recob::SpacePoint> const& sp,
		      art::FindManyP<recob::Hit> const& fmh) const;
  

  double SpacePointTime(art::Ptr<recob::SpacePoint> const& sp,
			art::FindManyP<recob::Hit> const& fmh) const; 
								
  double SpacePointCharge(art::Ptr<recob::SpacePoint> const& sp,
			  art::FindManyP<recob::Hit> const& fmh) const;
  
  double SpacePointProjection(const art::Ptr<recob::SpacePoint>&sp,
			      TVector3 const& vertex, 
			      TVector3 const& direction) const;
  
  TVector3 SpacePointPosition(art::Ptr<recob::SpacePoint> const& sp) const;

  double ShowerTrackWidthCal(std::vector<art::Ptr<recob::Hit> > const& showertrackhits
			     ,art::FindManyP<recob::SpacePoint> const& fmsphsp,
			     const TVector3& ShowerStart,
			     const TVector3& ShowerDirection) const;
  
  std::vector<int> ElectronNuDecays = {1,2,6,9};
  std::vector<int> MuonNuDecays = {3,4,5,7,8,10,11,12,13,14};

  std::clock_t c_start;
  std::clock_t c_end;

  TRandom3 rand;
  double fPOT;
  int MCTruthCounter; 
  bool dirtevent;

  TMVA::Reader *reader;
  std::map<std::string, float> MVAValues;

  cheat::ParticleInventory const *fpi_service;
  detinfo::DetectorPropertiesStandard const *fDetProp;
};

  }  // namespace SBNOsc
}  // namespace ana

#endif  // __sbnanalysis_ana_SBNOsc_NueSelectionReco__

