/**
 * @class SimpleFinder
 * @brief V0 particle (Lambda) reconstruction algorithm
 * @authors Oleksii Lubynets, Viktor Klochkov, Ilya Selyuzhenkov, Susanne Glaessel
 *
 * Simplified version of KFParticleFinder class. At the current moment it is developed to reconstruct 2- and 3-body-decays. \n
 * SimpleFinder is based on KFParticle package, it uses mathematical apparatus implemented in KFParticle.
 * The basic idea and reconstruction algorithm also reproduce the KFParticle, but SimpleFinder is free of overloading of too complicated procedure as in KFParticleFinder. \n
 * Also the advantage of SimpleFinder is that reconstruction procedure is clear and under full control of user, almost in the "hand mode".
 * It gives a possibility of detailed analysis of V0 reconstruction, in particular of decay parameters distributions in order to optomize selection criterias (cuts).\n
 */

#ifndef KFPARTICLESIMPLE_KFSIMPLE_SIMPLEFINDER_HPP_
#define KFPARTICLESIMPLE_KFSIMPLE_SIMPLEFINDER_HPP_

#include <deque>
#include <map>
#include <vector>

#include <KFPTrack.h>
#include <KFParticle.h>
#include <KFParticleSIMD.h>
#include <KFVertex.h>
#include <TVector3.h>

#include "InputContainer.hpp"
#include "OutputContainer.hpp"

#include "Constants.hpp"
#include "Decay.hpp"
#include "NonLinearCutBase.hpp"

class SimpleFinder {

  typedef std::array<float, 8> Param_t;
  typedef std::vector<Param_t> Parameters_t;

 public:
  SimpleFinder() = default;
  SimpleFinder(const SimpleFinder&) = default;
  SimpleFinder(SimpleFinder&&) = default;
  SimpleFinder& operator=(SimpleFinder&&) = default;
  SimpleFinder& operator=(const SimpleFinder&) = default;
  ~SimpleFinder() = default;

  void Init(std::vector<KFParticle>&& tracks, const KFVertex& pv);///< Initialize SimpleFinder object with PV and set of tracks of the current event
  void Init(const InputContainer& input);

  void InitTrackIndexesMap();

  void FindParticles() {
    for (const auto& decay : decays_) {
      InitTrackIndexesMap();
      ReconstructDecay(decay);
    }
    if (is_event_mixing_)
      FindMixedBackground();
  }

  void ReconstructDecay(const Decay& decay);
  void SetDecays(const std::vector<Decay>& decays) { decays_ = decays; }
  void AddDecay(const Decay& decay) { decays_.emplace_back(decay); }
  const std::vector<OutputContainer>& GetCandidates() const { return output_; }

  ///< Enable/disable background estimation via event mixing.
  void SetIsEventMixing(bool is = true) { is_event_mixing_ = is; }
  ///< PDG of the track species replaced by tracks from previous events (track-level mixing).
  void SetMixPdg(int pdg) { mix_pdg_ = pdg; }
  ///< Number of previous events whose tracks are combined with the current event to estimate the background.
  void SetMixBufferSize(int n) { mix_buffer_size_ = n; }
  ///< Background candidates obtained by combining tracks from different events (event mixing).
  const std::vector<OutputContainer>& GetMixedCandidates() const { return mixed_output_; }

 private:
  std::vector<KFParticle> tracks_;    ///< input information: vector of tracks
  std::vector<KFParticle> input_tracks_;///< snapshot of the current event input tracks (used for event mixing)
  KFVertex prim_vx_;                  ///< input information: primary vertex
  std::array<float, 3> sec_vx_;       ///< input information: primary vertex
  std::vector<Decay> decays_{};       ///< input information: list of decays to reconstruct
  NonLinearCutBase* ml_cuts_{nullptr};///< input information: non-linear cuts class (optional)

  std::map<Pdg_t, std::vector<int>> indexes_{};///< map of indexes for a given particle specie
  std::map<int, int> id2index_{};              ///< map of track Ids to indexes
  Parameters_t params_{};                      ///< vector of daughter parameters at current SV estimation
  SelectionValues values_{};                   ///< struct with mother and daughters properties used to apply cuts

  std::vector<OutputContainer> output_{};///< output information: vector of candidates
  std::vector<OutputContainer> mixed_output_{};///< output information: mixed-event background candidates
  int current_candidate_id_{0};
  int last_track_id_{0};

  // event mixing
  struct MixEvent {
    KFVertex pv;                      ///< primary vertex of the buffered event
    std::vector<KFParticle> mix_tracks;///< input tracks of the buffered event with PDG == mix_pdg_
  };
  ///< Number of centrality bins used for the event-mixing buffer (mirrors AFKF).
  static constexpr int kMixCentralityBins = 9;
  ///< Ring buffers of previous events' mix_pdg tracks, one per centrality bin.
  ///< An event only mixes with the last `mix_buffer_size_` events of the *same* bin.
  std::deque<MixEvent> mixed_events_buffer_[kMixCentralityBins];
  float current_centrality_{0.f};///< centrality of the current event (for bin selection)
  bool is_event_mixing_{false};
  int mix_pdg_{-211};///< PDG of the track species replaced by tracks from previous events
  int mix_buffer_size_{5};///< number of previous events kept in the mixing buffer (per bin)
  bool is_mixed_pass_{false};

  /**
  * Find indexes of good daughters
  * @param cuts daughter particle cuts container
  * @return vector of indexes
  */
  std::vector<int> GetGoodDaughtersIndexes(const Daughter& cuts);

  bool IsGoodDaughter(const KFParticle& track, const Daughter& cuts);
  bool IsGoodPair(const KFParticle& track1, const KFParticle& track2, const Decay& decay);
  bool IsGoodThree(const KFParticle& track1, const KFParticle& track2, const KFParticle& track3, const Decay& decay);
  bool IsGoodMother(const KFParticleSIMD& mother, const Mother& cuts, int id_mother);
  bool IsGoodMotherMass(const KFParticleSIMD& kf_mother, const Mother& mother);
  bool IsMotherFromPV(const KFParticleSIMD& mother, const Mother& cuts, int id_mother);
  bool IsGoodDecayLength(const KFParticleSIMD& mother, const Mother& cuts);
  bool IsGoodCos(const KFParticleSIMD& mother, const Decay& decay);
  bool IsGoodTrackIdsGenerations(const std::vector<KFParticle> tracks);

  KFParticleSIMD ConstructMother(const std::vector<KFParticle>& tracks, const std::vector<Pdg_t>& pdgs, const Decay& decay, bool print = false);
  void CalculateSecondaryVertex();

  void CalculateParamsInPCA(const KFParticle& track1, int pid1, const KFParticle& track2, int pid2);
  void CalculateParamsInSV(const std::array<KFParticle, 3>& tracks, std::vector<Pdg_t> pids);
  float CalculateChiToPrimaryVertex(const KFParticle& track, Pdg_t pid) const;
  float CalculateInvMassDiscrepancy(const KFParticle& track, const Mother& mother) const;
  float CalculateCosTopo(const KFParticleSIMD& mother) const;
  float CalculateDistanceBetweenParticles();
  float CalculateDistanceBetweenParticles(const KFParticleBaseSIMD& particle1, const KFParticleBaseSIMD& particle2);
  float CalculateDistanceToSV() const;
  float CalculateCosOpen(const int id_daughter_1, const int id_daughter_2) const;

  void FillDaughtersInfo(const std::vector<KFParticle>& tracks, const std::vector<Pdg_t>& pdgs);
  void SaveParticle(KFParticleSIMD& particle_simd, const Decay& decay);

  ///< Reconstruct background candidates by combining current event tracks with tracks of previous events.
  void FindMixedBackground();
  ///< Translate a buffered track from the source event PV to the target (current) event PV.
  KFParticle TranslateMixTrack(const KFParticle& track, const KFVertex& source_pv, const KFVertex& target_pv) const;
  ///< Check whether the given PDG can be used as a track-level mixing species.
  bool IsSupportedTrackPdg(int pdg) const;
  ///< Map an event `centrality_tracks` value to a 0-8 centrality bin (mirrors AFKF's
  ///< `Centrality_AT` with `centFull[9] = {80,70,60,50,40,30,20,10,5}`).
  static int CentralityAT(float gRefMult);

  bool ApplyNonLinearCut() const {
    return ml_cuts_ == nullptr || ml_cuts_->ApplyCut(values_);
  }

  KFParticle GetTrack(int i) {
    return tracks_.at(i);
  }

  static KFPTrack ToKFPTrack(const KFParticle& particle);
  void SetKFParticleEnergy(KFParticle& particle, int pdg) const;
};

#endif//KFPARTICLESIMPLE_KFSIMPLE_SIMPLEFINDER_HPP_
