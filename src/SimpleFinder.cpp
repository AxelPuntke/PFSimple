#include "SimpleFinder.hpp"
#include <KFParticleSIMD.h>

#include <algorithm>
#include <numeric>
#include <iostream>
using std::to_string;

void SimpleFinder::Init(std::vector<KFParticle>&& tracks, const KFVertex& pv) {
  output_.clear();
  mixed_output_.clear();
  current_candidate_id_ = 0;

  if (!tracks.empty())
    last_track_id_ = tracks.at(tracks.size() - 1).Id();
  else
    last_track_id_ = 0;

  tracks_ = tracks;
  input_tracks_ = tracks;
  prim_vx_ = pv;
}

void SimpleFinder::Init(const InputContainer& input) {// TODO know how to de-const without time lose for copying by value
  std::vector<KFParticle> tracks = input.GetTracks();
  Init(std::move(tracks), input.GetVertex());
  current_centrality_ = input.GetCentrality();
  current_n_good_tracks_ = input.GetNGoodTracks();
  centrality_source_ = input.GetCentralitySource();
}

int SimpleFinder::CentralityAT(float gRefMult) {
  // Mirrors AFKF's Centrality_AT with centFull[9] = {80,70,60,50,40,30,20,10,5}.
  // Bin 8 = most central (lowest centrality value), bin 0 = most peripheral.
  static const int centFull[kMixCentralityBins] = {80, 70, 60, 50, 40, 30, 20, 10, 5};
  int centrality;
  if (gRefMult <= centFull[8]) centrality = 8;
  else if (gRefMult <= centFull[7]) centrality = 7;
  else if (gRefMult <= centFull[6]) centrality = 6;
  else if (gRefMult <= centFull[5]) centrality = 5;
  else if (gRefMult <= centFull[4]) centrality = 4;
  else if (gRefMult <= centFull[3]) centrality = 3;
  else if (gRefMult <= centFull[2]) centrality = 2;
  else if (gRefMult <= centFull[1]) centrality = 1;
  else if (gRefMult <= centFull[0]) centrality = 0;
  else centrality = -1;
  return (centrality < 0 || centrality > 8) ? 0 : centrality;
}

int SimpleFinder::CentralityNTacks(int nTracks) {
  // Mirrors AFKF's Centrality(n_tracks) with centFull[9] = {6,12,22,38,62,95,145,212,257}.
  // Bin 8 = most central (highest track count), bin 0 = most peripheral.
  static const int centFull[kMixCentralityBins] = {6, 12, 22, 38, 62, 95, 145, 212, 257};
  int centrality;
  if (nTracks >= centFull[8]) centrality = 8;
  else if (nTracks >= centFull[7]) centrality = 7;
  else if (nTracks >= centFull[6]) centrality = 6;
  else if (nTracks >= centFull[5]) centrality = 5;
  else if (nTracks >= centFull[4]) centrality = 4;
  else if (nTracks >= centFull[3]) centrality = 3;
  else if (nTracks >= centFull[2]) centrality = 2;
  else if (nTracks >= centFull[1]) centrality = 1;
  else if (nTracks >= centFull[0]) centrality = 0;
  else centrality = -1;
  return (centrality < 0 || centrality > 8) ? 0 : centrality;
}

float SimpleFinder::CalculateDistanceBetweenParticles() {
  const float dx = params_.at(0).at(kX) - params_.at(1).at(kX);
  const float dy = params_.at(0).at(kY) - params_.at(1).at(kY);
  const float dz = params_.at(0).at(kZ) - params_.at(1).at(kZ);
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

float SimpleFinder::CalculateDistanceBetweenParticles(const KFParticleBaseSIMD& particle1, const KFParticleBaseSIMD& particle2) {
  return particle2.GetDistanceFromParticle(particle1)[0];
}

void SimpleFinder::CalculateParamsInPCA(const KFParticle& track1, int pid1, const KFParticle& track2, int pid2) {
  params_ = Parameters_t(2);

  KFParticle particle1(track1);
  KFParticle particle2(track2);
  particle1.SetPDG(pid1);
  particle2.SetPDG(pid2);
  KFParticleSIMD particleSIMD1(particle1);// the same particle is copied to each SIMD element
  KFParticleSIMD particleSIMD2(particle2);

  float_v dS[2];
  float_v params1[8], params2[8];
  particleSIMD1.GetDStoParticleFast(particleSIMD2, dS);
  particleSIMD1.TransportFast(dS[0], params1);
  particleSIMD2.TransportFast(dS[1], params2);

  for (int i = 0; i < 8; i++) {
    params_.at(0).at(i) = params1[i][0];
    params_.at(1).at(i) = params2[i][0];
  }
}

void SimpleFinder::CalculateParamsInSV(const std::array<KFParticle, 3>& tracks, std::vector<Pdg_t> pids) {

  std::vector<std::array<float, 3>> sec_vx;
  sec_vx.resize(3);

  int k = 0;
  for (int i = 0; i < tracks.size(); i++) {
    for (int j = i + 1; j < tracks.size(); j++) {
      CalculateParamsInPCA(tracks.at(i), pids.at(i), tracks.at(j), pids.at(j));
      CalculateSecondaryVertex();
      sec_vx.at(k) = sec_vx_;
      k++;
    }
  }

  for (int i = 0; i < 3; i++)
    sec_vx_.at(i) = (sec_vx.at(0).at(i) + sec_vx.at(1).at(i) + sec_vx.at(2).at(i)) / 3;

  params_.resize(3);

  for (int itrack = 0; itrack < 3; itrack++) {
    KFParticle particle(tracks.at(itrack));
    particle.SetPDG(pids.at(itrack));
    KFParticleSIMD particleSIMD(particle);

    float_v dS;
    float_v dsdr[6];
    float_v params[8];
    float_v sv[3];

    for (int i = 0; i < 3; i++)
      sv[i] = sec_vx_.at(i);

    dS = particleSIMD.GetDStoPoint(sv, dsdr);
    particleSIMD.TransportFast(dS, params);

    for (int i = 0; i < 8; i++)
      params_.at(itrack).at(i) = params[i][0];
  }
}

float SimpleFinder::CalculateDistanceToSV() const {
  const float dx = params_.at(2).at(0) - sec_vx_.at(0);
  const float dy = params_.at(2).at(1) - sec_vx_.at(1);
  const float dz = params_.at(2).at(2) - sec_vx_.at(2);
  const float dr = sqrt(dx * dx + dy * dy + dz * dz);
  return dr;
}

KFParticleSIMD SimpleFinder::ConstructMother(const std::vector<KFParticle>& tracks, const std::vector<Pdg_t>& pdgs, const Decay& decay, bool print) {
  const auto n = tracks.size();

  KFParticleSIMD mother;
  std::vector<KFParticle> particles = tracks;
  std::vector<KFParticleSIMD> particles_simd{};

  for (size_t i = 0; i < n; ++i) {
    if (tracks.at(i).Id() <= last_track_id_)
      SetKFParticleEnergy(particles.at(i), pdgs.at(i));
    particles.at(i).SetPDG(pdgs.at(i));
    particles.at(i).SetId(tracks.at(i).Id());// TODO rm obsolet ID copying?
    particles_simd.emplace_back(particles.at(i));
  }

  if (n == 2) {
    float_v ds[2] = {0.f, 0.f};
    float_v dsdr[4][6];
    particles_simd.at(0).GetDStoParticle(particles_simd.at(1), ds, dsdr);
    particles_simd.at(0).TransportToDS(ds[0], dsdr[0]);
    particles_simd.at(1).TransportToDS(ds[1], dsdr[3]);
  } else {
    float_v vertex[3] = {sec_vx_[0], sec_vx_[1], sec_vx_[2]};
    for (size_t i = 0; i < n; ++i) {
      float_v ds, dsdr[6];
      ds = particles_simd.at(i).GetDStoPoint(vertex, dsdr);
      particles_simd.at(i).TransportToDS(ds, dsdr);
    }
  }

  if (n == 2) {
    const KFParticleSIMD* vDaughtersPointer[2] = {&particles_simd.at(0), &particles_simd.at(1)};
    mother.Construct(vDaughtersPointer, 2, nullptr);

    if ((tracks.at(0).Id() <= last_track_id_ || tracks.at(1).Id() <= last_track_id_) && print) {
      std::cout << particles_simd.at(0).GetPx()[0] << "\t" << particles_simd.at(0).GetPy()[0] << "\t" << particles_simd.at(0).GetPz()[0] << "\t" << particles_simd.at(0).GetE()[0] << "\n";
      std::cout << particles_simd.at(1).GetPx()[0] << "\t" << particles_simd.at(1).GetPy()[0] << "\t" << particles_simd.at(1).GetPz()[0] << "\t" << particles_simd.at(1).GetE()[0] << "\n";
      std::cout << mother.GetPx()[0] << "\t" << mother.GetPy()[0] << "\t" << mother.GetPz()[0] << "\t" << mother.GetE()[0] << "\n";
    }

  } else if (n == 3) {
    const KFParticleSIMD* vDaughtersPointer[3] = {&particles_simd.at(0), &particles_simd.at(1), &particles_simd.at(2)};
    mother.Construct(vDaughtersPointer, 3, nullptr);
  }

  return mother;
}

float SimpleFinder::CalculateChiToPrimaryVertex(const KFParticle& track, Pdg_t pid) const {
  // Match AFKF ComputeChiToVertex (ATKFParticleFinder.cxx) exactly: construct a scalar
  // KFParticle from the track under the given PDG hypothesis and call
  // GetDeviationFromVertex once on the PV.  AFKF does a single transport inside
  // GetDeviationFromVertex; the previous SIMD path added an extra TransportToPoint
  // (and a SIMD round-trip) that lowered chi2_prim and undercounted secondary tracks.
  KFParticle tmpPart = track;
  SetKFParticleEnergy(tmpPart, pid);  // KFParticle(ToKFPTrack(particle), pid) -- same as AFKF
  return tmpPart.GetDeviationFromVertex(prim_vx_);
}

float SimpleFinder::CalculateInvMassDiscrepancy(const KFParticle& track, const Mother& mother) const {
  auto mothertmp = mother;
  float massPDG = mothertmp.GetMassPdg();
  float massPDGSigma = mothertmp.GetMassPdgSigma();
  return (track.GetMass() - massPDG) / massPDGSigma;
}

std::vector<int> SimpleFinder::GetGoodDaughtersIndexes(const Daughter& daughter) {
  std::vector<int> result{};
  for (auto pid : daughter.GetPids()) {
    auto it = indexes_.find(pid);
    if (it != indexes_.end()) {
      for (auto i_track : it->second) {
        auto track = GetTrack(i_track);
        if (IsGoodDaughter(track, daughter)) {
          result.emplace_back(i_track);
        }
      }
    }
  }
  // Accept nuclear fragments (|pdg|>=1e9) as proton-hypothesis positive daughters, mirroring
  // AFKF's Find2DaughterDecay active[1] (posPDG>1e9).  Fragments are always secondary, so
  // they bypass the chi2prim>=cut gate (see IsGoodDaughter is_fragment).
  if (std::abs(daughter.GetPdgHypo()) == 2212) {
    for (const auto& kv : indexes_) {
      if (std::abs(kv.first) > 1000000000) {
        for (auto i_track : kv.second) {
          auto track = GetTrack(i_track);
          if (IsGoodDaughter(track, daughter, true)) {
            result.emplace_back(i_track);
          }
        }
      }
    }
  }
  return result;
}

void SimpleFinder::InitTrackIndexesMap() {
  indexes_.clear();
  id2index_.clear();
  for (size_t i = 0; i < tracks_.size(); i++) {
    id2index_[tracks_.at(i).Id()] = i;
    auto pdg = tracks_.at(i).GetPDG();
    auto it = indexes_.find(pdg);
    if (it != indexes_.end()) {
      it->second.emplace_back(i);
    } else {
      indexes_[pdg] = {static_cast<int>(i)};
    }
  }
}

KFPTrack SimpleFinder::ToKFPTrack(const KFParticle& particle) {
  KFPTrack track;

  track.SetX(particle.GetX());
  track.SetY(particle.GetY());
  track.SetZ(particle.GetZ());
  track.SetPx(particle.GetPx());
  track.SetPy(particle.GetPy());
  track.SetPz(particle.GetPz());

  for (Int_t iC = 0; iC < NumberOfCovElements; iC++) {
    track.SetCovariance(iC, particle.GetCovariance(iC));
  }
  for (Int_t iF = 0; iF < 10; iF++) {
    track.SetFieldCoeff(particle.GetFieldCoeff()[iF], iF);
  }
  track.SetCharge(particle.GetQ());
  track.SetId(particle.Id());

  return track;
}

void SimpleFinder::SetKFParticleEnergy(KFParticle& particle, int pdg) const {
  const KFPTrack& track = ToKFPTrack(particle);
  particle = KFParticle(track, pdg);
}

bool SimpleFinder::IsGoodDaughter(const KFParticle& track, const Daughter& daughter, bool is_fragment) {
  int id = daughter.GetId();
  values_.chi2_prim.at(id) = CalculateChiToPrimaryVertex(track, daughter.GetPdgHypo());

  // AFKF: for |pdg|>=1e9 (fragments) fChiPrimaryCutFragments=0.0, so a fragment is always
  // secondary (pv_index=-1) regardless of chi_primary; only gate is chi_primary>0.0.
  // Ordinary particles must be secondary, i.e. chi_primary >= cut (chi2prim=5.0).
  if (is_fragment) {
    if (values_.chi2_prim.at(id) <= 0.0f || std::isnan(values_.chi2_prim.at(id))) { return false; }
  } else if (values_.chi2_prim.at(id) < daughter.GetCutChi2Prim() || std::isnan(values_.chi2_prim.at(id))) { return false; }

  if (!(track.Id() > last_track_id_)) {
    if (is_fragment) {
      // AFKF accepts nuclear fragments as positive daughters; He3/He4 carry Q=+2 (charge doubled
      // later in ReconstructDecay).  Only require positive charge here.
      if (track.Q() <= 0) { return false; }
    } else if (TMath::Abs(daughter.GetPdgHypo()) == 3112 || TMath::Abs(daughter.GetPdgHypo()) == 3312 || TMath::Abs(daughter.GetPdgHypo()) == 3334 || TMath::Abs(daughter.GetPdgHypo()) == 11 || TMath::Abs(daughter.GetPdgHypo()) == 13) {
      if ((TMath::Sign(1, daughter.GetPdgHypo()) > 0 && track.Q() != -1) || (TMath::Sign(1, daughter.GetPdgHypo()) < 0 && track.Q() != +1)) { return false; }
    } else {
      if ((TMath::Sign(1, daughter.GetPdgHypo()) > 0 && track.Q() != +1) || (TMath::Sign(1, daughter.GetPdgHypo()) < 0 && track.Q() != -1)) { return false; }
    }
  }

  return true;
}

bool SimpleFinder::IsGoodPair(const KFParticle& track1,
                              const KFParticle& track2,
                              const Decay& decay) {

  if (track1.Id() == track2.Id()) { return false; }

  std::vector<KFParticle> tracks = {track1, track2};
  if (IsGoodTrackIdsGenerations(tracks) == false) { return false; }

  int generation1 = track1.Id() > last_track_id_ ? 1 : 0;
  int generation2 = track2.Id() > last_track_id_ ? 1 : 0;

  std::vector<int> generation;
  for (int i = 0; i < tracks.size(); i++)
    generation.push_back(tracks[i].Id() > last_track_id_ ? 1 : 0);

  if (std::accumulate(generation.begin(), generation.end(), 0) > 0) {
    KFParticle track1_nonconst = track1;
    KFParticle track2_nonconst = track2;
    KFParticleSIMD track1_simd = KFParticleSIMD(track1_nonconst);
    KFParticleSIMD track2_simd = KFParticleSIMD(track2_nonconst);
    track1_simd.SetPDG(decay.GetDaughters().at(0).GetPdgHypo());
    track2_simd.SetPDG(decay.GetDaughters().at(1).GetPdgHypo());
    values_.distance = CalculateDistanceBetweenParticles(track1_simd, track2_simd);
  } else
    values_.distance = CalculateDistanceBetweenParticles();

  if (values_.distance > decay.GetMother().GetCutDistance() || std::isnan(values_.distance)) { return false; }

  const int id_mother = decay.GetNDaughters() == 2 ? 0 : 1;
  values_.cos_open[id_mother] = CalculateCosOpen(0, 1);
  if (values_.cos_open[id_mother] < decay.GetMother().GetCutCosOpen()[id_mother] || std::isnan(values_.cos_open[id_mother])) { return false; }

  return true;
}

bool SimpleFinder::IsGoodThree(const KFParticle& track1,
                               const KFParticle& track2,
                               const KFParticle& track3,
                               const Decay& decay) {

  if (track1.Id() == track2.Id()) { return false; }
  if (track1.Id() == track3.Id()) { return false; }
  if (track2.Id() == track3.Id()) { return false; }

  std::vector<KFParticle> tracks = {track1, track2, track3};
  if (IsGoodTrackIdsGenerations(tracks) == false) { return false; }

  values_.distance_sv = CalculateDistanceToSV();

  if (values_.distance_sv > decay.GetMother().GetCutDistanceToSV() || std::isnan(values_.distance_sv)) { return false; }

  for (int i = 0; i < 2; i++) {
    values_.cos_open[i + 2] = CalculateCosOpen(i, 2);
    if (values_.cos_open[i + 2] < decay.GetMother().GetCutCosOpen()[i + 2] || std::isnan(values_.cos_open[i + 2])) { return false; }
  }

  std::vector<float> cos_tmp;
  for (int i = 1; i < 4; i++)
    cos_tmp.push_back(values_.cos_open[i]);

  values_.cos_open[0] = *std::min_element(cos_tmp.begin(), cos_tmp.end());
  if (values_.cos_open[0] < decay.GetMother().GetCutCosOpen()[0] || std::isnan(values_.cos_open[0])) { return false; }

  return true;
}

bool SimpleFinder::IsGoodMother(const KFParticleSIMD& kf_mother, const Mother& mother, const int id_mother) {
  values_.chi2_geo[id_mother] = kf_mother.Chi2()[0] / simd_cast<float_v>(kf_mother.NDF())[0];
  if (values_.chi2_geo[id_mother] > mother.GetCutChi2Geo()[id_mother] || std::isnan(values_.chi2_geo[id_mother])) { return false; }
  if (values_.chi2_geo[id_mother] < 0. || !std::isfinite(values_.chi2_geo[id_mother])) { return false; }
  return true;
}

bool SimpleFinder::IsGoodMotherMass(const KFParticleSIMD& kf_mother, const Mother& mother) {
  KFParticle particle;
  KFParticleSIMD(kf_mother).GetKFParticle(particle, 0);
  values_.invmassdisc = CalculateInvMassDiscrepancy(particle, mother);
  if (std::fabs(values_.invmassdisc) > mother.GetCutInvMass()) { return false; }
  return true;
}

bool SimpleFinder::IsMotherFromPV(const KFParticleSIMD& kf_mother, const Mother& mother, const int id_mother) {
  KFVertex prim_vx_tmp = prim_vx_;
  const KFParticleSIMD prim_vx_Simd(prim_vx_tmp);
  KFParticleSIMD motherTopo = kf_mother;
  motherTopo.SetProductionVertex(prim_vx_Simd);

  values_.chi2_topo[id_mother] = motherTopo.GetChi2()[0] / simd_cast<float_v>(motherTopo.GetNDF())[0];
  values_.cos_topo[id_mother] = CalculateCosTopo(kf_mother);

  if (id_mother > 0) {
    std::array<bool, 2> is_from_pv = {true, true};
    if (values_.chi2_topo[id_mother] > mother.GetCutChi2Topo()[id_mother] || std::isnan(values_.chi2_topo[id_mother])) is_from_pv.at(0) = false;
    if (values_.cos_topo[id_mother] < mother.GetCutCosTopo()[id_mother] || std::isnan(values_.cos_topo[id_mother])) is_from_pv.at(1) = false;
    if ((is_from_pv.at(0) == false && is_from_pv.at(1) == false) || (is_from_pv.at(0) == false && mother.GetCutCosTopo()[id_mother] == huge_value) || (mother.GetCutChi2Topo()[id_mother] == -huge_value && is_from_pv.at(1) == false))
      return false;
  } else {
    if (values_.chi2_topo[id_mother] > mother.GetCutChi2Topo()[id_mother] || std::isnan(values_.chi2_topo[id_mother])) return false;
    if (values_.cos_topo[id_mother] < mother.GetCutCosTopo()[id_mother] || std::isnan(values_.cos_topo[id_mother])) return false;
  }

  KFParticle particle;
  KFParticleSIMD(kf_mother).GetKFParticle(particle, 0);
  values_.chi2_prim_mother = CalculateChiToPrimaryVertex(particle, mother.GetPdg());

  return true;
}

bool SimpleFinder::IsGoodDecayLength(const KFParticleSIMD& kf_mother, const Mother& mother) {

  float_v d_pv_Simd, dd_pv_Simd, l_Simd, dl_Simd;
  float_m isFromPV_Simd;
  KFVertex prim_vx_tmp = prim_vx_;
  const KFParticleSIMD prim_vx_Simd(prim_vx_tmp);

  kf_mother.GetDistanceToVertexLine(prim_vx_Simd, d_pv_Simd, dd_pv_Simd, &isFromPV_Simd);
  values_.distance_pv = d_pv_Simd[0];
  if (values_.distance_pv < mother.GetCutDistancePVLine() || std::isnan(values_.distance_pv)) { return false; }

  KFParticleSIMD motherTopo = kf_mother;
  motherTopo.SetProductionVertex(prim_vx_Simd);
  motherTopo.KFParticleBaseSIMD::GetDecayLength(l_Simd, dl_Simd);

  values_.l_over_dl = l_Simd[0] / dl_Simd[0];
  if (values_.l_over_dl < mother.GetCutLdL() || std::isnan(values_.l_over_dl)) { return false; }

  values_.l = l_Simd[0];
  if (values_.l < mother.GetCutDecayLength() || std::isnan(values_.l)) { return false; }

  values_.is_from_PV = isFromPV_Simd[0];

  return true;
}

void SimpleFinder::SaveParticle(KFParticleSIMD& particle_simd, const Decay& decay) {

  KFParticle particle;

  auto mothertmp = decay.GetMother();
  if (decay.GetIsApplyMassConstraint())
    particle_simd.SetNonlinearMassConstraint(float_v(mothertmp.GetMassPdg()));

  if (decay.GetIsTransportToPV()) {
    KFVertex prim_vx_tmp = prim_vx_;
    const KFParticleSIMD prim_vx_Simd(prim_vx_tmp);
    const float_v point[3] = {prim_vx_Simd.X(), prim_vx_Simd.Y(), prim_vx_Simd.Z()};
    particle_simd.TransportToPoint(point);
  }

  particle_simd.GetKFParticle(particle, 0);
  particle.SetPDG(decay.GetMother().GetPdg());

  int track_id = last_track_id_ + current_candidate_id_ + 1;
  particle.SetId(track_id);

  tracks_.emplace_back(particle);
  current_candidate_id_++;

  if (decay.GetIsDoNotWriteMother() == true) return;

  OutputContainer mother(particle);
  mother.SetSelectionValues(values_);
  mother.SetId(track_id);

  std::vector<int> daughters_generation;
  for (int i = 0; i < particle.NDaughters(); i++)
    if (particle.DaughterIds().at(i) > last_track_id_)
      daughters_generation.push_back(1);
    else
      daughters_generation.push_back(0);
  mother.SetDaughterGenerations(daughters_generation);

  if (is_mixed_pass_) {
    // All candidates built in the mixing pass are combinatorial background.
    mixed_output_.emplace_back(mother);
  } else {
    output_.emplace_back(mother);
  }
}

float SimpleFinder::CalculateCosTopo(const KFParticleSIMD& mother) const {
  const float px = mother.GetPx()[0];
  const float py = mother.GetPy()[0];
  const float pz = mother.GetPz()[0];

  const float delta_x = mother.GetX()[0] - prim_vx_.GetX();
  const float delta_y = mother.GetY()[0] - prim_vx_.GetY();
  const float delta_z = mother.GetZ()[0] - prim_vx_.GetZ();

  const float sp = delta_x * px + delta_y * py + delta_z * pz;
  const float norm = std::sqrt(delta_x * delta_x + delta_y * delta_y + delta_z * delta_z) * mother.GetP()[0];
  return sp / norm;
}

bool SimpleFinder::IsGoodCos(const KFParticleSIMD& mother, const Decay& decay) {
  for (int i = 0; i < decay.GetNDaughters(); ++i) {
    const auto cut = decay.GetDaughters()[i].GetCutCos();
    const auto& par = params_.at(i);
    const float norm = mother.GetP()[0] * std::sqrt(par[kPx] * par[kPx] + par[kPy] * par[kPy] + par[kPz] * par[kPz]);
    values_.cos[i] = (mother.GetPx()[0] * par[kPx] + mother.GetPy()[0] * par[kPy] + mother.GetPz()[0] * par[kPz]) / norm;
    if (values_.cos[i] < cut || std::isnan(values_.cos[i])) {
      return false;
    }
  }
  return true;
}

float SimpleFinder::CalculateCosOpen(const int id_daughter_1, const int id_daughter_2) const {
  const auto& par1 = params_.at(id_daughter_1);
  const auto& par2 = params_.at(id_daughter_2);
  const float sp = par1[kPx] * par2[kPx] + par1[kPy] * par2[kPy] + par1[kPz] * par2[kPz];
  const float norm = std::sqrt(par1[kPx] * par1[kPx] + par1[kPy] * par1[kPy] + par1[kPz] * par1[kPz]) * std::sqrt(par2[kPx] * par2[kPx] + par2[kPy] * par2[kPy] + par2[kPz] * par2[kPz]);

  return sp / norm;
}

void SimpleFinder::CalculateSecondaryVertex() {
  if (params_.size() < 2 || params_.size() > 3) {
    throw std::runtime_error("Daughter parameters size is wrong");
  }
  for (int i = 0; i < 3; ++i)
    sec_vx_.at(i) = (params_.at(0).at(kX + i) + params_.at(1).at(kX + i)) / 2;
}

bool SimpleFinder::IsGoodTrackIdsGenerations(const std::vector<KFParticle> tracks) {

  int n = tracks.size();

  std::vector<std::vector<int>> trackIdsGen0;
  std::vector<std::vector<KFParticle>> tracksLevel;

  tracksLevel.resize(1);
  for (int i = 0; i < n; i++)
    tracksLevel.at(0).push_back(tracks[i]);

  int ilevel = 0;
  while (tracksLevel.at(ilevel).size() > 0) {
    int ntracksL = tracksLevel.at(ilevel).size();

    trackIdsGen0.resize(ilevel + 1);
    tracksLevel.resize(ilevel + 2);

    for (int i = 0; i < ntracksL; i++) {
      if (!(tracksLevel.at(ilevel)[i].Id() > last_track_id_))
        trackIdsGen0.at(ilevel).push_back(tracksLevel.at(ilevel)[i].Id());
      else {
        for (int idaughter = 0; idaughter < tracksLevel.at(ilevel)[i].DaughterIds().size(); idaughter++) {
          int daughterId = tracksLevel.at(ilevel)[i].DaughterIds().at(idaughter);
          int index = id2index_.find(daughterId)->second;
          tracksLevel.at(ilevel + 1).push_back(tracks_[index]);
        }
      }
    }
    ilevel++;
  }

  ilevel = 0;
  while (trackIdsGen0.size() > ilevel + 1) {
    if (trackIdsGen0.at(ilevel + 1).size() == 0)
      return false;
    for (int itrack = 0; itrack < trackIdsGen0.at(ilevel).size(); itrack++) {
      int jlevel = ilevel + 1;
      while (trackIdsGen0.size() > jlevel) {
        for (int jtrack = 0; jtrack < trackIdsGen0.at(jlevel).size(); jtrack++) {
          if (trackIdsGen0.at(ilevel)[itrack] == trackIdsGen0.at(jlevel)[jtrack])
            return false;
        }
        jlevel++;
      }
    }
    ilevel++;
  }
  return true;
}

void SimpleFinder::ReconstructDecay(const Decay& decay) {

  std::vector<std::vector<int>> indexes{};
  std::vector<Pdg_t> pdgs{};
  for (const auto& daughter : decay.GetDaughters()) {
    indexes.emplace_back(GetGoodDaughtersIndexes(daughter));
    pdgs.emplace_back(daughter.GetPdgHypo());
  }

  for (auto index_1 : indexes.at(0)) {

    std::array<KFParticle, 3> track;

    track.at(0) = GetTrack(index_1);
    if (std::abs(pdgs.at(0)) == 1000020030 || std::abs(pdgs.at(0)) == 1000020040) {
      int charge = (int) track.at(0).Q();
      charge *= 2;
      track.at(0).Q() = charge;
      track.at(0).Px() *= std::abs(charge);
      track.at(0).Py() *= std::abs(charge);
      track.at(0).Pz() *= std::abs(charge);
    }

    for (auto index_2 : indexes.at(1)) {
      track.at(1) = GetTrack(index_2);
      if (std::abs(pdgs.at(1)) == 1000020030 || std::abs(pdgs.at(1)) == 1000020040
          || std::abs(track.at(1).GetPDG()) == 1000020030 || std::abs(track.at(1).GetPDG()) == 1000020040) {
        int charge = (int) track.at(1).Q();
        charge *= 2;
        track.at(1).Q() = charge;
        track.at(1).Px() *= std::abs(charge);
        track.at(1).Py() *= std::abs(charge);
        track.at(1).Pz() *= std::abs(charge);
      }

      CalculateParamsInPCA(track.at(0), pdgs.at(0), track.at(1), pdgs.at(1));
      if (!IsGoodPair(track.at(0), track.at(1), decay)) continue;

      KFParticleSIMD kf_mother = ConstructMother({track.at(0), track.at(1)}, pdgs, decay, false);
      if (decay.GetNDaughters() == 2) {

        int id_mother = 0;
        if (!IsGoodMother(kf_mother, decay.GetMother(), id_mother)) continue;
        if (!IsMotherFromPV(kf_mother, decay.GetMother(), id_mother)) continue;
        if (!IsGoodDecayLength(kf_mother, decay.GetMother())) continue;
        if (!IsGoodCos(kf_mother, decay)) continue;
        if (!IsGoodMotherMass(kf_mother, decay.GetMother())) continue;

        FillDaughtersInfo({track.at(0), track.at(1)}, pdgs);
        SaveParticle(kf_mother, decay);

      } else if (decay.GetNDaughters() == 3) {

        int id_mother = 1;

        if (!IsGoodMother(kf_mother, decay.GetMother(), id_mother)) continue;
        if (IsMotherFromPV(kf_mother, decay.GetMother(), id_mother)) continue;

        for (auto index_3 : indexes.at(2)) {

          track.at(2) = GetTrack(index_3);
          if (std::abs(pdgs.at(2)) == 1000020030 || std::abs(pdgs.at(2)) == 1000020040) {
            int charge = (int) track.at(2).Q();
            charge *= 2;
            track.at(2).Q() = charge;
            track.at(2).Px() *= std::abs(charge);
            track.at(2).Py() *= std::abs(charge);
            track.at(2).Pz() *= std::abs(charge);
          }

          CalculateParamsInSV(track, pdgs);

          if (!IsGoodThree(track.at(0), track.at(1), track.at(2), decay)) continue;

          id_mother = 2;
          kf_mother = ConstructMother({track.at(0), track.at(2)}, {pdgs.at(0), pdgs.at(2)}, decay);
          if (!IsGoodMother(kf_mother, decay.GetMother(), id_mother)) continue;
          if (IsMotherFromPV(kf_mother, decay.GetMother(), id_mother)) continue;

          id_mother = 3;
          kf_mother = ConstructMother({track.at(1), track.at(2)}, {pdgs.at(1), pdgs.at(2)}, decay);
          if (!IsGoodMother(kf_mother, decay.GetMother(), id_mother)) continue;
          if (IsMotherFromPV(kf_mother, decay.GetMother(), id_mother)) continue;

          id_mother = 0;
          kf_mother = ConstructMother({track.at(0), track.at(1), track.at(2)}, pdgs, decay);
          if (!IsGoodMother(kf_mother, decay.GetMother(), id_mother)) continue;
          if (!IsMotherFromPV(kf_mother, decay.GetMother(), id_mother)) continue;
          if (!IsGoodDecayLength(kf_mother, decay.GetMother())) continue;
          if (!IsGoodCos(kf_mother, decay)) continue;
          if (!IsGoodMotherMass(kf_mother, decay.GetMother())) continue;

          FillDaughtersInfo({track.at(0), track.at(1), track.at(2)}, pdgs);

          SaveParticle(kf_mother, decay);
        }
      } else {
        throw std::runtime_error("Number of daughters should be 2 or 3. Current number is " + std::to_string(decay.GetNDaughters()));
      }
    }
  }
}

void SimpleFinder::FillDaughtersInfo(const std::vector<KFParticle>& tracks, const std::vector<Pdg_t>& pdgs) {
  for (size_t i = 0; i < tracks.size(); ++i) {
    values_.chi2_prim[i] = CalculateChiToPrimaryVertex(tracks.at(i), pdgs.at(i));
  }
}

void SimpleFinder::FindMixedBackground() {
  // Track-level event mixing (afkf-style): replace the current event's `mix_pdg`
  // tracks with the `mix_pdg` tracks of previous events, translated to the current PV.
  // Every candidate built in this pass is combinatorial background.
  if (!IsSupportedTrackPdg(mix_pdg_) || mix_buffer_size_ <= 0)
    return;

  // Select the event-mixing bin for the current event (mirrors AFKF's
  // centrality-binned buffer). Events only mix with the last `mix_buffer_size_`
  // events of the *same* centrality bin.
  int cent_bin;
  switch (centrality_source_) {
    case InputContainer::CentSrc::kField:
      cent_bin = CentralityAT(current_centrality_);
      break;
    case InputContainer::CentSrc::kNTacks:
      cent_bin = CentralityNTacks(current_n_good_tracks_);
      break;
    case InputContainer::CentSrc::kFixed8:
      // AFKF's actual behavior for non-"aTree" inputs: fCentrality is hard-coded
      // to 8 at the start of every event, so the mixing buffer is never binned by
      // track count. All events mix within the single (bin-8) buffer.
      cent_bin = 8;
      break;
    case InputContainer::CentSrc::kNone:
    default:
      cent_bin = 0;
      break;
  }
  if (cent_bin < 0 || cent_bin >= kMixCentralityBins) cent_bin = 0;
  auto& buffer = mixed_events_buffer_[cent_bin];

  // Save the current event reconstruction state.
  const std::vector<KFParticle> saved_tracks = tracks_;
  const std::map<Pdg_t, std::vector<int>> saved_indexes = indexes_;
  const std::map<int, int> saved_id2index = id2index_;
  const int saved_last_track_id = last_track_id_;
  const int saved_current_candidate_id = current_candidate_id_;

  // Build the combined track set:
  //   - all current event input tracks except the `mix_pdg` species;
  //   - the buffered `mix_pdg` tracks of previous events (same centrality bin),
  //     translated to the current PV.
  // NOTE: `tracks_` has been extended with reconstructed mother candidates during the
  // normal pass, so the *input* tracks are taken from the `input_tracks_` snapshot.
  std::vector<KFParticle> mixed_tracks;
  mixed_tracks.reserve(input_tracks_.size());

  for (const auto& track : input_tracks_) {
    if (track.GetPDG() != mix_pdg_)
      mixed_tracks.push_back(track);
  }
  for (const MixEvent& mix_event : buffer) {
    for (const KFParticle& mix_track : mix_event.mix_tracks)
      mixed_tracks.push_back(TranslateMixTrack(mix_track, mix_event.pv, prim_vx_));
  }

  // Re-id contiguously so that the generation boundary (last_track_id_) marks all
  // input tracks as generation 0, independent of their original event.
  tracks_ = std::move(mixed_tracks);
  for (size_t i = 0; i < tracks_.size(); ++i)
    tracks_.at(i).SetId(static_cast<int>(i));
  last_track_id_ = static_cast<int>(tracks_.size()) - 1;
  current_candidate_id_ = 0;

  is_mixed_pass_ = true;
  for (const auto& decay : decays_) {
    InitTrackIndexesMap();
    ReconstructDecay(decay);
  }
  is_mixed_pass_ = false;

  // Restore the current event reconstruction state.
  tracks_ = saved_tracks;
  indexes_ = saved_indexes;
  id2index_ = saved_id2index;
  last_track_id_ = saved_last_track_id;
  current_candidate_id_ = saved_current_candidate_id;

  // Buffer the current event's `mix_pdg` tracks for future mixing (ring buffer of
  // size mix_buffer_size_ per centrality bin).
  MixEvent mix_event;
  mix_event.pv = prim_vx_;
  for (const auto& track : input_tracks_) {
    if (track.GetPDG() == mix_pdg_)
      mix_event.mix_tracks.push_back(track);
  }
  if (!mix_event.mix_tracks.empty()) {
    buffer.push_back(std::move(mix_event));
    while (static_cast<int>(buffer.size()) > mix_buffer_size_)
      buffer.pop_front();
  }
}

KFParticle SimpleFinder::TranslateMixTrack(const KFParticle& track,
                                           const KFVertex& source_pv,
                                           const KFVertex& target_pv) const {
  KFParticle translated = track;
  translated.X() += target_pv.GetX() - source_pv.GetX();
  translated.Y() += target_pv.GetY() - source_pv.GetY();
  translated.Z() += target_pv.GetZ() - source_pv.GetZ();
  translated.SetFieldCoeff(translated.GetFieldCoeff()[9] + target_pv.GetZ() - source_pv.GetZ(), 9);
  return translated;
}

bool SimpleFinder::IsSupportedTrackPdg(int pdg) const {
  switch (std::abs(pdg)) {
    case 1:   // charge-based PDG (pid_mode 0 sets the track PDG to the charge, e.g. -1 for a pion)
    case 11:
    case 13:
    case 211:
    case 321:
    case 2212:
    case 1000010020:
    case 1000010029:
    case 1000010030:
    case 1000020030:
    case 1000020040:
      return true;
    default:
      return false;
  }
}
