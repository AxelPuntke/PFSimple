#ifndef KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_PFSIMPLETASK_H_
#define KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_PFSIMPLETASK_H_

#include "Decay.hpp"
#include "SimpleFinder.hpp"

#include "AnalysisTree/Task.hpp"

class ConverterIn;

class PFSimpleTask : public AnalysisTree::Task {

 public:
  PFSimpleTask() = default;
  ~PFSimpleTask() override = default;

  void Init() override;
  void Exec() override;
  void Finish() override;
  void SetInTask(ConverterIn* in_task) { in_task_ = in_task; }
  SimpleFinder* GetSimpleFinder() { return pf_simple_; }
  void SetDecays(const std::vector<Decay>& decays) { decays_ = decays; }
  void AddDecay(const Decay& decay) { decays_.emplace_back(decay); }

  ///< Enable/disable background estimation via event mixing.
  void SetIsEventMixing(bool is = true) { is_event_mixing_ = is; }
  ///< PDG of the track species replaced by tracks from previous events (track-level mixing).
  void SetMixPdg(int pdg) { mix_pdg_ = pdg; }
  ///< Number of previous events whose tracks are combined with the current event to estimate the background.
  void SetMixBufferSize(int n) { mix_buffer_size_ = n; }
  ///< Select the event-mixing centrality source (mirrors AFKF's tree-name-dependent behavior).
  void SetCentralitySource(InputContainer::CentSrc src) { centrality_source_ = src; }
  ///< Name of the event-header field holding the centrality value (default "centrality_tracks").
  void SetCentralityFieldName(const std::string& name) { centrality_field_name_ = name; }
  ///< Background candidates obtained by combining tracks from different events.
  const std::vector<OutputContainer>& GetMixedCandidates() const;

 protected:
  ConverterIn* in_task_{nullptr};
  SimpleFinder* pf_simple_{nullptr};
  std::vector<Decay> decays_{};
  bool is_event_mixing_{false};
  int mix_pdg_{-211};
  int mix_buffer_size_{5};
  InputContainer::CentSrc centrality_source_{InputContainer::CentSrc::kField};
  std::string centrality_field_name_{"centrality_tracks"};
};

#endif//KFPARTICLESIMPLE_ANALYSISTREEINTERFACE_PFSIMPLETASK_H_
