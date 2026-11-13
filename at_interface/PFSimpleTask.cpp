#include "PFSimpleTask.hpp"

#include "ConverterIn.hpp"
#include "ConverterOut.hpp"

#include "AnalysisTree/TaskManager.hpp"

void PFSimpleTask::Init() {
  pf_simple_ = new SimpleFinder();
  pf_simple_->SetDecays(decays_);
  pf_simple_->SetIsEventMixing(is_event_mixing_);
  pf_simple_->SetMixPdg(mix_pdg_);
  pf_simple_->SetMixBufferSize(mix_buffer_size_);
  pf_simple_->SetChargeBasedPdg(charge_based_pdg_);
  pf_simple_->SetCentralitySource(centrality_source_);
  pf_simple_->SetAcceptFragments(accept_fragments_);
  pf_simple_->SetChi2PrimAfkf(chi2_prim_afkf_);
  in_task_->SetCentralitySource(centrality_source_);
  in_task_->SetCentralityFieldName(centrality_field_name_);
}

void PFSimpleTask::Exec() {
  // Record the entry number of the current event so the finder can tag mixed
  // candidates with the source event of each daughter.
  auto* man = AnalysisTree::TaskManager::GetInstance();
  auto* chain = man->GetChain();
  pf_simple_->SetCurrentEventIndex(static_cast<int>(chain->GetReadEntry()));
  pf_simple_->Init(in_task_->GetInputContainer());
  pf_simple_->FindParticles();
}

const std::vector<OutputContainer>& PFSimpleTask::GetMixedCandidates() const {
  return pf_simple_->GetMixedCandidates();
}

void PFSimpleTask::Finish() {
  delete pf_simple_;
}
