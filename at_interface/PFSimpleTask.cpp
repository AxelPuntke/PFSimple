#include "PFSimpleTask.hpp"

#include "ConverterIn.hpp"
#include "ConverterOut.hpp"

void PFSimpleTask::Init() {
  pf_simple_ = new SimpleFinder();
  pf_simple_->SetDecays(decays_);
  pf_simple_->SetIsEventMixing(is_event_mixing_);
  pf_simple_->SetMixPdg(mix_pdg_);
  pf_simple_->SetMixBufferSize(mix_buffer_size_);
}

void PFSimpleTask::Exec() {
  pf_simple_->Init(in_task_->GetInputContainer());
  pf_simple_->FindParticles();
}

const std::vector<OutputContainer>& PFSimpleTask::GetMixedCandidates() const {
  return pf_simple_->GetMixedCandidates();
}

void PFSimpleTask::Finish() {
  delete pf_simple_;
}
