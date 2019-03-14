#include "MuonTrackAnalyzer.C"

void runMuonTrackAnalyzer(bool Propagation, bool GenPropagation){
  MuonTrackAnalyzer t(Propagation, GenPropagation);
  t.Loop();
}
