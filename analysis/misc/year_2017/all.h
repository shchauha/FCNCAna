#ifndef YEAR2017ALL
#define YEAR2017ALL

namespace y2017 {
#include "fake_rates.h"
namespace allpt {
#include "new_fake_rates.h"
}
#include "flip_rates.h"
#include "lepton_sf.h"
#include "puw.h"
#include "trigger_sf.h"
#include "isr_weights.h"
#include "fastsim_lepton_sf.h"
float getLumi() { return 41.53; }
}

#endif
