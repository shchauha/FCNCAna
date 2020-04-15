#ifndef YEAR2016ALL
#define YEAR2016ALL

namespace y2016 {
#include "fake_rates.h"
namespace newbins {
#include "new_fake_rates.h"
}
#include "flip_rates.h"
#include "lepton_sf.h"
#include "puw.h"
#include "trigger_sf.h"
#include "isr_weights.h"
#include "fastsim_trigger_eff.h"
#include "fastsim_lepton_sf.h"
float getLumi() { return 35.922; }
}

#endif
