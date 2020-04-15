float isrWeight(int nisrmatch, int sample) {
    // see instructions in ../isr_weights/calc.py
    float scale = 1.;
    if (sample == 1) { // ttW -- 0 extra partons
        scale = 0.985;
        if (nisrmatch == 0) { return scale*1.056; }
        if (nisrmatch == 1) { return scale*0.886; }
        if (nisrmatch == 2) { return scale*1.007; }
        if (nisrmatch >= 3) { return scale*1.417; }
    } else if (sample == 2) { // ttZ -- 1 extra parton
        scale = 0.995;
        if (nisrmatch == 0) { return scale*1.080; }
        if (nisrmatch == 1) { return scale*0.902; }
        if (nisrmatch == 2) { return scale*0.934; }
        if (nisrmatch == 3) { return scale*1.165; }
        if (nisrmatch >= 4) { return scale*1.231; }
    } else if (sample == 10) { // amcatnlo ttbar sample
        scale = 0.985;
        if (nisrmatch == 0) { return scale*1.072; }
        if (nisrmatch == 1) { return scale*0.913; }
        if (nisrmatch == 2) { return scale*0.943; }
        if (nisrmatch == 3) { return scale*1.046; }
        if (nisrmatch >= 4) { return scale*1.366; }
    }
    return scale;
}
