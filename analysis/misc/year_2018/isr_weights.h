float isrWeight(int nisrmatch, int sample) {
    // see instructions in ../isr_weights/calc.py
    float scale = 1.;
    if (sample == 1) { // ttW -- 0 extra partons
        scale = 0.9933;
        if (nisrmatch == 0) { return scale*1.037; }
        if (nisrmatch == 1) { return scale*0.913; }
        if (nisrmatch == 2) { return scale*1.026; }
        if (nisrmatch >= 3) { return scale*1.236; }
    } else if (sample == 2) { // ttZ -- 1 extra parton
        scale = 0.9933;
        if (nisrmatch == 0) { return scale*1.061; }
        if (nisrmatch == 1) { return scale*0.904; }
        if (nisrmatch == 2) { return scale*0.999; }
        if (nisrmatch == 3) { return scale*1.180; }
        if (nisrmatch >= 4) { return scale*1.304; }
    } else if (sample == 10) { // amcatnlo ttbar sample
        scale = 0.9854;
        if (nisrmatch == 0) { return scale*1.073; }
        if (nisrmatch == 1) { return scale*0.912; }
        if (nisrmatch == 2) { return scale*0.944; }
        if (nisrmatch == 3) { return scale*1.056; }
        if (nisrmatch >= 4) { return scale*1.422; }
    }
    return scale;
}
