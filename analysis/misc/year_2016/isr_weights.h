float isrWeight(int nisrmatch, int sample) {
    // see instructions in ../isr_weights/calc.py
    float scale = 1.;
    if (sample == 1) { // ttW
        scale = 1.09063;
    } else if (sample == 2) { // ttZ
        scale = 1.10315;
    } else if (sample == 10) { // ttbar sample
        scale = 1.0643;
    }
    if (nisrmatch == 0) { return scale*1.00; }
    if (nisrmatch == 1) { return scale*0.86; }
    if (nisrmatch == 2) { return scale*0.79; }
    if (nisrmatch == 3) { return scale*0.75; }
    if (nisrmatch == 4) { return scale*0.77; }
    if (nisrmatch == 5) { return scale*0.74; }
    if (nisrmatch == 6) { return scale*0.73; }
    if (nisrmatch == 7) { return scale*0.56; }
    return scale;
}
