#ifndef XSEC
#define XSEC

struct Xsec { float xsec; float percErr; };
Xsec stop_xsec(int mStop);
Xsec go_xsec(int mgl);

#endif
