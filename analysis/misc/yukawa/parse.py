def get_stuff(fname):
    d = {}
    which = None
    for line in open(fname,"r").readlines():
        line = line.strip()
        if line.startswith("yukawa"):
            which = int(line.split()[-1])
            print which
            d[which] = {}

        if line.startswith("Obs UL:"):
            limit = float(line.split(":")[-1].split()[0])
            d[which]["limit"] = limit
        if line.startswith("Mu:"):
            parts = line.split(":")[-1]
            central = float(parts.split()[0])
            up = float(parts.split()[1].replace("+","").replace("(",""))
            down = float(parts.split()[2].replace("-","").replace(")",""))
            d[which]["central"] = 9.2*central
            d[which]["up"] = 9.2*up
            d[which]["down"] = 9.2*down


    # print d
    keys = sorted(d.keys())
    ktvals = sorted(map(lambda x: 0.1*x, keys))
    ul = map(lambda x: d[x]["limit"], keys)
    central = map(lambda x: d[x]["central"], keys)
    ups = map(lambda x: d[x]["up"], keys)
    downs = map(lambda x: d[x]["down"], keys)
    # print ktvals
    # print keys
    # print ul
    # print central, ups, downs

    return ktvals, ul, central, ups, downs

if __name__ == "__main__":
    get_stuff()
