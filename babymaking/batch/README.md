## Making babies
<img src="http://i.imgur.com/6BFuXOh.png" width="500">

### So you want to develop the babymaker?
```bash
# Copy over CORE and stuff
./copy_files.sh
# Edit the babymaker source files or whatever
# Compile with optimized flags
make opt -j10
# Test? Or don't, actually. We don't make mistakes.
./test_babymaker.sh
# Make a tar file to ship to the worker nodes
./make_tar.sh
```

### So you want to make some babies?
I know this guy -- great guy -- the best. He uses condor to deliver his babies instead of storks.
I'm going to let you in on a secret...

#### (new) ProjectMetis
```bash
# Edit submetis.py with samples and logic you want and also increment tag
# Open up a screen before proceeding
# Get environment
cd ProjectMetis
. setup.sh
cd -
# Submit
python submetis.py
```

#### (old) AutoTwopler
```bash
# Edit four_top.py with samples you want and increment tag near the bottom
# Edit ducks.py to build instructions with appropriate samples
# Open up a screen before proceeding
# Get environment
cd NtupleTools/AutoTwopler/
# edit `campaign` in `params.py` to be "FT_babymaker"
. setup.sh
cd ../../
# Submit
python ducks.py
```

### Tips and tricks
That pesky scale1fb text file makes it annoying to just run on an arbitrary file. Sometimes
we just want to check something simple. Luckily we can bypass all that with
```bash
noscale1fb=1 ./main.exe /hadoop/cms/store/user/namin/run2_moriond17_cms4/ProjectMetis/tttt-LOytscan_kt-2p1_v2_CMS4_CMS4_V00-00-02_2017Sep27/merged_ntuple_1.root
```
Also, can skip opposite sign events in MC with `noos=1`.
