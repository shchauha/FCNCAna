#!/bin/bash
echo "input dir:" $1  

python cards.py $1 sr_sr_fcnc_in fcnc_hut sr 2016
python cards.py $1 sr_sr_fcnc_in fcnc_hut sr 2017
python cards.py $1 sr_sr_fcnc_in fcnc_hut sr 2018
python cards.py $1 sr_sr_fcnc_in fcnc_hut sr run2


python cards.py $1 ssbr_sr_fcnc_in fcnc_hut ssbr 2016
python cards.py $1 ssbr_sr_fcnc_in fcnc_hut ssbr 2017
python cards.py $1 ssbr_sr_fcnc_in fcnc_hut ssbr 2018
python cards.py $1 ssbr_sr_fcnc_in fcnc_hut ssbr run2

python cards.py $1 mlbr_sr_fcnc_in fcnc_hut mlbr 2016
python cards.py $1 mlbr_sr_fcnc_in fcnc_hut mlbr 2017
python cards.py $1 mlbr_sr_fcnc_in fcnc_hut mlbr 2018
python cards.py $1 mlbr_sr_fcnc_in fcnc_hut mlbr run2


## hct

python cards.py $1 sr_sr_fcnc_in fcnc_hct sr 2016
python cards.py $1 sr_sr_fcnc_in fcnc_hct sr 2017
python cards.py $1 sr_sr_fcnc_in fcnc_hct sr 2018
python cards.py $1 sr_sr_fcnc_in fcnc_hct sr run2


python cards.py $1 ssbr_sr_fcnc_in fcnc_hct ssbr 2016
python cards.py $1 ssbr_sr_fcnc_in fcnc_hct ssbr 2017
python cards.py $1 ssbr_sr_fcnc_in fcnc_hct ssbr 2018
python cards.py $1 ssbr_sr_fcnc_in fcnc_hct ssbr run2

python cards.py $1 mlbr_sr_fcnc_in fcnc_hct mlbr 2016
python cards.py $1 mlbr_sr_fcnc_in fcnc_hct mlbr 2017
python cards.py $1 mlbr_sr_fcnc_in fcnc_hct mlbr 2018
python cards.py $1 mlbr_sr_fcnc_in fcnc_hct mlbr run2

