#!/usr/bin/env bash

branches="jets,lep2_closeJet,lep2_p4_gen,lep1_p4_gen,lep1_closeJet,jets_disc,jets_unc,jets_undoJEC,jets_JEC,filename,jets_flavor,btags,lep1_dZ,lep2_dZ,lep1_dxyPV,lep2_dxyPV,rawmetPhi,lep1_ip3d,lep2_ip3d,calometPhi,calomet,lep1_ip3d_err,lep2_d0_err,lep2_ip3d_err,lep1_d0_err,genht,weight,btags_cdisc,btags_unc,btags_undoJEC,btags_JEC,btags_disc,btags_flavor,lep1_iso,extragenbmoms,lep2_iso,weight_lepsf,all_leps_mc3idx,genps_id_mother,all_leps_veto,btags_effpt,all_leps_id,bjet_type,btags_eff,genps_id_grandma,genps_idx_mother,pdfweights,trijet_discs,trijet_frachigh,trijet_njetsnonb,trijet_leadingdisc,trijet_subleadingdisc,trijet_meandisc,trijet_numhigh,trijet_njetsb,jets_disc_dn,jets_disc_up,genps_status,genps_id,all_leps_p4,jets_jec_up,jets_jec_dn,bjets_jec_up,bjets_jec_dn,genps_p4,jet_type,jets_bsf,btags_sf,weight_lepsf2,weight_lepsf1,weight_pu,ngenjets,triggers,ngenjets30,lep2_mc3idx,lep1_mc3idx"
# cut="(hyp_class!=4 || is_real_data) && (njets >= 2)"
# skimdir="skim"

cut="(hyp_class!=4 || is_real_data) && (njets >= 2 || njets_unc_up >= 2 || njets_unc_dn >= 2 || njets_JER_up >= 2 || njets_JER_dn >= 2)"
skimdir="skimfix"

function make_skims() {
    basedir=$1
    mkdir -p ${basedir}/${skimdir}
    for sample in $(ls -1 ${basedir}/*.root); do
        sample=$(basename $sample);
        outname=${basedir}/${skimdir}/${sample}
        [ -e $outname ] && {
            echo "$outname already exists, skipping"
            continue
        }
        echo skim.py -t t -c "${cut}" "${basedir}/${sample}" -o ${outname} --flipbranches -b "${branches}"
        skim.py -t t -c "${cut}" "${basedir}/${sample}" -o ${outname} --flipbranches -b "${branches}"
    done

    # Now combine the two separate stop mass ranges into a single sample for convenience
    if [ -e ${basedir}/${skimdir}/T6TTHZ_m1.root -a -e ${basedir}/${skimdir}/T6TTHZ_m2.root ]; then
        if [ ! -e ${basedir}/${skimdir}/T6TTHZ_both.root ]; then
            echo hadd -k -f ${basedir}/${skimdir}/T6TTHZ_both.root ${basedir}/${skimdir}/T6TTHZ_m1.root ${basedir}/${skimdir}/T6TTHZ_m2.root 
            hadd -k -f ${basedir}/${skimdir}/T6TTHZ_both.root ${basedir}/${skimdir}/T6TTHZ_m1.root ${basedir}/${skimdir}/T6TTHZ_m2.root 
        else
            echo "${basedir}/${skimdir}/T6TTHZ_both.root already exists, not combining the two subsamples"
        fi
    fi
    if [ -e ${basedir}/${skimdir}/T5TTCC_m1.root -a -e ${basedir}/${skimdir}/T5TTCC_m2.root ]; then
        if [ ! -e ${basedir}/${skimdir}/T5TTCC_both.root ]; then
            echo hadd -k -f ${basedir}/${skimdir}/T5TTCC_both.root ${basedir}/${skimdir}/T5TTCC_m1.root ${basedir}/${skimdir}/T5TTCC_m2.root 
            hadd -k -f ${basedir}/${skimdir}/T5TTCC_both.root ${basedir}/${skimdir}/T5TTCC_m1.root ${basedir}/${skimdir}/T5TTCC_m2.root 
        else
            echo "${basedir}/${skimdir}/T5TTCC_both.root already exists, not combining the two subsamples"
        fi
    fi

}

mkdir -p logs

# make_skims "/nfs-7/userdata/namin/tupler_babies/merged/FT/v3.24/output/year_2016/" >& logs/skim_log_2016.txt &
# # make_skims "/nfs-7/userdata/namin/tupler_babies/merged/FT/v3.24/output/year_2017/" >& logs/skim_log_2017.txt &
# make_skims "/nfs-7/userdata/namin/tupler_babies/merged/FT/v3.24_newfixfsjecs/output/year_2017/" >& logs/skim_log_2017fix.txt &
# make_skims "/nfs-7/userdata/namin/tupler_babies/merged/FT/v3.28//output/year_2018/" >& logs/skim_log_2018.txt &

# make_skims "/nfs-7/userdata/namin/tupler_babies/merged/FT/v3.30//output/year_2016/" >& logs/skim_log_2016.txt &
# make_skims "/nfs-7/userdata/namin/tupler_babies/merged/FT/v3.30//output/year_2017/" >& logs/skim_log_2017.txt &
# make_skims "/nfs-7/userdata/namin/tupler_babies/merged/FT/v3.30//output/year_2018/" >& logs/skim_log_2018.txt &
# make_skims "/nfs-7/userdata/namin/tupler_babies/merged/FT/v3.29//output/year_2016_94x//" >& logs/skim_log_2016_94x.txt &

# # make_skims "/nfs-7/userdata/namin/tupler_babies/merged/FT/v3.30_hhat/output/year_2016/" >& logs/skim_log_2016.txt &
# make_skims "/nfs-7/userdata/namin/tupler_babies/merged/FT/v3.30_hhat/output/year_2017/" >& logs/skim_log_2017.txt &
# make_skims "/nfs-7/userdata/namin/tupler_babies/merged/FT/v3.30_hhat/output/year_2018/" >& logs/skim_log_2018.txt &

make_skims "/nfs-7/userdata/namin/tupler_babies/merged/FT/v3.31//output/year_2016_94x//" >& logs/skim_log_2016_94x.txt &
make_skims "/nfs-7/userdata/namin/tupler_babies/merged/FT/v3.31//output/year_2017/" >& logs/skim_log_2017.txt &
make_skims "/nfs-7/userdata/namin/tupler_babies/merged/FT/v3.31//output/year_2018/" >& logs/skim_log_2018.txt &
# make_skims "/nfs-7/userdata/namin/tupler_babies/merged/FT/v3.31_jec15//output/year_2018/" >& logs/skim_log_2018.txt &

wait
