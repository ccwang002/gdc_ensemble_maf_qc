# Quality check of GDC ensemble masked MAFs for CPTAC3 samples


## Data requirement
The repository requires the following files:

    external_data
    ├── Batch1through8_sample_attributes.csv.gz
    └── cptac3_aliquot_ensemble_masked_mafs.tar.gz


    mkdir external_data/disease_wg
    cd external_data/disease_wg
    rsync -a \
        denali:/diskmnt/Datasets/CPTAC/PGDAC/UCEC/Baylor_DataFreeze_V2.1/UCEC_somatic_mutation_site_level_V2.1.maf \
        ucec_hg38_all_samples.maf

    rsync -az --info=progress2 \
        denali:/diskmnt/Datasets/CPTAC/PGDAC/CCRCC/Umich_data_freeze_v_from20190125/mutation/somatic_hg38_v2.0_20181130/ \
        ccrcc_hg38/

    rsync -az --info=progress2 \
        vw3.gsc.wustl.edu:/gscmnt/gc2521/dinglab/mwyczalk/CPTAC3.Y1-submit/LUAD.all.Somatic/CPTAC3_WashU_hg38/LUAD/WXS_Somatic.LUAD.all.hg38_v1.3_20180125/ \
        luad_hg38/

    rsync -az \
        vw3.gsc.wustl.edu:/gscmnt/gc2521/dinglab/mwyczalk/CPTAC3.Y1-submit/LUAD.20190422/CPTAC3_WashU_hg38/LUAD/WXS_Somatic.LUAD.1-missing.v1.3.hg38.20190422/C3N-01024.Somatic.WXS.maf \
        luad_hg38/

    rsync -az \
        vw3.gsc.wustl.edu:/gscmnt/gc2518/dinglab/cptac3/hg38/luad/somatic/worklog/LUAD.Somatic.050919.mnp.annot.maf \
        luad_hg38_all_samples.song_wip.maf


## Setup

    python3 scripts/add_maf_to_db.py \
        processed_data/all_mutations.sqlite \
        disease_wg \
        external_data/disease_wg/luad_hg38_all_samples.song_wip.maf \
        external_data/disease_wg/ucec_hg38_all_samples.maf \
        external_data/disease_wg/ccrcc_hg38/*.maf
