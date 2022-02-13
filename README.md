# A comparison of topologically associating domain callers based on Hi-C data

This Github repository is used for reproducing the results in the manuscript *A comparison of topologically associating domain callers based on Hi-C data*. If you have any questions for this Github repository, please feel free to contact with me emre.sefer@ozyegin.edu.tr

## prepare Hi-C data (scripts are in the prepare_data folder)

### prepare the Hi-C data for Section 3.1

**step 1.**    Download merged sequencing reads from Hi-C experiments

	./download_raw_data.sh
	
**step 2.**    Down-sample reads from the pooled Hi-C data and generate KR normalized HI-C matrices using Juicertools 

	./downsample.sh

### prepare the Hi-C data for Section 3.2 and 3.3

**step 1.**    Download raw hic files for the human GM12878, IMR90, and K562 cell lines 

	./download_raw_hic.sh
	
**step 2.**    Generate KR normalized HI-C matrices using Juicertools 

	./dump_data_from_hic.sh

The synthetic Hi-C data used in this manuscript is deposited at work/simulate_data folder and work/benchmarks.

## call TADs (scripts are in the work folder)

The preprocessed data for Section 3.1 could be deposited at work/GM12878_downsample_diff_reso folder.

The preprocessed data for Section 3.2 and 3.3 could be deposited at work/Rao folder.
Directories structure in work/Rao folder

---Rao/

------HIC001/

------HIC002/

------...

Bash scripts in the folder of each TAD caller (such as work/Armatus) shall call TADs from the Hi-C data. The executive file path or command statement should be updated according to your actual situation. It will spend a lot of time for some methods to call TADs, such as TADtree. 

## reproduce the results in manuscript

**step 1.** Decompress the zip files in work folder

**step 2.** Convert the TADs identified by various TAD callers to TADs with a uniform format

	./extract_all.sh
	./extract_diff.sh

For convenience, the processed TADs with uniform format would be deposited at work/all_TADs/loci and work/all_TADs/bin folders, if the all_TADs.zip is decompressed.

**step 3.** Generate figures and tables in the manuscript 

Scripts for Figure1,2, Supplementary Figure1,2,3,4,7 and Supplementary 7,8,9 in work/Fig_Tab folder can directly plot corresponding figures and tables.

Run_TADadjRsquared.sh, Fig3.py, and Table5_Supp_Table7_8.py in work/Fig_Tab folder can generate corresponding figure and tables.

run_rep_size_contacts.sh and Fig4.py in work/Fig_Tab folder can generate Figure4.

run_enrichment.sh and Table2_3_4_Supp_3_4.py in work/Fig_Tab folder can generate corresponding tables.
