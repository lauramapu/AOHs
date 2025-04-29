## Brief explanation of scripts 

`/01-base_files`. This script should only be run once at the beggining because it generates the common files for all years: tropical forest mask (uploaded), tropical and terrestrial mammals SHP, habitat preferences and elevation ranges (uploaded), global SRTM, and translation table (uploaded).

`/02-process_per_year`. This script contains loops to process the desired years for all tropical and terrestrial mammals. With a vector named 'years' at the beggining of the script you can set the years of ESA-CCI to download and process base layers, and then filter AOHs.

`/02-in_slurm`. This script is the same as the previous ones but modified to run on SLURM (just setting the $STORE path to read/write big files).

`/gee_srtm_global`. This is the script used in Google Earth Engine to merge all SRTM (90m) and resample to 300m. It should return 3 files.

Other scripts are just trials.
