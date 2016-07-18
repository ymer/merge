How to use:


1. Install GWF

Grid WorkFlow (GWF) is created by Thomas Mailund at BiRC.

Go to an appropriate location for installation and run the following 3 commands:

gitproxy clone https://github.com/mailund/gwf.git
cd gwf
python2.7 setup.py install --user

GWF should now be installed in this folder, with the sub-path .local/bin/path.


2. Prepare for running

The newest version is found in IGDATA/faststorage/userdata/ijpal/merge_[date]. Copy this folder.

Load anaconda:
source /com/extra/Anaconda-Python/2.2.0-2.7/load.sh

If there is no file_locations file in the files folder (or not one with the same name as the data set you are interested in merging), create one using prep.py


3. Run the merging

The merging is run with run.py. For help on how to run the program, type run.py --help.

example command:
python run.py --gwf /home/ijpal/.local/bin/gwf --date 06072016 --vers phase3_woBLACKLIST_newQC --chunks one --genoprob 0.0 --run_check


4. Results

The merging results are found in output/[merging_folder]. The merging folder has aa name that includes the version, date, genotype cutoff and individuals included. Inside this folder there is a readme.txt file describing the content of the output folders.




