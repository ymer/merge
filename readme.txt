How to use:


1. Install GWF

Grid WorkFlow (GWF) is created by Thomas Mailund at BiRC.

Go to an appropriate location for installation and run the following commands:

gitproxy clone https://github.com/mailund/gwf.git
cd gwf
python2.7 setup.py install --user

GWF should now be installed in this folder, with the sub-path .local/bin/path.


2. Prepare for running

Load anaconda:
source /com/extra/Anaconda-Python/2.2.0-2.7/load.sh

Create new file_locations files using prep.py. This contains the paths of the input dosage and info files.


3. Run the merging

The merging is run with run.py. For help on how to run the program, type run.py --help.

example command:
python run.py --gwf /home/ijpal/.local/bin/gwf --date 06072016 --vers phase3 --chunks all --genoprob 0.0 --run_check


4. Results

The merging results are found in output/[merging_folder]. The merging folder has aa name that includes the version, date, genotype cutoff and individuals included. Inside this folder there is a readme.txt file describing the content of the output folders.




