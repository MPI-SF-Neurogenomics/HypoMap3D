# 05 Clustering

In this module, we perform our hierarchical multi-level clustering via iterative decomposition of higher clustering resolutions. 

Starting point are 'major features' that are derived from propagated ABC class id labels (--> GABA, Glut, NN classes).

Afterwards, each of neuronal and non-neuronal classes are decomposed separately to use decomposition parameters yielding best results.

After clustering, furter low quality cells/clusters were detected which were removed in 6_posthoc_QC.ipynb, after which some data needed updating. 
When running the pipeline yourself, you may want to choose omit 


### Notes on usage

For the scVI parts, the script with "submit" in its name will submit a job to SLURM while the other one will run scVI in the terminal. 

For the scripts calculating marker genes the same principle applies: the .py script will run a process in the terminal while the submit scripts will use SLURM.

Running either of the versions will work but using SLURM is heavily encouraged.