The files contained in this directory are used to take cleaned reads and convert them into an OTU table using Qiime and the Silva 132 database.
The whole process begins with the qiime-pipe.sh script, which starts the silva-qiime-pipe.sh script for each of the cleaned reads found in a sample. 
silva-qiime-pipe.sh then uses the pick_otus.py script from Qiime to build the OTU data which in then converted into a .biom file containing the OTU data for one sample. 
merge-biom-files.sh combinds the individual sample .biom files into one large .biom file containing all OTU information.