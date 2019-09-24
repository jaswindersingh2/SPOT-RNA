SPOT-RNA: RNA Secondary Structure Prediction using an Ensemble of Two-dimensional Recurrent and Residual Convolutional Neural Networks and Transfer Learning.
====

OVERVIEW
====
SPOT-RNA is a single sequence-based RNA secondary structure predictor. It used an ensemble of deep learning methods (ResNets and BiLSTM) (Figure 1) to infer the base-pair probability of nucleotides within the sequence. SPOT-RNA was initially trained on non-redundant bpRNA[1] dataset in which secondary structure was derived using comparative sequence analysis. After initial training, transfer learning was used to further train SPOT-RNA on high resolution non-redundant structured RNAs from Protein Data Bank (PDB)[2]. SPOT-RNA can predict all kind of base-pair including pseudoknots, non-canonical, lone, and triplets base-pairs. SPOT-RNA is also available as web-server at http://sparks-lab.org/jaswinder/server/SPOT-RNA/. The web-server provides an arc diagram and a 2D diagram of predicted RNA secondary the structure through Visualization Applet for RNA (VARNA)[3] tool along with a dot plot of SPOT-RNA predicted base-pair
probabilities.

|![](./SPOT-RNA-architecture.png)
|----|
| <p align="center"> <b>Figure 1:</b> The network layout of the SPOT-RNA, where L is the sequence length of a target RNA, Act. indicates the activation function, Norm. indicates the normalization function, and PreT indicates the pretrained (initial trained) models trained on the bpRNA dataset.|

SYSTEM REQUIREMENTS
====
Hardware Requirments:
----
SPOT-RNA predictor requires only a standard computer with around 16 GB RAM to support the in-memory operations for RNAs sequence length less than 500.

Software Requirments:
----
* [TensorFlow (v1.12) ](https://www.tensorflow.org/install/) 
* [Python3](https://docs.python-guide.org/starting/install3/linux/)
* pandas
* numpy
* tqdm
* argparse
* six
* virtualenv or Anaconda (optional)

SPOT-RNA has been tested on Ubuntu 16.04 operating system.

USAGE
====

Installation:
----
It is recommended to use a [virtual environment](http://virtualenvwrapper.readthedocs.io/en/latest/install.html) for installation.

To install:

1. `git clone https://github.com/jaswindersingh2/SPOT-RNA.git`
2. `cd SPOT-RNA`
3. `wget 'https://www.dropbox.com/s/dsrcf460nbjqpxa/SPOT-RNA-models.tar.gz' || wget 'https://app.nihaocloud.com/f/fbf3315a91d542c0bdc2/?dl=1'`
4. `tar -xvzf SPOT-RNA-models.tar.gz && rm SPOT-RNA-models.tar.gz`
5. `virtualenv -p python3 venv && source ./venv/bin/activate || conda create -n venv python=3.6 anaconda && source activate venv` (optional)
6. `pip install -r requirements.txt`

To run the SPOT-RNA Scripts
-----

**For single sequence:**
```
python3 SPOT-RNA.py  --inputs sample_inputs/single_seq.fasta  --outputs 'outputs/'
```
The outputcome of above command will be three files (.bpseq, .ct, and .prob) in 'outputs' folder. '.bpseq' and '.ct' file is standard format to represent RNA seondary structure. '.prob' file consists of base-pair probability of predicted secondary structure by SPOT-RNA which can be useful for plotting PR-curve and to check the confidence of predicted base-pair.

**For batch of sequences:**
```
python3 SPOT-RNA.py  --inputs sample_inputs/batch_seq.fasta  --outputs 'outputs/'
```

**To run on GPU:**
```
python3 SPOT-RNA.py  --inputs sample_inputs/batch_seq.fasta  --outputs 'outputs/' --GPU 0
```

**2D plots of predicted secondary structure:**

To get the 2D plots of SPOT-RNA output, VARNA[3] tool is used. To run this tool, please make sure Java plugin version 1.6 is installed in the system. Please refer to http://varna.lri.fr/ for more detailed information about this tool.
```
python3 SPOT-RNA.py  --inputs sample_inputs/single_seq.fasta  --outputs 'outputs/' --plots True
```
The output of above command will generate two additional files (arc plot and 2D plot of predicted secondary structure) along '.bpseq', '.ct', and '.prob' in 'outputs' folder.

**Secondary structure motifs form predicted structure:**

To get the secondary structure motifs like stem, helix, loops from predicted structure, SPOT-RNA used tool from bpRNA[1]. To run this script, please make sure 'Graph.pm' module (https://metacpan.org/pod/Graph) of perl is installed in the system.  Please refer to https://github.com/hendrixlab/bpRNA for more detailed information about this tool.
```
python3 SPOT-RNA.py  --inputs sample_inputs/single_seq.fasta  --outputs 'outputs/' --plots True --motifs True
```
The output of above command will generate one additional file '.st' in 'outputs' folder.
 
Datasets Used For Training, Validation, and Testing
====

The following datasets were used for Initial Training:
* [bpRNA](https://www.dropbox.com/s/w3kc4iro8ztbf3m/bpRNA_dataset.zip)[1]


The following datasets were used for Transfer Learning:
* [PDB](https://www.dropbox.com/s/rlr8n9r5mt456cd/PDB_dataset.zip)[2]

References
====
If you use SPOT-RNA for your research please cite the following papers:
----
Singh, J., Hanson, J., Paliwal, and Zhou, Y., 2019. SPOT-RNA: RNA Secondary Structure Prediction using an Ensemble of Two-dimensional Recurrent and Residual Convolutional Neural Networks and Transfer Learning (Under Review).

Other references:
----
[1] Padideh Danaee, Mason Rouches, Michelle Wiley, Dezhong Deng, Liang Huang, David Hendrix, bpRNA: large-scale automated annotation and analysis of RNA secondary structure, Nucleic Acids Research, Volume 46, Issue 11, 20 June 2018, Pages 5381–5394, https://doi.org/10.1093/nar/gky285

[2] H.M. Berman, J. Westbrook, Z. Feng, G. Gilliland, T.N. Bhat, H. Weissig, I.N. Shindyalov, P.E. Bourne.
(2000) The Protein Data Bank Nucleic Acids Research, 28: 235-242.

[3]  VARNA: Interactive drawing and editing of the RNA secondary structure Kévin Darty, Alain Denise and Yann Ponty Bioinformatics, pp. 1974-1975, Vol. 25, no. 15, 2009 


Licence
====
Mozilla Public License 2.0
