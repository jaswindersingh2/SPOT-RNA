SPOT-RNA: RNA Secondary Structure Prediction using an Ensemble of Two-dimensional Recurrent and Residual Convolutional Neural Networks and Transfer Learning.
====

OVERVIEW:
-----
SPOT-RNA is a single sequence-based RNA secondary structure predictor. It used an ensemble of deep learning methods (ResNets and BiLSTM) to infer the base-pair probability of nucleotides within the sequence. SPOT-RNA was initially trained on non-redundant bpRNA dataset in which secondary structure was derived using comparative sequence analysis. After initial training, transfer learning was used to further train SPOT-RNA on high resolution non-redundant structured RNAs from Protein Data Bank (PDB). SPOT-RNA is also available as web-server at http://sparks-lab.org/jaswinder/server/SPOT-RNA/. The web-server provides an arc diagram and a 2D diagram of predicted RNA secondary
the structure through Visualization Applet for RNA (VARNA) tool along with a dot plot of SPOT-RNA predicted base-pair
probabilities.

|![](./SPOT-RNA-architecture.png)
|----|
| <p align="center"> <b>Figure 1:</b> The network layout of the SPOT-RNA, where L is the sequence length of a target RNA, Act. indicates the activation function, Norm. indicates the normalization function, and PreT indicates the pretrained (initial trained) models trained on the bpRNA dataset.|

USAGE:
-----
SPOT-RNA uses RNA sequence in fasta format as input and provide output in ct format and base-pair probability in text file.

Installation:
----
It is recommended to use a [virtual environment](http://virtualenvwrapper.readthedocs.io/en/latest/install.html) for installation.

Prerequisites:

* [TensorFlow (v1.12) ](https://www.tensorflow.org/install/) 
* [Python3](https://docs.python-guide.org/starting/install3/linux/)
* pandas
* numpy
* tqdm
* argparse
* six
* virtualenv or Anaconda (optional)

To install:

1. `git clone https://github.com/jaswindersingh2/SPOT-RNA.git`
2. `cd SPOT-RNA`
3. `wget 'https://www.dropbox.com/s/dsrcf460nbjqpxa/SPOT-RNA-models.tar.gz' || wget 'https://app.nihaocloud.com/f/fbf3315a91d542c0bdc2/?dl=1'`
4. `tar -xvzf SPOT-RNA-models.tar.gz && rm SPOT-RNA-models.tar.gz`
5. `virtualenv -p python3 venv && source ./venv/bin/activate || conda create -n venv python=3.6 anaconda && source activate venv` (optional)
6. `pip install -r requirements.txt`

How to Use the SPOT-RNA Scripts
-----

For single sequence:
```
python3 SPOT-RNA.py  --input_seq sample_inputs/single_seq.fasta  --output_dir 'outputs/'
```

For batch of sequences:
```
python3 SPOT-RNA.py  --input_seq sample_inputs/batch_seq.fasta  --output_dir 'outputs/'
```

and check the outputs against the SPOT-RNA files in sample directory. To specify running on GPU please set the --gpu argument.

Datasets Used For Training, Validation, and Testing
-----

The following datasets were used for Initial Training:
* [bpRNA](https://www.dropbox.com/s/w3kc4iro8ztbf3m/bpRNA_dataset.zip)[1]


The following datasets were used for Transfer Learning:
* [PDB](https://www.dropbox.com/s/rlr8n9r5mt456cd/PDB_dataset.zip)[2]

References
-----
[1] Padideh Danaee, Mason Rouches, Michelle Wiley, Dezhong Deng, Liang Huang, David Hendrix, bpRNA: large-scale automated annotation and analysis of RNA secondary structure, Nucleic Acids Research, Volume 46, Issue 11, 20 June 2018, Pages 5381–5394, https://doi.org/10.1093/nar/gky285

[2] H.M. Berman, J. Westbrook, Z. Feng, G. Gilliland, T.N. Bhat, H. Weissig, I.N. Shindyalov, P.E. Bourne.
(2000) The Protein Data Bank Nucleic Acids Research, 28: 235-242.

If you use this code for your research please cite the following papers:
-----
Singh, J., Hanson, J., Paliwal, and Zhou, Y., 2019. RNA Secondary Structure Prediction using an Ensemble of Two-dimensional Recurrent and Residual Convolutional Neural Networks and Transfer Learning (Under Review).
