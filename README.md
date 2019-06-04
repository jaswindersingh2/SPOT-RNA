SPOT-RNA: RNA Secondary Structure Prediction using an Ensemble of Two-dimensional Recurrent and Residual Convolutional Neural Networks and Transfer Learning.
====

USAGE:
-----
SPOT-RNA uses RNA sequence as input:

Model Availability:
-----

The models used for SPOT-RNA can be downloaded from the Dropbox [here](https://www.dropbox.com/s/dsrcf460nbjqpxa/SPOT-RNA-models.tar.gz) or Nihao cloud service [here](https://app.nihaocloud.com/f/882db8caf4be43ddaa04/?dl=1)


Installation
-----

It is recommended to use a [virtual environment](http://virtualenvwrapper.readthedocs.io/en/latest/install.html) for installation.

Prerequisites:

* [TensorFlow (v1.12) ](https://www.tensorflow.org/install/) 
* [Python3](https://docs.python-guide.org/starting/install3/linux/)
* pandas
* numpy
* tqdm
* cPickle
* argparse
* six

To install:

1. `git clone git@github.com:jaswindersingh2/SPOT-RNA.git`
2. `cd SPOT-RNA`
2. `wget 'https://www.dropbox.com/s/dsrcf460nbjqpxa/SPOT-RNA-models.tar.gz' && tar -xvzf SPOT-RNA-models.tar.gz && rm SPOT-RNA-models.tar.gz`
3. `pip install -r requirements.txt`

How to Use the SPOT-RNA Scripts
-----

```
python3 SPOT-RNA.py  --input_seq samples/4qln-1-A  --output_dir 'outputs/'
```

```
and check the outputs against the SPOT-RNA files in sample directory. To specify running on GPU please set the --gpu argument. Furthermore, SPOT-RNA also accepts multiple sequences in one fasta file. 
Training:
```

Datasets Used For Training, Validation, and Testing
-----

The following datasets were used for Initial Training:
* bpRNA [here](https://www.dropbox.com/s/w3kc4iro8ztbf3m/bpRNA_dataset.zip?dl=0).


The following datasets were used for Transfer Learning:
* PDB [here](https://www.dropbox.com/s/rlr8n9r5mt456cd/PDB_dataset.zip?dl=0).

References
-----
```
Hanson, J., Paliwal, K. Yang, Y., and Zhou, Y., 2018. RNA Secondary Structure Prediction using an Ensemble of Two-dimensional Recurrent and Residual Convolutional Neural Networks and Transfer Learning.
```
