SPOT-RNA: RNA Secondary Structure Prediction using an Ensemble of Two-dimensional Recurrent and Residual Convolutional Neural Networks and Transfer Learning.
====

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
* [bpRNA](https://www.dropbox.com/s/w3kc4iro8ztbf3m/bpRNA_dataset.zip).


The following datasets were used for Transfer Learning:
* [PDB](https://www.dropbox.com/s/rlr8n9r5mt456cd/PDB_dataset.zip).

References
-----
```
Singh, J., Hanson, J., Paliwal, K. Yang, Y., and Zhou, Y., 2019. RNA Secondary Structure Prediction using an Ensemble of Two-dimensional Recurrent and Residual Convolutional Neural Networks and Transfer Learning.
```
