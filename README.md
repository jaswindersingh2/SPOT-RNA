;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;
;                           SPOT-RNA SECONDARY STRUCTURE PREDICTOR
;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

Hello,
Thank you for your interest in our predictor. SPOT-RNA requires 
several python packages, please verify that they are installed before use:

    tensorflow (v1.12.0)    ---> see https://www.tensorflow.org/install/
    pandas
    numpy
    tqdm
    cPickle
    argparse
    six

USAGE:
SPOT-RNA uses RNA sequence as input:


Model Availability:

The models used for SPOT-RNA are available from the Nihao cloud service [here](https://www.dropbox.com/s/g9ic16uz41zv9bp/SPOT-RNA-models.zip?dl=0)
The models used for SPOT-RNA are available from the Nihao cloud service [here](https://app.nihaocloud.com/f/882db8caf4be43ddaa04/?dl=1)


To run for a sample RNA, please  run:

python3 SPOT-RNA.py  --input_seq samples/4qln-1-A  --output_dir 'outputs/'
  
and check the outputs against the SPOT-RNA files in sample directory. To specify running on GPU please set the --gpu argument. Furthermore, SPOT-RNA also accepts multiple sequences in one fasta file. 

If you use this predictor in your research, please cite:

Hanson, J., Paliwal, K. Yang, Y., and Zhou, Y., 2018. Accurately Predicting Protein Contact Maps with Residual Two-Dimensional Bidirectional Long Short-Term Memory and Convolutional Neural Networks.


Thanks again!
Jaswinder Singh
