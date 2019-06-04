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

To run for a sample RNA, please  run:

python3 SPOT-RNA.py  --input_seq samples/4qln-1-A  --output_dir 'outputs/'
  
and check the outputs against the SPOT-RNA files in sample directory. To specify running on GPU please set the --gpu argument. Furthermore, SPOT-RNA also accepts multiple sequences in one fasta file. 

If you use this predictor in your research, please cite:

Hanson, J., Paliwal, K. Yang, Y., and Zhou, Y., 2018. Accurately Predicting Protein Contact Maps with Residual Two-Dimensional Bidirectional Long Short-Term Memory and Convolutional Neural Networks.


Thanks again!
Jaswinder Singh
