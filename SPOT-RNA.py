import tensorflow as tf
import numpy as np
import os
from tqdm import tqdm
import argparse
from utils.utils import create_tfr_files, prob_to_secondary_structure
from utils.FastaMLtoSL import FastaMLtoSL
import time
start = time.time()
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser()
parser.add_argument('--inputs', default='sample_inputs/2zzm-B.fasta', type=str, help='Path to input file in fasta format, accept multiple sequences as well in fasta format; default = ''sample_inputs/single_seq.fasta''\n', metavar='')
parser.add_argument('--outputs',default='outputs/', type=str, help='Path to output files; SPOT-RNA outputs at least three files .ct, .bpseq, and .prob files; default = ''outputs/\n', metavar='')
parser.add_argument('--gpu', default=-1, type=int, help='To run on GPU, specifiy GPU number. If only one GPU in computer specifiy 0; default = -1 (no GPU)\n', metavar='')
parser.add_argument('--plots',default=False, type=bool, help='Set this to "True" to get the 2D plots of predicted secondary structure by SPOT-RNA; default = False\n', metavar='')
parser.add_argument('--motifs',default=False, type=bool, help='Set this to "True" to get the motifs of predicted secondary structure by SPOT-RNA; default = False\n', metavar='')
parser.add_argument('--cpu',default=16, type=int, help='Specify number of cpu threads that SPOT-RNA can use; default = 16\n', metavar='')
#parser.add_argument('--NC',default=True, type=bool, help='Set this to "False" to predict only canonical pairs; default = True\n', metavar='')
args = parser.parse_args()

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)

FastaMLtoSL(args.inputs)

base_path = os.path.dirname(os.path.realpath(__file__))
input_file = os.path.basename(args.inputs)

create_tfr_files(args.inputs, base_path, input_file)

with open(args.inputs) as file:
    input_data = [line.strip() for line in file.read().splitlines() if line.strip()]

count = int(len(input_data)/2)

ids = [input_data[2*i].replace(">", "") for i in range(count)]
sequences = {}
for i,I in enumerate(ids):
    sequences[I] = input_data[2*i+1].replace(" ", "").upper().replace("T", "U")

os.environ["CUDA_VISIBLE_DEVICES"]= str(args.gpu)
#os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 
NUM_MODELS = 5

test_loc = [os.path.join(base_path, 'input_tfr_files', input_file+'.tfrecords')]

outputs = {}
mask = {}
def sigmoid(x):
    return 1/(1+np.exp(-np.array(x, dtype=np.float128)))

for MODEL in range(NUM_MODELS):

    if args.gpu==-1:
            config = tf.compat.v1.ConfigProto(intra_op_parallelism_threads=args.cpu, inter_op_parallelism_threads=args.cpu)
    else:
	    config = tf.compat.v1.ConfigProto()
	    config.allow_soft_placement=True
	    config.log_device_placement=False
	    config.gpu_options.allow_growth = True

    print('\nPredicting for SPOT-RNA model '+str(MODEL))
    with tf.compat.v1.Session(config=config) as sess:
        saver = tf.compat.v1.train.import_meta_graph(os.path.join(base_path, 'SPOT-RNA-models', 'model' + str(MODEL) + '.meta'))
        saver.restore(sess,os.path.join(base_path, 'SPOT-RNA-models', 'model' + str(MODEL)))
        graph = tf.compat.v1.get_default_graph()
        init_test =  graph.get_operation_by_name('make_initializer_2')
        tmp_out = graph.get_tensor_by_name('output_FC/fully_connected/BiasAdd:0')
        name_tensor = graph.get_tensor_by_name('tensors_2/component_0:0')
        RNA_name = graph.get_tensor_by_name('IteratorGetNext:0')
        label_mask = graph.get_tensor_by_name('IteratorGetNext:4')
        sess.run(init_test,feed_dict={name_tensor:test_loc})
        
        pbar = tqdm(total = count)
        while True:
            try:        
                out = sess.run([tmp_out,RNA_name,label_mask],feed_dict={'dropout:0':1})
                out[1] = out[1].decode()
                mask[out[1]] = out[2]
                
                if MODEL == 0:
                    outputs[out[1]] = [sigmoid(out[0])]
                else:
                    outputs[out[1]].append(sigmoid(out[0]))
                #print('RNA name: %s'%(out[1]))
                pbar.update(1)
            except tf.errors.OutOfRangeError:
                break
        pbar.close()
    tf.compat.v1.reset_default_graph()


RNA_ids = [i for i in list(outputs.keys())]
ensemble_outputs = {}

print('\nPost Processing and Saving Output')
for i in RNA_ids:
    ensemble_outputs[i] = np.mean(outputs[i],0)
    prob_to_secondary_structure(ensemble_outputs[i], mask[i], sequences[i], i, args, base_path)

print('\nFinished!')
end = time.time()
print('\nProcesssing Time {} seconds'.format(end - start))
