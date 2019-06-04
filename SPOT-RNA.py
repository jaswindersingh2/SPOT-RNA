import tensorflow as tf
import numpy as np
import os
from tqdm import tqdm
import pickle as pkl
import argparse
from utils import create_tfr_files, prob_to_secondary_structure
import time
start = time.time()

parser = argparse.ArgumentParser()
parser.add_argument('--input_seq', help='name of input file')
parser.add_argument('--non_canonical',default=True, help='whether to include non-canoical pairs or not')
parser.add_argument('--output_dir',default='outputs', type=str, help='specify path to save results')
parser.add_argument('--gpu', default=-1, type=int, help='Which GPU to use')
args = parser.parse_args()

create_tfr_files(args.input_seq)

with open(args.input_seq) as file:
    input_data = [line.strip() for line in file.read().splitlines() if line.strip()]

count = int(len(input_data)/2)

ids = [input_data[2*i].replace(">", "") for i in range(count)]
sequences = {}
for i,I in enumerate(ids):
    sequences[I] = input_data[2*i+1].replace(" ", "")

os.environ["CUDA_VISIBLE_DEVICES"]= str(args.gpu)
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 
NUM_MODELS = 5

test_loc = ["input_tfr_files/test_data.tfrecords"]

outputs = {}
mask = {}
def sigmoid(x):
    return 1/(1+np.exp(-np.array(x, dtype=np.float128)))

for MODEL in range(NUM_MODELS):
    config = tf.ConfigProto()
    #config.gpu_options.allow_growth = True
    config.allow_soft_placement=True
    config.log_device_placement=False
    #session_conf = tf.ConfigProto(intra_op_parallelism_threads=1, inter_op_parallelism_threads=1)
    #sess = tf.Session(config=session_conf)
    print('\nPredicting for SPOT-RNA model '+str(MODEL))
    with tf.Session(config=config) as sess:
        saver = tf.train.import_meta_graph('dat'+'/model'+str(MODEL)+'.meta')
        saver.restore(sess,'dat'+'/model'+str(MODEL))
        graph = tf.get_default_graph()
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
            except:
                break
        pbar.close()
    tf.reset_default_graph()


RNA_ids = [i for i in list(outputs.keys())]
ensemble_outputs = {}

print('\nPost Processing and Saving Output')
for i in tqdm(RNA_ids):
    ensemble_outputs[i] = np.mean(outputs[i],0)
    prob_to_secondary_structure(ensemble_outputs[i], mask[i], sequences[i], i, non_canonical=args.non_canonical, save_result_path=args.output_dir)
with open('ensemble_outputs', 'wb') as f:
    pkl.dump(ensemble_outputs, f, pkl.HIGHEST_PROTOCOL)
print('\nFinished!')
end = time.time()
print('\nProcesssing Time {} seconds'.format(end - start))
