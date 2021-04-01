import numpy as np
import os, six, sys, subprocess, time
import tensorflow as tf
import random
from tqdm import tqdm
from subprocess import PIPE
import pandas as pd

# ------------- one hot encoding of RNA sequences -----------------#
def one_hot(seq):
    RNN_seq = seq
    BASES = 'AUCG'
    bases = np.array([base for base in BASES])
    feat = np.concatenate(
        [[(bases == base.upper()).astype(int)] if str(base).upper() in BASES else np.array([[-1] * len(BASES)]) for base
         in RNN_seq])

    return feat


def z_mask(seq_len):
    mask = np.ones((seq_len, seq_len))
    return np.triu(mask, 2)


def l_mask(inp, seq_len):
    temp = []
    mask = np.ones((seq_len, seq_len))
    for k, K in enumerate(inp):
        if np.any(K == -1) == True:
            temp.append(k)
    mask[temp, :] = 0
    mask[:, temp] = 0
    return np.triu(mask, 2)


def get_data(seq):
    seq_len = len(seq)
    one_hot_feat = one_hot(seq)
    #print(one_hot_feat[-1])
    zero_mask = z_mask(seq_len)[None, :, :, None]
    label_mask = l_mask(one_hot_feat, seq_len)
    temp = one_hot_feat[None, :, :]
    temp = np.tile(temp, (temp.shape[1], 1, 1))
    feature = np.concatenate([temp, np.transpose(temp, [1, 0, 2])], 2)
    #out = true_output

    return seq_len, [i for i in (feature.astype(float)).flatten()], [i for i in zero_mask.flatten()], [i for i in label_mask.flatten()], [i for i in label_mask.flatten()]

def _int64_feature(value):
    if not isinstance(value, list) and not isinstance(value, np.ndarray):
        value = [value]

    return tf.train.Feature(int64_list=tf.train.Int64List(value=value))


def _float_feature(value):
    if not isinstance(value, list) and not isinstance(value, np.ndarray):
        value = [value]

    return tf.train.Feature(float_list=tf.train.FloatList(value=value))


def _bytes_feature(value):
    """Wrapper for inserting bytes features into Example proto."""
    if isinstance(value, six.string_types):
        value = six.binary_type(value, encoding='utf-8')
    return tf.train.Feature(bytes_list=tf.train.BytesList(value=[value]))

def create_tfr_files(all_seq, base_path, input_file):

    print('\nPreparing tfr records file for SPOT-RNA:')
    path_tfrecords = os.path.join(base_path, 'input_tfr_files', input_file+'.tfrecords')
    with open(all_seq) as file:
        input_data = [line.strip() for line in file.read().splitlines() if line.strip()]

    count = int(len(input_data)/2)

    ids = [input_data[2*i][1:].strip() for i in range(count)]
    
    with tf.io.TFRecordWriter(path_tfrecords) as writer:
        for i in tqdm(range(len(ids))):
            name     = input_data[2*i].replace(">", "") 
            sequence = input_data[2*i+1].replace(" ", "").upper().replace("T", "U")
            #print(sequence[-1])
            
            #print(len(sequence), name)                
            seq_len, feature, zero_mask, label_mask, true_label = get_data(sequence)

            example = tf.train.Example(features=tf.train.Features(feature={'rna_name': _bytes_feature(name),
                                                                           'seq_len': _int64_feature(seq_len),
                                                                           'feature': _float_feature(feature),
                                                                           'zero_mask': _float_feature(zero_mask),
                                                                           'label_mask': _float_feature(label_mask),
                                                                           'true_label': _float_feature(true_label)}))

            writer.write(example.SerializeToString())

    writer.close()

# ----------------------- hair pin loop assumption i - j < 2 --------------------------------#
def hair_pin_assumption(pred_pairs):
    pred_pairs_all = [i[:2] for i in pred_pairs]
    bad_pairs = []
    for i in pred_pairs_all:
        if abs(i[0] - i[1]) < 3:
            bad_pairs.append(i)
    return bad_pairs

def flatten(x):
    result = []
    for el in x:
        if hasattr(el, "__iter__") and not isinstance(el, str):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result

def type_pairs(pairs, sequence):
    sequence = [i.upper() for i in sequence]
    # seq_pairs = [[sequence[i[0]],sequence[i[1]]] for i in pairs]

    AU_pair = []
    GC_pair = []
    GU_pair = []
    other_pairs = []
    for i in pairs:
        if [sequence[i[0]],sequence[i[1]]] in [["A","U"], ["U","A"]]:
            AU_pair.append(i)
        elif [sequence[i[0]],sequence[i[1]]] in [["G","C"], ["C","G"]]:
            GC_pair.append(i)
        elif [sequence[i[0]],sequence[i[1]]] in [["G","U"], ["U","G"]]:
            GU_pair.append(i)
        else:
            other_pairs.append(i)
    watson_pairs_t = AU_pair + GC_pair
    wobble_pairs_t = GU_pair
    other_pairs_t = other_pairs
        # print(watson_pairs_t, wobble_pairs_t, other_pairs_t)
    return watson_pairs_t, wobble_pairs_t, other_pairs_t

# ----------------------- find multiplets pairs--------------------------------#
def multiplets_pairs(pred_pairs):

    pred_pair = [i[:2] for i in pred_pairs]
    temp_list = flatten(pred_pair)
    temp_list.sort()
    new_list = sorted(set(temp_list))
    dup_list = []
    for i in range(len(new_list)):
        if (temp_list.count(new_list[i]) > 1):
            dup_list.append(new_list[i])

    dub_pairs = []
    for e in pred_pair:
        if e[0] in dup_list:
            dub_pairs.append(e)
        elif e[1] in dup_list:
            dub_pairs.append(e)

    temp3 = []
    for i in dup_list:
        temp4 = []
        for k in dub_pairs:
            if i in k:
                temp4.append(k)
        temp3.append(temp4)
        
    return temp3

def multiplets_free_bp(pred_pairs, y_pred):
    L = len(pred_pairs)
    multiplets_bp = multiplets_pairs(pred_pairs)
    save_multiplets = []
    while len(multiplets_bp) > 0:
        remove_pairs = []
        for i in multiplets_bp:
            save_prob = []
            for j in i:
                save_prob.append(y_pred[j[0], j[1]])
            remove_pairs.append(i[save_prob.index(min(save_prob))])
            save_multiplets.append(i[save_prob.index(min(save_prob))])
        pred_pairs = [k for k in pred_pairs if k not in remove_pairs]
        multiplets_bp = multiplets_pairs(pred_pairs)
    save_multiplets = [list(x) for x in set(tuple(x) for x in save_multiplets)]
    assert L == len(pred_pairs)+len(save_multiplets)
    #print(L, len(pred_pairs), save_multiplets)
    return pred_pairs, save_multiplets
        
def output_mask(seq, NC=True):
    if NC:
        include_pairs = ['AU', 'UA', 'GC', 'CG', 'GU', 'UG', 'CC', 'GG', 'AG', 'CA', 'AC', 'UU', 'AA', 'CU', 'GA', 'UC']
    else:
        include_pairs = ['AU', 'UA', 'GC', 'CG', 'GU', 'UG']
    mask = np.zeros((len(seq), len(seq)))
    for i, I in enumerate(seq):
        for j, J in enumerate(seq):
            if str(I) + str(J) in include_pairs:
                mask[i, j] = 1
    return mask

def ct_file_output(pairs, seq, id, save_result_path):

    col1 = np.arange(1, len(seq) + 1, 1)
    col2 = np.array([i for i in seq])
    col3 = np.arange(0, len(seq), 1)
    col4 = np.append(np.delete(col1, 0), [0])
    col5 = np.zeros(len(seq), dtype=int)

    for i, I in enumerate(pairs):
        col5[I[0]] = int(I[1]) + 1
        col5[I[1]] = int(I[0]) + 1
    col6 = np.arange(1, len(seq) + 1, 1)
    temp = np.vstack((np.char.mod('%d', col1), col2, np.char.mod('%d', col3), np.char.mod('%d', col4),
                      np.char.mod('%d', col5), np.char.mod('%d', col6))).T
    #os.chdir(save_result_path)
    #print(os.path.join(save_result_path, str(id[0:-1]))+'.spotrna')
    np.savetxt(os.path.join(save_result_path, str(id))+'.ct', (temp), delimiter='\t\t', fmt="%s", header=str(len(seq)) + '\t\t' + str(id) + '\t\t' + 'SPOT-RNA output\n' , comments='')

    return

def bpseq_file_output(pairs, seq, id, save_result_path):

    col1 = np.arange(1, len(seq) + 1, 1)
    col2 = np.array([i for i in seq])
    #col3 = np.arange(0, len(seq), 1)
    #col4 = np.append(np.delete(col1, 0), [0])
    col5 = np.zeros(len(seq), dtype=int)

    for i, I in enumerate(pairs):
        col5[I[0]] = int(I[1]) + 1
        col5[I[1]] = int(I[0]) + 1
    #col6 = np.arange(1, len(seq) + 1, 1)
    temp = np.vstack((np.char.mod('%d', col1), col2, np.char.mod('%d', col5))).T
    #os.chdir(save_result_path)
    #print(os.path.join(save_result_path, str(id[0:-1]))+'.spotrna')
    np.savetxt(os.path.join(save_result_path, str(id))+'.bpseq', (temp), delimiter=' ', fmt="%s", header='#' + str(id) , comments='')

    return

def lone_pair(pairs):
    lone_pairs = []
    pairs.sort()
    for i, I in enumerate(pairs):
        if ([I[0] - 1, I[1] + 1] not in pairs) and ([I[0] + 1, I[1] - 1] not in pairs):
            lone_pairs.append(I)

    return lone_pairs

def prob_to_secondary_structure(ensemble_outputs, label_mask, seq, name, args, base_path):
    #save_result_path = 'outputs'
    Threshold = 0.335
    test_output = ensemble_outputs
    mask = output_mask(seq)
    inds = np.where(label_mask == 1)
    y_pred = np.zeros(label_mask.shape)
    for i in range(test_output.shape[0]):
        y_pred[inds[0][i], inds[1][i]] = test_output[i]
    y_pred = np.multiply(y_pred, mask)

    tri_inds = np.triu_indices(y_pred.shape[0], k=1)

    out_pred = y_pred[tri_inds]
    outputs = out_pred[:, None]
    seq_pairs = [[tri_inds[0][j], tri_inds[1][j], ''.join([seq[tri_inds[0][j]], seq[tri_inds[1][j]]])] for j in
                 range(tri_inds[0].shape[0])]

    outputs_T = np.greater_equal(outputs, Threshold)
    pred_pairs = [i for I, i in enumerate(seq_pairs) if outputs_T[I]]
    pred_pairs = [i[:2] for i in pred_pairs]
    pred_pairs, save_multiplets = multiplets_free_bp(pred_pairs, y_pred)
    
    watson_pairs, wobble_pairs, noncanonical_pairs = type_pairs(pred_pairs, seq)
    lone_bp = lone_pair(pred_pairs)

    tertiary_bp = save_multiplets + noncanonical_pairs + lone_bp
    tertiary_bp = [list(x) for x in set(tuple(x) for x in tertiary_bp)]

    str_tertiary = []
    for i,I in enumerate(tertiary_bp):
        if i==0: 
            str_tertiary += ('(' + str(I[0]+1) + ',' + str(I[1]+1) + '):color=""#FFFF00""')
        else:
            str_tertiary += (';(' + str(I[0]+1) + ',' + str(I[1]+1) + '):color=""#FFFF00""')   
        
    tertiary_bp = ''.join(str_tertiary) 

    if args.outputs=='outputs/':
        output_path = os.path.join(base_path, args.outputs)
    else:
        output_path = args.outputs

    ct_file_output(pred_pairs, seq, name, output_path)
    bpseq_file_output(pred_pairs, seq, name, output_path)
    np.savetxt(output_path + '/'+ name +'.prob', y_pred, delimiter='\t')
    
    if args.plots:
        try:
            subprocess.Popen(["java", "-cp", base_path + "/utils/VARNAv3-93.jar", "fr.orsay.lri.varna.applications.VARNAcmd", '-i', output_path + name + '.ct', '-o', output_path + name + '_radiate.png', '-algorithm', 'radiate', '-resolution', '8.0', '-bpStyle', 'lw', '-auxBPs', tertiary_bp], stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
            subprocess.Popen(["java", "-cp", base_path + "/utils/VARNAv3-93.jar", "fr.orsay.lri.varna.applications.VARNAcmd", '-i', output_path + name + '.ct', '-o', output_path + name + '_line.png', '-algorithm', 'line', '-resolution', '8.0', '-bpStyle', 'lw', '-auxBPs', tertiary_bp], stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
        except:
            print('\nUnable to generate 2D plots;\nplease refer to "http://varna.lri.fr/" for system requirments to use VARNA')	

    if args.motifs:
        try:
            os.chdir(output_path)
            p = subprocess.Popen(['perl', base_path + '/utils/bpRNA-master/bpRNA.pl', name + '.bpseq'])
            time.sleep(0.1)
            #print(os.path.exists(output_path + '/' + name + '.st'))
            with open(output_path + '/' + name + '.st') as f:
                df = pd.read_csv(f, comment='#', sep=";", header=None)
            np.savetxt(output_path + '/'+ name +'.dbn', np.array([df[0][0], df[0][1]]), fmt="%s", header='>' + name, comments='')
        except:
            print('\nUnable to run bpRNA script;\nplease refer to "https://github.com/hendrixlab/bpRNA/" for system requirments to use bpRNA')
        os.chdir(base_path)
    return
