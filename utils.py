import numpy as np
import os, six, sys
import tensorflow as tf
import random
from tqdm import tqdm

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

def create_tfr_files(all_seq):

    print('\nPreparing tfr records file for SPOT-RNA:')
    path_tfrecords = os.path.join('input_tfr_files',"test_data"+".tfrecords")
    with open(all_seq) as file:
        input_data = [line.strip() for line in file.read().splitlines() if line.strip()]

    count = int(len(input_data)/2)

    ids = [input_data[2*i][1:].strip() for i in range(count)]
    
    with tf.python_io.TFRecordWriter(path_tfrecords) as writer:
        for i in tqdm(range(len(ids))):
            name     = input_data[2*i].replace(">", "") 
            sequence = input_data[2*i+1].replace(" ", "")
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

# ----------------------- remove triplet pairs--------------------------------#
def find_triplet_pairs(pred_pairs, y_pred):

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

    delete_pair = []
    for i in temp3:
        if len(i) == 2:
            if y_pred[i[0][0], i[0][1]] > y_pred[i[1][0], i[1][1]]:
                delete_pair.append(i[1])
            else:
                delete_pair.append(i[0])
        elif len(i) == 3:
            if y_pred[i[0][0], i[0][1]] > y_pred[i[1][0], i[1][1]] and y_pred[i[0][0], i[0][1]] > y_pred[
                i[2][0], i[2][1]]:
                delete_pair.append(i[1])
                delete_pair.append(i[2])
            elif y_pred[i[1][0], i[1][1]] > y_pred[i[0][0], i[0][1]] and y_pred[i[1][0], i[1][1]] > y_pred[
                i[2][0], i[2][1]]:
                delete_pair.append(i[0])
                delete_pair.append(i[2])
            else:
                delete_pair.append(i[0])
                delete_pair.append(i[1])
        elif len(i) == 4:
            if y_pred[i[0][0], i[0][1]] > y_pred[i[1][0], i[1][1]] and y_pred[i[0][0], i[0][1]] > y_pred[i[2][0], i[2][1]] and y_pred[i[0][0], i[0][1]] > y_pred[i[3][0], i[3][1]]:
                delete_pair.append(i[1])
                delete_pair.append(i[2])
                delete_pair.append(i[3])
            elif y_pred[i[1][0], i[1][1]] > y_pred[i[0][0], i[0][1]] and y_pred[i[1][0], i[1][1]] > y_pred[i[2][0], i[2][1]] and y_pred[i[1][0], i[1][1]] > y_pred[i[3][0], i[3][1]]:
                delete_pair.append(i[0])
                delete_pair.append(i[2])
                delete_pair.append(i[3])              
            elif y_pred[i[2][0], i[2][1]] > y_pred[i[0][0], i[0][1]] and y_pred[i[2][0], i[2][1]] > y_pred[i[1][0], i[1][1]] and y_pred[i[2][0], i[2][1]] > y_pred[i[3][0], i[3][1]]:
                delete_pair.append(i[0])
                delete_pair.append(i[1])
                delete_pair.append(i[3])
        elif len(i) == 5:
            if y_pred[i[0][0], i[0][1]] > y_pred[i[1][0], i[1][1]] and y_pred[i[0][0], i[0][1]] > y_pred[i[2][0], i[2][1]] and y_pred[i[0][0], i[0][1]] > y_pred[i[3][0], i[3][1]] and y_pred[i[0][0], i[0][1]] > y_pred[i[4][0], i[4][1]]:
                delete_pair.append(i[1])
                delete_pair.append(i[2])
                delete_pair.append(i[3])
                delete_pair.append(i[4])
            elif y_pred[i[1][0], i[1][1]] > y_pred[i[0][0], i[0][1]] and y_pred[i[1][0], i[1][1]] > y_pred[i[2][0], i[2][1]] and y_pred[i[1][0], i[1][1]] > y_pred[i[3][0], i[3][1]] and y_pred[i[1][0], i[1][1]] > y_pred[i[4][0], i[4][1]]:
                delete_pair.append(i[1])
                delete_pair.append(i[2])
                delete_pair.append(i[3])
                delete_pair.append(i[4])              
            elif y_pred[i[2][0], i[2][1]] > y_pred[i[0][0], i[0][1]] and y_pred[i[2][0], i[2][1]] > y_pred[i[1][0], i[1][1]] and y_pred[i[2][0], i[2][1]] > y_pred[i[3][0], i[3][1]] and y_pred[i[2][0], i[2][1]] > y_pred[i[4][0], i[4][1]]:
                delete_pair.append(i[0])
                delete_pair.append(i[1])
                delete_pair.append(i[3])
                delete_pair.append(i[4])
            elif y_pred[i[3][0], i[3][1]] > y_pred[i[0][0], i[0][1]] and y_pred[i[3][0], i[3][1]] > y_pred[i[1][0], i[1][1]] and y_pred[i[3][0], i[3][1]] > y_pred[i[2][0], i[2][1]] and y_pred[i[3][0], i[3][1]] > y_pred[i[4][0], i[4][1]]:
                delete_pair.append(i[0])
                delete_pair.append(i[1])
                delete_pair.append(i[2])
                delete_pair.append(i[4])
            else:
                delete_pair.append(i[0])
                delete_pair.append(i[1])
                delete_pair.append(i[2])
                delete_pair.append(i[3])
        elif len(i) == 6:
            if y_pred[i[0][0], i[0][1]] > y_pred[i[1][0], i[1][1]] and y_pred[i[0][0], i[0][1]] > y_pred[i[2][0], i[2][1]] and y_pred[i[0][0], i[0][1]] > y_pred[i[3][0], i[3][1]] and y_pred[i[0][0], i[0][1]] > y_pred[i[4][0], i[4][1]] and y_pred[i[0][0], i[0][1]] > y_pred[i[5][0], i[5][1]]:
                delete_pair.append(i[1])
                delete_pair.append(i[2])
                delete_pair.append(i[3])
                delete_pair.append(i[4])
                delete_pair.append(i[5])
            elif y_pred[i[1][0], i[1][1]] > y_pred[i[0][0], i[0][1]] and y_pred[i[1][0], i[1][1]] > y_pred[i[2][0], i[2][1]] and y_pred[i[1][0], i[1][1]] > y_pred[i[3][0], i[3][1]] and y_pred[i[1][0], i[1][1]] > y_pred[i[4][0], i[4][1]] and y_pred[i[1][0], i[1][1]] > y_pred[i[5][0], i[5][1]]:
                delete_pair.append(i[1])
                delete_pair.append(i[2])
                delete_pair.append(i[3])
                delete_pair.append(i[4])              
                delete_pair.append(i[5])              
            elif y_pred[i[2][0], i[2][1]] > y_pred[i[0][0], i[0][1]] and y_pred[i[2][0], i[2][1]] > y_pred[i[1][0], i[1][1]] and y_pred[i[2][0], i[2][1]] > y_pred[i[3][0], i[3][1]] and y_pred[i[2][0], i[2][1]] > y_pred[i[4][0], i[4][1]] and y_pred[i[2][0], i[2][1]] > y_pred[i[5][0], i[5][1]]:
                delete_pair.append(i[0])
                delete_pair.append(i[1])
                delete_pair.append(i[3])
                delete_pair.append(i[4])
                delete_pair.append(i[5])
            elif y_pred[i[3][0], i[3][1]] > y_pred[i[0][0], i[0][1]] and y_pred[i[3][0], i[3][1]] > y_pred[i[1][0], i[1][1]] and y_pred[i[3][0], i[3][1]] > y_pred[i[2][0], i[2][1]] and y_pred[i[3][0], i[3][1]] > y_pred[i[4][0], i[4][1]] and y_pred[i[3][0], i[3][1]] > y_pred[i[5][0], i[5][1]]:
                delete_pair.append(i[0])
                delete_pair.append(i[1])
                delete_pair.append(i[2])
                delete_pair.append(i[4])
                delete_pair.append(i[5])
            elif y_pred[i[4][0], i[4][1]] > y_pred[i[0][0], i[0][1]] and y_pred[i[4][0], i[4][1]] > y_pred[i[1][0], i[1][1]] and y_pred[i[4][0], i[4][1]] > y_pred[i[2][0], i[2][1]] and y_pred[i[4][0], i[4][1]] > y_pred[i[3][0], i[3][1]] and y_pred[i[4][0], i[4][1]] > y_pred[i[5][0], i[5][1]]:
                delete_pair.append(i[0])
                delete_pair.append(i[1])
                delete_pair.append(i[2])
                delete_pair.append(i[3])
                delete_pair.append(i[5])
            else:
                delete_pair.append(i[0])
                delete_pair.append(i[1])
                delete_pair.append(i[2])
                delete_pair.append(i[3])
                delete_pair.append(i[4])
        elif len(i) > 6:
            print('one residue is in contact with more than 6 residues', len(i))

    return delete_pair

def output_mask(seq, non_canonical=True):
    if non_canonical:
        include_pairs = ['AU', 'UA', 'GC', 'CG', 'GU', 'UG', 'CC', 'GG', 'AG', 'CA', 'AC', 'UU', 'AA', 'CU', 'GA', 'UC']
    else:
        include_pairs = ['AU', 'UA', 'GC', 'CG', 'GU', 'UG']
    mask = np.zeros((len(seq), len(seq)))
    for i, I in enumerate(seq):
        for j, J in enumerate(seq):
            if str(I) + str(J) in include_pairs:
                mask[i, j] = 1
    return  mask

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
    np.savetxt(os.path.join(save_result_path, str(id))+'.spotrna', (temp), delimiter='\t\t\t\t', fmt="%s", header=str(len(seq)) + '\t\t\t\t' + str(id) + '\t\t\t\t' + 'SPOT-RNA output\n' , comments='')

    return

def prob_to_secondary_structure(ensemble_outputs, label_mask, seq, name, non_canonical, save_result_path):
    #save_result_path = 'outputs'
    Threshold = 0.335
    test_output = ensemble_outputs
    mask = output_mask(seq, non_canonical)
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
    delete_pairs = find_triplet_pairs(pred_pairs, y_pred)
    bad_pairs = hair_pin_assumption(pred_pairs)
    pred_pairs = [i for i in pred_pairs if i not in delete_pairs]   # remove triplets
    #pred_pairs = [i for i in pred_pairs if i not in bad_pairs]     # remove not follow hair pin assumption
    ct_file_output(pred_pairs, seq, name, save_result_path)
    np.savetxt('outputs/'+ name +'.prob', y_pred, delimiter='\t')

    
    
