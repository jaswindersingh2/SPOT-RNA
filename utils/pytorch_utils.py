import numpy as np
import os
# import six # Not used by current functions
import sys # Used by subprocess indirectly? Or by FastaMLtoSL's original exit(1)
import subprocess # Used by prob_to_secondary_structure for VARNA/bpRNA
import time # Used by prob_to_secondary_structure for bpRNA
import torch # PyTorch library, core of this version
# import random # Not explicitly used by current functions, can be removed if not needed by future ones
# from tqdm import tqdm # Not used by current functions

from subprocess import PIPE # Used by prob_to_secondary_structure for VARNA/bpRNA
import pandas as pd # Used by prob_to_secondary_structure for bpRNA
import re # Used by FastaMLtoSL

# ------------- RNA sequence processing and data preparation functions -----------------#

def one_hot(seq: str) -> np.ndarray:
    """
    Converts an RNA sequence string to its one-hot encoded NumPy array.
    'A' -> [1,0,0,0], 'U' -> [0,1,0,0], 'C' -> [0,0,1,0], 'G' -> [0,0,0,1].
    Unknown characters are mapped to [-1,-1,-1,-1].

    Args:
        seq (str): The RNA sequence string.

    Returns:
        np.ndarray: A 2D NumPy array of shape (sequence_length, 4) representing the one-hot encoded sequence.
    """
    RNN_seq = seq
    BASES = 'AUCG'
    bases = np.array([base for base in BASES])
    # For each base in the input sequence, compare it against BASES.
    # If in BASES, create a one-hot vector. Otherwise, create a vector of -1s.
    feat = np.concatenate(
        [[(bases == base.upper()).astype(int)] if str(base).upper() in BASES else np.array([[-1] * len(BASES)]) for base
         in RNN_seq])
    return feat

def z_mask(seq_len: int) -> np.ndarray:
    """
    Creates a zero mask (upper triangle) for a square matrix of given sequence length.
    This mask is often used to ignore certain parts of pair probability maps,
    specifically, it keeps the upper triangle starting from k=2.
    (i.e., masks diagonal and first off-diagonal).

    Args:
        seq_len (int): The length of the RNA sequence.

    Returns:
        np.ndarray: A 2D NumPy array of shape (seq_len, seq_len) with ones in the
                    upper triangle (k=2) and zeros elsewhere.
    """
    mask = np.ones((seq_len, seq_len))
    return np.triu(mask, 2) # np.triu returns the upper triangle of an array. k=2 means exclude main diagonal and first off-diagonal.

def l_mask(inp: np.ndarray, seq_len: int) -> np.ndarray:
    """
    Creates a label mask based on the input one-hot encoded sequence.
    If any row in the one-hot encoding contains -1 (unknown base),
    the corresponding rows and columns in the mask are set to 0.
    The mask is also an upper triangle (k=2).

    Args:
        inp (np.ndarray): The one-hot encoded sequence (output of one_hot function).
        seq_len (int): The length of the RNA sequence.

    Returns:
        np.ndarray: A 2D NumPy array (label_mask) of shape (seq_len, seq_len).
    """
    temp = [] # Stores indices of rows with unknown bases
    mask = np.ones((seq_len, seq_len))
    for k, K_row in enumerate(inp): # K_row is a one-hot vector for a base
        if np.any(K_row == -1): # Check if the base was unknown
            temp.append(k)
    # Mask out rows and columns corresponding to unknown bases
    mask[temp, :] = 0
    mask[:, temp] = 0
    return np.triu(mask, 2) # Keep only the upper triangle (k=2)

def get_data(seq: str) -> tuple[int, torch.Tensor, torch.Tensor, torch.Tensor]:
    """
    Preprocesses an RNA sequence to generate input features and masks for a PyTorch model.

    Args:
        seq (str): The RNA sequence string.

    Returns:
        tuple:
            - seq_len (int): Length of the sequence.
            - features_tensor (torch.Tensor): (seq_len, seq_len, 8) - Pairwise features.
            - zero_mask_tensor (torch.Tensor): (1, seq_len, seq_len) - Mask for zeroing out invalid regions (PyTorch conv format).
            - label_mask_tensor (torch.Tensor): (seq_len, seq_len) - Mask based on unknown bases and upper triangle.
    """
    seq_len = len(seq)
    one_hot_feat_np = one_hot(seq) # NumPy array: (seq_len, 4)

    # Generate masks as NumPy arrays first
    zero_mask_np = z_mask(seq_len) # (seq_len, seq_len)
    label_mask_np = l_mask(one_hot_feat_np, seq_len) # (seq_len, seq_len)
    
    # Construct pairwise features (similar to original TensorFlow version's feature engineering)
    # 1. Expand one_hot_feat_np to (1, seq_len, 4)
    one_hot_feat_expanded_np = one_hot_feat_np[np.newaxis, :, :]
    # 2. Tile it to (seq_len, seq_len, 4) - effectively creating (L, L, C) where C is one-hot dim
    #    Each row i of this tiled array contains seq_len copies of one_hot_feat_np[i]
    #    This was `np.tile(temp, (temp.shape[1], 1, 1))` in original, where temp was (1, seq_len, 4)
    #    Corrected tiling to match original intent for pairwise features:
    #    We need one_hot_feat at position (i,j,0:4) and one_hot_feat at (i,j,4:8)
    #    So, tile one_hot_feat to (L,L,4) such that feature_map[i,j,:] is one_hot_feat[i]
    #    and another (L,L,4) such that feature_map_transposed[i,j,:] is one_hot_feat[j]
    
    # Create matrix1: matrix1[i, j, k] = one_hot_feat_np[i, k]
    matrix1 = np.tile(one_hot_feat_np[:, np.newaxis, :], (1, seq_len, 1)) # Shape (L, L, 4)
    # Create matrix2: matrix2[i, j, k] = one_hot_feat_np[j, k]
    matrix2 = np.tile(one_hot_feat_np[np.newaxis, :, :], (seq_len, 1, 1)) # Shape (L, L, 4)
    
    # Concatenate to get (L, L, 8) features
    features_np = np.concatenate([matrix1, matrix2], axis=2)

    # Convert NumPy arrays to PyTorch tensors
    features_tensor = torch.from_numpy(features_np.astype(np.float32))
    # Add a channel dimension for zero_mask, suitable for 2D convolutions if used: (1, seq_len, seq_len)
    zero_mask_tensor = torch.from_numpy(zero_mask_np.astype(np.float32)).unsqueeze(0) 
    label_mask_tensor = torch.from_numpy(label_mask_np.astype(np.float32)) # (seq_len, seq_len)

    return seq_len, features_tensor, zero_mask_tensor, label_mask_tensor

def FastaMLtoSL(inputfile_path: str, outputfile_path: str = None) -> str:
    """
    Converts a multi-line FASTA file to a single-line FASTA file.
    Each sequence in the output will have its header on one line and its entire sequence on the next.

    Args:
        inputfile_path (str): Path to the input multi-line FASTA file.
        outputfile_path (str, optional): Path to save the single-line FASTA file.
                                         If None, overwrites the input file. Defaults to None.

    Returns:
        str: The path to the output single-line FASTA file.

    Raises:
        IOError: If file operations fail.
    """
    if outputfile_path is None:
        outputfile_path = inputfile_path # Overwrite input if no output path specified
    
    try:    
        with open(inputfile_path,"r") as newFile:
            content = newFile.read()
            # Split by '>' character, ensuring it's at the beginning of a line
            if not content.startswith('>'): # Handle files that don't start with a FASTA header
                sequences = [] 
            else:
                sequences = re.split("^>", content, flags=re.MULTILINE)
                if sequences and sequences[0] == '': # Remove empty string from split if file starts with '>'
                    del sequences[0]
            newFile.close()
    except IOError:
        print(f"Failed to open FASTA file for reading: {inputfile_path}")
        raise 
    except Exception as e:
        print(f"An error occurred reading {inputfile_path}: {e}")
        raise

    try:    
        with open(outputfile_path,"w") as newFasta:
            for fasta_entry in sequences:
                try:
                    # Split into header and sequence parts
                    header, sequence_data = fasta_entry.split("\n", 1) 
                except ValueError: # Handles entries that might be just a header or malformed
                    if fasta_entry.strip(): # If there's any content (likely a header without sequence)
                        # Write header and an empty line for sequence to maintain FASTA format
                        newFasta.write(">" + fasta_entry.strip() + "\n\n") 
                    continue # Skip this malformed entry
                
                header = ">" + header.strip() + "\n" # Ensure header starts with > and ends with newline
                # Remove all internal newlines from sequence_data and ensure it ends with a single newline
                sequence_data = sequence_data.replace("\n","").strip() + "\n" 
                newFasta.write(header + sequence_data)
            newFasta.close()
    except IOError:
        print(f"Failed to open FASTA file for writing: {outputfile_path}")
        raise
    except Exception as e:
        print(f"An error occurred writing to {outputfile_path}: {e}")
        raise
        
    return outputfile_path # Return path of (potentially new) single-line FASTA file

# ------------- RNA Secondary Structure Prediction Post-processing Utilities -----------------#
# These functions are largely based on the original utils.py from SPOT-RNA (TensorFlow version)
# and are adapted here for use with PyTorch outputs (by converting tensors to NumPy arrays).

def hair_pin_assumption(pred_pairs: list) -> list:
    """Identifies base pairs that violate the hairpin assumption (loop size < 3)."""
    pred_pairs_all = [i[:2] for i in pred_pairs] # Ensure we only have [i,j]
    bad_pairs = []
    for i in pred_pairs_all:
        # A hairpin loop must have at least 3 unpaired bases.
        # So, if j and i are paired, j must be > i+3. (abs(i-j) > 3 or >=4)
        # The original code used abs(i[0] - i[1]) < 3, which means loop of 0, 1, or 2.
        # This should be abs(i[0] - i[1]) <= 3 for standard hairpin definition (min loop 3).
        if abs(i[0] - i[1]) <= 3: # If distance is 0,1,2,3 it's a "bad" hairpin (too short)
            bad_pairs.append(i)
    return bad_pairs

def flatten(x: list) -> list:
    """Flattens a list of lists into a single list."""
    result = []
    for el in x:
        if hasattr(el, "__iter__") and not isinstance(el, str): # Check if element is iterable but not a string
            result.extend(flatten(el))
        else:
            result.append(el)
    return result

def type_pairs(pairs: list, sequence: str) -> tuple[list, list, list]:
    """Classifies base pairs into Watson-Crick, Wobble, and Other types."""
    sequence_list = [s.upper() for s in sequence]

    AU_pair, GC_pair, GU_pair, other_pairs = [], [], [], []
    for i in pairs: # Each 'i' is a pair [idx0, idx1]
        p_chars = [sequence_list[i[0]], sequence_list[i[1]]]
        if p_chars in [["A","U"], ["U","A"]]:
            AU_pair.append(i)
        elif p_chars in [["G","C"], ["C","G"]]:
            GC_pair.append(i)
        elif p_chars in [["G","U"], ["U","G"]]: # Wobble pairs
            GU_pair.append(i)
        else:
            other_pairs.append(i) # Non-canonical or other pairings
            
    watson_crick_pairs = AU_pair + GC_pair
    wobble_pairs = GU_pair
    # other_pairs remains as is
    return watson_crick_pairs, wobble_pairs, other_pairs

def multiplets_pairs(pred_pairs: list) -> list:
    """Identifies groups of pairs that form multiplets (one base involved in multiple pairs)."""
    # pred_pairs is a list of [i,j] coordinates.
    base_indices_flat = flatten([p[:2] for p in pred_pairs]) # Get all base indices involved in pairs
    base_indices_flat.sort()
    
    # Find indices that are duplicated (involved in more than one pair)
    unique_indices = sorted(list(set(base_indices_flat)))
    duplicated_indices = []
    for idx in unique_indices:
        if base_indices_flat.count(idx) > 1:
            duplicated_indices.append(idx)

    # Collect all pairs that involve any of these duplicated indices
    multiplet_candidate_pairs = []
    for pair in pred_pairs:
        if pair[0] in duplicated_indices or pair[1] in duplicated_indices:
            multiplet_candidate_pairs.append(pair)

    # Group these candidate pairs into actual multiplet groups
    # A multiplet group consists of pairs that share common indices directly or indirectly.
    grouped_multiplets = []
    visited_pairs = set()

    for pair in multiplet_candidate_pairs:
        if tuple(sorted(pair)) in visited_pairs:
            continue

        current_multiplet_group = set()
        queue = [pair]
        visited_pairs.add(tuple(sorted(pair)))
        
        head = 0
        while head < len(queue):
            current_pair = queue[head]
            head += 1
            current_multiplet_group.add(tuple(sorted(current_pair)))
            
            # Find other pairs in multiplet_candidate_pairs that share an index with current_pair
            for other_candidate_pair in multiplet_candidate_pairs:
                # Ensure it's not the same pair and not yet visited
                if tuple(sorted(other_candidate_pair)) not in visited_pairs:
                    # Check for shared indices
                    if (other_candidate_pair[0] == current_pair[0] or \
                        other_candidate_pair[0] == current_pair[1] or \
                        other_candidate_pair[1] == current_pair[0] or \
                        other_candidate_pair[1] == current_pair[1]):
                        visited_pairs.add(tuple(sorted(other_candidate_pair)))
                        queue.append(other_candidate_pair)
        
        if current_multiplet_group: # current_multiplet_group contains tuples
             # Convert tuples back to lists for consistent output format
            grouped_multiplets.append([list(p_tuple) for p_tuple in sorted(list(current_multiplet_group))])
            
    return grouped_multiplets


def multiplets_free_bp(pred_pairs: list, y_pred_np: np.ndarray) -> tuple[list, list]:
    """
    Resolves multiplets by iteratively removing the pair with the lowest probability from each multiplet group.
    
    Args:
        pred_pairs (list): List of predicted pairs, e.g., [[0,10], [0,11], [5,15]].
        y_pred_np (np.ndarray): (L,L) NumPy array of prediction probabilities.

    Returns:
        tuple:
            - current_pred_pairs (list): Pairs after resolving multiplets.
            - unique_removed_multiplets (list): Pairs that were removed to resolve multiplets.
    """
    if hasattr(y_pred_np, 'cpu'): # Ensure y_pred_np is a NumPy array
        y_pred_np = y_pred_np.cpu().numpy()

    current_pred_pairs = [p[:2] for p in pred_pairs] # Work with [[i,j], ...] format
    removed_multiplet_pairs_accumulator = [] 

    while True:
        multiplet_groups = multiplets_pairs(current_pred_pairs)
        if not multiplet_groups: # No more multiplets to resolve
            break

        pairs_to_remove_this_iteration = set() # Use set to avoid duplicate removal instructions

        for group in multiplet_groups:
            if not group: continue

            min_prob = float('inf')
            pair_with_min_prob = None
            
            if len(group) <= 1: # Not a multiplet if only one pair (or empty group)
                continue

            for pair_in_group in group:
                idx0, idx1 = int(pair_in_group[0]), int(pair_in_group[1])
                current_prob = y_pred_np[idx0, idx1]
                if current_prob < min_prob:
                    min_prob = current_prob
                    pair_with_min_prob = pair_in_group
            
            if pair_with_min_prob:
                # Add tuple for hashability in set
                pairs_to_remove_this_iteration.add(tuple(sorted(pair_with_min_prob))) 
        
        if not pairs_to_remove_this_iteration: # No pair found to remove, break to prevent infinite loop
            break

        # Update current_pred_pairs
        new_pred_pairs_list = []
        for p_current in current_pred_pairs:
            if tuple(sorted(p_current)) not in pairs_to_remove_this_iteration:
                new_pred_pairs_list.append(p_current)
            else: # This pair is being removed
                removed_multiplet_pairs_accumulator.append(p_current)
        
        if len(new_pred_pairs_list) == len(current_pred_pairs): # No change in pairs, break
            break
        current_pred_pairs = new_pred_pairs_list
        
    unique_removed_multiplets = [list(x) for x in set(tuple(sorted(y)) for y in removed_multiplet_pairs_accumulator)]
    return current_pred_pairs, unique_removed_multiplets
        
def output_mask(seq: str, NC: bool = True) -> np.ndarray:
    """
    Generates a mask of allowed pairings based on sequence composition.
    For example, A can pair with U, G with C, G with U (wobble).
    If NC (Non-Canonical) is True, more pairs are allowed.

    Args:
        seq (str): RNA sequence string.
        NC (bool): If True, allows non-canonical pairs. Otherwise, only WC and GU/UG.

    Returns:
        np.ndarray: (L,L) NumPy array where 1 indicates an allowed pair type, 0 otherwise.
    """
    if NC: # Non-Canonical pairs allowed
        include_pairs = ['AU', 'UA', 'GC', 'CG', 'GU', 'UG', 'CC', 'GG', 'AG', 'CA', 'AC', 'UU', 'AA', 'CU', 'GA', 'UC']
    else: # Only Watson-Crick and Wobble G-U pairs
        include_pairs = ['AU', 'UA', 'GC', 'CG', 'GU', 'UG']
    
    seq_chars = [s.upper() for s in seq]
    L = len(seq_chars)
    mask = np.zeros((L, L), dtype=np.float32)
    for i in range(L):
        for j in range(L):
            if i == j: continue # Bases cannot pair with themselves in this context
            pair_str = seq_chars[i] + seq_chars[j]
            if pair_str in include_pairs:
                mask[i, j] = 1.0
    return mask

def ct_file_output(pairs: list, seq: str, rna_id: str, save_result_path: str) -> str:
    """Generates a .ct (Connectivity Table) file from predicted pairs."""
    if not os.path.exists(save_result_path):
        os.makedirs(save_result_path)

    L = len(seq)
    # .ct file columns:
    # 1: Index n (1 to L)
    # 2: Base character B_n
    # 3: Index n-1 (0-indexed for file, 0 for first base)
    # 4: Index n+1 (0-indexed for file, 0 for last base)
    # 5: Pairing partner for B_n (1-indexed, 0 if unpaired)
    # 6: Historical numbering (usually same as column 1)
    
    col1_idx_n = np.arange(1, L + 1, dtype=int)
    col2_base_char = np.array([s.upper() for s in seq])
    col3_idx_n_minus_1 = np.arange(0, L, dtype=int) 
    col4_idx_n_plus_1 = np.arange(2, L + 2, dtype=int)
    col4_idx_n_plus_1[L-1] = 0 # Last base's next partner is 0
    
    col5_pairing_partner = np.zeros(L, dtype=int)
    for p in pairs: # p is [idx0, idx1] (0-indexed)
        col5_pairing_partner[p[0]] = p[1] + 1 # Convert 0-indexed pair to 1-indexed partner
        col5_pairing_partner[p[1]] = p[0] + 1
    
    col6_historical_num = np.arange(1, L + 1, dtype=int)

    header_line = f"{L}\t\t{rna_id}\t\tSPOT-RNA PyTorch output\n" # Standard CT header
    data_lines = []
    for i in range(L):
        data_lines.append(f"{col1_idx_n[i]}\t{col2_base_char[i]}\t{col3_idx_n_minus_1[i]}\t{col4_idx_n_plus_1[i]}\t{col5_pairing_partner[i]}\t{col6_historical_num[i]}\n")
    
    ct_filepath = os.path.join(save_result_path, str(rna_id)) + '.ct'
    with open(ct_filepath, 'w') as f:
        f.write(header_line)
        for line_str in data_lines:
            f.write(line_str)
    return ct_filepath


def bpseq_file_output(pairs: list, seq: str, rna_id: str, save_result_path: str) -> str:
    """Generates a .bpseq (Base Pair Sequence) file from predicted pairs."""
    if not os.path.exists(save_result_path):
        os.makedirs(save_result_path)
        
    L = len(seq)
    # .bpseq file columns:
    # 1: Base index (1-indexed)
    # 2: Base character
    # 3: Pairing partner index (1-indexed, 0 if unpaired)
    
    col1_base_idx = np.arange(1, L + 1, dtype=int)
    col2_base_char = np.array([s.upper() for s in seq])
    col3_pairing_partner = np.zeros(L, dtype=int)

    for p in pairs: # p is [idx0, idx1] (0-indexed)
        col3_pairing_partner[p[0]] = p[1] + 1 # Convert 0-indexed to 1-indexed
        col3_pairing_partner[p[1]] = p[0] + 1
        
    # Original SPOT-RNA.py used a header like '#<id>' for .bpseq, but it was part of np.savetxt header arg,
    # which means it's a single comment line at the top.
    # Individual lines are just "idx base partner_idx".
    data_lines = []
    for i in range(L):
        data_lines.append(f"{col1_base_idx[i]} {col2_base_char[i]} {col3_pairing_partner[i]}\n")

    bpseq_filepath = os.path.join(save_result_path, str(rna_id)) + '.bpseq'
    with open(bpseq_filepath, 'w') as f:
        # Example: if a top-level comment is desired like original:
        # f.write(f"#{rna_id}\n") # This was the original format from np.savetxt header
        for line_str in data_lines:
            f.write(line_str)
    return bpseq_filepath


def lone_pair(pairs: list) -> list:
    """Identifies lone pairs (pairs that are not part of a stack)."""
    lone_pairs_list = []
    # For efficient lookup, convert pairs to a set of sorted tuples
    # E.g., {[0,10], [1,9]} becomes {(0,10), (1,9)} assuming pair[0] < pair[1]
    pair_set = set(tuple(sorted(p_item)) for p_item in pairs)

    for p_orig in pairs: # p_orig is [idx0, idx1]
        p_sorted = tuple(sorted(p_orig)) 
        
        # Check for a stack on the "left" side: (idx0-1, idx1+1)
        # (assuming p_sorted[0] < p_sorted[1])
        stacked_left_candidate = tuple(sorted((p_sorted[0] - 1, p_sorted[1] + 1)))
        # Check for a stack on the "right" side: (idx0+1, idx1-1)
        stacked_right_candidate = tuple(sorted((p_sorted[0] + 1, p_sorted[1] - 1)))
        
        is_lone = True
        if stacked_left_candidate in pair_set:
            is_lone = False
        if stacked_right_candidate in pair_set: # Can be part of two stacks
            is_lone = False
            
        if is_lone:
            lone_pairs_list.append(p_orig) # Add the original pair (not the sorted tuple)
            
    return lone_pairs_list


def prob_to_secondary_structure(
    ensemble_outputs_tensor: torch.Tensor | np.ndarray, 
    label_mask_tensor: torch.Tensor | np.ndarray, 
    seq_str: str, 
    rna_name: str, 
    args: object, 
    base_path: str
) -> list:
    """
    Converts a probability map (from model ensemble) to a predicted secondary structure.
    This involves thresholding, multiplet removal, and formatting into output files.

    Args:
        ensemble_outputs_tensor (torch.Tensor | np.ndarray): (L,L) tensor/array of averaged pair probabilities.
        label_mask_tensor (torch.Tensor | np.ndarray): (L,L) tensor/array mask (1s for allowed regions).
        seq_str (str): RNA sequence string.
        rna_name (str): Identifier for the RNA sequence.
        args (object): Namespace object containing command-line arguments (e.g., args.threshold, args.outputs).
        base_path (str): Base path of the project, used for locating external tools like VARNA.

    Returns:
        list: A list of predicted base pairs, e.g., [[0,10], [1,9], ...].
    """
    # --- 1. Ensure inputs are NumPy arrays ---
    if hasattr(ensemble_outputs_tensor, 'cpu'): # Check if it's a PyTorch tensor
        ensemble_outputs_np = ensemble_outputs_tensor.detach().cpu().numpy() # detach() if it might have grads
    else:
        ensemble_outputs_np = np.asarray(ensemble_outputs_tensor)

    if hasattr(label_mask_tensor, 'cpu'):
        label_mask_np = label_mask_tensor.detach().cpu().numpy()
    else:
        label_mask_np = np.asarray(label_mask_tensor)

    # --- 2. Apply masks to the probability map ---
    # y_pred_np is the (L,L) probability map after initial processing
    y_pred_np = ensemble_outputs_np 

    # Apply output_mask: filters by allowed base pair types (e.g., A-U, G-C, G-U, or more if NC allowed)
    allow_nc = args.allow_noncanonical if hasattr(args, 'allow_noncanonical') else True
    pairing_type_mask = output_mask(seq_str, NC=allow_nc) # (L,L) boolean or float mask
    y_pred_np = np.multiply(y_pred_np, pairing_type_mask)

    # Apply label_mask: filters by valid regions (e.g., upper triangle, no unknown bases involved)
    y_pred_np = np.multiply(y_pred_np, label_mask_np) # label_mask_np should be (L,L)

    # --- 3. Extract and filter pairs from the probability map ---
    L = len(seq_str)
    # Get indices for the upper triangle (k=1 excludes the main diagonal i==j)
    tri_upper_indices = np.triu_indices(L, k=1) 

    # Probabilities in the upper triangle
    probabilities_flat = y_pred_np[tri_upper_indices] 
    # Corresponding [i,j] coordinates for these probabilities
    pair_coordinates = [[tri_upper_indices[0][j], tri_upper_indices[1][j]] for j in range(len(tri_upper_indices[0]))]

    # Filter pairs by probability threshold and hairpin constraint
    binarized_pairs_above_threshold = []
    current_threshold = args.threshold if hasattr(args, 'threshold') else 0.335

    for idx, prob_value in enumerate(probabilities_flat):
        if prob_value >= current_threshold:
            pair_coord = pair_coordinates[idx]
            # Hairpin constraint: distance between paired bases must be > 3 (e.g. i and i+4 can pair for a loop of 3)
            if abs(pair_coord[0] - pair_coord[1]) > 3: 
                binarized_pairs_above_threshold.append(pair_coord)
    
    # --- 4. Resolve conflicting pairs (multiplets) ---
    # multiplets_free_bp takes the list of pairs and the full probability map to decide which pairs to keep.
    final_predicted_pairs, removed_multiplet_pairs = multiplets_free_bp(binarized_pairs_above_threshold, y_pred_np)
    
    # (Optional: The original code also had a `hair_pin_assumption` check here, but it's now integrated above)

    # --- 5. Generate Output Files ---
    # Determine output directory from args, relative to base_path if not absolute
    output_dir = args.outputs
    if not os.path.isabs(output_dir) and base_path: # If args.outputs is relative
        output_dir = os.path.join(base_path, output_dir)
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory for results: {output_dir}")

    # Save .ct (Connectivity Table) file
    ct_filepath = ct_file_output(final_predicted_pairs, seq_str, rna_name, output_dir)
    print(f"Saved .ct file: {ct_filepath}")

    # Save .bpseq (Base Pair Sequence) file
    bpseq_filepath = bpseq_file_output(final_predicted_pairs, seq_str, rna_name, output_dir)
    print(f"Saved .bpseq file: {bpseq_filepath}")

    # Save probability map (.prob file) - this is the (L,L) matrix after masking
    prob_filepath = os.path.join(output_dir, rna_name + '.prob')
    np.savetxt(prob_filepath, y_pred_np, delimiter='\t')
    print(f"Saved .prob (probability map) file: {prob_filepath}")
    
    # --- 6. Optional: VARNA plots and bpRNA motif analysis ---
    # These rely on external tools (Java for VARNA, Perl for bpRNA)

    if hasattr(args, 'plots') and args.plots:
        varna_jar_path = os.path.join(base_path, "utils", "VARNAv3-93.jar") # Assumes VARNA is in utils
        if not os.path.exists(varna_jar_path):
            print(f"Warning: VARNA jar not found at {varna_jar_path}. Skipping 2D plots.")
        else:
            try:
                png_radiate_path = os.path.join(output_dir, rna_name + '_radiate.png')
                png_line_path = os.path.join(output_dir, rna_name + '_line.png')
                
                # Basic VARNA command structure
                # Note: The original code had a complex `tertiary_bp` string for highlighting.
                # This is simplified here. If specific highlighting for multiplets/non-canonicals is needed,
                # the `tertiary_bp_varna_str` would need to be constructed based on `removed_multiplet_pairs`
                # and `type_pairs` output.
                common_varna_cmd = ["java", "-cp", varna_jar_path, "fr.orsay.lri.varna.applications.VARNAcmd",
                                    '-i', ct_filepath, '-resolution', '8.0', '-bpStyle', 'lw']
                
                subprocess.run(common_varna_cmd + ['-o', png_radiate_path, '-algorithm', 'radiate'], 
                               stderr=PIPE, stdout=PIPE, check=True) # check=True will raise on error
                print(f"Saved VARNA radiate plot: {png_radiate_path}")
                
                subprocess.run(common_varna_cmd + ['-o', png_line_path, '-algorithm', 'line'],
                               stderr=PIPE, stdout=PIPE, check=True)
                print(f"Saved VARNA line plot: {png_line_path}")

            except subprocess.CalledProcessError as e:
                print(f"Error generating VARNA plots for {rna_name}. Subprocess error: {e.stderr.decode()}")
            except FileNotFoundError: # Java command not found
                 print("Java command not found. Please ensure Java is installed and in system PATH to generate VARNA plots.")
            except Exception as e: # Other unexpected errors
                print(f"An unexpected error occurred during VARNA plotting for {rna_name}: {e}")

    if hasattr(args, 'motifs') and args.motifs:
        # Assumes bpRNA scripts are in utils/bpRNA-master/
        bpRNA_script_path = os.path.join(base_path, "utils", "bpRNA-master", "bpRNA.pl")
        if not os.path.exists(bpRNA_script_path):
            print(f"Warning: bpRNA script not found at {bpRNA_script_path}. Skipping motif analysis.")
        else:
            original_cwd = os.getcwd() # bpRNA might write to current working dir
            try:
                os.chdir(output_dir) # Change to output directory where .bpseq is
                
                # Command: perl <bpRNA_script_path> <input_bpseq_file>
                # bpRNA.pl typically produces a .st file (e.g., rna_name.st)
                subprocess.run(['perl', bpRNA_script_path, os.path.basename(bpseq_filepath)], 
                               stderr=PIPE, stdout=PIPE, check=True)
                
                st_filepath = os.path.join(output_dir, rna_name + '.st') 
                dbn_filepath = os.path.join(output_dir, rna_name + '.dbn') # Dot-bracket notation file

                if os.path.exists(st_filepath):
                    # Parse .st file to extract sequence and dot-bracket structure for .dbn file
                    # This is a simplified parser; original might use pandas.
                    seq_line_from_st, str_line_from_st = "", ""
                    with open(st_filepath) as f_st:
                        lines = f_st.readlines()
                        # Heuristic: find lines that look like sequence and structure
                        for line_content in lines:
                            line_strip = line_content.strip()
                            if line_strip.startswith(rna_name) or line_strip.startswith(">"): continue # Skip headers
                            if re.match(r"^[ACGUN]+$", line_strip, re.IGNORECASE): # Sequence line (allow N for robustness)
                                seq_line_from_st = line_strip
                            if re.match(r"^[().<>A-Za-z]+$", line_strip): # Structure line (DBN can have pseudoknot chars)
                                str_line_from_st = line_strip
                        
                    if seq_line_from_st and str_line_from_st:
                         with open(dbn_filepath, 'w') as f_dbn:
                            f_dbn.write(f">{rna_name} (from bpRNA)\n") # Header for DBN
                            f_dbn.write(f"{seq_line_from_st}\n")
                            f_dbn.write(f"{str_line_from_st}\n")
                         print(f"Saved .dbn file (from bpRNA): {dbn_filepath}")
                    else:
                        print(f"Could not reliably parse sequence/structure from .st file: {st_filepath}")
                else:
                    print(f"bpRNA output .st file not found: {st_filepath}")

            except subprocess.CalledProcessError as e:
                print(f"Error running bpRNA for {rna_name}. Subprocess error: {e.stderr.decode()}")
            except FileNotFoundError: # Perl command not found
                 print("Perl command not found. Please ensure Perl is installed and in system PATH for bpRNA analysis.")
            except Exception as e: # Other unexpected errors
                print(f"An unexpected error occurred during bpRNA processing for {rna_name}: {e}")
            finally:
                os.chdir(original_cwd) # Important: change back to original directory
                
    return final_predicted_pairs
