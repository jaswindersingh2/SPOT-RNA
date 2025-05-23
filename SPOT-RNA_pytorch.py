import torch
import torch.nn as nn 
from torch.utils.data import Dataset, DataLoader
import os 
import numpy as np 
import argparse 
from argparse import RawTextHelpFormatter 

# Import PyTorch-specific utility functions
from utils.pytorch_utils import get_data, FastaMLtoSL, prob_to_secondary_structure

class RNADataset(Dataset):
    """
    PyTorch Dataset class for loading RNA sequences from a FASTA file.
    Handles FASTA file processing (converting multi-line to single-line) and
    preprocessing of sequences into tensors suitable for model input.
    """
    def __init__(self, fasta_file_path: str, processed_fasta_dir: str = "processed_fasta_temp"):
        """
        Initializes the RNADataset.

        Args:
            fasta_file_path (str): Path to the input FASTA file.
            processed_fasta_dir (str, optional): Directory to store intermediate single-line FASTA files.
                                                 Defaults to "processed_fasta_temp".
        """
        super().__init__()
        
        if not os.path.exists(processed_fasta_dir):
            os.makedirs(processed_fasta_dir)
            # print(f"Created temporary directory for processed FASTA: {processed_fasta_dir}")

        base_name = os.path.basename(fasta_file_path)
        name, ext = os.path.splitext(base_name)
        # Path for the (potentially) single-line version of the input FASTA
        processed_file_path = os.path.join(processed_fasta_dir, f"{name}_sl{ext}")

        try:
            # Convert input FASTA to single-line format for easier parsing
            actual_processed_path = FastaMLtoSL(fasta_file_path, processed_file_path)
            # print(f"FASTA file '{fasta_file_path}' processed to '{actual_processed_path}' for loading.")
        except Exception as e:
            print(f"Error processing FASTA file {fasta_file_path} with FastaMLtoSL: {e}")
            print(f"Warning: FastaMLtoSL failed. Attempting to read original file: {fasta_file_path}")
            actual_processed_path = fasta_file_path # Fallback to original path

        self.rna_data = []  # Stores tuples of (rna_id, sequence_string)
        self.rna_details_map = {} # Stores {rna_id: {'sequence': str, 'label_mask': Tensor, 'seq_len': int}}
                                  # This map is populated by __getitem__ for later use in post-processing.

        try:
            with open(actual_processed_path, 'r') as f:
                current_id = None
                current_seq_parts = []
                for line in f:
                    line = line.strip()
                    if line.startswith(">"): # Start of a new sequence entry
                        if current_id and current_seq_parts: # Save previous sequence if exists
                            full_sequence = "".join(current_seq_parts)
                            self.rna_data.append((current_id, full_sequence))
                        current_id = line[1:] # Store ID without ">"
                        current_seq_parts = []
                    elif current_id and line: # Sequence line for current_id
                        current_seq_parts.append(line.upper().replace("T", "U")) # Ensure uppercase and RNA bases
                if current_id and current_seq_parts: # Save the last sequence in the file
                    full_sequence = "".join(current_seq_parts)
                    self.rna_data.append((current_id, full_sequence))
            
            if not self.rna_data:
                print(f"Warning: No data read from {actual_processed_path}. Ensure it's a valid FASTA file.")

        except FileNotFoundError:
            print(f"Error: FASTA file not found: {actual_processed_path}. Cannot load dataset.")
            self.rna_data = [] 
        except Exception as e:
            print(f"An error occurred reading sequences from {actual_processed_path}: {e}")
            self.rna_data = []

    def __len__(self) -> int:
        """Returns the total number of RNA sequences in the dataset."""
        return len(self.rna_data)

    def __getitem__(self, idx: int) -> dict | None:
        """
        Retrieves and preprocesses a single RNA sequence from the dataset.

        Args:
            idx (int): Index of the RNA sequence to retrieve.

        Returns:
            dict | None: A dictionary containing the processed data for the sequence:
                         {'id', 'sequence', 'seq_len', 'features', 'zero_mask', 'label_mask'}.
                         Returns None if an error occurs during processing.
        """
        rna_id, sequence_str = self.rna_data[idx]
        
        try:
            # Preprocess sequence using get_data from pytorch_utils
            # Returns: seq_len, features_tensor, zero_mask_tensor, label_mask_tensor
            seq_len, features, zero_mask, label_mask = get_data(sequence_str)
            
            # Store details (original sequence, label_mask, seq_len) needed for post-processing (e.g., structure prediction)
            # This map is populated on-the-fly as items are requested.
            if rna_id not in self.rna_details_map:
                 self.rna_details_map[rna_id] = {
                    "sequence": sequence_str,    # Original sequence string
                    "label_mask": label_mask,    # (L,L) tensor from get_data, before padding
                    "seq_len": seq_len           # Original length
                }

            return {
                "id": rna_id,
                "sequence": sequence_str, # Original sequence string
                "seq_len": seq_len,       # Original length
                "features": features,     # (L, L, 8) tensor
                "zero_mask": zero_mask,   # (1, L, L) tensor
                "label_mask": label_mask  # (L, L) tensor (from get_data, before padding)
            }
        except Exception as e:
            print(f"Error processing sequence {rna_id} (index {idx}) in __getitem__: {e}")
            return None # DataLoader's collate_fn should handle None items

def rna_collate_fn(batch: list) -> dict | None:
    """
    Custom collate function for the DataLoader.
    Handles batches of RNA sequence data, which may have varying lengths.
    It pads sequences to the maximum length in the batch and stacks them into tensors.
    Filters out None items that may result from errors in RNADataset.__getitem__.

    Args:
        batch (list): A list of dictionaries, where each dictionary is an output from RNADataset.__getitem__.

    Returns:
        dict | None: A dictionary of batched and padded tensors:
                     {'ids', 'sequences_str', 'seq_lens', 'features', 'zero_masks', 'label_masks'}.
                     Returns None if the batch is empty after filtering None items.
    """
    batch = [item for item in batch if item is not None] # Filter out items that failed processing in __getitem__
    if not batch: # If all items in batch failed, return None
        return None 

    # Determine the maximum sequence length in the current batch for padding
    max_len = 0
    for item in batch:
        if item["seq_len"] > max_len:
            max_len = item["seq_len"]
            
    # Initialize lists to hold data for the batch
    ids = []
    sequences_str = [] 
    seq_lens = [] # Store original sequence lengths
    
    padded_features = []
    padded_zero_masks = []
    padded_label_masks = [] 

    # Pad each item's tensors to max_len and collect
    for item in batch:
        ids.append(item["id"])
        sequences_str.append(item["sequence"]) 
        seq_lens.append(item["seq_len"])
        
        current_len = item["seq_len"]
        
        # Padding for 'features' tensor (L, L, 8) to (max_L, max_L, 8)
        # Pad format: (pad_left, pad_right, pad_top, pad_bottom, pad_front, pad_back) for each dim.
        # For (L,L,8), dim0 (rows), dim1 (cols), dim2 (channels)
        # PyTorch pad is (pad_left, pad_right, pad_top, pad_bottom) for 4D input (N,C,H,W)
        # For 3D (C,H,W) or (H,W,C), need to be careful.
        # torch.nn.functional.pad expects ((pad_dimN_start, pad_dimN_end), ...) for last N dims
        # Features (L,L,8) -> pad last 3 dims. (Actually, only first 2 dims (H,W-like) for features)
        # Pad ( (0, max_len - current_len) for rows, (0, max_len - current_len) for columns, (0,0) for channels )
        pad_dims_features = (0, 0, 0, max_len - current_len, 0, max_len - current_len) # (pad_C_start,end, pad_W_start,end, pad_H_start,end)
        # Correcting pad for (H,W,C) like (L,L,8) - pad needs to be for last dim first
        # (pad_dim2_start, pad_dim2_end, pad_dim1_start, pad_dim1_end, pad_dim0_start, pad_dim0_end)
        # item["features"] is (L,L,8). We want to pad H and W dimensions.
        # So, (pad_cols_start, pad_cols_end, pad_rows_start, pad_rows_end) if thinking of (C,H,W)
        # If (H,W,C), pad is ((0,0) for C, (0,max_L-L) for W, (0,max_L-L) for H)
        # features are (current_len, current_len, 8)
        # padding for (dim -1), then (dim -2), etc.
        # No padding for channel dim (last one, index 2)
        # Pad width format for np.pad: ((before_axis0, after_axis0), (before_axis1, after_axis1), ...)
        # torch.nn.functional.pad input is (N, C, H, W) or (C, H, W) or (H, W)
        # item["features"] is (H, W, C) = (L, L, 8). Permute to (C, H, W) = (8, L, L) for pad, then permute back.
        features_permuted = item["features"].permute(2,0,1) # (8, L, L)
        pad_for_torch_L_L = (0, max_len - current_len, 0, max_len - current_len) # (pad_W, pad_H)
        padded_f_permuted = torch.nn.functional.pad(features_permuted, pad_for_torch_L_L, mode='constant', value=0)
        padded_f = padded_f_permuted.permute(1,2,0) # Back to (max_L, max_L, 8)
        padded_features.append(padded_f)

        # Padding for 'zero_mask' (1, L, L) to (1, max_L, max_L)
        # item["zero_mask"] is (C, H, W) = (1, L, L)
        # Pad H, W dimensions. Pad format is (pad_W_start,end, pad_H_start,end)
        padded_zm = torch.nn.functional.pad(item["zero_mask"], pad_for_torch_L_L, mode='constant', value=0) 
        padded_zero_masks.append(padded_zm)
        
        # Padding for 'label_mask' (L, L) to (max_L, max_L)
        # item["label_mask"] is (H,W)
        padded_lm = torch.nn.functional.pad(item["label_mask"], pad_for_torch_L_L, mode='constant', value=0)
        padded_label_masks.append(padded_lm)

    # Stack individual padded tensors to create batch tensors
    collated_batch = {
        "ids": ids,
        "sequences_str": sequences_str, # List of original sequence strings
        "seq_lens": torch.tensor(seq_lens, dtype=torch.long), # Tensor of original lengths
        "features": torch.stack(padded_features),       # Shape: (batch_size, max_L, max_L, 8)
        "zero_masks": torch.stack(padded_zero_masks),   # Shape: (batch_size, 1, max_L, max_L)
        "label_masks": torch.stack(padded_label_masks)  # Shape: (batch_size, max_L, max_L)
    }
    
    return collated_batch

class PlaceholderSPOTRNA(nn.Module):
    """
    A placeholder PyTorch model for SPOT-RNA.
    This model does not perform any actual learning or complex RNA structure prediction.
    Its purpose is to provide a structural stand-in for a real model within the pipeline,
    allowing the data loading, batching, and post-processing steps to be tested.
    The forward pass generates an output tensor of the expected shape based on input features.
    """
    def __init__(self):
        super(PlaceholderSPOTRNA, self).__init__()
        # No layers defined as this is a placeholder. A real model would define PyTorch layers here.
        # e.g., self.conv_layer = nn.Conv2d(...)

    def forward(self, batch_data: dict) -> torch.Tensor:
        """
        Placeholder forward pass.
        Generates a dummy output tensor with the shape expected by downstream processing.
        (batch_size, max_sequence_length, max_sequence_length).

        Args:
            batch_data (dict): A dictionary containing batched data, including 'features'.
                               'features' has shape (batch_size, max_L, max_L, 8).

        Returns:
            torch.Tensor: A placeholder output tensor, typically derived from input features
                          to match the expected output shape for pair probabilities.
                          Shape: (batch_size, max_L, max_L).
        """
        # Example: sum input features over the channel dimension (8) to get (batch_size, max_L, max_L)
        # This is an arbitrary operation to produce an output of the correct shape.
        placeholder_output = torch.sum(batch_data["features"], dim=3) 
        return placeholder_output


if __name__ == '__main__':
    # --- Argument Parser Setup ---
    # Defines command-line arguments for the script, similar to the original SPOT-RNA.py.
    parser = argparse.ArgumentParser(
        description="SPOT-RNA PyTorch: Predict RNA secondary structure.",
        formatter_class=RawTextHelpFormatter
    )
    parser.add_argument('--inputs', required=True, type=str, 
                        help='Path to input file in FASTA format (required). Accepts multiple sequences.', metavar='<file.fasta>')
    parser.add_argument('--outputs',default='outputs_pytorch/', type=str, 
                        help='Path to output files directory.\n(default: %(default)s)', metavar='<dir>')
    parser.add_argument('--gpu', default=-1, type=int, 
                        help='Specify GPU number (e.g., 0 or 1). Use -1 for CPU.\n(default: %(default)s, for CPU)', metavar='<gpu_id>')
    parser.add_argument('--plots',default=False, type=lambda x: (str(x).lower() == 'true'),
                        help='Set to "True" to generate 2D plots (requires Java & VARNA).\n(default: %(default)s)', metavar='<True|False>')
    parser.add_argument('--motifs',default=False, type=lambda x: (str(x).lower() == 'true'),
                        help='Set to "True" to generate motifs (requires Perl & bpRNA).\n(default: %(default)s)', metavar='<True|False>')
    parser.add_argument('--cpu',default=1, type=int, # Default to 1 CPU for PyTorch unless specified
                        help='Number of CPU threads for PyTorch.\n(default: %(default)s)', metavar='<threads>')
    parser.add_argument('--threshold', default=0.335, type=float, 
                        help='Probability threshold for selecting base pairs in post-processing.\n(default: %(default)s)', metavar='<float>')
    parser.add_argument('--allow_noncanonical', default=True, type=lambda x: (str(x).lower() == 'true'),
                        help='Allow non-canonical base pairs during structure prediction.\n(default: %(default)s)', metavar='<True|False>')
    
    args = parser.parse_args()

    # --- Setup base_path and output/device configurations ---
    base_path = os.path.dirname(os.path.realpath(__file__)) # Path to the directory containing this script
    
    if not os.path.exists(args.outputs):
        os.makedirs(args.outputs)
        print(f"Created output directory: {args.outputs}")

    print(f"SPOT-RNA PyTorch version starting...")
    print(f"  Base path: {base_path}")
    print(f"  Input FASTA file: {args.inputs}")
    print(f"  Output directory: {args.outputs}")

    # GPU/CPU settings
    if args.gpu != -1 and torch.cuda.is_available():
        device = torch.device(f"cuda:{args.gpu}")
        print(f"  Using GPU: {args.gpu}")
    else:
        if args.gpu != -1:
            print(f"  Warning: GPU {args.gpu} specified, but CUDA is not available. Using CPU.")
        device = torch.device("cpu")
        print("  Using CPU.")
    
    if args.cpu:
        torch.set_num_threads(args.cpu) # Set number of threads for CPU operations
        # print(f"  Attempting to set PyTorch CPU threads to: {args.cpu}")


    # Directory for intermediate processed FASTA files (e.g., single-line versions)
    processed_fasta_intermediate_dir = os.path.join(args.outputs, "processed_fasta_temp")

    # --- RNADataset and DataLoader instantiation ---
    print(f"\n--- Initializing RNADataset with input file: {args.inputs} ---")
    rna_dataset = RNADataset(fasta_file_path=args.inputs, processed_fasta_dir=processed_fasta_intermediate_dir)
    
    if not len(rna_dataset):
        print(f"RNADataset is empty. Please check the input file: {args.inputs}. Exiting.")
        if os.path.exists(processed_fasta_intermediate_dir) and not os.listdir(processed_fasta_intermediate_dir):
            try: os.rmdir(processed_fasta_intermediate_dir) # Clean up empty dir
            except OSError: pass 
        exit()
        
    print(f"RNADataset initialized. Number of sequences found: {len(rna_dataset)}")

    # --- Ensemble Model Setup ---
    # This section simulates an ensemble of models. In a real scenario, pre-trained model
    # weights would be loaded here. Currently, it uses multiple instances of PlaceholderSPOTRNA.
    num_models_ensemble = 5 
    ensemble_models = []
    for i in range(num_models_ensemble):
        model = PlaceholderSPOTRNA()
        # Placeholder for loading actual model weights:
        # model_path = os.path.join(base_path, 'models', f'spot_rna_pytorch_model_{i}.pth')
        # if os.path.exists(model_path):
        #     model.load_state_dict(torch.load(model_path, map_location=device))
        # else:
        #     print(f"Warning: Model weights not found at {model_path}. Using placeholder model {i}.")
        model.to(device) # Move model to the selected device (CPU or GPU)
        ensemble_models.append(model)
    
    print(f"\nInitialized {num_models_ensemble} placeholder models for the ensemble on device: {device}.")

    # Dictionary to store raw model outputs for each sequence from each model
    outputs_ensemble_dict = {} 
    # Dictionary to store sequence details (sequence string, label_mask, original length)
    # This is populated from rna_dataset.rna_details_map after data iteration.
    rna_details_for_postprocessing = {} 

    # --- Inference Loop for each model in the ensemble ---
    for model_idx, current_model in enumerate(ensemble_models):
        print(f"\n--- Processing with Model {model_idx + 1}/{num_models_ensemble} ---")
        current_model.eval() # Set model to evaluation mode

        # DataLoader for the current model. shuffle=False is typical for inference.
        # num_workers can be > 0 for performance, but 0 is simpler for debugging.
        rna_dataloader = DataLoader(rna_dataset, batch_size=2, shuffle=False, collate_fn=rna_collate_fn, num_workers=0) 
        
        # print(f"DataLoader ready for model {model_idx + 1}. Iterating through batches...")
        
        first_batch_processed_this_model = False
        for batch_num, batch in enumerate(rna_dataloader):
            if batch is None: # Skip if batch is None (e.g., all items failed in __getitem__)
                # print(f"Batch {batch_num + 1} was None. Skipping.")
                continue

            # Move batch data to the selected device
            features_batch = batch["features"].to(device)
            # zero_masks_batch = batch["zero_masks"].to(device) # If used by model
            # label_masks_batch = batch["label_masks"].to(device) # If used by model
            
            # Create a new batch dictionary for the model with tensors on the correct device
            model_input_batch = {
                "features": features_batch,
                "seq_lens": batch["seq_lens"] # seq_lens usually on CPU
                # Add other tensors if your model expects them:
                # "zero_masks": zero_masks_batch,
                # "label_masks": label_masks_batch,
            }

            with torch.no_grad(): # Disable gradient calculations during inference
                model_output = current_model(model_input_batch) # Output shape: (batch_size, max_L, max_L)
            
            # Apply sigmoid to get probabilities
            sigmoid_output_batch = torch.sigmoid(model_output) 

            if not first_batch_processed_this_model and model_idx == 0: # Print shape info only for the first batch of the first model
                print(f"  Output shape for first batch (Model 1, raw model output): {model_output.shape}")
                print(f"  Sigmoid output shape for first batch (Model 1): {sigmoid_output_batch.shape}")
                first_batch_processed_this_model = True 

            # Store sliced outputs (per sequence) for the current model
            batch_ids = batch["ids"]
            batch_seq_lens = batch["seq_lens"] # Original lengths

            for i in range(sigmoid_output_batch.shape[0]): # Iterate through items in the batch
                rna_id = batch_ids[i]
                actual_len = batch_seq_lens[i].item() 
                
                # Slice the output from (max_L, max_L) to actual (L, L)
                # Move to CPU for storage if it was on GPU, ensures NumPy compatibility later
                sliced_sigmoid_output = sigmoid_output_batch[i, :actual_len, :actual_len].cpu() 
                outputs_ensemble_dict.setdefault(rna_id, []).append(sliced_sigmoid_output)
        
    # Populate rna_details_for_postprocessing from the dataset's map
    # This map is filled by __getitem__ as data items are requested by DataLoader.
    rna_details_for_postprocessing = rna_dataset.rna_details_map

    # --- Averaging Ensemble Outputs and Generating Secondary Structures ---
    print("\n--- Averaging Ensemble Outputs and Generating Secondary Structures ---")
    if not outputs_ensemble_dict:
        print("No ensemble outputs were collected. Skipping structure prediction.")
    else:
        for rna_id, list_of_pred_tensors in outputs_ensemble_dict.items():
            if not list_of_pred_tensors:
                print(f"No predictions found for RNA ID: {rna_id}. Skipping.")
                continue
            
            # Stack all model predictions for this sequence and average them
            try:
                # Each tensor in list_of_pred_tensors is (L,L) and on CPU
                stacked_preds = torch.stack(list_of_pred_tensors, dim=0) # Shape: (num_models, L, L)
                averaged_preds_tensor = torch.mean(stacked_preds, dim=0) # Shape: (L, L)
            except Exception as e:
                print(f"Error stacking/averaging predictions for {rna_id}: {e}")
                continue

            # print(f"\nProcessing RNA ID: {rna_id}")
            # print(f"  Averaged predictions shape: {averaged_preds_tensor.shape}")

            # Retrieve the original sequence string and label_mask for this rna_id
            details = rna_details_for_postprocessing.get(rna_id)
            if not details:
                print(f"  Error: Could not find sequence details (sequence, label_mask) for {rna_id}. Skipping structure prediction.")
                continue

            sequence_str = details["sequence"]
            label_mask_tensor = details["label_mask"] # This is the (L,L) tensor from get_data (on CPU)
            
            try:
                # Call the main post-processing function
                predicted_pairs = prob_to_secondary_structure(
                    ensemble_outputs_tensor=averaged_preds_tensor, # (L,L) averaged probability map
                    label_mask_tensor=label_mask_tensor,           # (L,L) mask from get_data
                    seq_str=sequence_str,
                    rna_name=rna_id,
                    args=args, # Parsed command-line arguments
                    base_path=base_path 
                )
                # print(f"  Predicted {len(predicted_pairs)} pairs for {rna_id}.")

            except Exception as e:
                print(f"  Error during prob_to_secondary_structure for {rna_id}: {e}")
                import traceback
                traceback.print_exc()
    
    # --- Cleanup ---
    # Clean up the temporary directory used for processed FASTA files
    try:
        if os.path.exists(processed_fasta_intermediate_dir):
            import shutil
            shutil.rmtree(processed_fasta_intermediate_dir)
            # print(f"Cleaned up temporary processed FASTA directory: {processed_fasta_intermediate_dir}")
    except Exception as e:
        print(f"Warning: Could not clean up temporary directory {processed_fasta_intermediate_dir}: {e}")

    print("\nSPOT-RNA PyTorch script execution finished.")
