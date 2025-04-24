#!/usr/bin/env python3
"""
GPU-accelerated Gromark Cipher Solver for Kryptos K4
This script provides a CUDA implementation for testing large Gromark primer spaces.
Requires: CUDA toolkit and PyCUDA
"""

import argparse
import time
import string
import numpy as np
import math
import logging
from gromark_solver import score_text, print_results, K4_CIPHERTEXT, CRIB1, CRIB2, PROMISING_PRIMERS

try:
    import pycuda.driver as cuda
    import pycuda.autoinit
    from pycuda.compiler import SourceModule
    PYCUDA_AVAILABLE = True
except ImportError:
    PYCUDA_AVAILABLE = False
    print("Warning: PyCUDA not available. GPU acceleration will not work.")
    print("Install with: pip install pycuda")

# Optional tqdm import for progress bars
try:
    from tqdm import tqdm
    TQDM_AVAILABLE = True
except ImportError:
    TQDM_AVAILABLE = False
    
# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# CUDA kernel for Gromark cipher decryption and crib checking
CUDA_KERNEL = """
// Use unsigned char for clarity when dealing with text
typedef unsigned char uchar;

__device__ void expand_gromark_primer(int* primer, int primer_length, int base, int* key, int key_length) {
    // Copy primer to the beginning of the key
    for (int i = 0; i < primer_length; i++) {
        key[i] = primer[i];
    }
    
    // Expand the key
    for (int i = primer_length; i < key_length; i++) {
        key[i] = (key[i-1] + key[i-2]) % base;
    }
}

__device__ bool check_cribs(uchar* plaintext, 
                          const uchar* crib1, int crib1_start, int crib1_end,
                          const uchar* crib2, int crib2_start, int crib2_end) {
    // Check first crib
    for (int i = crib1_start; i <= crib1_end; i++) {
        if (plaintext[i] != crib1[i - crib1_start]) {
            return false;
        }
    }
    
    // Check second crib
    for (int i = crib2_start; i <= crib2_end; i++) {
        if (plaintext[i] != crib2[i - crib2_start]) {
            return false;
        }
    }
    
    return true;
}

__global__ void test_primers(int* primers, int num_primers, int primer_length, int base,
                           uchar* ciphertext, int ct_length,
                           const uchar* cipher_alphabet, int alphabet_length,
                           const uchar* crib1, int crib1_start, int crib1_end,
                           const uchar* crib2, int crib2_start, int crib2_end,
                           int* result_indices, int* num_results) {
    
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_primers) return;
    
    // Get the primer to test
    int primer[10];  // Assuming max primer length of 10
    for (int i = 0; i < primer_length; i++) {
        primer[i] = primers[idx * primer_length + i];
    }
    
    // Expand the primer into a key
    int key[500];  // Assuming max ciphertext length of 500
    expand_gromark_primer(primer, primer_length, base, key, ct_length);
    
    // Decrypt the ciphertext
    uchar plaintext[500];  // Assuming max ciphertext length of 500
    for (int i = 0; i < ct_length; i++) {
        // Find the index of the ciphertext character in the cipher alphabet
        int ct_index = -1;
        for (int j = 0; j < alphabet_length; j++) {
            if (ciphertext[i] == cipher_alphabet[j]) {
                ct_index = j;
                break;
            }
        }
        
        if (ct_index != -1) {
            // Decrypt the character: (ciphertext_index - key_digit) mod alphabet_length
            int pt_index = (ct_index - key[i] + alphabet_length) % alphabet_length;
            plaintext[i] = cipher_alphabet[pt_index];  // Using same alphabet for plaintext
        } else {
            plaintext[i] = ciphertext[i];  // Non-alphabet character
        }
    }
    
    // Check if the decrypted text contains the cribs
    if (check_cribs(plaintext, crib1, crib1_start, crib1_end, crib2, crib2_start, crib2_end)) {
        // Found a match, add the primer index to results
        int result_idx = atomicAdd(num_results, 1);
        if (result_idx < 1000) {  // Limit to 1000 results
            result_indices[result_idx] = idx;
        }
    }
}
"""

def generate_primers_batch(base, length, start_idx, batch_size):
    """
    Generate a batch of primers as a NumPy array for GPU processing.
    
    Args:
        base (int): The number base to use
        length (int): Length of primers to generate
        start_idx (int): Starting index in the primer space
        batch_size (int): Number of primers to generate in this batch
        
    Returns:
        numpy.ndarray: Array of primers, shape (batch_size, length)
    """
    # Generate a range of indices for this batch
    indices = np.arange(start_idx, start_idx + batch_size, dtype=np.int32)
    
    # Convert indices to base-N representation
    primers = np.zeros((batch_size, length), dtype=np.int32)
    for i in range(length):
        primers[:, length-1-i] = (indices // (base ** i)) % base
    
    return primers.reshape(-1)  # Flatten to 1D array in row-major order

def gpu_brute_force_gromark(base=10, primer_length=5, 
                          ciphertext=K4_CIPHERTEXT,
                          cipher_alphabet=string.ascii_uppercase,
                          max_batch_size=1000000,
                          specific_primers=None):
    """
    Use GPU to brute force Gromark primers and find those that match the cribs.
    
    Args:
        base (int): The number base to use
        primer_length (int): Length of primers to test
        ciphertext (str): The ciphertext to decrypt
        cipher_alphabet (str): Ciphertext alphabet
        max_batch_size (int): Maximum number of primers to test in a single batch
        specific_primers (list, optional): List of specific primers to try first
        
    Returns:
        list: List of matching primer indices
    """
    if not PYCUDA_AVAILABLE:
        logger.error("PyCUDA is required for GPU acceleration")
        return []
    
    # Prepare the CUDA kernel
    mod = SourceModule(CUDA_KERNEL)
    test_primers_kernel = mod.get_function("test_primers")
    
    # Prepare data for the kernel
    crib1_text, (crib1_start, crib1_end) = CRIB1
    crib2_text, (crib2_start, crib2_end) = CRIB2
    
    # Prepare ciphertext
    ciphertext_np = np.array([ord(c) for c in ciphertext], dtype=np.uint8)
    ciphertext_gpu = cuda.mem_alloc(ciphertext_np.nbytes)
    cuda.memcpy_htod(ciphertext_gpu, ciphertext_np)
    
    # Prepare cipher alphabet
    cipher_alphabet_np = np.array([ord(c) for c in cipher_alphabet], dtype=np.uint8)
    cipher_alphabet_gpu = cuda.mem_alloc(cipher_alphabet_np.nbytes)
    cuda.memcpy_htod(cipher_alphabet_gpu, cipher_alphabet_np)
    
    # Prepare cribs
    crib1_np = np.array([ord(c) for c in crib1_text], dtype=np.uint8)
    crib1_gpu = cuda.mem_alloc(crib1_np.nbytes)
    cuda.memcpy_htod(crib1_gpu, crib1_np)
    
    crib2_np = np.array([ord(c) for c in crib2_text], dtype=np.uint8)
    crib2_gpu = cuda.mem_alloc(crib2_np.nbytes)
    cuda.memcpy_htod(crib2_gpu, crib2_np)
    
    # Prepare result arrays
    max_results = 1000
    result_indices_np = np.zeros(max_results, dtype=np.int32)
    result_indices_gpu = cuda.mem_alloc(result_indices_np.nbytes)
    cuda.memcpy_htod(result_indices_gpu, result_indices_np)
    
    num_results_np = np.zeros(1, dtype=np.int32)
    num_results_gpu = cuda.mem_alloc(num_results_np.nbytes)
    cuda.memcpy_htod(num_results_gpu, num_results_np)
    
    all_result_indices = []
    
    # First, test any specific primers if provided
    if specific_primers:
        logger.info(f"Testing {len(specific_primers)} specific primers...")
        
        # Convert primers to numpy array
        primers_np = np.array([digit for primer in specific_primers for digit in primer], dtype=np.int32)
        primers_gpu = cuda.mem_alloc(primers_np.nbytes)
        cuda.memcpy_htod(primers_gpu, primers_np)
        
        # Reset result counters
        num_results_np.fill(0)
        cuda.memcpy_htod(num_results_gpu, num_results_np)
        
        # Calculate grid and block dimensions
        block_size = 256
        grid_size = (len(specific_primers) + block_size - 1) // block_size
        
        # Launch the kernel
        test_primers_kernel(
            primers_gpu, np.int32(len(specific_primers)), np.int32(primer_length), np.int32(base),
            ciphertext_gpu, np.int32(len(ciphertext)),
            cipher_alphabet_gpu, np.int32(len(cipher_alphabet)),
            crib1_gpu, np.int32(crib1_start), np.int32(crib1_end),
            crib2_gpu, np.int32(crib2_start), np.int32(crib2_end),
            result_indices_gpu, num_results_gpu,
            block=(block_size, 1, 1), grid=(grid_size, 1)
        )
        
        # Get the results
        cuda.memcpy_dtoh(num_results_np, num_results_gpu)
        cuda.memcpy_dtoh(result_indices_np, result_indices_gpu)
        
        num_matches = int(num_results_np[0])
        if num_matches > 0:
            logger.info(f"Found {num_matches} matches from specific primers")
            # These are indices into the specific_primers list
            for i in range(min(num_matches, max_results)):
                idx = result_indices_np[i]
                if idx < len(specific_primers):
                    logger.info(f"Match found for primer: {specific_primers[idx]}")
                    all_result_indices.append((specific_primers[idx], -1))  # -1 indicates a specific primer
    
    # Then test the entire primer space in batches
    total_primers = base ** primer_length
    logger.info(f"Testing up to {total_primers} primers (base {base}, length {primer_length})...")
    
    # Determine number of batches
    batch_size = min(max_batch_size, total_primers)
    num_batches = (total_primers + batch_size - 1) // batch_size
    
    # Prepare memory for a batch of primers
    primers_gpu = cuda.mem_alloc(batch_size * primer_length * np.dtype(np.int32).itemsize)
    
    overall_start_time = time.time()
    
    # Create batch processing loop
    batch_range = tqdm(range(num_batches)) if TQDM_AVAILABLE else range(num_batches)
    for batch_idx in batch_range:
        batch_start_time = time.time()
        
        # Calculate actual batch size (last batch may be smaller)
        start_primer_idx = batch_idx * batch_size
        actual_batch_size = min(batch_size, total_primers - start_primer_idx)
        
        # Generate batch of primers
        primers_np = generate_primers_batch(base, primer_length, start_primer_idx, actual_batch_size)
        cuda.memcpy_htod(primers_gpu, primers_np)
        
        # Reset result counters
        num_results_np.fill(0)
        cuda.memcpy_htod(num_results_gpu, num_results_np)
        
        # Calculate grid and block dimensions
        block_size = 256
        grid_size = (actual_batch_size + block_size - 1) // block_size
        
        # Launch the kernel
        test_primers_kernel(
            primers_gpu, np.int32(actual_batch_size), np.int32(primer_length), np.int32(base),
            ciphertext_gpu, np.int32(len(ciphertext)),
            cipher_alphabet_gpu, np.int32(len(cipher_alphabet)),
            crib1_gpu, np.int32(crib1_start), np.int32(crib1_end),
            crib2_gpu, np.int32(crib2_start), np.int32(crib2_end),
            result_indices_gpu, num_results_gpu,
            block=(block_size, 1, 1), grid=(grid_size, 1)
        )
        
        # Get the results
        cuda.memcpy_dtoh(num_results_np, num_results_gpu)
        cuda.memcpy_dtoh(result_indices_np, result_indices_gpu)
        
        num_matches = int(num_results_np[0])
        
        # Process any matches found in this batch
        if num_matches > 0:
            batch_end_time = time.time()
            logger.info(f"Batch {batch_idx+1}/{num_batches}: Found {num_matches} matches in {batch_end_time - batch_start_time:.2f} seconds")
            
            # Convert result indices to actual primer values and add to list
            for i in range(min(num_matches, max_results)):
                idx = result_indices_np[i]
                if idx < actual_batch_size:
                    primer_idx = start_primer_idx + idx
                    all_result_indices.append((None, primer_idx))  # None indicates we need to reconstruct from idx
        else:
            if TQDM_AVAILABLE:
                # Update progress bar description if using tqdm
                if isinstance(batch_range, tqdm):
                    batch_range.set_description(f"Batch {batch_idx+1}/{num_batches}")
            else:
                # Simple progress logging otherwise
                if batch_idx % 10 == 0 or batch_idx == num_batches - 1:
                    batch_end_time = time.time()
                    batch_time = batch_end_time - batch_start_time
                    estimated_remaining = batch_time * (num_batches - batch_idx - 1)
                    logger.info(f"Batch {batch_idx+1}/{num_batches} completed in {batch_time:.2f}s. " +
                               f"Est. remaining: {estimated_remaining:.2f}s")
    
    overall_end_time = time.time()
    logger.info(f"GPU search completed in {overall_end_time - overall_start_time:.2f} seconds")
    logger.info(f"Found {len(all_result_indices)} matches total")
    
    return all_result_indices

def get_primer_from_index(idx, base, length):
    """
    Convert a primer index back to primer values.
    
    Args:
        idx (int): Primer index
        base (int): The number base to use
        length (int): Length of primers
        
    Returns:
        list: Primer as a list of digits
    """
    primer = []
    value = idx
    for i in range(length):
        primer.insert(0, value % base)
        value //= base
    return primer

def expand_key_and_decrypt(primer, base, ciphertext, cipher_alphabet=string.ascii_uppercase):
    """
    Expand a primer into a key and decrypt the ciphertext.
    
    Args:
        primer (list): Primer to expand
        base (int): The number base to use
        ciphertext (str): The ciphertext to decrypt
        cipher_alphabet (str): Ciphertext alphabet
        
    Returns:
        str: Decrypted plaintext
    """
    # Expand the primer
    key = []
    for i in range(len(primer)):
        key.append(primer[i])
    
    while len(key) < len(ciphertext):
        key.append((key[-1] + key[-2]) % base)
    
    # Decrypt
    plaintext = []
    for i, c in enumerate(ciphertext):
        if c in cipher_alphabet:
            idx = cipher_alphabet.index(c)
            pt_idx = (idx - key[i]) % len(cipher_alphabet)
            plaintext.append(cipher_alphabet[pt_idx])
        else:
            plaintext.append(c)
    
    return ''.join(plaintext)

def main():
    """Main function to run the GPU accelerated Gromark solver."""
    if not PYCUDA_AVAILABLE:
        logger.error("PyCUDA is required for this script")
        logger.error("Please install with: pip install pycuda")
        return
    
    parser = argparse.ArgumentParser(description="GPU-accelerated Kryptos K4 Gromark Cipher Solver")
    parser.add_argument("--base", type=int, default=10, help="Number base to use (default: 10)")
    parser.add_argument("--length", type=int, default=5, help="Primer length (default: 5)")
    parser.add_argument("--batch-size", type=int, default=1000000, help="Maximum number of primers to test in a single batch (default: 1,000,000)")
    parser.add_argument("--keyword", type=str, default=None, help="Keyword for mixed alphabet")
    parser.add_argument("--specific-primers", type=str, help="Comma-separated list of specific primers to try")
    parser.add_argument("--bean-primers", action="store_true", help="Use primers identified in Bean's research paper")
    
    args = parser.parse_args()
    
    # Create cipher alphabet (mixed if keyword provided)
    cipher_alphabet = string.ascii_uppercase
    if args.keyword:
        seen = set()
        keyword = args.keyword.upper()
        unique_keyword = ''.join(c for c in keyword if c in cipher_alphabet and c not in seen and not seen.add(c))
        remaining = ''.join(c for c in cipher_alphabet if c not in unique_keyword)
        cipher_alphabet = unique_keyword + remaining
        logger.info(f"Using mixed alphabet: {cipher_alphabet}")
    
    # Prepare specific primers to test
    specific_primers = []
    
    # Add primers from Bean's paper if requested
    if args.bean_primers:
        key = (args.base, args.length)
        if key in PROMISING_PRIMERS:
            specific_primers.extend(PROMISING_PRIMERS[key])
            logger.info(f"Added {len(PROMISING_PRIMERS[key])} primers from Bean's paper")
    
    # Add user-specified primers if provided
    if args.specific_primers:
        for primer_str in args.specific_primers.split(','):
            try:
                primer = []
                for c in primer_str.strip():
                    if c.isdigit():
                        digit = int(c)
                    elif args.base > 10 and 'A' <= c.upper() <= 'Z':
                        digit = ord(c.upper()) - ord('A') + 10
                    else:
                        raise ValueError(f"Invalid character '{c}' for base {args.base}")
                    
                    if digit >= args.base:
                        raise ValueError(f"Digit {digit} is not valid for base {args.base}")
                    
                    primer.append(digit)
                
                if len(primer) != args.length:
                    logger.warning(f"Primer {primer_str} length doesn't match specified length {args.length}")
                specific_primers.append(primer)
            except ValueError as e:
                logger.error(f"Error parsing primer '{primer_str}': {e}")
    
    # Run the GPU brute force search
    result_indices = gpu_brute_force_gromark(
        base=args.base,
        primer_length=args.length,
        cipher_alphabet=cipher_alphabet,
        max_batch_size=args.batch_size,
        specific_primers=specific_primers
    )
    
    # Convert indices to primers
    logger.info("Processing results...")
    candidates = []
    for primer_or_none, idx in result_indices:
        if primer_or_none is not None:
            # This is a specific primer that matched
            primer = primer_or_none
        else:
            # This is an index from the brute force search
            primer = get_primer_from_index(idx, args.base, args.length)
        
        plaintext = expand_key_and_decrypt(primer, args.base, K4_CIPHERTEXT, cipher_alphabet)
        score = score_text(plaintext)
        candidates.append((primer, plaintext, score, 0.0))  # IC is set to 0.0 for compatibility
    
    # Sort by score
    candidates.sort(key=lambda x: x[2], reverse=True)
    
    # Print results
    print_results(candidates)

if __name__ == "__main__":
    main() 