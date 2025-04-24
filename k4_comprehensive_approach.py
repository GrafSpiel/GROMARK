#!/usr/bin/env python3
"""
K4 Comprehensive Approach - combines Gromark decryption with multiple post-processing
techniques including letter frequency analysis, Caesar shifts, substitution patterns,
and columnar transposition, with progress visualization.
"""
import string
import re
import time
import argparse
import multiprocessing
from collections import Counter
import os
import math
from datetime import datetime
import sys
from tqdm import tqdm
import logging  # Add logging import
import queue  # Import the standard queue module

# Constants
K4_CIPHERTEXT = "OBKRUOXOGHULBSOLIFBBWFLRVQQPRNGKSSOTWTQSJQSSEKZZWATJKLUDIAWINFBNYPVTTMZFPKWGDKZXTJCDIGKUHUAUEKCAR"
ENGLISH_FREQ = {
    'e': 0.1202, 't': 0.0910, 'a': 0.0812, 'o': 0.0768, 'i': 0.0731,
    'n': 0.0695, 's': 0.0628, 'r': 0.0602, 'h': 0.0592, 'd': 0.0432,
    'l': 0.0398, 'u': 0.0288, 'c': 0.0271, 'm': 0.0261, 'f': 0.0230,
    'y': 0.0211, 'w': 0.0209, 'g': 0.0203, 'p': 0.0182, 'b': 0.0149,
    'v': 0.0111, 'k': 0.0069, 'x': 0.0017, 'q': 0.0011, 'j': 0.0010, 'z': 0.0007
}

# Known plaintext segments from previous analysis
KNOWN_PLAINTEXT = {
    84: ('S', 'A'),
    85: ('O', 'T'),
    86: ('M', 'J'),
    87: ('E', 'K'),
    88: ('T', 'L'),
    89: ('H', 'U'),
    90: ('I', 'D'),
    91: ('N', 'I'),
    92: ('G', 'A'),
    93: ('W', 'W')
}

# Promising primers from Bean's research
BEAN_PRIMERS = [
    {"name": "Bean Primer 1", "primer": [2, 6, 7, 1, 7], "base": 10, "description": "Primary primer from Bean's paper"},
    {"name": "Bean Primer 2", "primer": [8, 4, 3, 9, 3], "base": 10, "description": "Secondary primer from Bean's paper"},
    {"name": "Bean Primer 3", "primer": [9, 8, 8, 0, 0], "base": 10, "description": "Tertiary primer from Bean's paper"}
]

def expand_key(primer, length, base=10):
    """
    Expand a primer into a key of the desired length using Fibonacci-style addition.
    """
    key = primer.copy()
    while len(key) < length:
        next_digit = (key[-1] + key[-2]) % base
        key.append(next_digit)
    return key

def create_keyword_alphabet(keyword):
    """
    Create a mixed alphabet starting with the given keyword.
    """
    # Remove duplicates while preserving order
    seen = set()
    keyword = ''.join([c for c in keyword.upper() if c not in seen and not seen.add(c)])
    
    # Standard alphabet without the letters from the keyword
    remaining = ''.join([c for c in string.ascii_uppercase if c not in keyword])
    
    # Keyword followed by remaining letters
    return keyword + remaining

def decrypt_gromark(ciphertext, key, pt_alphabet=string.ascii_uppercase, ct_alphabet=string.ascii_uppercase):
    """
    Decrypt using the Gromark cipher with the given key and alphabets.
    """
    plaintext = ""
    pt_len = len(pt_alphabet)
    ct_len = len(ct_alphabet)
    for i, char in enumerate(ciphertext):
        if char in ct_alphabet:
            # Find position in ciphertext alphabet
            pos = ct_alphabet.index(char)
            # Apply key digit (subtraction for decryption), modulo the *plaintext* alphabet length
            new_pos = (pos - key[i % len(key)]) % pt_len
            # Map to plaintext alphabet
            plaintext += pt_alphabet[new_pos]
        else:
            plaintext += char
    return plaintext

def calculate_ic(text):
    """
    Calculate the Index of Coincidence for the given text.
    """
    text = re.sub(r'[^A-Z]', '', text.upper())
    if len(text) <= 1:
        return 0
    
    counter = Counter(text)
    sum_freqs = sum(count * (count - 1) for count in counter.values())
    n = len(text)
    ic = sum_freqs / (n * (n - 1))
    return ic

def verify_known_plaintext(plaintext):
    """
    Check how many known plaintext characters match.
    Return percentage of matches.
    """
    matches = 0
    for pos, (pt_char, _) in KNOWN_PLAINTEXT.items():
        if pos < len(plaintext) and plaintext[pos] == pt_char:
            matches += 1
    
    match_percentage = (matches / len(KNOWN_PLAINTEXT)) * 100
    return match_percentage

def calculate_letter_frequency_fit(text):
    """
    Calculate how well the letter frequencies in the text match English.
    Returns a score (lower is better).
    """
    text = re.sub(r'[^A-Za-z]', '', text.upper())
    if len(text) == 0:
        return float('inf')
    
    # Calculate frequency distribution
    counter = Counter(text.lower())
    total = len(text)
    
    # Calculate chi-square statistic
    chi_square = 0
    for char in string.ascii_lowercase:
        observed = counter.get(char, 0) / total
        expected = ENGLISH_FREQ.get(char, 0.0001)
        chi_square += ((observed - expected) ** 2) / expected
    
    return chi_square

def apply_caesar_shift(text, shift):
    """
    Apply a Caesar shift to the text.
    """
    shifted = ""
    for char in text:
        if char in string.ascii_uppercase:
            shifted_index = (string.ascii_uppercase.index(char) + shift) % 26
            shifted += string.ascii_uppercase[shifted_index]
        else:
            shifted += char
    return shifted

def apply_substitution(text, pattern):
    """
    Apply a simple substitution based on the pattern.
    Pattern examples: +3 (increment by 3), -2 (decrement by 2), etc.
    """
    substituted = ""
    op = pattern[0]
    value = int(pattern[1:])
    
    for char in text:
        if char in string.ascii_uppercase:
            if op == '+':
                new_index = (string.ascii_uppercase.index(char) + value) % 26
            elif op == '-':
                new_index = (string.ascii_uppercase.index(char) - value) % 26
            else:
                new_index = string.ascii_uppercase.index(char)
            substituted += string.ascii_uppercase[new_index]
        else:
            substituted += char
            
    return substituted

def columnar_transposition(text, key_length):
    """
    Apply columnar transposition with the given key length.
    """
    if key_length <= 1 or len(text) <= key_length:
        return text
    
    # Calculate the number of rows needed
    rows = math.ceil(len(text) / key_length)
    
    # Pad the text if needed
    padded_text = text + 'X' * (rows * key_length - len(text))
    
    # Create the columns
    columns = []
    for i in range(key_length):
        column = ""
        for j in range(rows):
            index = j * key_length + i
            if index < len(padded_text):
                column += padded_text[index]
        columns.append(column)
    
    # Join columns to create the transposed text
    transposed = ''.join(columns)
    
    # Remove padding if added
    if len(transposed) > len(text):
        transposed = transposed[:len(text)]
    
    return transposed

def process_keyword_primer_combination(args):
    """
    Process a single keyword-primer combination.
    """
    keyword, primer_config, transformation_type, update_queue = args
    
    results = []
    keyword_upper = keyword.upper()
    base = primer_config["base"]
    primer = primer_config["primer"]
    
    # Create keyword-based alphabet
    pt_alphabet = string.ascii_uppercase  # Standard alphabet for plaintext
    ct_alphabet = create_keyword_alphabet(keyword_upper)
    
    # Expand the primer to match ciphertext length
    key = expand_key(primer, len(K4_CIPHERTEXT), base)
    
    # Decrypt using Gromark
    decrypted = decrypt_gromark(K4_CIPHERTEXT, key, pt_alphabet, ct_alphabet)
    
    # Evaluate base result
    ic = calculate_ic(decrypted)
    match_percentage = verify_known_plaintext(decrypted)
    freq_score = calculate_letter_frequency_fit(decrypted)
    
    base_result = {
        'keyword': keyword,
        'primer': primer_config["name"],
        'transformation': "Base",
        'transformation_param': "None",
        'decrypted': decrypted,
        'ic': ic,
        'match_percentage': match_percentage,
        'freq_score': freq_score
    }
    
    # Only proceed with transformations if there's some potential (IC or match percentage)
    if ic >= 0.035 or match_percentage > 0:
        results.append(base_result)
        
        # Apply transformations based on type
        if transformation_type in ["all", "freq"]:
            # Letter frequency analysis is already done
            freq_result = base_result.copy()
            freq_result['transformation'] = "Frequency"
            results.append(freq_result)
        
        if transformation_type in ["all", "caesar"]:
            # Test Caesar shifts
            for shift in range(1, 26):
                shifted = apply_caesar_shift(decrypted, shift)
                shifted_ic = calculate_ic(shifted)
                shifted_match = verify_known_plaintext(shifted)
                shifted_freq = calculate_letter_frequency_fit(shifted)
                
                if shifted_ic >= 0.035 or shifted_match > 0:
                    results.append({
                        'keyword': keyword,
                        'primer': primer_config["name"],
                        'transformation': "Caesar",
                        'transformation_param': str(shift),
                        'decrypted': shifted,
                        'ic': shifted_ic,
                        'match_percentage': shifted_match,
                        'freq_score': shifted_freq
                    })
        
        if transformation_type in ["all", "substitution"]:
            # Test simple substitution patterns
            for op in ['+', '-']:
                for val in [1, 2, 3]:
                    pattern = f"{op}{val}"
                    substituted = apply_substitution(decrypted, pattern)
                    sub_ic = calculate_ic(substituted)
                    sub_match = verify_known_plaintext(substituted)
                    sub_freq = calculate_letter_frequency_fit(substituted)
                    
                    if sub_ic >= 0.035 or sub_match > 0:
                        results.append({
                            'keyword': keyword,
                            'primer': primer_config["name"],
                            'transformation': "Substitution",
                            'transformation_param': pattern,
                            'decrypted': substituted,
                            'ic': sub_ic,
                            'match_percentage': sub_match,
                            'freq_score': sub_freq
                        })
        
        if transformation_type in ["all", "transposition"]:
            # Test columnar transposition
            for key_length in range(2, 11):
                transposed = columnar_transposition(decrypted, key_length)
                trans_ic = calculate_ic(transposed)
                trans_match = verify_known_plaintext(transposed)
                trans_freq = calculate_letter_frequency_fit(transposed)
                
                if trans_ic >= 0.035 or trans_match > 0:
                    results.append({
                        'keyword': keyword,
                        'primer': primer_config["name"],
                        'transformation': "Transposition",
                        'transformation_param': str(key_length),
                        'decrypted': transposed,
                        'ic': trans_ic,
                        'match_percentage': trans_match,
                        'freq_score': trans_freq
                    })
    
    # Update progress
    if update_queue:
        update_queue.put(1)
    
    return results

def load_wordlist(filename, min_length=3, max_length=20):
    """
    Load wordlist from a file, filtering words by length.
    """
    words = []
    try:
        with open(filename, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                word = line.strip()
                if min_length <= len(word) <= max_length:
                    words.append(word)
    except FileNotFoundError:
        print(f"Error: Wordlist file '{filename}' not found.")
        sys.exit(1)
    
    return words

def save_results(results, transformation_type, output_dir="output"):
    """
    Save results to files organized by transformation type.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Organize results by transformation
    transformation_files = {
        "all": f"{output_dir}/all_results_{timestamp}.txt",
        "freq": f"{output_dir}/frequency_results_{timestamp}.txt",
        "caesar": f"{output_dir}/caesar_results_{timestamp}.txt",
        "substitution": f"{output_dir}/substitution_results_{timestamp}.txt",
        "transposition": f"{output_dir}/transposition_results_{timestamp}.txt"
    }
    
    # Sort results by match percentage (primary) and IC (secondary)
    sorted_results = sorted(results, key=lambda x: (-x['match_percentage'], -x['ic']))
    
    # Write all results to a combined file
    with open(transformation_files["all"], 'w', encoding='utf-8') as f:
        f.write(f"Total results: {len(sorted_results)}\n")
        f.write("=" * 80 + "\n\n")
        
        for i, result in enumerate(sorted_results):
            f.write(f"Rank {i+1}\n")
            f.write(f"Keyword: {result['keyword']}\n")
            f.write(f"Primer: {result['primer']}\n")
            f.write(f"Transformation: {result['transformation']}\n")
            f.write(f"Transformation Parameter: {result['transformation_param']}\n")
            f.write(f"Match: {result['match_percentage']:.2f}%\n")
            f.write(f"IC: {result['ic']:.6f}\n")
            f.write(f"Frequency Score: {result['freq_score']:.6f}\n")
            f.write(f"Decrypted: {result['decrypted']}\n\n")
            f.write("-" * 80 + "\n\n")
    
    # Write results to individual transformation files
    for transformation in ["freq", "caesar", "substitution", "transposition"]:
        filtered_results = [r for r in sorted_results if r['transformation'].lower() == transformation or
                           (transformation == "freq" and r['transformation'] == "Base")]
        
        if filtered_results:
            with open(transformation_files[transformation], 'w', encoding='utf-8') as f:
                f.write(f"Total {transformation} results: {len(filtered_results)}\n")
                f.write("=" * 80 + "\n\n")
                
                for i, result in enumerate(filtered_results):
                    f.write(f"Rank {i+1}\n")
                    f.write(f"Keyword: {result['keyword']}\n")
                    f.write(f"Primer: {result['primer']}\n")
                    f.write(f"Transformation: {result['transformation']}\n")
                    f.write(f"Transformation Parameter: {result['transformation_param']}\n")
                    f.write(f"Match: {result['match_percentage']:.2f}%\n")
                    f.write(f"IC: {result['ic']:.6f}\n")
                    f.write(f"Frequency Score: {result['freq_score']:.6f}\n")
                    f.write(f"Decrypted: {result['decrypted']}\n\n")
                    f.write("-" * 80 + "\n\n")
    
    return transformation_files["all"]

def save_checkpoint(results, filename="search_checkpoint.txt"):
    """
    Save a checkpoint of the current top results.
    """
    # Sort results by match percentage (primary) and IC (secondary)
    sorted_results = sorted(results, key=lambda x: (-x['match_percentage'], -x['ic']))
    
    # Take top 100 or all if fewer
    top_results = sorted_results[:min(100, len(sorted_results))]
    
    with open(filename, 'w', encoding='utf-8') as f:
        f.write(f"Checkpoint saved at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Top {len(top_results)} results:\n")
        f.write("=" * 80 + "\n\n")
        
        for i, result in enumerate(top_results):
            f.write(f"Rank {i+1}\n")
            f.write(f"Keyword: {result['keyword']}\n")
            f.write(f"Primer: {result['primer']}\n")
            f.write(f"Transformation: {result['transformation']}\n")
            f.write(f"Transformation Parameter: {result['transformation_param']}\n")
            f.write(f"Match: {result['match_percentage']:.2f}%\n")
            f.write(f"IC: {result['ic']:.6f}\n")
            f.write(f"Frequency Score: {result['freq_score']:.6f}\n")
            f.write(f"Decrypted: {result['decrypted']}\n\n")
            f.write("-" * 40 + "\n\n")

def update_progress(queue_obj, total, checkpoint_interval=1000):
    """
    Update the progress bar based on queue messages.
    """
    progress_bar = tqdm(total=total, desc="Processing combinations")
    count = 0
    
    while count < total:
        try:
            msg = queue_obj.get(timeout=1)  # Use the passed queue object
            if msg == "DONE":
                break
            count += msg
            progress_bar.update(msg)
        except queue.Empty: # Use queue.Empty from the imported module
            # If the queue is empty after timeout, continue loop, maybe main process is just slow
            continue
        except (EOFError, BrokenPipeError) as e:
            logging.warning(f"IPC Error in progress tracking (Pipe closed?): {e}")
            break # Exit loop if pipe is broken
        except Exception as e:
            # Log other unexpected errors
            logging.exception(f"Unexpected error in progress tracking: {e}")
            # Depending on the error, you might want to break or continue
            break # Exit loop on other errors for safety
    
    progress_bar.close()

def main():
    parser = argparse.ArgumentParser(description='Comprehensive K4 Solver with multiple transformation strategies')
    parser.add_argument('--wordlist', type=str, required=True, help='Path to wordlist file')
    parser.add_argument('--min-length', type=int, default=3, help='Minimum word length to consider')
    parser.add_argument('--max-length', type=int, default=15, help='Maximum word length to consider')
    parser.add_argument('--transformation', type=str, default='all', choices=['all', 'freq', 'caesar', 'substitution', 'transposition'],
                        help='Type of transformation to apply')
    parser.add_argument('--processes', type=int, default=None, 
                        help='Number of processes to use (default: number of CPU cores)')
    parser.add_argument('--batch-size', type=int, default=1000, 
                        help='Number of words to process in each batch')
    parser.add_argument('--checkpoint-interval', type=int, default=5000, 
                        help='Save checkpoint after processing this many combinations')
    parser.add_argument('--output-dir', type=str, default='output', 
                        help='Directory to save results')
    parser.add_argument('--specific-keywords', type=str, 
                        help='Comma-separated list of specific keywords to test')
    
    args = parser.parse_args()
    
    # Start timing
    start_time = time.time()
    
    # Set up multiprocessing
    processes = args.processes if args.processes else multiprocessing.cpu_count()
    print(f"Using {processes} processes")
    
    # Load or generate wordlist
    if args.specific_keywords:
        words = [word.strip() for word in args.specific_keywords.split(',')]
        print(f"Testing {len(words)} specific keywords")
    else:
        words = load_wordlist(args.wordlist, args.min_length, args.max_length)
        print(f"Loaded {len(words)} words from {args.wordlist}")
    
    # Create tasks for all combinations
    all_tasks = []
    for word in words:
        for primer_config in BEAN_PRIMERS:
            all_tasks.append((word, primer_config, args.transformation, None))
    
    # Set up progress tracking
    total_combinations = len(all_tasks)
    print(f"Total combinations to test: {total_combinations}")
    
    # Process in batches with progress tracking
    all_results = []
    
    with multiprocessing.Manager() as manager:
        update_queue = manager.Queue()
        
        # Start the progress updater process
        progress_process = multiprocessing.Process(
            target=update_progress, 
            args=(update_queue, total_combinations, args.checkpoint_interval)
        )
        progress_process.start()
        
        # Process in batches
        batches = [all_tasks[i:i + args.batch_size] for i in range(0, len(all_tasks), args.batch_size)]
        
        for batch_idx, batch in enumerate(batches):
            print(f"\nProcessing batch {batch_idx + 1}/{len(batches)}")
            
            # Update tasks with queue
            batch_with_queue = [(word, primer, trans_type, update_queue) for word, primer, trans_type, _ in batch]
            
            # Process batch
            with multiprocessing.Pool(processes=processes) as pool:
                batch_results = pool.map(process_keyword_primer_combination, batch_with_queue)
            
            # Flatten results
            flat_results = [item for sublist in batch_results for item in sublist]
            all_results.extend(flat_results)
            
            # Save checkpoint if needed
            if (batch_idx + 1) * args.batch_size % args.checkpoint_interval < args.batch_size:
                save_checkpoint(all_results, f"{args.output_dir}/checkpoint_batch_{batch_idx + 1}.txt")
                print(f"Checkpoint saved after batch {batch_idx + 1}")
        
        # Signal progress process to finish
        try:
            update_queue.put("DONE")
        except Exception as e:
            logging.error(f"Error sending DONE signal to progress queue: {e}")

        progress_process.join(timeout=10) # Add a timeout for joining
        if progress_process.is_alive():
            logging.warning("Progress process did not terminate gracefully. Terminating.")
            progress_process.terminate()
            progress_process.join()
    
    # Save final results
    output_file = save_results(all_results, args.transformation, args.output_dir)
    print(f"Results saved to {output_file}")
    
    # Print execution time
    end_time = time.time()
    execution_time = end_time - start_time
    print(f"Execution time: {execution_time:.2f} seconds")

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO) # Configure basic logging
    main() 