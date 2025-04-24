#!/usr/bin/env python3
"""
Gromark Cipher Solver for Kryptos K4
This script implements a brute-force approach to solving Kryptos K4 using the Gromark cipher.

Based on findings from "Cryptodiagnosis of Kryptos K4" by Richard Bean, which identified
specific promising primers for different bases.
"""

import string
import time
from collections import Counter
import argparse
from multiprocessing import Pool, cpu_count
import os

# K4 ciphertext (97 letters)
K4_CIPHERTEXT = "OBKRUOXOGHULBSOLIFBBWFLRVQQPRNGKSSOTWTQSJQSSEKZZWATJKLUDIAWINFBNYPVTTMZFPKWGDKZXTJCDIGKUHUAUEKCAR"

# Known plaintext cribs
CRIB1 = ("EASTNORTHEAST", (20, 33))  # positions 21-34 (0-indexed)
CRIB2 = ("BERLINCLOCK", (48, 58))    # positions 49-59 (0-indexed)

# Default alphabets
DEFAULT_PLAIN_ALPHABET = string.ascii_uppercase
DEFAULT_CIPHER_ALPHABET = string.ascii_uppercase

# English language quadgram statistics for scoring
ENGLISH_QUADGRAMS = {}

# Promising primers identified in Bean's paper
PROMISING_PRIMERS = {
    # Format: (base, length): [list of primers]
    (10, 5): [
        # These are the 39 primers for base 10, length 5 from Bean's paper
        [2, 6, 7, 1, 7], [8, 4, 3, 9, 3],  # Equivalent primers using only 9 digits
        [9, 8, 8, 0, 0],  # Best IC value of 0.0625
        # Adding remaining unique primers from the paper
        [3, 2, 6, 8, 8], [4, 5, 9, 5, 1], [5, 4, 0, 5, 5], [6, 8, 5, 7, 8],
        [6, 9, 3, 2, 3], [3, 4, 5, 7, 0], [4, 1, 3, 5, 2], [5, 1, 6, 1, 4],
        [5, 3, 6, 9, 0], [5, 8, 6, 3, 6], [7, 0, 9, 7, 4], [7, 3, 5, 4, 8],
        [7, 5, 0, 5, 7], [8, 5, 2, 6, 9], [8, 5, 4, 2, 7], [9, 2, 1, 8, 5],
        [1, 3, 3, 6, 9], [1, 5, 0, 4, 3], [1, 6, 3, 9, 2], [1, 7, 4, 5, 4],
        [2, 4, 6, 4, 0], [3, 5, 9, 2, 9], [5, 9, 0, 2, 0], [5, 9, 1, 0, 0],
        [5, 9, 2, 0, 0], [5, 9, 3, 0, 0], [5, 9, 4, 0, 0], [5, 9, 5, 0, 0],
        [5, 9, 6, 0, 0], [5, 9, 7, 0, 0], [5, 9, 8, 0, 0], [5, 9, 9, 0, 0],
        [6, 0, 3, 5, 9], [6, 3, 0, 6, 2], [7, 4, 5, 3, 5], [7, 5, 7, 0, 3]
    ],
    (10, 4): [
        [3, 3, 0, 1], [6, 7, 4, 0], [9, 9, 0, 3]  # Base 10, length 4 primers
    ],
    (8, 5): [
        [0, 0, 3, 5, 1], [0, 0, 5, 3, 7],  # Base 8, length 5 primers with period 84
        [0, 0, 1, 3, 7], [0, 0, 7, 7, 3]   # Other two base 8, length 5 primers mentioned
    ],
    (5, 5): [
        # Base 5 primers suggested by the Berlin Clock reference
        [0, 0, 1, 2, 3], [0, 1, 2, 3, 4], [1, 2, 3, 4, 0]  # Suggested primers for base 5
    ],
    (12, 5): [
        # Base 12 primers suggested by the Berlin Clock reference
        [0, 1, 2, 3, 4], [1, 2, 3, 4, 5], [11, 0, 1, 2, 3]  # Suggested primers for base 12 (11=B)
    ]
}

# Deduplicate primers for each key in the PROMISING_PRIMERS dictionary
for key in PROMISING_PRIMERS:
    # Convert each list of primers to a set of tuples (for hashability), then back to list of lists
    unique_primers = {tuple(primer) for primer in PROMISING_PRIMERS[key]}
    PROMISING_PRIMERS[key] = [list(primer) for primer in unique_primers]
    print(f"Base {key[0]}, length {key[1]}: {len(PROMISING_PRIMERS[key])} unique primers")

def expand_gromark_primer(primer, base, length):
    """
    Expand a Gromark primer to generate a full period sequence.
    
    Args:
        primer (list): The initial primer values
        base (int): The number base (e.g., 10 for decimal)
        length (int): The desired length of the output key sequence
        
    Returns:
        Generator yielding digits of the full period sequence
    """
    # Create a working copy of the primer
    state = primer.copy()
    
    # Start yielding the initial primer values
    for digit in state:
        yield digit
    
    # Calculate how many more digits we need to generate
    remaining = length - len(primer)
    
    # Generate remaining digits
    for _ in range(remaining):
        # Calculate the new digit (sum of the two previous digits)
        new_digit = (state[-1] + state[-2]) % base
        
        # Yield the new digit
        yield new_digit
        
        # Shift the state and add the new digit
        state = state[1:] + [new_digit]

def fibonacci_expansion(primer, base, length):
    """
    Expand a primer using Fibonacci-style addition (a variant of Gromark).
    
    Args:
        primer (list): Initial primer digits
        base (int): The number base to use
        length (int): Desired length of the key sequence
        
    Returns:
        list: Expanded key sequence
    """
    key = list(primer)
    while len(key) < length:
        next_digit = 0
        # Sum the last n digits (where n is the primer length)
        for i in range(1, len(primer) + 1):
            if i <= len(key):
                next_digit = (next_digit + key[-i]) % base
        key.append(next_digit)
    return key[:length]

def expand_gromark_primer_full(primer, base, length):
    """
    Expand a Gromark primer to generate its full sequence as a list.
    This is a convenience function that collects all digits from the generator.
    
    Args:
        primer (list): The initial primer values
        base (int): The number base (e.g., 10 for decimal)
        length (int): The length of the primer
        
    Returns:
        list: The full sequence of digits
    """
    return list(expand_gromark_primer(primer, base, length))

def get_period(primer, base, length):
    """
    Calculate the period of a Gromark sequence.
    
    Args:
        primer (list): The initial primer values
        base (int): The number base (e.g., 10 for decimal)
        length (int): The length of the primer
        
    Returns:
        int: The period of the sequence
    """
    # Create a working copy of the primer
    state = primer.copy()
    
    # Track the states we've seen before to detect cycles
    seen_states = {tuple(state): 0}
    
    # Calculate the maximum period (theoretical)
    max_period = base ** length - 1
    
    # Start counting
    for i in range(max_period):
        # Calculate the new digit
        new_digit = sum(state) % base
        
        # Shift the state and add the new digit
        state = state[1:] + [new_digit]
        
        # Check if we've seen this state before
        state_tuple = tuple(state)
        if state_tuple in seen_states:
            # We've found a cycle
            return i + 1 - seen_states[state_tuple]
        
        # Store the current state
        seen_states[state_tuple] = i + 1
    
    # If we get here, the period is the maximum period
    return max_period

def decrypt_gromark(ct, key, plain_alphabet=DEFAULT_PLAIN_ALPHABET, cipher_alphabet=DEFAULT_CIPHER_ALPHABET):
    """
    Decrypt a ciphertext using the Gromark cipher.
    
    Args:
        ct (str): Ciphertext to decrypt
        key (list): Key sequence of numbers
        plain_alphabet (str): Plaintext alphabet
        cipher_alphabet (str): Ciphertext alphabet
        
    Returns:
        str: Decrypted plaintext
    """
    # Create inverse mapping for cipher alphabet
    invC = {c: i for i, c in enumerate(cipher_alphabet)}
    pt = []
    
    for i, c in enumerate(ct):
        cnum = invC[c]
        pnum = (cnum - key[i]) % len(plain_alphabet)
        pt.append(plain_alphabet[pnum])
    
    return ''.join(pt)

def check_cribs(plaintext, crib1=CRIB1, crib2=CRIB2, min_match_percentage=0.7):
    """
    Check if decrypted plaintext contains the expected cribs at the right positions.
    More lenient version that allows for partial matches.
    
    Args:
        plaintext (str): Decrypted plaintext
        crib1 (tuple): (text, (start, end)) for first crib
        crib2 (tuple): (text, (start, end)) for second crib
        min_match_percentage (float): Minimum percentage of characters that must match (0.0-1.0)
        
    Returns:
        bool: True if both cribs match sufficiently, False otherwise
    """
    crib1_text, (crib1_start, crib1_end) = crib1
    crib2_text, (crib2_start, crib2_end) = crib2
    
    # Check if the plaintext is long enough
    if len(plaintext) <= max(crib1_end, crib2_end):
        return False
    
    # Calculate matches for first crib
    crib1_matches = 0
    for i in range(len(crib1_text)):
        if crib1_start + i < len(plaintext) and plaintext[crib1_start + i] == crib1_text[i]:
            crib1_matches += 1
    
    # Calculate matches for second crib
    crib2_matches = 0
    for i in range(len(crib2_text)):
        if crib2_start + i < len(plaintext) and plaintext[crib2_start + i] == crib2_text[i]:
            crib2_matches += 1
    
    # Calculate match percentages
    crib1_percentage = crib1_matches / len(crib1_text) if len(crib1_text) > 0 else 0
    crib2_percentage = crib2_matches / len(crib2_text) if len(crib2_text) > 0 else 0
    
    # Both cribs must meet the minimum match percentage
    return crib1_percentage >= min_match_percentage and crib2_percentage >= min_match_percentage

def load_english_quadgrams(filename=None):
    """
    Load English quadgram statistics for scoring text.
    If filename is not provided, use a simple approximation based on common letter frequencies.
    
    Args:
        filename (str, optional): Path to quadgram statistics file
    """
    global ENGLISH_QUADGRAMS
    
    if filename:
        # Load from file
        try:
            with open(filename, 'r') as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        quadgram, count = parts[0], int(parts[1])
                        ENGLISH_QUADGRAMS[quadgram] = count
        except FileNotFoundError:
            print(f"Warning: Quadgram file {filename} not found. Using basic scoring.")
            ENGLISH_QUADGRAMS = {}
    
    # If no file or loading failed, use basic letter frequencies
    if not ENGLISH_QUADGRAMS:
        # Simple approximation using common letter frequencies
        letter_freq = {'E': 12.7, 'T': 9.1, 'A': 8.2, 'O': 7.5, 'I': 7.0, 'N': 6.7, 'S': 6.3, 
                      'H': 6.1, 'R': 6.0, 'D': 4.3, 'L': 4.0, 'U': 2.8, 'C': 2.8, 'M': 2.4,
                      'F': 2.2, 'W': 2.0, 'Y': 2.0, 'G': 1.9, 'P': 1.9, 'B': 1.5, 'V': 1.0,
                      'K': 0.8, 'J': 0.2, 'X': 0.2, 'Q': 0.1, 'Z': 0.1}
        
        # Generate basic quadgram frequencies based on letter frequencies
        for a in string.ascii_uppercase:
            for b in string.ascii_uppercase:
                for c in string.ascii_uppercase:
                    for d in string.ascii_uppercase:
                        quadgram = a + b + c + d
                        score = (letter_freq.get(a, 0.1) * letter_freq.get(b, 0.1) * 
                                letter_freq.get(c, 0.1) * letter_freq.get(d, 0.1))
                        ENGLISH_QUADGRAMS[quadgram] = int(score * 100)

def score_text(text):
    """
    Score a text for English-like qualities using quadgram statistics.
    
    Args:
        text (str): Text to score
        
    Returns:
        float: Score (higher is more English-like)
    """
    if not ENGLISH_QUADGRAMS:
        load_english_quadgrams()
    
    text = text.upper()
    score = 0
    
    for i in range(len(text) - 3):
        quadgram = text[i:i+4]
        if all(c in string.ascii_uppercase for c in quadgram):
            score += ENGLISH_QUADGRAMS.get(quadgram, 0)
    
    return score

def calculate_ic(text):
    """
    Calculate the Index of Coincidence for a text.
    
    Args:
        text (str): Text to analyze
        
    Returns:
        float: Index of Coincidence
    """
    # Count each letter
    counts = Counter(text.upper())
    
    # Calculate IC
    n = len(text)
    if n <= 1:
        return 0
    
    sum_freqs = 0
    for count in counts.values():
        sum_freqs += count * (count - 1)
    
    return sum_freqs / (n * (n - 1))

def test_primer(task):
    """
    Test a specific primer against the ciphertext.
    
    Args:
        task (tuple): (primer, base, ciphertext, plain_alphabet, cipher_alphabet, expansion_func, min_match_percentage)
        
    Returns:
        tuple: (primer, plaintext, score, ic) if matches cribs, None otherwise
    """
    primer, base, ciphertext, plain_alphabet, cipher_alphabet, expansion_func, min_match_percentage = task
    
    # Generate the key
    key = list(expansion_func(primer, base, len(ciphertext)))
    
    # Decrypt
    plaintext = decrypt_gromark(ciphertext, key, plain_alphabet, cipher_alphabet)
    
    # Check if it matches the cribs
    if check_cribs(plaintext, min_match_percentage=min_match_percentage):
        # Calculate index of coincidence
        ic = calculate_ic(plaintext)
        
        # Score based on IC and the presence of expected words
        score = calculate_ic(plaintext)
        
        # Return the candidate
        return (primer, plaintext, score, ic)
    
    return None

def generate_primers(base, length):
    """
    Generate all possible primers of a given length in a specified base.
    
    Args:
        base (int): The number base to use
        length (int): Length of primers to generate
        
    Yields:
        list: Each possible primer
    """
    # Generate primers from [0, 0, ..., 0] to [base-1, base-1, ..., base-1]
    if length == 0:
        yield []
        return
    
    for prefix in generate_primers(base, length - 1):
        for digit in range(base):
            yield prefix + [digit]

def check_vertical_bigrams(text, width=21):
    """
    Check the number of repeated vertical bigrams when written at specified width.
    
    Args:
        text (str): Text to analyze
        width (int): Width to use for formatting
        
    Returns:
        int: Number of repeated vertical bigrams
    """
    # Split the text into rows of specified width
    rows = [text[i:i+width] for i in range(0, len(text), width)]
    
    # Build vertical bigrams
    bigrams = []
    for col in range(min(width, len(rows[0]))):
        for row in range(len(rows) - 1):
            if col < len(rows[row]) and col < len(rows[row+1]):
                bigram = rows[row][col] + rows[row+1][col]
                bigrams.append(bigram)
    
    # Count repeated bigrams
    bigram_counts = Counter(bigrams)
    repeated = 0
    for bigram, count in bigram_counts.items():
        if count > 1:
            repeated += count - 1
    
    return repeated

def brute_force_gromark(base=10, primer_length=5, ciphertext=K4_CIPHERTEXT, 
                        plain_alphabet=DEFAULT_PLAIN_ALPHABET, 
                        cipher_alphabet=DEFAULT_CIPHER_ALPHABET,
                        processes=None, specific_primers=None,
                        expansion_func=expand_gromark_primer,
                        chunk_size=10000, min_match_percentage=0.7):
    """
    Brute force Gromark primers and find those that match the cribs.
    
    Args:
        base (int): The number base to use
        primer_length (int): Length of primers to test
        ciphertext (str): The ciphertext to decrypt
        plain_alphabet (str): Plaintext alphabet
        cipher_alphabet (str): Ciphertext alphabet
        processes (int, optional): Number of processes to use (None = use all available)
        specific_primers (list, optional): List of specific primers to try first
        expansion_func (function): Function to use for key expansion
        chunk_size (int): Size of chunks for batched processing
        min_match_percentage (float): Minimum percentage of characters that must match for cribs
        
    Returns:
        list: List of (primer, plaintext, score, ic) tuples for matching primers
    """
    start_time = time.time()
    candidates = []
    
    # First, try any specific primers if provided
    if specific_primers:
        print(f"Testing {len(specific_primers)} specific primers...")
        
        # Test each specific primer
        for primer in specific_primers:
            # Create task for multiprocessing
            task = (primer, base, ciphertext, plain_alphabet, cipher_alphabet, expansion_func, min_match_percentage)
            
            # Test the primer
            result = test_primer(task)
            
            if result:
                primer, plaintext, score, ic = result
                candidates.append(result)
                
                print(f"CANDIDATE FOUND: {primer}, Score: {score}, IC: {ic:.6f}")
                print(f"Plaintext: {plaintext}")
                
                # Check vertical bigrams at width 21
                vb_count = check_vertical_bigrams(plaintext, 21)
                print(f"Vertical bigrams at width 21: {vb_count}")
    
    # Count total primers for brute force search
    if specific_primers is None or len(candidates) == 0:
        total_primers = base ** primer_length
        print(f"Testing {total_primers} primers (base {base}, length {primer_length})...")
        
        # Use multiprocessing to speed up the brute force search
        if processes is None:
            processes = cpu_count()
        
        with Pool(processes=processes) as pool:
            # Process in chunks to avoid memory issues with large primer spaces
            total_chunks = (total_primers + chunk_size - 1) // chunk_size
            
            for chunk_idx in range(total_chunks):
                # Generate a chunk of primers
                chunk_start = chunk_idx * chunk_size
                chunk_end = min(chunk_start + chunk_size, total_primers)
                chunk_size_actual = chunk_end - chunk_start
                
                # Create generator for this chunk's primers
                def primers_for_chunk(chunk_start, chunk_size_actual):
                    # Generate digits for each position
                    for idx in range(chunk_start, chunk_start + chunk_size_actual):
                        primer = []
                        value = idx
                        for i in range(primer_length):
                            primer.insert(0, value % base)
                            value //= base
                        yield primer
                
                # Create tasks for this chunk
                tasks = [(primer, base, ciphertext, plain_alphabet, cipher_alphabet, expansion_func, min_match_percentage) 
                        for primer in primers_for_chunk(chunk_start, chunk_size_actual)]
                
                # Process this chunk
                print(f"Processing chunk {chunk_idx+1}/{total_chunks} ({chunk_start}-{chunk_end-1})...")
                chunk_start_time = time.time()
                
                for i, result in enumerate(pool.imap_unordered(test_primer, tasks)):
                    if result:
                        candidates.append(result)
                        
                    # Print progress every 1,000 primers or when a candidate is found
                    if i % 1000 == 0 or result:
                        elapsed = time.time() - start_time
                        global_idx = chunk_start + i
                        percent = global_idx / total_primers * 100
                        primers_per_sec = global_idx / elapsed if elapsed > 0 else 0
                        print(f"Progress: {global_idx}/{total_primers} ({percent:.2f}%) - {primers_per_sec:.2f} primers/sec")
                        
                        if result:
                            primer, plaintext, score, ic = result
                            print(f"CANDIDATE FOUND: {primer}, Score: {score}, IC: {ic:.6f}")
                            print(f"Plaintext: {plaintext}")
                            
                            # Check vertical bigrams at width 21
                            vb_count = check_vertical_bigrams(plaintext, 21)
                            print(f"Vertical bigrams at width 21: {vb_count}")
                
                chunk_end_time = time.time()
                chunk_time = chunk_end_time - chunk_start_time
                print(f"Chunk {chunk_idx+1}/{total_chunks} completed in {chunk_time:.2f}s")
    
    # Sort candidates by score (highest first)
    candidates.sort(key=lambda x: x[2], reverse=True)
    
    end_time = time.time()
    print(f"Brute force completed in {end_time - start_time:.2f} seconds")
    print(f"Found {len(candidates)} candidates")
    
    return candidates

def print_results(candidates, top_n=10):
    """
    Print the top candidates with their decrypted plaintext.
    
    Args:
        candidates (list): List of (primer, plaintext, score, ic) tuples
        top_n (int): Number of top candidates to print
    """
    if not candidates:
        print("No matching primers found.")
        return
    
    print("\n===== TOP CANDIDATES =====")
    for i, (primer, plaintext, score, ic) in enumerate(candidates[:top_n], 1):
        print(f"\n#{i}: Primer: {primer}, Score: {score}, IC: {ic:.6f}")
        print(f"Plaintext: {plaintext}")
        
        # Check vertical bigrams at width 21
        vb_count = check_vertical_bigrams(plaintext, 21)
        print(f"Vertical bigrams at width 21: {vb_count}")
        
        # Print the plaintext arranged at width 21
        print("Plaintext at width 21:")
        for j in range(0, len(plaintext), 21):
            print(plaintext[j:j+21])
        
        # Highlight the crib sections
        highlighted = list(plaintext)
        crib1_text, (crib1_start, crib1_end) = CRIB1
        crib2_text, (crib2_start, crib2_end) = CRIB2
        
        for j in range(crib1_start, crib1_end + 1):
            highlighted[j] = f"[{highlighted[j]}]"
        for j in range(crib2_start, crib2_end + 1):
            highlighted[j] = f"[{highlighted[j]}]"
        
        print(f"With cribs: {''.join(highlighted)}")

def primer_list_from_string(primer_str, base=10):
    """
    Convert a string representation of a primer to a list of digits.
    
    Args:
        primer_str (str): String representation of the primer (e.g., "26717")
        base (int): The number base to use
        
    Returns:
        list: Primer as a list of digits
    """
    # Convert each character to an integer, handling bases > 10
    primer = []
    for c in primer_str:
        if c.isdigit():
            digit = int(c)
        elif base > 10 and 'A' <= c.upper() <= 'Z':
            # For bases > 10, use A=10, B=11, etc.
            digit = ord(c.upper()) - ord('A') + 10
        else:
            raise ValueError(f"Invalid character '{c}' for base {base}")
        
        if digit >= base:
            raise ValueError(f"Digit {digit} is not valid for base {base}")
        
        primer.append(digit)
    
    return primer

def cpu_search_primers(ciphertext, primer_list, base, length, target_plaintext=None, target_substrings=None, 
                      min_alpha_ratio=0.7, min_ic=0.06, max_primers=None, expansion_fn=None,
                      plain_alphabet=DEFAULT_PLAIN_ALPHABET, cipher_alphabet=DEFAULT_CIPHER_ALPHABET):
    """
    Search through primers using CPU-based approach.
    
    Args:
        ciphertext (str): The Gromark ciphertext to decrypt
        primer_list (list): List of primers to test
        base (int): The base of the Gromark cipher (usually 10)
        length (int): The length of each primer
        target_plaintext (str, optional): Plaintext to look for
        target_substrings (list, optional): List of substrings to look for in decryption
        min_alpha_ratio (float, optional): Minimum alphabet ratio to consider
        min_ic (float, optional): Minimum index of coincidence to consider
        max_primers (int, optional): Maximum number of primers to test
        expansion_fn (function, optional): Function to use for expanding primers
        plain_alphabet (str, optional): Plain alphabet to use
        cipher_alphabet (str, optional): Cipher alphabet to use
        
    Returns:
        list: List of (primer, plaintext, score, period) tuples for promising candidates
    """
    # Default to standard Gromark expansion if none provided
    if expansion_fn is None:
        expansion_fn = expand_gromark_primer
    
    promising_candidates = []
    total_primers = len(primer_list)
    processed = 0
    
    # Limit the number of primers if specified
    if max_primers and max_primers < total_primers:
        primer_list = primer_list[:max_primers]
        total_primers = max_primers

    print(f"CPU searching through {total_primers} primers")
    
    # Process each primer
    for primer in primer_list:
        try:
            # Generate the key (convert generator to a list)
            key = list(expansion_fn(primer, base, len(ciphertext)))
            
            # Decrypt using the key
            plaintext = decrypt_gromark(ciphertext, key, plain_alphabet, cipher_alphabet)
            
            # Check if we have a promising candidate
            score = score_text(plaintext)
            ic = calculate_ic(plaintext)
            
            # Apply filtering criteria
            if (target_plaintext is None or target_plaintext in plaintext) and \
               (target_substrings is None or all(sub in plaintext for sub in target_substrings)) and \
               (min_alpha_ratio is None or sum(c.isalpha() for c in plaintext) / len(plaintext) >= min_alpha_ratio) and \
               (min_ic is None or ic >= min_ic):
                # Calculate the period for promising candidates
                period = get_period(primer, base, length)
                promising_candidates.append((primer, plaintext, score, period))
        except Exception as e:
            print(f"Error processing primer {primer}: {e}")
        
        processed += 1
        if processed % 100 == 0:
            print(f"Processed {processed}/{total_primers} primers")
    
    # Sort by score in descending order
    promising_candidates.sort(key=lambda x: x[2], reverse=True)
    
    return promising_candidates

def main():
    """Main function to run the Gromark solver."""
    parser = argparse.ArgumentParser(description="Kryptos K4 Gromark Cipher Solver")
    parser.add_argument("--base", type=int, default=10, help="Number base to use (default: 10)")
    parser.add_argument("--length", type=int, default=5, help="Primer length (default: 5)")
    parser.add_argument("--processes", type=int, default=None, help="Number of processes to use (default: all available)")
    parser.add_argument("--quadgrams", type=str, default=None, help="Path to quadgram statistics file")
    parser.add_argument("--top", type=int, default=10, help="Number of top candidates to display (default: 10)")
    parser.add_argument("--specific-primers", type=str, help="Comma-separated list of specific primers to try")
    parser.add_argument("--bean-primers", action="store_true", help="Use primers identified in Bean's research paper")
    parser.add_argument("--alt-expansion", action="store_true", help="Use alternative key expansion method")
    parser.add_argument("--min-match-percentage", type=float, default=0.7, help="Minimum percentage of characters that must match for cribs (default: 0.7)")
    
    args = parser.parse_args()
    
    # Load quadgram statistics if provided
    if args.quadgrams:
        load_english_quadgrams(args.quadgrams)
    else:
        load_english_quadgrams()
    
    # Prepare list of specific primers to try
    specific_primers = []
    
    # Add primers from Bean's paper if requested
    if args.bean_primers:
        key = (args.base, args.length)
        if key in PROMISING_PRIMERS:
            specific_primers.extend(PROMISING_PRIMERS[key])
    
    # Add user-specified primers if provided
    if args.specific_primers:
        for primer_str in args.specific_primers.split(','):
            try:
                primer = primer_list_from_string(primer_str.strip(), args.base)
                if len(primer) != args.length:
                    print(f"Warning: Primer {primer_str} length doesn't match specified length {args.length}.")
                specific_primers.append(primer)
            except ValueError as e:
                print(f"Error parsing primer '{primer_str}': {e}")
    
    # Choose expansion function
    expansion_func = fibonacci_expansion if args.alt_expansion else expand_gromark_primer
    
    # Run the brute force search
    candidates = brute_force_gromark(
        base=args.base,
        primer_length=args.length,
        processes=args.processes,
        specific_primers=specific_primers,
        expansion_func=expansion_func,
        min_match_percentage=args.min_match_percentage
    )
    
    # Print results
    print_results(candidates, top_n=args.top)

if __name__ == "__main__":
    main() 