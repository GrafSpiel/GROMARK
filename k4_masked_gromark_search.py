#!/usr/bin/env python3
"""
K4 Masked Gromark Search
Implements a comprehensive search approach combining Gromark decryption with additional transformations:
1. Analyzing letter frequency for English-like patterns
2. Testing alphabetical substitution (transposition)
3. Testing Caesar cipher shifts

This approach addresses the possibility that Sanborn may have applied additional masking
techniques to the plaintext after Gromark encryption.
"""

import string
import time
import argparse
import multiprocessing
import os
from collections import Counter, defaultdict
import itertools
import numpy as np
from tqdm import tqdm

# K4 ciphertext
K4_CIPHERTEXT = "OBKRUOXOGHULBSOLIFBBWFLRVQQPRNGKSSOTWTQSJQSSEKZZWATJKLUDIAWINFBNYPVTTMZFPKWGDKZXTJCDIGKUHUAUEKCAR"

# Known plaintext fragments
KNOWN_PLAINTEXT = {
    # "EASTNORTHEAST" at positions 21-33
    21: 'E',
    22: 'A',
    23: 'S',
    24: 'T',
    25: 'N',
    26: 'O',
    27: 'R',
    28: 'T',
    29: 'H',
    30: 'E',
    31: 'A',
    32: 'S',
    33: 'T',
    
    # "BERLINCLOCK" at positions 62-72
    62: 'B',
    63: 'E',
    64: 'R',
    65: 'L',
    66: 'I',
    67: 'N',
    68: 'C',
    69: 'L',
    70: 'O',
    71: 'C',
    72: 'K'
}

# Bean's promising primers
BEAN_PRIMERS = [
    {"name": "26717", "primer": [2, 6, 7, 1, 7], "base": 10, "desc": "Equivalent primer 1"},
    {"name": "98800", "primer": [9, 8, 8, 0, 0], "base": 10, "desc": "Highest IC value"},
    {"name": "84393", "primer": [8, 4, 3, 9, 3], "base": 10, "desc": "Equivalent primer 2"},
    {"name": "71688", "primer": [7, 1, 6, 8, 8], "base": 10, "desc": "Bean research primer 4"},
    {"name": "58664", "primer": [5, 8, 6, 6, 4], "base": 10, "desc": "Bean research primer 5"},
    {"name": "36833", "primer": [3, 6, 8, 3, 3], "base": 10, "desc": "Bean research primer 6"},
    {"name": "31861", "primer": [3, 1, 8, 6, 1], "base": 10, "desc": "Bean research primer 7"},
    {"name": "28088", "primer": [2, 8, 0, 8, 8], "base": 10, "desc": "Bean research primer 8"}
]

# English letter frequency (approximate percentages)
ENGLISH_LETTER_FREQ = {
    'E': 12.02, 'T': 9.10, 'A': 8.12, 'O': 7.68, 'I': 7.31, 'N': 6.95, 'S': 6.28, 'R': 6.02, 
    'H': 5.92, 'D': 4.32, 'L': 3.98, 'U': 2.88, 'C': 2.71, 'M': 2.61, 'F': 2.30, 'Y': 2.11, 
    'W': 2.09, 'G': 2.03, 'P': 1.82, 'B': 1.49, 'V': 1.11, 'K': 0.69, 'X': 0.17, 'Q': 0.11,
    'J': 0.10, 'Z': 0.07
}

def load_wordlist(filename, min_length=3, max_length=15):
    """Load words from a file."""
    words = []
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            for line in f:
                word = line.strip().upper()
                if min_length <= len(word) <= max_length and all(c in string.ascii_uppercase for c in word):
                    words.append(word)
        print(f"Loaded {len(words)} words from {filename}")
    except Exception as e:
        print(f"Error loading wordlist: {e}")
    
    return words

def create_keyword_alphabet(keyword):
    """Create a mixed alphabet based on a keyword."""
    # Remove duplicates while preserving order
    keyword = keyword.upper().replace(" ", "")
    unique_chars = []
    for c in keyword:
        if c in string.ascii_uppercase and c not in unique_chars:
            unique_chars.append(c)
    
    # Create the mixed alphabet
    remaining = [c for c in string.ascii_uppercase if c not in unique_chars]
    return ''.join(unique_chars + remaining)

def expand_key(primer, base, length, method="standard"):
    """Expand a primer to generate a key sequence."""
    key = list(primer.copy())
    
    while len(key) < length:
        if method == "standard":
            # Standard Gromark: add last two digits
            next_digit = (key[-1] + key[-2]) % base
            key.append(next_digit)
        elif method == "fibonacci":
            # Sum all digits in the primer
            next_digit = 0
            for i in range(1, len(primer) + 1):
                if i <= len(key):
                    next_digit = (next_digit + key[-i]) % base
            key.append(next_digit)
        else:
            # Default to standard
            next_digit = (key[-1] + key[-2]) % base
            key.append(next_digit)
            
    return key[:length]

def decrypt_gromark(ciphertext, key, plain_alphabet, cipher_alphabet):
    """Decrypt using the Gromark cipher."""
    plaintext = []
    
    for i, c in enumerate(ciphertext):
        if c in cipher_alphabet:
            # Find the position in the cipher alphabet
            pos = cipher_alphabet.index(c)
            # Subtract the key digit (modulo alphabet length)
            plain_pos = (pos - key[i]) % len(plain_alphabet)
            # Convert to plaintext
            plaintext.append(plain_alphabet[plain_pos])
        else:
            # Pass through non-alphabet characters
            plaintext.append(c)
    
    return ''.join(plaintext)

def calculate_ic(text):
    """Calculate the Index of Coincidence for a text."""
    counts = Counter(c for c in text.upper() if c in string.ascii_uppercase)
    n = sum(counts.values())
    if n <= 1:
        return 0
    
    sum_freqs = sum(count * (count - 1) for count in counts.values())
    return sum_freqs / (n * (n - 1))

def calculate_letter_frequency_fit(text):
    """
    Calculate how closely the text's letter frequency matches English.
    Returns a score from 0 to 1, where 1 is a perfect match.
    """
    # Count letter occurrences
    counts = Counter(c for c in text.upper() if c in string.ascii_uppercase)
    total = sum(counts.values())
    
    if total == 0:
        return 0
    
    # Calculate frequency for each letter
    freq = {letter: (count / total) * 100 for letter, count in counts.items()}
    
    # Calculate chi-squared statistic (lower is better)
    chi_squared = 0
    for letter in string.ascii_uppercase:
        observed = freq.get(letter, 0)
        expected = ENGLISH_LETTER_FREQ.get(letter, 0.01)  # Avoid division by zero
        chi_squared += ((observed - expected) ** 2) / expected
    
    # Convert to a score from 0 to 1 (using an exponential falloff)
    # A chi-squared of 0 means perfect match, higher values mean worse match
    score = np.exp(-chi_squared / 100)  # Adjust the divisor to scale appropriately
    
    return score

def apply_caesar_shift(text, shift):
    """Apply a Caesar cipher shift to the text."""
    result = []
    for char in text:
        if char in string.ascii_uppercase:
            # Convert to 0-25, apply shift, convert back to ASCII
            idx = (ord(char) - ord('A') + shift) % 26
            result.append(chr(idx + ord('A')))
        else:
            result.append(char)
    return ''.join(result)

def apply_alphabet_substitution(text, substitution_map):
    """Apply an alphabetical substitution to the text."""
    result = []
    for char in text:
        if char in string.ascii_uppercase:
            result.append(substitution_map.get(char, char))
        else:
            result.append(char)
    return ''.join(result)

def generate_column_transposition_key(length):
    """Generate a columnar transposition key of the specified length."""
    # Create all possible permutations of columns
    cols = list(range(length))
    # For practicality, limit to a reasonable number of permutations
    if length <= 7:  # Full permutation for small lengths
        return list(itertools.permutations(cols))
    else:
        # For longer keys, just use some basic permutations
        results = []
        # Add identity permutation
        results.append(tuple(cols))
        # Add reverse
        results.append(tuple(reversed(cols)))
        # Add some shifts
        for i in range(1, length):
            results.append(tuple(cols[i:] + cols[:i]))
        # Add some swaps
        for i in range(length - 1):
            new_cols = cols.copy()
            new_cols[i], new_cols[i+1] = new_cols[i+1], new_cols[i]
            results.append(tuple(new_cols))
        return results

def apply_columnar_transposition(text, key_permutation):
    """Apply a columnar transposition to the text using the given key permutation."""
    # Determine number of rows and columns
    key_length = len(key_permutation)
    num_rows = (len(text) + key_length - 1) // key_length
    
    # Pad the text if necessary
    padded_text = text.ljust(num_rows * key_length, 'X')
    
    # Fill the grid by columns according to the key
    grid = [[''] * key_length for _ in range(num_rows)]
    for i, char in enumerate(padded_text):
        row = i // key_length
        col = i % key_length
        grid[row][col] = char
    
    # Read out the transposed text
    result = []
    for col_idx in key_permutation:
        for row in range(num_rows):
            result.append(grid[row][col_idx])
    
    return ''.join(result)

def verify_known_plaintext(plaintext):
    """
    Check how many known plaintext characters match.
    Returns (match_percentage, matched_positions).
    """
    matches = []
    
    for pos, char in KNOWN_PLAINTEXT.items():
        if pos < len(plaintext) and plaintext[pos] == char:
            matches.append(pos)
    
    match_percentage = len(matches) / len(KNOWN_PLAINTEXT) * 100
    return match_percentage, matches

def test_gromark_with_word(args):
    """
    Test a word with Gromark decryption and various post-processing techniques.
    Used with multiprocessing.
    """
    word, primer_config, expansion_methods, post_process_options = args
    
    results = {
        "frequency": [],
        "caesar": [],
        "substitution": []
    }
    
    primer = primer_config["primer"]
    base = primer_config["base"]
    primer_desc = primer_config["desc"]
    
    # Create the keyword-mixed alphabet
    keyword_alphabet = create_keyword_alphabet(word)
    
    # Test with each expansion method
    for expansion_method in expansion_methods:
        # Expand the key
        key = expand_key(primer, base, len(K4_CIPHERTEXT), expansion_method)
        
        # Test with standard alphabet as plaintext, mixed as cipher
        plaintext1 = decrypt_gromark(
            K4_CIPHERTEXT,
            key,
            string.ascii_uppercase,  # Standard plaintext alphabet
            keyword_alphabet         # Mixed cipher alphabet
        )
        
        # Test with mixed alphabet as plaintext, standard as cipher
        plaintext2 = decrypt_gromark(
            K4_CIPHERTEXT,
            key,
            keyword_alphabet,        # Mixed plaintext alphabet
            string.ascii_uppercase   # Standard cipher alphabet
        )
        
        for plaintext_idx, plaintext in enumerate([plaintext1, plaintext2]):
            # Calculate base metrics
            ic = calculate_ic(plaintext)
            match_percentage, matched_positions = verify_known_plaintext(plaintext)
            
            # Skip if doesn't meet minimum threshold to save computation
            if match_percentage < post_process_options["min_match_percentage"] and ic < post_process_options["min_ic"]:
                continue
            
            # 1. Analyze letter frequency
            freq_score = calculate_letter_frequency_fit(plaintext)
            
            if freq_score >= post_process_options["min_freq_score"]:
                results["frequency"].append({
                    "word": word,
                    "primer": primer,
                    "primer_name": primer_config["name"],
                    "desc": primer_desc,
                    "expansion": expansion_method,
                    "alphabet_config": "standard-mixed" if plaintext_idx == 0 else "mixed-standard",
                    "plaintext": plaintext,
                    "ic": ic,
                    "match_percentage": match_percentage,
                    "matched_positions": matched_positions,
                    "freq_score": freq_score,
                    "transformation": "none"
                })
            
            # Only proceed with transformations if there's some promise in the plaintext
            if match_percentage >= 5 or ic >= 0.035:
                # 2. Test Caesar shifts
                for shift in range(1, 26):
                    shifted_text = apply_caesar_shift(plaintext, shift)
                    shift_ic = calculate_ic(shifted_text)
                    shift_match_pct, shift_matches = verify_known_plaintext(shifted_text)
                    
                    if shift_match_pct >= post_process_options["min_match_percentage"] or shift_ic >= post_process_options["min_ic"]:
                        results["caesar"].append({
                            "word": word,
                            "primer": primer,
                            "primer_name": primer_config["name"],
                            "desc": primer_desc,
                            "expansion": expansion_method,
                            "alphabet_config": "standard-mixed" if plaintext_idx == 0 else "mixed-standard",
                            "plaintext": shifted_text,
                            "original_plaintext": plaintext,
                            "ic": shift_ic,
                            "match_percentage": shift_match_pct,
                            "matched_positions": shift_matches,
                            "freq_score": calculate_letter_frequency_fit(shifted_text),
                            "transformation": f"caesar-{shift}"
                        })
                
                # 3. Test basic substitution patterns
                # For simplicity, just test a few key substitutions
                substitution_tests = [
                    # Reverse alphabet
                    {c: string.ascii_uppercase[25-i] for i, c in enumerate(string.ascii_uppercase)},
                    # A=Z, B=Y pattern
                    {c: chr(ord('Z') - (ord(c) - ord('A'))) for c in string.ascii_uppercase},
                    # Simple shifts of the keyword alphabet
                    {c: keyword_alphabet[(keyword_alphabet.index(c) + 13) % 26] 
                     if c in keyword_alphabet else c for c in string.ascii_uppercase},
                ]
                
                for sub_idx, sub_map in enumerate(substitution_tests):
                    sub_text = apply_alphabet_substitution(plaintext, sub_map)
                    sub_ic = calculate_ic(sub_text)
                    sub_match_pct, sub_matches = verify_known_plaintext(sub_text)
                    
                    if sub_match_pct >= post_process_options["min_match_percentage"] or sub_ic >= post_process_options["min_ic"]:
                        results["substitution"].append({
                            "word": word,
                            "primer": primer,
                            "primer_name": primer_config["name"],
                            "desc": primer_desc,
                            "expansion": expansion_method,
                            "alphabet_config": "standard-mixed" if plaintext_idx == 0 else "mixed-standard",
                            "plaintext": sub_text,
                            "original_plaintext": plaintext,
                            "ic": sub_ic,
                            "match_percentage": sub_match_pct,
                            "matched_positions": sub_matches,
                            "freq_score": calculate_letter_frequency_fit(sub_text),
                            "transformation": f"substitution-{sub_idx}"
                        })
    
    return results

def print_results(results, category, top_n=20):
    """Print top results for a category."""
    print(f"\nTOP {category.upper()} RESULTS")
    print("=" * 100)
    
    for i, result in enumerate(results[:top_n]):
        print(f"\nRank {i+1}:")
        print(f"  Word: {result['word']}")
        print(f"  Primer: {result['primer_name']} {result['primer']}")
        print(f"  Expansion: {result['expansion']}")
        print(f"  Config: {result['alphabet_config']}")
        print(f"  Transformation: {result.get('transformation', 'none')}")
        print(f"  Match: {result['match_percentage']:.2f}%")
        print(f"  IC: {result['ic']:.6f}")
        print(f"  Frequency Score: {result['freq_score']:.6f}")
        print(f"  Plaintext: {result['plaintext'][:50]}...")
        print("-" * 100)

def save_results(results, filename):
    """Save results to a file."""
    with open(filename, 'w', encoding='utf-8') as f:
        f.write(f"K4 Masked Gromark Search Results\n")
        f.write(f"Total results: {len(results)}\n")
        f.write("=" * 80 + "\n\n")
        
        for i, result in enumerate(results[:500]):  # Save top 500 results
            f.write(f"Rank {i+1}:\n")
            f.write(f"  Word: {result['word']}\n")
            f.write(f"  Primer: {result['primer_name']} {result['primer']}\n")
            f.write(f"  Expansion: {result['expansion']}\n")
            f.write(f"  Config: {result['alphabet_config']}\n")
            f.write(f"  Transformation: {result.get('transformation', 'none')}\n")
            f.write(f"  Match: {result['match_percentage']:.2f}%\n")
            f.write(f"  IC: {result['ic']:.6f}\n")
            f.write(f"  Frequency Score: {result['freq_score']:.6f}\n")
            f.write(f"  Plaintext: {result['plaintext']}\n\n")
    
    print(f"Results saved to {filename}")

def main():
    """Main function to run the masked Gromark search."""
    parser = argparse.ArgumentParser(description="K4 Masked Gromark Search")
    parser.add_argument("--wordlist", default="not.txt", help="Path to wordlist file")
    parser.add_argument("--min-length", type=int, default=4, help="Minimum word length")
    parser.add_argument("--max-length", type=int, default=15, help="Maximum word length")
    parser.add_argument("--processes", type=int, default=None, help="Number of processes to use")
    parser.add_argument("--words", type=str, help="Comma-separated list of specific words to test")
    parser.add_argument("--primers", type=str, help="Comma-separated list of specific primers to test")
    parser.add_argument("--batch-size", type=int, default=50, help="Batch size for processing")
    parser.add_argument("--min-match", type=float, default=10.0, help="Minimum match percentage to record")
    parser.add_argument("--min-ic", type=float, default=0.035, help="Minimum IC to record")
    parser.add_argument("--min-freq", type=float, default=0.2, help="Minimum frequency score to record")
    parser.add_argument("--output-dir", default="results", help="Directory for output files")
    
    args = parser.parse_args()
    
    print("K4 MASKED GROMARK SEARCH")
    print("=" * 50)
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Define output files
    freq_output = os.path.join(args.output_dir, "gromark_frequency_results.txt")
    caesar_output = os.path.join(args.output_dir, "gromark_caesar_results.txt")
    subst_output = os.path.join(args.output_dir, "gromark_substitution_results.txt")
    
    # Set up specific words if requested
    word_filter = None
    if args.words:
        word_filter = [w.strip().upper() for w in args.words.split(",")]
        print(f"Testing specific words: {word_filter}")
    
    # Load wordlist
    words = []
    if word_filter:
        words = word_filter
    else:
        words = load_wordlist(args.wordlist, args.min_length, args.max_length)
    
    if not words:
        print("No words to test. Exiting.")
        return
    
    # Set up specific primers if requested
    primers_to_use = BEAN_PRIMERS
    if args.primers:
        primer_names = [p.strip() for p in args.primers.split(",")]
        primers_to_use = [p for p in BEAN_PRIMERS if p["name"] in primer_names]
        print(f"Testing specific primers: {[p['name'] for p in primers_to_use]}")
    
    # Setup for multiprocessing
    num_processes = args.processes if args.processes else max(1, multiprocessing.cpu_count() - 1)
    print(f"Using {num_processes} processes")
    
    # Define expansion methods
    expansion_methods = ["standard", "fibonacci"]
    
    # Define post-processing options
    post_process_options = {
        "min_match_percentage": args.min_match,
        "min_ic": args.min_ic,
        "min_freq_score": args.min_freq
    }
    
    # Statistics
    all_results = {
        "frequency": [],
        "caesar": [],
        "substitution": []
    }
    
    start_time = time.time()
    
    # Calculate total tasks for progress bar
    total_tasks = len(words) * len(primers_to_use)
    
    # Process words in batches with progress bar
    batch_size = min(args.batch_size, len(words))
    num_batches = (len(words) + batch_size - 1) // batch_size
    
    print(f"\nTesting {len(words)} words with {len(primers_to_use)} primers")
    print(f"Total combinations to test: {total_tasks}")
    print(f"Processing in {num_batches} batches of up to {batch_size} words")
    
    progress_bar = tqdm(total=total_tasks, desc="Testing combinations", 
                      bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}{postfix}]")
    
    for batch_idx in range(num_batches):
        batch_start = batch_idx * batch_size
        batch_end = min(batch_start + batch_size, len(words))
        batch_words = words[batch_start:batch_end]
        
        # Prepare tasks for multiprocessing
        tasks = []
        for word in batch_words:
            for primer_config in primers_to_use:
                tasks.append((word, primer_config, expansion_methods, post_process_options))
        
        # Process batch using multiprocessing
        with multiprocessing.Pool(processes=num_processes) as pool:
            batch_results = pool.map(test_gromark_with_word, tasks)
        
        # Merge batch results with all results
        for category in all_results.keys():
            for result_dict in batch_results:
                all_results[category].extend(result_dict[category])
        
        # Update progress bar
        progress_bar.update(len(tasks))
        
        # Sort results by match percentage and IC
        for category in all_results.keys():
            all_results[category].sort(key=lambda x: (x["match_percentage"], x["ic"]), reverse=True)
        
        # Show interim summary
        batch_freq_count = sum(1 for r in batch_results if len(r["frequency"]) > 0)
        batch_caesar_count = sum(1 for r in batch_results if len(r["caesar"]) > 0)
        batch_subst_count = sum(1 for r in batch_results if len(r["substitution"]) > 0)
        
        tqdm.write(f"Batch {batch_idx+1}/{num_batches}: Found {batch_freq_count} frequency matches, "
                f"{batch_caesar_count} Caesar matches, {batch_subst_count} substitution matches")
    
    progress_bar.close()
    
    end_time = time.time()
    print(f"\nSearch completed in {end_time - start_time:.2f} seconds")
    print(f"Found {len(all_results['frequency'])} frequency matches, "
          f"{len(all_results['caesar'])} Caesar matches, "
          f"{len(all_results['substitution'])} substitution matches")
    
    # Print top results
    for category in all_results.keys():
        print_results(all_results[category], category)
    
    # Save results
    save_results(all_results["frequency"], freq_output)
    save_results(all_results["caesar"], caesar_output)
    save_results(all_results["substitution"], subst_output)
    
    # Print statistics for each category
    print("\nSTATISTICS")
    print("=" * 50)
    
    for category in all_results.keys():
        results = all_results[category]
        if results:
            top_words = Counter(r["word"] for r in results[:100]).most_common(5)
            top_primers = Counter(r["primer_name"] for r in results[:100]).most_common(3)
            top_expansions = Counter(r["expansion"] for r in results[:100]).most_common(2)
            
            print(f"\n{category.capitalize()} Analysis:")
            print(f"  Top 5 words: {', '.join(f'{word} ({count})' for word, count in top_words)}")
            print(f"  Top 3 primers: {', '.join(f'{primer} ({count})' for primer, count in top_primers)}")
            print(f"  Top 2 expansion methods: {', '.join(f'{exp} ({count})' for exp, count in top_expansions)}")
            
            max_match = max(results, key=lambda x: x["match_percentage"])
            max_ic = max(results, key=lambda x: x["ic"])
            max_freq = max(results, key=lambda x: x["freq_score"])
            
            print(f"  Best match: {max_match['word']} with {max_match['primer_name']} - {max_match['match_percentage']:.2f}%")
            print(f"  Highest IC: {max_ic['word']} with {max_ic['primer_name']} - {max_ic['ic']:.6f}")
            print(f"  Best frequency fit: {max_freq['word']} with {max_freq['primer_name']} - {max_freq['freq_score']:.6f}")

if __name__ == "__main__":
    main() 