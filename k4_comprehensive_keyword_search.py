#!/usr/bin/env python3
"""
K4 Comprehensive Keyword Search
Tests all of Bean's promising primers with each keyword in the wordlist,
focusing on the most efficient patterns to find optimal combinations.
"""

import string
import time
import argparse
import multiprocessing
from collections import Counter
import os
import sys

# K4 ciphertext
K4_CIPHERTEXT = "OBKRUOXOGHULBSOLIFBBWFLRVQQPRNGKSSOTWTQSJQSSEKZZWATJKLUDIAWINFBNYPVTTMZFPKWGDKZXTJCDIGKUHUAUEKCAR"

# Known plaintext segments and their positions (0-indexed)
KNOWN_PLAINTEXT = {
    # "EASTNORTHEAST" at positions 10-22
    10: ('E', 'F'),
    11: ('A', 'L'),
    12: ('S', 'R'),
    13: ('T', 'V'),
    14: ('N', 'Q'),
    15: ('O', 'Q'),
    16: ('R', 'P'),
    17: ('T', 'R'),
    18: ('H', 'N'),
    19: ('E', 'G'),
    20: ('A', 'K'),
    21: ('S', 'S'),
    22: ('T', 'S'),
    
    # "BERLINCLOCK" at positions 62-72
    62: ('B', 'N'),
    63: ('E', 'Y'),
    64: ('R', 'P'),
    65: ('L', 'V'),
    66: ('I', 'T'),
    67: ('N', 'T'),
    68: ('C', 'M'),
    69: ('L', 'Z'),
    70: ('O', 'F'),
    71: ('C', 'P'),
    72: ('K', 'K')
}

# All promising primers from Bean's research
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

# Top-performing patterns from previous testing
BEST_PATTERNS = [
    {"type": "standard", "desc": "Standard keyword alphabet"},
    {"type": "mod3", "shifts": {0: 6, 1: 3, 2: 20}, "desc": "High IC pattern (6,3,20)"},
    {"type": "mod3", "shifts": {0: 0, 1: 2, 2: 10}, "desc": "Minor differences analysis"},
    {"type": "mod3", "shifts": {0: 7, 1: 2, 2: 20}, "desc": "Modified best pattern (7,2,20)"}
]

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

def apply_mod3_pattern(cipher_alphabet, shifts):
    """
    Apply mod3 pattern to create position-specific alphabets.
    
    Args:
        cipher_alphabet (str): Base cipher alphabet
        shifts (dict): Dictionary mapping position mod 3 to shift values
        
    Returns:
        list: Position-specific alphabets
    """
    position_alphabets = []
    
    for i in range(len(K4_CIPHERTEXT)):
        mod = i % 3
        if mod in shifts:
            shift = shifts[mod]
            if shift == 0:
                position_alphabets.append(cipher_alphabet)
            else:
                shifted = ''.join(chr(((ord(c) - ord('A') + shift) % 26) + ord('A')) for c in cipher_alphabet)
                position_alphabets.append(shifted)
        else:
            position_alphabets.append(cipher_alphabet)
    
    return position_alphabets

def decrypt_gromark(ciphertext, key, plain_alphabets, cipher_alphabets):
    """Decrypt using the Gromark cipher with position-specific alphabets."""
    plaintext = []
    
    for i, c in enumerate(ciphertext):
        if i >= len(plain_alphabets) or i >= len(cipher_alphabets):
            plaintext.append(c)
            continue
            
        plain_alphabet = plain_alphabets[i]
        cipher_alphabet = cipher_alphabets[i]
        
        if c in cipher_alphabet:
            # Get the position in the cipher alphabet
            pos = cipher_alphabet.index(c)
            # Subtract the key digit (modulo alphabet length)
            plain_pos = (pos - key[i]) % len(plain_alphabet)
            # Convert back to plaintext
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

def verify_known_plaintext(plaintext):
    """Verify if the plaintext matches the known plaintext segments."""
    matches = []
    for pos, (expected, cipher) in KNOWN_PLAINTEXT.items():
        if pos < len(plaintext) and plaintext[pos] == expected:
            matches.append((pos, expected, cipher))
    
    match_percentage = len(matches) / len(KNOWN_PLAINTEXT) * 100
    return match_percentage, matches, len(KNOWN_PLAINTEXT)

def load_wordlist(filename, min_length=3, max_length=15, word_filter=None):
    """
    Load words from a wordlist file.
    
    Args:
        filename (str): Path to wordlist file
        min_length (int): Minimum word length to include
        max_length (int): Maximum word length to include
        word_filter (list): Optional list of specific words to include
        
    Returns:
        list: Filtered list of words
    """
    words = []
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            for line in f:
                word = line.strip().upper()
                if min_length <= len(word) <= max_length and all(c in string.ascii_uppercase for c in word):
                    if word_filter is None or word in word_filter:
                        words.append(word)
    except Exception as e:
        print(f"Error loading wordlist: {e}")
    
    print(f"Loaded {len(words)} words from {filename}")
    return words

def test_keyword_with_primer(args):
    """
    Test a keyword with a specific primer and patterns.
    Used with multiprocessing.
    
    Args:
        args (tuple): (word, primer_config, patterns, expansion_methods)
        
    Returns:
        list: Results for this keyword and primer
    """
    word, primer_config, patterns, expansion_methods = args
    
    results = []
    primer = primer_config["primer"]
    base = primer_config["base"]
    desc = primer_config["desc"]
    
    # Create the keyword-mixed alphabet
    keyword_alphabet = create_keyword_alphabet(word)
    
    # Test with each expansion method
    for expansion_method in expansion_methods:
        # Expand the key
        key = expand_key(primer, base, len(K4_CIPHERTEXT), expansion_method)
        
        # Test each pattern
        for pattern in patterns:
            if pattern["type"] == "standard":
                # Standard keyword alphabet (no pattern)
                cipher_alphabets = [keyword_alphabet] * len(K4_CIPHERTEXT)
                plain_alphabets = [string.ascii_uppercase] * len(K4_CIPHERTEXT)
                
                plaintext = decrypt_gromark(K4_CIPHERTEXT, key, plain_alphabets, cipher_alphabets)
                match_percentage, matches, total = verify_known_plaintext(plaintext)
                ic = calculate_ic(plaintext)
                
                if len(matches) >= 3:  # Only record if at least 3 matches
                    results.append({
                        "word": word,
                        "primer": primer,
                        "base": base,
                        "primer_name": primer_config["name"],
                        "desc": desc,
                        "pattern": pattern["desc"],
                        "expansion": expansion_method,
                        "plaintext": plaintext,
                        "match_percentage": match_percentage,
                        "matches": len(matches),
                        "matched_positions": [m[0] for m in matches],
                        "total": total,
                        "ic": ic
                    })
                
                # Also test with reverse standard (keyword as plaintext alphabet)
                cipher_alphabets = [string.ascii_uppercase] * len(K4_CIPHERTEXT)
                plain_alphabets = [keyword_alphabet] * len(K4_CIPHERTEXT)
                
                plaintext = decrypt_gromark(K4_CIPHERTEXT, key, plain_alphabets, cipher_alphabets)
                match_percentage, matches, total = verify_known_plaintext(plaintext)
                ic = calculate_ic(plaintext)
                
                if len(matches) >= 3:  # Only record if at least 3 matches
                    results.append({
                        "word": word,
                        "primer": primer,
                        "base": base,
                        "primer_name": primer_config["name"],
                        "desc": desc,
                        "pattern": "Keyword plaintext alphabet",
                        "expansion": expansion_method,
                        "plaintext": plaintext,
                        "match_percentage": match_percentage,
                        "matches": len(matches),
                        "matched_positions": [m[0] for m in matches],
                        "total": total,
                        "ic": ic
                    })
            
            elif pattern["type"] == "mod3":
                # Apply mod3 pattern to keyword alphabet
                cipher_alphabets = apply_mod3_pattern(keyword_alphabet, pattern["shifts"])
                plain_alphabets = [string.ascii_uppercase] * len(K4_CIPHERTEXT)
                
                plaintext = decrypt_gromark(K4_CIPHERTEXT, key, plain_alphabets, cipher_alphabets)
                match_percentage, matches, total = verify_known_plaintext(plaintext)
                ic = calculate_ic(plaintext)
                
                if len(matches) >= 3:  # Only record if at least 3 matches
                    results.append({
                        "word": word,
                        "primer": primer,
                        "base": base,
                        "primer_name": primer_config["name"],
                        "desc": desc,
                        "pattern": f"Keyword + {pattern['desc']}",
                        "expansion": expansion_method,
                        "plaintext": plaintext,
                        "match_percentage": match_percentage,
                        "matches": len(matches),
                        "matched_positions": [m[0] for m in matches],
                        "total": total,
                        "ic": ic
                    })
    
    return results

def print_results_at_width21(plaintext):
    """Print plaintext arranged at width 21."""
    for i in range(0, len(plaintext), 21):
        print(plaintext[i:i+21])

def print_top_results(results, top_n=20):
    """Print top results."""
    print("\nTOP RESULTS")
    print("=" * 100)
    
    for i, result in enumerate(results[:top_n]):
        print(f"\nRank {i+1}:")
        print(f"  Keyword: {result['word']}")
        print(f"  Primer: {result['primer_name']} {result['primer']} - {result['desc']}")
        print(f"  Pattern: {result['pattern']}")
        print(f"  Expansion: {result['expansion']}")
        print(f"  Match: {result['match_percentage']:.2f}% ({result['matches']}/{result['total']})")
        print(f"  Matched positions: {result['matched_positions']}")
        print(f"  IC: {result['ic']:.6f}")
        print(f"  Plaintext: {result['plaintext'][:50]}...")
        
        print("\n  Plaintext at width 21:")
        print_results_at_width21(result['plaintext'])
        print("-" * 100)

def save_results(results, filename="comprehensive_results.txt"):
    """Save results to a file."""
    with open(filename, 'w', encoding='utf-8') as f:
        f.write(f"K4 Comprehensive Keyword Search Results\n")
        f.write(f"Total results: {len(results)}\n")
        f.write("=" * 80 + "\n\n")
        
        for i, result in enumerate(results[:200]):  # Save top 200 results
            f.write(f"Rank {i+1}:\n")
            f.write(f"  Keyword: {result['word']}\n")
            f.write(f"  Primer: {result['primer_name']} {result['primer']} - {result['desc']}\n")
            f.write(f"  Pattern: {result['pattern']}\n")
            f.write(f"  Expansion: {result['expansion']}\n")
            f.write(f"  Match: {result['match_percentage']:.2f}% ({result['matches']}/{result['total']})\n")
            f.write(f"  Matched positions: {result['matched_positions']}\n")
            f.write(f"  IC: {result['ic']:.6f}\n")
            f.write(f"  Plaintext: {result['plaintext'][:100]}...\n\n")
    
    print(f"Results saved to {filename}")

def analyze_combinations(results, top_n=100):
    """Analyze which keyword-primer combinations perform best."""
    top_results = results[:top_n]
    
    # Count by primer
    primer_counts = Counter([r["primer_name"] for r in top_results])
    print("\nMost common primers in top results:")
    for primer, count in primer_counts.most_common():
        print(f"  {primer}: {count}")
    
    # Count by keyword
    keyword_counts = Counter([r["word"] for r in top_results])
    print("\nMost common keywords in top results:")
    for word, count in keyword_counts.most_common(10):
        print(f"  {word}: {count}")
    
    # Count by pattern
    pattern_counts = Counter([r["pattern"] for r in top_results])
    print("\nMost common patterns in top results:")
    for pattern, count in pattern_counts.most_common():
        print(f"  {pattern}: {count}")
    
    # Group by primer and find best matches
    primer_best = {}
    for result in results:
        primer_name = result["primer_name"]
        if primer_name not in primer_best or result["matches"] > primer_best[primer_name]["matches"]:
            primer_best[primer_name] = result
    
    print("\nBest match for each primer:")
    for primer_name in sorted(primer_best.keys()):
        result = primer_best[primer_name]
        print(f"  {primer_name}: {result['word']} - {result['matches']} matches ({result['match_percentage']:.2f}%), IC={result['ic']:.6f}")
    
    # Group by keyword and find best matches
    keyword_best = {}
    for result in results:
        word = result["word"]
        if word not in keyword_best or result["matches"] > keyword_best[word]["matches"]:
            keyword_best[word] = result
    
    # Sort by number of matches
    best_keywords = sorted(keyword_best.values(), key=lambda x: (x["matches"], x["ic"]), reverse=True)
    
    print("\nTop 10 keywords by best match:")
    for i, result in enumerate(best_keywords[:10]):
        print(f"  {i+1}. {result['word']} with {result['primer_name']}: {result['match_percentage']:.2f}% matches, IC={result['ic']:.6f}")

def save_checkpoint(results, filename="search_checkpoint.txt"):
    """Save a checkpoint with current top results."""
    with open(filename, 'w', encoding='utf-8') as f:
        f.write(f"Checkpoint - {len(results)} results found\n")
        f.write("=" * 80 + "\n\n")
        
        for i, result in enumerate(results[:20]):  # Save top 20 results
            f.write(f"Rank {i+1}: {result['word']} + {result['primer_name']} = {result['matches']} matches, IC={result['ic']:.6f}\n")
    
    print(f"Checkpoint saved to {filename}")

def main():
    """Main function to test all keywords with all primers."""
    parser = argparse.ArgumentParser(description="Comprehensive Keyword Search for K4")
    parser.add_argument("--wordlist", default="not.txt", help="Path to wordlist file")
    parser.add_argument("--min-length", type=int, default=4, help="Minimum word length")
    parser.add_argument("--max-length", type=int, default=15, help="Maximum word length")
    parser.add_argument("--processes", type=int, default=None, help="Number of processes to use")
    parser.add_argument("--words", type=str, help="Comma-separated list of specific words to test")
    parser.add_argument("--primers", type=str, help="Comma-separated list of specific primers to test")
    parser.add_argument("--batch-size", type=int, default=100, help="Batch size for processing")
    parser.add_argument("--min-matches", type=int, default=3, help="Minimum matches to record")
    parser.add_argument("--output", default="comprehensive_results.txt", help="Output file for results")
    parser.add_argument("--checkpoint-interval", type=int, default=10, help="Save checkpoint every N batches")
    
    args = parser.parse_args()
    
    print("K4 COMPREHENSIVE KEYWORD SEARCH")
    print("=" * 50)
    
    # Set up specific words if requested
    word_filter = None
    if args.words:
        word_filter = [w.strip().upper() for w in args.words.split(",")]
        print(f"Testing specific words: {word_filter}")
    
    # Set up specific primers if requested
    primers_to_use = BEAN_PRIMERS
    if args.primers:
        primer_names = [p.strip() for p in args.primers.split(",")]
        primers_to_use = [p for p in BEAN_PRIMERS if p["name"] in primer_names]
        print(f"Testing specific primers: {[p['name'] for p in primers_to_use]}")
    
    # Load wordlist
    words = load_wordlist(args.wordlist, args.min_length, args.max_length, word_filter)
    
    if not words:
        print("No words to test. Exiting.")
        return
    
    # Setup for multiprocessing
    num_processes = args.processes if args.processes else max(1, multiprocessing.cpu_count() - 1)
    print(f"Using {num_processes} processes")
    
    # Decide on expansion methods
    expansion_methods = ["standard", "fibonacci"]
    
    # Statistics
    all_results = []
    tested_combinations = 0
    total_combinations = len(words) * len(primers_to_use)
    start_time = time.time()
    
    # Process words in batches
    batch_size = min(args.batch_size, len(words))
    num_batches = (len(words) + batch_size - 1) // batch_size
    
    print(f"\nTesting {len(words)} words with {len(primers_to_use)} primers")
    print(f"Total combinations to test: {total_combinations}")
    print(f"Processing in {num_batches} batches of up to {batch_size} words")
    
    for batch_idx in range(num_batches):
        batch_start = batch_idx * batch_size
        batch_end = min(batch_start + batch_size, len(words))
        batch_words = words[batch_start:batch_end]
        
        batch_start_time = time.time()
        print(f"\nBatch {batch_idx+1}/{num_batches}: Testing words {batch_start+1}-{batch_end}")
        
        # Prepare tasks for multiprocessing
        tasks = []
        for word in batch_words:
            for primer_config in primers_to_use:
                tasks.append((word, primer_config, BEST_PATTERNS, expansion_methods))
        
        # Process batch using multiprocessing
        with multiprocessing.Pool(processes=num_processes) as pool:
            batch_results = pool.map(test_keyword_with_primer, tasks)
        
        # Flatten results
        batch_results = [item for sublist in batch_results for item in sublist]
        
        # Filter by minimum matches
        batch_results = [r for r in batch_results if r["matches"] >= args.min_matches]
        
        # Add to all results
        all_results.extend(batch_results)
        
        # Sort by matches then IC
        all_results.sort(key=lambda x: (x["matches"], x["ic"]), reverse=True)
        
        tested_combinations += len(tasks)
        batch_end_time = time.time()
        batch_time = batch_end_time - batch_start_time
        
        # Print batch summary
        print(f"Batch completed in {batch_time:.2f} seconds")
        print(f"Progress: {tested_combinations}/{total_combinations} combinations tested ({tested_combinations/total_combinations*100:.1f}%)")
        print(f"Found {len(batch_results)} results with {args.min_matches}+ matches in this batch")
        
        # Show best results in this batch
        if batch_results:
            best = batch_results[0]
            print(f"Best in batch: {best['word']} + {best['primer_name']} = {best['matches']} matches, IC={best['ic']:.6f}")
        
        # Save checkpoint periodically
        if batch_idx % args.checkpoint_interval == 0:
            save_checkpoint(all_results)
    
    end_time = time.time()
    print(f"\nSearch completed in {end_time - start_time:.2f} seconds")
    print(f"Tested {tested_combinations} keyword-primer combinations")
    print(f"Found {len(all_results)} results with {args.min_matches}+ matches")
    
    # Print top results
    print_top_results(all_results, top_n=20)
    
    # Analyze combinations
    analyze_combinations(all_results)
    
    # Save results
    save_results(all_results, args.output)

if __name__ == "__main__":
    main() 