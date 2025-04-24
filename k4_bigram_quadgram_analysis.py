#!/usr/bin/env python3
"""
K4 Bigram and Quadgram Analysis
Analyzes potential Kryptos K4 solutions by examining bigrams and quadgrams
rather than single letter matches. Uses statistics from English text to evaluate
the quality of decryptions.
"""

import string
import time
import argparse
import multiprocessing
from collections import Counter, defaultdict
import os
import sys
import math
import re

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

# Known plaintext bigrams (pairs of consecutive characters)
KNOWN_BIGRAMS = [
    ('EA', (10, 11)),
    ('AS', (11, 12)),
    ('ST', (12, 13)),
    ('TN', (13, 14)),
    ('NO', (14, 15)),
    ('OR', (15, 16)),
    ('RT', (16, 17)),
    ('TH', (17, 18)),
    ('HE', (18, 19)),
    ('EA', (19, 20)),
    ('AS', (20, 21)),
    ('ST', (21, 22)),
    
    ('BE', (62, 63)),
    ('ER', (63, 64)),
    ('RL', (64, 65)),
    ('LI', (65, 66)),
    ('IN', (66, 67)),
    ('NC', (67, 68)),
    ('CL', (68, 69)),
    ('LO', (69, 70)),
    ('OC', (70, 71)),
    ('CK', (71, 72))
]

# Known quadgrams (groups of 4 consecutive characters)
KNOWN_QUADGRAMS = [
    ('EAST', (10, 13)),
    ('STNO', (12, 15)),
    ('NORT', (14, 17)),
    ('RTHE', (16, 19)),
    ('HEAS', (18, 21)),
    ('EAST', (19, 22)),
    
    ('BERL', (62, 65)),
    ('RLIN', (64, 67)),
    ('INCL', (66, 69)),
    ('CLOC', (68, 71)),
    ('LOCK', (69, 72))
]

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

# Most common English bigrams with relative frequencies
ENGLISH_BIGRAMS = {
    'TH': 0.0356, 'HE': 0.0307, 'IN': 0.0243, 'ER': 0.0205, 'AN': 0.0199,
    'RE': 0.0185, 'ON': 0.0176, 'AT': 0.0149, 'EN': 0.0145, 'ND': 0.0135,
    'TI': 0.0134, 'ES': 0.0134, 'OR': 0.0128, 'TE': 0.0120, 'OF': 0.0115,
    'ED': 0.0117, 'IS': 0.0113, 'IT': 0.0112, 'AL': 0.0109, 'AR': 0.0107,
    'ST': 0.0105, 'TO': 0.0104, 'NT': 0.0104, 'NG': 0.0095, 'SE': 0.0093,
    'HA': 0.0093, 'AS': 0.0087, 'OU': 0.0087, 'IO': 0.0083, 'LE': 0.0083,
    'VE': 0.0083, 'CO': 0.0079, 'ME': 0.0079, 'DE': 0.0076, 'HI': 0.0076,
    'RI': 0.0073, 'RO': 0.0073, 'IC': 0.0070, 'NE': 0.0069, 'EA': 0.0069,
    'RA': 0.0069, 'CE': 0.0068, 'LI': 0.0062, 'CH': 0.0060, 'LL': 0.0058,
    'BE': 0.0058, 'MA': 0.0056, 'SI': 0.0055, 'OM': 0.0055, 'UR': 0.0054
}

# Top English quadgrams with log probabilities
# Load these from a file for better coverage
ENGLISH_QUADGRAMS = {}

def load_quadgram_stats(filename="english_quadgrams.txt"):
    """Load quadgram statistics from a file."""
    global ENGLISH_QUADGRAMS
    
    try:
        # If file doesn't exist, we'll use a small default set
        if not os.path.exists(filename):
            print(f"Quadgram file {filename} not found, using default set.")
            ENGLISH_QUADGRAMS = {
                'TION': -8.41, 'NTHE': -8.90, 'THER': -9.02, 'THAT': -9.05,
                'OFTH': -9.28, 'FTHE': -9.31, 'THES': -9.52, 'WITH': -9.55,
                'INTH': -9.59, 'ATIO': -9.82, 'THEC': -10.17, 'ETHE': -10.23
            }
            return
        
        total = 0
        counts = {}
        
        with open(filename, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    quadgram = parts[0].upper()
                    count = int(parts[1])
                    counts[quadgram] = count
                    total += count
        
        # Convert counts to log probabilities
        for quadgram, count in counts.items():
            ENGLISH_QUADGRAMS[quadgram] = math.log10(count / total)
        
        print(f"Loaded {len(ENGLISH_QUADGRAMS)} quadgrams from {filename}")
    
    except Exception as e:
        print(f"Error loading quadgram statistics: {e}")
        # Fallback to default set
        ENGLISH_QUADGRAMS = {
            'TION': -8.41, 'NTHE': -8.90, 'THER': -9.02, 'THAT': -9.05,
            'OFTH': -9.28, 'FTHE': -9.31, 'THES': -9.52, 'WITH': -9.55,
            'INTH': -9.59, 'ATIO': -9.82, 'THEC': -10.17, 'ETHE': -10.23
        }

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

def extract_bigrams(text):
    """Extract all bigrams (letter pairs) from text."""
    text = text.upper()
    return [text[i:i+2] for i in range(len(text)-1) if text[i:i+2].isalpha() and len(text[i:i+2]) == 2]

def extract_quadgrams(text):
    """Extract all quadgrams (4-letter sequences) from text."""
    text = text.upper()
    return [text[i:i+4] for i in range(len(text)-3) if text[i:i+4].isalpha() and len(text[i:i+4]) == 4]

def verify_known_bigrams(plaintext):
    """Check how many known bigrams match in the decrypted plaintext."""
    matches = []
    
    for bigram, (pos1, pos2) in KNOWN_BIGRAMS:
        if pos2 < len(plaintext) and plaintext[pos1:pos2+1] == bigram:
            matches.append((bigram, (pos1, pos2)))
    
    match_percentage = len(matches) / len(KNOWN_BIGRAMS) * 100
    return match_percentage, matches, len(KNOWN_BIGRAMS)

def verify_known_quadgrams(plaintext):
    """Check how many known quadgrams match in the decrypted plaintext."""
    matches = []
    
    for quadgram, (start, end) in KNOWN_QUADGRAMS:
        if end < len(plaintext) and plaintext[start:end+1] == quadgram:
            matches.append((quadgram, (start, end)))
    
    match_percentage = len(matches) / len(KNOWN_QUADGRAMS) * 100
    return match_percentage, matches, len(KNOWN_QUADGRAMS)

def score_bigrams(text):
    """
    Score text based on English bigram frequencies.
    Returns a normalized score between 0 and 1.
    """
    bigrams = extract_bigrams(text)
    if not bigrams:
        return 0
    
    bigram_counts = Counter(bigrams)
    total_bigrams = len(bigrams)
    
    # Calculate score based on how closely the text's bigram distribution
    # matches English bigram frequencies
    score = 0
    for bigram, count in bigram_counts.items():
        freq = count / total_bigrams
        if bigram in ENGLISH_BIGRAMS:
            # Reward common English bigrams
            score += freq * ENGLISH_BIGRAMS[bigram] * 10
    
    # Normalize the score (roughly)
    return min(1.0, score / 0.2)  # 0.2 is approximate max raw score for English text

def score_quadgrams(text):
    """
    Score text based on English quadgram log probabilities.
    Higher score indicates more English-like text.
    """
    quadgrams = extract_quadgrams(text)
    if not quadgrams:
        return -float('inf')
    
    # Use log probabilities for scoring
    score = 0
    floor = -15  # Default log probability for unseen quadgrams
    
    for quad in quadgrams:
        if quad in ENGLISH_QUADGRAMS:
            score += ENGLISH_QUADGRAMS[quad]
        else:
            score += floor
    
    return score

def test_keyword_with_primer(args):
    """
    Test a keyword with a specific primer, using bigram and quadgram analysis.
    Used with multiprocessing.
    
    Args:
        args (tuple): (word, primer_config, expansion_methods)
        
    Returns:
        list: Results for this keyword and primer
    """
    word, primer_config, expansion_methods = args
    
    results = []
    primer = primer_config["primer"]
    base = primer_config["base"]
    primer_desc = primer_config["desc"]
    
    # Create the keyword-mixed alphabet
    keyword_alphabet = create_keyword_alphabet(word)
    
    # Test with each expansion method
    for expansion_method in expansion_methods:
        # Expand the key
        key = expand_key(primer, base, len(K4_CIPHERTEXT), expansion_method)
        
        # Test standard configuration (mixed alphabet as cipher alphabet)
        plaintext = decrypt_gromark(
            K4_CIPHERTEXT, 
            key, 
            string.ascii_uppercase,  # Standard plaintext alphabet
            keyword_alphabet         # Mixed cipher alphabet
        )
        
        # Score the plaintext
        ic = calculate_ic(plaintext)
        bigram_score = score_bigrams(plaintext)
        quadgram_score = score_quadgrams(plaintext)
        
        # Check for known matches
        bigram_match_pct, bigram_matches, total_bigrams = verify_known_bigrams(plaintext)
        quadgram_match_pct, quadgram_matches, total_quadgrams = verify_known_quadgrams(plaintext)
        
        # Combined score (weighted combination of metrics)
        combined_score = (
            0.2 * ic + 
            0.3 * bigram_score +
            0.2 * (bigram_match_pct / 100) + 
            0.3 * (quadgram_match_pct / 100)
        )
        
        # Only record if we have sufficient matches or a good score
        if (len(bigram_matches) >= 3 or len(quadgram_matches) >= 2 or combined_score > 0.4):
            results.append({
                "word": word,
                "primer": primer,
                "primer_name": primer_config["name"],
                "desc": primer_desc,
                "expansion": expansion_method,
                "plaintext": plaintext,
                "ic": ic,
                "bigram_score": bigram_score,
                "quadgram_score": quadgram_score,
                "bigram_matches": len(bigram_matches),
                "bigram_match_pct": bigram_match_pct,
                "quadgram_matches": len(quadgram_matches),
                "quadgram_match_pct": quadgram_match_pct,
                "combined_score": combined_score,
                "match_details": {
                    "bigrams": bigram_matches,
                    "quadgrams": quadgram_matches
                }
            })
        
        # Test reverse configuration (standard as cipher alphabet, mixed as plaintext)
        plaintext = decrypt_gromark(
            K4_CIPHERTEXT, 
            key, 
            keyword_alphabet,       # Mixed plaintext alphabet
            string.ascii_uppercase  # Standard cipher alphabet
        )
        
        # Score the plaintext
        ic = calculate_ic(plaintext)
        bigram_score = score_bigrams(plaintext)
        quadgram_score = score_quadgrams(plaintext)
        
        # Check for known matches
        bigram_match_pct, bigram_matches, total_bigrams = verify_known_bigrams(plaintext)
        quadgram_match_pct, quadgram_matches, total_quadgrams = verify_known_quadgrams(plaintext)
        
        # Combined score
        combined_score = (
            0.2 * ic + 
            0.3 * bigram_score +
            0.2 * (bigram_match_pct / 100) + 
            0.3 * (quadgram_match_pct / 100)
        )
        
        if (len(bigram_matches) >= 3 or len(quadgram_matches) >= 2 or combined_score > 0.4):
            results.append({
                "word": word,
                "primer": primer,
                "primer_name": primer_config["name"],
                "desc": primer_desc,
                "expansion": expansion_method,
                "plaintext": plaintext,
                "ic": ic,
                "bigram_score": bigram_score,
                "quadgram_score": quadgram_score,
                "bigram_matches": len(bigram_matches),
                "bigram_match_pct": bigram_match_pct,
                "quadgram_matches": len(quadgram_matches),
                "quadgram_match_pct": quadgram_match_pct,
                "combined_score": combined_score,
                "match_details": {
                    "bigrams": bigram_matches,
                    "quadgrams": quadgram_matches
                }
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
        print(f"  Expansion: {result['expansion']}")
        print(f"  Combined Score: {result['combined_score']:.4f}")
        print(f"  Bigram Matches: {result['bigram_match_pct']:.2f}% ({result['bigram_matches']} of {len(KNOWN_BIGRAMS)})")
        print(f"  Quadgram Matches: {result['quadgram_match_pct']:.2f}% ({result['quadgram_matches']} of {len(KNOWN_QUADGRAMS)})")
        print(f"  IC: {result['ic']:.6f}")
        print(f"  Bigram Score: {result['bigram_score']:.6f}")
        print(f"  Plaintext: {result['plaintext'][:50]}...")
        
        print("\n  Plaintext at width 21:")
        print_results_at_width21(result['plaintext'])
        
        print("  Matched bigrams:")
        for bigram, positions in result['match_details']['bigrams']:
            print(f"    {bigram} at positions {positions}")
            
        print("  Matched quadgrams:")
        for quadgram, positions in result['match_details']['quadgrams']:
            print(f"    {quadgram} at positions {positions}")
            
        print("-" * 100)

def save_results(results, filename="bigram_quadgram_results.txt"):
    """Save results to a file."""
    with open(filename, 'w', encoding='utf-8') as f:
        f.write(f"K4 Bigram and Quadgram Analysis Results\n")
        f.write(f"Total results: {len(results)}\n")
        f.write("=" * 80 + "\n\n")
        
        for i, result in enumerate(results[:200]):  # Save top 200 results
            f.write(f"Rank {i+1}:\n")
            f.write(f"  Keyword: {result['word']}\n")
            f.write(f"  Primer: {result['primer_name']} {result['primer']} - {result['desc']}\n")
            f.write(f"  Expansion: {result['expansion']}\n")
            f.write(f"  Combined Score: {result['combined_score']:.4f}\n")
            f.write(f"  Bigram Matches: {result['bigram_match_pct']:.2f}% ({result['bigram_matches']} of {len(KNOWN_BIGRAMS)})\n")
            f.write(f"  Quadgram Matches: {result['quadgram_match_pct']:.2f}% ({result['quadgram_matches']} of {len(KNOWN_QUADGRAMS)})\n")
            f.write(f"  IC: {result['ic']:.6f}\n")
            f.write(f"  Bigram Score: {result['bigram_score']:.6f}\n")
            f.write(f"  Quadgram Score: {result['quadgram_score']:.6f}\n")
            f.write(f"  Plaintext: {result['plaintext'][:100]}...\n\n")
            
            f.write("  Matched bigrams:\n")
            for bigram, positions in result['match_details']['bigrams']:
                f.write(f"    {bigram} at positions {positions}\n")
                
            f.write("  Matched quadgrams:\n")
            for quadgram, positions in result['match_details']['quadgrams']:
                f.write(f"    {quadgram} at positions {positions}\n")
            
            f.write("\n")
    
    print(f"Results saved to {filename}")

def analyze_and_sort_positions(results):
    """
    Analyze which positions in the plaintext most commonly match bigrams/quadgrams.
    """
    position_counts = defaultdict(int)
    
    for result in results[:50]:  # Use top 50 results for analysis
        for _, positions in result['match_details']['bigrams']:
            position_counts[positions] += 1
        
        for _, positions in result['match_details']['quadgrams']:
            position_counts[positions] += 1
    
    # Sort by frequency
    sorted_positions = sorted(position_counts.items(), key=lambda x: x[1], reverse=True)
    
    print("\nMost common match positions:")
    for positions, count in sorted_positions[:10]:  # Show top 10
        print(f"  Positions {positions}: {count} occurrences")
    
    return sorted_positions

def main():
    """Main function to analyze K4 using bigram and quadgram matching."""
    parser = argparse.ArgumentParser(description="K4 Bigram and Quadgram Analysis")
    parser.add_argument("--wordlist", default="not.txt", help="Path to wordlist file")
    parser.add_argument("--min-length", type=int, default=4, help="Minimum word length")
    parser.add_argument("--max-length", type=int, default=15, help="Maximum word length")
    parser.add_argument("--processes", type=int, default=None, help="Number of processes to use")
    parser.add_argument("--words", type=str, help="Comma-separated list of specific words to test")
    parser.add_argument("--primers", type=str, help="Comma-separated list of specific primers to test")
    parser.add_argument("--batch-size", type=int, default=100, help="Batch size for processing")
    parser.add_argument("--quadgram-file", type=str, help="File with quadgram statistics")
    parser.add_argument("--output", default="bigram_quadgram_results.txt", help="Output file for results")
    
    args = parser.parse_args()
    
    print("K4 BIGRAM AND QUADGRAM ANALYSIS")
    print("=" * 50)
    
    # Load quadgram statistics if available
    if args.quadgram_file:
        load_quadgram_stats(args.quadgram_file)
    else:
        load_quadgram_stats()  # Use default
    
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
    words = []
    try:
        with open(args.wordlist, 'r', encoding='utf-8') as f:
            for line in f:
                word = line.strip().upper()
                if args.min_length <= len(word) <= args.max_length and all(c in string.ascii_uppercase for c in word):
                    if word_filter is None or word in word_filter:
                        words.append(word)
        print(f"Loaded {len(words)} words from {args.wordlist}")
    except Exception as e:
        print(f"Error loading wordlist: {e}")
    
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
                tasks.append((word, primer_config, expansion_methods))
        
        # Process batch using multiprocessing
        with multiprocessing.Pool(processes=num_processes) as pool:
            batch_results = pool.map(test_keyword_with_primer, tasks)
        
        # Flatten results
        batch_results = [item for sublist in batch_results for item in sublist]
        
        # Add to all results
        all_results.extend(batch_results)
        
        # Sort by combined score
        all_results.sort(key=lambda x: x["combined_score"], reverse=True)
        
        tested_combinations += len(tasks)
        batch_end_time = time.time()
        batch_time = batch_end_time - batch_start_time
        
        # Print batch summary
        print(f"Batch completed in {batch_time:.2f} seconds")
        print(f"Progress: {tested_combinations}/{total_combinations} combinations tested ({tested_combinations/total_combinations*100:.1f}%)")
        print(f"Found {len(batch_results)} notable results in this batch")
        
        # Show best result in this batch
        if batch_results:
            best = max(batch_results, key=lambda x: x["combined_score"])
            print(f"Best in batch: {best['word']} + {best['primer_name']} = Score {best['combined_score']:.4f}")
    
    end_time = time.time()
    print(f"\nAnalysis completed in {end_time - start_time:.2f} seconds")
    print(f"Tested {tested_combinations} keyword-primer combinations")
    print(f"Found {len(all_results)} promising results")
    
    # Print top results
    print_top_results(all_results)
    
    # Analyze position patterns
    analyze_and_sort_positions(all_results)
    
    # Save results
    save_results(all_results, args.output)

if __name__ == "__main__":
    main() 