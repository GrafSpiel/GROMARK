#!/usr/bin/env python3
"""
K4 Wordlist Explorer
Tests words from a wordlist as keywords for mixed alphabets,
using promising primers and patterns identified in previous analysis.
"""

import string
import time
import itertools
from collections import Counter

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

# Best primers from previous testing
PROMISING_PRIMERS = [
    {"name": "26717", "primer": [2, 6, 7, 1, 7], "base": 10, "desc": "Best performer - Bean's research"},
    {"name": "98800", "primer": [9, 8, 8, 0, 0], "base": 10, "desc": "High IC - Bean's research"},
    {"name": "84393", "primer": [8, 4, 3, 9, 3], "base": 10, "desc": "Alt primer - Bean's research"}
]

# Best patterns from minor difference analysis
PROMISING_PATTERNS = [
    # Position mod 3 patterns with different shifts
    {"type": "mod3", "shifts": {0: 0, 1: 2, 2: 10}, "desc": "Minor differences analysis"},
    {"type": "mod3", "shifts": {0: 6, 1: 3, 2: 20}, "desc": "High IC pattern"},
    {"type": "mod3", "shifts": {0: 3, 1: 6, 2: 20}, "desc": "Alternative pattern"},
    {"type": "mod3", "shifts": {0: 6, 1: 3, 2: 25}, "desc": "Variation with mod 2 = 25"},
    {"type": "mod3", "shifts": {0: 7, 1: 2, 2: 20}, "desc": "Modification of best pattern"}
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

def decrypt_gromark_with_position_alphabets(ciphertext, key, plain_alphabets, cipher_alphabets):
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

def load_wordlist(filename, min_length=3, max_length=15):
    """
    Load words from a wordlist file.
    
    Args:
        filename (str): Path to wordlist file
        min_length (int): Minimum word length to include
        max_length (int): Maximum word length to include
        
    Returns:
        list: Filtered list of words
    """
    words = []
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            for line in f:
                word = line.strip().upper()
                if min_length <= len(word) <= max_length and all(c in string.ascii_uppercase for c in word):
                    words.append(word)
    except Exception as e:
        print(f"Error loading wordlist: {e}")
    
    print(f"Loaded {len(words)} words from {filename}")
    return words

def test_word_with_patterns(word, primer_config, patterns, expansion_method="standard"):
    """
    Test a word as a keyword with different patterns.
    
    Args:
        word (str): Word to use as keyword
        primer_config (dict): Primer configuration
        patterns (list): List of patterns to test
        expansion_method (str): Key expansion method
        
    Returns:
        list: Results sorted by match percentage and IC
    """
    results = []
    primer = primer_config["primer"]
    base = primer_config["base"]
    desc = primer_config["desc"]
    
    # Create the keyword-mixed alphabet
    keyword_alphabet = create_keyword_alphabet(word)
    
    # Expand the key
    key = expand_key(primer, base, len(K4_CIPHERTEXT), expansion_method)
    
    # Test standard alphabet (no pattern)
    # 1. Standard plaintext, keyword ciphertext
    cipher_alphabets = [keyword_alphabet] * len(K4_CIPHERTEXT)
    plain_alphabets = [string.ascii_uppercase] * len(K4_CIPHERTEXT)
    
    plaintext = decrypt_gromark_with_position_alphabets(K4_CIPHERTEXT, key, plain_alphabets, cipher_alphabets)
    match_percentage, matches, total = verify_known_plaintext(plaintext)
    ic = calculate_ic(plaintext)
    
    if len(matches) >= 3:  # Only record if at least 3 matches
        results.append({
            "word": word,
            "primer": primer,
            "base": base,
            "desc": desc,
            "pattern": "Standard keyword alphabet",
            "expansion": expansion_method,
            "plaintext": plaintext,
            "match_percentage": match_percentage,
            "matches": len(matches),
            "matched_positions": [m[0] for m in matches],
            "total": total,
            "ic": ic
        })
    
    # 2. Keyword plaintext, standard ciphertext
    cipher_alphabets = [string.ascii_uppercase] * len(K4_CIPHERTEXT)
    plain_alphabets = [keyword_alphabet] * len(K4_CIPHERTEXT)
    
    plaintext = decrypt_gromark_with_position_alphabets(K4_CIPHERTEXT, key, plain_alphabets, cipher_alphabets)
    match_percentage, matches, total = verify_known_plaintext(plaintext)
    ic = calculate_ic(plaintext)
    
    if len(matches) >= 3:  # Only record if at least 3 matches
        results.append({
            "word": word,
            "primer": primer,
            "base": base,
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
    
    # Test with mod3 patterns
    for pattern in patterns:
        if pattern["type"] == "mod3":
            # Apply pattern to keyword alphabet
            cipher_alphabets = apply_mod3_pattern(keyword_alphabet, pattern["shifts"])
            plain_alphabets = [string.ascii_uppercase] * len(K4_CIPHERTEXT)
            
            plaintext = decrypt_gromark_with_position_alphabets(K4_CIPHERTEXT, key, plain_alphabets, cipher_alphabets)
            match_percentage, matches, total = verify_known_plaintext(plaintext)
            ic = calculate_ic(plaintext)
            
            if len(matches) >= 3:  # Only record if at least 3 matches
                pattern_desc = f"Keyword + {pattern['desc']}"
                results.append({
                    "word": word,
                    "primer": primer,
                    "base": base,
                    "desc": desc,
                    "pattern": pattern_desc,
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
        print(f"  Primer: {result['primer']} ({result['base']}) - {result['desc']}")
        print(f"  Pattern: {result['pattern']}")
        print(f"  Expansion: {result['expansion']}")
        print(f"  Match: {result['match_percentage']:.2f}% ({result['matches']}/{result['total']})")
        print(f"  Matched positions: {result['matched_positions']}")
        print(f"  IC: {result['ic']:.6f}")
        print(f"  Plaintext: {result['plaintext'][:50]}...")
        
        print("\n  Plaintext at width 21:")
        print_results_at_width21(result['plaintext'])
        print("-" * 100)

def save_results(results, filename="wordlist_results.txt"):
    """Save results to a file."""
    with open(filename, 'w', encoding='utf-8') as f:
        f.write(f"K4 Wordlist Explorer Results\n")
        f.write(f"Total results: {len(results)}\n")
        f.write("=" * 80 + "\n\n")
        
        for i, result in enumerate(results[:100]):  # Save top 100 results
            f.write(f"Rank {i+1}:\n")
            f.write(f"  Keyword: {result['word']}\n")
            f.write(f"  Primer: {result['primer']} ({result['base']}) - {result['desc']}\n")
            f.write(f"  Pattern: {result['pattern']}\n")
            f.write(f"  Expansion: {result['expansion']}\n")
            f.write(f"  Match: {result['match_percentage']:.2f}% ({result['matches']}/{result['total']})\n")
            f.write(f"  Matched positions: {result['matched_positions']}\n")
            f.write(f"  IC: {result['ic']:.6f}\n")
            f.write(f"  Plaintext: {result['plaintext'][:100]}...\n\n")
    
    print(f"Results saved to {filename}")

def analyze_keywords(results, top_n=100):
    """Analyze which keywords perform best."""
    top_results = results[:top_n]
    keyword_counts = Counter([r["word"] for r in top_results])
    
    print("\nMost common keywords in top results:")
    for word, count in keyword_counts.most_common(10):
        print(f"  {word}: {count}")
    
    # Group by keyword and find best match for each
    keyword_best = {}
    for result in results:
        word = result["word"]
        if word not in keyword_best or result["matches"] > keyword_best[word]["matches"]:
            keyword_best[word] = result
    
    # Sort by number of matches
    best_keywords = sorted(keyword_best.values(), key=lambda x: (x["matches"], x["ic"]), reverse=True)
    
    print("\nTop 10 keywords by best match:")
    for i, result in enumerate(best_keywords[:10]):
        print(f"  {i+1}. {result['word']}: {result['match_percentage']:.2f}% matches, IC={result['ic']:.6f}")

def main():
    """Main function to test wordlist."""
    print("K4 WORDLIST EXPLORER")
    print("=" * 50)
    
    # Load wordlist
    wordlist_file = "not.txt"
    min_word_length = 4
    max_word_length = 15
    words = load_wordlist(wordlist_file, min_word_length, max_word_length)
    
    # Statistics
    all_results = []
    tested_words = 0
    matches_found = 0
    
    start_time = time.time()
    
    # Use only the best primer and expansion method to speed up testing
    best_primer = PROMISING_PRIMERS[0]  # 26717
    expansion_method = "standard"
    
    print(f"\nTesting with primer {best_primer['name']} and {expansion_method} expansion")
    print(f"Processing {len(words)} words...")
    
    # Process in batches to show progress
    batch_size = 100
    num_batches = (len(words) + batch_size - 1) // batch_size
    
    for batch_idx in range(num_batches):
        batch_start = batch_idx * batch_size
        batch_end = min(batch_start + batch_size, len(words))
        batch_words = words[batch_start:batch_end]
        
        batch_results = []
        for word in batch_words:
            results = test_word_with_patterns(word, best_primer, PROMISING_PATTERNS, expansion_method)
            if results:
                batch_results.extend(results)
                matches_found += 1
            tested_words += 1
        
        all_results.extend(batch_results)
        
        # Sort by matches then IC
        all_results.sort(key=lambda x: (x["matches"], x["ic"]), reverse=True)
        
        # Keep only top results to manage memory
        if len(all_results) > 1000:
            all_results = all_results[:500]
        
        # Show progress
        print(f"Batch {batch_idx+1}/{num_batches}: Processed {tested_words}/{len(words)} words, found {matches_found} with matches")
        # Show top result in each batch
        if batch_results:
            batch_results.sort(key=lambda x: (x["matches"], x["ic"]), reverse=True)
            best = batch_results[0]
            print(f"  Best in batch: {best['word']} - {best['matches']} matches, IC={best['ic']:.6f}")
    
    # If we have enough promising results, test the top keywords with all primers and patterns
    if len(all_results) > 0:
        print("\nRefining top keyword results...")
        
        # Get unique top keywords
        top_keywords = set()
        for result in all_results[:50]:  # Take top 50 results
            top_keywords.add(result["word"])
        
        # Test these keywords more thoroughly
        detailed_results = []
        for word in top_keywords:
            print(f"  Detailed testing of keyword: {word}")
            for primer_config in PROMISING_PRIMERS:
                for expansion in ["standard", "fibonacci"]:
                    results = test_word_with_patterns(word, primer_config, PROMISING_PATTERNS, expansion)
                    detailed_results.extend(results)
        
        # Add to all results
        all_results.extend(detailed_results)
        
        # Sort by matches then IC
        all_results.sort(key=lambda x: (x["matches"], x["ic"]), reverse=True)
    
    end_time = time.time()
    print(f"\nExecution time: {end_time - start_time:.2f} seconds")
    print(f"Processed {tested_words} words, found {matches_found} with 3+ matches")
    
    # Print top results
    print_top_results(all_results, top_n=20)
    
    # Analyze keywords
    analyze_keywords(all_results)
    
    # Save results
    save_results(all_results)

if __name__ == "__main__":
    main() 