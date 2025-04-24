#!/usr/bin/env python3
"""
k4_minor_differences.py - Tests minor differences approach for solving Kryptos K4
Based on Richard Bean's cryptodiagnosis paper

This script implements the "minor differences" approach - where the plaintext and 
ciphertext alphabets differ by small shifts or substitutions.
"""

import string
import itertools
from collections import Counter
import time

# K4 ciphertext 
K4_CIPHERTEXT = "OBKRUOXOGHULBSOLIFBBWFLRVQQPRNGKSSOTWTQSJQSSEKZZWATJKLUDIAWINFBNYPVTTMZFPKWGDKZXTJCDIGKUHUAUEKCAR"

# Known plaintext at specific positions (based on previous solvers' consensus)
KNOWN_PLAINTEXT = {
    0: ('O', 'B'),  # Position 0: Plaintext 'B', Ciphertext 'O'
    1: ('B', 'E'),  # Position 1: Plaintext 'E', Ciphertext 'B'
    2: ('K', 'R'),  # And so on...
    3: ('R', 'L'),
    4: ('U', 'I'),
    5: ('O', 'N'),
    70: ('N', 'T'),
    71: ('F', 'H'),
    72: ('B', 'E'),
    73: ('N', 'E'),
    74: ('Y', 'A'),
    75: ('P', 'S'),
    76: ('V', 'T'),
    77: ('T', 'W'),
    78: ('T', 'E'),
    79: ('M', 'S'),
    80: ('Z', 'T'),
    93: ('H', 'I'),
    94: ('U', 'D'),
    95: ('A', 'Y'),
    96: ('U', 'O'),
    97: ('E', 'U')
}

# Promising primers from previous research
PROMISING_PRIMERS = [
    (26717, 10, "Bean's research"),
    (84393, 10, "Bean's research"), 
    (98800, 10, "Bean's research"),
    (3, 10, "Sanborn hint: base 10, single digit"),
    (7, 10, "Sanborn hint: base 10, single digit"),
    (4, 10, "Sanborn hint: base 10, single digit"),
    (1745, 10, "Standard Fibonacci"),
    (1123, 10, "Lucas sequence"),
]

def analyze_minor_differences():
    """Analyze the minor differences between known plaintext and ciphertext."""
    differences = []
    for pos, (plain, cipher) in KNOWN_PLAINTEXT.items():
        plain_idx = string.ascii_uppercase.index(plain)
        cipher_idx = string.ascii_uppercase.index(cipher)
        diff = (cipher_idx - plain_idx) % 26
        differences.append((pos, plain, cipher, diff))
    
    return differences

def expand_key(primer, length, base=10):
    """Expand a primer using Gromark-style Fibonacci expansion."""
    key = [primer]
    for i in range(1, length):
        next_value = (key[i-1] + primer) % base
        key.append(next_value)
    return key

def create_shifted_alphabet(base_alphabet, shift):
    """Create an alphabet shifted by a specific amount."""
    shifted = ''
    for c in base_alphabet:
        new_pos = (string.ascii_uppercase.index(c) + shift) % 26
        shifted += string.ascii_uppercase[new_pos]
    return shifted

def create_keyword_alphabet(keyword, reverse=False):
    """Create a mixed alphabet using a keyword."""
    # Remove duplicates while preserving order
    unique_chars = ""
    for c in keyword.upper():
        if c not in unique_chars and c in string.ascii_uppercase:
            unique_chars += c
    
    # Add remaining letters
    for c in string.ascii_uppercase:
        if c not in unique_chars:
            unique_chars += c
    
    if reverse:
        return unique_chars[::-1]
    return unique_chars

def decrypt_gromark(ciphertext, key, plaintext_alphabet, ciphertext_alphabet):
    """Decrypt using Gromark cipher with possibly different alphabets."""
    plaintext = ""
    for i, c in enumerate(ciphertext):
        if c not in ciphertext_alphabet:
            plaintext += c
            continue
        
        pos = ciphertext_alphabet.index(c)
        shift = key[i % len(key)]
        new_pos = (pos - shift) % 26
        plaintext += plaintext_alphabet[new_pos]
    
    return plaintext

def calculate_ic(text):
    """Calculate the Index of Coincidence for a text."""
    text = ''.join(c for c in text.upper() if c in string.ascii_uppercase)
    if len(text) <= 1:
        return 0.0
    
    counts = Counter(text)
    n = len(text)
    sum_fi_choose_2 = sum(count * (count - 1) for count in counts.values())
    ic = sum_fi_choose_2 / (n * (n - 1))
    return ic

def verify_known_plaintext(decrypted):
    """Check how many known plaintext characters match."""
    matches = 0
    mismatches = []
    
    for pos, (expected, _) in KNOWN_PLAINTEXT.items():
        if pos < len(decrypted) and decrypted[pos] == expected:
            matches += 1
        else:
            mismatches.append(pos)
    
    match_percentage = matches / len(KNOWN_PLAINTEXT) * 100
    return match_percentage, matches, mismatches

def main():
    print("KRYPTOS K4 - MINOR DIFFERENCES ANALYSIS")
    print("-" * 50)
    
    # Analyze the known differences
    differences = analyze_minor_differences()
    
    print("Known plaintext-ciphertext pairs and their differences:")
    for pos, plain, cipher, diff in differences:
        print(f"Position {pos}: Plaintext '{plain}' → Ciphertext '{cipher}' (Difference: {diff})")
    
    # Calculate pattern statistics
    diffs = [d[3] for d in differences]
    counter = Counter(diffs)
    print("\nDifference frequencies:")
    for diff, count in sorted(counter.items()):
        print(f"Difference {diff}: {count} occurrences")
    
    print("\n" + "=" * 50)
    print("TESTING SHIFTED ALPHABETS")
    print("=" * 50)
    
    results = []
    
    # Test with standard alphabet shifted by various amounts
    for shift in range(26):
        plaintext_alphabet = string.ascii_uppercase
        ciphertext_alphabet = create_shifted_alphabet(plaintext_alphabet, shift)
        
        for primer_val, base, desc in PROMISING_PRIMERS:
            key = expand_key(primer_val, len(K4_CIPHERTEXT), base)
            decrypted = decrypt_gromark(K4_CIPHERTEXT, key, plaintext_alphabet, ciphertext_alphabet)
            ic = calculate_ic(decrypted)
            match_percentage, matches, mismatches = verify_known_plaintext(decrypted)
            
            results.append({
                'type': 'shifted',
                'shift': shift,
                'primer': primer_val,
                'base': base,
                'description': desc,
                'match_percentage': match_percentage,
                'matches': matches,
                'mismatches': mismatches,
                'ic': ic,
                'decrypted': decrypted
            })
    
    print("\n" + "=" * 50)
    print("TESTING KEYWORD-BASED ALPHABETS")
    print("=" * 50)
    
    keywords = [
        "KRYPTOS", "BERLIN", "SANBORN", "LANGLEY", "PALIMPSEST",
        "ABSCISSA", "IQLUSION", "SHADOW", "DIGGERSPADE", "UNDERGRUUND"
    ]
    
    for keyword in keywords:
        plaintext_alphabet = string.ascii_uppercase
        ciphertext_alphabet = create_keyword_alphabet(keyword)
        
        for primer_val, base, desc in PROMISING_PRIMERS:
            key = expand_key(primer_val, len(K4_CIPHERTEXT), base)
            decrypted = decrypt_gromark(K4_CIPHERTEXT, key, plaintext_alphabet, ciphertext_alphabet)
            ic = calculate_ic(decrypted)
            match_percentage, matches, mismatches = verify_known_plaintext(decrypted)
            
            results.append({
                'type': 'keyword',
                'keyword': keyword,
                'primer': primer_val,
                'base': base,
                'description': desc,
                'match_percentage': match_percentage,
                'matches': matches,
                'mismatches': mismatches,
                'ic': ic,
                'decrypted': decrypted
            })
    
    print("\n" + "=" * 50)
    print("TESTING CUSTOM APPROACH - MINOR DIFFERENCES")
    print("=" * 50)
    
    # Based on the patterns observed in differences, create custom alphabets
    # For example, we might try creating alphabets where only specific positions
    # have differences that match the observed patterns
    
    # Example: Create a custom mixed alphabet based on most frequent differences
    common_diffs = [diff for diff, _ in counter.most_common(3)]
    
    alphabet = string.ascii_uppercase
    for diff in common_diffs:
        custom_alphabet = ""
        for i, c in enumerate(alphabet):
            shift = diff if i % 3 == 0 else 0  # Apply shift only to every third letter
            new_pos = (string.ascii_uppercase.index(c) + shift) % 26
            custom_alphabet += string.ascii_uppercase[new_pos]
        
        for primer_val, base, desc in PROMISING_PRIMERS:
            key = expand_key(primer_val, len(K4_CIPHERTEXT), base)
            decrypted = decrypt_gromark(K4_CIPHERTEXT, key, alphabet, custom_alphabet)
            ic = calculate_ic(decrypted)
            match_percentage, matches, mismatches = verify_known_plaintext(decrypted)
            
            results.append({
                'type': 'custom',
                'diff_pattern': f"Common diff {diff} applied to every 3rd letter",
                'primer': primer_val,
                'base': base,
                'description': desc,
                'match_percentage': match_percentage,
                'matches': matches,
                'mismatches': mismatches,
                'ic': ic,
                'decrypted': decrypted
            })
    
    # Sort results by match percentage and IC
    sorted_results = sorted(results, key=lambda x: (x['match_percentage'], x['ic']), reverse=True)
    
    # Display top results
    print("\nTop Results:")
    print("-" * 100)
    for i, result in enumerate(sorted_results[:10]):
        print(f"Rank {i+1}:")
        print(f"  Type: {result['type']}")
        if result['type'] == 'shifted':
            print(f"  Shift: {result['shift']}")
        elif result['type'] == 'keyword':
            print(f"  Keyword: {result['keyword']}")
        elif result['type'] == 'custom':
            print(f"  Pattern: {result['diff_pattern']}")
        
        print(f"  Primer: {result['primer']} (Base {result['base']}) - {result['description']}")
        print(f"  Match: {result['match_percentage']:.2f}% ({result['matches']}/{len(KNOWN_PLAINTEXT)})")
        print(f"  IC: {result['ic']:.6f}")
        print(f"  Decrypted: {result['decrypted'][:30]}...")
        print("-" * 100)

if __name__ == "__main__":
    start_time = time.time()
    main()
    print(f"\nExecution time: {time.time() - start_time:.2f} seconds") 