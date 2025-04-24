#!/usr/bin/env python3
"""
K4 Enhanced Pattern Explorer
Focuses on the most promising primers and patterns discovered previously, 
especially the 26717 primer with mod 3 patterns.
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

# Promising primers from previous exploration
PRIMER_CONFIGS = [
    {"name": "26717", "primer": [2, 6, 7, 1, 7], "base": 10, "desc": "Bean's research - equivalent primer 1"},
    {"name": "98800", "primer": [9, 8, 8, 0, 0], "base": 10, "desc": "Bean's research - highest IC"},
    {"name": "84393", "primer": [8, 4, 3, 9, 3], "base": 10, "desc": "Bean's research - equivalent primer 2"}
]

# Most promising pattern sets from previous exploration
PROMISING_PATTERNS = [
    # (modulo, remainder, shift) combinations that had good results
    [(3, 0, 6), (3, 1, 3), (3, 2, 20)],  # Best combination from previous run
    [(3, 0, 3), (3, 1, 6), (3, 2, 20)],  # Variation
    [(3, 0, 6), (3, 1, 3), (3, 2, 25)],  # Variation
    [(3, 0, 7), (3, 1, 2), (3, 2, 20)],  # Variation
    [(2, 0, 3), (2, 1, 6)],              # Simple even/odd pattern
    [(2, 0, 6), (2, 1, 3)]               # Opposite simple pattern
]

def expand_key(primer, base, length, method="standard"):
    """Expand a primer to generate a key sequence."""
    key = list(primer.copy())
    primer_length = len(primer)
    
    while len(key) < length:
        if method == "standard":
            # Standard Gromark: add last two digits
            next_digit = (key[-1] + key[-2]) % base
            key.append(next_digit)
        elif method == "fibonacci":
            # Sum all digits in the primer
            next_digit = 0
            for i in range(1, primer_length + 1):
                if i <= len(key):
                    next_digit = (next_digit + key[-i]) % base
            key.append(next_digit)
        else:
            # Default to standard
            next_digit = (key[-1] + key[-2]) % base
            key.append(next_digit)
            
    return key[:length]

def create_multi_pattern_alphabets(patterns):
    """
    Create alphabets based on multiple patterns.
    
    Args:
        patterns: List of (modulo, remainder, shift) tuples
        
    Returns:
        List of position-specific alphabets
    """
    alphabets = [string.ascii_uppercase] * len(K4_CIPHERTEXT)
    
    for mod, remainder, shift in patterns:
        for i in range(len(K4_CIPHERTEXT)):
            if i % mod == remainder:
                # Apply shift to positions that match the pattern
                alphabet = alphabets[i]
                alphabets[i] = ''.join(chr(((ord(c) - ord('A') + shift) % 26) + ord('A')) for c in alphabet)
    
    return alphabets

def create_width21_column_alphabets(shifts):
    """
    Create alphabets based on width 21 column patterns.
    
    Args:
        shifts: List of shifts for each column (length=21)
        
    Returns:
        List of position-specific alphabets
    """
    alphabets = []
    
    for i in range(len(K4_CIPHERTEXT)):
        col = i % 21
        shift = shifts[col]
        
        if shift == 0:
            alphabets.append(string.ascii_uppercase)
        else:
            shifted = ''.join(chr(((ord(c) - ord('A') + shift) % 26) + ord('A')) for c in string.ascii_uppercase)
            alphabets.append(shifted)
    
    return alphabets

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

def test_multi_pattern_combinations(primer_config, expansion_method="standard"):
    """
    Test combinations of different patterns applied together.
    
    Args:
        primer_config: Primer configuration
        expansion_method: Key expansion method
        
    Returns:
        List of results
    """
    results = []
    primer = primer_config["primer"]
    base = primer_config["base"]
    desc = primer_config["desc"]
    
    # Expand the key
    key = expand_key(primer, base, len(K4_CIPHERTEXT), expansion_method)
    
    # Test predefined promising pattern combinations
    for pattern_set in PROMISING_PATTERNS:
        # Create alphabets based on the pattern set
        cipher_alphabets = create_multi_pattern_alphabets(pattern_set)
        plain_alphabets = [string.ascii_uppercase] * len(K4_CIPHERTEXT)
        
        # Decrypt
        plaintext = decrypt_gromark_with_position_alphabets(K4_CIPHERTEXT, key, plain_alphabets, cipher_alphabets)
        
        # Analyze results
        match_percentage, matches, total = verify_known_plaintext(plaintext)
        ic = calculate_ic(plaintext)
        
        pattern_desc = " + ".join([f"(mod {m}, rem {r}, shift {s})" for m, r, s in pattern_set])
        
        results.append({
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
    
    # Test custom combinations of the most successful patterns
    successful_diffs = [3, 6, 20]
    for diff1, diff2, diff3 in itertools.product(successful_diffs, repeat=3):
        pattern_set = [(3, 0, diff1), (3, 1, diff2), (3, 2, diff3)]
        
        # Create alphabets based on the pattern set
        cipher_alphabets = create_multi_pattern_alphabets(pattern_set)
        plain_alphabets = [string.ascii_uppercase] * len(K4_CIPHERTEXT)
        
        # Decrypt
        plaintext = decrypt_gromark_with_position_alphabets(K4_CIPHERTEXT, key, plain_alphabets, cipher_alphabets)
        
        # Analyze results
        match_percentage, matches, total = verify_known_plaintext(plaintext)
        ic = calculate_ic(plaintext)
        
        if len(matches) >= 4:  # Only record if at least 4 matches
            pattern_desc = f"(mod 3, rem 0, shift {diff1}) + (mod 3, rem 1, shift {diff2}) + (mod 3, rem 2, shift {diff3})"
            
            results.append({
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

def test_width21_column_patterns(primer_config, expansion_method="standard"):
    """
    Test patterns based on width 21 columns.
    
    Args:
        primer_config: Primer configuration
        expansion_method: Key expansion method
        
    Returns:
        List of results
    """
    results = []
    primer = primer_config["primer"]
    base = primer_config["base"]
    desc = primer_config["desc"]
    
    # Expand the key
    key = expand_key(primer, base, len(K4_CIPHERTEXT), expansion_method)
    
    # Create different column shift patterns
    # 1. Alternating pattern: shifts for even and odd columns
    for even_shift, odd_shift in [(3, 6), (6, 3), (3, 20), (20, 3)]:
        shifts = [even_shift if col % 2 == 0 else odd_shift for col in range(21)]
        
        # Create alphabets
        cipher_alphabets = create_width21_column_alphabets(shifts)
        plain_alphabets = [string.ascii_uppercase] * len(K4_CIPHERTEXT)
        
        # Decrypt
        plaintext = decrypt_gromark_with_position_alphabets(K4_CIPHERTEXT, key, plain_alphabets, cipher_alphabets)
        
        # Analyze results
        match_percentage, matches, total = verify_known_plaintext(plaintext)
        ic = calculate_ic(plaintext)
        
        if len(matches) >= 4:  # Only record if at least 4 matches
            pattern_desc = f"Width 21: Even cols shift {even_shift}, Odd cols shift {odd_shift}"
            
            results.append({
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
    
    # 2. Three-pattern cycle for columns
    for shift1, shift2, shift3 in [(3, 6, 20), (6, 3, 20), (20, 3, 6), (20, 6, 3)]:
        shifts = [shift1 if col % 3 == 0 else shift2 if col % 3 == 1 else shift3 for col in range(21)]
        
        # Create alphabets
        cipher_alphabets = create_width21_column_alphabets(shifts)
        plain_alphabets = [string.ascii_uppercase] * len(K4_CIPHERTEXT)
        
        # Decrypt
        plaintext = decrypt_gromark_with_position_alphabets(K4_CIPHERTEXT, key, plain_alphabets, cipher_alphabets)
        
        # Analyze results
        match_percentage, matches, total = verify_known_plaintext(plaintext)
        ic = calculate_ic(plaintext)
        
        if len(matches) >= 4:  # Only record if at least 4 matches
            pattern_desc = f"Width 21: Cols mod 3 = 0: shift {shift1}, mod 3 = 1: shift {shift2}, mod 3 = 2: shift {shift3}"
            
            results.append({
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

def test_targeted_position_patterns(primer_config, expansion_method="standard"):
    """
    Test patterns specifically targeted at known plaintext positions.
    
    Args:
        primer_config: Primer configuration
        expansion_method: Key expansion method
        
    Returns:
        List of results
    """
    results = []
    primer = primer_config["primer"]
    base = primer_config["base"]
    desc = primer_config["desc"]
    
    # Expand the key
    key = expand_key(primer, base, len(K4_CIPHERTEXT), expansion_method)
    
    # Group known plaintext by remainder mod 3
    by_mod3 = {0: [], 1: [], 2: []}
    for pos, (plain, cipher) in KNOWN_PLAINTEXT.items():
        by_mod3[pos % 3].append((pos, plain, cipher))
    
    # Try different shifts for each remainder group
    for diff0, diff1, diff2 in itertools.product(range(1, 26), range(1, 26), range(1, 26)):
        # Create alphabets
        cipher_alphabets = [string.ascii_uppercase] * len(K4_CIPHERTEXT)
        
        # Apply different shifts to positions based on their mod 3 remainder
        for pos, _, _ in by_mod3[0]:
            shifted = ''.join(chr(((ord(c) - ord('A') + diff0) % 26) + ord('A')) for c in string.ascii_uppercase)
            cipher_alphabets[pos] = shifted
        
        for pos, _, _ in by_mod3[1]:
            shifted = ''.join(chr(((ord(c) - ord('A') + diff1) % 26) + ord('A')) for c in string.ascii_uppercase)
            cipher_alphabets[pos] = shifted
        
        for pos, _, _ in by_mod3[2]:
            shifted = ''.join(chr(((ord(c) - ord('A') + diff2) % 26) + ord('A')) for c in string.ascii_uppercase)
            cipher_alphabets[pos] = shifted
        
        plain_alphabets = [string.ascii_uppercase] * len(K4_CIPHERTEXT)
        
        # Decrypt
        plaintext = decrypt_gromark_with_position_alphabets(K4_CIPHERTEXT, key, plain_alphabets, cipher_alphabets)
        
        # Analyze results
        match_percentage, matches, total = verify_known_plaintext(plaintext)
        
        # Only record if we achieve at least 50% match
        if len(matches) >= 12:  # 12 = 50% of 24
            ic = calculate_ic(plaintext)
            pattern_desc = f"Targeted: mod3=0 shift {diff0}, mod3=1 shift {diff1}, mod3=2 shift {diff2}"
            
            results.append({
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

def print_top_results(results, top_n=10):
    """Print top results."""
    print("\nTOP RESULTS")
    print("=" * 100)
    
    for i, result in enumerate(results[:top_n]):
        print(f"\nRank {i+1}:")
        print(f"  Primer: {result['primer']} ({result['base']})")
        print(f"  Description: {result['desc']}")
        print(f"  Pattern: {result['pattern']}")
        print(f"  Expansion: {result['expansion']}")
        print(f"  Match: {result['match_percentage']:.2f}% ({result['matches']}/{result['total']})")
        print(f"  Matched positions: {result['matched_positions']}")
        print(f"  IC: {result['ic']:.6f}")
        print(f"  Plaintext: {result['plaintext'][:50]}...")
        
        print("\n  Plaintext at width 21:")
        print_results_at_width21(result['plaintext'])
        print("-" * 100)

def analyze_match_positions(results, min_matches=4):
    """
    Analyze which positions are most frequently matched across results.
    
    Args:
        results: List of results
        min_matches: Minimum number of matches to include in analysis
        
    Returns:
        Dict mapping position to match frequency
    """
    position_matches = Counter()
    qualifying_results = [r for r in results if r["matches"] >= min_matches]
    
    for result in qualifying_results:
        for pos in result["matched_positions"]:
            position_matches[pos] += 1
    
    return position_matches

def main():
    """Main function to explore patterns."""
    print("K4 ENHANCED PATTERN EXPLORER")
    print("=" * 50)
    
    all_results = []
    start_time = time.time()
    
    # For the enhanced explorer, focus primarily on the 26717 primer
    for primer_config in PRIMER_CONFIGS:
        print(f"\nTesting primer {primer_config['name']}: {primer_config['desc']}")
        
        # Test multi-pattern combinations
        print(f"  Testing multi-pattern combinations...")
        for expansion in ["standard", "fibonacci"]:
            results = test_multi_pattern_combinations(primer_config, expansion)
            all_results.extend(results)
        
        # Test width 21 column patterns
        print(f"  Testing width 21 column patterns...")
        results = test_width21_column_patterns(primer_config)
        all_results.extend(results)
        
        # Only do targeted position testing for most promising primer (26717)
        if primer_config["name"] == "26717":
            print(f"  Testing targeted position patterns...")
            results = test_targeted_position_patterns(primer_config)
            all_results.extend(results)
    
    # Sort all results
    all_results.sort(key=lambda x: (x["matches"], x["ic"]), reverse=True)
    
    # Print top results
    print_top_results(all_results, top_n=15)
    
    # Analyze which positions are most frequently matched
    match_frequency = analyze_match_positions(all_results)
    print("\nFrequency of matched positions:")
    for pos, count in match_frequency.most_common(10):
        expected, cipher = KNOWN_PLAINTEXT[pos]
        print(f"  Position {pos}: '{expected}' → '{cipher}' - Matched in {count} results")
    
    end_time = time.time()
    print(f"\nExecution time: {end_time - start_time:.2f} seconds")

if __name__ == "__main__":
    main() 