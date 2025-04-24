#!/usr/bin/env python3
"""
K4 Enhanced Pattern Explorer V2
An improved version focusing on the most promising primers and patterns,
incorporating minor differences analysis and additional pattern testing.
"""

import string
import time
import itertools
import random
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

# Promising minor differences observed
MINOR_DIFFERENCES = [
    # Position-specific shift patterns inferred from analysis
    {"desc": "Mod 3 with increasing shifts", "shifts": {0: 6, 1: 3, 2: 20}},
    {"desc": "Alternate shift based on position content", "shifts": {}}
]

# Keyword alphabets to test
KEYWORDS = [
    "KRYPTOS", "SANBORN", "BERLIN", "CLOCK", "PALIMPSEST", "ABSCISSA", 
    "NORTHEAST", "BERLIN CLOCK", "BERLIN WALL", "EASTBERLIN",
    "DIGKUHOPFWEXVBNYZRCLSJTQAM",  # Mixed alphabet from minor differences
]

def analyze_minor_differences():
    """Analyze differences between known plaintext and ciphertext."""
    differences = {}
    position_differences = {}
    
    for pos, (plain, cipher) in KNOWN_PLAINTEXT.items():
        # Calculate the shift from plain to cipher
        plain_num = ord(plain) - ord('A')
        cipher_num = ord(cipher) - ord('A')
        shift = (cipher_num - plain_num) % 26
        
        # Store the difference
        differences[pos] = shift
        
        # Group by position mod 3
        pos_mod = pos % 3
        if pos_mod not in position_differences:
            position_differences[pos_mod] = []
        position_differences[pos_mod].append((pos, shift))
    
    # Analyze differences by position mod 3
    mod3_patterns = {}
    for mod, shifts in position_differences.items():
        shift_counts = Counter([s for _, s in shifts])
        most_common = shift_counts.most_common(1)[0][0]
        mod3_patterns[mod] = most_common
    
    return differences, mod3_patterns

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
        elif method == "variable_modulo":
            # Use modulo based on position in the sequence
            if len(key) % 3 == 0:
                next_digit = (key[-1] + key[-2]) % base
            elif len(key) % 3 == 1:
                next_digit = (key[-1] + key[-3]) % base if len(key) >= 3 else (key[-1] + key[-2]) % base
            else:
                next_digit = (key[-2] + key[-3]) % base if len(key) >= 3 else (key[-1] + key[-2]) % base
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

def create_keyword_alphabet(keyword):
    """Create a mixed alphabet based on a keyword."""
    # Remove duplicates while preserving order
    unique_chars = []
    for c in keyword.upper():
        if c in string.ascii_uppercase and c not in unique_chars:
            unique_chars.append(c)
    
    # Create the mixed alphabet
    remaining = [c for c in string.ascii_uppercase if c not in unique_chars]
    return ''.join(unique_chars + remaining)

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
        
        if len(matches) >= 4:  # Only record if at least 4 matches
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
    
    return results

def test_width21_column_patterns(primer_config, expansion_method="standard"):
    """
    Test patterns for a text displayed at width 21.
    
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
    
    # 1. Same shift for all positions in each column
    for shift in [3, 6, 19, 20, 25]:
        shifts = [shift] * 21  # Same shift for all columns
        
        # Create alphabets
        cipher_alphabets = create_width21_column_alphabets(shifts)
        plain_alphabets = [string.ascii_uppercase] * len(K4_CIPHERTEXT)
        
        # Decrypt
        plaintext = decrypt_gromark_with_position_alphabets(K4_CIPHERTEXT, key, plain_alphabets, cipher_alphabets)
        
        # Analyze results
        match_percentage, matches, total = verify_known_plaintext(plaintext)
        ic = calculate_ic(plaintext)
        
        if len(matches) >= 4:  # Only record if at least 4 matches
            pattern_desc = f"Width 21: All columns shift {shift}"
            
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
    
    # Get differences analysis
    differences, mod3_patterns = analyze_minor_differences()
    
    # Apply the most common differences by mod 3
    cipher_alphabets = [string.ascii_uppercase] * len(K4_CIPHERTEXT)
    for i in range(len(K4_CIPHERTEXT)):
        mod = i % 3
        if mod in mod3_patterns:
            shift = mod3_patterns[mod]
            shifted = ''.join(chr(((ord(c) - ord('A') + shift) % 26) + ord('A')) for c in string.ascii_uppercase)
            cipher_alphabets[i] = shifted
    
    plain_alphabets = [string.ascii_uppercase] * len(K4_CIPHERTEXT)
    
    # Decrypt
    plaintext = decrypt_gromark_with_position_alphabets(K4_CIPHERTEXT, key, plain_alphabets, cipher_alphabets)
    
    # Analyze results
    match_percentage, matches, total = verify_known_plaintext(plaintext)
    ic = calculate_ic(plaintext)
    
    # Record the result
    pattern_desc = f"Minor differences: Mod 3 patterns {mod3_patterns}"
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
    
    # Test specific targeted patterns based on grouping in mod 3
    for diff0, diff1, diff2 in [(6, 3, 20), (6, 3, 25), (3, 6, 20), (7, 2, 20)]:
        # Create alphabets
        cipher_alphabets = [string.ascii_uppercase] * len(K4_CIPHERTEXT)
        
        # Apply different shifts based on position mod 3
        for i in range(len(K4_CIPHERTEXT)):
            if i % 3 == 0:
                shifted = ''.join(chr(((ord(c) - ord('A') + diff0) % 26) + ord('A')) for c in string.ascii_uppercase)
                cipher_alphabets[i] = shifted
            elif i % 3 == 1:
                shifted = ''.join(chr(((ord(c) - ord('A') + diff1) % 26) + ord('A')) for c in string.ascii_uppercase)
                cipher_alphabets[i] = shifted
            else:  # i % 3 == 2
                shifted = ''.join(chr(((ord(c) - ord('A') + diff2) % 26) + ord('A')) for c in string.ascii_uppercase)
                cipher_alphabets[i] = shifted
        
        plain_alphabets = [string.ascii_uppercase] * len(K4_CIPHERTEXT)
        
        # Decrypt
        plaintext = decrypt_gromark_with_position_alphabets(K4_CIPHERTEXT, key, plain_alphabets, cipher_alphabets)
        
        # Analyze results
        match_percentage, matches, total = verify_known_plaintext(plaintext)
        ic = calculate_ic(plaintext)
        
        # Record the result
        pattern_desc = f"Mod 3 patterns: 0→{diff0}, 1→{diff1}, 2→{diff2}"
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

def test_keyword_mixed_alphabets(primer_config, expansion_method="standard"):
    """
    Test keyword-mixed alphabets.
    
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
    
    # Test each keyword
    for keyword in KEYWORDS:
        # Create mixed alphabet
        mixed_alphabet = create_keyword_alphabet(keyword)
        
        # Test standard plaintext to mixed ciphertext
        cipher_alphabets = [mixed_alphabet] * len(K4_CIPHERTEXT)
        plain_alphabets = [string.ascii_uppercase] * len(K4_CIPHERTEXT)
        
        # Decrypt
        plaintext = decrypt_gromark_with_position_alphabets(K4_CIPHERTEXT, key, plain_alphabets, cipher_alphabets)
        
        # Analyze results
        match_percentage, matches, total = verify_known_plaintext(plaintext)
        ic = calculate_ic(plaintext)
        
        if len(matches) >= 4:  # Only record if at least 4 matches
            pattern_desc = f"Keyword alphabet: {keyword}"
            
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
        
        # Test mixed plaintext to standard ciphertext
        cipher_alphabets = [string.ascii_uppercase] * len(K4_CIPHERTEXT)
        plain_alphabets = [mixed_alphabet] * len(K4_CIPHERTEXT)
        
        # Decrypt
        plaintext = decrypt_gromark_with_position_alphabets(K4_CIPHERTEXT, key, plain_alphabets, cipher_alphabets)
        
        # Analyze results
        match_percentage, matches, total = verify_known_plaintext(plaintext)
        ic = calculate_ic(plaintext)
        
        if len(matches) >= 4:  # Only record if at least 4 matches
            pattern_desc = f"Mixed plaintext alphabet: {keyword}"
            
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
        
        # For the most promising primer (26717), test combined with mod3 patterns
        if primer_config["name"] == "26717":
            # Get differences analysis
            _, mod3_patterns = analyze_minor_differences()
            
            # Create position-specific alphabets
            cipher_alphabets = []
            for i in range(len(K4_CIPHERTEXT)):
                mod = i % 3
                if mod in mod3_patterns:
                    shift = mod3_patterns[mod]
                    # Shift the mixed alphabet
                    shifted = ''.join(chr(((ord(c) - ord('A') + shift) % 26) + ord('A')) for c in mixed_alphabet)
                    cipher_alphabets.append(shifted)
                else:
                    cipher_alphabets.append(mixed_alphabet)
            
            plain_alphabets = [string.ascii_uppercase] * len(K4_CIPHERTEXT)
            
            # Decrypt
            plaintext = decrypt_gromark_with_position_alphabets(K4_CIPHERTEXT, key, plain_alphabets, cipher_alphabets)
            
            # Analyze results
            match_percentage, matches, total = verify_known_plaintext(plaintext)
            ic = calculate_ic(plaintext)
            
            if len(matches) >= 4:  # Only record if at least 4 matches
                pattern_desc = f"Keyword + Mod3: {keyword}"
                
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

def test_position_specific_offsets(primer_config, expansion_method="standard"):
    """
    Test position-specific offsets based on minor differences analysis.
    
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
    
    # Get the differences
    differences, _ = analyze_minor_differences()
    
    # Create alphabets using the exact differences
    cipher_alphabets = [string.ascii_uppercase] * len(K4_CIPHERTEXT)
    for pos, shift in differences.items():
        shifted = ''.join(chr(((ord(c) - ord('A') + shift) % 26) + ord('A')) for c in string.ascii_uppercase)
        cipher_alphabets[pos] = shifted
    
    plain_alphabets = [string.ascii_uppercase] * len(K4_CIPHERTEXT)
    
    # Decrypt
    plaintext = decrypt_gromark_with_position_alphabets(K4_CIPHERTEXT, key, plain_alphabets, cipher_alphabets)
    
    # Analyze results
    match_percentage, matches, total = verify_known_plaintext(plaintext)
    ic = calculate_ic(plaintext)
    
    pattern_desc = "Exact position-specific differences"
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
    
    # Try extending the pattern to neighboring positions
    extended_cipher_alphabets = [string.ascii_uppercase] * len(K4_CIPHERTEXT)
    for pos, shift in differences.items():
        shifted = ''.join(chr(((ord(c) - ord('A') + shift) % 26) + ord('A')) for c in string.ascii_uppercase)
        extended_cipher_alphabets[pos] = shifted
        
        # Apply to neighboring positions when they don't have a specific shift
        for offset in [-1, 1]:
            neighbor_pos = pos + offset
            if 0 <= neighbor_pos < len(K4_CIPHERTEXT) and neighbor_pos not in differences:
                extended_cipher_alphabets[neighbor_pos] = shifted
    
    # Decrypt with extended pattern
    plaintext = decrypt_gromark_with_position_alphabets(K4_CIPHERTEXT, key, plain_alphabets, extended_cipher_alphabets)
    
    # Analyze results
    match_percentage, matches, total = verify_known_plaintext(plaintext)
    ic = calculate_ic(plaintext)
    
    pattern_desc = "Extended position-specific differences"
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

def print_top_results(results, top_n=15):
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
    print("K4 ENHANCED PATTERN EXPLORER V2")
    print("=" * 50)
    
    all_results = []
    start_time = time.time()
    
    # First, analyze the minor differences
    differences, mod3_patterns = analyze_minor_differences()
    print("\nKNOWN PLAINTEXT ANALYSIS:")
    print(f"Position differences: {differences}")
    print(f"Mod3 patterns: {mod3_patterns}")
    
    # New expansion methods to test
    expansion_methods = ["standard", "fibonacci", "variable_modulo"]
    
    # For each primer configuration
    for primer_config in PRIMER_CONFIGS:
        print(f"\nTesting primer {primer_config['name']}: {primer_config['desc']}")
        
        # Test each expansion method
        for expansion in expansion_methods:
            print(f"  Testing with {expansion} expansion...")
            
            # Test multi-pattern combinations
            results = test_multi_pattern_combinations(primer_config, expansion)
            all_results.extend(results)
            
            # Test width 21 column patterns
            results = test_width21_column_patterns(primer_config, expansion)
            all_results.extend(results)
            
            # Test keyword mixed alphabets
            results = test_keyword_mixed_alphabets(primer_config, expansion)
            all_results.extend(results)
            
            # Test targeted position patterns
            results = test_targeted_position_patterns(primer_config, expansion)
            all_results.extend(results)
            
            # Test position-specific offsets from minor differences
            results = test_position_specific_offsets(primer_config, expansion)
            all_results.extend(results)
    
    # Sort all results
    all_results.sort(key=lambda x: (x["matches"], x["ic"]), reverse=True)
    
    # Print top results
    print_top_results(all_results, top_n=20)
    
    # Analyze which positions are most frequently matched
    match_frequency = analyze_match_positions(all_results)
    print("\nFrequency of matched positions:")
    for pos, count in match_frequency.most_common(10):
        expected, cipher = KNOWN_PLAINTEXT[pos]
        print(f"  Position {pos}: '{expected}' → '{cipher}' - Matched in {count} results")
    
    # Analysis of top results
    top_results = all_results[:10]
    print("\nAnalysis of top results:")
    
    # Count patterns
    pattern_counts = Counter([r["pattern"] for r in top_results])
    print("\nMost common patterns in top results:")
    for pattern, count in pattern_counts.most_common(5):
        print(f"  {pattern}: {count}")
    
    # Count primers
    primer_counts = Counter([r["primer"][0] for r in top_results])
    print("\nMost common primers in top results:")
    for primer, count in primer_counts.most_common(3):
        print(f"  {primer}: {count}")
    
    # Count expansion methods
    expansion_counts = Counter([r["expansion"] for r in top_results])
    print("\nMost common expansion methods in top results:")
    for expansion, count in expansion_counts.most_common(3):
        print(f"  {expansion}: {count}")
    
    end_time = time.time()
    print(f"\nExecution time: {end_time - start_time:.2f} seconds")

if __name__ == "__main__":
    main() 