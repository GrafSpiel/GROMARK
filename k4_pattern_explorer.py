#!/usr/bin/env python3
"""
K4 Pattern Explorer
Explores promising patterns from the minor differences approach for Kryptos K4.
Focuses on analyzing primer 98800 and patterns like "every 3rd letter differences".
"""

import string
import time
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

# Promising primers
PRIMER_CONFIGS = [
    {"name": "98800", "primer": [9, 8, 8, 0, 0], "base": 10, "desc": "Bean's research - highest IC"},
    {"name": "26717", "primer": [2, 6, 7, 1, 7], "base": 10, "desc": "Bean's research - equivalent primer 1"},
    {"name": "84393", "primer": [8, 4, 3, 9, 3], "base": 10, "desc": "Bean's research - equivalent primer 2"},
    {"name": "1745", "primer": [1, 7, 4, 5], "base": 10, "desc": "Standard Fibonacci"}
]

def expand_key(primer, base, length, method="standard"):
    """Expand a primer to generate a full key sequence."""
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

def create_pattern_alphabet(base_alphabet, pattern, diff):
    """
    Create an alphabet with differences applied according to a pattern.
    
    Args:
        base_alphabet (str): Base alphabet (usually A-Z)
        pattern (str): Pattern descriptor (e.g., 'every3rd', 'every2nd', 'column1')
        diff (int): Difference to apply to the shifted letters
        
    Returns:
        list: Custom alphabets for each position
    """
    alphabets = []
    
    for i in range(len(K4_CIPHERTEXT)):
        if pattern == 'every3rd' and i % 3 == 0:
            # Apply shift to every 3rd letter
            shifted = ''.join(chr(((ord(c) - ord('A') + diff) % 26) + ord('A')) for c in base_alphabet)
            alphabets.append(shifted)
        elif pattern == 'every2nd' and i % 2 == 0:
            # Apply shift to every 2nd letter
            shifted = ''.join(chr(((ord(c) - ord('A') + diff) % 26) + ord('A')) for c in base_alphabet)
            alphabets.append(shifted)
        elif pattern == 'every4th' and i % 4 == 0:
            # Apply shift to every 4th letter
            shifted = ''.join(chr(((ord(c) - ord('A') + diff) % 26) + ord('A')) for c in base_alphabet)
            alphabets.append(shifted)
        elif pattern.startswith('column') and len(pattern) > 6:
            # Apply shift to specific column in width 21 layout
            col = int(pattern[6:])
            if (i % 21) == col:
                shifted = ''.join(chr(((ord(c) - ord('A') + diff) % 26) + ord('A')) for c in base_alphabet)
                alphabets.append(shifted)
            else:
                alphabets.append(base_alphabet)
        else:
            alphabets.append(base_alphabet)
    
    return alphabets

def decrypt_gromark_with_position_alphabets(ciphertext, key, plain_alphabets, cipher_alphabets):
    """
    Decrypt using the Gromark cipher with position-specific alphabets.
    
    Args:
        ciphertext (str): Ciphertext to decrypt
        key (list): Key sequence
        plain_alphabets (list): List of plaintext alphabets for each position
        cipher_alphabets (list): List of ciphertext alphabets for each position
        
    Returns:
        str: Decrypted plaintext
    """
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
    matches = 0
    for pos, (expected, _) in KNOWN_PLAINTEXT.items():
        if pos < len(plaintext) and plaintext[pos] == expected:
            matches += 1
    
    match_percentage = matches / len(KNOWN_PLAINTEXT) * 100
    return match_percentage, matches, len(KNOWN_PLAINTEXT)

def test_every_difference_pattern(primer_config, pattern_type, expansion_method="standard"):
    """
    Test every possible difference (1-25) with a specific pattern.
    
    Args:
        primer_config (dict): Primer configuration
        pattern_type (str): Pattern type to test
        expansion_method (str): Key expansion method
        
    Returns:
        list: Results sorted by match percentage
    """
    results = []
    primer = primer_config["primer"]
    base = primer_config["base"]
    desc = primer_config["desc"]
    
    # Expand the key
    key = expand_key(primer, base, len(K4_CIPHERTEXT), expansion_method)
    
    for diff in range(1, 26):
        # Create position-specific alphabets
        cipher_alphabets = create_pattern_alphabet(string.ascii_uppercase, pattern_type, diff)
        plain_alphabets = [string.ascii_uppercase] * len(K4_CIPHERTEXT)
        
        # Decrypt
        plaintext = decrypt_gromark_with_position_alphabets(K4_CIPHERTEXT, key, plain_alphabets, cipher_alphabets)
        
        # Analyze results
        match_percentage, matches, total = verify_known_plaintext(plaintext)
        ic = calculate_ic(plaintext)
        
        if matches >= 5:  # Only record if at least 5 matches (modify threshold as needed)
            results.append({
                "primer": primer,
                "base": base,
                "desc": desc,
                "pattern": pattern_type,
                "diff": diff,
                "expansion": expansion_method,
                "plaintext": plaintext,
                "match_percentage": match_percentage,
                "matches": matches,
                "total": total,
                "ic": ic
            })
    
    # Sort by match percentage
    results.sort(key=lambda x: (x["matches"], x["ic"]), reverse=True)
    return results

def test_width21_columns(primer_config, expansion_method="standard"):
    """
    Test applying different shifts to each column in width 21 layout.
    
    Args:
        primer_config (dict): Primer configuration
        expansion_method (str): Key expansion method
        
    Returns:
        list: Results sorted by match percentage
    """
    results = []
    primer = primer_config["primer"]
    base = primer_config["base"]
    desc = primer_config["desc"]
    
    # Expand the key
    key = expand_key(primer, base, len(K4_CIPHERTEXT), expansion_method)
    
    for col in range(21):
        for diff in [3, 6, 20]:  # Common differences from minor_differences.py
            pattern = f"column{col}"
            
            # Create position-specific alphabets
            cipher_alphabets = create_pattern_alphabet(string.ascii_uppercase, pattern, diff)
            plain_alphabets = [string.ascii_uppercase] * len(K4_CIPHERTEXT)
            
            # Decrypt
            plaintext = decrypt_gromark_with_position_alphabets(K4_CIPHERTEXT, key, plain_alphabets, cipher_alphabets)
            
            # Analyze results
            match_percentage, matches, total = verify_known_plaintext(plaintext)
            ic = calculate_ic(plaintext)
            
            if matches >= 5:  # Only record if at least 5 matches
                results.append({
                    "primer": primer,
                    "base": base,
                    "desc": desc,
                    "pattern": f"Column {col} shifted by {diff}",
                    "expansion": expansion_method,
                    "plaintext": plaintext,
                    "match_percentage": match_percentage,
                    "matches": matches,
                    "total": total,
                    "ic": ic
                })
    
    # Sort by match percentage
    results.sort(key=lambda x: (x["matches"], x["ic"]), reverse=True)
    return results

def test_hybrid_patterns(primer_config, expansion_method="standard"):
    """
    Test hybrid patterns combining different shifts.
    
    Args:
        primer_config (dict): Primer configuration
        expansion_method (str): Key expansion method
        
    Returns:
        list: Results sorted by match percentage
    """
    results = []
    primer = primer_config["primer"]
    base = primer_config["base"]
    desc = primer_config["desc"]
    
    # Expand the key
    key = expand_key(primer, base, len(K4_CIPHERTEXT), expansion_method)
    
    # Define promising combinations
    combinations = [
        # (Modulo, Remainder, Difference)
        (3, 0, 6),  # Every 3rd letter (remainder 0) shifted by 6
        (3, 1, 3),  # Every 3rd letter (remainder 1) shifted by 3
        (3, 2, 20), # Every 3rd letter (remainder 2) shifted by 20
        (2, 0, 3),  # Every 2nd letter shifted by 3
        (2, 1, 6),  # Every odd position shifted by 6
    ]
    
    for mod, remainder, diff in combinations:
        # Create position-specific alphabets
        cipher_alphabets = []
        
        for i in range(len(K4_CIPHERTEXT)):
            if i % mod == remainder:
                # Apply shift to positions that match the pattern
                shifted = ''.join(chr(((ord(c) - ord('A') + diff) % 26) + ord('A')) for c in string.ascii_uppercase)
                cipher_alphabets.append(shifted)
            else:
                cipher_alphabets.append(string.ascii_uppercase)
        
        plain_alphabets = [string.ascii_uppercase] * len(K4_CIPHERTEXT)
        
        # Decrypt
        plaintext = decrypt_gromark_with_position_alphabets(K4_CIPHERTEXT, key, plain_alphabets, cipher_alphabets)
        
        # Analyze results
        match_percentage, matches, total = verify_known_plaintext(plaintext)
        ic = calculate_ic(plaintext)
        
        results.append({
            "primer": primer,
            "base": base,
            "desc": desc,
            "pattern": f"Mod {mod}, Remainder {remainder}, Shift {diff}",
            "expansion": expansion_method,
            "plaintext": plaintext,
            "match_percentage": match_percentage,
            "matches": matches,
            "total": total,
            "ic": ic
        })
    
    # Sort by match percentage
    results.sort(key=lambda x: (x["matches"], x["ic"]), reverse=True)
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
        print(f"  IC: {result['ic']:.6f}")
        print(f"  Plaintext: {result['plaintext'][:50]}...")
        
        print("\n  Plaintext at width 21:")
        print_results_at_width21(result['plaintext'])
        print("-" * 100)

def main():
    """Main function to explore patterns."""
    print("K4 PATTERN EXPLORER")
    print("=" * 50)
    
    all_results = []
    start_time = time.time()
    
    # Test configurations
    for primer_config in PRIMER_CONFIGS:
        print(f"\nTesting primer {primer_config['name']}: {primer_config['desc']}")
        
        # Test basic patterns
        for pattern in ["every3rd", "every2nd", "every4th"]:
            print(f"  Testing {pattern} pattern...")
            for expansion in ["standard", "fibonacci"]:
                results = test_every_difference_pattern(primer_config, pattern, expansion)
                all_results.extend(results)
        
        # Test width 21 columns (only for the most promising primers)
        if primer_config["name"] in ["98800", "1745"]:
            print(f"  Testing width 21 columns...")
            results = test_width21_columns(primer_config)
            all_results.extend(results)
        
        # Test hybrid patterns
        print(f"  Testing hybrid patterns...")
        results = test_hybrid_patterns(primer_config)
        all_results.extend(results)
    
    # Sort all results
    all_results.sort(key=lambda x: (x["matches"], x["ic"]), reverse=True)
    
    # Print top results
    print_top_results(all_results, top_n=15)
    
    end_time = time.time()
    print(f"\nExecution time: {end_time - start_time:.2f} seconds")

if __name__ == "__main__":
    main() 