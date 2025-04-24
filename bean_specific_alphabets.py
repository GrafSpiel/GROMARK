#!/usr/bin/env python3
"""
Kryptos K4 Solver - Bean Specific Alphabet Analysis
Based on Richard Bean's cryptodiagnosis paper

This script focuses specifically on exploring variations of the alphabet offsets
described in Richard Bean's paper, with a focus on the [2,6,7,1,7] primer which
showed the highest match percentage in prior testing.
"""

import string
import time
import argparse
from collections import Counter

# K4 ciphertext (97 letters)
K4_CIPHERTEXT = "OBKRUOXOGHULBSOLIFBBWFLRVQQPRNGKSSOTWTQSJQSSEKZZWATJKLUDIAWINFBNYPVTTMZFPKWGDKZXTJCDIGKUHUAUEKCAR"

# Known plaintext mapping
KNOWN_PLAINTEXT = {
    # position: (plaintext char, ciphertext char)
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

# All observed differences from the known plaintext analysis
OBSERVED_DIFFERENCES = {
    'E': [1, 2, 20],      # E→F, E→G, E→Y
    'A': [10, 11],        # A→K, A→L
    'S': [0, 25],         # S→S, S→R
    'T': [2, 24, 25],     # T→V, T→R, T→S
    'N': [3, 6],          # N→Q, N→T
    'O': [2, 17],         # O→Q, O→F
    'R': [24],            # R→P
    'H': [6],             # H→N
    'B': [12],            # B→N
    'I': [11],            # I→T
    'L': [10, 14],        # L→V, L→Z
    'C': [10, 13],        # C→M, C→P
    'K': [0]              # K→K
}

def expand_key(primer, base, length, expansion_func='standard'):
    """
    Expand a primer to generate a full key sequence.
    
    Args:
        primer: List of digits in the primer
        base: Base of the number system
        length: Desired key length
        expansion_func: 'standard' or 'fibonacci'
    
    Returns:
        List of expanded key digits
    """
    if expansion_func == 'standard':
        # Standard Gromark expansion (add last two digits)
        key = list(primer)
        while len(key) < length:
            next_digit = (key[-1] + key[-2]) % base
            key.append(next_digit)
        return key
    
    elif expansion_func == 'fibonacci':
        # Sum all digits in the primer
        key = list(primer)
        primer_length = len(primer)
        while len(key) < length:
            next_digit = 0
            for i in range(1, primer_length + 1):
                if i <= len(key):
                    next_digit = (next_digit + key[-i]) % base
            key.append(next_digit)
        return key
    
    else:
        raise ValueError(f"Unknown expansion function: {expansion_func}")

def create_bean_specific_alphabet(variation="exact_offsets", offset=0, custom_rules=None):
    """
    Create specific alphabet variations based on observations from Bean's paper.
    
    Args:
        variation: Type of variation to create:
            - "exact_offsets": Use exactly the observed offsets from known plaintext
            - "primary_offsets": Use the most common offset for each letter
            - "consistent_offsets": Use offsets that maintain consistency in the alphabet
            - "positional_offsets": Create offsets based on position in the alphabet
            - "progressive_offsets": Gradually increase/decrease offsets through alphabet
            - "custom": Use custom rules provided in custom_rules parameter
        offset: Additional global offset to apply to the cipher alphabet
        custom_rules: Dictionary of custom mapping rules (only used with "custom" variation)
    
    Returns:
        Tuple of (plaintext_alphabet, ciphertext_alphabet)
    """
    standard = string.ascii_uppercase
    
    if variation == "exact_offsets":
        # Create cipher alphabet with specific offsets from the known plaintext
        # E→F (1 right), A→L (11 right), S→R (25 left or -1), etc.
        plain_to_cipher = {
            'E': 'F',      # +1
            'A': 'L',      # +11
            'S': 'R',      # -1 (or +25)
            'T': 'V',      # +2
            'N': 'Q',      # +3
            'O': 'Q',      # +2
            'R': 'P',      # -2 (or +24)
            'H': 'N',      # +6
            'B': 'N',      # +12
            'I': 'T',      # +11
            'L': 'V',      # +10
            'C': 'M',      # +10
            'K': 'K',      # +0
        }
        
        # Create cipher alphabet by replacing letters according to the mapping
        cipher_list = list(standard)
        for p, c in plain_to_cipher.items():
            p_idx = ord(p) - ord('A')
            c_idx = ord(c) - ord('A')
            cipher_list[p_idx] = standard[c_idx]
        
        # Fill in any unmapped letters with their standard positions
        for i in range(len(cipher_list)):
            if i >= len(standard) or cipher_list[i] not in standard:
                cipher_list[i] = standard[i]
    
    elif variation == "primary_offsets":
        # Create cipher alphabet using the most common offset for each letter
        # Where there are multiple offsets, use the one that occurs most often
        cipher_list = [''] * 26
        
        for p, offsets in OBSERVED_DIFFERENCES.items():
            p_idx = ord(p) - ord('A')
            # Use the most frequent offset (first in the list)
            most_common_offset = offsets[0]
            c_idx = (p_idx + most_common_offset) % 26
            cipher_list[p_idx] = standard[c_idx]
        
        # Fill in any unmapped letters with their standard positions
        for i in range(len(cipher_list)):
            if i < len(standard) and (i >= len(cipher_list) or not cipher_list[i]):
                cipher_list[i] = standard[i]
    
    elif variation == "consistent_offsets":
        # Create a cipher alphabet with offsets that maintain consistency
        # Use either a consistent +1 or +2 shift pattern with a few exceptions
        cipher_list = []
        
        base_offset = 1  # Most common shift seems to be +1 or +2
        
        for p in standard:
            p_idx = ord(p) - ord('A')
            if p in OBSERVED_DIFFERENCES:
                # Use the observed offset if available
                offsets = OBSERVED_DIFFERENCES[p]
                # Try to find an offset close to base_offset
                best_offset = min(offsets, key=lambda x: abs(x - base_offset))
                c_idx = (p_idx + best_offset) % 26
            else:
                # Use the base offset for unmapped letters
                c_idx = (p_idx + base_offset) % 26
            
            cipher_list.append(standard[c_idx])
    
    elif variation == "positional_offsets":
        # Create offsets based on position in the alphabet
        # Alternate between different offsets based on position
        cipher_list = []
        
        for i, p in enumerate(standard):
            position_offset = i % 3  # Alternate between 0, 1, 2
            
            if p in OBSERVED_DIFFERENCES:
                # Try to find an offset consistent with the position pattern
                offsets = OBSERVED_DIFFERENCES[p]
                # Find the offset closest to position_offset
                closest_offset = min(offsets, key=lambda x: abs(x - position_offset) % 26)
                c_idx = (ord(p) - ord('A') + closest_offset) % 26
            else:
                # Use position-based offset for unmapped letters
                c_idx = (ord(p) - ord('A') + position_offset) % 26
            
            cipher_list.append(standard[c_idx])
    
    elif variation == "progressive_offsets":
        # Gradually increase/decrease offsets through the alphabet
        cipher_list = []
        
        for i, p in enumerate(standard):
            # Gradually increase offset from 0 to 10
            progressive_offset = min(i // 2, 10)
            
            if p in OBSERVED_DIFFERENCES:
                # Try to use an observed offset close to the progressive offset
                offsets = OBSERVED_DIFFERENCES[p]
                closest_offset = min(offsets, key=lambda x: abs(x - progressive_offset) % 26)
                c_idx = (ord(p) - ord('A') + closest_offset) % 26
            else:
                # Use progressive offset for unmapped letters
                c_idx = (ord(p) - ord('A') + progressive_offset) % 26
            
            cipher_list.append(standard[c_idx])
    
    elif variation == "custom" and custom_rules:
        # Use custom mapping rules
        cipher_list = list(standard)  # Start with standard alphabet
        
        for p, c in custom_rules.items():
            if p in standard and c in standard:
                p_idx = ord(p) - ord('A')
                c_idx = ord(c) - ord('A')
                cipher_list[p_idx] = standard[c_idx]
    
    else:
        # Default to standard alphabet
        cipher_list = list(standard)
    
    # Create the cipher alphabet
    cipher_alphabet = ''.join(cipher_list)
    
    # Apply additional offset if needed
    if offset != 0:
        cipher_alphabet = ''.join(chr(((ord(c) - ord('A') + offset) % 26) + ord('A')) 
                              for c in cipher_alphabet)
    
    return standard, cipher_alphabet

def decrypt_gromark(ciphertext, key, plain_alphabet, cipher_alphabet):
    """
    Decrypt using the Gromark cipher.
    
    Args:
        ciphertext: The ciphertext to decrypt
        key: List of key digits
        plain_alphabet: The plaintext alphabet
        cipher_alphabet: The ciphertext alphabet
    
    Returns:
        Decrypted plaintext
    """
    # Create mapping from cipher to indices
    cipher_to_indices = {c: i for i, c in enumerate(cipher_alphabet)}
    
    plaintext = []
    for i, c in enumerate(ciphertext):
        if c in cipher_to_indices:
            # Get the position in the cipher alphabet
            pos = cipher_to_indices[c]
            # Subtract the key digit (modulo alphabet length)
            plain_pos = (pos - key[i]) % len(plain_alphabet)
            # Convert back to plaintext
            plaintext.append(plain_alphabet[plain_pos])
        else:
            # Pass through non-alphabet characters
            plaintext.append(c)
    
    return ''.join(plaintext)

def calculate_ic(text):
    """
    Calculate the Index of Coincidence for a text.
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

def evaluate_plaintext(plaintext, known_plaintext=KNOWN_PLAINTEXT):
    """
    Evaluate the quality of a decryption by checking against known plaintext.
    
    Args:
        plaintext: The decrypted plaintext to evaluate
        known_plaintext: Dictionary of known plaintext positions and values
    
    Returns:
        Tuple of (match_percentage, total_matches, total_chars)
    """
    total_matches = 0
    total_chars = 0
    
    for pos, (p, c) in known_plaintext.items():
        if pos < len(plaintext):
            if plaintext[pos] == p:
                total_matches += 1
            total_chars += 1
    
    match_percentage = total_matches / total_chars if total_chars > 0 else 0
    return match_percentage, total_matches, total_chars

def analyze_result_patterns(plaintext):
    """
    Analyze patterns in the decrypted plaintext.
    
    Args:
        plaintext: The decrypted plaintext to analyze
    
    Returns:
        Dictionary with pattern analysis results
    """
    # Common English patterns to check for
    common_bigrams = ['TH', 'HE', 'IN', 'ER', 'AN', 'RE', 'ON', 'AT', 'EN', 'ND']
    common_trigrams = ['THE', 'AND', 'ING', 'ENT', 'ION', 'FOR', 'OUR', 'THA', 'NTH', 'INT']
    common_words = ['THE', 'AND', 'THAT', 'HAVE', 'FOR', 'NOT', 'WITH', 'FROM', 'THIS', 'THEY']
    
    # Count occurrences
    found_bigrams = []
    for i in range(len(plaintext) - 1):
        bigram = plaintext[i:i+2]
        if bigram in common_bigrams:
            found_bigrams.append(bigram)
    
    found_trigrams = []
    for i in range(len(plaintext) - 2):
        trigram = plaintext[i:i+3]
        if trigram in common_trigrams:
            found_trigrams.append(trigram)
    
    found_words = []
    for word in common_words:
        if word in plaintext:
            found_words.append(word)
    
    return {
        'bigrams': found_bigrams,
        'trigrams': found_trigrams,
        'words': found_words
    }

def print_analysis_result(result, rank=None):
    """
    Print a formatted analysis result.
    
    Args:
        result: Dictionary containing analysis result
        rank: Optional rank number to display
    """
    prefix = f"#{rank}: " if rank is not None else ""
    print(f"\n{prefix}Match: {result['match_percentage']:.2%}, IC: {result['ic']:.6f}")
    print(f"Variation: {result['variation']}, Offset: {result['offset']}")
    print(f"Primer: {result['primer']} (Base 10), Expansion: {result['expansion']}")
    print(f"Plaintext: {result['plaintext']}")
    
    # Print plaintext in rows of 21 characters
    print("\nPlaintext at width 21:")
    for i in range(0, len(result['plaintext']), 21):
        print(result['plaintext'][i:i+21])
    
    # Print any interesting patterns found
    patterns = result['patterns']
    if patterns['bigrams'] or patterns['trigrams'] or patterns['words']:
        print("\nInteresting patterns found:")
        if patterns['bigrams']:
            print(f"  Common bigrams: {', '.join(patterns['bigrams'])}")
        if patterns['trigrams']:
            print(f"  Common trigrams: {', '.join(patterns['trigrams'])}")
        if patterns['words']:
            print(f"  Common words: {', '.join(patterns['words'])}")

def analyze_bean_variations():
    """
    Analyze various alphabet variations based on Bean's paper findings.
    Focus on the primer [2,6,7,1,7] which showed the best match.
    
    Returns:
        List of result dictionaries
    """
    print("=== ANALYZING BEAN-SPECIFIC ALPHABET VARIATIONS ===")
    
    # Focus on the primer that gave the best match
    primer = [2, 6, 7, 1, 7]  # Base 10
    base = 10
    
    # Define variations to test
    variations = [
        "exact_offsets",
        "primary_offsets",
        "consistent_offsets",
        "positional_offsets",
        "progressive_offsets"
    ]
    
    # Additional custom variations
    custom_variations = [
        {
            "name": "custom_bean1",
            "rules": {
                'E': 'F',      # +1
                'A': 'L',      # +11
                'S': 'R',      # -1
                'T': 'V',      # +2
                'N': 'Q',      # +3
                'O': 'Q',      # +2
                'R': 'P',      # -2
                'H': 'N',      # +6
                'B': 'N',      # +12
                'I': 'T',      # +11
                'L': 'V',      # +10
                'C': 'M',      # +10
                'K': 'K',      # +0
                # Add consistent +1 shift for remaining letters
                'D': 'E',      # +1
                'F': 'G',      # +1
                'G': 'H',      # +1
                'J': 'K',      # +1
                'M': 'N',      # +1
                'P': 'Q',      # +1
                'Q': 'R',      # +1
                'U': 'V',      # +1
                'V': 'W',      # +1
                'W': 'X',      # +1
                'X': 'Y',      # +1
                'Y': 'Z',      # +1
                'Z': 'A',      # +1
            }
        },
        {
            "name": "custom_bean2",
            "rules": {
                # Try to create a more consistent pattern based on position % 21
                'E': 'F',      # +1 (position 10 % 21 = 10)
                'A': 'L',      # +11 (position 11 % 21 = 11)
                'S': 'R',      # -1 (position 12 % 21 = 12)
                'T': 'V',      # +2 (position 13 % 21 = 13)
                'N': 'Q',      # +3 (position 14 % 21 = 14)
                'O': 'Q',      # +2 (position 15 % 21 = 15)
                'R': 'P',      # -2 (position 16 % 21 = 16)
                'H': 'N',      # +6 (position 18 % 21 = 18)
                # Use consistent +1 for other letters
                'B': 'C',      # +1
                'C': 'D',      # +1
                'D': 'E',      # +1
                'F': 'G',      # +1
                'G': 'H',      # +1
                'I': 'J',      # +1
                'J': 'K',      # +1
                'K': 'L',      # +1
                'L': 'M',      # +1
                'M': 'N',      # +1
                'P': 'Q',      # +1
                'Q': 'R',      # +1
                'U': 'V',      # +1
                'V': 'W',      # +1
                'W': 'X',      # +1
                'X': 'Y',      # +1
                'Y': 'Z',      # +1
                'Z': 'A',      # +1
            }
        }
    ]
    
    # Add custom variations to the list
    for custom in custom_variations:
        variations.append((custom["name"], custom["rules"]))
    
    results = []
    
    # Test each variation with different offsets and key expansion methods
    for variation in variations:
        print(f"\nTesting variation: {variation}")
        
        # Determine if this is a standard variation or custom
        is_custom = isinstance(variation, tuple)
        var_name = variation[0] if is_custom else variation
        custom_rules = variation[1] if is_custom else None
        
        # Try different offsets
        for offset in range(0, 26, 3):  # Try offsets 0, 3, 6, 9, 12, 15, 18, 21, 24
            # Try both standard and Fibonacci key expansion
            for expansion_func in ['standard', 'fibonacci']:
                # Expand the primer
                key = expand_key(primer, base, len(K4_CIPHERTEXT), expansion_func)
                
                # Create the alphabet variation
                if is_custom:
                    plain_alphabet, cipher_alphabet = create_bean_specific_alphabet(
                        variation="custom", offset=offset, custom_rules=custom_rules
                    )
                else:
                    plain_alphabet, cipher_alphabet = create_bean_specific_alphabet(
                        variation=var_name, offset=offset
                    )
                
                # Decrypt
                plaintext = decrypt_gromark(K4_CIPHERTEXT, key, plain_alphabet, cipher_alphabet)
                
                # Evaluate the decryption
                match_percentage, matches, total = evaluate_plaintext(plaintext)
                
                # Calculate Index of Coincidence
                ic = calculate_ic(plaintext)
                
                # Analyze patterns in the plaintext
                patterns = analyze_result_patterns(plaintext)
                
                # Store the result
                result = {
                    'variation': var_name,
                    'offset': offset,
                    'primer': primer,
                    'base': base,
                    'expansion': expansion_func,
                    'plain_alphabet': plain_alphabet,
                    'cipher_alphabet': cipher_alphabet,
                    'plaintext': plaintext,
                    'match_percentage': match_percentage,
                    'matches': matches,
                    'total': total,
                    'ic': ic,
                    'patterns': patterns
                }
                
                results.append(result)
                
                # Print good matches immediately
                if match_percentage >= 0.25:  # Using a lower threshold to catch potentially interesting results
                    print(f"Good match ({match_percentage:.2%}, {matches}/{total}) with offset {offset}, {expansion_func} expansion")
                    print(f"IC: {ic:.6f}")
                    print(f"Plaintext: {plaintext[:50]}...")
    
    # Sort results by match percentage and IC
    results.sort(key=lambda x: (x['match_percentage'], x['ic']), reverse=True)
    
    # Print the best results
    print("\n=== TOP RESULTS ===")
    for i, result in enumerate(results[:10], 1):
        print_analysis_result(result, i)
    
    return results

def main():
    """Main function to run the analysis."""
    parser = argparse.ArgumentParser(
        description="Analyze variations of Bean's approach for Kryptos K4"
    )
    parser.add_argument(
        '--custom-only', action='store_true',
        help='Only test custom variations'
    )
    args = parser.parse_args()
    
    start_time = time.time()
    
    # Run the analysis
    results = analyze_bean_variations()
    
    end_time = time.time()
    print(f"\nTotal execution time: {end_time - start_time:.2f} seconds")

if __name__ == "__main__":
    main() 