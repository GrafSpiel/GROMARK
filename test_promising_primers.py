#!/usr/bin/env python3
"""
Kryptos K4 Promising Primers Tester
Tests the most promising primers identified in Richard Bean's cryptodiagnosis paper.
"""

import string
import itertools
from collections import Counter
import time

# K4 ciphertext (97 letters)
K4_CIPHERTEXT = "OBKRUOXOGHULBSOLIFBBWFLRVQQPRNGKSSOTWTQSJQSSEKZZWATJKLUDIAWINFBNYPVTTMZFPKWGDKZXTJCDIGKUHUAUEKCAR"

# Known plaintext segments with their positions (0-indexed)
KNOWN_PLAINTEXT = {
    10: "EASTNORTHEAST",
    62: "BERLINCLOCK"
}

# Most promising primers from the paper
PROMISING_PRIMERS = [
    # Base 10, length 5
    (10, 5, [2, 6, 7, 1, 7], "Bean's paper - identified as one of two equivalent primers"),
    (10, 5, [8, 4, 3, 9, 3], "Bean's paper - identified as one of two equivalent primers"),
    (10, 5, [9, 8, 8, 0, 0], "Bean's paper - gives highest IC of 0.0625"),
    
    # Base 8, length 5
    (8, 5, [0, 0, 3, 5, 1], "Base 8 primer with period 84"),
    (8, 5, [0, 0, 5, 3, 7], "Base 8 primer with period 84"),
    
    # Base 10, length 4
    (10, 4, [3, 3, 0, 1], "Base 10, length 4 primer"),
    (10, 4, [6, 7, 4, 0], "Base 10, length 4 primer"),
    (10, 4, [9, 9, 0, 3], "Base 10, length 4 primer")
]

# Potential keywords for cipher alphabet mixing
KEYWORDS = [
    "KRYPTOS", "PALIMPSEST", "ABSCISSA", "SHADOW", "BERLIN", 
    "CLOCK", "LANGLEY", "SANBORN", "SCHEIDT", "NYPVTTMZFPK",
    "FLRVQQPRNGKSS", "CIA", "NORTHEAST"
]

# Key expansion functions
def standard_expansion(primer, base, length):
    """Standard Gromark expansion (add last two digits)"""
    key = list(primer)
    while len(key) < length:
        next_digit = (key[-1] + key[-2]) % base
        key.append(next_digit)
    return key

def fibonacci_expansion(primer, base, length):
    """Sum all digits in the primer"""
    key = list(primer)
    primer_length = len(primer)
    while len(key) < length:
        next_digit = 0
        for i in range(1, primer_length + 1):
            if i <= len(key):
                next_digit = (next_digit + key[-i]) % base
        key.append(next_digit)
    return key

def alternating_expansion(primer, base, length):
    """Alternating addition and subtraction"""
    key = list(primer)
    while len(key) < length:
        if len(key) % 2 == 0:
            next_digit = (key[-1] + key[-2]) % base
        else:
            next_digit = (key[-1] - key[-2]) % base
            if next_digit < 0:
                next_digit += base
        key.append(next_digit)
    return key

def weighted_expansion(primer, base, length):
    """Weighted expansion - double the last digit plus previous digit"""
    key = list(primer)
    while len(key) < length:
        next_digit = ((2 * key[-1]) + key[-2]) % base
        key.append(next_digit)
    return key

def create_keyword_alphabet(keyword, reverse=False):
    """Create a keyword-mixed alphabet."""
    # Convert to uppercase and remove duplicates while preserving order
    keyword = keyword.upper()
    seen = set()
    unique_keyword = ''.join(c for c in keyword if c in string.ascii_uppercase 
                          and c not in seen and not seen.add(c))
    
    # Remove keyword letters from the remaining alphabet
    remaining = ''.join(c for c in string.ascii_uppercase if c not in unique_keyword)
    
    # Create the mixed alphabet
    mixed_alphabet = unique_keyword + remaining
    
    # Reverse if requested
    if reverse:
        mixed_alphabet = mixed_alphabet[::-1]
    
    return mixed_alphabet

def decrypt_gromark(ciphertext, key, plain_alphabet, cipher_alphabet):
    """Decrypt using the Gromark cipher."""
    # Create mapping from cipher to plain
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
    """Calculate the Index of Coincidence for a text."""
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

def check_decryption_quality(plaintext):
    """Evaluate the quality of decryption."""
    # Count how many valid English letters
    english_letters = sum(1 for c in plaintext if c in string.ascii_uppercase)
    english_ratio = english_letters / len(plaintext)
    
    # Calculate IC (closer to 0.067 is more like English)
    ic = calculate_ic(plaintext)
    ic_quality = 1 - abs(ic - 0.067) / 0.067
    
    # Simple word check - look for common English words
    common_words = ["THE", "AND", "FOR", "THAT", "WITH", "HAVE", "THIS", "FROM", "THEY", "WILL"]
    word_score = sum(1 for word in common_words if word in plaintext)
    
    # Check if known plaintext segments are reasonably present
    known_match = False
    for pos, known_text in KNOWN_PLAINTEXT.items():
        segment = plaintext[pos:pos+len(known_text)]
        match_percentage = sum(1 for a, b in zip(known_text, segment) if a == b) / len(known_text)
        if match_percentage > 0.7:  # 70% match threshold
            known_match = True
            break
    
    # Combined quality score
    quality = (0.4 * english_ratio) + (0.4 * ic_quality) + (0.2 * (word_score / 10))
    
    # If known plaintext is not found, quality is very low
    if not known_match:
        quality *= 0.1
    
    return {
        "english_ratio": english_ratio,
        "ic": ic,
        "ic_quality": ic_quality,
        "word_score": word_score,
        "known_match": known_match,
        "quality": quality
    }

def test_primer_configuration(base, primer, expansion_func, plain_alphabet, cipher_alphabet):
    """Test a specific primer configuration and return details about the result."""
    # Generate the key
    key = expansion_func(primer, base, len(K4_CIPHERTEXT))
    
    # Decrypt the ciphertext
    plaintext = decrypt_gromark(K4_CIPHERTEXT, key, plain_alphabet, cipher_alphabet)
    
    # Evaluate the quality
    quality = check_decryption_quality(plaintext)
    
    return {
        "primer": primer,
        "plaintext": plaintext,
        "key": key[:10],  # First 10 key digits
        "quality": quality
    }

def main():
    """Test the most promising primers with various configurations."""
    results = []
    
    # Expansion functions to try
    expansion_functions = [
        ("standard", standard_expansion),
        ("fibonacci", fibonacci_expansion),
        ("alternating", alternating_expansion),
        ("weighted", weighted_expansion)
    ]
    
    print("Testing promising primers from Bean's paper...")
    
    # Track start time
    start_time = time.time()
    
    # Try each primer with various configurations
    for base, length, primer, description in PROMISING_PRIMERS:
        print(f"\nTesting primer {primer} (Base {base}, Length {length}): {description}")
        
        # Try each expansion function
        for exp_name, exp_func in expansion_functions:
            print(f"  Using {exp_name} expansion...")
            
            # Try standard alphabets first
            result = test_primer_configuration(
                base, primer, exp_func, 
                string.ascii_uppercase, string.ascii_uppercase
            )
            
            results.append({
                "base": base,
                "primer": primer,
                "expansion": exp_name,
                "plain_alphabet": "STANDARD",
                "cipher_alphabet": "STANDARD",
                "plaintext": result["plaintext"],
                "key_sample": result["key"],
                "quality": result["quality"]
            })
            
            # Only test keywords if the standard alphabet shows some promise
            if result["quality"]["quality"] > 0.3:
                print("    Promising result with standard alphabets, trying keywords...")
                
                # Try with keyword mixed alphabets
                for keyword in KEYWORDS:
                    # Create keyword alphabet
                    mixed_alphabet = create_keyword_alphabet(keyword)
                    
                    # Try as cipher alphabet
                    result = test_primer_configuration(
                        base, primer, exp_func, 
                        string.ascii_uppercase, mixed_alphabet
                    )
                    
                    results.append({
                        "base": base,
                        "primer": primer,
                        "expansion": exp_name,
                        "plain_alphabet": "STANDARD",
                        "cipher_alphabet": f"KEYWORD:{keyword}",
                        "plaintext": result["plaintext"],
                        "key_sample": result["key"],
                        "quality": result["quality"]
                    })
                    
                    # Try as plain alphabet
                    result = test_primer_configuration(
                        base, primer, exp_func, 
                        mixed_alphabet, string.ascii_uppercase
                    )
                    
                    results.append({
                        "base": base,
                        "primer": primer,
                        "expansion": exp_name,
                        "plain_alphabet": f"KEYWORD:{keyword}",
                        "cipher_alphabet": "STANDARD",
                        "plaintext": result["plaintext"],
                        "key_sample": result["key"],
                        "quality": result["quality"]
                    })
                    
                    # Try as both alphabets
                    result = test_primer_configuration(
                        base, primer, exp_func, 
                        mixed_alphabet, mixed_alphabet
                    )
                    
                    results.append({
                        "base": base,
                        "primer": primer,
                        "expansion": exp_name,
                        "plain_alphabet": f"KEYWORD:{keyword}",
                        "cipher_alphabet": f"KEYWORD:{keyword}",
                        "plaintext": result["plaintext"],
                        "key_sample": result["key"],
                        "quality": result["quality"]
                    })
    
    # Track end time
    end_time = time.time()
    
    # Sort results by quality
    results.sort(key=lambda x: x["quality"]["quality"], reverse=True)
    
    # Print top results
    print("\n\n===== TOP RESULTS =====")
    for i, result in enumerate(results[:10], 1):
        print(f"\n#{i}: Quality Score: {result['quality']['quality']:.4f}")
        print(f"Primer: {result['primer']} (Base {result['base']})")
        print(f"Expansion: {result['expansion']}")
        print(f"Plain Alphabet: {result['plain_alphabet']}")
        print(f"Cipher Alphabet: {result['cipher_alphabet']}")
        print(f"Key (first 10 digits): {result['key_sample']}")
        print(f"IC: {result['quality']['ic']:.6f}")
        print(f"Plaintext: {result['plaintext']}")
        
        # Print the plaintext at width 21
        print("Plaintext at width 21:")
        for j in range(0, len(result['plaintext']), 21):
            print(result['plaintext'][j:j+21])
    
    print(f"\nTotal configurations tested: {len(results)}")
    print(f"Time taken: {end_time - start_time:.2f} seconds")

if __name__ == "__main__":
    main() 