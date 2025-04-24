#!/usr/bin/env python3
"""
This script explores position-dependent approaches for the Kryptos K4 cipher,
focusing on how alphabet shifts might vary based on position in the text.
It tests various approaches including exact shifts observed in known plaintext,
row patterns, and other variations based on Richard Bean's findings.
"""

import string
import time
import argparse
from collections import Counter

# Constants
K4_CIPHERTEXT = "OBKRUOXOGHULBSOLIFBBWFLRVQQPRNGKSSOTWTQSJQSSEKZZWATJKLUDIAWINFBNYPVTTMZFPKWGDKZXTJCDIGKUHUAUEKCAR"
ROW_LENGTH = 21

# Known plaintext mappings based on K1-K3 techniques
KNOWN_PLAINTEXT = {
    10: ('E', 'F'),  # Position 10: Expected 'E', Observed 'F', Difference: 1
    11: ('A', 'L'),  # Position 11: Expected 'A', Observed 'L', Difference: 11
    12: ('S', 'R'),  # Position 12: Expected 'S', Observed 'R', Difference: 25
    13: ('T', 'V'),  # Position 13: Expected 'T', Observed 'V', Difference: 2
    14: ('N', 'Q'),  # Position 14: Expected 'N', Observed 'Q', Difference: 3
    15: ('O', 'Q'),  # Position 15: Expected 'O', Observed 'Q', Difference: 2
    16: ('R', 'P'),  # Position 16: Expected 'R', Observed 'P', Difference: 25
    17: ('T', 'R'),  # Position 17: Expected 'T', Observed 'R', Difference: 24
    18: ('H', 'N'),  # Position 18: Expected 'H', Observed 'N', Difference: 6
    19: ('E', 'G'),  # Position 19: Expected 'E', Observed 'G', Difference: 2
    20: ('A', 'K'),  # Position 20: Expected 'A', Observed 'K', Difference: 10
    21: ('S', 'S'),  # Position 21: Expected 'S', Observed 'S', Difference: 0
    22: ('T', 'S'),  # Position 22: Expected 'T', Observed 'S', Difference: 25
    
    62: ('B', 'N'),  # Position 62: Expected 'B', Observed 'N', Difference: 12
    63: ('E', 'Y'),  # Position 63: Expected 'E', Observed 'Y', Difference: 20
    64: ('R', 'P'),  # Position 64: Expected 'R', Observed 'P', Difference: 25
    65: ('L', 'V'),  # Position 65: Expected 'L', Observed 'V', Difference: 10
    66: ('I', 'T'),  # Position 66: Expected 'I', Observed 'T', Difference: 11
    67: ('N', 'T'),  # Position 67: Expected 'N', Observed 'T', Difference: 6
    68: ('C', 'M'),  # Position 68: Expected 'C', Observed 'M', Difference: 10
    69: ('L', 'Z'),  # Position 69: Expected 'L', Observed 'Z', Difference: 14
    70: ('O', 'F'),  # Position 70: Expected 'O', Observed 'F', Difference: 13
    71: ('C', 'P'),  # Position 71: Expected 'C', Observed 'P', Difference: 13
    72: ('K', 'K'),  # Position 72: Expected 'K', Observed 'K', Difference: 0
}

# Bean's best primer based on previous analyses
BEST_PRIMER = [2, 6, 7, 1, 7]

def expand_key(primer, length, method="standard"):
    """
    Expand a key primer to the required length.
    
    Args:
        primer: List of digits representing the primer.
        length: Desired length of the expanded key.
        method: Expansion method (standard or fibonacci).
        
    Returns:
        List of expanded key digits.
    """
    expanded_key = primer.copy()
    
    while len(expanded_key) < length:
        if method == "standard":
            # Standard Gromark expansion: next digit is sum of all digits mod base
            next_digit = sum(expanded_key) % 10
            expanded_key.append(next_digit)
        elif method == "fibonacci":
            # Fibonacci-style expansion: next digit is sum of last two digits mod base
            next_digit = (expanded_key[-1] + expanded_key[-2]) % 10
            expanded_key.append(next_digit)
    
    return expanded_key[:length]

def create_shifted_alphabet(shift):
    """
    Create an alphabet shifted by the specified amount.
    
    Args:
        shift: Amount to shift the alphabet (0-25).
        
    Returns:
        Shifted alphabet string.
    """
    alphabet = string.ascii_uppercase
    return alphabet[shift:] + alphabet[:shift]

def create_keyword_alphabet(keyword):
    """
    Create an alphabet based on a keyword.
    
    Args:
        keyword: Keyword for alphabet creation.
        
    Returns:
        Keyword-mixed alphabet.
    """
    # Ensure uppercase and remove duplicates while maintaining order
    keyword = keyword.upper()
    seen = set()
    unique_chars = [c for c in keyword if c not in seen and not seen.add(c)]
    
    # Add remaining alphabet characters
    alphabet = string.ascii_uppercase
    for c in alphabet:
        if c not in unique_chars:
            unique_chars.append(c)
    
    return "".join(unique_chars)

def create_position_based_alphabet(position, variation, plain_shift=0, cipher_shift=0):
    """
    Generate an alphabet based on the position in the text.
    
    Args:
        position: Position in the text.
        variation: Type of variation to use (exact, row_pattern, progressive, etc.).
        plain_shift: Base shift for plaintext alphabet.
        cipher_shift: Base shift for ciphertext alphabet.
        
    Returns:
        Tuple of (plaintext alphabet, ciphertext alphabet).
    """
    plain_alphabet = string.ascii_uppercase
    cipher_alphabet = string.ascii_uppercase
    
    if variation == "exact":
        # Use exact shifts as observed in known plaintext positions
        plain_alphabet = create_shifted_alphabet(plain_shift)
        cipher_alphabet = create_shifted_alphabet(cipher_shift)
        
    elif variation == "row_pattern":
        # Pattern based on row position (mod ROW_LENGTH)
        row_position = position % ROW_LENGTH
        p_shift = (plain_shift + row_position) % 26
        c_shift = (cipher_shift + row_position) % 26
        plain_alphabet = create_shifted_alphabet(p_shift)
        cipher_alphabet = create_shifted_alphabet(c_shift)
        
    elif variation == "progressive":
        # Progressive shift based on position
        p_shift = (plain_shift + position//10) % 26
        c_shift = (cipher_shift + position//10) % 26
        plain_alphabet = create_shifted_alphabet(p_shift)
        cipher_alphabet = create_shifted_alphabet(c_shift)
        
    elif variation == "alternating":
        # Alternating pattern based on odd/even position
        if position % 2 == 0:
            p_shift = plain_shift
            c_shift = cipher_shift
        else:
            p_shift = (plain_shift + 13) % 26  # Halfway shift
            c_shift = (cipher_shift + 13) % 26
        plain_alphabet = create_shifted_alphabet(p_shift)
        cipher_alphabet = create_shifted_alphabet(c_shift)
        
    elif variation == "bean_pattern":
        # Use a pattern described in Bean's paper
        mod_pos = position % 21
        if mod_pos < 7:
            # First third of row
            offset = 0
        elif mod_pos < 14:
            # Middle third of row
            offset = 5
        else:
            # Last third of row
            offset = 10
        
        p_shift = (plain_shift + offset) % 26
        c_shift = (cipher_shift + offset) % 26
        plain_alphabet = create_shifted_alphabet(p_shift)
        cipher_alphabet = create_shifted_alphabet(c_shift)
        
    elif variation == "keyword":
        # Use keyword-based alphabet
        keywords = [
            "KRYPTOS", 
            "BERLIN", 
            "SANBORN",
            "PALIMPSEST", 
            "ABSCISSA"
        ]
        # Choose keyword based on row number
        row_num = position // ROW_LENGTH
        keyword_idx = row_num % len(keywords)
        plain_alphabet = create_keyword_alphabet(keywords[keyword_idx])
        cipher_alphabet = string.ascii_uppercase
        
    return plain_alphabet, cipher_alphabet

def decrypt_gromark_position_based(ciphertext, key, variation, plain_shift=0, cipher_shift=0, expansion="standard"):
    """
    Decrypt Gromark cipher with position-based alphabets.
    
    Args:
        ciphertext: The ciphertext to decrypt.
        key: The key primer.
        variation: Type of position-based variation to use.
        plain_shift: Base shift for plaintext alphabet.
        cipher_shift: Base shift for ciphertext alphabet.
        expansion: Key expansion method.
        
    Returns:
        Decrypted plaintext.
    """
    # Expand key to match ciphertext length
    expanded_key = expand_key(key, len(ciphertext), expansion)
    
    # Decrypt using position-based alphabets
    plaintext = ""
    for i, char in enumerate(ciphertext):
        plain_alphabet, cipher_alphabet = create_position_based_alphabet(
            i, variation, plain_shift, cipher_shift
        )
        
        # Find position in cipher alphabet
        if char in cipher_alphabet:
            char_pos = cipher_alphabet.index(char)
            
            # Apply key shift and wrap around
            plain_pos = (char_pos - expanded_key[i]) % 26
            
            # Get plaintext character
            plaintext += plain_alphabet[plain_pos]
        else:
            plaintext += char
    
    return plaintext

def evaluate_plaintext(plaintext, expected_mappings):
    """
    Evaluate plaintext against known mappings.
    
    Args:
        plaintext: The plaintext to evaluate.
        expected_mappings: Dictionary of expected character mappings.
        
    Returns:
        Tuple of (match percentage, match count, mismatches).
    """
    matches = []
    mismatches = []
    
    for pos, (expected_plain, expected_cipher) in expected_mappings.items():
        if pos < len(plaintext):
            if plaintext[pos] == expected_plain:
                matches.append(pos)
            else:
                mismatches.append((pos, expected_plain, plaintext[pos], expected_cipher))
    
    match_count = len(matches)
    if match_count + len(mismatches) > 0:
        match_percentage = match_count / (match_count + len(mismatches)) * 100
    else:
        match_percentage = 0
    
    return match_percentage, matches, mismatches

def calculate_ic(text):
    """
    Calculate the Index of Coincidence for a given text.
    
    Args:
        text: Text to calculate IC for.
        
    Returns:
        Index of Coincidence value.
    """
    # Count character frequencies
    counts = Counter(text)
    n = len(text)
    
    # Calculate IC
    if n <= 1:
        return 0
    
    ic_sum = sum(counts[c] * (counts[c] - 1) for c in counts)
    ic = ic_sum / (n * (n - 1))
    
    return ic

def analyze_result_patterns(plaintext):
    """
    Analyze patterns in the plaintext for common English n-grams and words.
    
    Args:
        plaintext: The plaintext to analyze.
        
    Returns:
        Dictionary of analysis results.
    """
    results = {}
    
    # Check for common English bigrams
    common_bigrams = ["TH", "HE", "IN", "ER", "AN", "RE", "ON", "AT", "EN", "ND", "TI", "ES", "OR", "TE"]
    found_bigrams = []
    
    for bigram in common_bigrams:
        if bigram in plaintext:
            found_bigrams.append(bigram)
    
    # Check for common English trigrams
    common_trigrams = ["THE", "AND", "ING", "ENT", "ION", "TIO", "FOR", "OUR", "THI", "STH"]
    found_trigrams = []
    
    for trigram in common_trigrams:
        if trigram in plaintext:
            found_trigrams.append(trigram)
    
    # Check for common English words
    common_words = ["THE", "AND", "THAT", "HAVE", "FOR", "NOT", "WITH", "THIS", "FROM", "THEY"]
    found_words = []
    
    for word in common_words:
        if word in plaintext:
            found_words.append(word)
    
    results["common_bigrams"] = found_bigrams
    results["common_trigrams"] = found_trigrams
    results["common_words"] = found_words
    
    return results

def display_result(variation, expansion, plain_shift, cipher_shift, plaintext, match_percentage, matches, mismatches, ic):
    """
    Display detailed results of a decryption attempt.
    
    Args:
        variation: Type of variation used.
        expansion: Key expansion method used.
        plain_shift: Plaintext alphabet shift.
        cipher_shift: Ciphertext alphabet shift.
        plaintext: Decrypted plaintext.
        match_percentage: Percentage of characters matching expected.
        matches: List of positions with matching characters.
        mismatches: List of positions with mismatching characters.
        ic: Index of Coincidence value.
    """
    print(f"\nMatch: {match_percentage:.2f}%, IC: {ic:.6f}")
    print(f"Variation: {variation}")
    print(f"Plain shift: {plain_shift}, Cipher shift: {cipher_shift}")
    print(f"Expansion: {expansion}")
    print(f"Plaintext: {plaintext}")
    
    # Print plaintext in rows of ROW_LENGTH
    print("\nPlaintext at width 21:")
    for i in range(0, len(plaintext), ROW_LENGTH):
        print(plaintext[i:i+ROW_LENGTH])
    
    # Print pattern analysis
    patterns = analyze_result_patterns(plaintext)
    
    if patterns["common_bigrams"] or patterns["common_trigrams"] or patterns["common_words"]:
        print("\nInteresting patterns found:")
        
        if patterns["common_bigrams"]:
            print(f"  Common bigrams: {', '.join(patterns['common_bigrams'])}")
        
        if patterns["common_trigrams"]:
            print(f"  Common trigrams: {', '.join(patterns['common_trigrams'])}")
        
        if patterns["common_words"]:
            print(f"  Common words: {', '.join(patterns['common_words'])}")
    
    # Print matched positions
    if matches:
        print("\nMatched positions:")
        for pos in matches:
            expected_plain, expected_cipher = KNOWN_PLAINTEXT[pos]
            print(f"  Position {pos}: Expected '{expected_plain}', Got '{plaintext[pos]}', Ciphertext '{expected_cipher}'")
    
    # Print mismatched positions
    if mismatches:
        print("\nMismatched positions:")
        for pos, expected, actual, cipher in mismatches:
            print(f"  Position {pos}: Expected '{expected}', Got '{actual}', Ciphertext '{cipher}'")

def test_custom_keyword_approach():
    """
    Test a custom approach using keyword-based alphabets with position dependency.
    
    This approach combines keyword-based alphabets with position-based shifts,
    inspired by the findings from previous analyses.
    """
    results = []
    
    # Test with promising keywords from Kryptos context
    keywords = [
        "KRYPTOS",
        "BERLIN",
        "SANBORN",
        "PALIMPSEST",
        "ABSCISSA",
        "BERLINCLOCK",
        "NYPVTTMZFPK",  # Segment from K4
        "EASTBERLIN",
        "WESTBERLIN",
        "BELOUSOV"  # Referenced in Kryptos
    ]
    
    primer = BEST_PRIMER
    
    print("\n=== TESTING CUSTOM KEYWORD APPROACH ===")
    
    for keyword in keywords:
        print(f"Testing keyword: {keyword}")
        
        # Create plaintext alphabet from keyword
        plain_alphabet = create_keyword_alphabet(keyword)
        
        # Try various position-dependent strategies for cipher alphabet
        for position_pattern in ["fixed", "row_based", "progressive"]:
            plaintext = ""
            expanded_key = expand_key(primer, len(K4_CIPHERTEXT), "standard")
            
            for i, char in enumerate(K4_CIPHERTEXT):
                # Use keyword-based alphabet for plaintext
                if position_pattern == "fixed":
                    cipher_alphabet = string.ascii_uppercase
                elif position_pattern == "row_based":
                    # Shift based on row position
                    row_pos = i % ROW_LENGTH
                    cipher_alphabet = create_shifted_alphabet(row_pos)
                elif position_pattern == "progressive":
                    # Progressive shift based on key digit
                    cipher_alphabet = create_shifted_alphabet(expanded_key[i])
                
                # Find position in cipher alphabet
                if char in cipher_alphabet:
                    char_pos = cipher_alphabet.index(char)
                    
                    # Apply key shift and wrap
                    plain_pos = (char_pos - expanded_key[i]) % 26
                    
                    # Get plaintext character
                    if plain_pos < len(plain_alphabet):
                        plaintext += plain_alphabet[plain_pos]
                    else:
                        plaintext += "?"
                else:
                    plaintext += char
            
            # Evaluate the result
            match_percentage, matches, mismatches = evaluate_plaintext(plaintext, KNOWN_PLAINTEXT)
            ic = calculate_ic(plaintext)
            
            results.append((match_percentage, ic, keyword, position_pattern, plaintext, matches, mismatches))
    
    # Sort results by match percentage and IC
    results.sort(key=lambda x: (x[0], x[1]), reverse=True)
    
    # Display top results
    print("\n=== TOP KEYWORD RESULTS ===")
    for i, (match_percentage, ic, keyword, pattern, plaintext, matches, mismatches) in enumerate(results[:10]):
        print(f"\n#{i+1}: Match: {match_percentage:.2f}%, IC: {ic:.6f}")
        print(f"Keyword: {keyword}")
        print(f"Pattern: {pattern}")
        print(f"Plaintext: {plaintext}")
        
        # Print plaintext in rows of ROW_LENGTH
        print("\nPlaintext at width 21:")
        for j in range(0, len(plaintext), ROW_LENGTH):
            print(plaintext[j:j+ROW_LENGTH])
        
        # Print pattern analysis
        patterns = analyze_result_patterns(plaintext)
        
        if patterns["common_bigrams"] or patterns["common_trigrams"] or patterns["common_words"]:
            print("\nInteresting patterns found:")
            
            if patterns["common_bigrams"]:
                print(f"  Common bigrams: {', '.join(patterns['common_bigrams'])}")
            
            if patterns["common_trigrams"]:
                print(f"  Common trigrams: {', '.join(patterns['common_trigrams'])}")
            
            if patterns["common_words"]:
                print(f"  Common words: {', '.join(patterns['common_words'])}")
        
        # Print matched positions
        if matches:
            print("\nMatched positions:")
            for pos in matches:
                expected_plain, expected_cipher = KNOWN_PLAINTEXT[pos]
                print(f"  Position {pos}: Expected '{expected_plain}', Got '{plaintext[pos]}', Ciphertext '{expected_cipher}'")

def main():
    """Main function to test position-dependent approaches."""
    start_time = time.time()
    
    parser = argparse.ArgumentParser(description="Test position-dependent approaches for Kryptos K4")
    parser.add_argument("--custom", action="store_true", help="Test custom keyword approach")
    args = parser.parse_args()
    
    if args.custom:
        test_custom_keyword_approach()
    else:
        # List of variations to test
        variations = ["exact", "row_pattern", "progressive", "alternating", "bean_pattern"]
        
        # Track results for comparison
        results = []
        
        print("=== ANALYZING POSITION-DEPENDENT APPROACHES ===")
        
        # Test each variation with different shifts
        for variation in variations:
            print(f"\nTesting variation: {variation}")
            
            # Test with different shifts
            for plain_shift in range(0, 26, 3):  # Try different shifts
                for cipher_shift in range(0, 26, 3):
                    for expansion in ["standard"]:  # Could also try "fibonacci"
                        # Decrypt using current configuration
                        plaintext = decrypt_gromark_position_based(
                            K4_CIPHERTEXT, BEST_PRIMER, variation, 
                            plain_shift, cipher_shift, expansion
                        )
                        
                        # Evaluate against known plaintext
                        match_percentage, matches, mismatches = evaluate_plaintext(plaintext, KNOWN_PLAINTEXT)
                        
                        # Calculate Index of Coincidence
                        ic = calculate_ic(plaintext)
                        
                        # Store result if it's promising
                        if match_percentage > 0:
                            results.append((
                                match_percentage, ic, variation, expansion, 
                                plain_shift, cipher_shift, plaintext, matches, mismatches
                            ))
        
        # Sort results by match percentage and IC
        results.sort(key=lambda x: (x[0], x[1]), reverse=True)
        
        # Display top results
        print("\n=== TOP RESULTS ===")
        for i, result in enumerate(results[:10]):
            match_percentage, ic, variation, expansion, plain_shift, cipher_shift, plaintext, matches, mismatches = result
            print(f"\n#{i+1}: Match: {match_percentage:.2f}%, IC: {ic:.6f}")
            print(f"Variation: {variation}")
            print(f"Plain shift: {plain_shift}, Cipher shift: {cipher_shift}")
            print(f"Expansion: {expansion}")
            print(f"Plaintext: {plaintext}")
            
            # Print plaintext in rows of ROW_LENGTH
            print("\nPlaintext at width 21:")
            for j in range(0, len(plaintext), ROW_LENGTH):
                print(plaintext[j:j+ROW_LENGTH])
            
            # Print pattern analysis
            patterns = analyze_result_patterns(plaintext)
            
            if patterns["common_bigrams"] or patterns["common_trigrams"] or patterns["common_words"]:
                print("\nInteresting patterns found:")
                
                if patterns["common_bigrams"]:
                    print(f"  Common bigrams: {', '.join(patterns['common_bigrams'])}")
                
                if patterns["common_trigrams"]:
                    print(f"  Common trigrams: {', '.join(patterns['common_trigrams'])}")
                
                if patterns["common_words"]:
                    print(f"  Common words: {', '.join(patterns['common_words'])}")
            
            # Print matched positions
            if matches:
                print("\nMatched positions:")
                for pos in matches:
                    expected_plain, expected_cipher = KNOWN_PLAINTEXT[pos]
                    print(f"  Position {pos}: Expected '{expected_plain}', Got '{plaintext[pos]}', Ciphertext '{expected_cipher}'")
            
            # Print mismatched positions
            if mismatches:
                print("\nMismatched positions:")
                for pos, expected, actual, cipher in mismatches:
                    print(f"  Position {pos}: Expected '{expected}', Got '{actual}', Ciphertext '{cipher}'")
    
    end_time = time.time()
    print(f"\nTotal execution time: {end_time - start_time:.2f} seconds")

if __name__ == "__main__":
    main() 