#!/usr/bin/env python3
"""
Gromark Cipher Variations for Kryptos K4
This script provides utilities to test different variations of the Gromark cipher,
including keyword-mixed alphabets and different key expansion rules.
"""

import string
import argparse
from gromark_solver import brute_force_gromark, print_results, K4_CIPHERTEXT, expand_gromark_primer, cpu_search_primers

def create_keyword_alphabet(keyword):
    """
    Create a keyword-mixed alphabet.
    
    Args:
        keyword (str): Keyword to use for alphabet mixing
        
    Returns:
        str: Mixed alphabet
    """
    # Convert to uppercase and remove duplicates while preserving order
    keyword = keyword.upper()
    seen = set()
    unique_keyword = ''.join(c for c in keyword if c not in seen and not seen.add(c))
     
    # Remove keyword letters from the remaining alphabet
    remaining = ''.join(c for c in string.ascii_uppercase if c not in unique_keyword)
    
    # Return the mixed alphabet
    return unique_keyword + remaining

def fibonacci_expansion(primer, base, length):
    """
    Expand a primer using Fibonacci-style addition (a variant of Gromark).
    
    Args:
        primer (list): Initial primer digits
        base (int): The number base to use
        length (int): Desired length of the key sequence
        
    Returns:
        list: Expanded key sequence
    """
    key = list(primer)
    while len(key) < length:
        next_digit = 0
        # Sum the last n digits (where n is the primer length)
        for i in range(1, len(primer) + 1):
            if i <= len(key):
                next_digit = (next_digit + key[-i]) % base
        key.append(next_digit)
    return key[:length]

def test_keyword_alphabets(keywords, base=10, primer_length=5, processes=None):
    """
    Test a list of keywords as mixed alphabets for both plaintext and ciphertext.
    
    Args:
        keywords (list): List of keywords to try
        base (int): The number base to use
        primer_length (int): Length of primers to test
        processes (int): Number of processes to use
        
    Returns:
        dict: Dictionary mapping keywords to candidate lists
    """
    results = {}
    
    print(f"Testing {len(keywords)} keyword alphabets...")
    
    for keyword in keywords:
        mixed_alphabet = create_keyword_alphabet(keyword)
        print(f"\nTesting keyword: {keyword}")
        print(f"Mixed alphabet: {mixed_alphabet}")
        
        # Test with mixed alphabet as both plaintext and ciphertext alphabet
        candidates = brute_force_gromark(
            base=base,
            primer_length=primer_length,
            ciphertext=K4_CIPHERTEXT,
            plain_alphabet=string.ascii_uppercase,  # Standard plaintext alphabet
            cipher_alphabet=mixed_alphabet,        # Mixed ciphertext alphabet
            processes=processes
        )
        
        if candidates:
            results[keyword] = candidates
            print(f"Found {len(candidates)} candidates with keyword '{keyword}'")
            print_results(candidates, top_n=3)  # Show top 3 for each keyword
        else:
            print(f"No candidates found with keyword '{keyword}'")
    
    return results

def custom_key_expansion(base=10, primer_length=5, processes=None, expansion_rules=None):
    """
    Test Gromark variations with different key expansion rules.
    
    Args:
        base (int): The number base to use
        primer_length (int): Length of primers to test
        processes (int): Number of processes to use
        expansion_rules (list): List of tuples (name, function) for expansion rules to test
    """
    if expansion_rules is None:
        expansion_rules = [
            ("Standard", expand_gromark_primer),
            ("Fibonacci", fibonacci_expansion)
        ]
    
    # Use a subset of good primers from Bean's paper for testing
    test_primers = [
        [2, 6, 7, 1, 7], [8, 4, 3, 9, 3], [9, 8, 8, 0, 0]
    ]
    
    print(f"Testing {len(test_primers)} promising primers with {len(expansion_rules)} expansion rules...")
    
    for name, expansion_fn in expansion_rules:
        print(f"\n=== Testing {name} expansion rule ===")
        
        candidates = cpu_search_primers(
            ciphertext=K4_CIPHERTEXT,
            primer_list=test_primers,
            base=base,
            length=primer_length,
            expansion_fn=expansion_fn,
            max_primers=None,
            min_ic=0.06
        )
        
        if candidates:
            print(f"Found {len(candidates)} candidates with {name} expansion rule")
            print_results(candidates, top_n=3)
        else:
            print(f"No candidates found with {name} expansion rule")
            
    return expansion_rules

def try_different_expansion_rules(base=10, primer_length=5, processes=None):
    """
    Test different key expansion rules.
    
    Args:
        base (int): The number base to use
        primer_length (int): Length of primers to test
        processes (int): Number of processes to use
    """
    # Define custom expansion rules
    def running_sum_expansion(primer, base, length):
        """Expansion that uses a running sum of the last primer_length digits"""
        return fibonacci_expansion(primer, base, length)
    
    def alternating_expansion(primer, base, length):
        """Expansion that alternates between addition and subtraction"""
        key = list(primer)
        while len(key) < length:
            # If at even position, add; if at odd, subtract
            if len(key) % 2 == 0:
                next_digit = (key[-1] + key[-2]) % base
            else:
                next_digit = (key[-1] - key[-2]) % base
            key.append(next_digit)
        return key[:length]
    
    def weighted_expansion(primer, base, length):
        """Expansion that weights the last two digits differently"""
        key = list(primer)
        while len(key) < length:
            next_digit = (2 * key[-1] + key[-2]) % base
            key.append(next_digit)
        return key[:length]
    
    # List of expansion rules to test
    expansion_rules = [
        ("Standard", expand_gromark_primer),
        ("Fibonacci", fibonacci_expansion),
        ("Running Sum", running_sum_expansion),
        ("Alternating", alternating_expansion),
        ("Weighted", weighted_expansion)
    ]
    
    return custom_key_expansion(base, primer_length, processes, expansion_rules)

def main():
    parser = argparse.ArgumentParser(description='Test variations of the Gromark cipher.')
    parser.add_argument('--mode', choices=['keywords', 'expansions', 'all'], default='all',
                       help='Test mode: keywords, expansions, or all')
    parser.add_argument('--base', type=int, default=10, help='Number base (default: 10)')
    parser.add_argument('--length', type=int, default=5, help='Primer length (default: 5)')
    parser.add_argument('--processes', type=int, default=None, help='Number of processes to use')
    
    args = parser.parse_args()
    
    # Default keywords to test
    test_keywords = ["KRYPTOS", "PALIMPSEST", "ABSCISSA", "SHADOW", "BERLIN", 
                    "CLOCK", "NOISE", "DIGETAL", "LUCID", "UNDERDOG"]
    
    if args.mode in ['keywords', 'all']:
        print("\n=== Testing Different Keyword Alphabets ===")
        test_keyword_alphabets(test_keywords, args.base, args.length, args.processes)
    
    if args.mode in ['expansions', 'all']:
        print("\n=== Testing Different Expansion Rules ===")
        try_different_expansion_rules(args.base, args.length, args.processes)
    
    print("\nAll tests completed.")

if __name__ == "__main__":
    main() 