#!/usr/bin/env python3
"""
Kryptos K4 Primer Constraint Solver
Based on Richard Bean's cryptodiagnosis paper

This script implements the constraint-based approach to identify valid
Gromark primers that could produce the K4 ciphertext with the known plaintext.
"""

import string
import itertools
from collections import defaultdict

# K4 ciphertext (97 letters)
K4_CIPHERTEXT = "OBKRUOXOGHULBSOLIFBBWFLRVQQPRNGKSSOTWTQSJQSSEKZZWATJKLUDIAWINFBNYPVTTMZFPKWGDKZXTJCDIGKUHUAUEKCAR"

# Known plaintext mapping
# "EASTNORTHEAST" corresponds to "FLRVQQPRNGKSS" at position 10 (0-indexed)
# "BERLINCLOCK" corresponds to "NYPVTTMZFPK" at position 62 (0-indexed)
PLAINTEXT_MAP = {
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

def expand_primer(primer, base, length):
    """Expand a primer to generate a full key sequence."""
    key = list(primer)
    while len(key) < length:
        next_digit = (key[-1] + key[-2]) % base
        key.append(next_digit)
    return key

def get_key_constraints(base=10, key_length=5):
    """
    Derive constraints on the key digits based on known plaintext/ciphertext pairs.
    
    Returns a dictionary of constraints:
    - Equalities: key positions that must have the same value
    - Inequalities: key positions that must have different values
    """
    # Dictionary to track constraint types
    equalities = set()
    inequalities = set()
    
    # Derive constraints from plaintext/ciphertext mapping
    positions = sorted(PLAINTEXT_MAP.keys())
    
    # Compare each pair of positions
    for i, pos1 in enumerate(positions):
        plain1, cipher1 = PLAINTEXT_MAP[pos1]
        
        for pos2 in positions[i+1:]:
            plain2, cipher2 = PLAINTEXT_MAP[pos2]
            
            # If same plaintext maps to same ciphertext, key values must be equal
            if plain1 == plain2 and cipher1 == cipher2:
                equalities.add((pos1, pos2))
            
            # If same plaintext maps to different ciphertext, key values must be different
            elif plain1 == plain2 and cipher1 != cipher2:
                inequalities.add((pos1, pos2))
            
            # If different plaintext maps to same ciphertext, related constraint on key values
            elif plain1 != plain2 and cipher1 == cipher2:
                # The difference in key values must equal the difference in plaintext positions
                # Not directly expressible as simple equality/inequality
                pass
    
    # For standard alphabets, calculate key value from plaintext and ciphertext positions
    required_key_values = {}
    for pos, (plain, cipher) in PLAINTEXT_MAP.items():
        p_val = ord(plain) - ord('A')
        c_val = ord(cipher) - ord('A')
        # Key digit at this position = (cipher_pos - plain_pos) mod 26
        key_val = (c_val - p_val) % 26
        # For base 10, we need to ensure the key value is < 10
        # This is an additional constraint
        if base <= key_val:
            pass  # Can't use standard A-Z mapping with this base
        else:
            required_key_values[pos] = key_val
    
    # Add derived equalities from known key values
    for pos1, val1 in required_key_values.items():
        for pos2, val2 in required_key_values.items():
            if pos1 < pos2:
                if val1 == val2:
                    equalities.add((pos1, pos2))
                else:
                    inequalities.add((pos1, pos2))
    
    # Make the constraints canonical (smaller position first)
    equalities = {(min(a, b), max(a, b)) for a, b in equalities}
    inequalities = {(min(a, b), max(a, b)) for a, b in inequalities}
    
    return {
        "equalities": equalities,
        "inequalities": inequalities,
        "key_values": required_key_values
    }

def print_constraints(constraints):
    """Print the derived constraints in a readable format."""
    print("=== KEY CONSTRAINTS ===")
    
    print("\nEqualities (key positions that must have the same value):")
    for pos1, pos2 in sorted(constraints["equalities"]):
        print(f"k{pos1} = k{pos2}")
    
    print("\nInequalities (key positions that must have different values):")
    for pos1, pos2 in sorted(constraints["inequalities"]):
        print(f"k{pos1} ≠ k{pos2}")
    
    print("\nRequired key values:")
    for pos, val in sorted(constraints["key_values"].items()):
        print(f"k{pos} = {val}")

def check_key_validity(key, constraints):
    """Check if a key satisfies all the constraints."""
    # Check equalities
    for pos1, pos2 in constraints["equalities"]:
        if pos1 < len(key) and pos2 < len(key) and key[pos1] != key[pos2]:
            return False
    
    # Check inequalities
    for pos1, pos2 in constraints["inequalities"]:
        if pos1 < len(key) and pos2 < len(key) and key[pos1] == key[pos2]:
            return False
    
    # Check required values
    for pos, val in constraints["key_values"].items():
        if pos < len(key) and key[pos] != val:
            return False
    
    return True

def find_valid_primers(base, primer_length):
    """Find all valid primers that satisfy the constraints."""
    # Generate constraints
    constraints = get_key_constraints(base, primer_length)
    print_constraints(constraints)
    
    # Derive constraints specifically for the first 'primer_length' positions
    # This would require analyzing how the key expansion affects later positions
    
    # As a brute force approach, we can try all possible primers
    # and check if the expanded key satisfies the constraints
    valid_primers = []
    total_primers = base ** primer_length
    
    print(f"\nSearching through {total_primers} possible primers for base {base}, length {primer_length}...")
    
    # Generate all possible primers
    for primer_idx in range(total_primers):
        # Convert index to primer digits
        primer = []
        value = primer_idx
        for i in range(primer_length):
            primer.insert(0, value % base)
            value //= base
        
        # Expand primer to full key
        key = expand_primer(primer, base, max(PLAINTEXT_MAP.keys()) + 1)
        
        # Check if key satisfies constraints
        if check_key_validity(key, constraints):
            valid_primers.append(primer)
            print(f"Found valid primer: {primer}")
    
    return valid_primers

def main():
    """Main function to search for valid primers."""
    # Check a range of bases and primer lengths
    bases_to_check = [10, 8, 5, 12]
    lengths_to_check = [4, 5]
    
    all_valid_primers = {}
    
    for base in bases_to_check:
        for length in lengths_to_check:
            print(f"\n\n=== SEARCHING BASE {base}, LENGTH {length} ===\n")
            valid_primers = find_valid_primers(base, length)
            all_valid_primers[(base, length)] = valid_primers
            print(f"Found {len(valid_primers)} valid primers for base {base}, length {length}")
    
    # Print summary
    print("\n=== SUMMARY OF VALID PRIMERS ===")
    for (base, length), primers in all_valid_primers.items():
        print(f"Base {base}, Length {length}: {len(primers)} valid primers found")
        if primers:
            print(f"Examples: {primers[:5]}")

if __name__ == "__main__":
    main() 