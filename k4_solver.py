#!/usr/bin/env python3
"""
Kryptos K4 Solver
Implementing the Gromark cipher approach based on Richard Bean's cryptodiagnosis paper
"""

import string
import itertools
from collections import Counter

# K4 ciphertext (97 letters)
K4_CIPHERTEXT = "OBKRUOXOGHULBSOLIFBBWFLRVQQPRNGKSSOTWTQSJQSSEKZZWATJKLUDIAWINFBNYPVTTMZFPKWGDKZXTJCDIGKUHUAUEKCAR"

# Known plaintext segments
KNOWN_PLAIN = ["EASTNORTHEAST", "BERLINCLOCK"]
KNOWN_CIPHER = ["FLRVQQPRNGKSS", "NYPVTTMZFPK"]

# Starting positions of known segments (0-indexed)
KNOWN_POSITIONS = [10, 62]

# Most promising primers from the cryptodiagnosis paper
PROMISING_PRIMERS = {
    (10, 5): [
        [2, 6, 7, 1, 7],  # equivalent to 84393
        [8, 4, 3, 9, 3],  # equivalent to 26717
        [9, 8, 8, 0, 0],  # Best IC value of 0.0625
    ],
    (10, 4): [
        [3, 3, 0, 1], 
        [6, 7, 4, 0], 
        [9, 9, 0, 3]
    ],
    (8, 5): [
        [0, 0, 3, 5, 1], 
        [0, 0, 5, 3, 7],
        [0, 0, 1, 3, 7], 
        [0, 0, 7, 7, 3]
    ],
    (5, 5): [
        [0, 0, 1, 2, 3], 
        [0, 1, 2, 3, 4], 
        [1, 2, 3, 4, 0]
    ],
    (12, 5): [
        [0, 1, 2, 3, 4], 
        [1, 2, 3, 4, 5], 
        [11, 0, 1, 2, 3]
    ]
}

def expand_gromark_primer(primer, base, length, expansion_func='standard'):
    """
    Expand a Gromark primer to generate a full key sequence.
    
    Args:
        primer (list): The initial primer values
        base (int): The number base (e.g., 10 for decimal)
        length (int): The desired length of the output key sequence
        expansion_func (str): The expansion function to use
        
    Returns:
        list: Expanded key sequence
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
    
    elif expansion_func == 'alternating':
        # Alternating addition and subtraction
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
    
    elif expansion_func == 'weighted':
        # Weighted expansion - double the last digit plus previous digit
        key = list(primer)
        while len(key) < length:
            next_digit = ((2 * key[-1]) + key[-2]) % base
            key.append(next_digit)
        return key
    
    else:
        raise ValueError(f"Unknown expansion function: {expansion_func}")

def decrypt_gromark(ciphertext, key, plain_alphabet, cipher_alphabet):
    """
    Decrypt a ciphertext using the Gromark cipher.
    
    Args:
        ciphertext (str): Ciphertext to decrypt
        key (list): Key sequence of numbers
        plain_alphabet (str): Plaintext alphabet
        cipher_alphabet (str): Ciphertext alphabet
        
    Returns:
        str: Decrypted plaintext
    """
    # Create mapping from cipher to plain
    cipher_to_plain = {c: i for i, c in enumerate(cipher_alphabet)}
    
    plaintext = []
    for i, c in enumerate(ciphertext):
        if c in cipher_to_plain:
            # Get the position in the cipher alphabet
            pos = cipher_to_plain[c]
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
    
    Args:
        text (str): Text to analyze
        
    Returns:
        float: Index of Coincidence
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

def score_english_text(text):
    """
    Score a text for English-like qualities.
    
    Args:
        text (str): Text to score
        
    Returns:
        float: Score (higher is more English-like)
    """
    # Common English letter frequencies
    english_freqs = {
        'E': 12.7, 'T': 9.1, 'A': 8.2, 'O': 7.5, 'I': 7.0, 'N': 6.7, 'S': 6.3, 
        'H': 6.1, 'R': 6.0, 'D': 4.3, 'L': 4.0, 'U': 2.8, 'C': 2.8, 'M': 2.4,
        'F': 2.2, 'W': 2.0, 'Y': 2.0, 'G': 1.9, 'P': 1.9, 'B': 1.5, 'V': 1.0,
        'K': 0.8, 'J': 0.2, 'X': 0.2, 'Q': 0.1, 'Z': 0.1
    }
    
    # Calculate letter frequencies in the text
    text = text.upper()
    text_counts = Counter(text)
    text_length = len(text)
    
    # Calculate chi-squared statistic
    chi2 = 0
    for letter, expected_freq in english_freqs.items():
        observed = text_counts.get(letter, 0)
        expected = expected_freq * text_length / 100
        chi2 += ((observed - expected) ** 2) / expected if expected > 0 else 0
    
    # Lower chi-squared means better fit to English
    return -chi2

def create_keyword_alphabet(keyword, reverse=False):
    """
    Create a keyword-mixed alphabet.
    
    Args:
        keyword (str): Keyword to use for alphabet mixing
        reverse (bool): Whether to reverse the resulting alphabet
        
    Returns:
        str: Mixed alphabet
    """
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

def verify_known_plaintext(decrypted, min_match_percentage=0.7):
    """
    Verify if the decrypted text contains the known plaintext segments.
    
    Args:
        decrypted (str): Decrypted text to check
        min_match_percentage (float): Minimum percentage of characters that must match
        
    Returns:
        bool: True if both segments match sufficiently, False otherwise
    """
    matches = []
    
    for i, (known_plain, position) in enumerate(zip(KNOWN_PLAIN, KNOWN_POSITIONS)):
        segment = decrypted[position:position+len(known_plain)]
        matches_count = sum(p == d for p, d in zip(known_plain, segment))
        match_percentage = matches_count / len(known_plain)
        matches.append(match_percentage >= min_match_percentage)
    
    return all(matches)

def try_primer(primer, base, expansion_func='standard', plain_alphabet=string.ascii_uppercase, 
               cipher_alphabet=string.ascii_uppercase, min_match_percentage=0.7):
    """
    Try a specific primer with the given parameters.
    
    Args:
        primer (list): Primer to test
        base (int): Number base to use
        expansion_func (str): Key expansion function to use
        plain_alphabet (str): Plaintext alphabet
        cipher_alphabet (str): Ciphertext alphabet
        min_match_percentage (float): Minimum match percentage for known plaintext
        
    Returns:
        tuple: (primer, plaintext, score, ic) if matching, None otherwise
    """
    # Expand the primer to get the key
    key = expand_gromark_primer(primer, base, len(K4_CIPHERTEXT), expansion_func)
    
    # Decrypt the ciphertext
    plaintext = decrypt_gromark(K4_CIPHERTEXT, key, plain_alphabet, cipher_alphabet)
    
    # Verify if it matches the known plaintext segments
    if verify_known_plaintext(plaintext, min_match_percentage):
        # Calculate index of coincidence
        ic = calculate_ic(plaintext)
        
        # Score for English-like qualities
        score = score_english_text(plaintext)
        
        return (primer, plaintext, score, ic)
    
    return None

def solve_k4():
    """
    Attempt to solve Kryptos K4 by trying various primers and parameters.
    """
    candidates = []
    
    # List of potential keywords for mixed alphabets
    keywords = [
        "KRYPTOS", "PALIMPSEST", "ABSCISSA", "SHADOW", "BERLIN", 
        "CLOCK", "NOISE", "DIGETAL", "LUCID", "UNDERDOG",
        # Add additional potential keywords
        "CIA", "LANGLEY", "SANBORN", "SCHEIDT", "ILLUSION",
        "SECRET", "HIDDEN", "ENCRYPTION", "CIPHER", "NORTHWEST"
    ]
    
    # List of expansion functions to try
    expansion_functions = [
        'standard', 'fibonacci', 'alternating', 'weighted'
    ]
    
    # First try the promising primers with standard alphabet
    for (base, length), primers in PROMISING_PRIMERS.items():
        print(f"Testing {len(primers)} primers with base {base}, length {length}")
        
        for primer in primers:
            for expansion_func in expansion_functions:
                # Try standard alphabets first
                result = try_primer(
                    primer, 
                    base, 
                    expansion_func=expansion_func,
                    plain_alphabet=string.ascii_uppercase,
                    cipher_alphabet=string.ascii_uppercase,
                    min_match_percentage=0.5  # More permissive matching
                )
                
                if result:
                    candidates.append(result + (expansion_func, 'STANDARD', 'STANDARD'))
                    print(f"CANDIDATE: {primer}, Base: {base}, Expansion: {expansion_func}")
                    print(f"Plaintext: {result[1]}")
    
    # If no candidates found, try with keyword alphabets
    if not candidates:
        print("No candidates with standard alphabets, trying keyword alphabets...")
        
        # Try most promising primers with keyword alphabets
        for (base, length), primers in PROMISING_PRIMERS.items():
            if (base, length) not in [(10, 5), (8, 5)]:  # Focus on most promising bases
                continue
                
            for primer in primers:
                for keyword in keywords:
                    # Try each keyword for plain, cipher, or both alphabets
                    plain_alphabets = [string.ascii_uppercase, create_keyword_alphabet(keyword)]
                    
                    for plain_alphabet in plain_alphabets:
                        # Try regular and reversed keyword alphabets for cipher
                        for rev in [False, True]:
                            cipher_alphabet = create_keyword_alphabet(keyword, reverse=rev)
                            
                            for expansion_func in expansion_functions:
                                result = try_primer(
                                    primer,
                                    base,
                                    expansion_func=expansion_func,
                                    plain_alphabet=plain_alphabet,
                                    cipher_alphabet=cipher_alphabet,
                                    min_match_percentage=0.5
                                )
                                
                                if result:
                                    plain_type = 'STANDARD' if plain_alphabet == string.ascii_uppercase else f'KEY:{keyword}'
                                    cipher_type = f'KEY:{keyword}' + (':REV' if rev else '')
                                    candidates.append(result + (expansion_func, plain_type, cipher_type))
                                    print(f"CANDIDATE: {primer}, Base: {base}, Expansion: {expansion_func}")
                                    print(f"Plain alphabet: {plain_type}")
                                    print(f"Cipher alphabet: {cipher_type}")
                                    print(f"Plaintext: {result[1]}")
    
    # Sort candidates by score
    candidates.sort(key=lambda x: x[2], reverse=True)
    
    # Display the top candidates
    print("\n===== TOP CANDIDATES =====")
    for i, (primer, plaintext, score, ic, exp_func, plain_type, cipher_type) in enumerate(candidates[:10], 1):
        print(f"\n#{i}: Primer: {primer}, Score: {score:.2f}, IC: {ic:.6f}")
        print(f"Expansion: {exp_func}, Plain: {plain_type}, Cipher: {cipher_type}")
        print(f"Plaintext: {plaintext}")
        
        # Print the plaintext arranged at width 21 (as discussed in the paper)
        print("Plaintext at width 21:")
        for j in range(0, len(plaintext), 21):
            print(plaintext[j:j+21])
    
    return candidates

if __name__ == "__main__":
    print("Attempting to solve Kryptos K4...")
    candidates = solve_k4()
    
    if not candidates:
        print("No viable candidates found. Try adjusting parameters or implementing additional cipher variations.") 