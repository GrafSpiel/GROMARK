#!/usr/bin/env python3
"""
K4 Gromark Comprehensive Solver
This script implements a systematic approach to solve Kryptos K4 using the Gromark cipher,
incorporating known plaintext, statistical properties, and multiple key expansion methods.
Based on Richard Bean's cryptodiagnosis paper and known constraints.
"""

import string
import itertools
import time
import random
import math
from collections import Counter, defaultdict

# K4 ciphertext (97 letters)
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

# Most promising primers from Bean's paper
PROMISING_PRIMERS = [
    # Format: (primer, base, description)
    ([2, 6, 7, 1, 7], 10, "Bean's paper - identified as one of two equivalent primers"),
    ([8, 4, 3, 9, 3], 10, "Bean's paper - identified as one of two equivalent primers"),
    ([9, 8, 8, 0, 0], 10, "Bean's paper - gives highest IC of 0.0625"),
    ([0, 0, 3, 5, 1], 8, "Base 8 primer with period 84"),
    ([0, 0, 5, 3, 7], 8, "Base 8 primer with period 84"),
    ([3, 3, 0, 1], 10, "Base 10, length 4 primer"),
    ([6, 7, 4, 0], 10, "Base 10, length 4 primer"),
    ([9, 9, 0, 3], 10, "Base 10, length 4 primer")
]

# Common English words for scoring
COMMON_WORDS = set([
    "THE", "AND", "THAT", "HAVE", "FOR", "NOT", "WITH", "YOU", "THIS", "BUT",
    "HIS", "FROM", "THEY", "SAY", "SHE", "WILL", "ONE", "ALL", "WOULD", "THERE",
    "THEIR", "WHAT", "OUT", "ABOUT", "WHO", "GET", "WHICH", "WHEN", "MAKE", "CAN",
    "LIKE", "TIME", "JUST", "HIM", "KNOW", "TAKE", "PEOPLE", "INTO", "YEAR", "YOUR",
    "GOOD", "SOME", "COULD", "THEM", "SEE", "OTHER", "THAN", "THEN", "NOW", "LOOK",
    "ONLY", "COME", "ITS", "OVER", "THINK", "ALSO", "BACK", "AFTER", "USE", "TWO",
    "HOW", "OUR", "WORK", "FIRST", "WELL", "WAY", "EVEN", "NEW", "WANT", "BECAUSE",
    "ANY", "THESE", "GIVE", "DAY", "MOST", "BERLIN", "CLOCK", "EAST", "NORTH", "UNDERGROUND",
    "LANGLEY", "CIPHER", "CODE", "SECRET", "HIDDEN", "SPY", "AGENT", "CIA", "KRYPTOS"
])

# Kryptos-related keywords for alphabet mixing
KEYWORDS = [
    "KRYPTOS", "PALIMPSEST", "ABSCISSA", "SHADOW", "BERLIN", 
    "CLOCK", "LANGLEY", "SANBORN", "SCHEIDT", "CIA", 
    "BERLINCLOCK", "EASTNORTHEAST", "UNDERGROUND", "HIDDEN"
]

def expand_key(primer, base, length, method="standard"):
    """
    Expand a primer to generate a full key sequence based on different methods.
    
    Args:
        primer (list): Initial primer digits
        base (int): Number base to use
        length (int): Desired length of expanded key
        method (str): Key expansion method to use
        
    Returns:
        list: Expanded key sequence
    """
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
            
        elif method == "alternating":
            # Alternating addition and subtraction
            if len(key) % 2 == 0:
                next_digit = (key[-1] + key[-2]) % base
            else:
                next_digit = (key[-1] - key[-2]) % base
                if next_digit < 0:
                    next_digit += base
            key.append(next_digit)
            
        elif method == "weighted":
            # Double the last digit plus previous digit
            next_digit = ((2 * key[-1]) + key[-2]) % base
            key.append(next_digit)
            
        elif method == "m-sequence":
            # Implementation of Golomb's m-sequence (period 63)
            if base == 4:  # Specific to F₂²
                # Map 0,1,2,3 to field elements in F₂²
                field_map = {0: [0, 0], 1: [0, 1], 2: [1, 0], 3: [1, 1]}
                inv_map = {tuple(v): k for k, v in field_map.items()}
                
                # Get the last two elements in F₂²
                a = field_map[key[-2]]
                b = field_map[key[-1]]
                
                # Compute the next element using the polynomial x² + x + 1
                # This is a simplified version of field multiplication
                c0 = (a[1] ^ b[0]) % 2
                c1 = (a[0] ^ a[1] ^ b[1]) % 2
                
                next_digit = inv_map[(c0, c1)]
                key.append(next_digit)
            else:
                # Fall back to standard method for other bases
                next_digit = (key[-1] + key[-2]) % base
                key.append(next_digit)
                
        else:
            # Default to standard method
            next_digit = (key[-1] + key[-2]) % base
            key.append(next_digit)
            
    return key[:length]

def check_key_periodicity(key, max_period=100):
    """
    Check if a key has a periodic pattern and return the period if found.
    
    Args:
        key (list): Key sequence to check
        max_period (int): Maximum period to check
        
    Returns:
        int or None: Period length if found, None otherwise
    """
    # Start with period lengths that might be interesting
    for period in [21, 42, 63, 84]:
        if period >= len(key):
            continue
            
        is_periodic = True
        for i in range(len(key) - period):
            if key[i] != key[i + period]:
                is_periodic = False
                break
                
        if is_periodic:
            return period
    
    # Check other periods
    for period in range(1, min(max_period, len(key) // 2)):
        is_periodic = True
        for i in range(len(key) - period):
            if key[i] != key[i + period]:
                is_periodic = False
                break
                
        if is_periodic:
            return period
    
    return None

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

def create_caesar_shifted_alphabet(shift):
    """
    Create an alphabet shifted by a given amount (Caesar shift).
    
    Args:
        shift (int): Amount to shift the alphabet
        
    Returns:
        str: Shifted alphabet
    """
    return ''.join(chr(((ord(c) - ord('A') + shift) % 26) + ord('A')) for c in string.ascii_uppercase)

def decrypt_gromark(ciphertext, key, plain_alphabet, cipher_alphabet):
    """
    Decrypt using the Gromark cipher.
    
    Args:
        ciphertext (str): Ciphertext to decrypt
        key (list): Key sequence
        plain_alphabet (str): Plaintext alphabet
        cipher_alphabet (str): Ciphertext alphabet
        
    Returns:
        str: Decrypted plaintext
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

def undo_masking(text, shift=0, method="caesar"):
    """
    Undo the masking step applied to the plaintext.
    
    Args:
        text (str): Text to unmask
        shift (int): Shift value for Caesar shift
        method (str): Masking method (caesar or custom)
        
    Returns:
        str: Unmasked text
    """
    if method == "caesar":
        # Simple Caesar shift
        return ''.join(chr(((ord(c) - ord('A') - shift) % 26) + ord('A')) for c in text)
    else:
        # For future custom masking methods
        return text

def calculate_ic(text):
    """
    Calculate the Index of Coincidence for a text.
    
    Args:
        text (str): Text to analyze
        
    Returns:
        float: Index of Coincidence
    """
    # Count each letter
    counts = Counter(c for c in text.upper() if c in string.ascii_uppercase)
    
    # Calculate IC
    n = sum(counts.values())
    if n <= 1:
        return 0
    
    sum_freqs = sum(count * (count - 1) for count in counts.values())
    
    return sum_freqs / (n * (n - 1))

def count_width21_bigrams(text):
    """
    Count repeated bigrams when text is arranged at width 21.
    
    Args:
        text (str): Text to analyze
        
    Returns:
        int: Number of repeated vertical bigrams
    """
    bigrams = []
    for col in range(21):
        for row in range(len(text) // 21):
            if row + 1 < len(text) // 21:
                pos1 = row * 21 + col
                pos2 = (row + 1) * 21 + col
                if pos1 < len(text) and pos2 < len(text):
                    bigrams.append(text[pos1] + text[pos2])
    
    # Count repeated bigrams
    bigram_counts = Counter(bigrams)
    repeated = sum(count - 1 for count in bigram_counts.values() if count > 1)
    
    return repeated

def verify_known_plaintext(plaintext, min_match_percentage=0.5):
    """
    Verify if the plaintext matches the known plaintext segments.
    
    Args:
        plaintext (str): Plaintext to check
        min_match_percentage (float): Minimum percentage of matches required
        
    Returns:
        tuple: (match_percentage, matches, total)
    """
    matches = 0
    for pos, (expected, _) in KNOWN_PLAINTEXT.items():
        if pos < len(plaintext) and plaintext[pos] == expected:
            matches += 1
    
    match_percentage = matches / len(KNOWN_PLAINTEXT)
    
    return match_percentage, matches, len(KNOWN_PLAINTEXT)

def score_plaintext(text, known_match_percentage):
    """
    Score plaintext for English-like qualities.
    
    Args:
        text (str): Text to score
        known_match_percentage (float): Percentage of known plaintext matched
        
    Returns:
        float: Score (higher is better)
    """
    # If known plaintext doesn't match well, give very low score
    if known_match_percentage < 0.5:
        return -1000
    
    # Calculate Index of Coincidence (closer to 0.067 is better)
    ic = calculate_ic(text)
    ic_quality = 1 - abs(ic - 0.067) / 0.067
    
    # Count English words
    word_count = 0
    for length in range(3, 10):  # Check words of different lengths
        for i in range(len(text) - length + 1):
            if text[i:i+length] in COMMON_WORDS:
                word_count += length  # Give more weight to longer words
    
    # Calculate bigram and trigram frequencies compared to English
    english_bigrams = {
        'TH': 3.56, 'HE': 3.07, 'IN': 2.43, 'ER': 2.05, 'AN': 1.99,
        'RE': 1.85, 'ON': 1.76, 'AT': 1.49, 'EN': 1.45, 'ND': 1.35
    }
    
    bigram_score = 0
    for bigram, expected_freq in english_bigrams.items():
        count = text.count(bigram)
        bigram_score += min(count * expected_freq * 0.1, 5)  # Cap the impact
    
    # Count width 21 bigrams (should be close to 11 for K4)
    width21_count = count_width21_bigrams(text)
    width21_score = 5 - abs(width21_count - 11)  # Penalize deviation from 11
    
    # Combined score
    score = (
        known_match_percentage * 50 +  # Weight matches heavily
        ic_quality * 20 +               # Weight IC quality
        word_count * 0.5 +              # Weight word matches
        bigram_score * 0.2 +            # Weight bigram frequency
        width21_score * 2                # Weight width 21 property
    )
    
    return score

def analyze_minor_differences(plaintext, ciphertext):
    """
    Analyze the minor differences between plaintext and ciphertext.
    
    Args:
        plaintext (str): Plaintext
        ciphertext (str): Corresponding ciphertext
        
    Returns:
        tuple: (mean_diff, count_small_diffs)
    """
    diffs = []
    for p, c in zip(plaintext, ciphertext):
        p_idx = string.ascii_uppercase.index(p) if p in string.ascii_uppercase else -1
        c_idx = string.ascii_uppercase.index(c) if c in string.ascii_uppercase else -1
        
        if p_idx >= 0 and c_idx >= 0:
            diff = min((c_idx - p_idx) % 26, (p_idx - c_idx) % 26)
            diffs.append(diff)
    
    if not diffs:
        return 13, 0  # Default to middle value and no small diffs
    
    mean_diff = sum(diffs) / len(diffs)
    count_small_diffs = sum(1 for d in diffs if d <= 5)
    
    return mean_diff, count_small_diffs

def try_key_and_alphabets(primer, base, key_method, plain_alphabet, cipher_alphabet, masking_shift=0):
    """
    Try a specific key and alphabet combination.
    
    Args:
        primer (list): Primer to test
        base (int): Number base
        key_method (str): Key expansion method
        plain_alphabet (str): Plaintext alphabet
        cipher_alphabet (str): Ciphertext alphabet
        masking_shift (int): Caesar shift for masking
        
    Returns:
        dict: Results including plaintext, score, and analysis
    """
    # Expand the key
    key = expand_key(primer, base, len(K4_CIPHERTEXT), key_method)
    
    # Decrypt
    decrypted = decrypt_gromark(K4_CIPHERTEXT, key, plain_alphabet, cipher_alphabet)
    
    # Undo masking if applied
    if masking_shift > 0:
        plaintext = undo_masking(decrypted, masking_shift, "caesar")
    else:
        plaintext = decrypted
    
    # Verify known plaintext
    match_percentage, matches, total = verify_known_plaintext(plaintext)
    
    # Only proceed with full scoring if we have a reasonable match
    if match_percentage < 0.3:
        return {
            "primer": primer,
            "base": base,
            "key_method": key_method,
            "plain_alphabet": plain_alphabet[:10] + "..." if len(plain_alphabet) > 10 else plain_alphabet,
            "cipher_alphabet": cipher_alphabet[:10] + "..." if len(cipher_alphabet) > 10 else cipher_alphabet,
            "masking_shift": masking_shift,
            "matches": matches,
            "total": total,
            "match_percentage": match_percentage,
            "plaintext": plaintext[:30] + "..." if len(plaintext) > 30 else plaintext,
            "score": -1000,
            "ic": calculate_ic(plaintext),
            "period": check_key_periodicity(key)
        }
    
    # Full scoring for promising candidates
    score = score_plaintext(plaintext, match_percentage)
    ic = calculate_ic(plaintext)
    period = check_key_periodicity(key)
    
    # Analyze minor differences
    mean_diff, small_diffs = analyze_minor_differences(plaintext, K4_CIPHERTEXT)
    
    return {
        "primer": primer,
        "base": base,
        "key_method": key_method,
        "plain_alphabet": plain_alphabet[:10] + "..." if len(plain_alphabet) > 10 else plain_alphabet,
        "cipher_alphabet": cipher_alphabet[:10] + "..." if len(cipher_alphabet) > 10 else cipher_alphabet,
        "masking_shift": masking_shift,
        "matches": matches,
        "total": total,
        "match_percentage": match_percentage,
        "plaintext": plaintext,
        "score": score,
        "ic": ic,
        "period": period,
        "mean_diff": mean_diff,
        "small_diffs": small_diffs,
        "key": key[:10]  # First 10 digits of the key
    }

def search_promising_primers():
    """
    Search promising primers from Bean's paper.
    
    Returns:
        list: Sorted results
    """
    results = []
    
    # Key expansion methods to try
    key_methods = ["standard", "fibonacci", "alternating", "weighted", "m-sequence"]
    
    # Masking shifts to try (0 = no masking)
    masking_shifts = [0, 1, 2, 3, 5, 7, 13]
    
    # Standard alphabets
    print("Testing promising primers with standard alphabets...")
    for primer, base, desc in PROMISING_PRIMERS:
        print(f"Testing primer {primer} (Base {base}): {desc}")
        
        for key_method in key_methods:
            for masking_shift in masking_shifts:
                result = try_key_and_alphabets(
                    primer, base, key_method,
                    string.ascii_uppercase, string.ascii_uppercase,
                    masking_shift
                )
                
                # Only add promising results
                if result["match_percentage"] >= 0.3:
                    results.append(result)
                    print(f"  Found promising result: Match {result['match_percentage']:.2f}, Score {result['score']:.2f}")
    
    # Keyword-based alphabets (only for most promising primers)
    print("\nTesting promising primers with keyword alphabets...")
    top_primers = [(primer, base) for primer, base, _ in PROMISING_PRIMERS[:3]]
    
    for primer, base in top_primers:
        for key_method in key_methods[:2]:  # Only use standard and fibonacci methods
            for keyword in KEYWORDS:
                # Create keyword alphabets
                mixed_alphabet = create_keyword_alphabet(keyword)
                
                # Try plain=standard, cipher=keyword
                result = try_key_and_alphabets(
                    primer, base, key_method,
                    string.ascii_uppercase, mixed_alphabet,
                    0  # No masking shift when using keyword alphabets
                )
                
                if result["match_percentage"] >= 0.3:
                    results.append(result)
                    print(f"  Found promising result with keyword {keyword}: Match {result['match_percentage']:.2f}")
                
                # Try plain=keyword, cipher=standard
                result = try_key_and_alphabets(
                    primer, base, key_method,
                    mixed_alphabet, string.ascii_uppercase,
                    0
                )
                
                if result["match_percentage"] >= 0.3:
                    results.append(result)
                    print(f"  Found promising result with keyword {keyword}: Match {result['match_percentage']:.2f}")
                
                # Try both alphabets as keywords
                result = try_key_and_alphabets(
                    primer, base, key_method,
                    mixed_alphabet, mixed_alphabet,
                    0
                )
                
                if result["match_percentage"] >= 0.3:
                    results.append(result)
                    print(f"  Found promising result with keyword {keyword}: Match {result['match_percentage']:.2f}")
    
    # Sort results by score
    results.sort(key=lambda x: x["score"], reverse=True)
    
    return results

def test_constraint_satisfied(key):
    """
    Test if a key satisfies the known constraints from plaintext-ciphertext pairs.
    
    Args:
        key (list): Key to test
        
    Returns:
        bool: Whether the key satisfies all constraints
    """
    # Equalities and inequalities derived from known plaintext
    # k27 = k65 (from plaintext 'R' at these positions)
    if key[27] != key[65]:
        return False
    
    # Sample of key inequalities
    inequalities = [
        (24, 28), (28, 33), (24, 33),  # T → different ciphertext
        (21, 30), (21, 64), (30, 64),  # E → different ciphertext
        (22, 31),                     # A → different ciphertext
        (66, 70),                     # L → different ciphertext
        (26, 71),                     # O → different ciphertext
        (69, 72),                     # C → different ciphertext
        (23, 32),                     # S → different ciphertext
        (25, 26),                     # N,O → same ciphertext
        (24, 66),                     # T,L → same ciphertext
        (29, 63),                     # H,B → same ciphertext
        (32, 33),                     # S,T → same ciphertext
        (67, 68),                     # I,N → same ciphertext
        (27, 72)                      # R,C → same ciphertext
    ]
    
    for a, b in inequalities:
        if a < len(key) and b < len(key) and key[a] == key[b]:
            return False
    
    return True

def generate_constrained_primer(base, length):
    """
    Generate a random primer that satisfies the known constraints.
    
    Args:
        base (int): Number base
        length (int): Primer length
        
    Returns:
        list: Generated primer
    """
    max_attempts = 1000
    for _ in range(max_attempts):
        # Generate random primer
        primer = [random.randint(0, base-1) for _ in range(length)]
        
        # Expand to full key
        key = expand_key(primer, base, len(K4_CIPHERTEXT))
        
        # Check constraints
        if test_constraint_satisfied(key):
            return primer
    
    # If we couldn't find a constrained primer, return a random one
    return [random.randint(0, base-1) for _ in range(length)]

class SimulatedAnnealing:
    """Simulated annealing implementation for K4 solution."""
    
    def __init__(self, primer, base, key_method, plain_shift=0, cipher_shift=0):
        """Initialize the annealing process."""
        self.primer = primer.copy()
        self.base = base
        self.key_method = key_method
        self.plain_shift = plain_shift
        self.cipher_shift = cipher_shift
        self.best_score = float('-inf')
        self.best_state = None
        self.best_plaintext = None
        
    def energy(self):
        """Calculate the energy (negative score) of the current state."""
        plain_alphabet = create_caesar_shifted_alphabet(self.plain_shift)
        cipher_alphabet = create_caesar_shifted_alphabet(self.cipher_shift)
        
        key = expand_key(self.primer, self.base, len(K4_CIPHERTEXT), self.key_method)
        decrypted = decrypt_gromark(K4_CIPHERTEXT, key, plain_alphabet, cipher_alphabet)
        
        match_percentage, matches, total = verify_known_plaintext(decrypted)
        score = score_plaintext(decrypted, match_percentage)
        
        # Update best state if found
        if score > self.best_score:
            self.best_score = score
            self.best_state = (self.primer.copy(), self.plain_shift, self.cipher_shift)
            self.best_plaintext = decrypted
            
        return -score  # Negative because we're minimizing energy
        
    def move(self):
        """Make a random move in the state space."""
        choice = random.random()
        
        if choice < 0.7:  # 70% chance to change primer
            idx = random.randint(0, len(self.primer) - 1)
            old_value = self.primer[idx]
            while self.primer[idx] == old_value:
                self.primer[idx] = random.randint(0, self.base - 1)
        elif choice < 0.85:  # 15% chance to change plain shift
            self.plain_shift = (self.plain_shift + random.randint(1, 25)) % 26
        else:  # 15% chance to change cipher shift
            self.cipher_shift = (self.cipher_shift + random.randint(1, 25)) % 26
            
    def anneal(self, steps=10000, temp_start=10.0, temp_end=0.1):
        """Run the simulated annealing process."""
        temp = temp_start
        e_curr = self.energy()
        
        for step in range(steps):
            # Save old state
            old_primer = self.primer.copy()
            old_plain_shift = self.plain_shift
            old_cipher_shift = self.cipher_shift
            
            # Make a move
            self.move()
            
            # Calculate new energy
            e_new = self.energy()
            
            # Accept or reject the move
            temp = temp_start * (temp_end / temp_start) ** (step / steps)
            delta_e = e_new - e_curr
            
            if delta_e <= 0 or random.random() < math.exp(-delta_e / temp):
                # Accept the move
                e_curr = e_new
            else:
                # Reject the move, restore old state
                self.primer = old_primer
                self.plain_shift = old_plain_shift
                self.cipher_shift = old_cipher_shift
                
            # Print progress every 1000 steps
            if (step + 1) % 1000 == 0:
                print(f"Step {step+1}/{steps}, Temp: {temp:.4f}, Best score: {self.best_score:.2f}")
                
        return self.best_state, self.best_score, self.best_plaintext

def run_simulated_annealing():
    """
    Run simulated annealing to find optimal parameters.
    
    Returns:
        list: Results from annealing
    """
    results = []
    
    print("\nRunning simulated annealing optimization...")
    
    # Bases and primer lengths to try
    base_lengths = [(10, 5), (8, 5), (5, 5), (12, 5), (10, 4)]
    key_methods = ["standard", "fibonacci"]
    
    for base, length in base_lengths:
        for key_method in key_methods:
            print(f"\nOptimizing for base {base}, length {length}, method {key_method}")
            
            # Generate constrained primer
            primer = generate_constrained_primer(base, length)
            print(f"Starting with primer: {primer}")
            
            # Run annealing
            annealer = SimulatedAnnealing(primer, base, key_method)
            best_state, best_score, best_plaintext = annealer.anneal(steps=5000)
            
            if best_state:
                best_primer, plain_shift, cipher_shift = best_state
                
                # Calculate IC and other metrics
                ic = calculate_ic(best_plaintext)
                match_percentage, matches, total = verify_known_plaintext(best_plaintext)
                
                result = {
                    "primer": best_primer,
                    "base": base,
                    "key_method": key_method,
                    "plain_shift": plain_shift,
                    "cipher_shift": cipher_shift,
                    "plaintext": best_plaintext,
                    "score": -best_score,  # Convert back to positive
                    "ic": ic,
                    "match_percentage": match_percentage,
                    "matches": matches,
                    "total": total
                }
                
                results.append(result)
                
                print(f"Best result for this configuration:")
                print(f"  Primer: {best_primer}")
                print(f"  Score: {-best_score:.2f}, IC: {ic:.6f}")
                print(f"  Matches: {matches}/{total} ({match_percentage:.2f})")
                print(f"  Plaintext: {best_plaintext[:50]}...")
    
    # Sort results by score
    results.sort(key=lambda x: x["score"], reverse=True)
    
    return results

def print_results_at_width21(plaintext):
    """
    Print plaintext arranged at width 21, as in the Kryptos sculpture.
    
    Args:
        plaintext (str): Plaintext to display
    """
    for i in range(0, len(plaintext), 21):
        print(plaintext[i:i+21])

def main():
    """Main function to run the K4 solver."""
    print("KRYPTOS K4 COMPREHENSIVE SOLVER")
    print("=" * 50)
    print(f"Ciphertext: {K4_CIPHERTEXT}")
    print(f"Length: {len(K4_CIPHERTEXT)} characters")
    print("Known plaintext positions:")
    for pos, (plain, cipher) in KNOWN_PLAINTEXT.items():
        print(f"  Position {pos}: '{plain}' → '{cipher}'")
    print("=" * 50)
    
    start_time = time.time()
    
    # Step 1: Initial search with promising primers
    print("\nSTEP 1: SEARCHING PROMISING PRIMERS")
    print("-" * 50)
    promising_results = search_promising_primers()
    
    # Step 2: Run simulated annealing
    print("\nSTEP 2: RUNNING SIMULATED ANNEALING")
    print("-" * 50)
    annealing_results = run_simulated_annealing()
    
    # Combine results
    all_results = promising_results + annealing_results
    all_results.sort(key=lambda x: x["score"], reverse=True)
    
    # Print top results
    print("\nTOP RESULTS")
    print("=" * 100)
    for i, result in enumerate(all_results[:10]):
        print(f"\nRank {i+1}:")
        print(f"  Primer: {result['primer']}")
        print(f"  Base: {result['base']}, Method: {result['key_method']}")
        
        if "plain_shift" in result:
            print(f"  Plain shift: {result['plain_shift']}, Cipher shift: {result['cipher_shift']}")
        elif "masking_shift" in result:
            print(f"  Masking shift: {result['masking_shift']}")
            
        print(f"  Plain alphabet: {result['plain_alphabet'] if 'plain_alphabet' in result else 'Standard'}")
        print(f"  Cipher alphabet: {result['cipher_alphabet'] if 'cipher_alphabet' in result else 'Standard'}")
        print(f"  Score: {result['score']:.2f}, IC: {result['ic']:.6f}")
        print(f"  Matches: {result['matches']}/{result['total']} ({result['match_percentage']*100:.2f}%)")
        
        if "period" in result and result["period"]:
            print(f"  Key periodicity: {result['period']}")
        
        if "mean_diff" in result:
            print(f"  Mean letter difference: {result['mean_diff']:.2f}, Small diffs: {result['small_diffs']}")
        
        print(f"  Plaintext: {result['plaintext'][:50]}...")
        
        print("\n  Plaintext at width 21:")
        print_results_at_width21(result['plaintext'])
        
        print("-" * 100)
    
    # Print execution time
    end_time = time.time()
    print(f"\nExecution time: {end_time - start_time:.2f} seconds")

if __name__ == "__main__":
    main() 