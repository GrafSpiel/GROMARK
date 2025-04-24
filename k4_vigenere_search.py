#!/usr/bin/env python3
"""
Kryptos K4 Vigenere Wordlist Search
Implements a brute-force search using wordlist keys with various Vigenere cipher variants.
Focuses on finding adjacent character matches at known positions.
"""

import string
import time
import argparse
import multiprocessing
from collections import Counter, defaultdict

# K4 ciphertext
K4_CIPHERTEXT = "OBKRUOXOGHULBSOLIFBBWFLRVQQPRNGKSSOTWTQSJQSSEKZZWATJKLUDIAWINFBNYPVTTMZFPKWGDKZXTJCDIGKUHUAUEKCAR"

# Target plaintext with known fragments
TARGET_PLAINTEXT = "?????????????????????" + "EASTNORTHEAST" + "???????????????????????" + "BERLINCLOCK" + "??????????????????"

# Make sure the target plaintext length matches ciphertext
if len(TARGET_PLAINTEXT) != len(K4_CIPHERTEXT):
    # Adjust target plaintext if needed
    TARGET_PLAINTEXT = TARGET_PLAINTEXT[:len(K4_CIPHERTEXT)].ljust(len(K4_CIPHERTEXT), '?')

def load_wordlist(filename, min_length=3, max_length=15):
    """Load words from a file."""
    words = []
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            for line in f:
                word = line.strip().upper()
                if min_length <= len(word) <= max_length and all(c in string.ascii_uppercase for c in word):
                    words.append(word)
        print(f"Loaded {len(words)} words from {filename}")
    except Exception as e:
        print(f"Error loading wordlist: {e}")
    
    return words

def vigenere_decrypt(ciphertext, key):
    """Standard Vigenere decryption."""
    plaintext = []
    key_length = len(key)
    
    for i, char in enumerate(ciphertext):
        if char in string.ascii_uppercase:
            # Convert both to 0-25 range
            char_idx = ord(char) - ord('A')
            key_char = key[i % key_length]
            key_idx = ord(key_char) - ord('A')
            
            # Decrypt: (ciphertext_char - key_char) mod 26
            plain_idx = (char_idx - key_idx) % 26
            plaintext.append(chr(plain_idx + ord('A')))
        else:
            plaintext.append(char)
    
    return ''.join(plaintext)

def vigenere_periodic_rotate_forward_decrypt(ciphertext, key):
    """Vigenere with key that rotates forward after each use."""
    plaintext = []
    key_length = len(key)
    rotated_key = key
    
    for i, char in enumerate(ciphertext):
        if char in string.ascii_uppercase:
            # Use the current rotated key
            key_char = rotated_key[i % key_length]
            char_idx = ord(char) - ord('A')
            key_idx = ord(key_char) - ord('A')
            
            # Decrypt: (ciphertext_char - key_char) mod 26
            plain_idx = (char_idx - key_idx) % 26
            plaintext.append(chr(plain_idx + ord('A')))
            
            # Rotate the key forward by 1 position after using it
            if (i + 1) % key_length == 0:
                rotated_key = rotated_key[1:] + rotated_key[0]
        else:
            plaintext.append(char)
    
    return ''.join(plaintext)

def vigenere_periodic_rotate_reverse_decrypt(ciphertext, key):
    """Vigenere with key that rotates backward after each use."""
    plaintext = []
    key_length = len(key)
    rotated_key = key
    
    for i, char in enumerate(ciphertext):
        if char in string.ascii_uppercase:
            # Use the current rotated key
            key_char = rotated_key[i % key_length]
            char_idx = ord(char) - ord('A')
            key_idx = ord(key_char) - ord('A')
            
            # Decrypt: (ciphertext_char - key_char) mod 26
            plain_idx = (char_idx - key_idx) % 26
            plaintext.append(chr(plain_idx + ord('A')))
            
            # Rotate the key backward by 1 position after using it
            if (i + 1) % key_length == 0:
                rotated_key = rotated_key[-1] + rotated_key[:-1]
        else:
            plaintext.append(char)
    
    return ''.join(plaintext)

def calculate_ic(text):
    """Calculate the Index of Coincidence for a text."""
    counts = Counter(c for c in text.upper() if c in string.ascii_uppercase)
    n = sum(counts.values())
    if n <= 1:
        return 0
    
    sum_freqs = sum(count * (count - 1) for count in counts.values())
    return sum_freqs / (n * (n - 1))

def check_adjacent_matches(plaintext, target=TARGET_PLAINTEXT, min_adjacent=4):
    """
    Check for adjacent character matches at known positions.
    Return (max_consecutive_matches, matched_positions)
    """
    max_consecutive = 0
    current_consecutive = 0
    all_match_positions = []
    
    for i in range(min(len(plaintext), len(target))):
        if target[i] != '?' and plaintext[i] == target[i]:
            current_consecutive += 1
            all_match_positions.append(i)
        else:
            # End of a run or no match
            max_consecutive = max(max_consecutive, current_consecutive)
            current_consecutive = 0
    
    # Final check in case the longest run is at the end
    max_consecutive = max(max_consecutive, current_consecutive)
    
    return max_consecutive, all_match_positions

def format_match_ranges(positions):
    """Convert a list of positions to a concise range representation."""
    if not positions:
        return "None"
    
    ranges = []
    start = positions[0]
    prev = start
    
    for pos in positions[1:]:
        if pos != prev + 1:
            if start == prev:
                ranges.append(f"{start}")
            else:
                ranges.append(f"{start}-{prev}")
            start = pos
        prev = pos
    
    # Add the last range
    if start == prev:
        ranges.append(f"{start}")
    else:
        ranges.append(f"{start}-{prev}")
    
    return ", ".join(ranges)

def test_key(args):
    """Test a key with multiple Vigenere variants."""
    word, min_adjacent = args
    results = []
    
    # Dictionary of cipher functions
    cipher_functions = {
        "Standard Vigenere": vigenere_decrypt,
        "Periodic Forward Rotating": vigenere_periodic_rotate_forward_decrypt,
        "Periodic Reverse Rotating": vigenere_periodic_rotate_reverse_decrypt
    }
    
    for cipher_name, cipher_func in cipher_functions.items():
        plaintext = cipher_func(K4_CIPHERTEXT, word)
        max_adjacent, match_positions = check_adjacent_matches(plaintext, TARGET_PLAINTEXT, min_adjacent)
        ic = calculate_ic(plaintext)
        
        if max_adjacent >= min_adjacent:
            results.append({
                "word": word,
                "cipher": cipher_name,
                "max_adjacent": max_adjacent,
                "match_positions": match_positions,
                "plaintext": plaintext,
                "ic": ic
            })
    
    return results

def main():
    """Main function to run the Vigenere wordlist search."""
    parser = argparse.ArgumentParser(description="K4 Vigenere Wordlist Search")
    parser.add_argument("--wordlist", default="not.txt", help="Path to wordlist file")
    parser.add_argument("--min-length", type=int, default=3, help="Minimum word length")
    parser.add_argument("--max-length", type=int, default=15, help="Maximum word length")
    parser.add_argument("--processes", type=int, default=None, help="Number of processes to use")
    parser.add_argument("--min-adjacent", type=int, default=4, help="Minimum adjacent matches required")
    parser.add_argument("--words", type=str, help="Comma-separated list of specific words to test")
    parser.add_argument("--output", default="vigenere_results.txt", help="Output file for results")
    parser.add_argument("--batch-size", type=int, default=100, help="Batch size for processing")
    
    args = parser.parse_args()
    
    print("K4 VIGENERE WORDLIST SEARCH")
    print("=" * 50)
    
    # Set up specific words if requested
    word_filter = None
    if args.words:
        word_filter = [w.strip().upper() for w in args.words.split(",")]
        print(f"Testing specific words: {word_filter}")
    
    # Load wordlist
    words = []
    if word_filter:
        words = word_filter
    else:
        words = load_wordlist(args.wordlist, args.min_length, args.max_length)
    
    if not words:
        print("No words to test. Exiting.")
        return
    
    # Setup for multiprocessing
    num_processes = args.processes if args.processes else max(1, multiprocessing.cpu_count() - 1)
    print(f"Using {num_processes} processes")
    
    # Statistics
    all_results = []
    start_time = time.time()
    
    # Process words in batches
    batch_size = min(args.batch_size, len(words))
    num_batches = (len(words) + batch_size - 1) // batch_size
    
    print(f"\nTesting {len(words)} words with 3 Vigenere variants")
    print(f"Processing in {num_batches} batches of up to {batch_size} words")
    
    for batch_idx in range(num_batches):
        batch_start = batch_idx * batch_size
        batch_end = min(batch_start + batch_size, len(words))
        batch_words = words[batch_start:batch_end]
        
        batch_start_time = time.time()
        print(f"\nBatch {batch_idx+1}/{num_batches}: Testing words {batch_start+1}-{batch_end}")
        
        # Prepare tasks for multiprocessing
        tasks = [(word, args.min_adjacent) for word in batch_words]
        
        # Process batch using multiprocessing
        with multiprocessing.Pool(processes=num_processes) as pool:
            batch_results = pool.map(test_key, tasks)
        
        # Flatten results
        batch_results = [item for sublist in batch_results for item in sublist]
        
        # Add to all results
        all_results.extend(batch_results)
        
        # Sort by max_adjacent and then by IC
        all_results.sort(key=lambda x: (x["max_adjacent"], x["ic"]), reverse=True)
        
        batch_end_time = time.time()
        batch_time = batch_end_time - batch_start_time
        
        # Print batch summary
        print(f"Batch completed in {batch_time:.2f} seconds")
        print(f"Found {len(batch_results)} results with {args.min_adjacent}+ adjacent matches in this batch")
        
        # Show best result in this batch
        if batch_results:
            best = max(batch_results, key=lambda x: (x["max_adjacent"], x["ic"]))
            print(f"Best in batch: {best['word']} + {best['cipher']} = {best['max_adjacent']} adjacent matches, IC={best['ic']:.6f}")
    
    end_time = time.time()
    print(f"\nSearch completed in {end_time - start_time:.2f} seconds")
    print(f"Found {len(all_results)} results with {args.min_adjacent}+ adjacent matches")
    
    # Print top results
    print("\nTOP RESULTS")
    print("=" * 100)
    
    for i, result in enumerate(all_results[:20]):  # Show top 20
        print(f"\nRank {i+1}:")
        print(f"  Word: {result['word']}")
        print(f"  Cipher: {result['cipher']}")
        print(f"  Max Adjacent Matches: {result['max_adjacent']}")
        print(f"  Match Positions: {format_match_ranges(result['match_positions'])}")
        print(f"  IC: {result['ic']:.6f}")
        print(f"  Plaintext: {result['plaintext'][:50]}...")
        print("-" * 100)
    
    # Save results
    with open(args.output, 'w', encoding='utf-8') as f:
        f.write(f"K4 Vigenere Wordlist Search Results\n")
        f.write(f"Total results: {len(all_results)}\n")
        f.write("=" * 80 + "\n\n")
        
        for i, result in enumerate(all_results[:200]):  # Save top 200 results
            f.write(f"Rank {i+1}:\n")
            f.write(f"  Word: {result['word']}\n")
            f.write(f"  Cipher: {result['cipher']}\n")
            f.write(f"  Max Adjacent Matches: {result['max_adjacent']}\n")
            f.write(f"  Match Positions: {format_match_ranges(result['match_positions'])}\n")
            f.write(f"  IC: {result['ic']:.6f}\n")
            f.write(f"  Plaintext: {result['plaintext']}\n\n")
    
    print(f"Results saved to {args.output}")

if __name__ == "__main__":
    main() 