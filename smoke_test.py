#!/usr/bin/env python3
"""
Smoke test for the Gromark cipher solver.
This script verifies that the CPU and GPU implementations produce consistent results
before launching a full search.
"""

import argparse
import time
import string
import random
import numpy as np
from gromark_solver import (
    decrypt_gromark,
    expand_gromark_primer,
    check_cribs,
    K4_CIPHERTEXT,
    test_primer
)

try:
    from gromark_gpu import (
        get_primer_from_index,
        expand_key_and_decrypt
    )
    GPU_AVAILABLE = True
except ImportError:
    GPU_AVAILABLE = False
    print("Warning: GPU functions not available. Only CPU tests will be run.")

def smoke_test_cpu(known_primers, base=10, size=1000, ciphertext=K4_CIPHERTEXT):
    """
    Test the CPU implementation with known primers and a sample of random primers.
    
    Args:
        known_primers (list): List of known primers that should be tested
        base (int): Number base to use
        size (int): Number of random primers to test
        ciphertext (str): Ciphertext to test with
    
    Returns:
        set: Set of indices of primers that match the cribs
    """
    print("\n== CPU Smoke Test ==")
    print(f"Testing known primers: {known_primers}")
    
    cipher_alphabet = string.ascii_uppercase
    plain_alphabet = string.ascii_uppercase
    expansion_func = expand_gromark_primer
    min_match_percentage = 0.7
    
    # Test known primers
    matching_primers = set()
    for primer in known_primers:
        result = test_primer((primer, base, ciphertext, plain_alphabet, cipher_alphabet, expansion_func, min_match_percentage))
        if result:
            primer, plaintext, score, ic = result
            print(f"MATCH: Primer {primer} -> {plaintext[:50]}...")
            matching_primers.add(tuple(primer))
        else:
            print(f"NO MATCH: Primer {primer}")
    
    # Test random sample of primers
    print(f"\nTesting {size} random primers...")
    start_time = time.time()
    random_matches = 0
    
    # Generate random primers
    for _ in range(size):
        primer = [random.randint(0, base-1) for _ in range(len(known_primers[0]))]
        result = test_primer((primer, base, ciphertext, plain_alphabet, cipher_alphabet, expansion_func, min_match_percentage))
        if result:
            primer, plaintext, score, ic = result
            print(f"RANDOM MATCH: Primer {primer} -> {plaintext[:50]}...")
            matching_primers.add(tuple(primer))
            random_matches += 1
    
    end_time = time.time()
    print(f"CPU test completed in {end_time - start_time:.2f} seconds")
    print(f"Found {len(matching_primers)} matching primers total, {random_matches} from random sample")
    
    return matching_primers

def smoke_test_gpu(known_primers, base=10, size=1000, ciphertext=K4_CIPHERTEXT):
    """
    Test the GPU implementation with known primers and a sample of random primers.
    
    Args:
        known_primers (list): List of known primers that should be tested
        base (int): Number base to use
        size (int): Number of random primers to test
        ciphertext (str): Ciphertext to test with
    
    Returns:
        set: Set of indices of primers that match the cribs
    """
    if not GPU_AVAILABLE:
        print("GPU testing not available")
        return set()
    
    from gromark_gpu import gpu_brute_force_gromark
    
    print("\n== GPU Smoke Test ==")
    print(f"Testing known primers: {known_primers}")
    
    # Test known primers
    cipher_alphabet = string.ascii_uppercase
    
    # Run GPU test with known primers
    result_indices = gpu_brute_force_gromark(
        base=base,
        primer_length=len(known_primers[0]),
        ciphertext=ciphertext,
        cipher_alphabet=cipher_alphabet,
        specific_primers=known_primers
    )
    
    matching_primers = set()
    for primer_or_none, idx in result_indices:
        if primer_or_none is not None:  # This is from specific_primers
            matching_primers.add(tuple(primer_or_none))
    
    # Generate random primers for testing
    print(f"\nTesting {size} random primers...")
    
    random_primers = []
    for _ in range(size):
        primer = [random.randint(0, base-1) for _ in range(len(known_primers[0]))]
        random_primers.append(primer)
    
    # Run GPU test with random primers
    start_time = time.time()
    result_indices = gpu_brute_force_gromark(
        base=base,
        primer_length=len(known_primers[0]),
        ciphertext=ciphertext,
        cipher_alphabet=cipher_alphabet,
        specific_primers=random_primers
    )
    
    random_matches = 0
    for primer_or_none, idx in result_indices:
        if primer_or_none is not None:  # This is from specific_primers
            print(f"RANDOM MATCH: Primer {primer_or_none}")
            matching_primers.add(tuple(primer_or_none))
            random_matches += 1
    
    end_time = time.time()
    print(f"GPU test completed in {end_time - start_time:.2f} seconds")
    print(f"Found {len(matching_primers)} matching primers total, {random_matches} from random sample")
    
    return matching_primers

def compare_results(cpu_matches, gpu_matches):
    """
    Compare the results from CPU and GPU implementations.
    
    Args:
        cpu_matches (set): Set of matching primers from CPU test
        gpu_matches (set): Set of matching primers from GPU test
    
    Returns:
        bool: True if results match, False otherwise
    """
    if not GPU_AVAILABLE:
        print("\nGPU not available, skipping comparison")
        return True
    
    print("\n== Comparison of CPU and GPU Results ==")
    
    if cpu_matches == gpu_matches:
        print("SUCCESS: CPU and GPU implementations return the same results")
        return True
    else:
        print("ERROR: CPU and GPU implementations returned different results")
        
        # Show differences
        cpu_only = cpu_matches - gpu_matches
        gpu_only = gpu_matches - cpu_matches
        
        if cpu_only:
            print(f"Primers found by CPU but not GPU: {cpu_only}")
        if gpu_only:
            print(f"Primers found by GPU but not CPU: {gpu_only}")
        
        return False

def main():
    """Main function to run the smoke test."""
    parser = argparse.ArgumentParser(description="Smoke test for Gromark cipher solver")
    parser.add_argument("--base", type=int, default=10, help="Number base to use (default: 10)")
    parser.add_argument("--size", type=int, default=1000, help="Number of random primers to test (default: 1000)")
    parser.add_argument("--cpu-only", action="store_true", help="Run only CPU tests")
    
    args = parser.parse_args()
    
    # Known primers that should be tested
    known_primers = [
        [2, 6, 7, 1, 7],  # The two equivalent primers from Bean's paper
        [8, 4, 3, 9, 3],
        [9, 8, 8, 0, 0]   # The primer with best IC value
    ]
    
    # Run CPU tests
    cpu_matches = smoke_test_cpu(known_primers, args.base, args.size)
    
    # Run GPU tests if available and not disabled
    gpu_matches = set()
    if GPU_AVAILABLE and not args.cpu_only:
        gpu_matches = smoke_test_gpu(known_primers, args.base, args.size)
        
        # Compare results
        match = compare_results(cpu_matches, gpu_matches)
        
        if match:
            print("\nSMOKE TEST PASSED: CPU and GPU implementations are consistent")
        else:
            print("\nSMOKE TEST FAILED: CPU and GPU implementations are not consistent")
    else:
        print("\nSMOKE TEST COMPLETED: CPU only")

if __name__ == "__main__":
    main() 