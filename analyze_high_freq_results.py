#!/usr/bin/env python3
"""
Analyze K4 Gromark results with high frequency scores.
This script focuses on applying various transposition methods to promising decryption candidates.
"""

import os
import argparse
import time
from k4_transpositions import process_promising_results

def main():
    parser = argparse.ArgumentParser(description='Analyze high frequency score results with advanced transpositions')
    parser.add_argument('--input', type=str, default='results/all_results_latest.txt', help='Input results file')
    parser.add_argument('--output', type=str, default='results/transposition_results.txt', help='Output file')
    parser.add_argument('--min-freq-score', type=float, default=15.0, help='Minimum frequency score to consider')
    parser.add_argument('--max-results', type=int, default=50, help='Maximum number of results to process')
    
    args = parser.parse_args()
    
    # Find the most recent results file if not specified
    if args.input == 'results/all_results_latest.txt':
        result_files = []
        
        # Check if directory exists
        if os.path.exists('results'):
            for file in os.listdir('results'):
                if file.startswith('all_results_') and file.endswith('.txt'):
                    result_files.append(os.path.join('results', file))
        
        if result_files:
            # Get the most recent file
            latest_file = max(result_files, key=os.path.getmtime)
            args.input = latest_file
            print(f"Using most recent results file: {latest_file}")
        else:
            print("No results files found in 'results' directory. Please specify an input file.")
            return
    
    print(f"Analyzing results with frequency score >= {args.min_freq_score}")
    print(f"Input file: {args.input}")
    print(f"Output file: {args.output}")
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Measure execution time
    start_time = time.time()
    
    # Process the results
    best_results = process_promising_results(
        args.input, 
        args.output, 
        args.min_freq_score, 
        args.max_results
    )
    
    # Print execution time
    execution_time = time.time() - start_time
    print(f"Execution time: {execution_time:.2f} seconds")
    
    # Print top 3 results
    print("\nTop 3 Results after Transposition Analysis:")
    print("="*80)
    
    for i, result in enumerate(best_results[:3]):
        print(f"Rank {i+1}")
        print(f"Original Keyword: {result['original_keyword']}")
        print(f"Original Primer: {result['original_primer']}")
        print(f"Transposition: {result['transposition']} ({result['trans_param']})")
        print(f"IC: {result['ic']:.6f}")
        print(f"Text: {result['text']}")
        print("-"*80)

if __name__ == "__main__":
    main() 