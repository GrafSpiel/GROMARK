#!/usr/bin/env python3
"""
K4 Advanced Transposition Methods

This module provides various transposition techniques to apply to promising decryption
results from the Gromark analysis. These techniques include:

1. Keyword-based Columnar Transposition
2. Double Columnar Transposition
3. Myszkowski Transposition
4. Route Cipher
5. Rail Fence Cipher
6. Scytale (Cylinder) Transposition

Each technique includes both encryption and decryption methods for testing.
"""

import string
import math
import random
import itertools
from collections import defaultdict

def create_columnar_key_order(keyword):
    """
    Create a column order based on keyword alphabetical ordering.
    Returns a list of indices representing the order to read columns.
    """
    # Convert keyword to uppercase
    keyword = keyword.upper()
    
    # Create a list of (character, position) tuples
    char_positions = [(char, i) for i, char in enumerate(keyword)]
    
    # Sort by character (alphabetically)
    sorted_positions = sorted(char_positions)
    
    # Extract the column order
    column_order = [pos for _, pos in sorted_positions]
    
    return column_order

def columnar_transposition_encrypt(text, keyword):
    """
    Encrypt using columnar transposition with a keyword.
    
    Args:
        text (str): Text to encrypt
        keyword (str): Keyword to determine column order
        
    Returns:
        str: Encrypted text
    """
    # Get column order from keyword
    col_order = create_columnar_key_order(keyword)
    num_cols = len(keyword)
    
    # Calculate number of rows needed
    num_rows = math.ceil(len(text) / num_cols)
    
    # Pad the text if necessary
    padded_text = text + 'X' * (num_rows * num_cols - len(text))
    
    # Create grid (row by row)
    grid = []
    for i in range(0, len(padded_text), num_cols):
        grid.append(padded_text[i:i+num_cols])
    
    # Read out by column in keyword order
    result = ''
    for col_idx in col_order:
        for row_idx in range(num_rows):
            if col_idx < len(grid[row_idx]):
                result += grid[row_idx][col_idx]
    
    return result

def columnar_transposition_decrypt(text, keyword):
    """
    Decrypt using columnar transposition with a keyword.
    
    Args:
        text (str): Text to decrypt
        keyword (str): Keyword to determine column order
        
    Returns:
        str: Decrypted text
    """
    # Get column order from keyword
    col_order = create_columnar_key_order(keyword)
    num_cols = len(keyword)
    
    # Calculate number of rows needed
    num_rows = math.ceil(len(text) / num_cols)
    
    # Calculate lengths of columns (last row might be incomplete)
    last_row_cols = len(text) % num_cols
    if last_row_cols == 0:
        last_row_cols = num_cols
    
    # Create empty grid
    grid = [[''] * num_cols for _ in range(num_rows)]
    
    # Fill in the grid by columns in the given order
    idx = 0
    for col_pos in range(num_cols):
        # Find the actual column index from the order
        col_idx = col_order.index(col_pos)
        
        # Determine number of rows in this column
        rows_in_col = num_rows if col_idx < last_row_cols or last_row_cols == 0 else num_rows - 1
        
        # Fill the column
        for row_idx in range(rows_in_col):
            grid[row_idx][col_idx] = text[idx]
            idx += 1
    
    # Read out by rows
    result = ''
    for row in grid:
        result += ''.join(row)
    
    return result

def double_columnar_transposition(text, keyword1, keyword2=None):
    """
    Apply double columnar transposition.
    
    Args:
        text (str): Text to transform
        keyword1 (str): First keyword
        keyword2 (str, optional): Second keyword. If None, use the same as first.
        
    Returns:
        str: Transformed text
    """
    if not keyword2:
        keyword2 = keyword1
    
    # Apply first transposition
    intermediate = columnar_transposition_encrypt(text, keyword1)
    
    # Apply second transposition
    result = columnar_transposition_encrypt(intermediate, keyword2)
    
    return result

def myszkowski_transposition_encrypt(text, keyword):
    """
    Encrypt using Myszkowski transposition.
    
    In this variant, repeated letters in the keyword create columns that
    are read in same order, but consecutively.
    
    Args:
        text (str): Text to encrypt
        keyword (str): Keyword with possible repeated letters
        
    Returns:
        str: Encrypted text
    """
    # Create mapping of keyword characters to positions
    char_positions = defaultdict(list)
    for i, char in enumerate(keyword.upper()):
        char_positions[char].append(i)
    
    # Create priority order for characters (alphabetical)
    priority_chars = sorted(char_positions.keys())
    
    # Calculate number of rows
    num_cols = len(keyword)
    num_rows = math.ceil(len(text) / num_cols)
    
    # Pad the text if necessary
    padded_text = text + 'X' * (num_rows * num_cols - len(text))
    
    # Create the grid
    grid = []
    for i in range(0, len(padded_text), num_cols):
        grid.append(list(padded_text[i:i+num_cols]))
    
    # Read out by columns according to Myszkowski rules
    result = ''
    for char in priority_chars:
        for position in char_positions[char]:
            for row in range(min(num_rows, len(grid))):
                if position < len(grid[row]):
                    result += grid[row][position]
    
    return result

def rail_fence_encrypt(text, rails):
    """
    Encrypt using Rail Fence cipher.
    
    Args:
        text (str): Text to encrypt
        rails (int): Number of rails to use
        
    Returns:
        str: Encrypted text
    """
    # Create the fence pattern
    fence = [[] for _ in range(rails)]
    rail = 0
    direction = 1  # 1 for down, -1 for up
    
    # Place each character in the fence
    for char in text:
        fence[rail].append(char)
        
        # Change direction if at the top or bottom rail
        if rail == 0:
            direction = 1
        elif rail == rails - 1:
            direction = -1
        
        # Move to next rail
        rail += direction
    
    # Read off the fence
    result = ''
    for rail_chars in fence:
        result += ''.join(rail_chars)
    
    return result

def rail_fence_decrypt(text, rails):
    """
    Decrypt using Rail Fence cipher.
    
    Args:
        text (str): Text to decrypt
        rails (int): Number of rails used in encryption
        
    Returns:
        str: Decrypted text
    """
    # Calculate the pattern of the fence
    fence_pattern = []
    rail = 0
    direction = 1
    
    for i in range(len(text)):
        fence_pattern.append(rail)
        if rail == 0:
            direction = 1
        elif rail == rails - 1:
            direction = -1
        rail += direction
    
    # Create empty fence
    fence = [[] for _ in range(rails)]
    
    # Sort pattern and text by rail
    sorted_pairs = sorted(zip(fence_pattern, range(len(text))))
    
    # Distribute text to fence based on sorted order
    idx = 0
    for _, pos in sorted_pairs:
        fence_pattern[pos] = text[idx]
        idx += 1
    
    # Read off in zigzag pattern
    result = ''.join(fence_pattern)
    return result

def route_cipher(text, rows, columns, pattern="spiral"):
    """
    Apply a route cipher with various patterns.
    
    Args:
        text (str): Text to transform
        rows (int): Number of rows in the grid
        columns (int): Number of columns in the grid
        pattern (str): Pattern to follow ("spiral", "snake", "diagonal")
        
    Returns:
        str: Transformed text
    """
    # Pad the text if necessary
    padded_text = text + 'X' * (rows * columns - len(text))
    
    # Create the grid by filling it row by row
    grid = []
    for i in range(0, len(padded_text), columns):
        grid.append(list(padded_text[i:i+columns]))
    
    result = ''
    
    if pattern == "spiral":
        # Spiral pattern (clockwise from outside)
        top, bottom = 0, rows - 1
        left, right = 0, columns - 1
        
        while top <= bottom and left <= right:
            # Top row
            for i in range(left, right + 1):
                if top < len(grid) and i < len(grid[top]):
                    result += grid[top][i]
            top += 1
            
            # Right column
            for i in range(top, bottom + 1):
                if i < len(grid) and right < len(grid[i]):
                    result += grid[i][right]
            right -= 1
            
            # Bottom row
            if top <= bottom:
                for i in range(right, left - 1, -1):
                    if bottom < len(grid) and i < len(grid[bottom]):
                        result += grid[bottom][i]
                bottom -= 1
            
            # Left column
            if left <= right:
                for i in range(bottom, top - 1, -1):
                    if i < len(grid) and left < len(grid[i]):
                        result += grid[i][left]
                left += 1
    
    elif pattern == "snake":
        # Snake pattern (alternating left-to-right and right-to-left)
        for i in range(rows):
            if i % 2 == 0:  # Left to right
                for j in range(columns):
                    if i < len(grid) and j < len(grid[i]):
                        result += grid[i][j]
            else:  # Right to left
                for j in range(columns - 1, -1, -1):
                    if i < len(grid) and j < len(grid[i]):
                        result += grid[i][j]
    
    elif pattern == "diagonal":
        # Diagonal pattern (top-left to bottom-right, then next diagonal)
        for s in range(rows + columns - 1):
            for i in range(rows):
                j = s - i
                if 0 <= j < columns and i < len(grid) and j < len(grid[i]):
                    result += grid[i][j]
    
    return result

def scytale_encrypt(text, diameter):
    """
    Encrypt using a Scytale (cylinder) cipher.
    
    Args:
        text (str): Text to encrypt
        diameter (int): Diameter (circumference) of the cylinder
        
    Returns:
        str: Encrypted text
    """
    # Calculate number of rows needed
    rows = math.ceil(len(text) / diameter)
    
    # Pad the text if necessary
    padded_text = text + 'X' * (rows * diameter - len(text))
    
    # Create the grid (row by row)
    grid = []
    for i in range(0, len(padded_text), diameter):
        grid.append(padded_text[i:i+diameter])
    
    # Read off by columns
    result = ''
    for col in range(diameter):
        for row in range(rows):
            if col < len(grid[row]):
                result += grid[row][col]
    
    return result

def scytale_decrypt(text, diameter):
    """
    Decrypt using a Scytale (cylinder) cipher.
    
    Args:
        text (str): Text to decrypt
        diameter (int): Diameter (circumference) of the cylinder
        
    Returns:
        str: Decrypted text
    """
    # Calculate number of rows
    rows = math.ceil(len(text) / diameter)
    
    # Create the empty grid
    grid = [[''] * diameter for _ in range(rows)]
    
    # Fill the grid by columns
    idx = 0
    for col in range(diameter):
        for row in range(rows):
            if idx < len(text):
                grid[row][col] = text[idx]
                idx += 1
    
    # Read off by rows
    result = ''
    for row in grid:
        result += ''.join(row)
    
    return result

def turning_grille(text, mask_size=4, rotation_steps=4):
    """
    Apply a turning grille cipher.
    This is a simplified implementation where we create a random mask.
    
    Args:
        text (str): Text to transform
        mask_size (int): Size of the grille (4 for 4x4, etc.)
        rotation_steps (int): Number of rotation steps (usually 4 for 90° rotations)
        
    Returns:
        tuple: (transformed text, mask used)
    """
    # Determine total cells in mask
    total_cells = mask_size * mask_size
    
    # Ensure text is long enough
    if len(text) < total_cells:
        padded_text = text + 'X' * (total_cells - len(text))
    else:
        padded_text = text[:total_cells]
    
    # Create the grid
    grid = [[''] * mask_size for _ in range(mask_size)]
    
    # Create a mask (which cells are revealed in each rotation)
    mask = [[False] * mask_size for _ in range(mask_size)]
    
    # Each rotation should reveal exactly 1/4 of the cells (for 4 rotations)
    cells_per_rotation = total_cells // rotation_steps
    
    # Track which positions have been used
    used_positions = set()
    
    # For each rotation
    for rot in range(rotation_steps):
        # Add cells_per_rotation positions for this rotation
        count = 0
        while count < cells_per_rotation:
            row = random.randint(0, mask_size - 1)
            col = random.randint(0, mask_size - 1)
            
            # Calculate where this position would be in all rotations
            positions = []
            for r in range(rotation_steps):
                new_row, new_col = rotate_position(row, col, mask_size, r)
                positions.append((new_row, new_col))
            
            # If none of the rotated positions are used yet
            if not any((pos in used_positions) for pos in positions):
                used_positions.add((row, col))
                mask[row][col] = (rot == 0)  # Mark in mask if this is rotation 0
                count += 1
    
    # Apply the mask through all rotations
    text_idx = 0
    for rot in range(rotation_steps):
        for r in range(mask_size):
            for c in range(mask_size):
                # If this is a hole in the original mask
                if rot == 0 and mask[r][c]:
                    grid[r][c] = padded_text[text_idx]
                    text_idx += 1
                # If this is a hole after rotation
                elif rot > 0:
                    # Find where this position was before rotation
                    orig_r, orig_c = rotate_position(r, c, mask_size, rotation_steps - rot)
                    if mask[orig_r][orig_c]:
                        grid[r][c] = padded_text[text_idx]
                        text_idx += 1
    
    # Read off the grid row by row
    result = ''
    for row in grid:
        result += ''.join(row)
    
    return result, mask

def rotate_position(row, col, size, rotations):
    """Helper function to rotate a position in a grid."""
    for _ in range(rotations % 4):
        row, col = col, size - 1 - row
    return row, col

def process_promising_results(input_file, output_file, min_freq_score=15, max_results=100):
    """
    Process promising results from the Gromark solver, applying various transpositions.
    
    Args:
        input_file (str): Input file containing Gromark solver results
        output_file (str): Output file to write transposition results
        min_freq_score (float): Minimum frequency score to consider
        max_results (int): Maximum number of results to process
        
    Returns:
        list: Best results after transpositions
    """
    promising_results = []
    current_result = {}
    parsing_result = False
    
    # Read the input file
    with open(input_file, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            
            if line.startswith("Rank "):
                # Start of a new result
                if parsing_result and current_result:
                    promising_results.append(current_result.copy())
                current_result = {}
                parsing_result = True
            
            elif parsing_result:
                if line.startswith("Keyword:"):
                    current_result["keyword"] = line.split(":", 1)[1].strip()
                elif line.startswith("Primer:"):
                    current_result["primer"] = line.split(":", 1)[1].strip()
                elif line.startswith("Frequency Score:"):
                    current_result["freq_score"] = float(line.split(":", 1)[1].strip())
                elif line.startswith("IC:"):
                    current_result["ic"] = float(line.split(":", 1)[1].strip())
                elif line.startswith("Decrypted:"):
                    current_result["decrypted"] = line.split(":", 1)[1].strip()
                elif line.startswith("Match:"):
                    current_result["match"] = float(line.split(":", 1)[1].strip().strip('%'))
    
    # Add the last result if there is one
    if parsing_result and current_result:
        promising_results.append(current_result.copy())
    
    # Filter by frequency score and sort by match percentage and IC
    promising_results = [r for r in promising_results if r.get("freq_score", 0) >= min_freq_score]
    promising_results.sort(key=lambda x: (-x.get("match", 0), -x.get("ic", 0)))
    
    # Limit the number of results to process
    promising_results = promising_results[:max_results]
    
    print(f"Processing {len(promising_results)} promising results...")
    
    # Common keywords to try
    test_keywords = ["KRYPTOS", "BERLIN", "PALIMPSEST", "ABSCISSA", "UNDERDOG", "SANBORN"]
    
    # Apply various transpositions
    all_transposition_results = []
    
    for idx, result in enumerate(promising_results):
        print(f"Processing result {idx+1}/{len(promising_results)}")
        decrypted_text = result.get("decrypted", "")
        
        if not decrypted_text:
            continue
        
        # Try columnar transposition with various keywords
        for keyword in test_keywords:
            trans_text = columnar_transposition_encrypt(decrypted_text, keyword)
            all_transposition_results.append({
                "original_rank": idx + 1,
                "original_keyword": result.get("keyword", ""),
                "original_primer": result.get("primer", ""),
                "original_match": result.get("match", 0),
                "original_ic": result.get("ic", 0),
                "original_freq_score": result.get("freq_score", 0),
                "transposition": "Columnar",
                "trans_param": keyword,
                "text": trans_text
            })
        
        # Try double columnar transposition
        for kw1, kw2 in itertools.combinations(test_keywords, 2):
            trans_text = double_columnar_transposition(decrypted_text, kw1, kw2)
            all_transposition_results.append({
                "original_rank": idx + 1,
                "original_keyword": result.get("keyword", ""),
                "original_primer": result.get("primer", ""),
                "original_match": result.get("match", 0),
                "original_ic": result.get("ic", 0),
                "original_freq_score": result.get("freq_score", 0),
                "transposition": "Double Columnar",
                "trans_param": f"{kw1}-{kw2}",
                "text": trans_text
            })
        
        # Try Myszkowski transposition
        for keyword in test_keywords:
            trans_text = myszkowski_transposition_encrypt(decrypted_text, keyword)
            all_transposition_results.append({
                "original_rank": idx + 1,
                "original_keyword": result.get("keyword", ""),
                "original_primer": result.get("primer", ""),
                "original_match": result.get("match", 0),
                "original_ic": result.get("ic", 0),
                "original_freq_score": result.get("freq_score", 0),
                "transposition": "Myszkowski",
                "trans_param": keyword,
                "text": trans_text
            })
        
        # Try Rail Fence
        for rails in range(2, 8):
            trans_text = rail_fence_encrypt(decrypted_text, rails)
            all_transposition_results.append({
                "original_rank": idx + 1,
                "original_keyword": result.get("keyword", ""),
                "original_primer": result.get("primer", ""),
                "original_match": result.get("match", 0),
                "original_ic": result.get("ic", 0),
                "original_freq_score": result.get("freq_score", 0),
                "transposition": "Rail Fence",
                "trans_param": str(rails),
                "text": trans_text
            })
        
        # Try Route Cipher with different patterns
        for pattern in ["spiral", "snake", "diagonal"]:
            for rows in [4, 5, 7, 8, 9, 10, 11]:
                cols = math.ceil(len(decrypted_text) / rows)
                trans_text = route_cipher(decrypted_text, rows, cols, pattern)
                all_transposition_results.append({
                    "original_rank": idx + 1,
                    "original_keyword": result.get("keyword", ""),
                    "original_primer": result.get("primer", ""),
                    "original_match": result.get("match", 0),
                    "original_ic": result.get("ic", 0),
                    "original_freq_score": result.get("freq_score", 0),
                    "transposition": f"Route ({pattern})",
                    "trans_param": f"{rows}x{cols}",
                    "text": trans_text
                })
        
        # Try Scytale
        for diameter in range(3, 12):
            trans_text = scytale_encrypt(decrypted_text, diameter)
            all_transposition_results.append({
                "original_rank": idx + 1,
                "original_keyword": result.get("keyword", ""),
                "original_primer": result.get("primer", ""),
                "original_match": result.get("match", 0),
                "original_ic": result.get("ic", 0),
                "original_freq_score": result.get("freq_score", 0),
                "transposition": "Scytale",
                "trans_param": str(diameter),
                "text": trans_text
            })
    
    # Calculate IC for all transposition results
    for result in all_transposition_results:
        text = result["text"]
        # Calculate IC
        text_upper = ''.join(c for c in text.upper() if c in string.ascii_uppercase)
        if len(text_upper) <= 1:
            result["ic"] = 0
            continue
            
        char_count = {}
        for c in text_upper:
            char_count[c] = char_count.get(c, 0) + 1
            
        n = len(text_upper)
        sum_freqs = sum(count * (count - 1) for count in char_count.values())
        ic = sum_freqs / (n * (n - 1))
        
        result["ic"] = ic
    
    # Sort by IC (higher is better)
    all_transposition_results.sort(key=lambda x: -x["ic"])
    
    # Write results to output file
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(f"Total transposition results: {len(all_transposition_results)}\n")
        f.write("=" * 80 + "\n\n")
        
        for i, result in enumerate(all_transposition_results[:500]):  # Limit to top 500 results
            f.write(f"Rank {i+1}\n")
            f.write(f"Original Rank: {result['original_rank']}\n")
            f.write(f"Original Keyword: {result['original_keyword']}\n")
            f.write(f"Original Primer: {result['original_primer']}\n")
            f.write(f"Original Match: {result['original_match']:.2f}%\n")
            f.write(f"Original IC: {result['original_ic']:.6f}\n")
            f.write(f"Original Frequency Score: {result['original_freq_score']:.6f}\n")
            f.write(f"Transposition: {result['transposition']}\n")
            f.write(f"Transposition Parameter: {result['trans_param']}\n")
            f.write(f"IC: {result['ic']:.6f}\n")
            f.write(f"Text: {result['text']}\n\n")
            f.write("-" * 80 + "\n\n")
    
    print(f"Results written to {output_file}")
    return all_transposition_results[:50]  # Return top 50 results

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Apply advanced transposition techniques to promising Gromark results')
    parser.add_argument('--input', type=str, required=True, help='Input results file')
    parser.add_argument('--output', type=str, default='transposition_results.txt', help='Output file for transposition results')
    parser.add_argument('--min-freq-score', type=float, default=15.0, help='Minimum frequency score to consider')
    parser.add_argument('--max-results', type=int, default=100, help='Maximum number of results to process')
    
    args = parser.parse_args()
    
    process_promising_results(args.input, args.output, args.min_freq_score, args.max_results) 