# Kryptos K4 Gromark Solver

A Python implementation of a brute-force Gromark cipher solver targeting the Kryptos K4 segment, using known plaintext cribs.

## About Kryptos K4

Kryptos is a sculpture located at CIA headquarters in Langley, Virginia, created by artist Jim Sanborn. It contains four encrypted messages (K1-K4), with the fourth section (K4) remaining unsolved since its installation in 1990.

Two known plaintext "cribs" have been identified in K4:
- "EASTNORTHEAST" at positions 21-34
- "BERLINCLOCK" at positions 49-59

This tool uses these cribs to brute-force a Gromark cipher solution.

## The Gromark Cipher

The Gromark cipher is a numerical key-based encryption system that:
1. Starts with a primer (a small set of digits)
2. Expands the primer into a full key by adding consecutive pairs of digits (modulo the base)
3. Uses the expanded key to perform simple additive shift operations on each letter

## Bean's Cryptodiagnosis

This implementation incorporates findings from Richard Bean's paper "Cryptodiagnosis of Kryptos K4," which performed a detailed statistical analysis of K4's properties. The paper identified:

- 39 specific primers for base 10, length 5 that satisfy all constraints from known plaintexts
- 3 specific primers for base 10, length 4
- 4 specific primers for base 8, length 5
- Observations about the minor differences between plaintext and ciphertext letters
- Statistical properties of vertical bigrams when the ciphertext is arranged at width 21

The solver automatically tests these promising primers first.

## Usage

```
python gromark_solver.py [options]
```

### Options

- `--base INT`: Number base to use for the key (default: 10)
- `--length INT`: Primer length to test (default: 5)
- `--processes INT`: Number of CPU processes to use (default: all available)
- `--quadgrams PATH`: Path to quadgram statistics file (optional)
- `--top INT`: Number of top candidates to display (default: 10)
- `--specific-primers LIST`: Comma-separated list of specific primers to try
- `--bean-primers`: Use the primers identified in Bean's research paper
- `--alt-expansion`: Use alternative key expansion method (sum all primer digits)

### Examples

Basic usage with default settings (base 10, primer length 5):
```
python gromark_solver.py
```

Try a different base and primer length:
```
python gromark_solver.py --base 8 --length 4
```

Use specific primers identified in Bean's paper:
```
python gromark_solver.py --bean-primers
```

Try a specific primer:
```
python gromark_solver.py --specific-primers 26717
```

## Performance

The script uses multiprocessing to leverage all available CPU cores for testing primers in parallel. For reference:
- Base 10, length 5 = 100,000 primers (completes in seconds on a modern CPU)
- Base 8, length 6 = ~260,000 primers
- Base 12, length 5 = ~250,000 primers

## Variations and Extensions

Use `gromark_variations.py` to test keyword-mixed alphabets:
```
python gromark_variations.py --keywords KRYPTOS,BERLIN,CLOCK
```

For GPU acceleration with large primer spaces, use `gromark_gpu.py`:
```
python gromark_gpu.py --base 10 --length 5
```

## Extending the Solver

To test variations of the Gromark cipher:
1. Use `--alt-expansion` to try the alternative key expansion method
2. Try different alphabets by passing custom keywords to `gromark_variations.py`
3. Explore alternative bases like 5 or 12 as suggested by the "Berlin Clock" clue

## License

This code is provided for educational and research purposes. 