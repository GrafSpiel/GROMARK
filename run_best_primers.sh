#!/bin/bash

echo "Kryptos K4 Gromark Solver - Testing Bean's Most Promising Primers"
echo "================================================================="

# First, test the base 10, length 5 primers
echo
echo "Testing base 10, length 5 primers..."
python3 gromark_solver.py --base 10 --length 5 --bean-primers

# Then test equivalent primers with only 9 digits (26717 and 84393)
echo
echo "Testing the two equivalent primers with only 9 digits..."
python3 gromark_solver.py --specific-primers 26717,84393

# Test the primer with best IC value
echo
echo "Testing primer with best IC value (98800)..."
python3 gromark_solver.py --specific-primers 98800

# Test base 8, length 5 primers with period 84
echo
echo "Testing base 8, length 5 primers..."
python3 gromark_solver.py --base 8 --length 5 --bean-primers

# Test base 10, length 4 primers
echo
echo "Testing base 10, length 4 primers..."
python3 gromark_solver.py --base 10 --length 4 --bean-primers

# Try some alternative suggestions from the Berlin Clock reference
echo
echo "Testing base 5 primers (Berlin Clock reference)..."
python3 gromark_solver.py --base 5 --length 5 --bean-primers

echo
echo "Testing base 12 primers (Berlin Clock reference)..."
python3 gromark_solver.py --base 12 --length 5 --bean-primers

# Try the most promising primers with alternative expansion
echo
echo "Testing most promising primers with alternative key expansion..."
python3 gromark_solver.py --specific-primers 26717,84393,98800 --alt-expansion

echo
echo "Testing complete!" 