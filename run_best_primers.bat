@echo off
echo Kryptos K4 Gromark Solver - Testing Bean's Most Promising Primers
echo =================================================================

REM First, test the base 10, length 5 primers
echo.
echo Testing base 10, length 5 primers...
python gromark_solver.py --base 10 --length 5 --bean-primers

REM Then test equivalent primers with only 9 digits (26717 and 84393)
echo.
echo Testing the two equivalent primers with only 9 digits...
python gromark_solver.py --specific-primers 26717,84393

REM Test the primer with best IC value
echo.
echo Testing primer with best IC value (98800)...
python gromark_solver.py --specific-primers 98800

REM Test base 8, length 5 primers with period 84
echo.
echo Testing base 8, length 5 primers...
python gromark_solver.py --base 8 --length 5 --bean-primers

REM Test base 10, length 4 primers
echo.
echo Testing base 10, length 4 primers...
python gromark_solver.py --base 10 --length 4 --bean-primers

REM Try some alternative suggestions from the Berlin Clock reference
echo.
echo Testing base 5 primers (Berlin Clock reference)...
python gromark_solver.py --base 5 --length 5 --bean-primers

echo.
echo Testing base 12 primers (Berlin Clock reference)...
python gromark_solver.py --base 12 --length 5 --bean-primers

REM Try the most promising primers with alternative expansion
echo.
echo Testing most promising primers with alternative key expansion...
python gromark_solver.py --specific-primers 26717,84393,98800 --alt-expansion

echo.
echo Testing complete!
pause 