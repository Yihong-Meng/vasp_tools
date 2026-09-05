#!/bin/bash

set -e

echo "========================================="
echo "  Step 1: Running wannier90.x"
echo "========================================="
wannier90.x wannier90.win
if [ $? -ne 0 ]; then
    echo "ERROR: wannier90.x failed. Exiting."
    exit 1
fi
echo "wannier90.x finished successfully."

echo ""
echo "========================================="
echo "  Step 2: Running python3 plot_comparar.py"
echo "========================================="
python3 plot_comparar.py
if [ $? -ne 0 ]; then
    echo "ERROR: plot_comparar.py failed. Exiting."
    exit 1
fi
echo "plot_comparar.py finished successfully."

echo ""
echo "========================================="
echo "  Extracting Final Spread from wannier90.wout"
echo "========================================="
grep -i "Final Spread (Ang^2)" wannier90.wout
if [ $? -ne 0 ]; then
    echo "Warning: Could not find 'Final Spread (Ang^2)' in wannier90.wout"
fi

echo ""
echo "========================================="
echo "  Checking for generated comparison figure"
echo "========================================="
if [ -f "band_comparison.png" ]; then
    echo "  >> band_comparison.png  (size: $(du -h band_comparison.png | cut -f1))"
    echo "  Transferring file to local via sz..."
    sz band_comparison.png
else
    echo "  WARNING: band_comparison.png not found! Please check plot_comparar.py."
fi

echo ""
echo "All tasks completed."
