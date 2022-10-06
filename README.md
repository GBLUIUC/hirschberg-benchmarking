## hirschberg-benchmarking
Implementation and Benchmarking of Hirschberg Algorithm against Needleman-Wunsch

# Dependencies
This code requires

Python 3 

Python memory_profiler package

# Usage
You can simply run the comparison between Needleman-Wunsch and Hirschberg with command

`python3 hirschberg.py`

This will run both methods on the three biological datasets in ./data/cmp*.txt. The format of these input data files is simply a line for each sequence. After running the algorithms, it will report the scores, alignments, running times, and peak memory usages in the corresponding ./results/out*.txt file. Please see the project report for analysis and discussion of these results.

As as sanity check, this program will also print to the terminal the results of two alignment problems provided in HW1.

