import numpy as np 
from memory_profiler import memory_usage
import time

# Constants 
UP = (-1,0)
LEFT = (0, -1)
TOPLEFT = (-1, -1)
ORIGIN = (0, 0)
PATHS = [UP, LEFT, TOPLEFT]

# From HW1
def read_blosum62(path):
    """
    Reads in the ncbi's BLOSUM62.txt file and loads the scoring matrix
    into a dictionary.
    
    :param: path is the full path in the local filesystem at which the .txt file is located
    :return: a dictionary of dictionaries which will hold the cost of various amino acid
    substitutions as defined in BLOSUM62.
    """
    delta = {}
    with open(path, 'r') as f:
        lines = f.readlines()[6:]
        keys = lines[0].split()
        keys[-1] = '-'
        for i, line in enumerate(lines[1:]):
            delta[keys[i]] = {k : int(v) for (k,v) in zip(keys, line.split()[1:])}  
    return delta


# From HW1
def traceback_global(v, w, pointers):
    i,j = len(v), len(w)
    new_v = []
    new_w = []
    while True:
        #print(str(i) + ", " + str(j))

        di, dj = pointers[i][j]
        if (di,dj) == LEFT:
            new_v.append('-')
            new_w.append(w[j-1])
        elif (di,dj) == UP:
            new_v.append(v[i-1])
            new_w.append('-')
        elif (di,dj) == TOPLEFT:
            new_v.append(v[i-1])
            new_w.append(w[j-1])
        i, j = i + di, j + dj
        if (i <= 0 and j <= 0):
            break
    return ''.join(new_v[::-1])+'\n'+''.join(new_w[::-1])

# From HW1
def global_align(v, w, delta):
    """
    Returns the score of the maximum scoring alignment of the strings v and w, as well as the actual alignment as 
    computed by traceback_global. 
    
    :param: v
    :param: w
    :param: delta
    """
    M = [[0 for j in range(len(w)+1)] for i in range(len(v)+1)]
    pointers = [[ORIGIN for j in range(len(w)+1)] for i in range(len(v)+1)]
    score, alignment = None, None
    # YOUR CODE HERE

    for i in range(len(v)+1):
      for j in range(len(w)+1):

        # Set the initial values for "only gaps"
        if i == 0 and j == 0:
          M[i][j] = 0
          pointers[i][j] = ORIGIN 
        elif i == 0:
          M[i][j] = M[i][j - 1] + delta['-'][w[j - 1]]
          pointers[i][j] = LEFT 
        elif j == 0:
          M[i][j] = M[i - 1][j] + delta[v[i - 1]]['-']
          pointers[i][j] = UP 
        else:
          char_v = v[i - 1]
          char_w = w[j - 1]

          # Calculate the scores coming from up, left, or topleft 
          up = M[i - 1][j] + delta[char_v]['-']
          left = M[i][j - 1] + delta['-'][char_w]
          topleft = M[i - 1][j - 1] + delta[char_v][char_w]

          scores = [up, left, topleft]
          pos = np.argmax(scores)

          M[i][j] = scores[pos]
          pointers[i][j] = PATHS[pos]

    score = M[len(v)][len(w)]
    
    alignment = traceback_global(v,w, pointers)
    return score, alignment

# Linear space global alignment function
def hirschberg(v, w, delta):
    """
    Returns the score of the maximum scoring alignment of the strings v and w, as well as the actual alignment.
    
    :param: v
    :param: w
    :param: delta
    """
    score = 0
    alignment = None

    n = len(v)
    m = len(w)

    # Base case
    if m < 2 or n < 2:
        return global_align(v, w, delta)

    # Find the optimal i that must pass through the midpoint
    midpoint = m // 2
    weights = weight(v, w, delta, midpoint)
    i_star = np.argmax(weights)

    # Now recurse on the prefix and suffix subproblems 
    prefix_score, prefix_alignment = hirschberg(v[:i_star], w[:midpoint], delta)
    suffix_score, suffix_alignment = hirschberg(v[i_star:], w[midpoint:], delta)

    # Combine the prefix and suffix alignments
    prefix_split = prefix_alignment.split('\n')
    suffix_split = suffix_alignment.split('\n')
    alignment = prefix_split[0] + suffix_split[0] + '\n' + prefix_split[1] + suffix_split[1]

    return prefix_score + suffix_score, alignment 

# Weight is the sum of the prefix and the suffix
def weight(v, w, delta, split):
    
    n = len(v)
    prefix_w = w[:split+1]
    suffix_w = v[split:]

    prefix = get_prefix(v, prefix_w, delta)
    suffix = get_suffix(v, suffix_w, delta)

    return [prefix[i] + suffix[n - i] for i in range(n + 1)]

# Find the score of the prefix, keeping two columns in memory at a time
def get_prefix(v, w, delta):

    n = len(v)
    m = len(w)

    # Go two columns at a time 
    curr = [0 for i in range(n + 1)]
    prev = [0 for i in range(n + 1)]

    # For each column j calculate all of the M[i][j] values by using the previous column
    for j in range(m + 1):
        for i in range(n + 1):

            if i == 0 and j == 0: # Base Case
                curr[i] = 0
            elif i == 0: # First row case, just add gap panelty to previous column's start 
                curr[i] = prev[0] + delta['-'][w[j - 1]]
            elif j == 0: # First column just adds up gap penalties
                curr[i] = curr[i-1] + delta[v[i - 1]]['-']
            else: # Normal case, check all three positions before
                char_v = v[i - 1]
                char_w = w[j - 1]

                # Calculate the scores coming from up, left, or topleft 
                up = curr[i - 1] + delta[char_v]['-']
                left = prev[i] + delta['-'][char_w]
                topleft = prev[i-1] + delta[char_v][char_w]

                scores = [up, left, topleft]
                curr[i] = np.max(scores)

        # Now current is previous and we clear current    
        prev = curr 
        curr = [0 for i in range(n + 1)]

    # Returm the whole last column
    return prev

# Find the score of the suffix, keeping two columns in memory at a time
def get_suffix(v, w, delta):
    n = len(v)
    m = len(w)

    # Go two columns at a time 
    curr = [0 for i in range(n + 1)]
    prev = [0 for i in range(n + 1)]

    # For each column j calculate all of the M[i][j] values by using the previous column
    for j in range(m + 1):
        for i in range(n + 1):

            if i == 0 and j == 0: # Base Case
                curr[i] = 0
            elif i == 0: # First row case, just add gap panelty to previous column's start 
                curr[i] = prev[0] + delta['-'][w[j - 1]]
            elif j == 0: # First column just adds up gap penalties
                curr[i] = curr[i-1] + delta[v[i - 1]]['-']
            else: # Normal case, check all three positions before
                char_v = v[n - i]
                char_w = w[m - j]

                # Calculate the scores coming from up, left, or topleft 
                up = curr[i - 1] + delta[char_v]['-']
                left = prev[i] + delta['-'][char_w]
                topleft = prev[i-1] + delta[char_v][char_w]

                scores = [up, left, topleft]
                curr[i] = np.max(scores)

        # Now current is previous and we clear current    
        prev = curr 
        curr = [0 for i in range(n + 1)]

    # Return the whole last column
    return prev


def main():
    blosum = read_blosum62('./BLOSUM62.txt')

    # Local Sanity test from HW1 
    scoreA = None
    scoreB = None
    alignmentA = None
    alignmentB = None

    # Human Insulin from https://www.rcsb.org/fasta/entry/4F1C/display
    humanA = 'GIVEQCCTSICSLYQLENYCN'
    humanB = 'FVNQHLCGSHLVEALYLVCGERGFFYTPKT'

    # Bovine Insulin from https://www.rcsb.org/fasta/entry/2ZP6/display
    bovineA = 'GIVEQCCASVCSLYQLENYCN'
    bovineB = 'FVNQHLCGSHLVEALYLVCGERGFFYTPKA'

    scoreA, alignmentA = global_align(humanA, bovineA, blosum)
    scoreB, alignmentB = global_align(humanB, bovineB, blosum)

    print("Needleman-Wunsch Score: ", scoreA)
    print(alignmentA)
    print()

    h_scoreA, h_alignmentA = hirschberg(humanA, bovineA, blosum)
    h_scoreB, h_alignmentB = hirschberg(humanB, bovineB, blosum)
    print("Hirshberg Score: ", h_scoreA)
    print(h_alignmentA)
    print('--------------------------------------------------------------------')
    
    print("Needleman-Wunsch Score: ", scoreB)
    print(alignmentB)
    print()

    print("Hirshberg Score: ", h_scoreB)
    print(h_alignmentB)

    # Run both Needleman-Wunsch and Hirschberg for the three datasets 
    for i in range(1, 4):
        infile = "cmp" + str(i) + ".txt"
        outfile = "out" + str(i) + ".txt"

        f = open(infile, "r")
        seqA = f.readline().rstrip()
        seqB = f.readline().rstrip()
        f.close()

        start = time.time()
        nw_score, nw_alignment = global_align(seqA, seqB, blosum)
        end = time.time()
        nw_time = end - start 
        nw_mem = memory_usage((global_align, (seqA, seqB, blosum)))

        start = time.time()
        h_score, h_alignment = hirschberg(seqA, seqB, blosum)
        end = time.time()
        h_time = end - start 
        h_mem = memory_usage((hirschberg, (seqA, seqB, blosum)))

        try:
            f = open(outfile, "x")
            f.write("Needleman-Wunsch Score: " + str(nw_score) + '\n')
            f.write("Needleman-Wunsch Time: " + str(nw_time) + '\n')
            f.write("Needleman-Wunsch Memory: " + str(max(nw_mem)) + '\n')
            f.write(nw_alignment)
            f.write('\n')
            f.write('\n')
            f.write("Hirschberg Score: " + str(nw_score) + '\n')
            f.write("Hirschberg Time: " + str(h_time) + '\n')
            f.write("Hirschberg Memory: " + str(max(h_mem)) + '\n')
            f.write(h_alignment)
            f.close()
        except:
            pass


if __name__ == '__main__':
    main()