import numpy as np
import math
import random
import itertools
from collections import Counter


# ----------------------------------------------Bioinformatics I----------------------------------------------
""" Interactive Text for Week 1 (Coursera I) """

def PatternCount(Text, Pattern):
    count = 0
    k = len(Pattern)
    for i in range(len(Text) - k + 1):
        if Text[i:i+k] == Pattern:
            count += 1
    return count


def FrequentWords(Text, k):
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            words.append(key)
    return words


def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern] = 0
        
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern] += 1
    return freq


def ReverseComplement(Pattern):
    rev_complement = ""
    for char in Pattern:
        if char == 'A': rev_complement += 'T'
        elif char == 'T': rev_complement += 'A'
        elif char == 'C': rev_complement += 'G'
        elif char == 'G': rev_complement += 'C'
    return rev_complement[::-1]


def PatternMatching(Pattern, Genome):
    positions = []
    k, n = len(Pattern), len(Genome)
    for i in range(n-k+1):
        if Genome[i:i+k] == Pattern:
            positions.append(i)
    return positions


def FindClump(Genome, k, L, t):
    res = set()
    freq = FrequencyMap(Genome[0:L], k)
    for key, value in freq.items():
        if value >= t: res.add(key)
    for i in range(1, len(Genome)-L+1):
        freq[Genome[i-1:i-1+k]] -= 1
        new_pattern = Genome[i+L-k:i+L]
        freq[new_pattern] = freq.get(new_pattern, 0) + 1
        if freq[new_pattern] >= t: res.add(new_pattern)
    return res


""" ---------------------------------------- """
""" Interactive Text for Week 2 (Coursera I) """

def MinimumSkew(Genome):
    positions = []
    skew = SkewArray(Genome)
    min_skew = min(skew.values())
    for key, value in skew.items():
        if value == min_skew:
            positions.append(key)
    return positions


def SkewArray(Genome):
    skew = {0:0}
    n = len(Genome)
    for i in range(1, n+1):
        skew[i] = skew[i-1]
        if Genome[i-1] == 'G': skew[i] += 1
        elif Genome[i-1] == 'C': skew[i] -= 1
    return skew


def HammingDistance(p, q):
    k = len(p)
    hamming_dist = 0
    for i in range(k):
        if p[i] != q[i]: hamming_dist += 1
    return hamming_dist


def ApproximatePatternMatching(Text, Pattern, d):
    positions = []
    k, n = len(Pattern), len(Text)
    for i in range(n-k+1):
        hamming_dist = HammingDistance(Text[i:i+k], Pattern)
        if hamming_dist <= d: positions.append(i)
    return positions


def Neighbors(Pattern, d):
    if d == 0:
        return {Pattern}
    if len(Pattern) == 0:
        return {''}
    first_char = Pattern[0]
    suffix = Pattern[1:]
    neighborhood = set()
    suffix_neighbors = Neighbors(suffix, d)
    for neighbor in suffix_neighbors:
        neighborhood.add(first_char + neighbor)
    suffix_neighbors_d_minus_1 = Neighbors(suffix, d - 1)
    for neighbor in suffix_neighbors_d_minus_1:
        for char in ['A', 'C', 'G', 'T']:
            if char != first_char:
                neighborhood.add(char + neighbor)
    return neighborhood


def FrequentWordsWithMismatches(Text, k, d):
    Patterns = set()
    freqMap = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        neighborhood = Neighbors(Pattern, d)
        rev_neighborhood = Neighbors(ReverseComplement(Pattern), d)
        for neighbor in neighborhood:
            freqMap[neighbor] = freqMap.get(neighbor, 0) + 1
        for rev_neighbor in rev_neighborhood:
            freqMap[rev_neighbor] = freqMap.get(rev_neighbor, 0) + 1
    m = max(freqMap.values())
    for Pattern, freq in freqMap.items():
        if freq == m:
            Patterns.add(Pattern)
    return Patterns


""" ---------------------------------------- """
""" Interactive Text for Week 3 (Coursera I) """

def MotifEnumeration(Dna, k, d):
    Patterns = set()
    all_kmers = set()
    for text in Dna:
        for i in range(len(text) - k + 1):
            all_kmers.add(text[i:i+k])
    for kmer in all_kmers:
        neighborhood = Neighbors(kmer, d)
        for neighbor in neighborhood:
            is_motif = True
            for text in Dna:
                found_in_text = False
                for i in range(len(text) - k + 1):
                    if HammingDistance(neighbor, text[i:i+k]) <= d:
                        found_in_text = True
                        break
                if not found_in_text:
                    is_motif = False
                    break
            if is_motif:
                Patterns.add(neighbor)
    return Patterns


def DistanceBetweenPatternAndStrings(Pattern, Dna):
    n, k = len(Dna[0]), len(Pattern)
    dist = 0
    for Text in Dna:
        hamming_dist = 1e9
        for i in range(n-k+1):
            hamming_dist = min(hamming_dist, HammingDistance(Text[i:i+k], Pattern))
        dist += hamming_dist
    return dist


# Median String
def MedianString(Dna, k):
    distance = 1e9
    ori = ''
    for _ in range(k): ori += 'A'
    Median = ''
    Patterns = Neighbors(ori, k)
    for Pattern in Patterns:
        if DistanceBetweenPatternAndStrings(Pattern, Dna) < distance:
            distance = DistanceBetweenPatternAndStrings(Pattern, Dna)
            Median = Pattern
    return Median


def Pr(Pattern, Profile):
    k = len(Profile['A'])
    prob = 1
    for j in range(k):
        prob *= Profile[Pattern[j]][j]
    return prob


def ProfileMostProbableKmer(Text, k, Profile):
    n = len(Text)
    max_prob, most_kmer = -1, ""
    for i in range(n-k+1):
        prob = Pr(Text[i:i+k], Profile)
        if prob > max_prob: max_prob, most_kmer = prob, Text[i:i+k]
    return most_kmer


def Consensus(Motifs):
    n, k = len(Motifs), len(Motifs[0])
    count = CountWithPseudocounts(Motifs)
    consensus = ""
    for j in range(k):
        freq_sym, max_freq = '', -1
        for char in "ACGT":
            if count[char][j] > max_freq: max_freq, freq_sym = count[char][j], char
        consensus += freq_sym
    return consensus


def Count(Motifs):
    m, n = len(Motifs), len(Motifs[0])
    count = {}
    for char in "ACGT":
        count[char] = []
        for j in range(n):
            count[char].append(0)
    for i in range(m):
        for j in range(n):
            count[Motifs[i][j]][j] += 1
    return count


def CountWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    count = {}
    for char in "ACGT":
        count[char] = [1 for _ in range(k)]
    for i in range(t):
        for j in range(k):
            count[Motifs[i][j]][j] += 1
    return count


def Profile(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    count = Count(Motifs)
    profile = {}
    for char in "ACGT":
        profile[char] = []
        for j in range(k):
            profile[char].append(count[char][j] / t)
    return profile


def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}
    for char in "ACGT":
        profile[char] = []
    count = CountWithPseudocounts(Motifs)
    for j in range(k):
        for char in "ACGT":
            profile[char].append(count[char][j] / (t+4))
    return profile


def Score(Motifs):
    n, k = len(Motifs), len(Motifs[0])
    count = Count(Motifs)
    consensus = Consensus(Motifs)
    score = 0
    for j in range(k):
        score += n - count[consensus[j]][j]
    return score


# Greedy Motif Search 
def GreedyMotifSearch(Dna, k, t):
    BestMotif = []
    for i in range(t):
        BestMotif.append(Dna[i][0:k])
        
    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = ProfileWithPseudocounts(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotif):
            BestMotif = Motifs
    return BestMotif


""" ---------------------------------------- """
""" Interactive Text for Week 4 (Coursera I) """

def RandomMotifs(Dna, k, t):
    random_motifs = []
    for i in range(t):
        rd = random.randint(0, len(Dna[0])-k)
        random_motifs.append(Dna[i][rd:rd+k])
    return random_motifs


def Motifs(Profile, Dna):
    motifs = []
    n, k = len(Dna), len(Profile['A'])
    for i in range(n):
        motifs.append(ProfileMostProbableKmer(Dna[i], k, Profile))
    return motifs


# Randomized Motif Search
def RandomizedMotifSearch(Dna, k, t):
    BestMotifs = RandomMotifs(Dna, k, t)
    M = RandomMotifs(Dna, k, t)
    while True:
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile, Dna)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifs


def Normalize(Probabilities):
    sum = 0
    for _, value in Probabilities.items():
        sum += value
    normalized = {}
    for key, value in Probabilities.items():
        normalized[key] = value / sum
    return normalized


def WeightedDie(Probabilities):
    kmer = '' # output variable
    p = random.uniform(0, 1)
    for kmer, prob in Probabilities.items():
        if p <= prob: return kmer
        else: p -= prob


def ProfileGeneratedString(Text, profile, k):
    n = len(Text)
    Probabilities = {}
    for i in range(n-k+1):
        Probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)
    probabilities = Normalize(Probabilities)
    return WeightedDie(probabilities)


# Gibbs Sampling
def GibbsSampler(Dna, k, t, N):
    Motifs = RandomMotifs(Dna, k, t)
    BestMotifs = Motifs.copy()
    for _ in range(N):
        i = random.randint(0, t-1)
        profile = ProfileWithPseudocounts(Motifs[:i] + Motifs[i+1:])
        new_kmer = ProfileGeneratedString(Dna[i], profile, k)
        Motifs[i] = new_kmer
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs.copy()
    return BestMotifs


""" ---------------------------------------- """
""" Additional Functions """

def FrequencyArray(Text, k):
    frequency_array = [0] * (4**k)
    for i in range(len(Text)-k+1):
        Pattern = Text[i:i+k]
        frequency_array[PatternToNumber(Pattern)] += 1
    return frequency_array


def PatternToNumber(Pattern):
    symbol_to_number = {
        'A': 0,
        'C': 1,
        'G': 2,
        'T': 3,
    }
    number = 0
    for char in Pattern:
        number = 4*number + symbol_to_number[char]
    return number


def NumberToPattern(number, k):
    pattern = ''
    number_to_symbol = {
        0: 'A',
        1: 'C',
        2: 'G',
        3: 'T',
    }
    tmp = number
    while k > 0:
        pattern += number_to_symbol[tmp % 4]
        tmp = tmp // 4
        k -= 1
    return pattern[::-1]



# ----------------------------------------------Bioinformatics II----------------------------------------------
"""Interactive Text for Week 1 (Coursera II)"""

def Composition(Text, k):
    s = []
    for i in range(len(Text)-k+1):
        s.append(Text[i:i+k])
    return s


def PathToGenome(path):
    n, k = len(path), len(path[0])
    genome = path[0]
    for i in range(1, n):
        genome += path[i][-1]
    return genome


def Overlap(Patterns):
    dict = {}
    n = len(Patterns)
    for i in range(n):
        node = Patterns[i]
        dict[node] = []
    for i in range(n):
        node = Patterns[i]
        for j in range(n):
            if j == i: continue
            adj_node = Patterns[j]
            if node[1:] == adj_node[:-1]: dict[node].append(adj_node)
    return dict


def DeBruijn(Text, k):
    dict = {}
    n = len(Text)
    for i in range(n-k+2):
        node = Text[i:i+k-1]
        dict[node] = []
    for i in range(n-k+1):
        node = Text[i:i+k]
        start_node = node[:-1]
        end_node = node[1:]
        dict[start_node].append(end_node)
    return dict


def DeBruijn(Patterns):
    dict = {}        
    for Pattern in Patterns:
        start_node = Pattern[:-1]
        end_node = Pattern[1:]
        if start_node not in dict: dict[start_node] = [end_node]
        else: dict[start_node].append(end_node)
        if end_node not in dict: dict[end_node] = []
    return dict


""" ----------------------------------------- """
""" Interactive Text for Week 2 (Coursera II) """

def EulerianCycle(Graph):
    start_node = list(Graph.keys())[0]
    cycle = []
    current_path = [start_node]
    while current_path:
        current_node = current_path[-1]
        if current_node in Graph and Graph[current_node]:
            next_node = Graph[current_node].pop()
            current_path.append(next_node)
        else:
            cycle.append(current_path.pop())
    cycle = cycle[::-1]
    return cycle


def EulerianPath(Graph):
    dict = {}
    for node in Graph: 
        dict[node] = {'out': 0, 'in': 0}
    for start_node, end_nodes in Graph.items():
        dict[start_node]['out'] += len(end_nodes)
        for end_node in end_nodes:
            dict[end_node]['in'] += 1
    start_node, end_node = None, None
    for node in dict:
        if dict[node]['out'] ==  dict[node]['in'] + 1: start_node = node
        elif dict[node]['in'] ==  dict[node]['out'] + 1: end_node = node
    Graph[end_node].append(start_node)
    cycle = []
    current_path = [start_node]
    while current_path:
        current_node = current_path[-1]
        if current_node in Graph and Graph[current_node]:
            next_node = Graph[current_node].pop()
            current_path.append(next_node)
        else:
            cycle.append(current_path.pop())
    cycle.reverse()
    path = []
    for i in range(len(cycle) - 1):
        if cycle[i] == end_node and cycle[i+1] == start_node:
            path = cycle[i+1:-1] + cycle[:i+1]
    return path


def StringRescontruction(Patterns):
    dB = DeBruijn(Patterns)
    path = EulerianPath(dB)
    Text = PathToGenome(path)
    return Text


def BinaryStrings(k):
    combinations = itertools.product("01", repeat=k)
    binary_strings = ["".join(combo) for combo in combinations]
    return binary_strings


def kUniversalCircularString(k):
    binary_kmers = BinaryStrings(k)
    graph = DeBruijn(binary_kmers)
    eulerian_cycle = EulerianCycle(graph)
    circular_str = PathToGenome(eulerian_cycle[:2**k - k + 2])
    return circular_str


def StringSpelledByGappedPatterns(GappedPatterns, k, d):
    prefixs, suffixs = [], []
    n = len(GappedPatterns)
    for GappedPattern in GappedPatterns:
        prefix, suffix = GappedPattern.split("|")
        prefixs.append(prefix)
        suffixs.append(suffix)
    prefix_str = PathToGenome(prefixs)
    suffix_str = PathToGenome(suffixs)
    for i in range(k+d, len(prefix_str)):
        if prefix_str[i] != suffix_str[i - (k+d)]:
            return
    Text = prefix_str + suffix_str[n-d-2:]
    return Text


def DeBruijnReadPairs(GappedPatterns):
    prefixs, suffixs = [], []
    n = len(GappedPatterns)
    for GappedPattern in GappedPatterns:
        prefix, suffix = GappedPattern.split("|")
        prefixs.append(prefix)
        suffixs.append(suffix)
    dict = {}
    for i in range(n):
        start_node = prefixs[i][:-1] + '|' + suffixs[i][:-1]
        end_node = prefixs[i][1:] + '|' + suffixs[i][1:]
        dict[start_node] = []
        dict[end_node] = []
    for i in range(n):
        start_node = prefixs[i][:-1] + '|' + suffixs[i][:-1]
        end_node = prefixs[i][1:] + '|' + suffixs[i][1:]
        dict[start_node].append(end_node)
    return dict


def StringRescontructionReadPairs(GappedPatterns, k, d):
    dB = DeBruijnReadPairs(GappedPatterns)
    path = EulerianPath(dB)
    Text = StringSpelledByGappedPatterns(path, k, d)
    return Text


def MaximalNonBranchingPaths(Graph):
    dict = {}
    for node in Graph: 
        dict[node] = {'out': 0, 'in': 0}
    for start_node, end_nodes in Graph.items():
        dict[start_node]['out'] += len(end_nodes)
        for end_node in end_nodes:
            dict[end_node]['in'] += 1
    paths = []
    visited = {}
    for node in Graph.keys(): visited[node] = False
    for start_node in Graph:
        if (dict[start_node]['in'] != 1 or dict[start_node]['out'] != 1) and dict[start_node]['out'] >= 1:
            visited[start_node] = True
            for current_node in Graph[start_node]:
                visited[current_node] = True
                path = [start_node, current_node]
                while dict[current_node]['in'] == dict[current_node]['out'] == 1:
                    current_node = Graph[current_node][0]
                    path.append(current_node)
                    visited[current_node] = True
                paths.append(path)
    for node in Graph:
        if visited[node] == False:
            visited[node] = True
            path = [node]
            current_node = Graph[node][0]
            while current_node != node:
                visited[current_node] = True
                path.append(current_node)
                current_node = Graph[current_node][0]
            path.append(current_node)
            paths.append(path)
    return paths


def ContigGeneration(Patterns):
    Graph = DeBruijn(Patterns)
    paths = MaximalNonBranchingPaths(Graph)
    contigs = []
    for path in paths:
        contigs.append(PathToGenome(path))
    return contigs


""" ----------------------------------------- """
""" Interactive Text for Week 3 (Coursera II) """

with open('RNA_codon_table_1.txt', 'r') as file:
    RNA_codon_table = {}
    lines = file.readlines()
    for line in lines:
        line = line.strip()
        codon = line[:3]
        if len(line) == 5: RNA_codon_table[codon] = line[4]
        else: RNA_codon_table[codon] = '!'


amino_acid_mass = {
    'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101,
    'C': 103, 'I': 113, 'L': 113, 'N': 114, 'D': 115, 'K': 128,
    'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156,
    'Y': 163, 'W': 186
}


def ProteinTranslation(RNA, RNA_codon_table):
    protein = ''
    for i in range(0, len(RNA), 3):
        codon = RNA[i:i+3]
        if RNA_codon_table[codon] != '!': protein += RNA_codon_table[codon]
        else: return protein
    return protein


def DNAtoRNA(DNA):
    RNA = ''
    for char in DNA:
        if char == 'T': RNA += 'U'
        else: RNA += char
    return RNA


def PeptideEncoding(Text, Peptide, RNA_codon_table):
    substrs = []
    protein = [''] * 3
    rev_protein = [''] * 3
    k = len(Peptide)
    for i in range(0, 3):
        substr = Text[i:i+3*k]
        rev_substr = ReverseComplement(substr)
        for j in range(0, len(substr), 3):
            codon = substr[j:j+3]
            protein[i] += RNA_codon_table[DNAtoRNA(codon)]
            codon = rev_substr[j:j+3]
            rev_protein[i] += RNA_codon_table[DNAtoRNA(codon)]
        if Peptide in [protein[i], rev_protein[i]]:
            substrs.append(substr)
    for i in range(3, len(Text)-3*k+1):
        substr = Text[i:i+3*k]
        protein[i%3] = protein[i%3][1:] + RNA_codon_table[DNAtoRNA(substr[-3:])]
        rev_protein[i%3] = RNA_codon_table[DNAtoRNA(ReverseComplement(substr[-3:]))] + rev_protein[i%3][:-1]
        if Peptide in [protein[i%3], rev_protein[i%3]]:
            substrs.append(substr)
    return substrs


def LinearSpectrum(Peptide):
    n = len(Peptide)
    prefix_mass = [0] * (n+1)
    for i in range(1, n+1):
        prefix_mass[i] = prefix_mass[i-1] + Peptide[i-1]
    linear_spectrum = []
    linear_spectrum.append(0)
    for i in range(0, n):
        for j in range(i+1, n+1):
            mass = prefix_mass[j] - prefix_mass[i]
            linear_spectrum.append(mass)
    linear_spectrum.sort()
    return linear_spectrum


def CycloSpectrum(Peptide):
    n = len(Peptide)
    prefix_mass = [0] * (n+1)
    for i in range(1, n+1):
        prefix_mass[i] = prefix_mass[i-1] + Peptide[i-1]
    cyclo_spectrum = []
    cyclo_spectrum.append(0)
    for i in range(0, n):
        for j in range(i+1, n+1):
            mass = prefix_mass[j] - prefix_mass[i]
            cyclo_spectrum.append(mass)
            if i > 0 and j < n:
                cyclo_spectrum.append(prefix_mass[n] - mass)
    cyclo_spectrum.sort()
    return cyclo_spectrum


def CountPeptidesWithMass(target_mass, amino_acid_masses):
    dp = [0] * (target_mass + 1)
    dp[0] = 1
    for i in range(1, target_mass + 1):
        for mass in amino_acid_masses:
            if i >= mass: dp[i] += dp[i - mass]
    return dp[target_mass]


def Mass(peptide):
    return sum(peptide)


def ParentMass(spectrum):
    return max(spectrum) if spectrum else 0


def Expand(peptides, amino_acid_masses):
    expanded = []
    for peptide in peptides:
        for mass in amino_acid_masses:
            new_peptide = peptide + [mass]
            expanded.append(new_peptide)
    return expanded


def IsConsistent(peptide, spectrum):
    linear_spec = LinearSpectrum(peptide)
    spec_counts = {}
    for mass in spectrum:
        spec_counts[mass] = spec_counts.get(mass, 0) + 1
    for mass in linear_spec:
        if spec_counts.get(mass, 0) == 0:
            return False
        spec_counts[mass] -= 1
    return True


def CyclopeptideSequencing(spectrum):
    candidates_peptides = [[]]
    final_peptides = []
    while candidates_peptides:
        candidates_peptides = Expand(candidates_peptides)
        surviving_candidates = []
        for peptide in candidates_peptides:
            if Mass(peptide) == ParentMass(spectrum):
                if CycloSpectrum(peptide) == spectrum and peptide not in final_peptides:
                    final_peptides.append(peptide)
            elif IsConsistent(peptide, spectrum):
                surviving_candidates.append(peptide)
        candidates_peptides = surviving_candidates
    return final_peptides


""" ----------------------------------------- """
""" Interactive Text for Week 4 (Coursera II) """

def CycloScore(peptide, spectrum):
    peptide_spectrum = CycloSpectrum(peptide)
    spec_dict = {}
    for mass in spectrum:
        spec_dict[mass] = spec_dict.get(mass, 0) + 1
    score = 0
    for mass in peptide_spectrum:
        if mass in spec_dict and spec_dict[mass] > 0:
            score += 1
            spec_dict[mass] -= 1
    return score


def LinearScore(peptide, spectrum):
    peptide_spectrum = LinearSpectrum(peptide)
    spec_dict = {}
    for mass in spectrum:
        spec_dict[mass] = spec_dict.get(mass, 0) + 1
    score = 0
    for mass in peptide_spectrum:
        if mass in spec_dict and spec_dict[mass] > 0:
            score += 1
            spec_dict[mass] -= 1
    return score


def Trim(leaderboard, spectrum, N):
    if len(leaderboard) <= N:
        return leaderboard
    scores = [LinearScore(peptide, spectrum) for peptide in leaderboard]
    sorted_pairs = sorted(zip(leaderboard, scores), key=lambda item: item[1], reverse=True)
    trimmed_leaderboard = []
    cutoff_score = sorted_pairs[N-1][1]
    for peptide, score in sorted_pairs:
        if score >= cutoff_score:
            trimmed_leaderboard.append(peptide)
        else: break
    return trimmed_leaderboard


def LeaderboardCyclopeptideSequencing(spectrum, N):
    leaderboard = [[]]
    leader_peptides = []
    best_score = 0
    while leaderboard:
        leaderboard = Expand(leaderboard)
        surviving_peptides = []
        for peptide in leaderboard:
            if Mass(peptide) == ParentMass(spectrum):
                current_score = CycloScore(peptide, spectrum)
                if current_score > best_score:
                    best_score = current_score
                    leader_peptides = [peptide]
                elif current_score == best_score:
                    leader_peptides.append(peptide)
            if Mass(peptide) <= ParentMass(spectrum):
                surviving_peptides.append(peptide)
        leaderboard = Trim(surviving_peptides, spectrum, N)
    return leader_peptides


def SpectralConvolution(spectrum):
    n = len(spectrum)
    convolution = []
    for i in range (1, n):
        for j in range(0, i):
            if spectrum[j] != spectrum[i]:
                convolution.append(abs(spectrum[i] - spectrum[j]))
    return convolution


def TrimConvolution(spectrum, M):
    convolution = SpectralConvolution(spectrum)
    filtered_convo = list(filter(lambda x: 57 <= x <= 200, convolution))
    mass_counts = Counter(filtered_convo)
    sorted_masses = sorted(mass_counts.items(), key=lambda item: item[1], reverse=True)
    trimmed_leaderboard = []
    cutoff_count = sorted_masses[M-1][1]
    for mass, count in sorted_masses:
        if count >= cutoff_count:
            trimmed_leaderboard.append(mass)
        else: break
    return trimmed_leaderboard


def ConvolutionCyclopeptideSequencing(spectrum, M, N):
    leaderboard = [[]]
    leader_peptide = ''
    # leader_peptides = []
    # best_score = 0
    amino_acid_masses = TrimConvolution(spectrum, M)
    while leaderboard:
        leaderboard = Expand(leaderboard, amino_acid_masses)
        surviving_peptides = []
        for peptide in leaderboard:
            if Mass(peptide) == ParentMass(spectrum):
                current_score = CycloScore(peptide, spectrum)
                if current_score > CycloScore(leader_peptide, spectrum):
                    leader_peptide = peptide
                # if current_score > best_score:
                #     best_score = current_score
                #     leader_peptides = [peptide]
                # elif current_score == best_score:
                #     leader_peptides.append(peptide)
            if Mass(peptide) <= ParentMass(spectrum):
                surviving_peptides.append(peptide)
        leaderboard = Trim(surviving_peptides, spectrum, N)
    return leader_peptide


def recur_Turnpike(points, dist_counts):
    if not any(dist_counts.values()):
        return points
    y = 0
    for key, value in dist_counts.items():
        if value > 0 and key > y:
            y = key
    # 1. New point = y
    x = y
    temp_counts = dist_counts.copy()
    is_consistent = True
    for point in points:
        dist = abs(x - point)
        if temp_counts.get(dist, 0) == 0: 
            is_consistent = False
            break
        else:
            temp_counts[dist] -= 1
    if is_consistent:
        new_points = points + [x]
        solution = recur_Turnpike(new_points, temp_counts)
        if solution:
            return solution
    # 2. New point = max(points) - y
    x = max(points) - y
    temp_counts = dist_counts.copy()
    is_consistent = True
    for point in points:
        dist = abs(x - point)
        if temp_counts.get(dist, 0) == 0: 
            is_consistent = False
            break
        else:
            temp_counts[dist] -= 1
    if is_consistent:
        new_points = points + [x]
        solution = recur_Turnpike(new_points, temp_counts)
        if solution:
            return solution
    # No solution for both path
    return None


def Turnpike(pairwise_dist):
    if not pairwise_dist: return []
    pairwise_dist.sort()
    sz = int(np.sqrt(len(pairwise_dist)))
    pairwise_dist = pairwise_dist[(sz**2 + sz)//2:]
    dist_counts = Counter(pairwise_dist)
    max_width = max(pairwise_dist)
    points = [0, max_width]
    dist_counts[max_width] -= 1
    solution = recur_Turnpike(points, dist_counts)
    return sorted(solution) if solution else None



# ----------------------------------------------Bioinformatics III---------------------------------------------
""" Interactive Text for Week 1 (Coursera III) """

def DPChange(money, coins):
    dp = [0] * (money + 1)
    for i in range(1, money+1):
        dp[i] = 1e9
        for coin in coins:
            if i >= coin: dp[i] = min(dp[i-coin] + 1, dp[i])
    return dp[money]


def ManhattanTourist(n, m, down, right):
    s = [[0 for _ in range(m+1)] for _ in range(n+1)]
    for j in range(1, m+1): s[0][j] = s[0][j-1] + right[0][j-1]
    for i in range(1, n+1): s[i][0] = s[i-1][0] + down[i-1][0]
    for i in range(1, n+1):
        for j in range(1, m+1):
            s[i][j] = max(s[i-1][j] + down[i-1][j], s[i][j-1] + right[i][j-1])
    return s[n][m]


def LongestCommonSubsequence(v, w):
    n, m = len(v), len(w)
    s = [[0 for _ in range(m+1)] for _ in range(n+1)]
    backtrack = [['' for _ in range(m+1)] for _ in range(n+1)]
    for i in range(1, n+1):
        for j in range(1, m+1):
            match = 0
            if v[i-1] == w[j-1]: match = 1
            s[i][j] = max(s[i-1][j], s[i][j-1], s[i-1][j-1] + match)
            if s[i][j] == s[i-1][j]: backtrack[i][j] = 'S'
            elif s[i][j] == s[i][j-1]: backtrack[i][j] = 'E'
            elif s[i][j] == s[i-1][j-1] + match: backtrack[i][j] = 'SE'
    lcs = ''
    i, j = n, m
    while i > 0 and j > 0:
        if backtrack[i][j] == 'S': i -= 1
        elif backtrack[i][j] == 'E': j -= 1
        elif backtrack[i][j] == 'SE':
            if v[i-1] == w[j-1]: lcs += v[i-1]
            i -= 1
            j -= 1
    lcs = lcs[::-1]
    return lcs


def DAG_LongestPath(DAG_graph, start_node, end_node):
    if start_node == end_node: return 0, [start_node]
    max_length = -1e9
    longest_path = [start_node]
    for node, weight in DAG_graph[start_node]:
        length, path = DAG_LongestPath(DAG_graph, node, end_node)
        if max_length < length  + weight:
            max_length = length + weight
            longest_path = [start_node] + path
    return max_length, longest_path


""" ------------------------------------------ """
""" Interactive Text for Week 2 (Coursera III) """

def ParseScoringMatrix(filepath):
    with open(filepath, 'r') as file:
        matrix = {}
        lines = file.readlines()
        proteins = lines[0].strip().split()
        for i in range(1, len(lines)):
            scores = lines[i].strip().split()
            protein = scores[0]
            matrix[protein] = {}
            for j in range(1, len(scores)):
                matrix[protein][proteins[j-1]] = int(scores[j])
    return matrix


def GlobalAlignment(score_matrix_path, indel_pen, nucl_1, nucl_2):
    score_matrix = ParseScoringMatrix(score_matrix_path)
    n, m = len(nucl_1), len(nucl_2)
    s = [[0 for _ in range(m+1)] for _ in range(n+1)]
    backtrack = [['' for _ in range(m+1)] for _ in range(n+1)]
    for j in range(1, m+1):
        s[0][j] = s[0][j-1] - indel_pen
        backtrack[0][j] = 'E'
    for i in range(1, n+1):
        s[i][0] = s[i-1][0] - indel_pen
        backtrack[i][0] = 'S'
    for i in range(1, n+1):
        for j in range(1, m+1):
            match = score_matrix[nucl_1[i-1]][nucl_2[j-1]]
            s[i][j] = max(s[i-1][j] - indel_pen, s[i][j-1] - indel_pen, s[i-1][j-1] + match)
            if s[i][j] == s[i-1][j] - indel_pen: backtrack[i][j] = 'S'
            elif s[i][j] == s[i][j-1] - indel_pen: backtrack[i][j] = 'E'
            elif s[i][j] == s[i-1][j-1] + match: backtrack[i][j] = 'SE'
    align_nucl_1, align_nucl_2 = backtrack_alignment(n, m, backtrack, nucl_1, nucl_2)
    return s[n][m], align_nucl_1, align_nucl_2


def LocalAlignment(score_matrix_path, indel_pen, nucl_1, nucl_2):
    score_matrix = ParseScoringMatrix(score_matrix_path)
    n, m = len(nucl_1), len(nucl_2)
    s = [[0 for _ in range(m+1)] for _ in range(n+1)]
    backtrack = [['' for _ in range(m+1)] for _ in range(n+1)]
    for i in range(1, n+1):
        for j in range(1, m+1):
            match = score_matrix[nucl_1[i-1]][nucl_2[j-1]]
            s[i][j] = max(0, s[i-1][j] + indel_pen, s[i][j-1] + indel_pen, s[i-1][j-1] + match)
            if s[i][j] == s[i-1][j] + indel_pen: backtrack[i][j] = 'S'
            elif s[i][j] == s[i][j-1] + indel_pen: backtrack[i][j] = 'E'
            elif s[i][j] == s[i-1][j-1] + match: backtrack[i][j] = 'SE'
    arr = np.array(s)
    row_idx, col_idx = np.unravel_index(np.argmax(arr), arr.shape)
    max_align_score = arr[row_idx, col_idx]
    align_nucl_1, align_nucl_2 = backtrack_alignment(row_idx, col_idx, backtrack, nucl_1, nucl_2)
    return max_align_score, align_nucl_1, align_nucl_2


def EditDistanceAlignment(nucl_1, nucl_2):
    n, m = len(nucl_1), len(nucl_2)
    dp = [[0] * (m+1) for _ in range(n+1)]
    backtrack = [[''] * (m+1) for _ in range(n+1)]
    for i in range(n+1):
        dp[i][0] = i
        backtrack[i][0] = 'S'
    for j in range(m+1):
        dp[0][j] = j
        backtrack[0][j] = 'E'
    for i in range(1, n+1):
        for j in range(1, m+1):
            cost = 0 if nucl_1[i-1] == nucl_2[j-1] else 1
            dp[i][j] = min(dp[i-1][j] + 1, dp[i][j-1] + 1, dp[i-1][j-1] + cost)
            if dp[i][j] == dp[i-1][j] + 1: backtrack[i][j] = 'S'
            elif dp[i][j] == dp[i][j-1] + 1: backtrack[i][j] = 'E'
            elif dp[i][j] == dp[i-1][j-1] + cost: backtrack[i][j] = 'SE'       
    align_nucl_1, align_nucl_2 = backtrack_alignment(n, m, backtrack, nucl_1, nucl_2)
    return dp[n][m], align_nucl_1, align_nucl_2


def FittingAlignment(score_matrix_path, indel_pen, nucl_1, nucl_2):
    score_matrix = ParseScoringMatrix(score_matrix_path)
    n, m = len(nucl_1), len(nucl_2)
    s = [[0 for _ in range(m+1)] for _ in range(n+1)]
    backtrack = [['' for _ in range(m+1)] for _ in range(n+1)]
    for j in range(1, m+1):
        s[0][j] = s[0][j-1] - indel_pen
        backtrack[0][j] = 'E'
    for i in range(1, n+1):
        for j in range(1, m+1):
            # match = score_matrix[nucl_1[i-1]][nucl_2[j-1]]
            match = 1 if nucl_1[i-1] == nucl_2[j-1] else -1
            s[i][j] = max(s[i-1][j] - indel_pen, s[i][j-1] - indel_pen, s[i-1][j-1] + match)
            if s[i][j] == s[i-1][j] - indel_pen: backtrack[i][j] = 'S'
            elif s[i][j] == s[i][j-1] - indel_pen: backtrack[i][j] = 'E'
            elif s[i][j] == s[i-1][j-1] + match: backtrack[i][j] = 'SE'
    arr = np.array(s)
    row_idx = np.argmax(arr[:, m])
    max_align_score = arr[row_idx][m]
    align_nucl_1, align_nucl_2 = backtrack_alignment(row_idx, m, backtrack, nucl_1, nucl_2)
    return max_align_score, align_nucl_1, align_nucl_2


def OverlapAlignment(score_matrix_path, indel_pen, nucl_1, nucl_2):
    score_matrix = ParseScoringMatrix(score_matrix_path)
    n, m = len(nucl_1), len(nucl_2)
    s = [[0 for _ in range(m+1)] for _ in range(n+1)]
    backtrack = [['' for _ in range(m+1)] for _ in range(n+1)]
    for j in range(1, m+1):
        s[0][j] = s[0][j-1] - indel_pen
        backtrack[0][j] = 'E'
    for i in range(1, n+1):
        for j in range(1, m+1):
            # match = score_matrix[nucl_1[i-1]][nucl_2[j-1]]
            match = 1 if nucl_1[i-1] == nucl_2[j-1] else -2
            s[i][j] = max(s[i-1][j] - indel_pen, s[i][j-1] - indel_pen, s[i-1][j-1] + match)
            if s[i][j] == s[i-1][j] - indel_pen: backtrack[i][j] = 'S'
            elif s[i][j] == s[i][j-1] - indel_pen: backtrack[i][j] = 'E'
            elif s[i][j] == s[i-1][j-1] + match: backtrack[i][j] = 'SE'
    arr = np.array(s)
    col_idx = np.argmax(arr[n])
    max_align_score = arr[n][col_idx]
    align_nucl_1, align_nucl_2 = backtrack_alignment(n, col_idx, backtrack, nucl_1, nucl_2)
    return max_align_score, align_nucl_1, align_nucl_2


def backtrack_alignment(start_row, start_col, backtrack, nucl_1, nucl_2):
    align_nucl_1, align_nucl_2 = '', ''
    i, j = start_row, start_col
    while i > 0 or j > 0:
        if backtrack[i][j] == 'S':
            align_nucl_1 += nucl_1[i-1]
            align_nucl_2 += '-'
            i -= 1
        elif backtrack[i][j] == 'E':
            align_nucl_1 += '-' 
            align_nucl_2 += nucl_2[j-1]
            j -= 1
        elif backtrack[i][j] == 'SE':
            align_nucl_1 += nucl_1[i-1]
            align_nucl_2 += nucl_2[j-1]
            i -= 1
            j -= 1
        else:
            break
    return align_nucl_1[::-1], align_nucl_2[::-1]


""" ------------------------------------------ """
""" Interactive Text for Week 3 (Coursera III) """

def AffineGapPenaltiesAlignment(score_matrix_path, gap_opening_pen, gap_extension_pen, nucl_1, nucl_2):
    score_matrix = ParseScoringMatrix(score_matrix_path)
    n, m = len(nucl_1), len(nucl_2)
    lower = [[0] * (m + 1) for _ in range(n+1)]
    middle = [[0] * (m + 1) for _ in range(n+1)]
    upper = [[0] * (m + 1) for _ in range(n+1)]
    backtrack_lower = [[''] * (m + 1) for _ in range(n+1)]
    backtrack_middle = [[''] * (m + 1) for _ in range(n+1)]
    backtrack_upper = [[''] * (m + 1) for _ in range(n+1)]
    lower[0] = [-1e9 for _ in range(m+1)]
    for i in range(1, n+1):
        middle[i][0] = -float('inf')
        upper[i][0] = -float('inf')
        lower[i][0] = gap_opening_pen + (i-1) * gap_extension_pen
    for j in range(1, m+1):
        middle[0][j] = -float('inf')
        lower[0][j] = -float('inf')
        upper[0][j] = gap_opening_pen + (j-1) * gap_extension_pen
    for i in range(1, n+1):
        for j in range(1, m+1):
            score = score_matrix[nucl_1[i-1]][nucl_2[j-1]]
            # score = 1 if nucl_1[i-1] == nucl_2[j-1] else -5
            # Lower matrix
            lower_scores = [lower[i-1][j] + gap_extension_pen, middle[i-1][j] + gap_opening_pen]
            lower[i][j] = max(lower_scores)
            backtrack_lower[i][j] = 'lower' if lower_scores[0] >= lower_scores[1] else 'middle'
            # Upper matrix
            upper_scores = [upper[i][j-1] + gap_extension_pen, middle[i][j-1] + gap_opening_pen]
            upper[i][j] = max(upper_scores)
            backtrack_upper[i][j] = 'upper' if upper_scores[0] >= upper_scores[1] else 'middle'
            # Middle matrix
            middle_scores = [lower[i][j], upper[i][j], middle[i-1][j-1] + score]
            middle[i][j] = max(middle_scores)
            if middle[i][j] == middle_scores[0]: backtrack_middle[i][j] = 'lower'
            elif middle[i][j] == middle_scores[1]: backtrack_middle[i][j] = 'upper'
            else: backtrack_middle[i][j] = 'middle'
    final_scores = [lower[n][m], middle[n][m], upper[n][m]]
    max_align_score = max(final_scores)
    start_matrix = ''
    if max_align_score == final_scores[0]: start_matrix = 'lower'
    elif max_align_score == final_scores[1]: start_matrix = 'middle'
    else: start_matrix = 'upper'
    align_nucl_1, align_nucl_2 = backtrack_AffineAlignment(n, m, backtrack_lower, backtrack_middle, backtrack_upper, nucl_1, nucl_2, start_matrix)
    return max_align_score, align_nucl_1, align_nucl_2


def backtrack_AffineAlignment(start_row, start_col, bt_lower, bt_middle, bt_upper, nucl_1, nucl_2, current_matrix):
    align_nucl_1, align_nucl_2 = '', ''
    i, j = start_row, start_col
    while i > 0 or j > 0:
        if current_matrix == 'middle':
            prev_matrix = bt_middle[i][j]
            if prev_matrix == 'middle':
                align_nucl_1 += nucl_1[i-1]
                align_nucl_2 += nucl_2[j-1]
                i -= 1
                j -= 1
            current_matrix = prev_matrix
        elif current_matrix == 'lower':
            align_nucl_1 += nucl_1[i-1]
            align_nucl_2 += '-'
            current_matrix = bt_lower[i][j]
            i -= 1
        elif current_matrix == 'upper':
            align_nucl_1 += '-' 
            align_nucl_2 += nucl_2[j-1]
            current_matrix = bt_upper[i][j]
            j -= 1
    return align_nucl_1[::-1], align_nucl_2[::-1]


def check_score(score_matrix_path, gap_opening_pen, gap_extension_pen, nucl_1, nucl_2):
    align_score = 0
    score_matrix = ParseScoringMatrix(score_matrix_path)
    for i in range(len(nucl_1)):
        if nucl_1[i] == '-':
            if i >= 1 and nucl_1[i-1] == '-': align_score += gap_extension_pen
            else: align_score += gap_opening_pen
        elif nucl_2[i] == '-':
            if i >= 1 and nucl_2[i-1] == '-': align_score += gap_extension_pen
            else: align_score += gap_opening_pen
        else:
            score = score_matrix[nucl_1[i]][nucl_2[i]]
            # score = 1 if nucl_1[i] == nucl_2[i] else -5
            align_score += score
    return align_score


def MiddleEdgeLinearSpace(score_matrix_path, indel_pen, nucl_1, nucl_2):
    score_matrix = ParseScoringMatrix(score_matrix_path)
    n, m = len(nucl_1), len(nucl_2)
    mid = m // 2

    # Forward source to mid
    col = [i * indel_pen for i in range(n + 1)]
    for j in range(1, mid + 1):
        new_col = [col[0] + indel_pen]
        for i in range(1, n + 1):
            # match_score = 1 if nucl_2[i-1] == nucl_1[j-1] else -1
            match_score = score_matrix[nucl_1[i-1]][nucl_2[j-1]]
            score_diag = col[i-1] + match_score
            score_down = new_col[i-1] + indel_pen
            score_right = col[i] + indel_pen
            new_col.append(max(score_diag, score_down, score_right))
        col = new_col
    from_source = col

    # Backward sink to mid
    rev_col = [(n - i) * indel_pen for i in range(n + 1)]
    to_sink_mid_plus_1, to_sink_mid = rev_col, rev_col
    for j in range(m - 1, mid - 1, -1):
        new_rev = [0] * (n + 1)
        new_rev[n] = rev_col[n] + indel_pen
        for i in range(n - 1, -1, -1):
            # match_score = 1 if nucl_2[i] == nucl_1[j] else -1
            match_score = score_matrix[nucl_1[i]][nucl_2[j]]
            score_diag = rev_col[i+1] + match_score
            score_down = new_rev[i+1] + indel_pen
            score_right = rev_col[i] + indel_pen
            new_rev[i] = max(score_diag, score_down, score_right)
        
        rev_col = new_rev
        if j == mid + 1: to_sink_mid_plus_1 = rev_col
        if j == mid: to_sink_mid = rev_col

    # Find middle node
    max_node_score = -float('inf')
    middle_node = -1
    for i in range(n + 1):
        score = from_source[i] + to_sink_mid[i]
        if score > max_node_score:
            max_node_score = score
            middle_node = i
    print(f"{middle_node} {max_node_score}")

    # Find middle edge
    i = middle_node
    best_edge = ()
    score_right = from_source[i] + indel_pen + to_sink_mid_plus_1[i] # Right edge
    if score_right == max_node_score:
        best_edge = ((i, mid), (i, mid + 1))
    
    if not best_edge and i < n:
        score_down = from_source[i] + indel_pen + to_sink_mid[i+1] # Down edge
        if score_down == max_node_score:
            best_edge = ((i, mid), (i + 1, mid))
    
    if not best_edge and i < n:
        # match_score = 1 if nucl_2[i] == nucl_1[j] else -1
        match_score = score_matrix[nucl_1[i]][nucl_2[j]]
        score_diag = from_source[i] + match_score + to_sink_mid_plus_1[i + 1] # Diagonal edge
        if score_diag == max_node_score:
            best_edge = ((i, mid), (i + 1, mid + 1))
    
    print(f"{score_right} {score_down} {score_diag}")
    return best_edge


def get_scores_forward(score_matrix_path, indel_pen, nucl_1, nucl_2, top, bottom, left, right):
    score_matrix = ParseScoringMatrix(score_matrix_path)
    n_sub, m_sub = bottom - top, right - left

    # Forward source to mid
    col = [i * indel_pen for i in range(n_sub + 1)]
    for j in range(1, m_sub + 1):
        new_col = [col[0] + indel_pen]
        for i in range(1, n_sub + 1):
            # match_score = 1 if nucl_1[i - 1 + top] == nucl_2[j - 1 + left] else -1
            match_score = score_matrix[nucl_1[i - 1 + top]][nucl_2[j - 1 + left]]
            score_diag = col[i-1] + match_score
            score_down = new_col[i-1] + indel_pen
            score_right = col[i] + indel_pen
            new_col.append(max(score_diag, score_down, score_right))
        col = new_col
    return col


def get_scores_backward(score_matrix_path, indel_pen, nucl_1, nucl_2, top, bottom, left, right):
    score_matrix = ParseScoringMatrix(score_matrix_path)
    n_sub, m_sub = bottom - top, right - left
    rev_col = [(n_sub - i) * indel_pen for i in range(n_sub + 1)]
    col, col_second_last = rev_col, rev_col
    for j in range(m_sub - 1, - 1, -1):
        new_rev = [0] * (n_sub + 1)
        new_rev[n_sub] = rev_col[n_sub] + indel_pen
        for i in range(n_sub - 1, -1, -1):
            # match_score = 1 if nucl_1[i + top] == nucl_2[j + left] else -1
            match_score = score_matrix[nucl_1[i + top]][nucl_2[j + left]]
            score_diag = rev_col[i+1] + match_score
            score_down = new_rev[i+1] + indel_pen
            score_right = rev_col[i] + indel_pen
            new_rev[i] = max(score_diag, score_down, score_right)
        
        rev_col = new_rev
        if j == 1: col_second_last = rev_col
        if j == 0: col = rev_col
    return col, col_second_last


def MiddleEdge(score_matrix_path, indel_pen, nucl_1, nucl_2, top, bottom, left, right):
    score_matrix = ParseScoringMatrix(score_matrix_path)
    n_sub, m_sub = bottom - top, right - left
    mid = (left + right) // 2

    from_source = get_scores_forward(score_matrix_path, indel_pen, nucl_1, nucl_2, top, bottom, left, mid)
    to_sink_mid, to_sink_mid_plus_1 = get_scores_backward(score_matrix_path, indel_pen, nucl_1, nucl_2, top, bottom, mid, right)

    # Find middle node
    max_node_score = -float('inf')
    middle_node = -1
    for i in range(n_sub + 1):
        score = from_source[i] + to_sink_mid[i]
        if score > max_node_score:
            max_node_score = score
            middle_node = i

    # Find middle edge
    i = middle_node
    best_edge = ()
    score_right = from_source[i] + indel_pen + to_sink_mid_plus_1[i] # Right edge
    if score_right == max_node_score:
        best_edge = ((i, mid - left), (i, mid + 1 - left), 'R')
    
    if not best_edge and i < n_sub:
        score_down = from_source[i] + indel_pen + to_sink_mid[i+1] # Down edge
        if score_down == max_node_score:
            best_edge = ((i, mid - left), (i + 1, mid - left), 'Dwn')
    
    if not best_edge and i < n_sub:
        # match_score = 1 if nucl_1[i + top] == nucl_2[mid] else -1
        match_score = score_matrix[nucl_1[i + top]][nucl_2[mid]]
        score_diag = from_source[i] + match_score + to_sink_mid_plus_1[i + 1] # Diagonal edge
        if score_diag == max_node_score:
            best_edge = ((i, mid - left), (i + 1, mid + 1 - left), 'D')

    return best_edge, max_node_score


def LinearSpaceAlignment(score_matrix_path, indel_pen, nucl_1, nucl_2, top, bottom, left, right):
    # Base cases
    if left == right: return ['Dwn'] * (bottom - top)
    if top == bottom: return ['R'] * (right - left)

    edge_tuple, _ = MiddleEdge(score_matrix_path, indel_pen, nucl_1, nucl_2, top, bottom, left, right)
    mid_node_start, mid_node_end, edge_type = edge_tuple
    mid_i_start = top + mid_node_start[0]
    mid_j_start = left + mid_node_start[1]
    mid_i_end = top + mid_node_end[0]
    mid_j_end = left + mid_node_end[1]

    path_left = LinearSpaceAlignment(score_matrix_path, indel_pen, nucl_1, nucl_2, top, mid_i_start, left, mid_j_start)
    path_right = LinearSpaceAlignment(score_matrix_path, indel_pen, nucl_1, nucl_2, mid_i_end, bottom, mid_j_end, right)

    return path_left + [edge_type] + path_right


def Multiple_LCS(score_matrix_path, indel_pen, nucl_1, nucl_2, nucl_3):
    score_matrix = ParseScoringMatrix(score_matrix_path)
    m, n, p = len(nucl_1), len(nucl_2), len(nucl_3)
    s = [[[0 for _ in range(p+1)] for _ in range(n+1)] for _ in range(m+1)]
    backtrack = [[['' for _ in range(p+1)] for _ in range(n+1)] for _ in range(m+1)]
    for i in range(1, m+1):
        s[i][0][0] = s[i-1][0][0] + indel_pen
        backtrack[i][0][0] = 'S'
    for j in range(1, n+1):
        s[0][j][0] = s[0][j-1][0] + indel_pen
        backtrack[0][j][0] = 'E'
    for k in range(1, p+1):
        s[0][0][k] = s[0][0][k-1] + indel_pen
        backtrack[0][0][k] = 'F'
    for i in range(1, m+1):
        for j in range(1, n+1):
            s[i][j][0] = max(s[i-1][j][0] + indel_pen, s[i][j-1][0] + indel_pen, s[i-1][j-1][0] + indel_pen)
            if s[i][j][0] == s[i-1][j][0] + indel_pen: backtrack[i][j][0] = 'S'
            elif s[i][j][0] == s[i][j-1][0] + indel_pen: backtrack[i][j][0] = 'E'
            else: backtrack[i][j][0] = 'SE'
    for j in range(1, n+1):
        for k in range(1, p+1):
            s[0][j][k] = max(s[0][j-1][k] + indel_pen, s[0][j][k-1] + indel_pen, s[0][j-1][k-1] + indel_pen)
            if s[0][j][k] == s[0][j-1][k] + indel_pen: backtrack[0][j][k] = 'E'
            elif s[0][j][k] == s[0][j][k-1] + indel_pen: backtrack[0][j][k] = 'F'
            else: backtrack[0][j][k] = 'FE'
    for i in range(1, m+1):
        for k in range(1, p+1):
            s[i][0][k] = max(s[i-1][0][k] + indel_pen, s[i][0][k-1] + indel_pen, s[i-1][0][k-1] + indel_pen)
            if s[i][0][k] == s[i-1][0][k] + indel_pen: backtrack[i][0][k] = 'S'
            elif s[i][0][k] == s[i][0][k-1] + indel_pen: backtrack[i][0][k] = 'F'
            else: backtrack[i][0][k] = 'FS'

    for i in range(1, m+1):
        for j in range(1, n+1):
            for k in range(1, p+1):
                # match = score_matrix[nucl_1[i-1]][nucl_2[j-1]]
                match = 1 if (nucl_1[i-1] == nucl_2[j-1] == nucl_3[k-1]) else 0
                S = s[i-1][j][k] + indel_pen
                E = s[i][j-1][k] + indel_pen
                F = s[i][j][k-1] + indel_pen
                FE = s[i][j-1][k-1] + indel_pen
                FS = s[i-1][j][k-1] + indel_pen
                SE = s[i-1][j-1][k] + indel_pen
                FSE = s[i-1][j-1][k-1] + match

                s[i][j][k] = max(S, E, F, FE, FS, SE, FSE)
                if s[i][j][k] == FSE: backtrack[i][j][k] = 'FSE'
                elif s[i][j][k] == FE: backtrack[i][j][k] = 'FE'
                elif s[i][j][k] == FS: backtrack[i][j][k] = 'FS'
                elif s[i][j][k] == SE: backtrack[i][j][k] = 'SE'
                elif s[i][j][k] == S: backtrack[i][j][k] = 'S'
                elif s[i][j][k] == F: backtrack[i][j][k] = 'F'
                elif s[i][j][k] == E: backtrack[i][j][k] = 'E'
    align_nucl_1, align_nucl_2, align_nucl_3 = backtrack_multi_alignment(m, n, p, backtrack, nucl_1, nucl_2, nucl_3)
    return s[m][n][p], align_nucl_1, align_nucl_2, align_nucl_3


def backtrack_multi_alignment(start_row, start_col, start_hgt, backtrack, nucl_1, nucl_2, nucl_3):
    align_nucl_1, align_nucl_2, align_nucl_3 = '', '', ''
    i, j, k = start_row, start_col, start_hgt
    while i > 0 or j > 0 or k > 0:
        if backtrack[i][j][k] == 'S':
            align_nucl_1 += nucl_1[i-1]
            align_nucl_2 += '-'
            align_nucl_3 += '-'
            i -= 1
        elif backtrack[i][j][k] == 'E':
            align_nucl_1 += '-'
            align_nucl_2 += nucl_2[j-1]
            align_nucl_3 += '-'
            j -= 1
        elif backtrack[i][j][k] == 'F':
            align_nucl_1 += '-'
            align_nucl_2 += '-'
            align_nucl_3 += nucl_3[k-1]
            k -= 1
        elif backtrack[i][j][k] == 'FE':
            align_nucl_1 += '-'
            align_nucl_2 += nucl_2[j-1]
            align_nucl_3 += nucl_3[k-1]
            j -= 1
            k -= 1
        elif backtrack[i][j][k] == 'FS':
            align_nucl_1 += nucl_1[i-1]
            align_nucl_2 += '-'
            align_nucl_3 += nucl_3[k-1]
            i -= 1
            k -= 1
        elif backtrack[i][j][k] == 'SE':
            align_nucl_1 += nucl_1[i-1]
            align_nucl_2 += nucl_2[j-1]
            align_nucl_3 += '-'
            i -= 1
            j -= 1
        elif backtrack[i][j][k] == 'FSE':
            align_nucl_1 += nucl_1[i-1]
            align_nucl_2 += nucl_2[j-1]
            align_nucl_3 += nucl_3[k-1]
            i -= 1
            j -= 1
            k -= 1
        else:
            break
    return align_nucl_1[::-1], align_nucl_2[::-1], align_nucl_3[::-1]


def parse_adjacency_list(filepath):
    graph = {}
    with open(filepath, 'r') as file:
        lines = file.readlines()
        for line in lines:
            source_str, dest_str = line.strip().split(' -> ')
            src, dests = int(source_str), [int(dest) for dest in dest_str.split(',')]
            if src not in graph: graph[src] = dests
            else: graph[src] += dests
            for dest in dests:
                if dest not in graph: graph[dest] = []
    return graph


def TopologicalOrdering(graph):
    topological_list = []
    in_degree = {}
    for node, adj_lst in graph.items():
        for dst in adj_lst:
            if dst not in in_degree: in_degree[dst] = 1
            else: in_degree[dst] += 1
        if node not in in_degree: in_degree[node] = 0
    candidates = [node for node, value in in_degree.items() if value == 0]
    while candidates:
        a = candidates.pop()
        topological_list.append(a)
        while graph[a]:
            b = graph[a].pop()
            in_degree[b] -= 1
            if in_degree[b] == 0: candidates.append(b)
    for node, adj_lst in graph.items():
        if adj_lst: return "The input graph is not a DAG."
    return topological_list


""" ------------------------------------------ """
""" Interactive Text for Week 4 (Coursera III) """

def GreedySorting(permutation):
    permu_seq = []
    n = len(permutation)
    current_permu = permutation.copy()
    for i in range(1, n+1):
        if current_permu[i-1] == i: continue
        elif current_permu[i-1] == -i:
            current_permu[i-1] = i
            permu_seq.append(current_permu.copy())
        else:
            if i in current_permu:
                index1 = current_permu.index(i)
                current_permu[i-1: index1+1] = current_permu[i-1: index1+1][::-1]
                for j in range (i-1, index1+1): current_permu[j] = -current_permu[j]
                permu_seq.append(current_permu.copy())
                current_permu[i-1] = i
                permu_seq.append(current_permu.copy())
            else:
                index2 = current_permu.index(-i)
                current_permu[i-1: index2+1] = current_permu[i-1: index2+1][::-1]
                for j in range (i-1, index2+1): current_permu[j] = -current_permu[j]
                permu_seq.append(current_permu.copy())
    return permu_seq


def CountBreakpoints(permutation):
    n = len(permutation)
    permutation = [0] + permutation + [n+1]
    bp = 0
    for i in range(1, n+2):
        if permutation[i] - permutation[i-1] != 1: bp += 1
    return bp


""" ------------------------------------------ """
""" Interactive Text for Week 5 (Coursera III) """

def ChromosomeToCycle(Chromosome):
    n = len(Chromosome)
    Nodes = [0] * (2*n)
    for j in range(0, n):
        i = Chromosome[j]
        if i > 0:
            Nodes[2*j] = 2*i-1
            Nodes[2*j+1] = 2*i
        else:
            Nodes[2*j] = -2*i
            Nodes[2*j+1] = -2*i-1
    return Nodes


def CycleToChromosome(Nodes):
    n = len(Nodes)
    Chromosome = [0] * (n//2)
    for j in range(0, n//2):
        if Nodes[2*j] < Nodes[2*j+1]: Chromosome[j] = Nodes[2*j+1] // 2
        else: Chromosome[j] = -Nodes[2*j] // 2 
    return Chromosome


def ColoredEdges(P):
    edges = []
    for chromosome in P:
        nodes = ChromosomeToCycle(chromosome)
        n = len(chromosome)
        for j in range(0, n-1):
            edges.append((nodes[2*j+1], nodes[2*j+2]))
        edges.append((nodes[2*n-1], nodes[0]))
    return edges


def GraphToGenome(genome_graph):
    P = []
    nodes = []
    for edge in genome_graph:
        nodes.append(edge[0])
        if edge[0] < edge[1]: nodes.append(edge[1])
        else:
            nodes = [edge[1]] + nodes
            chromosome = CycleToChromosome(nodes)
            P.append(chromosome)
            nodes = []
    return P


def CountCycles(graph):
    n = (len(graph) - 1) // 2
    visited = [False for _ in range(2*n+1)]
    cycles = 0
    for node in range(1, 2*n+1):
        if visited[node] or not graph[node]: continue
        while not visited[node]:
            visited[node] = True
            for dst in graph[node]:
                if not visited[dst]:
                    node = dst
                    break
        cycles += 1
    return cycles


def TwoBreakDistance(P, Q):
    edges_P = ColoredEdges(P)
    edges_Q = ColoredEdges(Q)
    n = len(edges_P)
    graph = [[] for _ in range(2*n+1)]
    for src, dst in edges_P:
        graph[src].append(dst)
        graph[dst].append(src)
    for src, dst in edges_Q:
        graph[src].append(dst)
        graph[dst].append(src)
    cycles_num = CountCycles(graph)
    two_break_dist = n - cycles_num
    return two_break_dist


def TwoBreakOnGenomeGraph(edges, i1, i2, i3, i4):
    n = len(edges)
    graph = [[] for _ in range(2*n+1)]
    for src, dst in edges:
        graph[src].append(dst)
    graph[i1].remove(i2) if i2 in graph[i1] else graph[i2].remove(i1)
    graph[i1].append(i3)
    graph[i3].remove(i4) if i4 in graph[i3] else graph[i4].remove(i3)
    graph[i2].append(i4)
    two_break_edges = [(i, j) for i, neighbors in enumerate(graph) for j in neighbors]
    return two_break_edges


def TwoBreakOnGenome(genome, i1, i2, i3, i4):
    edges = ColoredEdges(genome)
    print(edges)
    n = len(edges)
    graph = [[] for _ in range(2*n+1)]
    for src, dst in edges:
        graph[src].append(dst)
    if i1 in graph[i2]: i1, i2 = i2, i1
    if i3 in graph[i4]: i3, i4 = i4, i3
    graph[i1].remove(i2)
    graph[i1].append(i4)
    graph[i3].remove(i4)
    graph[i3].append(i2)
    two_break_edges = [(i, j) for i, neighbors in enumerate(graph) for j in neighbors]
    print(two_break_edges)
    two_break_genome = GraphToGenome(two_break_edges)
    return two_break_genome