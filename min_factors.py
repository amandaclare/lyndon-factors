import Lyndon
import re
import itertools
import os

class Mapping:
    def __init__(self, n):
        self.alphabet_loc = ord('a') + n
        self.mapping = {}

    def assign(self,letter):
        if letter not in self.mapping:
            self.mapping[letter] = self.alphabet_loc
            self.alphabet_loc += 1

    def unassign_up_to(self, letter):
        n = self.look_up(letter)
        # remove items bigger than n
        self.mapping = { k:v for item in self.mapping.items if v <= n }
        self.alphabet_loc = n + 1

    def clear(self):
        self.mapping.clear()
        self.reset(0)
        
    def reset(self, n):
        '''reset the alphabet_loc to position n in alphabet (n is 0-indexed)'''
        self.alphabet_loc = ord('a') + n
            
    def look_up(self,letter):
        return self.mapping.get(letter, -1)

    def map_string(self, string):
        out = []
        for c in string:
            if c not in self.mapping:
                self.assign(c)
            out.append(chr(self.mapping[c]))
        return ''.join(out)

    def assign_all(self,string):
        for letter in string:
            self.assign(letter)
    
    def __str__(self):
        out = []
        out.append("The remapping of characters:\n")
        for (k,v) in self.mapping.items():
            out.append(k + "-->" + chr(v)+"\n")
        return(''.join(out))



    
def read_fasta(filename):
    seq = ""
    with open(filename) as f:
        header = next(f)
        for line in f:
            if line.startswith('>'):
                break
            else:
                seq += line.strip().upper().translate({ord(c):'' for c in 'KWNRSMY'})
    return seq



def exp_parikh_vector(s):
    """get a list of exponent parikh vectors, in order of first occurrence of chars.
    return a list of pairs: unique char, corresponding exponent vector"""
    uniques = []
    counts = {}
    last_seen = None
    new_count = True
    for c in s:
        if c != last_seen:
            new_count = True
            last_seen = c
        else:
            new_count = False
        if c not in counts:
            uniques.append(c)
            counts[c] = []
        if new_count:
            counts[c].append(1)
        else:
            counts[c][-1] += 1
    eps = [(c,''.join(map(str, counts[c]))) for c in uniques]
    return eps



def ep_lyndon_factors(u_ep):
    """Find Lyndon factors of the exponent parikh vector components"""
    for (u,p) in u_ep:
        factors = list(Lyndon.ChenFoxLyndon(p, co_lex=True))
        yield (u,factors)


# -------------------------------------------------------
# printing

def print_eps(eps):
    print("Exponent parikhs:")
    for u,p in eps:
        print(u,p)
    print()


def print_fL(fL_prs):
    print("Lyndon factors of eps:")
    for (letter, facs) in fL_prs:
        print(letter, " >= ".join(facs))
    print()

# ---------------------------------------------------------    

def sort_pis(fL_prs):
    return sorted(((len(facs), char, facs) for (char, facs) in fL_prs), key=lambda fl: fl[0])


def choose_pi(fL_prs):
    """Choose the leftmost pr with the minimal number of factors"""
    # assumes fL_prs is not empty
    (orig, facs) = fL_prs[0]
    min_num = len(facs)
    if len(fL_prs) == 1:
        return (min_num, orig, facs)
    else: # check the others
        for (o, fs) in fL_prs[1:]:
            num = len(fs)
            if num < min_num:
                min_num = num
                orig = o
                facs = fs
        return (min_num, orig, facs)



def get_Xs_after_char(s, char):
    """Get all pieces of string that exist between the runs of the char provided
    (ignoring any prefix to the first occurrence of char)"""
    return re.split(char+'+', s)[1:]



        
def assign_to_Xs(m, facs, Xs):
    """
    We need to assign ordering to the letters in order:
    by considering each exponent group (largest to smallest) and by leftmost first
    Params:
    - a mapping object m for storing the new ordering
    - the list of Lyndon factors for the exponents for this letter
    - the list of string components occurring after each run of the letter in the original string
    """
    for fac in facs:
        lf = len(fac)
        Xs_block = Xs[0:lf]
        Xs = Xs[lf:]

        # we map each exponent to their set of X blocks.
        exp_X_dict = {}
        for (e,X) in zip(map(int, fac), Xs_block):
            if e in exp_X_dict:
                exp_X_dict[e].append(X)
            else:
                exp_X_dict[e] = [X]

        # in reverse numeric order of exponents (greatest first)
        sorted_exponents = iter(sorted(exp_X_dict.items(), reverse=True))
        (e, X_vals) = next(sorted_exponents)
        # for (e, X_vals) in
        for later_X_val in X_vals[1:]:
            good = True
            if len(later_X_val) > len(X_vals[0]):
                # problem because no longer Lyndon
                good = False
                break
            else:
                for (x,y) in zip(X_vals[0], later_X_val):
                    if x == y:
                        m.assign(x)
                    else: # different
                        m.assign(x)
                        m.assign(y)
                        if m.look_up(x) > m.look_up(y):
                            good = False  # inconsistent
                        break
            if not good:
                # undo m assignments
                # break # from the later X_vals with same exponent
                return False

    return True

    # now assign any leftover letters in the remaining X components arbitrarily
    #all_letters = ''.join(X_vals)
    #m.assign_all(all_letters)

        

def before_char(s, c):
    # this could be optimised - we can work out from parikhs
    n = s.find(c)
    return s[:n]


def alg1(s, m, no_backtrack = False):
    # make exponent parikh vectors
    eps = exp_parikh_vector(s)
    #print_eps(eps)

    # factorise each vector (reverse order)
    fL_prs = list(ep_lyndon_factors(eps))
    #print_fL(fL_prs)

    good = False
    pis = sort_pis(fL_prs)

    for (num, unique_char, facs) in pis: #choose_pi(fL_prs)
        string_prefix = before_char(s, unique_char)
        #print('prefix',string_prefix)
        m.reset(len(set(string_prefix))) # allow lower nums for prefix
        m.assign(unique_char)
        Xs = get_Xs_after_char(s, unique_char)
        good = assign_to_Xs(m, facs, Xs)
        #print(m)
        if good or no_backtrack:
            break
        else:
            m.clear()

    if good or no_backtrack:
        # sort prefix
        if string_prefix:
            n = m.alphabet_loc - ord('a')
            m.reset(0)
            alg1(string_prefix, m) # not caring about result
            m.reset(n)
        m.assign_all(s)
        reordered = m.map_string(s)
        return reordered
    else:
        
        # failed to make consistent - just use first one
        (num, unique_char, facs) = pis[0]
        m.clear()
        string_prefix = before_char(s, unique_char)
        m.reset(len(set(string_prefix))) # allow lower nums for prefix
        m.assign(unique_char)
        Xs = get_Xs_after_char(s, unique_char)
        assign_to_Xs(m, facs, Xs)
        if string_prefix:
            n = m.alphabet_loc - ord('a')
            m.reset(0)
            alg1(string_prefix, m) # not caring about result
            m.reset(n)
        m.assign_all(s)
        reordered = m.map_string(s)
        return reordered
        
# --------------------------------------------------------------------

def all_reorders(perms):
    for cs in perms:
        m = Mapping(0)
        for c in cs:
            m.assign(c)
        yield m

        

def do_all_reorderings(s):
    """Make all possible reorderings so that we can compare our chosen reordering"""
    all_lens = []
    smallest = len(s)
    smallest_m = None
    #print('\nAll possible reorderings')
    unique_letters = set(s)
    perms = itertools.permutations(unique_letters)
    ms = all_reorders(perms)
    for mapping in ms:
        r = mapping.map_string(s)
        fs = list(Lyndon.ChenFoxLyndon(r))
        #print(len(fs))
        #print(" >= ".join(fs))
        if len(fs) < smallest:
            smallest = len(fs)
            smallest_m = mapping
        #print()
        all_lens.append(len(fs))
    #print(smallest)
    #print(smallest_m)
    r = smallest_m.map_string(s)
    fs = list(Lyndon.ChenFoxLyndon(r))
    #print(len(fs))
    #print(" >= ".join(fs))
    #print()
    print(sorted(all_lens))
    #print()
    


# --------------------------------------------------------------------

if __name__ == "__main__":
    
    #s = "bbbffbbcf"
    #s = "ATGGTACTGACGATTTATCCTGACGAACTCGTACAAATAGTGTCTGATAAAATTGCTTCAAATAAGGGAAGTATGTTCA"
    #s = "ATGGTACTGACGATTTATCCTGACGAACTCGTACAAATAGTGTCTGATAAAATTGCTTCAAATAAGGGAAGTATGTTCATGTCTCATTCTCCTTTTCGGCTCCGTTTAGGTGATAAACGTACTATATTGTGAAAGATTATTTACTAACGACACATTGAAGAAATCACTTTGAATCAGCTGTGGGATATATCTGGTAAATATTTTGATTTGTCTGATAAAAAAGTTAAACAGTTCGTGCTTTCATGCGTGATATTGAAAAAGGACATTGAGGTGTATTGTGATGGTGCTATAACAACTAAAAATGTGACTGATATTATAGGCGACGCTAATCATTCATACTCGGTTGGGATTACTGAGGACAGCCTATGGACATTATTAACGGGATACACAAAAAAGGAGTCAACTATTGGAAATTCTGCATTTGAACTACTTCTCGAAGTTGCCAAATCAGGAGAAAAAGGGATCAATACTATGGATTTGGCGCAGGTAACTGGGCAAGATCCTAGAAGTGTGACTGGACGTATCAAGAAAATAAACCACCTGTTAACAAGTTCACAACTGATTTATAAGGGACACGTCGTGAAGCAATTGAAGCTAAAAAAATTCAGCCATGACGGGGTGGATAGTAATCCCTATATTAATATTAGGGATCATTTAGCAACAATAGTTGAGGTGGTAAAACGATCAAAAAATGGTATTCGCCAGATAATTGATTTAAAGCGTGAATTGAAATTTGACAAAGAGAAAAGACTTTCTAAAGCTTTTATTGCAGCTATTGCATGGTTAGATGAAAAGGAGTACTTAAAGAAAGTGCTTGTAGTATCACCCAAGAATCCTGCCATTAAAATCAGATGTGTAAAATACGTGAAAGATATTCCAGACTCTAAAGGCTCGCCTTCATTTGAGTATGATAGCAATAGCGCGGATGAAGATTCTGTATCAGATAGCAAGGCAGCTTTCGAAGATGAAGACTTAGTCGAAGGTTTAGATAATTTCAATGCGACTGATTTATTACAAAATCAAGGCCTTGTTATGGAAGAGAAAGAGGATGCTGTAAAGAATGAAGTTCTTCTTAATCGATTTTATCCACTTCAAAATCAGACTTATGACATTGCAGATAAGTCTGGCCTTAAAGGAATTTCAACTATGGATGTTGTAAATCGAATTACCGGAAAAGAATTTCAGCGAGCTTTTACCAAATCAAGCGAATATTATTTAGAAAGTGTGGATAAGCAAAAAGAAAATACAGGGGGGTATAGGCTTTTTCGCATATACGATTTTGAGGGAAAGAAGAAGTTTTTTAGGCTGTTCACAGCTCAGAACTTTCAAAAGTTAACAAATGCGGAAGACGAAATATCCGTTCCAAAAGGGTTTGATGAGCTAGGCAAATCTCGTACCGATTTGAAAACTCTCAACGAGGATAATTTCGTCGCACTCAACAACACTGTTAGATTTACAACGGACAGCGATGGACAGGATATATTCTTCTGGCACGGTGAATTAAAAATTCCCCCAAACTCAAAAAAAACTCCGAATAAAAACAAACGGAAGAGGCAGGTTAAAAACAGTACTAATGCTTCTGTTGCAGGAAACATTTCGAATCCCAAAAGGATTAAGCTAGAGCAGCATGTCAGCACTGCACAGGAGCCGAAATCTGCTGAAGATAGTCCAAGTTCAAACGGAGGCACTGTTGTCAAAGGCAAGGTGGTTAACTTCGGCGGCTTTTCTGCCCGCTCTTTGCGTTCACTACAGAGACAGAGAGCCATTTTGAAAGTTATGAATACGATTGGTGGGGTAGCATACCTGAGAGAACAATTTTACGAAAGCGTTTCTAAATATATGGGCTCCACAACGACATTAGATAAAAAGACTGTCCGTGGTGATGTTGATTTGATGGTAGAAAGCGAAAAATTAGGAGCCAGAACAGAGCCTGTATCAGGAAGAAAAATTATTTTTTTGCCCACTGTTGGAGAGGACGCTATCCAAAGGTACATCCTGAAAGAAAAAGATAGTAAAAAAGCAACCTTTACTGATGTTATACATGATACGGAAATATACTTCTTTGACCAAACGGAAAAAAATAGGTTTCACAGAGGAAAGAAATCAGTTGAAAGAATTCGTAAGTTTCAGAACCGCCAAAAGAATGCTAAGATCAAAGCTTCAGATGACGCTATCTCTAAGAAGAGTACGTCGGTCAACGTATCAGATGGAAAGATCAAAAGGAGAGACAAAAAAGTGTCTGCTGGTAGGACAACGGTGGTCGTGGAAAATACTAAAGAAGACAAAACTGTCTATCATGCAGGCACTAAAGATGGTGTTCAGGCTTTAATCAGAGCTGTTGTAGTTACTAAAAGTATTAAAAATGAAATAATGTGGGACAAAATAACAAAATTATTTCCTAATAATTCTTTAGATAACCTAAAAAAGAAATGGACGGCACGGCGAGTAAGAATGGGTCATAGTGGTTGGAGGGCATATGTCGATAAGTGGAAAAAAATGCTCGTTCTAGCCATTAAAAGTGAAAAGATTTCACTGAGGGATGTTGAAGAACTAGATCTTATCAAATTGCTTGATATTTGGACCTCTTTTGATGAAAAGGAAATAAAAAGGCCGCTCTTTCTTTATAAGAACTACGAAGAGAATAGAAAAAAATTTACTCTGGTACGTGATGACACACTTACACATTCTGGCAACGATCTGGCCATGTCTTCTATGATTCAAAGAGAGATCTCTTCTTTAAAAAAAACTTACACTAGAAAGATTTCCGCTTCTACTAAGGACTTATCGAAGAGTCAAAGCGACGATTATATTCGCACAGTGATCCGGTCCATATTAATAGAAAGTCCTTCGACCACTAGAAATGAAATAGAGGCGTTGAAGAACGTTGGAAACGAATCAATAGATAACGTCATCATGGATATGGCTAAGGAAAAGCAAATTTATCTCCATGGCTCAAAACTTGAATGTACTGATACTTTACCAGACATTTTGGAAAATAGAGGAAATTATAAAGATTTTGGTGTAGCTTTTCAGTATAGATGTAAGGTTAATGAATTATTGGAGGCCGGAAACGCTATTGTTATCAATCAAGAGCCGTCCGATATATCCTCTTGGGTTTTAATTGATTTGATTTCGGGAGAGCTATTGAATATGGATGTAATTCCAATGGTGAGAAATGTTCGACCTTTAACGTATACTTCAAGGAGATTTGAAATACGAACATTAACTCCCCCTCTGATTATATATGCCAATTCTCAGACAAAATTGAATACAGCAAGGAAGTCTGCTGTCAAAGTTCCACTGGGCAAACCATTTTCTCGTTTATGGGTGAATGGATCTGGTTCCATTAGGCCAAACATATGGAAGCAGGTAGTTACTATGGTCGTTAACGAAATAATATTTCATCCAGGGATAACATTGAGTAGATTGCAATCTAGGTGTCGTGAAGTACTTTCGCTTCATGAAATATCAGAAATATGCAAATGGCTCCTAGAAAGACAAGTATTAATAACTACTGATTTTGATGGCTATTGGGTCAATCATAATTGGTATTCTATATATGAATCTACATAA"
    #s = "303308882768827988825829665"
    #s = read_fasta("Escherichia_coli_k_12.ASM80076v1.dna.toplevel.fa")
    #g = read_fasta("ecoli_aceb.fa")

    #print(s)

    directory = 'refseq_prok/ncbi-genomes-2018-04-30'
    for fasta_file in os.listdir(directory):
        if fasta_file.endswith('.fna'):
            print(fasta_file)
            s = read_fasta(directory + '/' + fasta_file)
            do_all_reorderings(s)

            m = Mapping(0)
            r = alg1(s, m, no_backtrack=False)
            fs = list(Lyndon.ChenFoxLyndon(r))
            #print("Original string:")
            #print(s)
            #print()
            #print("Reordered alphabet string:")
            #print(r)
            #print()
            #print("Lyndon factorisation after algorithm:")
            #print(" >= ".join(fs))
            #print(m)
            print(len(fs))
            #for f in fs:
            #    print(f[:40])
            print()

            #g_remapped = m.map_string(g)
            #reads = [g_remapped[i:i+100] for i in range(0,len(g_remapped)-100,100)]

            #for read in reads:
            #    print(read[:15])

        
# lyndon factors -> suffix array -> lcp array 
# star algorithm
