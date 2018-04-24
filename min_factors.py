import Lyndon
import re
import itertools

class Mapping:
    def __init__(self):
        self.alphabet_loc = ord('a')
        self.mapping = {}

    def assign(self,letter):
        if letter not in self.mapping:
            self.mapping[letter] = self.alphabet_loc
            self.alphabet_loc += 1

    def look_up(self,letter):
        return self.mapping.get(letter, -1)

    def map_string(self, string):
        out = []
        for c in string:
            if c not in self.mapping:
                self.assign(c)
            out.append(chr(self.mapping[c]))
        return ''.join(out)
                

    def __str__(self):
        out = []
        out.append("The remapping of characters:\n")
        for (k,v) in self.mapping.items():
            out.append(k + "-->" + chr(v)+"\n")
        return(''.join(out))


def read_fasta(filename):
    with open(filename) as f:
        header = next(f)
        seq = ''.join(map(str.strip,f.readlines()))
    return seq


def exp_parikh_vector(s):
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
    ep = [''.join(map(str, counts[c])) for c in uniques]
    return (uniques, ep)



def ep_lyndon_factors(u_ep):
    """Find Lyndon factors of the exponent parikh vector components"""
    for (u,p) in u_ep:
        factors = list(Lyndon.ChenFoxLyndon(p, co_lex=True))
        yield (u,factors)


# -------------------------------------------------------
# printing

def print_eps(us,ep):
    print("Exponent parikhs:")
    for u,p in zip(us,ep):
        print(u,p)
    print()


def print_fL(fL_prs):
    print("Lyndon factors of eps:")
    for (letter, facs) in fL_prs:
        print(letter, " >= ".join(facs))
    print()

# ---------------------------------------------------------    


def choose_pi(fL_prs):
    """Choose the leftmost pr with the minimal number of factors"""
    (orig, facs) = fL_prs[0]
    min_num = len(facs)
    if len(fL_prs) < 2:
        return (min_num, orig, facs)
    else:
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
        for (e, X_vals) in sorted(exp_X_dict.items(), reverse=True):
            for later_X_val in X_vals[1:]:
                finished = True
                for (x,y) in zip(X_vals[0], later_X_val):
                    if x == y:
                        m.assign(x)
                    else:
                        m.assign(x)
                        m.assign(y)
                        finished = False
                        break
                if not finished:
                    break

            # print(m)
            # now assign any leftover letters in the remaining X components arbitrarily
            all_letters = ''.join(X_vals)
            for letter in all_letters:
                m.assign(letter)

        



def alg1(s, m):
    (us,ep) = exp_parikh_vector(s)
    #print_eps(us,ep)
    fL_prs = list(ep_lyndon_factors(zip(us,ep)))
    #print_fL(fL_prs)

    flag = True
    #while flag: # still more to go
    (num, orig, facs) = choose_pi(fL_prs)
    m.assign(orig)
    Xs = get_Xs_after_char(s, orig)
    assign_to_Xs(m, facs, Xs)
    #print(m)
    reordered = m.map_string(s)
    return reordered


def all_reorders(perms):
    for cs in perms:
        m = Mapping()
        for c in cs:
            m.assign(c)
        yield m

def do_all_reorderings(s):
    smallest = len(s)
    smallest_m = None
    print('\nAll possible reorderings')
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
    print(smallest)
    print(smallest_m)
    r = smallest_m.map_string(s)
    fs = list(Lyndon.ChenFoxLyndon(r))
    #print(len(fs))
    #print(" >= ".join(fs))
    print()

        

if __name__ == "__main__":
    
    #s = "bbbffbbcf"
    #s = "ATGGTACTGACGATTTATCCTGACGAACTCGTACAAATAGTGTCTGATAAAATTGCTTCAAATAAGGGAAGTATGTTCA"
    #s = "ATGGTACTGACGATTTATCCTGACGAACTCGTACAAATAGTGTCTGATAAAATTGCTTCAAATAAGGGAAGTATGTTCATGTCTCATTCTCCTTTTCGGCTCCGTTTAGGTGATAAACGTACTATATTGTGAAAGATTATTTACTAACGACACATTGAAGAAATCACTTTGAATCAGCTGTGGGATATATCTGGTAAATATTTTGATTTGTCTGATAAAAAAGTTAAACAGTTCGTGCTTTCATGCGTGATATTGAAAAAGGACATTGAGGTGTATTGTGATGGTGCTATAACAACTAAAAATGTGACTGATATTATAGGCGACGCTAATCATTCATACTCGGTTGGGATTACTGAGGACAGCCTATGGACATTATTAACGGGATACACAAAAAAGGAGTCAACTATTGGAAATTCTGCATTTGAACTACTTCTCGAAGTTGCCAAATCAGGAGAAAAAGGGATCAATACTATGGATTTGGCGCAGGTAACTGGGCAAGATCCTAGAAGTGTGACTGGACGTATCAAGAAAATAAACCACCTGTTAACAAGTTCACAACTGATTTATAAGGGACACGTCGTGAAGCAATTGAAGCTAAAAAAATTCAGCCATGACGGGGTGGATAGTAATCCCTATATTAATATTAGGGATCATTTAGCAACAATAGTTGAGGTGGTAAAACGATCAAAAAATGGTATTCGCCAGATAATTGATTTAAAGCGTGAATTGAAATTTGACAAAGAGAAAAGACTTTCTAAAGCTTTTATTGCAGCTATTGCATGGTTAGATGAAAAGGAGTACTTAAAGAAAGTGCTTGTAGTATCACCCAAGAATCCTGCCATTAAAATCAGATGTGTAAAATACGTGAAAGATATTCCAGACTCTAAAGGCTCGCCTTCATTTGAGTATGATAGCAATAGCGCGGATGAAGATTCTGTATCAGATAGCAAGGCAGCTTTCGAAGATGAAGACTTAGTCGAAGGTTTAGATAATTTCAATGCGACTGATTTATTACAAAATCAAGGCCTTGTTATGGAAGAGAAAGAGGATGCTGTAAAGAATGAAGTTCTTCTTAATCGATTTTATCCACTTCAAAATCAGACTTATGACATTGCAGATAAGTCTGGCCTTAAAGGAATTTCAACTATGGATGTTGTAAATCGAATTACCGGAAAAGAATTTCAGCGAGCTTTTACCAAATCAAGCGAATATTATTTAGAAAGTGTGGATAAGCAAAAAGAAAATACAGGGGGGTATAGGCTTTTTCGCATATACGATTTTGAGGGAAAGAAGAAGTTTTTTAGGCTGTTCACAGCTCAGAACTTTCAAAAGTTAACAAATGCGGAAGACGAAATATCCGTTCCAAAAGGGTTTGATGAGCTAGGCAAATCTCGTACCGATTTGAAAACTCTCAACGAGGATAATTTCGTCGCACTCAACAACACTGTTAGATTTACAACGGACAGCGATGGACAGGATATATTCTTCTGGCACGGTGAATTAAAAATTCCCCCAAACTCAAAAAAAACTCCGAATAAAAACAAACGGAAGAGGCAGGTTAAAAACAGTACTAATGCTTCTGTTGCAGGAAACATTTCGAATCCCAAAAGGATTAAGCTAGAGCAGCATGTCAGCACTGCACAGGAGCCGAAATCTGCTGAAGATAGTCCAAGTTCAAACGGAGGCACTGTTGTCAAAGGCAAGGTGGTTAACTTCGGCGGCTTTTCTGCCCGCTCTTTGCGTTCACTACAGAGACAGAGAGCCATTTTGAAAGTTATGAATACGATTGGTGGGGTAGCATACCTGAGAGAACAATTTTACGAAAGCGTTTCTAAATATATGGGCTCCACAACGACATTAGATAAAAAGACTGTCCGTGGTGATGTTGATTTGATGGTAGAAAGCGAAAAATTAGGAGCCAGAACAGAGCCTGTATCAGGAAGAAAAATTATTTTTTTGCCCACTGTTGGAGAGGACGCTATCCAAAGGTACATCCTGAAAGAAAAAGATAGTAAAAAAGCAACCTTTACTGATGTTATACATGATACGGAAATATACTTCTTTGACCAAACGGAAAAAAATAGGTTTCACAGAGGAAAGAAATCAGTTGAAAGAATTCGTAAGTTTCAGAACCGCCAAAAGAATGCTAAGATCAAAGCTTCAGATGACGCTATCTCTAAGAAGAGTACGTCGGTCAACGTATCAGATGGAAAGATCAAAAGGAGAGACAAAAAAGTGTCTGCTGGTAGGACAACGGTGGTCGTGGAAAATACTAAAGAAGACAAAACTGTCTATCATGCAGGCACTAAAGATGGTGTTCAGGCTTTAATCAGAGCTGTTGTAGTTACTAAAAGTATTAAAAATGAAATAATGTGGGACAAAATAACAAAATTATTTCCTAATAATTCTTTAGATAACCTAAAAAAGAAATGGACGGCACGGCGAGTAAGAATGGGTCATAGTGGTTGGAGGGCATATGTCGATAAGTGGAAAAAAATGCTCGTTCTAGCCATTAAAAGTGAAAAGATTTCACTGAGGGATGTTGAAGAACTAGATCTTATCAAATTGCTTGATATTTGGACCTCTTTTGATGAAAAGGAAATAAAAAGGCCGCTCTTTCTTTATAAGAACTACGAAGAGAATAGAAAAAAATTTACTCTGGTACGTGATGACACACTTACACATTCTGGCAACGATCTGGCCATGTCTTCTATGATTCAAAGAGAGATCTCTTCTTTAAAAAAAACTTACACTAGAAAGATTTCCGCTTCTACTAAGGACTTATCGAAGAGTCAAAGCGACGATTATATTCGCACAGTGATCCGGTCCATATTAATAGAAAGTCCTTCGACCACTAGAAATGAAATAGAGGCGTTGAAGAACGTTGGAAACGAATCAATAGATAACGTCATCATGGATATGGCTAAGGAAAAGCAAATTTATCTCCATGGCTCAAAACTTGAATGTACTGATACTTTACCAGACATTTTGGAAAATAGAGGAAATTATAAAGATTTTGGTGTAGCTTTTCAGTATAGATGTAAGGTTAATGAATTATTGGAGGCCGGAAACGCTATTGTTATCAATCAAGAGCCGTCCGATATATCCTCTTGGGTTTTAATTGATTTGATTTCGGGAGAGCTATTGAATATGGATGTAATTCCAATGGTGAGAAATGTTCGACCTTTAACGTATACTTCAAGGAGATTTGAAATACGAACATTAACTCCCCCTCTGATTATATATGCCAATTCTCAGACAAAATTGAATACAGCAAGGAAGTCTGCTGTCAAAGTTCCACTGGGCAAACCATTTTCTCGTTTATGGGTGAATGGATCTGGTTCCATTAGGCCAAACATATGGAAGCAGGTAGTTACTATGGTCGTTAACGAAATAATATTTCATCCAGGGATAACATTGAGTAGATTGCAATCTAGGTGTCGTGAAGTACTTTCGCTTCATGAAATATCAGAAATATGCAAATGGCTCCTAGAAAGACAAGTATTAATAACTACTGATTTTGATGGCTATTGGGTCAATCATAATTGGTATTCTATATATGAATCTACATAA"
    #s = "303308882768827988825829665"
    s = read_fasta("Escherichia_coli_k_12.ASM80076v1.dna.toplevel.fa")
    g = read_fasta("ecoli_aceb.fa")
    
    #do_all_reorderings(s)

    m = Mapping()
    r = alg1(s, m)
    fs = list(Lyndon.ChenFoxLyndon(r))
    #print("Original string:")
    #print(s)
    #print()
    #print("Reordered alphabet string:")
    #print(r)
    #print()
    print("Lyndon factorisation after algorithm:")
    #print(" >= ".join(fs))
    print(m)
    print(len(fs))
    for f in fs:
        print(f[:100])
    print()

    g_remapped = m.map_string(g)
    reads = [g_remapped[i:i+100] for i in range(0,len(g_remapped)-100,100)]

    for read in reads:
        print(read[:15])

        
# lyndon factors -> suffix array -> lcp array 
# star algorithm
