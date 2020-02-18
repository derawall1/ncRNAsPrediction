
import sys
import re
import time
import os

from features_extraction_utility.util import frequency
from features_extraction_utility.util import get_data
from features_extraction_utility.data import index_list


def make_kmer_list(k, sequence_type):
    if k < 0:
        print("Error, k must be an inter and larger than 0.")

    kmers = []
    for i in range(1, k + 1):
        if len(kmers) == 0:
            kmers = list(sequence_type)
        else:
            new_kmers = []
            for kmer in kmers:
                for c in sequence_type:
                    new_kmers.append(kmer + c)
            kmers = new_kmers

    return kmers


def find_revcomp(sequence, revcomp_dictionary):
    # Save time by storing reverse complements in a hash.
    if sequence in revcomp_dictionary:
        return revcomp_dictionary[sequence]

    # Make a reversed version of the string.
    rev_sequence = list(sequence)
    rev_sequence.reverse()
    rev_sequence = ''.join(rev_sequence)

    return_value = ""
    for letter in rev_sequence:
        if letter == "A":
            return_value += "T"
        elif letter == "C":
            return_value += "G"
        elif letter == "G":
            return_value += "C"
        elif letter == "T":
            return_value += "A"
        elif letter == "N":
            return_value += "N"
        else:
            error_info = ("Unknown DNA character (%s)\n" % letter)
            sys.exit(error_info)

    # Store this value for future use.
    revcomp_dictionary[sequence] = return_value

    return return_value


def _cmp(a, b):
    return (a > b) - (a < b)


def make_revcomp_kmer_list(kmer_list):
    revcomp_dictionary = {}
    new_kmer_list = [kmer for kmer in kmer_list if _cmp(kmer, find_revcomp(kmer, revcomp_dictionary)) <= 0]
    return new_kmer_list


def make_kmer_vector(k, sequence_type, filename, revcomp=False):
    """Generate kmer vector."""
    #print(filename)
    with open(filename) as f:
        seq_list = get_data(f, sequence_type=sequence_type)

        if revcomp and re.search(r'[^acgtACGT]', ''.join(sequence_type)) is not None:
            sys.exit("Error, Only DNA sequence can be reverse compliment.")

        vector = []
        kmer_list = make_kmer_list(k, sequence_type)
        for seq in seq_list:
            count_sum = 0

            # Generate the kmer frequency dict.
            kmer_count = {}
            for kmer in kmer_list:
                temp_count = frequency(seq, kmer)
                if not revcomp:
                    if kmer not in kmer_count:
                        kmer_count[kmer] = 0
                    kmer_count[kmer] += temp_count
                else:
                    rev_kmer = find_revcomp(kmer, {})
                    if kmer <= rev_kmer:
                        if kmer not in kmer_count:
                            kmer_count[kmer] = 0
                        kmer_count[kmer] += temp_count
                    else:
                        if rev_kmer not in kmer_count:
                            kmer_count[rev_kmer] = 0
                        kmer_count[rev_kmer] += temp_count

                count_sum += temp_count
            #print(count_sum)
            # Normalize.
            if not revcomp:
                count_vec = [kmer_count[kmer] for kmer in kmer_list]
            else:
                revc_kmer_list = make_revcomp_kmer_list(kmer_list)
                count_vec = [kmer_count[kmer] for kmer in revc_kmer_list]
            count_vec = [round(float(e)/count_sum, 8) for e in count_vec]

            vector.append(count_vec)

    return vector
def make_kmer_list_from_seq(k, sequence_type, seq, revcomp=False):
    """Generate kmer list."""
    if sequence_type == 'DNA':
        sequence_type = index_list.DNA
    elif sequence_type == 'RNA':
        sequence_type = index_list.RNA
    elif sequence_type == 'Protein':
        sequence_type = index_list.PROTEIN
    if revcomp and re.search(r'[^acgtACGT]', ''.join(sequence_type)) is not None:
        sys.exit("Error, Only DNA sequence can be reverse compliment.")
    
    vector = []
    kmer_list = make_kmer_list(k, sequence_type)
    #print(kmer_list)
    count_sum = 0

    # Generate the kmer frequency dict.
    kmer_count = {}
    for kmer in kmer_list:
        temp_count = frequency(seq, kmer)
        if not revcomp:
            if kmer not in kmer_count:
                kmer_count[kmer] = 0
            kmer_count[kmer] += temp_count
        else:
            rev_kmer = find_revcomp(kmer, {})
            if kmer <= rev_kmer:
                if kmer not in kmer_count:
                    kmer_count[kmer] = 0
                kmer_count[kmer] += temp_count
            else:
                if rev_kmer not in kmer_count:
                    kmer_count[rev_kmer] = 0
                kmer_count[rev_kmer] += temp_count

        count_sum += temp_count
    #print(count_sum)
    # Normalize.
    if not revcomp:
        count_vec = [kmer_count[kmer] for kmer in kmer_list]
    else:
        revc_kmer_list = make_revcomp_kmer_list(kmer_list)
        count_vec = [kmer_count[kmer] for kmer in revc_kmer_list]
    #print(count_vec)
    count_vec = [round(float(e)/count_sum, 8) for e in count_vec]

    vector.append(count_vec)

    return vector

def main(args):
    # Set revcomp parameter.
    if args.r != 1:
        args.r = False
    elif args.r == 1 and args.sequence_type != 'DNA':
        print("Error, the -r parameter can only be used in DNA.")
    elif args.r == 1 and args.sequence_type == 'DNA':
        args.r = True

    # Set sequence_type parameter.
    if args.sequence_type == 'DNA':
        args.sequence_type = index_list.DNA
    elif args.sequence_type == 'RNA':
        args.sequence_type = index_list.RNA
    elif args.sequence_type == 'Protein':
        args.sequence_type = index_list.PROTEIN

    res = make_kmer_vector(k=args.k, sequence_type=args.sequence_type, filename=args.inputfile, revcomp=args.r)

    # Write correspond res file.
    if args.f == 'svm':
        from util import write_libsvm
        write_libsvm(res, [args.l] * len(res), args.outputfile)
    elif args.f == 'tab':
        from util import write_tab
        write_tab(res, args.outputfile)
    elif args.f == 'csv':
        from util import write_csv
        write_csv(res, args.outputfile)


if __name__ == '__main__':
    import argparse
    from argparse import RawTextHelpFormatter

    parse = argparse.ArgumentParser(description="This is kmer module for generate kmer vector.",
                                    formatter_class=RawTextHelpFormatter)
    parse.add_argument('inputfile',
                       help="The input file in FASTA format.")
    parse.add_argument('outputfile',
                       help="The output file stored results.")
    parse.add_argument('sequence_type', choices=['DNA', 'RNA', 'Protein'],
                       help="The sequence type.")

    parse.add_argument('-k', type=int, default=2,
                       help="The k value of kmer.")
    parse.add_argument('-r', default=0, type=int, choices=[1, 0],
                       help="Whether consider the reverse complement or not.\n"
                            "1 means True, 0 means False. (default = 0)")
    parse.add_argument('-f', default='tab', choices=['tab', 'svm', 'csv'],
                       help="The output format (default = tab).\n"
                            "tab -- Simple format, delimited by TAB.\n"
                            "svm -- The libSVM training data format.\n"
                            "csv -- The format that can be loaded into a spreadsheet program.")
    parse.add_argument('-l', default='+1', choices=['+1', '-1'],
                       help="The libSVM output file label.")

    args = parse.parse_args()

    print("Calculating...")
    start_time = time.time()
    main(args)
    print("Done.")
    print("Used time: %ss" % (time.time() - start_time))