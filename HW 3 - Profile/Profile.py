import numpy as np
import math
from copy import deepcopy

PSEUDOCOUNT = 2


def find_all_amino_acids_in_seqs(seqs):
    amino_acids = []
    for seq in seqs:
        amino_acids_in_seq = list(seq)
        for amino_acid in amino_acids_in_seq:
            if amino_acid not in amino_acids:
                amino_acids.append(amino_acid)

    if '-' not in amino_acids:
        amino_acids.append('-')

    return amino_acids


def create_np_2d_array_from_seqs(seqs):
    seqs_2d_arr = []
    for seq in seqs:
        amino_acids_in_seq = list(seq)
        seqs_2d_arr.append(amino_acids_in_seq)

    return np.array(seqs_2d_arr)


def count_amino_acid_in_column(seqs_2d_arr, amino_acid, column_number):
    # print(seqs_2d_arr)
    return np.count_nonzero(seqs_2d_arr[:, column_number] == amino_acid)


def count_amino_acid_in_MSA(seqs_2d_arr, amino_acid):
    return np.count_nonzero(seqs_2d_arr == amino_acid)


def create_PSSM(all_amino_acids, seqs_2d_arr):
    number_of_columns = len(seqs_2d_arr[0]) + 1
    number_of_rows = len(all_amino_acids) + 1
    N = len(seqs_2d_arr)  # number_of_seqs
    B = len(all_amino_acids)

    PSSM = np.zeros((number_of_rows, number_of_columns))
    for j in range(number_of_columns):
        PSSM[0][j] = j

    # Convert multiple alignment to a raw frequency table
    for i in range(1, number_of_rows):
        for j in range(number_of_columns):
            # amino_acid = all_amino_acids[i - 1]
            amino_acid = all_amino_acids[i - 1]

            if j == number_of_columns - 1:
                # score = (count_amino_acid_in_MSA(seqs_2d_arr, amino_acid)) / (
                #         N * len(seqs_2d_arr[0]))
                row_sum = PSSM.sum(axis=1)
                score = row_sum[i] / j

            else:
                frequency = count_amino_acid_in_column(seqs_2d_arr, amino_acid, j)
                score = (frequency + PSEUDOCOUNT) / (N + B * PSEUDOCOUNT)

            PSSM[i, j] = score

    # print(PSSM)
    # Normalize the values by dividing them by overall frequency and
    # then convert the values to log to base of 2
    for i in range(1, number_of_rows):
        for j in range(number_of_columns - 1):
            overall_freq = PSSM[i, number_of_columns - 1]
            PSSM[i, j] = PSSM[i, j] / overall_freq
            PSSM[i, j] = math.log2(PSSM[i, j])

    return PSSM


def calculate_score_seq(PSSM, seq, all_amino_acids_in_seqs):
    seq_length = len(seq)
    seq_score = 0
    for i in range(seq_length):
        amino_acid = seq[i]
        amino_acid_index = all_amino_acids_in_seqs.index(amino_acid)
        seq_score += PSSM[amino_acid_index + 1, i]

    return seq_score


def add_gap_to_seq(seq):
    arr = []
    for i in range(len(seq) + 1):
        seq_copy = deepcopy(seq)
        seq_with_gap = seq_copy[:i] + '-' + seq_copy[i:]
        arr.append(seq_with_gap)

    return arr


def generate_all_possible_seq_with_gap(seq, max_length):
    generated_seq_with_gap_arr = add_gap_to_seq(seq)
    generated_seq_with_gap_arr = generated_seq_with_gap_arr

    while len(generated_seq_with_gap_arr[0]) < max_length:
        new_generated_seq_with_gap = []
        for new_seq in generated_seq_with_gap_arr:
            generated_seq_with_gap = add_gap_to_seq(new_seq)
            new_generated_seq_with_gap = new_generated_seq_with_gap + generated_seq_with_gap

        generated_seq_with_gap_arr = list(set(new_generated_seq_with_gap))

    return generated_seq_with_gap_arr


def search_in_long_seq(long_seq, PSSM, all_amino_acids_in_seqs):
    length_of_seq = len(PSSM[0]) - 1
    long_seq_arr = list(long_seq)
    length_of_long_seq = len(long_seq_arr)
    # print(length_of_long_seq)

    best_sub_seq = long_seq[: length_of_seq]
    max_score = calculate_score_seq(PSSM, best_sub_seq, all_amino_acids_in_seqs)

    for l in range(length_of_seq, 2, -1):
        seqs_without_gap = []
        calculated_seq = set([])
        for i in range(0, length_of_long_seq - l + 1):
            # print(calculated_seq)
            # don't have gap
            if l == length_of_seq:
                seq = long_seq[i: i + l]
                if seq not in calculated_seq:
                    calculated_seq.add(seq)
                    seq_score = calculate_score_seq(PSSM, seq, all_amino_acids_in_seqs)
                    # check score
                    if seq_score > max_score:
                        max_score = seq_score
                        best_sub_seq = seq

            else:
                seq_without_gap = long_seq[i: i + l]
                if seq_without_gap not in seqs_without_gap:
                    seqs_without_gap.append(seq_without_gap)
                    all_possible_seq_with_gap = generate_all_possible_seq_with_gap(seq_without_gap, length_of_seq)
                    # all_possible_seq_with_gap = add_gaps(seq_without_gap, length_of_seq - len(seq_without_gap))

                    for possible_seq in set(all_possible_seq_with_gap):
                        if possible_seq not in calculated_seq:
                            seq_score = calculate_score_seq(PSSM, possible_seq, all_amino_acids_in_seqs)
                            if seq_score > max_score and len(possible_seq) == length_of_seq:
                                max_score = seq_score
                                best_sub_seq = possible_seq

                    calculated_seq = set.union(calculated_seq, all_possible_seq_with_gap)
    return best_sub_seq, max_score


if __name__ == '__main__':
    number_of_seqs = int(input())
    seqs = []
    for i in range(number_of_seqs):
        seq = input()
        seqs.append(seq)

    long_seq = input()

    all_amino_acids_in_seqs = find_all_amino_acids_in_seqs(seqs)
    seqs_2d_arr = create_np_2d_array_from_seqs(seqs)
    PSSM = create_PSSM(all_amino_acids_in_seqs, seqs_2d_arr)
    # print(PSSM)
    # print(all_amino_acids_in_seqs)

    print(search_in_long_seq(long_seq, PSSM, all_amino_acids_in_seqs)[0])
