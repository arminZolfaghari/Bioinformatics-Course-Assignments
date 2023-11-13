import numpy as np
from copy import deepcopy


def global_align(x, y, s_match=3, s_mismatch=-1, s_gap=-2):
    A = []
    for i in range(len(y) + 1):
        A.append([0] * (len(x) + 1))
    for i in range(len(y) + 1):
        A[i][0] = s_gap * i
    for i in range(len(x) + 1):
        A[0][i] = s_gap * i
    for i in range(1, len(y) + 1):
        for j in range(1, len(x) + 1):
            A[i][j] = max(
                A[i][j - 1] + s_gap,
                A[i - 1][j] + s_gap,
                A[i - 1][j - 1] + (s_match if (y[i - 1] == x[j - 1] and y[i - 1] != '-') else 0) + (
                    s_mismatch if (y[i - 1] != x[j - 1] and y[i - 1] != '-' and x[j - 1] != '-') else 0) + (
                    s_gap if (y[i - 1] == '-' or x[j - 1] == '-') else 0)
            )
    align_X = ""
    align_Y = ""
    i = len(x)
    j = len(y)
    while i > 0 or j > 0:
        current_score = A[j][i]
        if i > 0 and j > 0 and (
                ((x[i - 1] == y[j - 1] and y[j - 1] != '-') and current_score == A[j - 1][i - 1] + s_match) or
                ((y[j - 1] != x[i - 1] and y[j - 1] != '-' and x[i - 1] != '-') and current_score == A[j - 1][
                    i - 1] + s_mismatch) or
                ((y[j - 1] == '-' or x[i - 1] == '-') and current_score == A[j - 1][i - 1] + s_gap)
        ):
            align_X = x[i - 1] + align_X
            align_Y = y[j - 1] + align_Y
            i = i - 1
            j = j - 1
        elif i > 0 and (current_score == A[j][i - 1] + s_gap):
            align_X = x[i - 1] + align_X
            align_Y = "-" + align_Y
            i = i - 1
        else:
            align_X = "-" + align_X
            align_Y = y[j - 1] + align_Y
            j = j - 1
    return (align_X, align_Y, A[len(y)][len(x)])


def calculate_distance_matrix(seqs):
    distance_matrix = [['-' for i in range(len(seqs))] for j in range(len(seqs))]
    for i in range(len(seqs)):
        for j in range(len(seqs)):
            if i != j:
                distance_matrix[i][j] = global_align(seqs[i], seqs[j])[2]
            else:
                distance_matrix[i][j] = 0

    return distance_matrix


def find_order_sequence(seqs):
    distance_matrix = calculate_distance_matrix(seqs)
    # print(distance_matrix)
    new_distance_matrix = np.array(distance_matrix)
    sum_rows = new_distance_matrix.sum(axis=1)
    center_index = (-sum_rows).argsort()[0]
    order_sequences = new_distance_matrix[center_index]
    # print(order_sequences)
    order_sequences = (-order_sequences).argsort()

    return center_index, order_sequences


def find_difference(old_center_seq, new_center_seq):
    gaps_index = []
    if len(old_center_seq) != len(new_center_seq):
        number_of_gaps = len(new_center_seq) - len(old_center_seq)
        i, j = 0, 0
        while j < len(new_center_seq):
            # print(i, j)
            # print(old_center_seq)
            # print(new_center_seq)
            if i >= len(old_center_seq) and new_center_seq[j] == '-':
                gaps_index.append(j)
                j += 1
            elif old_center_seq[i] != new_center_seq[j] and new_center_seq[j] == '-':
                gaps_index.append(j)
                j += 1
            else:
                i += 1
                j += 1

    return gaps_index


def update_MSA_with_gap(MSA, old_center_seq, new_center_seq):
    MSA_copy = deepcopy(MSA)
    gaps_index = find_difference(old_center_seq, new_center_seq)
    for sequence_number in MSA_copy:
        sequence_before_update = list(MSA_copy[sequence_number])
        for index in gaps_index:
            sequence_before_update.insert(index, '-')

        sequence_after_update = ''.join(sequence_before_update)
        MSA_copy[sequence_number] = sequence_after_update

    return MSA_copy


def cstar_algorithm(sequences):
    center_index, order_sequences = find_order_sequence(sequences)
    MSA = {center_index: sequences[center_index]}
    for seq_index in order_sequences:
        aligned_seq, new_center_seq, score = global_align(sequences[seq_index], MSA[center_index])

        MSA = update_MSA_with_gap(MSA, MSA[center_index], new_center_seq)
        MSA[seq_index] = aligned_seq
        MSA[center_index] = new_center_seq

    return MSA, center_index


def sort_final_MSA(MSA):
    for i in range(len(MSA)):
        seq = ''.join(MSA[i])
        print(seq)


def calculate_score(MSA):
    length_seq = len(MSA[0])
    final_score = 0
    number_of_sequences = len(MSA)
    for seq1_index in range(number_of_sequences):
        seq1 = list(MSA[seq1_index])
        for seq2_index in range(seq1_index + 1, number_of_sequences):
            if seq1_index != seq2_index:
                seq2 = MSA[seq2_index]
                for i in range(length_seq):
                    if seq1[i] == seq2[i] and seq1[i] != '-':  # match
                        final_score += 3
                    elif seq1[i] == seq2[i] and seq1[i] == '-':  # gap, gap
                        final_score += 0
                    elif (seq1[i] == '-' and seq2[i] != '-') or (seq1[i] != '-' and seq2[i] == '-'):  # gap
                        final_score += -2
                    elif seq1[i] != '-' and seq2[i] != '-' and seq1[i] != seq2[i]:  # mismatch
                        final_score += -1

    return final_score


def find_column_for_block(MSA):
    number_of_sequences = len(MSA)
    length_seq = len(MSA[0])

    MSA_arr = []
    for i in range(number_of_sequences):
        MSA_arr.append(list(MSA[i]))

    MSA_arr = np.array(MSA_arr)

    status = []

    for column in MSA_arr.T:
        first_char = column[0]
        number_of_this_char = np.count_nonzero(column == first_char)
        if number_of_this_char != number_of_sequences:
            status.append(1)
        else:
            status.append(-1)

    if status[0] == 1 and status[1] == -1:
        status[0] = -1
    if status[-1] == 1 and status[-2] == -1:
        status[-1] = -1
    for i in range(1, len(status) - 1):
        if status[i] == 1 and status[i - 1] == status[i + 1] == -1:
            status[i] = -1

    block_index = []
    start = -1
    end = -1
    for i in range(len(status)):
        if status[i] == 1:
            if start == -1:
                start = i
            if i + 1 == len(status):
                end = i
                if end - start >= 1:
                    block_index.append((start, end))
        if status[i] == -1:
            if start != -1:
                end = i - 1
                if end - start >= 1:
                    block_index.append((start, end))
                start = -1

    return block_index, MSA_arr


def remove_gap_in_seqs(MSA_arr):
    number_of_seqs = len(MSA_arr)
    length_seq = len(MSA_arr[0])

    sequences = []
    for i in range(number_of_seqs):
        seq = ""
        for j in range(length_seq):
            if MSA_arr[i, j] != '-':
                seq += MSA_arr[i, j]
        sequences.append(seq)

    return sequences


def create_arr_from_dict(MSA_dict):
    number_of_seqs = len(MSA_dict)
    arr = [[] for i in range(number_of_seqs)]

    for seq_index in MSA_dict:
        arr[seq_index] = list(MSA_dict[seq_index])

    MSA_arr = np.array(arr)
    return MSA_arr


def get_start_end(r):
    (start, end) = r
    # print(start, end)
    distance = end - start
    possible_start_end = []
    for d in range(distance, 0, -1):
        for s in range(start, end + 1, 1):
            e = s + d
            if e > end:
                break
            possible_start_end.append((s, e))

    return possible_start_end


# print(get_start_end((10, 11)))
# exit()


if __name__ == '__main__':
    arr = []
    n = input()  # number of sequences
    sequences = []
    for i in range(int(n)):
        sequence = input()
        sequences.append(sequence)

    # TODO remove index
    MSA, center_index = cstar_algorithm(sequences)
    pre_score = calculate_score(create_arr_from_dict(MSA))
    print(pre_score)
    sort_final_MSA(MSA)

    # exit()

    flag = True
    while True:
        block_index, MSA_arr = find_column_for_block(MSA)
        if len(block_index) == 0:
            flag = False
            break

        # print(block_index)

        for i in range(len(block_index)):
            start, end = block_index[i][0], block_index[i][1]

            block_seq = MSA_arr[:, start: end + 1]
            sequences = remove_gap_in_seqs(block_seq)
            sub_MSA, sub_center_index = cstar_algorithm(sequences)
            current_score = calculate_score(sub_MSA)
            if start != 0:
                part1 = MSA_arr[:, : start]
                current_score += calculate_score(part1)

            if end + 1 != len(MSA_arr[0]):
                part2 = MSA_arr[:, end + 1:]
                current_score += calculate_score(part2)

            sub_MSA_arr = create_arr_from_dict(sub_MSA)

            # print(current_score, pre_score)
            if current_score > pre_score:
                if start != 0 and end + 1 != len(MSA_arr[0]):
                    MSA = np.concatenate((part1, sub_MSA_arr, part2), axis=1)
                elif start == 0:
                    MSA = np.concatenate((sub_MSA_arr, part2), axis=1)
                elif end + 1 == len(MSA_arr[0]):
                    MSA = np.concatenate((part1, sub_MSA_arr), axis=1)
                pre_score = current_score
                break
                # print(pre_score)
                # print(MSA)

        if i + 1 == len(block_index):
            flag = False
            break

if flag == False:
    print(pre_score)
    sort_final_MSA(MSA)
