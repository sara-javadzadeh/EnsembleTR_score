import numpy as np
import sys
import copy
import math
import itertools
import pickle
from pprint import pprint


# Function used for filtering VNTRs
def get_normalized_sequence_similarity(a, b):
    # Compute Hamming distance between two short sequences.
    sequences_len = len(a)
    if len(a) != len(b):
        print("Warning: Computing sequence similarity between two sequences " +\
              "of different lengths: {} and {}".format(len(a), len(b)))
        sequences_len = min(len(a), len(b))
    similarity_score = 0
    for idx in range(sequences_len):
        # "M" character represent the masked character that matches
        # any other character for the purpose of this filter.
        if a[idx] == b[idx] or a[idx] == "M" or b[idx] == "M":
            similarity_score += 1
    return similarity_score / sequences_len

# Function used for filtering VNTRs
def get_motif_complexity_score(motif):
    self_match_score = 0
    for idx_in_mask in range(len(motif)):
        # Masking a single character at a time and computing the
        # max similarity score among all masked motifs.
        # Masked character matches any other character.
        masked_motif = motif[:idx_in_mask] + "M" + motif[idx_in_mask + 1:]
        # Creating a rolling window to compare the masked motif with itself.
        motif_window = masked_motif + masked_motif
        for idx in range(1, len(masked_motif)):
            end_idx_in_window = idx + len(masked_motif)
            # Compute the max score among all possible positions
            # when sliding the motif within the motif window to compare.
            self_match_score = max(self_match_score,
                                   get_normalized_sequence_similarity(masked_motif,
                                                         motif_window[idx:end_idx_in_window]))

    return self_match_score


# Function used for filtering VNTRs with variable number of masked characters based on motif length
def get_motif_complexity_score_flexible_masking(motif, num_masked_positions=None):
    self_match_score = 0
    if len(motif) > 40:
        # Long motif. These computations will take along time.
        # Compute a single character masking instead.
        score = get_motif_complexity_score(motif)
        print("motif map for > 40 for motif {} len {} score {}".format(motif, len(motif), score))
        return score
    if num_masked_positions is None:
        num_masked_positions = math.ceil(len(motif)/10.0)
    all_combinations = itertools.combinations(range(len(motif)), num_masked_positions)
    for combination in all_combinations:
        masked_motif = motif
        for idx in combination:
            masked_motif = masked_motif[:int(idx)] + "M" + masked_motif[int(idx)+1:]
        motif_window = masked_motif + masked_motif
        for idx in range(1, len(masked_motif)):
            end_idx_in_window = idx + len(masked_motif)
            # Compute the max score among all possible positions
            # when sliding the motif within the motif window to compare.
            self_match_score = max(self_match_score,
                                   get_normalized_sequence_similarity(masked_motif,
                                                         motif_window[idx:end_idx_in_window]))
            if self_match_score == 1.0:
                # This is the maximum possible score. No need to search for other combinations.
                return self_match_score

    return self_match_score


# Function used for filtering VNTRs
def isValidVNTR(motif, motif_complexity_threshold, load_motif_complexity_scores, motif_complexity_score_map):
    if load_motif_complexity_scores and motif in motif_complexity_score_map:
        score = motif_complexity_score_map[motif]
    elif motif not in motif_complexity_score_map:
        score = get_motif_complexity_score_flexible_masking(motif)
        motif_complexity_score_map[motif] = score
    else:
        score = get_motif_complexity_score_flexible_masking(motif)
        motif_complexity_score_map[motif] = score
    if score > motif_complexity_threshold:
        # VNTR is very much STR like. Skip this VNTR.
        return False, score
    return True, score

if __name__ == "__main__":
    motif = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTAA"
    print("for motif {} with len {}".format(motif, len(motif)))
    print("normalized score is ", get_motif_complexity_score_flexible_masking(motif, num_masked_positions=1))

    motif = "ATGATG"
    print("for motif {} with len {}".format(motif, len(motif)))
    print("normalized score is ", get_motif_complexity_score_flexible_masking(motif))
