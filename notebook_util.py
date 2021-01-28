# Functions
import align.algs
from align import algs
import pickle
import numpy as np
from sklearn.metrics import confusion_matrix


def paired_file_reader(file):
    pair = []
    amino = []
    with open(file) as file:
        for n, line in enumerate(file):
            pair = []
            for name in line.split():
                pair.append(name)
            amino.append(pair)
    return amino


def align_score_gen(pair_list, sim_matrix, gap, align_type, base_path):
    scores = []
    if align_type == 'nw':
        for seq_a, seq_b in pair_list:
            sw = algs.NeedlemanWunsch(sim_matrix, gap)
            sw.align(base_path + seq_a, base_path + seq_b)
            scores.append(sw.score_[0])

    elif align_type == 'sw':
        for seq_a, seq_b in pair_list:
            sw = algs.SmithWaterman(sim_matrix, gap)
            sw.align(base_path + seq_a, base_path + seq_b)
            scores.append(sw.score_[0])
    else:
        raise Exception("align_type must be either 'nw' or 'sw'")
    return scores


def pred_gen(scores):
    pred_list = []
    for thresh in np.sort(scores):
        pred = []
        for value in scores:
            if value >= thresh:
                pred.append(1)
            else:
                pred.append(0)
        pred_list.append(pred)
    return pred_list


def pr_calc(actual, prediction_list):
    tpr, fpr = [], []
    for prediction in prediction_list:
        cm = confusion_matrix(actual, prediction)
        tn, fp, fn, tp = cm.ravel()
        tpr.append(tp/(tp + fn))
        fpr.append(fp/(fp + tn))
    return tpr, fpr


def save_obj(obj, name):
    with open(name + '.pkl', 'wb+') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)