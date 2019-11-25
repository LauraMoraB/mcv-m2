import os
import matplotlib.pyplot as plt
import numpy as np

from pandas import ExcelFile

from pystruct.models import ChainCRF, MultiClassClf
from pystruct.learners import OneSlackSSVM, NSlackSSVM, FrankWolfeSSVM
from sklearn.model_selection import KFold
from sklearn.svm import LinearSVC

from plot_segments import plot_segments

num_segments_per_jacket = 40
add_gaussian_noise_to_features = False
sigma_noise = 0.1
plot_labeling = True
plot_coefficients = True

measures_path = 'more_samples/man_jacket_hand_measures.xls'
segments_path = 'more_samples/segments'

""" 
Load the segments and the groundtruth for all jackets
"""
xl = ExcelFile(measures_path)
sheet = xl.parse(xl.sheet_names[0])
""" be careful, parse() just reads literals, does not execute formulas """
xl.close()

it = sheet.iterrows()
labels_segments = []
segments = []
for row in it:
    ide = row[1]['ide']
    segments.append(np.load(os.path.join(segments_path, ide+'_front.npy'), allow_pickle=True, encoding='latin1'))
    labels_segments.append(list(row[1].values[-num_segments_per_jacket:]))

labels_segments = np.array(labels_segments).astype(int)

"""
Show groundtruth for 3rd jacket
"""
n = 2
plot_segments(segments[n], sheet.ide[n], labels_segments[n])

""" 
Make matrices X of shape (number of jackets, number of features) 
and Y of shape (number of jackets, number of segments) where, 
for all jackets,
    X = select the features for each segment 
    Y = the grountruth label for each segment
"""
Y = labels_segments
num_jackets = labels_segments.shape[0]
num_labels = np.unique(np.ravel(labels_segments)).size

""" CHANGE THIS IF YOU CHANGE NUMBER OF FEATURES """
num_features = 7
X = np.zeros((num_jackets, num_segments_per_jacket, num_features))

for jacket_segments, i in zip(segments, range(num_jackets)):
    for s, j in zip(jacket_segments, range(num_segments_per_jacket)):
        """ set the features """
        X[i, j, 0:num_features] = \
            s.x0norm, s.y0norm, s.x1norm, s.y1norm, \
            (s.x0norm + s.x1norm) / 2., (s.y0norm + s.y1norm) / 2., \
            s.angle / (2 * np.pi)
        """ all possible features at present, see segment.py """
        # s.x0, s.y0, s.x1, s.y1, \
        # s.x0norm, s.y0norm, s.x1norm, s.y1norm, \
        # (s.x0norm+s.x1norm)/2., (s.y0norm+s.y1norm)/2., \
        # np.sqrt((s.x0norm-s.x1norm)**2 + (s.y0norm-s.y1norm)**2), \
        # s.angle, \

print('X, Y done')

""" you can add some noise to the features """
if add_gaussian_noise_to_features:
    print('Noise sigma {}'.format(sigma_noise))
    X = X + np.random.normal(0.0, sigma_noise, size=X.size).reshape(np.shape(X))

"""
DEFINE HERE YOUR GRAPHICAL MODEL AND CHOOSE ONE LEARNING METHOD
(OneSlackSSVM, NSlackSSVM, FrankWolfeSSVM)
"""
model = ChainCRF()
ssvm = FrankWolfeSSVM(model, verbose=1)

svm = LinearSVC()

""" 
Compare SVM with S-SVM doing k-fold cross validation, k=5, see scikit-learn.org 
"""
n_folds = 5
""" with 5 in each fold we have 4 jackets for testing, 19 for training, 
with 23 we have leave one out : 22 for training, 1 for testing"""
scores_crf = np.zeros(n_folds)
scores_svm = np.zeros(n_folds)
wrong_segments_crf = []
wrong_segments_svm = []

kf = KFold(n_splits=n_folds)
fold = 0
for train_index, test_index in kf.split(X):
    print(' ')
    print('train index {}'.format(train_index))
    print('test index {}'.format(test_index))
    print('{} jackets for training, {} for testing'.format(len(train_index), len(test_index)))
    X_train = X[train_index]
    Y_train = Y[train_index]
    X_test = X[test_index]
    Y_test = Y[test_index]

    """ YOUR S-SVM TRAINING CODE HERE """
    ssvm.fit(X_train, Y_train)

    """ LABEL THE TESTING SET AND PRINT RESULTS """
    Y_pred = ssvm.predict(X_test)
    wrong_segments_crf.append(np.sum(Y_pred != Y_test))
    score = ssvm.score(X_test, Y_test)
    scores_crf[fold] = score

    """ figure showing the result of classification of segments for
    each jacket in the testing part of present fold """
    if plot_labeling:
        for ti, pred in zip(test_index, Y_pred):
            print(ti)
            print(pred)
            s = segments[ti]
            plot_segments(s, caption='SSVM predictions for jacket ' + str(ti + 1), labels_segments=pred)

    """ YOUR LINEAR SVM TRAINING CODE HERE """
    svm.fit(X_train.reshape((-1, num_features)), Y_train.reshape((-1)))

    """ LABEL THE TESTING SET AND PRINT RESULTS """
    Y_pred = svm.predict(X_test.reshape((-1, num_features))).reshape((-1, num_segments_per_jacket))
    wrong_segments_svm.append(np.sum(Y_pred != Y_test))
    score = svm.score(X_test.reshape((-1, num_features)), Y_test.reshape((-1)))
    scores_svm[fold] = score

    fold += 1

"""
Global results
"""
total_segments = num_jackets * num_segments_per_jacket
wrong_segments_crf = np.array(wrong_segments_crf)
wrong_segments_svm = np.array(wrong_segments_svm)
print('Results per fold ')
print('Scores CRF : {}'.format(scores_crf))
print('Scores SVM : {}'.format(scores_svm))
print('Wrongs CRF : {}'.format(wrong_segments_crf))
print('Wrongs SVM : {}'.format(wrong_segments_svm))
print(' ')
print('Final score CRF: {}, {} wrong labels in total out of {}'.
      format(1.0 - wrong_segments_crf.sum() / float(total_segments),
             wrong_segments_crf.sum(),
             total_segments))
print('Final score SVM: {}, {} wrong labels in total out of {}'.
      format(1.0 - wrong_segments_svm.sum() / float(total_segments),
             wrong_segments_svm.sum(),
             total_segments))

if plot_coefficients:
    name_of_labels = [
        'neck',
        'left shoulder',
        'outer left sleeve',
        'left wrist',
        'inner left sleeve',
        'left chest',
        'waist',
        'right chest',
        'inner right sleeve',
        'right wrist',
        'outer right sleeve',
        'right shoulder',
    ]

    """ SHOW IMAGE OF THE LEARNED UNARY COEFFICIENTS, size (num_labels, num_features)"""
    """ use matshow() and colorbar()"""
    fig, ax = plt.subplots()
    mat = ax.matshow(np.array(ssvm.w[:num_labels*num_features]).reshape((num_labels, num_features)))
    fig.colorbar(mat, ax=ax)
    ax.set_yticks(np.arange(num_labels))
    ax.set_yticklabels(name_of_labels)
    plt.show(block=False)

    """ SHOW IMAGE OF PAIRWISE COEFFICIENTS size (num_labels, num_labels)"""
    fig, ax = plt.subplots()
    mat = ax.matshow(np.array(ssvm.w[-num_labels**2:]).reshape((num_labels, num_labels)))
    fig.colorbar(mat, ax=ax)
    ax.set_yticks(np.arange(num_labels))
    ax.set_xticks(np.arange(num_labels))
    ax.set_yticklabels(name_of_labels)
    ax.set_xticklabels(name_of_labels, rotation=45, ha='left')
    plt.show(block=False)
