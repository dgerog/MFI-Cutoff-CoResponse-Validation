from itertools import permutations
import math

from mfiStudy import *

# in file
FILE_INPUT = [
    #'../../datasets/EVANGELISMOS-LS1A04-10-and-11_UNIQUE_ANONYMOUS_WITHOUT-HLA-TYPING.xlsx',
    #'../../datasets/Gennhmatas-LS1_LOT-10-11_2017-2018-CLEANED-Anonymous.xls',
    '../../datasets/ippokratio-LS1_ALL.xls'
]

STORAGE_DIR = [
    #'results/evangelismos-ls1/',
    #'results/gennimatas-ls1/',
    'results/ippokratio-ls1'
]

# Cross Validation Iterations
ITERATIONS = 10
TRAIN_PCT  = .8

# skip these rows (if any)
LINES_TO_SKIP = [
    1, 
    1, 
    1
]

# which columns to use in order to extract the input vectors
COLS_TO_USE = [
    "AB:DT",
    "N:DF",
    "U:DM"
]

# which column stores the patient ID
LABELS_ID_COL = None

# which column is used as control point (thresolding)
CONTROL_ID_COL = None

# characters to remove from header title
REMOVE_CHARACTERS = [";", ".", "-", ","]

# 
# Start parsing data 
#
study = mfiStudy()

for i in range(0, len(STORAGE_DIR)):
    #load data
    study.readTable(FILE_INPUT[i], COLS_TO_USE[i], LINES_TO_SKIP[i], True, LABELS_ID_COL, CONTROL_ID_COL, _cleanCharacters=REMOVE_CHARACTERS)
    
    # Co-Response Validation
    score = study.coResponsesValidation(_ColBasis=1, _ColControl=[15,28,34], _allInd=None)

    print('%s: %.2d%%' % (FILE_INPUT[i], 100*score))