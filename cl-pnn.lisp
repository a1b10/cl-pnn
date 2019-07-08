;;;; cl-pnn.lisp

(in-package #:cl-pnn)

#|

import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


class PNN:
    """
    This follows R's pnn package! Tested!
    """
    def __init__(self, x_data = None, y_data = None):
        if not type(x_data) == pd.core.frame.DataFrame \
           and not type(y_data) == pd.core.frame.DataFrame:
            self.N_features = 0
            self.N_samples  = 0
            self.classes    = []
        else:
            self.N_features = x_data.shape[1]
            self.N_samples  = x_data.shape[0]
            self.classes = self.unique(str(x) for x in y_data)
        self.x_data = x_data   # either data itself or pointer to data
        self.y_data = y_data
    
    def unique(self, seq):
        seen = set()
        seen_add = seen.add
        return [x for x in seq if not (x in seen or seen_add(x))]
    
    def diff_square_sum(self, x, w_j):
        """With x and w_j being np.arrays."""
        return sum(diff ** 2 for diff in x - w_j)
    
    def gaussian_prbdf(self, x, w_j, sigma): # prbdf = probability density function
        """With x and w_j being np.arrays."""
        return math.exp( -(self.diff_square_sum(x, w_j)) / (2 * sigma ** 2))
    
    def gaussian_prbdf_group(self, x_df, w_j, sigma):
        """With x_df being a data frame of x values,
        with w_j being np.array."""
        m, p = x_df.shape
        return 1 / m * 1 / ((2 * np.pi) ** (p/2) * (sigma ** p))  * sum([self.gaussian_prbdf(row, w_j, sigma) \
            for row in [x_df.iloc[x, :] for x in range(x_df.shape[0])]])
    
    def predict_probs_single(self, x, sigma):
        df_list = self.x_data.groupby(self.y_data)
        probs = [(name, self.gaussian_prbdf_group(df, x, sigma)) for name, df in df_list ]
        total = sum([x[1] for x in probs])
        probs = [(name, x/total) for name, x in probs]
        sorted_probs = sorted( probs, key=lambda x: -x[1]) # decreasing sort by values
        return sorted_probs
    
    def predict_single(self, x, sigma, value = False):
        sorted_probs = self.predict_probs_single(x, sigma)
        return sorted_probs[0][0] if not value else sorted_probs[0][1]

|#

(ql:quickload :computable-reals)
(use-package :computable-reals)
(setq *PRINT-PREC* 50)
(defconstant +e+ (exp-r 1))

(defun unique (l &optional (acc '()))
  "Return unique list."
  (cond ((null l) (nreverse acc))
	((member (car l) acc) (unique (cdr l) acc))
	(t (unique (cdr l) (cons (car l) acc)))))

(defun unique (l)
  (remove-duplicates l :from-end t))

(defun difference-square-sum (l1 l2)
  "Return sum of squared differences."
  (loop for x in l1
	for y in l2
	sum (expt (-r x y) 2) into res
	finally (return res)))
  
(defun gaussian-probability-density-function (l1 l2 sigma)
  "Return gaussian-probability-density-function of probe l1 and row l2."
  (/r (expt +e+ (- (difference-square-sum l1 l2))) (*r 2 (expt-r sigma 2))))

(defun dim (lol)
  "Return dimensions of a list of list or matrix."
  (values (length lol) (length (elt lol 0))))

#|
        m, p = x_df.shape
        return 1 / m * 1 / ((2 * np.pi) ** (p/2) * (sigma ** p))  * sum([self.gaussian_prbdf(row, w_j, sigma) \
            for row in [x_df.iloc[x, :] for x in range(x_df.shape[0])]])
|#
(defun gaussian-probability-density-function-matrix (l1 x-train sigma)
  "Return probability value for a sample and matrix of train values."
  (multiple-value-bind (m p) (dim x-train)
    (let ((expression (*r m (/r 1 (*r (expt (*r 2 pi) (/r p 2)) (expt sigma p))))))
      (/r 1 (*r expression
		(mapcar #'(lambda (l2) (gaussian-probability-density-function l1 l2 sigma))
			x-train))))))

#|
    def predict_probs_single(self, x, sigma):
        df_list = self.x_data.groupby(self.y_data)
        probs = [(name, self.gaussian_prbdf_group(df, x, sigma)) for name, df in df_list ]
        total = sum([x[1] for x in probs])
        probs = [(name, x/total) for name, x in probs]
        sorted_probs = sorted( probs, key=lambda x: -x[1]) # decreasing sort by values
        return sorted_probs
|#

(defun predict-single-probe-probabilities (x x-train y sigma)
  "Return probabilities for a single probe given train data."
  (multiple-value-bind (x-grouped names) (split-by-group x-train y)
    (let* ((probs (loop for group in x-grouped
				for name in names
				collect (list name (gaussian-probability-density-function-matrix x group sigma))))
	   (total (+ (lambda (x) (elt x 1)) probs))
	   (probabilities (loop for x in probs
				for n in names
				collect (list n (/ x total))))
	   (sorted-probabilities (sort probabilities #'> :key #'cadr)))
      sorted-probabilities)))

(defun split-by-group (data labels &optional (acc '()))
  "Return data (list of list) split by different kind of labels.
   First return value is data split by labels second value is unique labels."
  (let* ((unique-labels (unique labels))
	 (positions     (mapcar #'(lambda (x) (position x unique-labels))
				labels))
	 (result        (make-list (length unique-labels))))
    (loop for line in data
	  for pos in positions
	  do (setf (elt result pos) (cons line (elt result pos))))
    (values result unique-labels))) ;; works as expected!

#|
    def predict_single(self, x, sigma, value = False):
        sorted_probs = self.predict_probs_single(x, sigma)
        return sorted_probs[0][0] if not value else sorted_probs[0][1]

|#

(defun predict-single (sample x-train &optional (valuep nil))
  "Return class or probabilities of sample when running over x-train as PNN."
  (let ((sorted-probabilities (predict-single-probe-probabilities x x-train y sigma)))
    (elt (elt sorted-probabilites 0) (if valuep 1 0))))
