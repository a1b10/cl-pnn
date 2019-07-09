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

;; (defun unique (l &optional (acc '()))
;;   "Return unique list."
;;   (cond ((null l) (nreverse acc))
;; 	((member (car l) acc) (unique (cdr l) acc))
;; 	(t (unique (cdr l) (cons (car l) acc)))))

(defun unique (l &key (test #'string=))
  "Return unique list."
  (remove-duplicates l :from-end t :test test))

(defun difference-square-sum (l1 l2)
  "Return sum of squared differences."
  (loop for x in l1
	for y in l2
	sum (expt (-r x y) 2) into res
	finally (return res)))
  
(defun gaussian-probability-density-function (l1 l2 sigma)
  "Return gaussian-probability-density-function of probe l1 and row l2."
  (/r (expt +e+ (- (difference-square-sum l1 l2))) (*r 2 (expt-r sigma 2))))

(defun dim (lol &key (contains-header-p t) (contains-row-names-p t))
  "Return dimensions of a list of list or matrix."
  (values (if contains-header-p
              (1- (length lol))
	      (length lol))
	  (if contains-row-names-p
              (1- (length (elt lol 0)))
	      (length (elt lol 0)))
	  contains-header-p
	  contains-row-names-p))

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
	   (total (apply #'+ (mapcar (lambda (x) (elt x 1)) probs)))
	   (probabilities (loop for x in probs
				for n in names
				collect (list n (/ x total))))
	   (sorted-probabilities (sort probabilities #'> :key #'cadr)))
      sorted-probabilities)))

(defun split-by-group (data labels &key (test #'string=))
  "Return data (list of list) split by different kind of labels.
   First return value is data split by labels second value is unique labels."
  (let* ((unique-labels (unique labels :test test))
	 (positions     (mapcar #'(lambda (x) (position x unique-labels :test test))
				labels))
	 (result        (make-list (length unique-labels))))
    (loop for line in data
	  for pos in positions
	  do (setf (elt result pos) (cons line (elt result pos))))
    (values result unique-labels))) ;; works as expected!

(defun split-by-group-indexes (labels &key (test #'string=) (levels nil) (alphabetical-p t) (alpha-test #'string-lessp))
  "Return indexes of the groups in labels. If not levels manually given, and if 
   alphabetical-p is not t, then the levels would be sorted alphabetically, otherwise
   in order of their occurrence."
  (let* ((levels (if (null levels)
		     (unique labels :test test)
		     levels))
         (levels (if alphabetical-p
		     (sort levels alpha-test)
		     levels))
	 (positions (mapcar #'(lambda (x) (position x levels :test test))
			    labels))
	 (result    (make-list (length levels))))
    (loop for i from 0 to (length labels)
	  for pos in positions
	  do (setf (elt result pos) (cons i (elt result pos))))
    (values result levels)))
#|
    def predict_single(self, x, sigma, value = False):
        sorted_probs = self.predict_probs_single(x, sigma)
        return sorted_probs[0][0] if not value else sorted_probs[0][1]

|#

(defun predict-single (sample x-train y sigma &optional (valuep nil))
  "Return class or probabilities of sample when running over x-train as PNN."
  (let ((sorted-probabilities (predict-single-probe-probabilities sample x-train y sigma)))
    (elt (elt sorted-probabilities 0) (if valuep 1 0))))






;; read tab delimited files
;; https://quickref.common-lisp.net/cl-csv.html#Introduction

(ql:quickload :cl-csv)
(use-package :cl-csv)

(defun get-row (matrix i)
  "Return matrix row as numbers and row name as second value."
  (let* ((row (elt matrix i))
	 (counts (mapcar #'parse-integer (cdr row)))
	 (name (car row)))
    (values counts name)))


(defun extract-column-by-position (lol idx)
  "Return content of column in lol by idx position."
  (mapcar (lambda (l) (elt l idx)) lol))

(defun get-column-by-header-name (matrix col-name)
  "Return content of col-name; first row of matrix is header!"
  ;; without header this wouldn't be possible so header is given!
  (let* ((colnames (elt matrix 0))
	 (col-idx  (position col-name colnames :test #'string=)))
    (extract-column-by-position (cdr matrix) (1+ col-idx))))

(defun replace-na (l &key (by "other"))
  "Return l with \"NA\" by 'by' value."
  (substitute by "NA" l :test #'string=)) ;; works!


(defun categorical-to-nums (l &key (levels nil) (alphabetical-p t))
  "Return categories by numbers - if levels not given by unique categories."
  (let* ((categories (if (null levels)
			 (unique l)
			 levels))
	 (categories (if alphabetical-p
			 (sort categories #'string-lessp) ;; case-insensitive alphabetical sort
			 categories))
	 (nums (mapcar (lambda (x) (position x categories :test #'string=))
		       l)))
    (values nums categories)))

(defun categorical-to-one-hot (l &key (levels nil) (alphabetical-p t))
  "Return categories as one-hot encoding."
  (multiple-value-bind (nums categories) (categorical-to-nums l :levels levels :alphabetical-p alphabetical-p)
    (let* ((length-categories (length categories))
	   (one-hot-list (loop for x in nums
			       collect (let ((l (loop for i from 0 to length-categories
						      collect 0)))
					 (setf (elt l x) 1)
					 l))))
      (values one-hot-list categories))))

(defun arg-max (one-hot-list)
  "Return positions of 1 in one-hot-list."
  (mapcar (lambda (l) (position 1 l)) one-hot-list))


(defun one-hot-to-categories (one-hot-list categories)
  "Return categories from one-hot-encoded list."
  (labels ((choose (one-hot-row category-list)
	     (cond ((null one-hot-row) nil)
		   ((= (car one-hot-row) 1) (car category-list))
		   (t (choose (cdr one-hot-row) (cdr category-list))))))
    (mapcar (lambda (l) (choose l categories))
	    one-hot-list)))

(multiple-value-bind (one-hot categories) (categorical-to-one-hot (replace-na *labels*))
  (one-hot-to-categories one-hot categories)) ;; works!

;; you have train data with labels (lol with labels)
;; these are split into a list of lols with labels
;; it wouldn't be bad to write them out into separate files
;; which one can read line by line to keep memory footprint low.

;; one general problem will be to normalize such data while processing line by line.

;; maybe it will be good to have the mathematical mean plus n
;; and using the mean plus n information and the new line one can generate a new mean plus new n
;; and is able to process the next line. By this, one can process and calculate a mean
;; just by processing one line a time row-wise.
;; At the end of the table, one has the means and the n's. And also the maxes and the mins for each column.
;; And also the standard deviation of each column. This information of the first run-through one could save somewhere.
;; And using this information, one can normalize row by row in a second run.

;; this might be the best strategy.



(defun extract-rows-by-index (lol index-list)
  "Return subset of lol selected by index-list for rows."
  (let ((index-list (sort (copy-seq index-list) #'<))) ;; sort has side effect on index-list! Thus copy!
    (loop for i from 0 to (length lol)
	  when (member i index-list)
	    collect (progn
		      (setf index-list (cdr index-list))
		      (elt lol i))))) ;; this works!
;; and it would work much faster, if using vectors or arrays!



(defparameter *data* (cl-csv:read-csv #P"/media/josephus/My Book/sc_atlas/data/atlas_data/atlas/e7_50_molecules_genes_transposed_r1000.tab" :separator #\Tab))
(length *data*)
    
(defparameter *t* (get-row *data* 10))

(defparameter *anno* (cl-csv:read-csv #P"/media/josephus/My Book/sc_atlas/data/atlas_data/atlas/anno_e7_50.tab" :separator #\Tab))

(defparameter *labels* (get-column-by-header-name *anno* "celltype" :one-correct-p t))

(defparameter *labels-replaced* (replace-na *labels*))

;; (multiple-value-bind (nums categories) (categorical-to-nums (replace-na *labels*))
;;   (defparameter *labels-nums* nums)
;;   (defparameter *categories* categories))

;; (multiple-value-bind (lols categories) (split-by-group (cdr *data*) *labels-replaced*)
;;   (defparameter *train-data* lols)
;;   (defparameter *categories* categories))

(multiple-value-bind (list-of-indexes categories) (split-by-group-indexes *labels-replaced*)
  (defparameter *train-indexes* list-of-indexes)
  (defparameter *categories* categories))

(defparameter *test* (extract-rows-by-index (cdr *data*) (elt *train-indexes* 0)))





;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; utilities
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun transpose (lol)
  (apply #'mapcar #'list lol))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; read-in data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun table-read (fpath)
  (cl-csv:read-csv fpath :separator #\Tab))

(defun molecules-read (fpath)
  (transpose (table-read fpath)))

(defun annotations-read (fpath)
  (let* ((anno (table-read fpath))
	 (y    (replace-na (get-column-by-header-name anno "celltype")))
	 (y-1h (categorical-to-one-hot y))
	 (y-nums (categorical-to-nums y)))
    (values y y-1h y-nums)))

(ql:quickload :alexandria)

        


;; ;; (macroexpand-1 (mutli-assign (*train-indexes* *categories*) (list-of-indexes categories) (split-by-group-indexes *labels-replaced*)))

;; ;; (multiple-value-bind (list-of-indexes categories) (split-by-group-indexes *labels-replaced*)
;; ;;   (defparameter *train-indexes* list-of-indexes)
;; ;;   (defparameter *categories* categories))

;; ;; (defmacro multi-assign (global-names local-names multi-value-expression)
;; ;;     `(multiple-value-bind (,@local-names) (,@multi-value-expression)
;; ;;        `(loop for g in (list ,@global-names)
;; ;; 	     for l in (list ,@local-names)
;; ;; 	      do `(defparameter ,g ,l))))

;; (defmacro multi-assign ((&rest globals) (&rest values))
;;   (destructuring-bind
;;       ((globals values))
;;       `((,globals) (,values))
;;   `(loop for g in ,globals
;; 	for v in ,values
;; 	do `(defparameter ,g ,v))))

;; (defmacro single-assign (var value)
;;   `(defparameter ,var ,value))

;; (macroexpand-1 (single-assign *b* 2))
;; (macroexpand-1 '(single-assign *c* 3))
;; (single-assign *c* 3) ;; works

;; (defmacro double-assign (vars values)
;;   `(loop for g in ,vars
;; 	for v in ,values
;; 	do `(defparameter ,g ,v)))
;; (macroexpand-1 '(double-assign (*a* *b* *c*) (5 6 7)))
;; ;; => (LOOP FOR G IN (*A* *B* *C*)
;; ;;          FOR V IN (5 6 7)
;; ;;          DO `(DEFPARAMETER ,G ,V)), T

;; ;; (double-assign (*a* *b* *c*) (5 6 7))
;; ;; gives error

;; (defmacro double-assign (vars values)
;;   `(loop for g in ',vars
;; 	for v in ',values
;; 	 do `(defparameter ,g ,v)))
;; (macroexpand-1 '(double-assign (*a* *b* *c*) (5 6 7)))
;; ;; => (LOOP FOR G IN '(*A* *B* *C*)
;; ;;          FOR V IN '(5 6 7)
;; ;;          DO `(DEFPARAMETER ,G ,V)), T

;; (double-assign (*a* *b* *c*) (5 6 7))

;; ;; (LOOP FOR G IN '(*A* *B* *C*)
;; ;;          FOR V IN '(5 6 7)
;; ;;       collect `(DEFPARAMETER ,G ,V))
;; ;; ;; ((DEFPARAMETER *A* 5) (DEFPARAMETER *B* 6) (DEFPARAMETER *C* 7))

;; (defmacro multi-assign (vars values)
;;   ``(progn ,@,`(loop for g in ',vars
;; 	for v in ',values
;; 	collect `(defparameter ,g ,v))))
;; (macroexpand-1 '(multi-assign (*a* *b* *c*) (5 6 7)))
;; (eval (multi-assign (*a* *b* *c*) (5 6 7))) ;; how can I get rid of eval?

;; (defmacro multi-assign (vars values)
;;   `(progn ,@`(loop for g in ',vars
;; 	for v in ',values
;; 		   collect `(defparameter ,g ,v)))) ;; wrong!

;; (defun multi-assign (vars values)
;;   (loop for g in vars
;; 	for v in values
;; 	collect `(defparameter ,g ,v)))
;; (multi-assign '(*a* *b* *c*) '(10 11 12))
;; ;; => ((DEFPARAMETER *A* 10) (DEFPARAMETER *B* 11) (DEFPARAMETER *C* 12))

;; ;; so if I make a macro out of it, then it gets the extra eval step ... (this was correct!)

;; (defmacro multi-assign (vars values)
;;   (cons 'progn (loop for g in vars
;; 		     for v in values
;; 		     collect `(defparameter ,g ,v))))
;; (multi-assign (*a* *b* *c*) (10 11 12))  ;; this worked!
;; (macroexpand-1 '(multi-assign (*a* *b* *c*) (10 11 12)))
;; ;; => (PROGN (DEFPARAMETER *A* 10) (DEFPARAMETER *B* 11) (DEFPARAMETER *C* 12)), T

(defmacro multi-assign (vars values)
  `(progn ,@(loop for g in vars
		  for v in values
		  collect `(defparameter ,g ,v))))
;; (macroexpand-1 '(multi-assign (*a* *b* *c*) (10 11 12)))
;; ;; => (PROGN (DEFPARAMETER *A* 10) (DEFPARAMETER *B* 11) (DEFPARAMETER *C* 12)), T
;; (multi-assign (*a* *b* *c*) (20 21 22)) ;; works!

;; so my conclusion is: get as much as possible done within a function!
;; think in conses and really pasting together the syntax-tree!


;; (macroexpand-1 (mutli-assign (*train-indexes* *categories*) (list-of-indexes categories) (split-by-group-indexes *labels-replaced*)))

;; (multiple-value-bind (list-of-indexes categories) (split-by-group-indexes *labels-replaced*)
;;   (defparameter *train-indexes* list-of-indexes)
;;   (defparameter *categories* categories))

;; so what I want to do is:

(macroexpand-1 (multi-bind (*train-indexes* *categories*) (list-of-indexes categories) (split-by-group-indexes *labels-replaced*)))
;; to expand to:
;; (multiple-value-bind (list-of-indexes categories) (split-by-group-indexes *labels-replaced*)
;;   (multi-assign (*train-indexes* *categories*) (list-of-indexes categories)))

(defmacro multi-bind (global-vars inner-vars multi-value-expression)
  `(multiple-value-bind ,inner-vars ,multi-value-expression
     (multi-assign ,global-vars ,inner-vars)))

(multi-bind (*a* *b* *c*) (a b c) (values 30 31 32)) ;; works! Just at first try!

;; now is only the question how to make inner-vars safe - gensym

(defmacro multi-bind (global-vars multi-value-expression)
  `(let ((inner-vars (loop for x in ,global-vars collect (gensym))))
     `(multiple-value-bind ,inner-vars ,multi-value-expression
	(multi-assign ,global-vars ,inner-vars))))
(macroexpand-1 '(multi-bind (*a* *b* *c*) (41 42 43)))
;; => (LET ((INNER-VARS
;;           (LOOP FOR X IN (*A* *B* *C*)
;;                 COLLECT (GENSYM))))
;;      `(MULTIPLE-VALUE-BIND ,INNER-VARS
;;           ,MULTI-VALUE-EXPRESSION
;;         (MULTI-ASSIGN ,GLOBAL-VARS ,INNER-VARS))), T
(multi-bind (*a* *b* *c*) (41 42 43))
;; one evaluation level too little, so I thought get rid of outer '`'
;; and I learned that '`,` is NOT equivalent with leaving it out!
;; because it places there the chunk of code without execution
;; and protects it when evaluated.
;; it is '<place it here unevaluated> or `(quote ,<place it here>)
;; and it NOT equivalent with <place it here>!!

;; generally, '`' says: the following is a template!
;; no matter how deeply nested the list is, one can replace things
;; and then evaluate the whole code as a program!

;; by putting let-clauses or other clauses go in front of
;; the code-template part, one can e.g. prepare gensyms
;; which one inside the '`' template part genuinely ','ify them,
;; so that the gensym can come to power!

;; one has to think in paper-formal level when replacing with '`' and ',' and ',@'.
;; objectify the code.

;; the '`,'-ed parts are litterally replaced.
;; and at the end this code will be evaluated!

;; and it is totally worth to break down parts of replacements into other macros,
;; which makes thinking a lot easier!

;; this is the correct version!!
(defmacro multi-bind (global-vars multi-value-expression)
  (let ((inner-vars (loop for x in `,global-vars collect (gensym))))
     `(multiple-value-bind ,inner-vars ,multi-value-expression
	(multi-assign ,global-vars ,inner-vars))))
(macroexpand-1 '(multi-bind (*a* *b* *c*) (values 41 42 43)))
;; (MULTIPLE-VALUE-BIND (#:G754 #:G755 #:G756)
;;        (VALUES 41 42 43)
;;   (MULTI-ASSIGN (*A* *B* *C*) (#:G754 #:G755 #:G756))), T
(multi-bind (*a* *b* *c*) (values 41 42 43)) ;; works!!!
