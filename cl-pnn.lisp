;;;; cl-pnn.lisp

(in-package #:cl-pnn)

(ql:quickload :alexandria)

(ql:quickload :computable-reals)
(use-package :computable-reals)
(setq *PRINT-PREC* 50)
(defconstant +e+ (exp-r 1))

(defun unique (l &key (test #'string=))
  "Return unique list."
  (remove-duplicates l :from-end t :test test))

(defun difference-square-sum (l1 l2)
  "Return sum of squared differences."
  (let ((res 0))
    (loop for x in l1
	  for y in l2
	  do (setf res (+r res (expt-r (-r x y) 2)))
	  finally (return res))))
  
(defun gaussian-probability-density-function (l1 l2 sigma)
  "Return gaussian-probability-density-function of probe l1 and row l2."
  (/r (expt-r +e+ (-r (difference-square-sum l1 l2))) (*r 2 (expt-r sigma 2))))

(defun dim (lol &optional (contains-row-names-p nil) (contains-header-p nil))
  "Return dimensions of a list of list or matrix."
  (values (if contains-header-p
              (1- (length lol))
	      (length lol))
	  (if contains-row-names-p
              (1- (length (elt lol 0)))
	      (length (elt lol 0)))
	  contains-header-p
	  contains-row-names-p))

(defun gaussian-probability-density-function-matrix (l1 x-train sigma &key (contains-row-names-p nil) (contains-header-p nil))
  "Return probability value for a sample and matrix of train values."
  (multiple-value-bind (m p) (dim x-train contains-row-names-p contains-header-p)
    (let ((expression (*r m (/r 1 (*r (expt-r (*r 2 +pi-r+) (/r p 2)) (expt-r sigma p))))))
      (/r 1 (*r expression
		(mapcar #'(lambda (l2) (gaussian-probability-density-function l1 l2 sigma))
			x-train))))))

;; new approach

(defun sum-r (l)
  (reduce #'+r l))

(defun gaussian-exponent (l1 l2 sigma)
  (expt-r +e+ (/r (difference-square-sum l1 l2) (*r 2 (expt-r sigma 2)))))

(defun gaussian-density (l1 x-train sigma)
  "density function using gaussian distribution after ch1_3-ch1_3.pdf"
  (multiple-value-bind (m p) (dim x-train)
    (let* ((coefficient (/r 1 (*r (expt-r (*r 2 +pi-r+) (/r p 2)) m))))
      (*r coefficient (sum-r (mapcar (lambda (l2) (gaussian-exponent l1 l2 sigma))
				     x-train))))))

(defparameter *l1* (elt *test-matrix* 0))
(defparameter *x-train* *train-matrix*)
(defparameter *sigma* 1)
(multi-bind (*m* *p*) (dim *x-train* nil nil))
(defparameter *expr* (*r *m* (/r 1 (*r (expt-r (*r 2 +pi-r+) (/r *p* 2)) (expt-r *sigma* *p*)))))
(defparameter *t2* (gaussian-density *l1* *x-train* *sigma*)) ;; takes quite long!

#|
The question would be by which strategy to reduce (prune) the train-data set
without any information loss.

Here not only feature reduction but also train-data reduction is needed!

Actually the point is that one can randomly select - a minimal core of train datasets.
They can be taken as seed and allt he remaining data will be added only if they improve the outcome.
However, this is very validation dataset-dependent.
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
    (values (mapcar (lambda (x) (nreverse x)) result)
	    levels)))

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

;; PNN brain image 100-80-73
;; http://www.rroij.com/open-access/probabilistic-neural-network-forbrain-tumor-classification.php?aid=41337

;; pnn mathworks https://de.mathworks.com/help/deeplearning/ug/probabilistic-neural-networks.html


;; keras ensembl with softmax https://towardsdatascience.com/the-softmax-function-neural-net-outputs-as-probabilities-and-ensemble-classifiers-9bd94d75932


;; PNN with complex exponential activation function image recognition
;; https://arxiv.org/pdf/1708.02733.pdf


;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; label conversions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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

(defun one-hot-to-categories (one-hot-list categories)
  "Return categories from one-hot-encoded list."
  (labels ((choose (one-hot-row category-list)
	     (cond ((null one-hot-row) nil)
		   ((= (car one-hot-row) 1) (car category-list))
		   (t (choose (cdr one-hot-row) (cdr category-list))))))
    (mapcar (lambda (l) (choose l categories))
	    one-hot-list)))

(defun arg-max (one-hot-list)
  "Return positions of 1 in one-hot-list."
  (mapcar (lambda (l) (position 1 l)) one-hot-list))

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



(defun index-extract (lol index-list)
  "Return subset of lol or l selected by index-list for rows."
  (let ((index-list (sort (copy-seq index-list) #'<))) ;; sort has side effect on index-list! Thus copy!
    (loop for i from 0 to (length lol)
	  when (member i index-list)
	    collect (progn
		      (setf index-list (cdr index-list))
		      (elt lol i))))) ;; this works!
;; and it would work much faster, if using vectors or arrays!


#|

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

(defparameter *test* (index-extract (cdr *data*) (elt *train-indexes* 0))) ;; yes it extracts rows by a given index!

|#



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


;; well, the working version of multi-assign had a 'progn' which is implicit inside multiple-value-bind.
;; let's get rid of it!
;; (defmacro multi-assign (vars values)
;;   `(progn ,@(loop for g in vars
;; 		  for v in values
;; 		  collect `(defparameter ,g ,v))))

;; let's think very flat

(defmacro multi-bind (global-vars multi-value-expression)
  (let ((inner-vars (loop for _ in `,global-vars collect (gensym))))
     `(multiple-value-bind ,inner-vars ,multi-value-expression
	,@(loop for g in `,global-vars
		for v in `,inner-vars
		collect `(defparameter ,g ,v)))))

(macroexpand-1 '(multi-bind (*a* *b* *c*) (values 41 42 43)))
;; => (MULTIPLE-VALUE-BIND (#:G761 #:G762 #:G763)
;;        (VALUES 41 42 43)
;;      (DEFPARAMETER *A* #:G761)
;;      (DEFPARAMETER *B* #:G762)
;;      (DEFPARAMETER *C* #:G763)), T
(macroexpand-1 '(multi-bind (*a* *b* *c* *d* *e*) (values 41 42 43 44 45)))
;; => (MULTIPLE-VALUE-BIND (#:G764 #:G765 #:G766 #:G767 #:G768)
;;        (VALUES 41 42 43 44 45)
;;      (DEFPARAMETER *A* #:G764)
;;      (DEFPARAMETER *B* #:G765)
;;      (DEFPARAMETER *C* #:G766)
;;      (DEFPARAMETER *D* #:G767)
;;      (DEFPARAMETER *E* #:G768)), T   ;; yes, that looks very correct!
(multi-bind (*a* *b* *c* *d* *e*) (values 41 42 43 44 45)) ;; works!





;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; improved train test split
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defparameter *random-state* (make-random-state t))

(defun integer-like-p (n)
  "Return t if modulo rest is zero else nil."
  (multiple-value-bind (divider modulo) (floor n)
    (zerop modulo)))

(defun divisor-zero-p (n)
  "Return t if divider is zero othervise nil."
  (zerop (floor n)))

(defun nshuffle (sequence &key (random-state *random-state*))
  (loop for i from (length sequence) downto 2
        do (rotatef (elt sequence (random i *random-state*))
                    (elt sequence (1- i))))
  sequence) ;; rosettac code 

(defun split (list count)
  "split list by position count (count is length of first list)."
  (let ((count (round count)))
    (values (subseq list 0 count) (subseq list count))))

(defun split-list (list count)
  "Return lol with inner list being split by psoition count."
  (multiple-value-list (split list count)))

(defun train-test-index-split (l &key (test-size 0.3) (random-state *random-state*))
  "Return indexes split given list by ratio test-size."
  (let* ((length (length l))
         (ratio-full-number (* test-size length))
	 (indexes (loop for i from 0 below length
			collect i)))
    (if (or (> 1 ratio-full-number) (not (integer-like-p ratio-full-number)))
      (setf ratio-full-number (floor ratio-full-number))
      (format t "~%Not even divider. Using next lower number:~A" ratio-full-number))
    (let ((indexes (nshuffle indexes :random-state random-state))
	  (train-number (- length ratio-full-number)))
      (split-list indexes train-number)))) ;; works

(defun train-test-split (l &key (test-size 0.3) (random-state *random-state*))
  "Return lol split given list of list by ratio test-size."
  (let ((split-indexes (train-test-index-split l :test-size test-size :random-state random-state)))
    (flet ((select (idxs) (mapcar (lambda (i) (elt l i)) idxs)))
      (list (select (elt split-indexes 0)) (select (elt split-indexes 1))))))

;; (defparameter *l* '("a" "b" "c" "d" "e" "f" "g" "h" "i" "j" "k" "l" "m"))
;; (length *l*)
;; (* 0.3 13) ;; 3.9
;; (loop for i from 0 to 13 collect i)
;; (floor 3.9) ;; 3
;; (train-test-index-split '(a b c d e f g h i j k l m) :test-size 0.3)
;; (train-test-split *l* :test-size 0.3)

;; (defparameter *length* (length *lr*))
;; (defparameter *rfn* (* 0.3 *length*))
;; (defparameter *indexes* (loop for i from 0 below *length* collect i))
;; (> 1 *rfn*) ;; nil
;; (integer-like-p *rfn*) ;; T -> null
;; (defparameter *indexes* (nshuffle *indexes*))
;; (defparameter *train-number* (- *length* *rfn*))
;; (split-list *indexes* *train-number*)

(defun head (l &optional (n 5))
  (subseq l 0 n))

(defun tail (l &optional (n 5))
  (let ((len (length l)))
    (subseq l (- len n) len)))

(defun stratified-train-test-index-split (labels &key (test-size 0.3) (random-state *random-state*) (verbose t))
  "Return split indexes (a much more memory saving method)."
  (multiple-value-bind (split-indexes categories) (split-by-group-indexes labels)
    (let* ((nshuffled-stratified-indexes (mapcar (lambda (x) (train-test-split
							      (nshuffle x :random-state random-state)
							      :test-size test-size
							      :random-state random-state))
						 split-indexes))
	   (train-indexes (alexandria:flatten (mapcar #'first nshuffled-stratified-indexes)))
	   (test-indexes  (alexandria:flatten (mapcar #'second nshuffled-stratified-indexes))))
      
      (when verbose
	(format t "Category counts:~%")
	(loop for si in split-indexes
	      for c in categories
	      do (format t "~%~A~A: ~A" c #\Tab (length si)))
	(format t "~%~%Only train no test (too few members):~%")
	(loop for si in split-indexes
	      for c in categories
	      when (< (* (length si) test-size) 1)
		do (format t "~%~A~A: ~A" c #\Tab (length si))))
      (list train-indexes test-indexes categories test-size)))) ;; it works!!

;; (multi-bind (*t1* *cat*) (split-by-group-indexes *labels-replaced*))
;; (defparameter *t2* (mapcar (lambda (x) (train-test-split (nshuffle x) :test-size 0.3)) *t1*))

;; (stratified-train-test-index-split '(a a a b b b b b b c c c d d d d d d d d d) :test-size 0.4)
      
;; (defparameter *ls* '(a a a b b b b b b c c c d d d d d d d d d))
;; (defparameter *ls* (mapcar #'string *ls*))
;; (defparameter *ixs* (split-by-group-indexes *ls*))
#| 
(let ((split-indexes *ixs*)
      (random-state *random-state*)
      (test-size 0.34))
(mapcar (lambda (x) (train-test-split
                       (nshuffle x :random-state random-state)
                       :test-size test-size
                       :random-state random-state))
        split-indexes))
|#
;; (stratified-train-test-index-split *ls* :test-size 0.3)

;; (defparameter *lr* (head *labels-replaced* 1000))
;; (stratified-train-test-index-split *lr* :test-size 0.3)
;; (untrace stratified-train-test-index-split)


(defun stratified-train-test-split (lol labels &key (test-size 0.3) (random-state *random-state*))
  "Return data lol split in a stratified manner."
  (destructuring-bind (train-indexes test-indexes categories tsize)
      (stratified-train-test-index-split labels :test-size test-size :random-state random-state)
    (flet ((select (idxs) (mapcar (lambda (i) (elt lol i)) idxs)))
      (list (select train-indexes) (select test-indexes) train-indexes test-indexes categories tsize))))

;; (stratified-train-test-split *ls* *ls* :test-size 0.4) ;; works!!

;; arrays and reshape
;; https://lispcookbook.github.io/cl-cookbook/arrays.html
;; (ql:quickload :numcl)
;; (in-package :numcl)
;; https://github.com/numcl/numcl
;; this is really good for reshaping around!



;; python text classification
;; https://towardsdatascience.com/text-classification-in-python-dd95d264c802

(defun string-to-number (s)
  "Return number from string."
  (with-input-from-string (in s)
    (read in)))

(defun numbers-extract (train-test-data)
  "Return numeric lol from train-data or test-data."
  (mapcar (lambda (row)
	    (mapcar #'string-to-number
		    (cdr row)))
	  train-test-data))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; test the big data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defparameter *data* (cl-csv:read-csv #P"/media/josephus/My Book/sc_atlas/data/atlas_data/atlas/e7_50_molecules_genes_transposed_r1000.tab" :separator #\Tab))
(length *data*)
    
(defparameter *t* (get-row *data* 10))

(defparameter *anno* (cl-csv:read-csv #P"/media/josephus/My Book/sc_atlas/data/atlas_data/atlas/anno_e7_50.tab" :separator #\Tab))

(defparameter *labels* (get-column-by-header-name *anno* "celltype" :one-correct-p t))

(defparameter *labels-replaced* (replace-na *labels*))

(defparameter *result* (split-by-group-indexes *labels-replaced*))


;; split data in a stratified manner

(defparameter *split-result* (stratified-train-test-index-split *labels-replaced* :test-size 0.1))
(defparameter *train-indexes* (elt *split-result* 0))
(defparameter *test-indexes* (elt *split-result* 1))
(defparameter *categories* (elt *split-result* 2))
(defparameter *split* (elt *split-result* 3))
(length *split-result*)

;; extract train test data

(defparameter *train-data* (index-extract (cdr *data*) *train-indexes*))
(defparameter *test-data*  (index-extract (cdr *data*) *test-indexes*))
(defparameter *train-labels* (index-extract *labels-replaced* *train-indexes*))
(defparameter *test-labels* (index-extract *labels-replaced* *train-indexes*))

;; extract row names
;; convert to numbers
(defparameter *train-row-names* (mapcar #'first *train-data*))
(defparameter *test-row-names*  (mapcar #'first *test-data*))

(defparameter *train-matrix* (numbers-extract *train-data*))
(defparameter *test-matrix*  (numbers-extract *test-data*))

(difference-square-sum (elt *train-matrix* 0) (elt *test-matrix* 0))
(gaussian-probability-density-function (elt *train-matrix* 0) (elt *test-matrix* 0) 1)
(defparameter *t1* (gaussian-probability-density-function-matrix (elt *test-matrix* 0) *train-matrix* 1))
