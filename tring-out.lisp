;; (defun index-extract (lol index-list)
;;   "Return subset of lol or l selected by index-list for rows."
;;   (let ((index-list (sort (copy-seq index-list) #'<))) ;; sort has side effect on index-list! Thus copy!
;;     (loop for i from 0 to (length lol)
;; 	  when (member i index-list)
;; 	    collect (progn
;; 		      (setf index-list (cdr index-list))
;; 		      (elt lol i))))) ;; this works!
;; and it would work much faster, if using vectors or arrays!


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; test new approach
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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

The other possibility is to cluster the train data and see which data are similar
and to remove same data points (or highly similar data points).

Ah! or similar to xgboost or random-forest, one should build 
small classifiers who are used in ensemble!!
|#


#|
(defparameter *mydata*
	   "
row	length 	area 	label
0 	0.5 	0.7 	BLUE
1 	0.2 	0.5 	BLUE
2 	0.8 	0.8 	RED
3 	0.4 	0.5 	RED
4 	0.8 	0.5 	GREEN
5 	0.6 	0.3 	GREEN
6 	0.3 	0.2 	GREEN
")

(defparameter *mydata* (parse-whitespace-delimited-string *mydata*))

((ROW LENGTH AREA LABEL) (0 0.5 0.7 BLUE) (1 0.2 0.5 BLUE) (2 0.8 0.8 RED)
 (3 0.4 0.5 RED) (4 0.8 0.5 GREEN) (5 0.6 0.3 GREEN) (6 0.3 0.2 GREEN))
|#


;; PNN brain image 100-80-73
;; http://www.rroij.com/open-access/probabilistic-neural-network-forbrain-tumor-classification.php?aid=41337

;; pnn mathworks https://de.mathworks.com/help/deeplearning/ug/probabilistic-neural-networks.html

;; keras ensembl with softmax https://towardsdatascience.com/the-softmax-function-neural-net-outputs-as-probabilities-and-ensemble-classifiers-9bd94d75932

;; PNN with complex exponential activation function image recognition
;; https://arxiv.org/pdf/1708.02733.pdf



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


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; multi-bind macro making of
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
         
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

  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; train-test split tests
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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

;; (train-test-split-limited *l* :test-size 0.3 :max-n 5)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; stratified train-test
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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

;; (stratified-train-test-split *ls* *ls* :test-size 0.4) ;; works!!

;; arrays and reshape
;; https://lispcookbook.github.io/cl-cookbook/arrays.html
;; (ql:quickload :numcl)
;; (in-package :numcl)
;; https://github.com/numcl/numcl
;; this is really good for reshaping around!

;; python text classification
;; https://towardsdatascience.com/text-classification-in-python-dd95d264c802

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; (defparameter *t3* '("a" "a" "a" "b" "b" "b" "b" "c" "c"))
;; (stratified-train-test-index-split-limited '("a" "a" "a" "b" "b" "b" "b" "c" "c") :test-size (/ 1 3) :max-n 2)
;; (multi-bind (*si* *c*) (split-by-group-indexes *t3*))
;; (defparameter *test-size* (/ 1 3))
;; (defparameter *nsi* (mapcar (lambda (x) (train-test-split-limited (nshuffle x :random-state *random-state*)
;; 								  :test-size *test-size*
;; 								  :max-n 2
;; 								  :random-state *random-state*))
;; 			    *si*))
;; (defparameter *tis* (alexandria:flatten (mapcar #'first *nsi*)))
;; (defparameter *dis* (alexandria:flatten (mapcar #'second *nsi*)))
;; (defparameter *sis* (alexandria:flatten (mapcar #'third *nsi*)))



(multiple-value-bind (one-hot categories) (categorical-to-one-hot (replace-na *labels*))
  (one-hot-to-categories one-hot categories)) ;; works!
