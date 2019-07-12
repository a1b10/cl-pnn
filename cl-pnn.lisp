;;;; cl-pnn.lisp

(defpackage :cl-pnn
  (:use :cl :alexandria :cl-csv :computable-reals)
  (:export :remove-empty-string))

(in-package #:cl-pnn)

;; (ql:quickload :alexandria)
;; (ql:quickload :cl-csv)
;; (use-package :cl-csv)

;; (ql:quickload :computable-reals)
;; (use-package :computable-reals)

(setq *PRINT-PREC* 50)
(defconstant +e+ (exp-r 1))

(defun head (l &optional (n 5))
  (subseq l 0 n))

(defun tail (l &optional (n 5))
  (let ((len (length l)))
    (subseq l (- len n) len)))

(defun index-extract (x index-list &key (exclude-nil-p t))
  (mapcar (lambda (i) (elt x i)) (if exclude-nil-p
				     (remove-if #'null index-list)
				     index-list)))

(defun unique (l &key (test #'string=))
  "Return unique list."
  (remove-duplicates l :from-end t :test test))

(defun dim (lol &optional (contains-row-names-p nil) (contains-header-p nil))
  "Return dimensions of a list of list or matrix."
  (let ((nrows (if contains-header-p
		   (1- (length lol))
		   (length lol)))
	(ncols (if contains-row-names-p
		   (1- (length (elt lol 0)))
		   (length (elt lol 0)))))
    (values nrows ncols contains-header-p contains-row-names-p)))


  
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

(defun remove-empty-string (string-list)
  (remove-if (lambda (s) (string= "" s))
	     string-list))

(ql:quickload :cl-ppcre)
(defun parse-whitespace-delimited-string (multiline-string)
  "Return whitespace delimited multiline string as lol whitespace-only lines ignoring."
  (mapcar (lambda (line)
	    (read-from-string (concatenate 'string "(" line ")")))
	  (remove-empty-string (cl-ppcre:split "\\n" multiline-string))))

(defun difference-square-sum (l1 l2)
  "Return sum of squared differences."
  (let ((res 0))
    (loop for x in l1
	  for y in l2
	  do (setf res (+r res (expt-r (-r x y) 2)))
	  finally (return res))))

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


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; read tab delimited files
;; https://quickref.common-lisp.net/cl-csv.html#Introduction
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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




(defun gaussian-probability-density-function (l1 l2 sigma)
  "Return gaussian-probability-density-function of probe l1 and row l2."
  (/r (expt-r +e+ (-r (difference-square-sum l1 l2))) (*r 2 (expt-r sigma 2))))

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


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; multi-bind macro
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmacro multi-bind (global-vars multi-value-expression)
  (let ((inner-vars (loop for _ in `,global-vars collect (gensym))))
     `(multiple-value-bind ,inner-vars ,multi-value-expression
	,@(loop for g in `,global-vars
		for v in `,inner-vars
		collect `(defparameter ,g ,v)))))

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
    (list (index-extract l (elt split-indexes 0))
	  (index-extract l (elt split-indexes 1)))))

(defun train-test-index-split-limited (l &key (test-size 0.3) (max-n 100) (random-state *random-state*))
  "Return l split indexes given ratio-test-size and max train number."
  (let* ((split-indexes (train-test-index-split l :test-size test-size :random-state random-state))
	 (train-indexes (elt split-indexes 0))
	 (test-indexes  (elt split-indexes 1))
	 (len           (length train-indexes))
	 (train-sel-indexes nil)
	 (train-dumped-indexes nil))
    (setf train-sel-indexes (if (> len max-n)
				(subseq train-indexes 0 max-n)
				train-indexes))
    (setf train-dumped (if (> len max-n)
			   (subseq train-indexes max-n len)
			   ()))
    (list train-sel-indexes train-dumped test-indexes)))

(defun train-test-split-limited (l &key (test-size 0.3) (max-n 100) (random-state *random-state*))
  "Return l split given ratio-test-size and max train number."
  (let ((split-indexes (train-test-index-split-limited l
						       :test-size test-size
						       :max-n max-n
						       :random-state random-state)))
    (list (index-extract l (elt split-indexes 0))
	  (index-extract l (elt split-indexes 1))
	  (index-extract l (elt split-indexes 2)))))

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

(defun stratified-train-test-split (lol labels &key (test-size 0.3) (random-state *random-state*))
  "Return data lol split in a stratified manner."
  (destructuring-bind (train-indexes test-indexes categories tsize)
      (stratified-train-test-index-split labels :test-size test-size :random-state random-state)
    (flet ((select (idxs) (mapcar (lambda (i) (elt lol i)) idxs)))
      (list (select train-indexes) (select test-indexes) train-indexes test-indexes categories tsize))))

(defun stratified-train-test-index-split-limited (labels &key (test-size 0.3) (max-n 100) (random-state *random-state*) (verbose t))
  "Return split indexes randomly selecting maximally max-n sample indexes 
   (a much more memory saving method); limit for calculation saving."
  (multiple-value-bind (split-indexes categories) (split-by-group-indexes labels)
    (let* ((nshuffled-stratified-indexes (mapcar (lambda (x) (train-test-split-limited (nshuffle x :random-state random-state)
										       :test-size test-size
										       :max-n max-n
										       :random-state random-state))
						 split-indexes))
	   (train-indexes (alexandria:flatten (mapcar #'first nshuffled-stratified-indexes)))
	   (dumped-indexes  (alexandria:flatten (mapcar #'second nshuffled-stratified-indexes)))
	   (test-indexes (alexandria:flatten (mapcar #'third nshuffled-stratified-indexes))))
      (when verbose
	(format t "Category counts:~%")
	(loop for si in split-indexes
	      for c in categories
	      do (format t "~%~A~A: ~A" c #\Tab (length si)))
	(format t "~%~%Only train no test (too few members):~%")
	(loop for si in split-indexes
	      for c in categories
	      when (< (* (length si) test-size) 1)
		do (format t "~%~A~A: ~A" c #\Tab (length si)))
	(format t "~%~%Limits:~%")
	(loop for nsi in nshuffled-stratified-indexes
	      for c in categories
	      when (not (zerop (length (cadr nsi))))
		do (format t "~%Dumped ~A samples for ~A" (length (cadr nsi)) c)))
      (list train-indexes dumped-indexes test-indexes categories test-size))))
