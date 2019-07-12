;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; test the big data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(load #P"cl-pnn.lisp")
(in-package :cl-pnn)

(defparameter *data* (cl-csv:read-csv #P"/media/josephus/My Book/sc_atlas/data/atlas_data/atlas/e7_50_molecules_genes_transposed_r1000.tab" :separator #\Tab))
(defparameter *data* cl-user::*data*)
(makunbound cl-user::*data*)
(length *data*)
    
(defparameter *t* (get-row *data* 10))

(defparameter *sigma* 1)

(defparameter *anno* (cl-csv:read-csv #P"/media/josephus/My Book/sc_atlas/data/atlas_data/atlas/anno_e7_50.tab" :separator #\Tab))

(defparameter *labels* (get-column-by-header-name *anno* "celltype"))

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


;; split and limit data in stratified manner
(defparameter *split-limited-result* (stratified-train-test-index-split-limited *labels-replaced* :test-size 0.3 :max-n 50))
(defparameter *train-lim-indexes* (elt *split-limited-result* 0))
(defparameter *test-lim-indexes* (elt *split-limited-result* 1))
(defparameter *categories-lim* (elt *split-limited-result* 2))
(defparameter *split-lim* (elt *split-limited-result* 3))

(defparameter *train-lim-data* (index-extract (cdr *data*) *train-lim-indexes*))
(defparameter *test-lim-data*  (index-extract (cdr *data*) *test-lim-indexes*))
(defparameter *train-lim-labels* (index-extract *labels-replaced* *train-lim-indexes*))
(defparameter *test-lim-labels* (index-extract *labels-replaced* *train-lim-indexes*))

;; extract row names
;; convert to numbers
(defparameter *train-lim-row-names* (mapcar #'first *train-lim-data*))
(defparameter *test-lim-row-names*  (mapcar #'first *test-lim-data*))

(defparameter *train-lim-matrix* (numbers-extract *train-lim-data*))
(defparameter *test-lim-matrix*  (numbers-extract *test-lim-data*))

(difference-square-sum (elt *train-lim-matrix* 0) (elt *test-lim-matrix* 0))
(gaussian-probability-density-function (elt *train-lim-matrix* 0) (elt *test-lim-matrix* 0) 1)
(defparameter *t4* (gaussian-density (elt *test-lim-matrix* 0) (head *train-lim-matrix* 9) 1))
(defparameter *t5* (mapcar (lambda (l2) (gaussian-exponent (elt *test-lim-matrix* 0) l2 *sigma*))
				     (head *train-lim-matrix* 9)))
(defparameter *t6* (sum-r *t5*))
(defparameter *t1* (gaussian-density (elt *test-lim-matrix* 0) (head *train-lim-matrix* 9) 1))
;; this is strange (it gives a huuuuge number)... what is correct at the end?  


#|

#################
# kernel.R
#################

# D square
# @param Xa A vector describing a new observation
# @param X The training set of observation
ds <- function(Xa, X) {
    value <- (X - Xa) %*% t(X - Xa)
    return(as.numeric(value))
}

Xa <- c(1, 2, 3, 4, 5, 6)
X <- c(10, 10, 10, 10, 10, 10)


;; this is actually 
;;;; (defun as-creal (x) (*r 1 x))

(defun ds (x-train-vec x-observed-vec)
  (let ((x-train-vec (mapcar #'as-creal x-train-vec))
        (x-observed-vec (mapcar #'as-creal x-observed-vec)))
    (sum-r (mapcar #'*r x-train-vec x-observed-vec))))



# Exponential kernel
# @param Xa A vector describing a new observation
# @param X The training set of observation
# @param sigma The smooth parameter
pattern <- function(Xa, X, sigma) {
    res <- exp( - ds(Xa, X) / (2 * sigma ^ 2) )
    return(as.numeric(res))
}

(defun exponent-expression (x-train-vec x-observed-vec sigma)
  (exp-r (/r (-r (ds x-train-vec x-observed-vec)) (*r 2 sigma sigma))))

;; (exp-r (/r (-r (ds (elt *tlm* 0) *xov*)) (*r 2 1 1)))

# Apply kernel over all patterns from category A
# @param Xa A vector describing a new observation
# @param X The training set of observation
# @param sigma The smooth parameter
patterns <- function(Xa, X, sigma)
    apply(Xa, 1, pattern, X, sigma)

(defun exponent-over-train-matrix (x-train-matrix x-observed-vec sigma)
  (mapcar (lambda (x-vec) (exponent-expression x-vec x-observed-vec sigma))
          x-train-matrix))

;; (mapcar (lambda (x-vec) (head x-vec)) *tlm*)
;; (exponent-expression (elt *tlm* 0) *xov* 1)

# Sum the results of applying the kernel over all patterns
# @param X Pattern from which we have to decide a category. It is a set of measurements represented by a p-dimensional vector
# @param Xa One of the training patterns from category A
# @param sigma Smoothing parameter
fA <- function(Xa, X, sigma) {
    if(missing(Xa)) stop("Xa is missing")
    if(missing(X)) stop("X is missing")
    if(missing(sigma)) stop("sigma is missing")
    p <- length(X) # Dimensionality of measurement space # seems wrong to me!
    m <- length(Xa[,1]) # Total number of training patterns from category A
    f <- 1 /((2 * pi) ^ (p / 2) * sigma ^ p) / m * sum(patterns(Xa, X, sigma)) # Probability density function
    return(f)
}

(defun gaussian-density (x-train-matrix x-observed-vec sigma)
  (multiple-value-bind (m p) (dim x-train-matrix)
    (let ((coeff (/r 1 (*r (expt-r (*r 2 +pi-r+) (/r p 2)) m))))
      (sum-r (mapcar (lambda (exponent) (*r coeff exponent))
                     (exponent-over-train-matrix x-train-matrix x-observed-vec sigma))))))

(defparameter *tlm* (head *train-lim-matrix* 9))
(defparameter *xov* (elt *test-lim-matrix* 0))
(gaussian-density *tlm* *xov* 1)                   ;; this was tested

(multiple-value-bind (m p) (dim *tlm*)
  (defparameter *m* m)
  (defparameter *p* p))

(exponent-over-train-matrix *tlm* *xov* 1) ;; error

53 x 30 = 1590 => 1600 => (setq *PRINT-PREC* 1600)


#################
# guess-probabilities.R
#################

# Predict the probabilities of each category given X
# @param nn An already trained probabilistic neural network
# @param X Pattern from which we have to decide a category. It is a set of measurements represented by a p-dimensional vector
guess.probabilities.of.each.category <- function(nn, X) {
    results <- vector()
    for(category in nn$categories) {
        Xa <- nn$set[nn$set[,nn$category.column] == category,]
        Xa <- as.matrix(Xa[,-nn$category.column])
        results <- c(results, fA(Xa, X, nn$sigma))
    }
    probs <- results / sum(results)
    names(probs) <- nn$categories
    return(probs)
}

(defun creal-to-double-float (creal-num)
  "Return double float from creal representation."
  (let* ((s (with-output-to-string (out)
	      (print-r creal-num 20 :stream out)))
	 (s (string-trim "." s)))
    (let ((*read-default-float-format* 'double-float))
      (with-input-from-string (in s)
	(read in)))))

(defun predict (x-train-matrix y-train x-observed-vec sigma)
  (multiple-value-bind (split-x-trains categories) (split-by-group x-train-matrix y-train)
    (let* ((results (mapcar (lambda (x-train-mat) 
                              (gaussian-density x-train-mat x-observed-vec sigma))
                            split-x-trains))
           (total   (sum-r results))
           (probabilities (mapcar (lambda (x) (/r x total))
                                  results))
           (probabilities (mapcar #'creal-to-double-float probabilities)) ;; for comparison
           (prob-cat-pairs (zip probabilities categories))
           (sorted-prob-cat-pairs (sort (copy-seq prob-cat-pairs) #'> :key #'car))
           (prediction (car sorted-prob-cat-pairs)))
      (values prediction
              sorted-prob-cat-pairs))))


(predict *train-lim-matrix* *train-lim-labels* (elt *test-lim-matrix* 0) 1)


(multiple-value-bind (split-x-trains categories)
                     (split-by-group *train-lim-matrix* *train-lim-labels*)
  (mapcar (lambda (x-train-mat) (gaussian-density x-train-mat (elt *test-lim-matrix* 0) 1))
          split-x-trains))

(setq *PRINT-PREC* 3000)

(multiple-value-bind (split-x-trains categories)
                     (split-by-group *train-lim-matrix* *train-lim-labels*)
  (mapcar (lambda (x-train-mat) (gaussian-density x-train-mat (elt *test-lim-matrix* 0) 1))
          split-x-trains))

;; with 3000 positions still many are below radar

;; would the numbers be smaller, then the numbers would be better


;; the sum of TPMs in different sampels always add up to the same number
;; denominator required to calculate proportions is the same
;; regardless of the sample you are looking at.

(defparameter *train-lim-matrix-tpm* (tpm-pseudo-r *train-lim-matrix*)) 

(defparameter *cols* (transpose *train-lim-matrix*))
(defparameter *col-totals* (mapcar (lambda (col) (reduce #'+r col)) *cols*))
(defparameter *cols-normed* (mapcar (lambda (col total)
                                     (mapcar (lambda (x) (if (and (zerop x)
                                                                  (zerop total))
                                                             0
                                                             (/r x total))) 
                                             col))
                                    *cols* *col-totals*))
(defparameter *rows* (transpose *cols-normed*))
(defparameter *cols-rows-normed* (mapcar #'norm-list-r *rows*)) ;; works!
(defparameter *train-lim-matrix-tpm* (tpm-pseudo-r *train-lim-matrix*))


(setq *PRINT-PREC* 500) ;; actually 500 arleady good for this
(defparameter *tlm* (head *train-lim-matrix-tpm* 9))
(defparameter *xov* (elt *test-lim-matrix* 0))
(gaussian-density *tlm* *xov* 1)                   ;; this was tested


(setq *PRINT-PREC* 1000)
(multi-bind (*train-lim-matrix-tpm* *col-totals*) (tpm-pseudo-r *train-lim-matrix*))
(defparameter *test-lim-matrix-tpm* (tpm-pseudo-r *test-lim-matrix* *col-totals*))
(predict *train-lim-matrix-tpm* *train-lim-labels* (elt *test-lim-matrix* 0) 1)
;; zero division error

;; one has to bring rigorous rules
;; like any division by zero -> zero!

;; also a special treatment of the sparse data!

(defun any-in (x l &key (test #'eql))
  (cond ((null l) nil)
        ((funcall test x (car l)) l)
        (t (any-in x (cdr l) :test test))))
  

(defun //r (x &rest rest)
  (if (any-in 0 rest)
      0
      (/r )))
  



################
# guess.r
################

#' Guess
#' 
#' Infers the category of a new observation.
#' 
#' Given an already trained and smoothed Probabilistic neural network, the function \code{guess} gives the category with the highest probability, and the probabilities of each category.
#' 
#' @param nn A trained and smoothed Probabilistic neural network.
#' @param X A vector describing a new observation.
#' 
#' @seealso \code{\link{pnn-package}}, \code{\link{learn}}, \code{\link{smooth}}, \code{\link{perf}}, \code{\link{norms}}
#' 
#' @return A \code{list} of the guessed category and the probabilities of each category.
#' 
#' @examples
#' library(pnn)
#' data(norms)
#' pnn <- learn(norms)
#' pnn <- smooth(pnn, sigma=0.8)
#' guess(pnn, c(1,1))
#' guess(pnn, c(1,1))$category
#' guess(pnn, c(1,1))$probabilities
#' guess(pnn, c(2,1))
#' guess(pnn, c(1.5,1))
#' @export
guess <- function(nn, X) {
    X <- matrix(X, ncol=nn$k)
    probs <- guess.probabilities.of.each.category(nn, X)
    if(is.na(probs[1])) return(NA)
    category <- names(probs[probs == max(probs)])
    results <- list(category=category, probabilities=probs)
    return(results)
}


|#









;; at least test with small data

#|

(defparameter *dstring* 
"
row     sepal_length  sepal_width  petal_length  petal_width
0             5.1          3.5           1.4          0.2
1             4.9          3.0           1.4          0.2
2             4.7          3.2           1.3          0.2
3             4.6          3.1           1.5          0.2
4             5.0          3.6           1.4          0.2
5             5.4          3.9           1.7          0.4
6             4.6          3.4           1.4          0.3
7             5.0          3.4           1.5          0.2
8             4.4          2.9           1.4          0.2
9             4.9          3.1           1.5          0.1
10            5.4          3.7           1.5          0.2
11            4.8          3.4           1.6          0.2
12            4.8          3.0           1.4          0.1
13            4.3          3.0           1.1          0.1
14            5.8          4.0           1.2          0.2
15            5.7          4.4           1.5          0.4
16            5.4          3.9           1.3          0.4
17            5.1          3.5           1.4          0.3
18            5.7          3.8           1.7          0.3
19            5.1          3.8           1.5          0.3
20            5.4          3.4           1.7          0.2
21            5.1          3.7           1.5          0.4
22            4.6          3.6           1.0          0.2
23            5.1          3.3           1.7          0.5
24            4.8          3.4           1.9          0.2
25            5.0          3.0           1.6          0.2
26            5.0          3.4           1.6          0.4
27            5.2          3.5           1.5          0.2
28            5.2          3.4           1.4          0.2
29            4.7          3.2           1.6          0.2
30            4.8          3.1           1.6          0.2
31            5.4          3.4           1.5          0.4
32            5.2          4.1           1.5          0.1
33            5.5          4.2           1.4          0.2
34            4.9          3.1           1.5          0.2
35            5.0          3.2           1.2          0.2
36            5.5          3.5           1.3          0.2
37            4.9          3.6           1.4          0.1
38            4.4          3.0           1.3          0.2
39            5.1          3.4           1.5          0.2
40            5.0          3.5           1.3          0.3
41            4.5          2.3           1.3          0.3
42            4.4          3.2           1.3          0.2
43            5.0          3.5           1.6          0.6
44            5.1          3.8           1.9          0.4
45            4.8          3.0           1.4          0.3
46            5.1          3.8           1.6          0.2
47            4.6          3.2           1.4          0.2
48            5.3          3.7           1.5          0.2
49            5.0          3.3           1.4          0.2
50            7.0          3.2           4.7          1.4
51            6.4          3.2           4.5          1.5
52            6.9          3.1           4.9          1.5
53            5.5          2.3           4.0          1.3
54            6.5          2.8           4.6          1.5
55            5.7          2.8           4.5          1.3
56            6.3          3.3           4.7          1.6
57            4.9          2.4           3.3          1.0
58            6.6          2.9           4.6          1.3
59            5.2          2.7           3.9          1.4
60           5.0          2.0           3.5          1.0
61           5.9          3.0           4.2          1.5
62           6.0          2.2           4.0          1.0
63           6.1          2.9           4.7          1.4
64           5.6          2.9           3.6          1.3
65           6.7          3.1           4.4          1.4
66           5.6          3.0           4.5          1.5
67           5.8          2.7           4.1          1.0
68           6.2          2.2           4.5          1.5
69           5.6          2.5           3.9          1.1
70           5.9          3.2           4.8          1.8
71           6.1          2.8           4.0          1.3
72           6.3          2.5           4.9          1.5
73           6.1          2.8           4.7          1.2
74           6.4          2.9           4.3          1.3
75           6.6          3.0           4.4          1.4
76           6.8          2.8           4.8          1.4
77           6.7          3.0           5.0          1.7
78           6.0          2.9           4.5          1.5
79           5.7          2.6           3.5          1.0
80           5.5          2.4           3.8          1.1
81           5.5          2.4           3.7          1.0
82           5.8          2.7           3.9          1.2
83           6.0          2.7           5.1          1.6
84           5.4          3.0           4.5          1.5
85           6.0          3.4           4.5          1.6
86           6.7          3.1           4.7          1.5
87           6.3          2.3           4.4          1.3
88           5.6          3.0           4.1          1.3
89           5.5          2.5           4.0          1.3
90            5.5          2.6           4.4          1.2
91            6.1          3.0           4.6          1.4
92            5.8          2.6           4.0          1.2
93            5.0          2.3           3.3          1.0
94            5.6          2.7           4.2          1.3
95            5.7          3.0           4.2          1.2
96            5.7          2.9           4.2          1.3
97            6.2          2.9           4.3          1.3
98            5.1          2.5           3.0          1.1
99            5.7          2.8           4.1          1.3
100           6.3          3.3           6.0          2.5
101           5.8          2.7           5.1          1.9
102           7.1          3.0           5.9          2.1
103           6.3          2.9           5.6          1.8
104           6.5          3.0           5.8          2.2
105           7.6          3.0           6.6          2.1
106           4.9          2.5           4.5          1.7
107           7.3          2.9           6.3          1.8
108           6.7          2.5           5.8          1.8
109           7.2          3.6           6.1          2.5
110           6.5          3.2           5.1          2.0
111           6.4          2.7           5.3          1.9
112           6.8          3.0           5.5          2.1
113           5.7          2.5           5.0          2.0
114           5.8          2.8           5.1          2.4
115           6.4          3.2           5.3          2.3
116           6.5          3.0           5.5          1.8
117           7.7          3.8           6.7          2.2
118           7.7          2.6           6.9          2.3
119           6.0          2.2           5.0          1.5
120           6.9          3.2           5.7          2.3
121           5.6          2.8           4.9          2.0
122           7.7          2.8           6.7          2.0
123           6.3          2.7           4.9          1.8
124           6.7          3.3           5.7          2.1
125           7.2          3.2           6.0          1.8
126           6.2          2.8           4.8          1.8
127           6.1          3.0           4.9          1.8
128           6.4          2.8           5.6          2.1
129           7.2          3.0           5.8          1.6
130           7.4          2.8           6.1          1.9
131           7.9          3.8           6.4          2.0
132           6.4          2.8           5.6          2.2
133           6.3          2.8           5.1          1.5
134           6.1          2.6           5.6          1.4
135           7.7          3.0           6.1          2.3
136           6.3          3.4           5.6          2.4
137           6.4          3.1           5.5          1.8
138           6.0          3.0           4.8          1.8
139           6.9          3.1           5.4          2.1
140           6.7          3.1           5.6          2.4
141           6.9          3.1           5.1          2.3
142           5.8          2.7           5.1          1.9
143           6.8          3.2           5.9          2.3
144           6.7          3.3           5.7          2.5
145           6.7          3.0           5.2          2.3
146           6.3          2.5           5.0          1.9
147           6.5          3.0           5.2          2.0
148           6.2          3.4           5.4          2.3
149           5.9          3.0           5.1          1.8
")

(defun remove-header-row-names (mat)
  (let* ((row-names (cdr (mapcar #'first mat)))
         (header    (elt mat 0))
         (matrix    (mapcar #'cdr (cdr mat))))
    (values matrix header row-names)))

(defparameter *iris* (parse-whitespace-delimited-string *dstring*))
(defparameter *iris-label* '(0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
                             0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                             0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1
                             1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
                             1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
                             2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
                             2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
                             2 2 2 2 2 2 2 2 2 2))
(multi-bind (*iris-matrix* *iris-header* *iris-rownames*) 
            (remove-header-row-names *iris*))

(defparameter *obs* '(5.9 3.0 5.1 1.8)) ;; should be 2

(defparameter *y-iris* (mapcar (lambda (x) (format nil "~A" x)) 
                       *iris-label*)) ;; as strings

;; (general-stratified-train-test-split-indexes *y-iris* :test-size 0.1)
(multi-bind-destructuring (*train-matrix-iris*
             *test-matrix-iris*
             *train-labels-iris*
             *test-labels-iris*
             *train-indexes-iris*
             *test-indexes-iris*
             *dumped-indexes-iris*
             *categories-iris*
             *split-iris*
             *max-n-iris*
             *train-row-names-iris*
             *test-row-names-iris*
             *header-iris*)
  (general-stratified-train-test-split-from-named-matrix *iris* *y-iris*
  :test-size 0.1))

;;  *test-labels-iris*
;; ("0" "0" "0" "0" "0" "1" "1" "1" "1" "1" "2" "2" "2" "2" "2")

(predict *train-matrix-iris* *train-labels-iris* (elt *test-matrix-iris* 0) 1)
(predict *train-matrix-iris* *train-labels-iris* (elt *test-matrix-iris* 1) 1)
(predict *train-matrix-iris* *train-labels-iris* (elt *test-matrix-iris* 8) 1)
(predict *train-matrix-iris* *train-labels-iris* (elt *test-matrix-iris* 13) (as-creal 0.1))

;; it brings worse problems ... horrible ...


;; but maybe probabilistic neural network is in the last step
;; that thing what could make the last inch ...
;; after hefty pruning of features using the visualization with CNN ...


|#
