;;;; cl-pnn.asd

(asdf:defsystem #:cl-pnn
  :description "Probabilistic Neural Network (PNN) package for common lisp for big data"
  :author "Gwang-Jin Kim"
  :license  "free one"
  :version "0.0.1"
  :serial t
  :components ((:file "package")
               (:file "cl-pnn")))
