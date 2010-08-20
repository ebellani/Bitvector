#lang racket

(require rackunit
         rackunit/text-ui
         "bitvector-genealogy.scm")

(define utilities-suite
  (test-suite
   "Testing utilities"
   (test-case
    "Testing retrieving the genetic information from a file"
    (check-equal? (vector-length
                   (filename->bitvector
                    (open-input-file "bitvectors-genes.data.small")))
                  500)
    (check-equal? (filename->bitvector
                   (open-input-file "bitvectors-genes.data.tiny"))
                  #(#b10100011
                    #b10000000
                    #b11110001
                    #b10011110))
    (check-equal? (binary-difference #b11111001 #b1001) 4)
    (check-equal? (binary-difference #b01101 #b10001) 3)
    (check-equal? (binary-difference #b01101 #b00000001101) 0))))

(define dna-suite
  (test-suite
   "Testing everything that depends on the size of the DNA"
   #:before (lambda () (set-dna! 10))
   (test-case
    "Testing the % of difference between 2 individuals"
    
    (check-equal? (exact->inexact (%-difference #b0110101010
                                                #b1110100010))
                  0.2)
    (check-equal? (%-difference #b0110101010
                                #b1001010101)
                  1)
    (check-equal? (%-difference #b0110101010
                                #b0110101010)
                  0)
    (check-equal? (exact->inexact (%-difference #b0110101010
                                                #b0101010101))
                  0.8))
   (test-case
    "Test if there is kinship between 2 individuals"
    (check-equal? (has-kinship? #b0110101010
                                #b1001100101)
                  #f)
    (check-equal? (has-kinship? #b0110101010
                                #b0111101110)
                  #t))
   
   (test-case
    "Test if there is we can create a table of relationships from a vector of individuals"
    (check-equal? (make-relationship-hash #(#b0110101010
                                            #b1001100101
                                            #b1001100100
                                            #b1001101100
                                            #b0110101011))
                  (make-hash (list '(0 . (4))
                                   '(4 . (0))
                                   '(3 . (2 1))
                                   '(2 . (3 1))
                                   '(1 . (3 2)))))
    
    (check-equal? (make-relationship-hash
                   (filename->bitvector
                    (open-input-file "bitvectors-genes.data.little.txt")))
                  (make-hash (list '(0 . ())
                                   '(5 . ())
                                   '(4 . ())
                                   '(3 . ())
                                   '(2 . ())
                                   '(1 . ())))))
   
   (test-case
    "Test if there is we can create a table of kinships, mapping individual to its parent"
    (check-equal? (make-kinship-hash
                   (make-hash (list '(0 . (1 2))
                                    '(1 . (0))
                                    '(2 . (0 3))
                                    '(3 . (2)))))
                  (make-hash (list '(0 . -1)
                                   '(3 . 2)
                                   '(2 . 0)
                                   '(1 . 0)))))
   (test-case
    "Test to see if I can compare results"
    (check-equal? (compare-results
                   (make-hash (list '(0 . -1)
                                    '(3 . 2)
                                    '(2 . 0)
                                    '(1 . 0)))
                   #(-1 0 0 2))
                  #t)
    (check-equal? (compare-results
                   (make-hash (list '(0 . -1)
                                    '(3 . -1)
                                    '(2 . 0)
                                    '(1 . 2)))
                   #(-1 0 0 2))
                  (list '(3 2 -1)
                        '(1 0 2))))))

(define small-dna-suite
  (test-suite
   "Testing including the small sample data"
   #:before (lambda () (set-dna! 500))
   
   (test-case
    "Check what the result is, compared to mine."
    ;(list '(347 119 58)
    ;      '(284 -1 125)
    ;      '(211 28 -1)
    ;      '(125 284 -1)
    ;      '(58 347 -1))
    
    ;    (list '(125 284 399)
    ;          '(211 28 -1)
    ;          '(284 -1 125)
    ;          '(347 119 -1)
    ;          '(399 125 -1))
    (check-equal? (compare-results
                   (make-kinship-hash
                    (make-relationship-hash
                     (filename->bitvector
                      (open-input-file "bitvectors-genes.data.small"))))
                   (filename->bitvector
                    (open-input-file "bitvectors-parents.data.small.txt") 10))
                  #f))))



;(run-tests utilities-suite 'verbose)
;(run-tests dna-suite 'verbose)
;(run-tests small-dna-suite 'verbose)
