#lang racket

(require rackunit)

;;  The BitVectors are an ancient and immortal race of 10,000, each with a 10,000
;;  bit genome. The race evolved from a single individual by the following process:
;;  9,999 times a BitVector chosen at random from amongst the population was cloned
;;  using an error-prone process that replicates each bit of the genome with 80%
;;  fidelity (any particular bit is flipped from parent to child 20% of the time,
;;  independent of other bits).  
;;  
;;  Write a program to guess the reproductive history of BitVectors from their
;;  genetic material. The randomly-ordered file bitvectors-genes.data.gz contains a
;;  10,000 bit line for each individual. Your program's output should be, for each
;;  input line, the 0-based line number of that individual's parent, or -1 if it is
;;  the progenitor. Balance performance against probability of mistakes as you see
;;  fit.
(define ROOT-SYMBOL -1)

(define MUTATION-RATE 0.2)

(define DNA-LENGTH 0)

;; http://mathworld.wolfram.com/BinomialDistribution.html

(define DEVIATION 3) ;; 3 = 99.73 percent

(define (get-dna)
  DNA-LENGTH)

(define (set-dna! number)
  (set! DNA-LENGTH number))

;; string->bitvector : string -> binary-number
;; turns a long bit symbol into a binary-number.
(define (string->binary-number bit-sequence)
  (when (> (string-length bit-sequence) 0)
    (string->number bit-sequence 2)))


;; filename->bitvector : input-port -> (vectorof binary-number)
;; read everything from a given file and returns a 
;; vector containing binary numbers. 
;; Side Effect : Sets the DNA-LENGTH.
(define (filename->bitvector file-path)
  (local [(define (read-accumulative index)
            (let ([line-readed (read-line file-path)])
              (cond
                [(eof-object? line-readed) index]
                [else
                 (read-accumulative
                  (begin
                    (set-dna! (add1 (get-dna)))
                    (vector-append
                     index
                     (vector (string->binary-number line-readed)))))])))]
    (cond
      [(equal? "" file-path)]
      [else
       (read-accumulative (vector))])))


;;  binary-difference : binary-number binary-number -> integer
;; returns the number of bits that are different between 2 binary numbers.
(define (binary-difference bits1 bits2)
  (let ([differences 0])
    (for ((i (in-string (number->string (bitwise-xor bits1 bits2) 2))))
      (when (equal? i #\1)
        (set! differences (add1 differences))))
    differences))

;; %-difference : binary-number binary-number -> number
;; returns a number between 0 and 1 representing the % of difference between 
;; 2 binary numbers and the current global dna value.
(define (%-difference individual-1
                      individual-2
                      (dna-length (get-dna)))
  (/ (binary-difference individual-1 individual-2)
     dna-length))

;; binomial-standard-deviation : number -> number
;; Compute the expected binomial standard deviation for N events of probability P.
(define (binomial-standard-deviation n p)
  (sqrt (* n p (- 1 p))))

;; has-kinship? : integer integer -> boolean
;; verify if 2 individuals share a kinship. 
;; The relationship is assumed if the bit difference between the two is
;; less than DEVIATIONS standard deviations from the expected average for
;; the MUTATION-RATE.
(define (has-kinship? individual-1
                      individual-2
                      (mutation-rate MUTATION-RATE)
                      (deviations DEVIATION)
                      (dna-length (get-dna)))
  (let* ([expected-difference (* dna-length mutation-rate)]
         [standard-deviation
          (binomial-standard-deviation dna-length
                                       mutation-rate)]
         [max-difference (+ expected-difference
                            (* deviations standard-deviation))])
    (< (binary-difference individual-1
                          individual-2)
       max-difference)))



;; create-relationship-hash : (vectorof number) -> hash
;; creates a hash mapping an individual to it's list of relationships.
;; Both parents and children are mapped to an individual.
(define (make-relationship-hash bitvector)
  (let ([relationship-hash (make-hash)]
        [dna-range (in-range 0 (vector-length bitvector))])
    (for* ((individual-1-index dna-range)
           (individual-2-index dna-range))
      (begin (unless (hash-has-key? relationship-hash individual-1-index)
               (hash-set! relationship-hash individual-1-index empty))
             (let ([individual-1 (vector-ref bitvector
                                             individual-1-index)]
                   [individual-2 (vector-ref bitvector
                                             individual-2-index)])
               (when (and (has-kinship? individual-1
                                        individual-2)
                          (not (= individual-1 individual-2)))
                 (hash-update! relationship-hash
                               individual-1-index
                               (λ (relationships)
                                 (cons individual-2-index relationships)))))))
    relationship-hash))

;Write a program to guess the reproductive history of BitVectors from their genetic material. The randomly-ordered file bitvectors-genes.data.gz contains a 10,000 bit line for each individual. Your program's output should be, for each input line, the 0-based line number of that individual's parent, or -1 if it is the progenitor. Balance performance against probability of mistakes as you see fit. 

;; make-kinship-hash : hash -> hash
;; creates a hash mapping an individual only to it's parent.
;    * Find all the bit vectors with only a single parent-child relationship
;    * Assume that those vectors are the child nodes
;    * Record that information
;    * Remove the child node from the parent's list of relationships
;    * Repeat
(define (make-kinship-hash relationship-hash-initial)
  (local [(define (make-kinships kinship-hash
                                 relationship-hash)
            (local [(define (process-possible-kinship individual-index
                                                      relationships-indexes)
                      (let ([number-of-relationships
                             (length relationships-indexes)])
                        (unless (> number-of-relationships 1)
                          (begin
                            (hash-remove! relationship-hash
                                          individual-index)
                            (if (= number-of-relationships 1)
                                ; leaf node
                                (begin
                                  (hash-set! kinship-hash
                                             individual-index
                                             (first relationships-indexes))
                                  (hash-update!
                                   relationship-hash
                                   (first relationships-indexes)
                                   (λ (relationships)
                                     (remove individual-index
                                             relationships))))
                                ; root node
                                (begin (hash-remove! relationship-hash
                                                     individual-index)
                                       (hash-set! kinship-hash
                                                  individual-index
                                                  ROOT-SYMBOL)))))))]
              (if (false? (hash-iterate-first relationship-hash))
                  kinship-hash
                  (begin
                    (hash-for-each relationship-hash
                                   process-possible-kinship)
                    (make-kinships kinship-hash
                                   relationship-hash)))))]
    (make-kinships (make-hash) relationship-hash-initial)))


(provide filename->bitvector
         binary-difference
         get-dna
         set-dna!
         has-kinship?
         %-difference
         make-relationship-hash
         make-kinship-hash
         string->binary-number)