#lang racket
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

(require rnrs/arithmetic/bitwise-6)
(require (planet dherman/memoize:3:1))

(define ROOT-SYMBOL -1)

(define MUTATION-RATE 0.2)

;; http://mathworld.wolfram.com/BinomialDistribution.html

(define MAXIMUM-DIFERENCE 0)

(define (set-conditions! dna-value
                         (deviation 4));; 3 = 99.73 percent
  (let* ([expected-difference
          (* dna-value
             MUTATION-RATE)]
         [standard-deviation
          (sqrt (* dna-value
                   MUTATION-RATE
                   (- 1 MUTATION-RATE)))])
    (set! MAXIMUM-DIFERENCE
          (+ expected-difference
             (* deviation
                standard-deviation)))))

;; filename->bitvector : input-port -> (vectorof binary-number)
;; read everything from a given file and returns a 
;; vector containing binary numbers. 
(define (filename->bitvector file-path (radix 2))
  (local [(define (read-accumulative bitvector)
            (let ([line-read (read-line file-path 'any)])
              (if (eof-object? line-read)
                  bitvector
                  (let ([possible-number (string->number line-read radix)])
                    (read-accumulative
                     (if (number? possible-number)
                         (vector-append bitvector (vector possible-number))
                         bitvector))))))]
    (cond
      [(equal? "" file-path)]
      [else
       (read-accumulative (vector))])))


;;  binary-difference : binary-number binary-number -> integer
;; returns the number of bits that are different between 2 binary numbers.
(define/memo (binary-difference bits1 bits2)
  (let ([differences 0])
    (for ((i (in-string (number->string (bitwise-xor bits1 bits2)
                                        2)))
          #:when (equal? i #\1))
      (set! differences (add1 differences)))
    differences))

;; has-kinship? : integer integer -> boolean
;; verify if 2 individuals share a kinship. 
;; The relationship is assumed if the bit difference between the two is
;; less than DEVIATIONS standard deviations from the expected average for
;; the MUTATION-RATE.
;; You must use set-conditions! before using this function, as it relies on
;; external variables.
(define (has-kinship? individual-1 individual-2)
  (< (binary-difference individual-1
                        individual-2)
     MAXIMUM-DIFERENCE))

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
               (when (and (not (= individual-1 individual-2))
                          (has-kinship? individual-1
                                        individual-2))
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
                                  (hash-update! relationship-hash
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

;; compare-results : hash vector -> (listof (listof number)) or true
;; returns a list with this format: 
;; (index-of-element correct-answer found-answer)
(define (compare-results found-hash answers-vector)
  (let ([results empty])
    (for ((individual-index (in-range 0 (vector-length answers-vector))))
      (let ([answer (vector-ref answers-vector individual-index)]
            [found-result (hash-ref found-hash individual-index)])
        (unless (equal? answer found-result)
          (set! results
                (cons (list individual-index
                            answer
                            found-result)
                      results)))))
    (if (empty? results)
        #t
        results)))

;; profiling...
(set-conditions! 500)
(define population (filename->bitvector
                    (open-input-file "small-data")))
(define relationships (make-relationship-hash
                       population))

;(define results empty)
;
;(time (set! results (make-kinship-hash relationships)))
;
;(compare-results
; results
; (filename->bitvector
;  (open-input-file "bitvectors-parents.data.small.txt") 10))


;(make-relationship-hash
; (filename->bitvector
;  (open-input-file "bitvectors-genes.data.smaller")))

(provide filename->bitvector
         binary-difference
         set-conditions!
         has-kinship?
         compare-results
         make-relationship-hash
         make-kinship-hash)