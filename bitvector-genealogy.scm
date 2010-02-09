#lang scheme

(require test-engine/scheme-tests)

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

;; string->bit-vector : symbol -> (vectorof char)
;; turns a long bit symbol into a vector of bits.
(define (string->bit-vector the-string)
  (local [(define (string->bit-vector-acc bits accumulator)
            (cond
              [(>= accumulator (string-length the-string)) bits]
              [else
               (string->bit-vector-acc
                (vector-append bits
                               (vector (string-ref the-string
                                                   accumulator)))
                (add1 accumulator))]))]
    (string->bit-vector-acc (vector) 0)))

(check-expect (string->bit-vector "01101")
              (vector #\0 #\1 #\1 #\0 #\1))

(check-expect (string->bit-vector "")
              (vector))


;; filename->bit-vector : input-port -> (vectorof (vectorof char))
;; read everything from a given file and returns a 
;; vector containing vectors of bits.
(define (filename->bit-vector file-path)
  (local [(define (read-accumulative accumulator)
            (let ([line-readed (read-line file-path)])
              (cond
                [(eof-object? line-readed) accumulator]
                [else
                 (read-accumulative
                  (vector-append
                   accumulator
                   (vector (string->bit-vector line-readed))))])))]
    (cond
      [(equal? "" file-path)]
      [else
       (read-accumulative (vector))])))


(test)