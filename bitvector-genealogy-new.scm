#lang racket

;; starting over. This time with new references.
;; 2 main objectives now:
;; * faster by doing less things and caching (memoization?) functions,
;;   using for in the faster way.

;; * get the roots right
;;  http://www.reddit.com/r/programming/comments/5z74c/new_ita_software_puzzle_bitvector_genealogy/
;; http://web.itu.edu.tr/~yazicivo/files/bitvector-genealogy.lisp.txt

(provide filename->population
         kinship-probability
         population->probabilities-sum
         set-conditions!
         hamming-weigth)

(require rnrs/arithmetic/bitwise-6)
(require (planet dherman/memoize:3:1))

(define ROOT-SYMBOL -1)

(define MUTATION-RATE 0.2)

(define FIDELITY-RATE
  (- 1 MUTATION-RATE))

(define GENOME-LENGTH 0)
(define POPULATION-SIZE 0)

(define MAXIMUM-DIFFERENCE 0)

;(define-syntax-rule (while condition body)
;  (let loop ()
;    (when condition
;      body (loop))))


(define (set-conditions! genome-length)
  (set! GENOME-LENGTH genome-length)
  (set! MAXIMUM-DIFFERENCE
        (let* ([deviation 3] ;; 3 = 99.73 percent
               [expected-difference
                (* GENOME-LENGTH MUTATION-RATE)]
               [standard-deviation
                (sqrt (* GENOME-LENGTH
                         MUTATION-RATE
                         (- 1 MUTATION-RATE)))])
          (+ expected-difference
             (* deviation standard-deviation)))))

;; filename->population : input-port -> (hashof integer integer)
;; read everything from a given file and maps the index to 
;; the content of the file.
(define (filename->population file-path (radix 2))
  (local [(define population (make-hash))
          (define (read-accumulative index)
            (let ([line-read (read-line file-path 'any)])
              (if (eof-object? line-read)
                  (begin
                    (set! POPULATION-SIZE index) population)
                  (begin (hash-set! population
                                    index
                                    (string->number line-read radix))
                         (read-accumulative (add1 index))))))]
    (read-accumulative 0)))

;; hamming-weight : string string -> integer
;; http://en.wikipedia.org/wiki/Hamming_weight
;; counts the number of 1's in the XOR of 2 binary numbers.
(define/memo (hamming-weigth bits-1 bits-2)
  (for/fold ([differences 0])
    ((i (in-string (number->string (bitwise-xor bits-1 bits-2)
                                   2)))
     #:when (equal? i #\1))
    (add1 differences)))

; kinship-probability : number number hash -> number
; verify the probability that 2 individuals share a kinship.
(define (kinship-probability genome-1
                             genome-2
                             population)
  (let ([differences (/ (hamming-weigth (hash-ref population
                                                  genome-1)
                                        (hash-ref population
                                                  genome-2))
                        GENOME-LENGTH)])
    (* (expt FIDELITY-RATE
             (- 1 differences))
       (expt MUTATION-RATE
             differences))))

;; population->relationships : hash -> hash
;; creates a hash mapping an individual to it's list of relationships.
;; Both parents and children are mapped to an individual.
(define (population->relationships population)
  (let ([relationships (make-hash)])
    (for* ((genome-1-index (in-range 0 POPULATION-SIZE))
           (genome-2-index (in-range 0 POPULATION-SIZE))
           #:when (and (not (= genome-1-index
                               genome-2-index))))
      (if (not (hash-has-key? relationships genome-1-index))
          (hash-set! relationships genome-1-index empty)
          (when (has-kinship? genome-1-index
                              genome-2-index
                              population)
            (hash-update! relationships
                          genome-1-index
                          (λ (current-relationships)
                            (cons genome-2-index
                                  current-relationships))
                          empty))))
    relationships))

;; has-kinship? : integer integer hash -> boolean
;; verify if 2 individuals share a kinship. 
;; The relationship is assumed if the bit difference between the two is
;; less than DEVIATIONS standard deviations from the expected average for
;; the MUTATION-RATE.
;; You must use set-conditions! before using this function, as it relies on
;; external variables.
(define (has-kinship? genome-1
                      genome-2
                      population)
  (< (hamming-weigth (hash-ref population
                               genome-1)
                     (hash-ref population
                               genome-2))
     MAXIMUM-DIFFERENCE))

;; is-source? : integer integer hash -> boolean
;; Figures if a genome is the source for another.
(define (is-source? possible-source
                    possible-child
                    probabilities-sum)
  (> (hash-ref probabilities-sum
               possible-source)
     (hash-ref probabilities-sum
               possible-child)))


;; population->probabilities-sum  : hash -> integer
;; sums up all probabilities and finds the item with the 
;; largest probabilities to all elements
(define (population->probabilities-sum population)
  (let ([population-probabilities
         (for*/fold ([probabilities (hash)])
           ([individual-1-key (in-range 0 POPULATION-SIZE)]
            [individual-2-key (in-range 0 POPULATION-SIZE)]
            #:when (not (= individual-1-key individual-2-key)))
           (hash-update probabilities
                        individual-1-key
                        (λ (probability-sum)
                          (+ (kinship-probability individual-1-key
                                                  individual-2-key
                                                  population)
                             probability-sum))
                        0))])
    population-probabilities))



;; compare-results : hash hash -> (listof (listof number)) or true
;; returns a list with this format:
;; (index-of-element correct-answer found-answer)
(define (compare-results found-results answers)
  (let ([results empty])
    (for ((individual-index (in-range 0 POPULATION-SIZE)))
      (let* ([answer (hash-ref answers individual-index)]
             [found-result (hash-ref found-results
                                     individual-index
                                     (λ ()
                                       (printf "missing ~s.~n"
                                               individual-index)))])
        (unless (equal? answer found-result)
          (set! results (cons (list individual-index
                                    answer
                                    found-result)
                              results)))))
    results))

;; find-source-in-relationships : integer (listof integer) hash ->
;;                                number or false
;; finds who is the source of the genome from all its relationships,
;; returns false if the element is the source of all its relatinships (root).
(define (find-source-in-relationships genome
                                      relationships
                                      probabilities-sum)
  (let ([sources (filter (λ (another-genome)
                           (is-source? another-genome
                                       genome
                                       probabilities-sum))
                         relationships)])
    (if (empty? sources)
        #f
        (if (= (length sources) 1)
            (first sources)
            (first (sort sources
                         (λ (genome1 genome2)
                           (is-source? genome1 genome2 probabilities-sum))))))))

;; find-genealogy : hash hash -> hash
;; creates a mapping from an individual only to it's parent.
(define (find-genealogy relationships
                        probabilities-sum)
  (let ([genealogy (make-hash)])
    (for ([genome-key (in-range 0 POPULATION-SIZE)])
      (let ([source (find-source-in-relationships
                     genome-key
                     (hash-ref relationships genome-key)
                     probabilities-sum)])
        (if (false? source)
            (hash-set! genealogy
                       genome-key
                       ROOT-SYMBOL)
            (hash-set! genealogy
                       genome-key
                       source))))
    genealogy))

(define (find-genealogy2 population probabilities-sum (passes 3))
  (let ([genealogy (make-hash)])
    (for* ([genome-1-key (in-range 0 POPULATION-SIZE)]
           [genome-2-key (in-range 0 POPULATION-SIZE)]
           #:when (not (= genome-1-key
                          genome-2-key)))
      (let* ([found-possible-source (hash-ref genealogy
                                              genome-1-key
                                              #f)]
             [updater
              (if (has-kinship? genome-1-key
                                genome-2-key
                                population)
                  (if (not (false? found-possible-source))
                      (if (is-source? genome-2-key
                                      found-possible-source
                                      probabilities-sum)
                          genome-2-key
                          found-possible-source)
                      genome-2-key)
                  #f)])
        (when (not (false? updater))
          (hash-update! genealogy
                        genome-1-key
                        (λ (x)
                          updater)
                        #f))))
    genealogy))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; TESTS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(set-conditions! 500)
(define population (filename->population
                    (open-input-file "small-data")))

(time (population->relationships population))

(define relationships (population->relationships population))
(define probabilities-sum (population->probabilities-sum population))
(define results empty)
(define results2 empty)

(time (set! results (find-genealogy relationships
                                    probabilities-sum)))
(time (set! results2 (find-genealogy2 population
                                      probabilities-sum)))

(define answers (filename->population
                 (open-input-file "bitvectors-parents.data.small.txt")
                 10))


(display "results 1 = ")
(compare-results results answers)

(display "results 2 = ")
(compare-results results2 answers)

;; find-root : hash -> integer
;; given the probabilities, find the one that has the biggest 
;; chance of being the root
;(define (find-root population-probabilities)
;  (car (first (sort (hash-map population-probabilities
;                              (λ (x y)
;                                (cons x y)))
;                    (λ (x y)
;                      (> (cdr x)
;                         (cdr y)))))))


