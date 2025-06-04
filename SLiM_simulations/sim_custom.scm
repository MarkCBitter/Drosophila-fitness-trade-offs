(define (r-from-ps pAB pA pB) (/ (- pAB (* pA pB))
                                (sqrt (* pA pB (- 1 pA) (- 1 pB)))))

;; from VanLiere & Rosenberg (2008)
;; returns (list rmin rmax)
(define (r-bounds pA pB) (cond 
        ;; segment S1 (pA=0.4, pB=0.9)
        [(and (<= (- 1 pA) pB) (<= pA 0.5) (> pB 0.5)) 
            ;; on S1 r is minimised at pAB = pA + pB - 1
            ;;        and maximized at pAB = pA (= min(pA, pB))
            (list 
                (r-from-ps (+ pA pB -1) pA pB)
                (r-from-ps pA pA pB))]
        ;; segment S2 (pA=0.6, pB=0.9)
        [(and (<= pA pB) (> pA 0.5) (> pB 0.5)) 
            ;; on S2 r is minimised at pAB = pA + pB - 1
            ;;        and maximized at pAB = pA
            (list 
                (r-from-ps (+ pA pB -1) pA pB)
                (r-from-ps pA pA pB))]
        ;; segment S3 (pA=0.9, pB=0.6)
        [(and (> pA pB) (> pA 0.5) (> pB 0.5)) 
            ;; on S3 r is minimised at pAB = pA + pB - 1
            ;;        and maximized at pAB = pB
            (list 
                (r-from-ps (+ pA pB -1) pA pB)
                (r-from-ps pB pA pB))]
        ;; segment S4 (pA=0.8, pB=0.4)
        [(and (<= (- 1 pA) pB) (> pA 0.5) (<= pB 0.5)) 
            ;; on S3 r is minimised at pAB = pA + pB - 1
            ;;        and maximized at pAB = pB
            (list 
                (r-from-ps (+ pA pB -1) pA pB)
                (r-from-ps pB pA pB))]
        ;; segment S5 (pA=0.7, pB=0.1)
        [(and (> (- 1 pA) pB)  (> pA 0.5) (<= pB 0.5)) 
            ;; on S5 r is minimised at pAB = 0
            ;;        and maximized at pAB = pB
            (list 
                (r-from-ps 0.0 pA pB)
                (r-from-ps pB pA pB))]
        ;; segment S6 (pA=0.4, pB=0.3)
        [(and (> pA pB) (<= pA 0.5) (<= pB 0.5)) 
            ;; on S6 r is minimised at pAB = 0
            ;;        and maximized at pAB = pB
            (list 
                (r-from-ps 0.0 pA pB)
                (r-from-ps pB pA pB))]
        ;; segment S7 (pA=0.2 pB=0.4)
        [(and (<= pA pB) (<= pA 0.5) (<= pB 0.5)) 
            ;; on S7 r is minimised at pAB = 0
            ;;        and maximized at pAB = pA
            (list 
                (r-from-ps 0.0 pA pB)
                (r-from-ps pA pA pB))]
        ;; segment S8 (pA=0.2, pB=0.6)
        [(and (> (- 1 pA) pB)  (<= pA 0.5) (> pB 0.5)) 
            ;; on S8 r is minimised at pAB = 0
            ;;        and maximized at pAB = pA
            (list 
                (r-from-ps 0.0 pA pB)
                (r-from-ps pA pA pB))]
))

;; pAB defined based on r rather than r-squared
;; allows to produce different values for attraction and repulsion
(define (pAB-from-r pA pB r)
    (+ (* pA pB) 
       (* r (sqrt (* pA pB (- 1 pA) (- 1 pB))))
    )
)

; ;; sample initial frequencies
(define pA (random-range-float f0min f0max))
(define pB (random-range-float f0min f0max))
(define bounds (r-bounds pA pB))

; ;; define limits for r
;; unlinked is when r-squared < 0.01
(define r-unlinked-min -0.1)
(define r-unlinked-max 0.1)
;; attraction is when r-squared > 0.05, with r > 0
(define r-attraction-min 0.2236)
(define r-attraction-max (cadr bounds))
;; repulsion is when r-squared > 0.05, with r < 0
(define r-repulsion-max -0.2236)
(define r-repulsion-min (car bounds))

;; define r based on input category
(define r 
    (cond 
        ;; sample r in repulsion
        [(< rsqCategory -0.5) 
            (random-range-float r-repulsion-min r-repulsion-max)]
        ;; sample r in attraction
        [(> rsqCategory  0.5) 
            (random-range-float r-attraction-min r-attraction-max)]
        ;; sample unlinked r
        [else 
            (random-range-float r-unlinked-min r-unlinked-max)]))

(define nA  (ceiling (* 2 nSeed pA)))
(define nB  (ceiling (* 2 nSeed pB)))
(define nAB (floor (* 2 nSeed (pAB-from-r pA pB r))))

(parameters nA nB nAB r)