(define colors '(blue green red))

(define samples
  (mh-query
   200 100

   (define phi (dirichlet '(1 1 1)))
   (define alpha 0.1)
   (define prototype (map (lambda (w) (* alpha w)) phi))

   (define bag->prototype (mem (lambda (bag) (dirichlet prototype))))

   (define obs->bag
     (mem (lambda (obs-name)
            (uniform-draw '(bag1 bag2 bag3)))))

   (define draw-marble
     (mem (lambda (obs-name)
            (multinomial colors (bag->prototype (obs->bag obs-name))))))

   (list (equal? (obs->bag 'obs1) (obs->bag 'obs2))
         (equal? (obs->bag 'obs1) (obs->bag 'obs3)))

   (and
    (equal? 'red (draw-marble 'obs1))
    (equal? 'red (draw-marble 'obs2))
    (equal? 'blue (draw-marble 'obs3))
    (equal? 'blue (draw-marble 'obs4))
    (equal? 'red (draw-marble 'obs5))
    (equal? 'blue (draw-marble 'obs6))
    )))