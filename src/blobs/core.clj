(ns blobs.core
  (:require [quil.core :as q]
            [quil.middleware :as mw]
            [thi.ng.geom.core :as g]
            [thi.ng.geom.core.vector :as v :refer [vec2 vec3]]
            [thi.ng.geom.core.matrix :as mat]
            [thi.ng.geom.circle :as c]
            [thi.ng.geom.rect :as rect]
            [thi.ng.geom.spatialtree :as accel]
            [thi.ng.geom.physics.core :as phys])
  (:import (thi.ng.geom.physics.core VerletParticle)))

(declare aux-control-points take-nth-closed)

(defmacro draw-poly
  "Macro to draw a filled polygon with verticies given by a sequence of position vectors"
  [verts & {:keys [show-points]
            :or   {:show-points nil}
            :as   args}]
  `(let [verts# ~verts]
     (q/begin-shape)
     (doseq [v# verts#]
       (q/vertex (v# 0) (v# 1)))
     (q/end-shape :close)
     (if (:show-points ~args)
       (do
         ;(q/ellipse (aux2# 0) (aux2# 1) 10 10)
         (doseq [v# verts#]
           (q/ellipse (v# 0) (v# 1) 5 5))
         ;(q/ellipse (aux1# 0) (aux1# 1) 10 10)
         ))))

(defmacro draw-curve
  [verts & {:keys [show-points step]
            :or   {:show-points nil :step 1}
            :as   args}]
  "Macro to interpolate all given position vectors with a closed catmull-rom curve"
  `(let [verts# (take-nth-closed (:step ~args) ~verts)
         [aux1# aux2#] (aux-control-points verts#)
         ]
     (q/begin-shape)
     (q/curve-vertex (aux1# 0) (aux1# 1))
     (doseq [v# verts#]
       (q/curve-vertex (v# 0) (v# 1)))
     (q/curve-vertex (aux2# 0) (aux2# 1))
     (q/end-shape :close)
     (if (:show-points ~args)
       (do
         (doseq [v# ~verts]
           (q/ellipse (v# 0) (v# 1) 5 5))))))

(defn take-nth-closed [n coll]
  "Returns list of every nth element, with the final element being a duplicate of the first"
  (let [coll' (take-nth n coll)]
    (if (= (first coll') (last coll'))
      coll'
      (concat coll' (list (first coll'))))))

(defn aux-control-points
  "Takes a closed path, where the first and last verticies are the same,
  and returns two vectors that represent the beginning and end control
  points of a smooth closed catmull-rom curve that interpolates those vectors"
  [verts]
  (let [p0          (first verts)                           ;first/nth knot
        p1          (second verts)                          ;second knot
        pn-1        ((comp second reverse) verts)           ;(n-1)th knot
        unit-chord1 (g/* (g/- pn-1 p0) (/ 1 (g/mag (g/- pn-1 p0))))
        unit-chord2 (g/* (g/- p1 p0) (/ 1 (g/mag (g/- p0 p1))))
        k1          (g/+ p0 (g/* unit-chord1 (g/dist p0 p1)))
        k2          (g/+ p0 (g/* unit-chord2 (g/dist p0 pn-1)))]
    (list k1 k2)))

(defn attract!
  ;Adapted from https://github.com/thi-ng/demos/blob/master/geom/src/physics_demos/strands.cljs
  [p q rsq strength delta]
  (let [d (g/- p (phys/position q))
        l (+ (g/mag-squared d) 1e-6)]
    (if (< l rsq)
      (phys/add-force q
        (g/* d
          (/ (* (- 1.0 (/ 1 rsq)) (* strength delta)) (Math/sqrt l)))))))

(defn accelerated-force-field
  ;Adapted from https://github.com/thi-ng/demos/blob/master/geom/src/physics_demos/strands.cljs
  [accel r strength]
  (fn [p delta]
    (let [p' (phys/position p)]
      (loop [neighbors (accel/select-with-circle accel p' r)]
        (when-let [n (first neighbors)]
          (if-not (= p n) (attract! p' n (* r r) strength delta))
          (recur (next neighbors)))))))

(defn update-accelerator
  ;Adapted from https://github.com/thi-ng/demos/blob/master/geom/src/physics_demos/strands.cljs
  [accel]
  (fn [physics _]
    (reduce #(g/add-point % (phys/position %2) %2)
      (g/clear! accel)
      (:particles physics))))

(defn particle [pos w lock?]
  "Alternate constructor for a better-behaved particle on initialization"
  (VerletParticle. pos pos (g/clear* pos) lock? nil nil (/ 1.0 w) nil))

(defn make-blob [num-particles spr-len spr-str]
  "Creates a cyclic chain of new particles connected by springs"
  (let [particles (map #(particle
                          (vec2 (+ (/ (q/width) 2) (* 100 (Math/cos (* (/ 2 100) Math/PI %)))) (+ (/ (q/height) 2) (* 100 (Math/sin (* (/ 2 100) Math/PI %))))) 1 false)
                    (range num-particles))
        springs   (cons (phys/spring (last particles) (first particles) spr-len spr-str)
                    (map (fn [[a b]] (phys/spring a b spr-len spr-str)) (partition 2 1 particles)))]
    {:blob-particles particles
     :blob-springs   springs}))

(defn make-init-state
  "Returns state map containing physics information"
  [num-particles spr-len spr-str drag grav ff-dist ff-str]
  (let [blob1         (make-blob num-particles spr-len spr-str)
        blob2         (make-blob num-particles spr-len spr-str)
        all-particles (concat (:blob-particles blob1) (:blob-particles blob2))
        all-springs   (concat (:blob-springs blob1) (:blob-springs blob2))
        accel         (accel/quadtree 0 0 (* 2 (q/width)) (* 2 (q/height)))
        r             {:c (phys/shape-constraint-inside
                            (c/circle (/ (q/width) 2) (/ (q/height) 2) (/ (q/height) 2.3)))}]

    (doseq [p all-particles] (phys/add-constraints p r))
    {:physics           (phys/physics
                          {:blobs     (list blob1 blob2)
                           :particles all-particles
                           :springs   all-springs
                           :drag      drag
                           :behaviors {:g (phys/gravity grav)
                                       :f (accelerated-force-field accel ff-dist ff-str)}
                           :listeners {:iter (update-accelerator accel)}
                           })
     :constraint        r
     :blobs             (list blob1 blob2)
     :quadtree          accel
     :selected-particle nil
     :draw-points       false
     :step              1}))

(defn add-to-chain
  "Creates a new particle at position, inserts it into blob before
  the target particle, and returns a new state where this is reflected in
  :blobs and :physics"
  [state position target-particle]
  (let [blob          (first (filter #(some #{target-particle} (:blob-particles %)) (:blobs state)))
        rest-blobs    (remove #(= blob %) (:blobs state))
        [s-h s-t] (split-with #(not= target-particle (:b %)) (:blob-springs blob))
        [p-h p-t] (split-with #(not= target-particle %) (:blob-particles blob))
        old-spr       (first s-t)
        new-part      (particle position 1 false)
        updated-spr   (assoc old-spr :b new-part)
        new-spr       (phys/spring new-part target-particle (:rest-len old-spr) (:strength old-spr))
        new-blob      {:blob-particles (concat p-h (cons new-part p-t))
                       :blob-springs   (concat s-h (list updated-spr new-spr) (rest s-t))}
        new-particles (cons new-part (-> state :physics :particles))
        new-springs   (cons new-spr (replace {old-spr updated-spr} (-> state :physics :springs)))
        new-physics   (assoc (:physics state) :particles new-particles
                                              :springs new-springs)]
    ;(println "found blob" blob)
    (phys/add-constraints new-part (:constraint state))
    (assoc state :physics new-physics
                 :blobs (cons new-blob rest-blobs))))

(defn setup []
  (q/frame-rate 30)
  (make-init-state 30 5 1.5 0.01 (vec2 0 0) 25 -1.5))

(defn update-state [state]
  (let [state (if (nil? state) (setup) state)
        s     (atom state)]
    (phys/timestep (:physics @s) 2)
    @s))

(defn draw-state [state]
  (q/background 22 131 131)
  (q/fill 255)
  (q/stroke-weight 3)
  ;(println (:blobs state))
  (doseq [blob (:blobs state)]
    ;(println blob)
    (let [b-particles (:blob-particles blob)
          b-verts     (map phys/position b-particles)]
      (draw-curve b-verts
        :show-points (:draw-points state)
        :step (:step state)))))

(defn key-pressed [old-state event]
  ;(println event)
  (case (:key-code event)
    32 (update-in old-state [:draw-points] not)             ;space
    37 (update-in old-state [:step] #(if (> % 1) (dec %) %)) ;left
    39 (update-in old-state [:step] inc)                    ;right
    old-state))

(defn mouse-released [old-state event]
  (when (:selected-particle old-state)
    (phys/unlock (:selected-particle old-state)))
  (assoc-in old-state [:selected-particle] nil))

(defn mouse-pressed [old-state event]
  ;(println event)
  (let [pos (vec2 (:x event) (:y event))
        particle-of-interest
            (first (accel/select-with-circle (:quadtree old-state) pos 10))]
    (if (#{:left :right} (:button event))
      (assoc old-state :selected-particle particle-of-interest)
      (if particle-of-interest
        (add-to-chain old-state pos particle-of-interest)
        (assoc old-state :selected-particle particle-of-interest)))))

(defn mouse-dragged [old-state event]
  (let [mx                   (:x event)
        my                   (:y event)
        particle-of-interest (:selected-particle old-state)]
    (when particle-of-interest
      (phys/set-position particle-of-interest (vec2 mx my))
      (if (= (:button event) :left)
        (phys/unlock particle-of-interest)
        (phys/lock particle-of-interest)))
    old-state))

(q/defsketch blobs
  :title "Blobs-- running"
  :size [500 500]
  :setup setup
  :update update-state
  :draw draw-state
  :key-pressed key-pressed
  :mouse-dragged mouse-dragged
  :mouse-released mouse-released
  :mouse-pressed mouse-pressed
  :features [:resizable :keep-on-top]
  :middleware [mw/fun-mode mw/pause-on-error])



