# Locating the integer slit heights from the four working intervals

Let the stopband be \(A\), and write the three prescribed passbands as
\(B'_1,B'_2,B'_3\) in increasing real order. Compute the canonical SC
numerator and the three positive passband heights \(H'_j\) directly from
these eight endpoints, treating the supplied bands as full bands for this
auxiliary calculation. Put

\[
q_1=n\frac{H'_1}{H'_1+H'_2+H'_3},\qquad
q_2=n\frac{H'_1+H'_2}{H'_1+H'_2+H'_3}.
\]

Suppose that the prescribed bands do arise by two allowed one-tile trims
of a solution with indices \(0<m_1<m_2<n\), period \(T>0\), and trim lengths
\(\delta_1,\delta_2\in[0,1]\) in units of \(T\). Then

\[
 |q_j-m_j|\leq\delta_1+\delta_2\leq 2,\qquad j=1,2.
\]

Thus in exact arithmetic each closed interval \([q_j-2,q_j+2]\) contains
at most five candidate integers. Removing unordered and out-of-range pairs
leaves at most 25 pairs and at most 100 assignments after both trim sides
are varied. The bound holds before solving for the hidden prevertices and
uses no estimate of the conditioning of that solve.

## Proof by comparison of the two condenser potentials

Write \(B_j\) for the unknown full passbands and set
\(B=\bigcup_j B_j\), \(B'=\bigcup_j B'_j\). Each \(B'_j\subseteq B_j\), and
the stopband \(A\) is common. Let \(u\) be the real part of the normalized
SC primitive for \((A,B)\), and let \(u'\) be the corresponding function for
\((A,B')\). Use the normalization \(-1\) on the stopband and \(+1\) on
each passband in both problems.

For either set of eight ordered endpoints the SC differential is real along
the open gaps and the two exterior rays. The imaginary part of the primitive
is consequently constant on each such segment. Subtracting this local
imaginary constant gives Schwarz reflection across the segment. In particular
the real part extends by \(u(\bar z)=u(z)\). It is harmonic on the sphere
cut along the four conducting intervals, continuous up to both banks and
at their endpoints. At a finite endpoint the SC primitive has a finite
limit with a square-root local coordinate. At infinity its differential
is \(O(z^{-2})\), so the bounded reflected real part extends harmonically
there as well. The maximum principle now gives \(-1\leq u,u'\leq1\).
These facts also apply to the auxiliary SC primitive for the supplied
endpoints without assuming any integer-height condition on it.

Compare \(w=u-u'\) on the common domain
\(D=\widehat{\mathbb C}\setminus(A\cup B)\). The function is harmonic on
\(D\). Its limiting boundary values are zero on \(A\cup B'\) and are
\(1-u'\geq0\) on the deleted passband arcs. Thus \(w\geq0\) on \(D\).
At an interior point of \(A\) or a retained passband, \(w(x)=0\) and
\(w(x+iy)\geq0\) for \(y>0\), so \(\partial_y w(x,0^+)\geq0\).
Both primitives have regular one-sided derivatives there. This yields

\[
 -\partial_y u'\geq-\partial_y u\quad\text{on the retained passbands},
 \qquad
 \partial_y u'\leq\partial_y u\quad\text{on the stopband}.
\]

By Cauchy–Riemann, \(-\partial_y u=\partial_x\Im U\). Integrating this
positive density along a full passband gives its SC height. On the stopband
the sign reverses. Height closure identifies the integral of \(\partial_y u\)
over the stopband with the sum \(C\) of the three positive full heights.
The same identity holds for the auxiliary geometry with total \(C'\).
Consequently \(0<C'\leq C\). All these are upper-bank integrals. The possible
inverse-square-root endpoint densities are integrable and contribute no
point masses.

Define probability measures \(\nu\) and \(\nu'\) on the real line by
integrating \(-\partial_y u/C\) on \(B\) and \(-\partial_y u'/C'\) on \(B'\),
respectively, extending each by zero elsewhere. On retained arcs the density
inequality and \(C'\leq C\) imply \(\nu'\geq\nu\). The old mass of the deleted
arcs is exactly

\[
 d=\nu(B\setminus B')=\frac{(\delta_1+\delta_2)T}{nT}.
\]

Both measures have total mass one. The signed measure \(\nu'-\nu\) therefore
has positive mass \(d\) on the retained set and negative mass of magnitude
\(d\) on the deleted set. For every measurable \(S\),
\(|\nu'(S)-\nu(S)|\leq d\). Take \(S=B_1\cup\cdots\cup B_j\), \(j=1,2\).
Its old mass is \(m_j/n\) by quantization and its new mass is \(q_j/n\).
Multiplying by \(n\) gives the stated bound.

## What the numerical search establishes

The exact bound concerns the SC integrals of the working endpoints.
Evaluating these integrals in floating point also introduces error into
\(q_j\), which is multiplied by \(n\). A numerical guard expands the list of
integer candidates. Guaranteed inclusion under rounding would require a
certified enclosure of \(q_j\) or a verified error smaller than the guard.
Each trial also requires a two-variable numerical inverse. Its failure can
mean insufficient accuracy, poor conditioning, an unsuitable initial point,
or absence of a solution in that mode. Exhausting the finite trial list in
floating point does not prove nonexistence. The window applies to the marked
one-stopband/three-passband family with its two allowed trims. The comparison
of harmonic normal derivatives is an analytical proof. The subsequent mass bound, tile rescaling and five-integer count are checked in the companion Lean project (Stiefel/ModeWindow.lean).
