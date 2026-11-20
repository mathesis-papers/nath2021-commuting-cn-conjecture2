# Common-neighborhood energy: groups, rings, surfaces and cliques

Haobo Ma — ChronoAI Pte Ltd, Singapore  
Wenlin Zhang — National University of Singapore; ChronoAI Pte Ltd, Singapore

## Results

The note proves the planar/toroidal and commuting-graph conclusions in
Nath–Fasfous–Das–Shang, *Symmetry* **13** (2021), 1651, Section 4,
Conjectures 1–2, printed page 11, DOI `10.3390/sym13091651`.
For a finite nonabelian group, and also for a finite noncommutative ring
without an identity assumption, put `n = |S \ Z(S)|` and `z = |Z(S)|`.
The noncentral commuting graph satisfies

`E_CN <= n(n-z-2) < 2(n-1)(n-2)`.

Vertices are individual elements, adjacency means commuting, and the
common-neighborhood matrix has zero diagonal. Its energy counts every
eigenvalue with multiplicity. The note also proves centralizer-profile and
least-prime-divisor refinements.

The clique-union conclusion, Conjecture 3, has a prior proof by
Tak–Dutta–Nath, `arXiv:2605.22160v1` (2026), Theorem 5.4, pages 14–15.
The alternative proof here gives its non-strict form and all equality cases:
one clique, or two isolated vertices. The same authors' Corollary 5.5
proves Fasfous–Nath's CC-ring Conjecture 3.12, *Palestine Journal of
Mathematics* **13**(1) (2024), page 75; the ring corollary removes the
CC hypothesis. The coverage table distinguishes these statements from
questions about other graph-energy matrices.

## Build

Build `note.tex` with Tectonic **0.17.0**:

```sh
tectonic --untrusted note.tex
```

Alternatively run `pdflatex -interaction=nonstopmode note.tex` twice.
The bibliography is contained in the main source; BibTeX and external figures
are unnecessary. Normal TeX packages and a first-use network connection may
be needed by Tectonic.

## Mathematical coverage

The universal conclusions have complete written proofs in `note.tex`.
The surface proof uses `E_CN <= 4m-2r`, where `r` counts positive-degree
vertices, Euler's edge bound, and explicitly displayed quotient matrices for
the small dense cases. It does not extrapolate from enumerated examples.

The formal declaration `NoncentralBound.result` states the non-strict
complete-graph bound for all finite nonabelian groups. The following pinned
files specify that statement and the zero-diagonal, full-spectrum definitions:

- https://github.com/the-omega-institute/trureturing/blob/b3bd002bc644799939f299fe6c6c26195e3bb998/D5/S3/Combinatorics/Graph/CommonNeighborhood/NoncentralBound.lean
- https://github.com/the-omega-institute/trureturing/blob/b3bd002bc644799939f299fe6c6c26195e3bb998/D5/S3/Combinatorics/Graph/SubgroupGenerating/CyclicExtensionRefutation.lean

The strict comparison, ring and least-prime-divisor results, surface theorem
and clique equality characterization have separate written proofs.

## Exact finite checks

Python **3.10 or later**, standard library only:

```sh
python3 verify.py --output verification.json
python3 verify_extensions.py --output extension-verification.json
```

`verify.py` checks the matrix identity on 33,868 labelled graphs of orders
0–6, clique energies and equality on 7,337 partitions of orders 1–24,
and 14 named finite group models, including the exact `S4` energy
`17 + sqrt(145)`.

`verify_extensions.py` checks the neighborhood-clique decomposition on the same
33,868 labelled graphs, five noncomplete dense quotient matrices and the five
complete graphs of orders 3–7. It checks exact multiplication, addition,
centralizers and spectra for four finite matrix rings: the four-element row
ring over `F2`, `UT2(F2)`, `M2(F2)` and `UT2(F3)`. The last has CN-energy 160;
the row ring has no multiplicative identity.
All arithmetic is integer or rational. Counts and polynomial coefficients
appear in the two JSON output files.

`CHECKSUMS.sha256` lists every other member of the source archive.
License: **CC BY 4.0** (`cc-by-4.0`); see `LICENSE.txt`.
