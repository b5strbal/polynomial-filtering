# polynomial-filtering

Code for searching for polynomials that can possibly arise as characteristic polynomials of the action on homology for pseudo-Anosov maps with an orientable invariant foliation and with small stretch factor.

The code works in three cases: non-orientable surfaces, and orientable surfaces with orientation-preserving and -reserving maps. The orientation-preserving version is essentially identical to a Mathematica program used by Erwan Lanneau and Jean-Luc Thiffeault (On the minimum dilatation of pseudo-Anosov homeromorphisms on surfaces of small genus) to determine the minimal stretch factor in some cases. The other two versions are also modeled on similar ideas. 

There is also a file (`lefschetz.sage`) implementing the Lefschetz tests used by Lanneau and Thiffeault, but with slight modifications so it works for non-orientable surfaces.

## Installation

Just download the files. 

## Requirements

Sage (http://www.sagemath.org) is needed to run the code.

## Usage

The easiest way to run experiments is by running the `demo.py` script from the shell. For example, for listing all polynomials of degree 5 whose largest root can possibly be the minimal stretch factor on a nonorientable surface, run

```
$ sage demo.py 'nonor' 5.
```

It is possible to specify an upper bound for the largest root of the polynomials:

```
$ sage demo.py nonor 3 --largest_root_bound 2
```

For large degree cases, computation case be sped up by separating the cases depending on the trace of the companion matrix of the polynomial. (This trace is the negative of the ``x^{d-1}``-coefficient of a degree `d` polynomial.) For example,

```
$ sage demo.py nonor 3 --first_trace -1
```

only searches for polynomials where the this trace is -1. The command

```
$ sage demo.py reversing 4 --separate_by_first_trace
```

searches for degree 4 polynomials for orientation-reversing pseudo-Anosov maps and invokes several commands of the form

```
$  sage demo.py reversing 4 --first_trace T
```
in order to test various cases at the same time.

More extensive documentation is in demo.py.
