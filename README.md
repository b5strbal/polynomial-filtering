# polynomial-filtering

Code for searching for polynomials that can possibly arise as characteristic polynomials of the action on homology for pseudo-Anosov maps with an orientable invariant foliation and with small stretch factor.

The code works in three cases: non-orientable surfaces, and orientable surfaces with orientation-preserving and -reserving maps. The orientation-preserving version is essentially identical to a Mathematica program used by Erwan Lanneau and Jean-Luc Thiffeault (On the minimum dilatation of pseudo-Anosov homeromorphisms on surfaces of small genus) to determine the minimal stretch factor in some cases. The other two versions are also modeled on similar ideas. 

There is also a file (`lefschetz.sage`) implementing the Lefschetz tests used by Lanneau and Thiffeault, but with slight modifications so it works for non-orientable surfaces.

## Installation

Download the files. 

## Usage

Sage is needed to run the code. The easiest thing to do is to run the `demo.py` script from the shell, such as

```
$ sage demo.py 'nonor' 5
```

More extensive documentation is in demo.py.
