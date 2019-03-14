# Parareal for Problems with Discontinuous Sources

This repository contains the main source code of the algorithm documented in the paper "A New Parareal Algorithm for Problems with Discontinuous Sources", see https://arxiv.org/abs/1803.05503.

## Getting Started

* [GetDP](http://getdp.info) - Finite element solver
* [Octave](https://www.octave.org) - Interpreter
* [Parallel](https://octave.sourceforge.io/parallel/) - Parallel execution package for Octave

### Prerequisites

The scripts can be executed in Matlab or Octave. However, Octave with an installed parallel execution package is recommended to obtain acceleration due to parallelization. In Octave, the parallel package can be automatically installed by `pkg install -forge struct` and `pkg install -forge parallel`. The execution of the numerical example ("im_3kW") requires the installation of the finite element package GetDP, e.g. in `~/bin`.

## Running the example

Execute in Matlab/Octave the m-script `runme_para_im3kw`.

## Authors

* I. Kulchytska-Ruchka, S. Schöps, Technische Universitaet Darmstadt
* I. Niyonzima, Université Grenoble Alpes
* M. J. Gander, Université de Genève

## License

This project is licensed under the terms of the GNU General Public License (GPL) (version 2 or later).

## Acknowledgments

* the authors like to thank J. Corno for several code beautifications
* the software GetDP was mainly written by P. Dular and C. Geuzaine, University of Liege
* the model im_3kW is modified but originally part of the GetDP distribution. The authors are J. Gyselinck and R.V. Sabariego
