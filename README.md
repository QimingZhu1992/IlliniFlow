# IlliniFlow

A C++ library for Navier-Stokes problem with Arbitrary Lagrange Eulerian method using FEniCS. Variational multiscale method is used in our code. The following article outline our code:
```
@article{zhu2019moving,
  title={A moving-domain CFD solver in FEniCS with applications to tidal turbine simulations in turbulent flows},
  author={Zhu, Qiming and Yan, Jinhui},
  journal={Computers and Mathematics with Applications},
  year={2019},
  publisher={Elsevier}
}
```

### Installing

For local personal computer, following code is used to install FEniCS:  
```bash
sudo apt-get install software-properties-common  
sudo add-apt-repository ppa:fenics-packages/fenics  
sudo apt-get update  
sudo apt-get install --no-install-recommends fenics  
```
For supercomputer, shell code is provided for compiling FEniCS with source code. Be patient with the compiling process.  
```bash
build-pre.sh	is used to compile prerequisites of FEniCS library  
build-dolfin.sh is used to compile FEniCS library  
preuse.sh	is used to find the path of all needed package  
```
Keything to mention:
The preconditioning of our code is based on fieldsplit from PETSC. This part does not fit well with ubuntu 18.04, code works well with ubuntu 16.04. If you want to use ubuntu 18.04, please change your preconditioner to other preconditioner without fieldsplit.

## Running the tests

Two test demo is provided in our code. One case is for flow past stationary sphere with independent rotating mesh. This case should have same result with flow past stationary sphere. The second case is for flow past rotating turbine. The mesh file is too large to include in github. Please generate your own mesh and change the parameter to fit to your case. Please check the details in our paper.
```bash
run.sh    is used to clean the old file and recompile the whole code  
cl.sh     is used to clean the old file  
ns.ufl    is formulation for variational multiscale method of Navier-Stokes problem  
calcd.ufl is formulation for calculating force coefficient  
calQ.ufl  is formulation for calculating Q criterion by L2 Projection  
```
## Authors

Qiming Zhu, phd student in UIUC,          qiming2@illinois.edu  
Jinhui Yan, Assistant Professor in UIUC,  yjh@illinois.edu  

Feel free to send me email, if you have problems with using this code.  

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details  


