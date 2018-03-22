# TransformationOptics

This series of MATLAB scripts allows the user to create cloaks for arbitrary shaped objects using the numerical transformation optics method.

## Included Files
1. TransformationOptics.m: Main simulation file, object and cloak are defined here. Script outputs transformed coordinate space along with dielectric tensors for transofming said space.

2. fdder.m: Finite Difference Derivative Matrix function. Given grid, resolution, and boundary condition parameters, the function outputs 4 derivative matrices (D<sub>x</sub>, D<sub>x</sub><sup>2</sup>, D<sub>y</sub>, D<sub>y</sub><sup>2</sup>).

## Sample Ouput
In the case of a simple squared cloak enclosing a triangular object, the output of the script is as follows:

Name        | Output                                                                                                  |
------------|:-------------------------------------------------------------------------------------------------------:|
CLOAK       | ![CLOAK](https://github.com/martinezManuelF/TransformationOptics/blob/master/Graphics/CLOAK.png)        |
ECLOAK      | ![ECLK](https://github.com/martinezManuelF/TransformationOptics/blob/master/Graphics/CLOAK2ECLK.png)    |
OBJECT      | ![OBJ](https://github.com/martinezManuelF/TransformationOptics/blob/master/Graphics/OBJECT.png)         |
EOBJ        | ![EOBJ](https://github.com/martinezManuelF/TransformationOptics/blob/master/Graphics/OBJECT2EOBJ.png)   |
FORCE       | ![FORCE](https://github.com/martinezManuelF/TransformationOptics/blob/master/Graphics/FMAP.png)         |
INTERMEDIATE| ![INT](https://github.com/martinezManuelF/TransformationOptics/blob/master/Graphics/XAYAXFYF.png)       |
TRANSFORMED | ![XBYB](https://github.com/martinezManuelF/TransformationOptics/blob/master/Graphics/XBYB.png)          |
ER          | ![ER](https://github.com/martinezManuelF/TransformationOptics/blob/master/Graphics/ER.png)              |
