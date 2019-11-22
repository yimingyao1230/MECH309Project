```
MECH 309: Numerical Methods in Mechanical Engineering
```
# MECH 309: Numerical Methods in Mechanical Engineering

# Department of Mechanical Engineering, McGill University

# Final Project: Murman-Cole Scheme for the Transonic Small Disturbance

# Equation

# Due 6pm, 4th. December, 2019

Instructions: Write a report (LateX, LibreOffice, Word...) with comments, equations, figures, and algo-
rithms. Save it as a .pdf file which should neither exceed 2 Mb nor 12 pages. Fontsize should be 11 pt. You
are free to organize the plots in an effective manner. Submit all your Matlab scripts in a .zip archive into
which you will also include your report. The name of the archive should be a concatenation of the three
family names of the members of your team. Submit only one archive per team. Unless stated otherwise, the
use of advanced Matlab commands is prohibited. Equally important is the prohibited use of the symbolic
toolbox in Matlab. All plots must have bothx- andy-axis labels, a legend clearly describing the various
lines, and Figure captions. Plots generated with MS Excel are not acceptable. Provide every step of your
derivations. Please submit a single PDF file.

Grading scheme:Below is the grading scheme used only for the final report to be submitted on April 7th:

1. Comprehension of information and concepts of fundamental engineering sciences (10%)
2. Critical evaluation of the validity and accuracy of solution methods (45%)
3. Appropriate selection of solution techniques and resources (30%)
4. High quality written engineering report (15%)

Solve the transonic small disturbance (TSD) theory over a circular arc airfoil at various Mach numbers using
the Murman-Cole method. The TSD equation is simplified from the two-dimensional full potential equations
and can be written as [
(
1 −M∞^2

## )

```
−(γ+ 1)M∞^2
```
```
φx
U∞
```
## ]

```
φxx+φyy= 0
```
whereγ= 1.4 is the specific heat ratio for air, the gas constant, R = 287.058 J kg−^1 K−^1 , the freestream
static temperature and pressures areT∞= 293 K andp∞= 100 kN/m^2 , and lastly, (x,y)∈[0, 50 .0]^2 spans
the two-dimensional domain, and the following boundary conditions hold:

```
φ(x,y) = 0, ∀(x,50),(j,0),and (j,50) (1)
∂φ
∂n
= 0 ∀(x,0)6∈ 20 ≤x≤ 21 (2)
∂φ
∂n
```
## =U∞

```
dy
dx
```
```
∀(x,0)∈ 20 ≤x≤ 21. (3)
```
The airfoil is defined by the following equation for a circular-arc:

```
y(x) = (t/c)(− 2 x^2 + 82x−840), ∀(x,0)∈ 20 ≤x≤ 21
```
where,t/c= 0.08 is the thickness ratio, andxcis the position of the mid chord. Employ a constant grid
spacing in thexdirection over the airfoil but an exponential or polynomial stretching of the grid along the
xandydirections. You may initialize the flow with,φ= 0. Either specify the freestream Mach number or
velocity as an input parameter. Use the isentropic equations to evaluate Mach, velocity, or pressure in the
domain.
Provide the following in a written report:

1. Solve the TSD using the Murman-Cole method on the computational grid and use at a minimum 60
    points in thex-direction, with at least 60 grid points in they-direction. The first grid spacing adjacent
    to the airfoil in they-direction should bedy=t/c 2. Use the Gauss-Seidel method to solve the equation
    along each column (y-direction).

## 1


```
MECH 309: Numerical Methods in Mechanical Engineering
```
2. Provide a plot of the pressure coefficient along the airfoil surface with the negative pointing upwards.
    Vary the freestream Mach number between [0. 80 , 0 .90] with 0.02 increments. For each case, provide a
    convergence plot of theL∞-norm, surface pressure coefficient as a function ofx, and pressure contour
    forx∈[20,21] andy∈[0,1]. Discuss your findings. A four order reduction in the residual is sufficient
    for the Gauss-Seidel method. Ensure that the residual reaches at least 1× 10 −^4.
3. Vary the grid size, by doubling it in each direction as well as the number of points on the airfoil surface.
    Produce a coarse, medium, and fine grid. Plot the surface coefficient of pressure as a function ofxon
    the airfoil surface for the three grids on the same plot at Mach 0.88. Discuss your findings. Does the
    shock location change.
4. For a select grid size, plot the coefficient of pressure along the airfoil surface for Mach number between
    [0. 80 , 0 .90] with 0.02 increments on the same plot. Discuss your observations. What is the effect of
    increasing the Mach number.

## 2



This is a offline tool, your data stays locally and is not send to any server!
Feedback & Bug Reports