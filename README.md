# LieAlgebraRankCondition.jl


`LieAlgebraRankCondition` is a small tool (implemented in Julia Language) to determine a subset of L. I. elements of the set of matrices generated by Lie products obtained from two given arrays `A` and `B`. If you are interested in knowing more about the Lie Algebra Rank Condition, we recommend [1]. The process of installation and use are simple and are described below. 

## Installation

First, you need to have the `Julia language` (v1.0.5 or later) installed. [Click here](https://julialang.org/) to get and install. Now, in Julia REPL, you need to acess the package manager. It can be done by pressing the `]` key and typing

```julia
Pkg> add https://github.com/evcastelani/LieAlgebraRankCondition.jl
```

## Usage

First, you need to load the package. It can be done typing

```julia
julia> using LieAlgebraRankCondition
```

There are two importants functions implemented. They are: `isli` and `rankcondition`. The first one is used to determine if a set of arrays is L.I or not. For example,

```julia
julia> C = [[0 2 0 -1;2 0 2 0;0 2 0 2;-1 0 2 0.0],[1.0 0 0 0;0 4 0 0; 0 0 -2 0;0 0 0 -3]]

julia> isli(C)
```
returns `true`. 

Function `rankcondition` provides, if it exists, a L.I set of matrices formed from convenient Lie products from given matrices. For example

```julia
julia> A = [0 2 0 -1;2 0 2 0;0 2 0 2;-1 0 2 0.0]

julia> B=[1.0 0 0 0;0 4 0 0; 0 0 -2 0;0 0 0 -3]

julia> rankcondition(A,B)

```
if there is a L. I set as mentioned before, the output given by the function describes the set and the operations (`path`) used to get the mentioned set. For example, using `A` and `B` given before, the output will be:


```julia
🎈 A set of L.I arrays with 15 elements was found.

 👉 Here are the set and the path:  

 ↪  Element 1: A = 
4×4 Array{Float64,2}:
  1.0  1.0   0.0   0.0
 -1.0  1.0   0.0   0.0
  0.0  0.0  -1.0   0.5
  0.0  0.0  -0.5  -1.0

 ↪  Element 2: B = 
4×4 Array{Float64,2}:
 0.789769  0.404268   0.315184  0.132109
 0.828743  0.293395   0.330319  0.878407
 0.715575  0.0538233  0.413952  0.212379
 0.936372  0.614032   0.240568  0.304594

 ↪  Element 3: [A,B] = 
4×4 Array{Float64,2}:
  1.23301   -0.496374   1.02674     0.985033 
 -0.496374  -1.23301    0.784658    1.45955  
 -0.909142  -0.516206   0.226473   -0.0546791
 -1.6165    -2.19135   -0.0546791  -0.226473 

 ↪  Element 4: [A,[A,B]] = 
4×4 Array{Float64,2}:
 -0.992747  -2.46602    3.33066     2.91624  
 -2.46602    0.992747   1.27235     1.54173  
  0.493827   0.84588   -0.0546791  -0.226473 
  1.49622    6.2573    -0.226473    0.0546791

 ↪  Element 5: [A,[A,[A,B]]] = 
4×4 Array{Float64,2}:
 -4.93204     1.98549    9.39178     5.70888  
  1.98549     4.93204   -0.0150999  -0.468959 
  0.606337    0.943062  -0.226473    0.0546791
  3.01794   -14.4338     0.0546791   0.226473 

 ↪  Element 6: [A,[A,[A,[A,B]]]] = 
4×4 Array{Float64,2}:
   3.97099   9.86409  21.6229      6.25291  
   9.86409  -3.97099  -9.65646    -6.63925  
   1.23936  -9.70934   0.0546791   0.226473 
 -20.7728   25.378     0.226473   -0.0546791

 ↪  Element 7: [A,[A,[A,[A,[A,B]]]]] = 
4×4 Array{Float64,2}:
  19.7282    -7.94198   36.7158      -4.94488  
  -7.94198  -19.7282   -44.2555     -14.7032   
 -22.5745    30.8683     0.226473    -0.0546791
  66.304    -25.1286    -0.0546791   -0.226473 

 ↪  Element 8: [A,[A,[A,[A,[A,[A,B]]]]]] = 
4×4 Array{Float64,2}:
  -15.884   -39.4564    26.7037     -42.9508   
  -39.4564   15.884   -132.578       -2.33373  
  109.169   -51.7265    -0.0546791   -0.226473 
 -146.449   -31.4809    -0.226473     0.0546791

 ↪  Element 9: [A,[A,[A,[A,[A,[A,[A,B]]]]]]] = 
4×4 Array{Float64,2}:
  -78.9127   31.7679  -100.646      -101.587    
   31.7679   78.9127  -293.027       104.573    
 -343.29    -21.4566    -0.226473      0.0546791
  206.833   235.274      0.0546791     0.226473 

 ↪  Element 10: [A,[A,[A,[A,[A,[A,[A,[A,B]]]]]]]] = 
4×4 Array{Float64,2}:
  63.5358    157.825   -545.113      -48.2788   
 157.825     -63.5358  -433.122      457.246    
 768.539     503.84       0.0546791    0.226473 
  -6.74703  -666.654      0.226473    -0.0546791

 ↪  Element 11: [A,[A,[A,[A,[A,[A,[A,[A,[A,B]]]]]]]]] = 
4×4 Array{Float64,2}:
   315.651   -127.072  -1547.49        633.245    
  -127.072   -315.651    -92.5073     1179.33     
 -1036.61   -2109.55       0.226473     -0.0546791
 -1037.43    1088.13      -0.0546791    -0.226473 

 ↪  Element 12: [A,[A,[A,[A,[A,[A,[A,[A,[A,[A,B]]]]]]]]]] = 
4×4 Array{Float64,2}:
 -254.143  -631.302   -2870.86       3219.57     
 -631.302   254.143    1952.14       1771.67     
 -555.037  5799.77       -0.0546791    -0.226473 
 3681.3     -84.0662     -0.226473      0.0546791

 ↪  Element 13: [A,[A,[A,[A,[A,[A,[A,[A,[A,[A,[A,B]]]]]]]]]]] = 
4×4 Array{Float64,2}:
 -1262.6       508.287  -2179.8        9646.23     
   508.287    1262.6     7660.98       -652.291    
  8750.5    -11086.5       -0.226473      0.0546791
 -7169.15    -6413.05       0.0546791     0.226473 

 ↪  Element 14: [B,[A,[A,[A,[A,[A,[A,[A,[A,[A,[A,[A,B]]]]]]]]]]]] = 
4×4 Array{Float64,2}:
  -5877.5   -8874.13    187.301   4599.72
 -10623.2   -9091.43  -3150.45    5198.51
   3500.71  -5804.45   -243.39   15450.1 
  10027.7    1410.46   7040.82   15212.3 

 ↪  Element 15: [B,[B,[A,[A,[A,[A,[A,[A,[A,[A,[A,[A,[A,B]]]]]]]]]]]]] = 
4×4 Array{Float64,2}:
   1046.79  -10181.8    3327.46  19744.4 
  15288.0    -6760.59  11740.1   32278.3 
 -13446.0   -18128.6   -1443.05  13179.5 
 -36501.1   -29046.2   -9873.69   7156.86

```

Alternatively, you can store the answer in a variable and access each field separately. For example,

```julia
julia> s = rankcondition(A,B)

julia> s.elements

julia> s.path
``` 

Finally, a secondary function implemented is the `lieprod`, which returns the Lie product between `A` and `B`.  


## Credits

1. Alexandre Jose Santana
1. Eduardo Celso Viscovini
1. Emerson Vitor Castelani
1. João Cossich

## Reference

[1] Do Rocio, O. G., & Santana, A. J. (2006). Invariant cones and convex sets for bilinear control systems and parabolic type of semigroups. Journal of dynamical and control systems, 12(3), 419-432.
