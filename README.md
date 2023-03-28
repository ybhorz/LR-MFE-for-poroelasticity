# Reproducibility Information 

*A Robust and Partially Superconvergent Finite Element Method for Poroelasticity Problems*

**Authors** 	
- Bohan Yang <ybh.orz@gmail.com>  
- Hongpeng Li <lihongpeng_sd@163.com>

**Corresponding Author**  
- Hongxing Rui <hxrui@sdu.edu.cn>

## Numerical Experiment

### Experiment 1

#### Rate of Error

`test1_1.m`

`test1_result\eRate_output.m`

**Parameter Setting:**

```matlab
domn = [0,1,0,1]; nSubs = [8,16,32,64,128];
EgDire = [2,3,4]; EgNeum = 1;
lmd = 1; mu = 1; alph = 1; c0 = 0; k = 1;
vu = 0.1; vp = [];
test = 'test1'; method = []; index = [];

ux = 'vu*sin(pi*x)*(1-cos(pi*y))';
uy = 'vu*(1-cos(2*pi*x))*sin(pi/2*y)';
p = 'vp*cos(pi*x)*cos(pi*y)';
```

##### case 1

```matlab
vp = 1; method = 'LR'; index = '1';
```

result: `test1_result\eRate_LR_1.mat`

formated output: `test1_result\eRate_LR_1.txt` 

##### case 2

```matlab
vp = 1e4; method = 'LR'; index = '2';
```

result: `test1_result\eRate_LR_2.mat`

formated output: `test1_result\eRate_LR_2.txt` 

##### case 3

```matlab
vp = 1e4; method = 'BR'; index = '2';
```

result: `test1_result\eRate_BR_2.mat`

formated output: `test1_result\eRate_BR_2.txt` 

#### Evolution of Error with Pressure

`test1_2.m`

`test1_result\pRob_output.m`

**Parameter Setting:**

```matlab
domn = [0,1,0,1]; nSub = 100;
EgDire = [2,3,4]; EgNeum = 1;
lmd = 1; mu = 1; alph = 1; c0 = 0; k = 1;
vu = 0.1; vps = 10.^(-1:6);
test = 'test1'; method = []; index = '1';

ux = sprintf('%.1g*sin(pi*x)*(1-cos(pi*y))',vu);
uy = sprintf('%.1g*(1-cos(2*pi*x))*sin(pi/2*y)',vu);
p = 'vp*cos(pi*x)*cos(pi*y)';
```

##### LR-MFE

```matlab
method = 'LR';
```

result: `test1_result\pRob_LR_1.mat`

figure: `test1_result\pRob_LR_1.fig`

##### BR-MFE

```matlab
method = 'BR';	
```

result: `test1_result\pRob_BR_1.mat`

figure: `test1_result\pRob_BR_1.fig`

### Experiment 2

#### Rate of Error

`test2_1.m`

`test2_result\eRate_output.m`

**Parameter Setting:**

```matlab
domn = [0,1,0,1]; nSubs = [8,16,32,64,128];
EgDire = [2,3,4]; EgNeum = 1;
lmd = 1e6; mu = 1; alph = 1; c0 = 0; k = 1;
test = 'test2'; method = []; index = '1';

ux = '0.1*(sin(2*pi*y)*(-1+cos(2*pi*x))+1/lmd*sin(pi*x)*sin(pi*y))';
uy = '0.1*(sin(2*pi*x)*(1-cos(2*pi*y))+1/lmd*sin(pi*x)*sin(pi*y))';
p = 'cos(pi*x)*(1-cos(pi*y))';
```

##### LR-MFE

```matlab
method = 'LR';
```

result: `test2_result\eRate_LR_1.mat`

formated output: `test2_result\eRate_LR_1.txt`

##### CG-MFE

```matlab
method = 'P1';
```

result: `test2_result\eRate_P1_1.mat`

formated output: `test2_result\eRate_P1_1.txt`

#### Evolution of Error with Lambda

`test2_2.m`

`test2_result\lamRob_output.m`

**Parameter Setting:**

```matlab
domn = [0,1,0,1]; nSub = 32;
EgDire = [2,3,4]; EgNeum = 1;
lmds = 10.^(0:7); 
mu = 1; alph = 1; c0 = 0; k = 1;
test = 'test2'; method = []; index = '1';

ux = '0.1*(sin(2*pi*y)*(-1+cos(2*pi*x))+1/lmd*sin(pi*x)*sin(pi*y))';
uy = '0.1*(sin(2*pi*x)*(1-cos(2*pi*y))+1/lmd*sin(pi*x)*sin(pi*y))';
p = 'cos(pi*x)*(1-cos(pi*y))';
```

##### LR-MFE

```matlab
method = 'LR';
```

result: `test2_result\lmdRob_LR_1.mat`

figure: `test2_result\lmdRob_LR_1.fig`

##### CG-MFE

```matlab
mathod = 'P1';
```

result: `test2_result\lmdRob_P1_1.mat`

figure: `test2_result\lmdRob_P1_1.fig`

### Experiment 3

#### cantilever bracket problem

##### LR-MFE

`test3_1_LR.m`

figure: `test3_1_result\deformation_LR.fig`, `test3_1_result\pressure_LR.fig`

##### CG-MFE

`test3_1_P1.m`

figure: `test3_1_result\deformation_P1.fig`, `test3_1_result\pressure_P1.fig`

#### Barry and Mercer problem

##### LR-MFE

`test3_2_LR.m`

figure: `test3_2_result\pressure_LR.fig`

##### CG-MFE

`test3_2_P1.m`

figure: `test3_2_result\pressure_P1.fig`

### Experiment 4

`test4.m`

`test4_result\eRate_output.m`

#### Case 1

**Parameter Setting:**

```matlab
domn = [0,1,0,1]; nSubs = [8,16,32,64,128];
EgDire = [2,3,4]; EgNeum = 1;
lmd = 1; mu = 1; alph = 1; c0 = 0; k = 1;
vu = 0.1; vp = 1e4;
test = 'test4'; method = 'LR'; post = []; index = '1'; 

ux = 'vu*sin(pi*x)*(1-cos(pi*y))';
uy = 'vu*(1-cos(2*pi*x))*sin(pi/2*y)';
p = 'vp*cos(pi*x)*cos(pi*y)';
```

##### Post-Processing 1

```matlab
post = 1;
```

result: `test4_result\eRate_LR_post1_1.mat`

formated output: `test4_result\eRate_LR_post1_1.txt`

##### Port-Processing 2

```
post = 2;
```

result: `test4_result\eRate_LR_post2_1.mat`

formated output: `test4_result\eRate_LR_post2_1.txt`

#### Case 2

**Parameter Setting:**

```matlab
domn = [0,1,0,1]; nSubs = [8,16,32,64,128];
EgDire = [2,3,4]; EgNeum = 1;
lmd = 1; mu = 1; alph = 1; c0 = 0; k = 1e-6;
vu = 0.1; vp = 1e4;
test = 'test4'; method = 'LR'; post = []; index = 2; 

ux = 'vu*sin(pi*x)*(1-cos(pi*y))';
uy = 'vu*(1-cos(2*pi*x))*sin(pi/2*y)';
p = 'vp*cos(pi*x)*cos(pi*y)';
```

##### LR-MFE

```matlab
post = [];
```

result: `test4_result\eRate_LR_2.mat`

formated output: `test4_result\eRate_LR_2.txt`

##### Post-Processing 2

```matlab
post = 2;
```

result: `test4_result\eRate_LR_post2_2.mat`

formated output: `test4_result\eRate_LR_post2_2.txt`

## Appendix

`FEM_Package`:  basic components for the finite element method.

`mesh2D_[0,1,0,1]`, `mesh2D_[-1,1,-1,1]`: mesh informations (already generated and stored). Not used.

`Biot_LR.m`, `Biot_BR.m`, `Biot_P1.m`:  Solve the Biot's equation by using LR-MFE, BR-MFE and CG-MFE methods.

**Parameters:**  
- `domn`: `[xmin,xmax,ymin,ymax]`  
- `nSub`:  Controls the size of subdivisions  
- `EgDire, EgNeum`: Dirichlet boundary, Neumann boundary. Integer 1-4 for each side  
	- 1: top side  
	- 2: right side  
	- 3: bottom side  
	- 4: left side  
- `lmd`: lame constant  
- `mu`: lame constant  
- `alph`: biot-willis parameter  
- `c0`: constrained specific storage coefficient  
- `k`: dynamic permeability  
- `ux`: displacement in x-direction  
- `uy`: displacement in y-direction  
- `p`: pressure  
- `post`: post-processing. Integer 1, 2 for post-processing 1, 2. Not used in `Biot_BR.m, Biot_P1.m`.  
