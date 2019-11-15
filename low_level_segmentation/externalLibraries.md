# External Libraries

## Linear Programming

Follow instructions below to install GLPK on windows. Instead of using 2.8 (not compiled for win64) use version **2.11**

http://asbidyarthy.blogspot.com/2012/09/how-to-install-glpk-on-matlab-in-windows.html

## Graph Cuts

Download build from:
https://es.mathworks.com/matlabcentral/fileexchange/21310-maxflow

Follow instructions on Readme.txt:

1. Unzip code from Mathworks into
`<lib_folder>`

2. Download src code from: https://vision.cs.uwaterloo.ca/code/ and unzip it in: 
`<lib_folder>/maxflow-v3.0` (see version in Readme)

3. In <lib_folder>: `run make.m`
```
>> run make.m
Building with 'Microsoft Visual C++ 2015'.
MEX completed successfully.
```

To test the correct installation:
```
>> run test1.m
flow = 3
labels =
  2×1 int32 column vector

   1
   1
```