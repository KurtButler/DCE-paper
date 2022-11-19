# Using Gaussian Process Regression to Measure the Strength of Causation
In this repo, we provide code to reproduce the figures from our paper "[A Differential Measure of the Strength of Causation](https://ieeexplore.ieee.org/abstract/document/9924562)," published to IEEE Signal Processing Letters. 

> **Abstract:** The ability to quantify the strength of an interaction between random variables is important in many applications such as medicine and environmental science. We present the problem of measuring the strength of a causal interaction, starting from the linear perspective and generalizing to a nonlinear measure of causal influence, using a differential calculus approach. The proposed measure of causal strength is interpretable and may be estimated efficiently using Gaussian process regression. We validate our estimation approach on several synthesized data sets, considering both static variables and time series.


This code creates four figures when executed, corresponding to figures from the paper:
- Figure 1: A demonstration of Simpson's paradox, which explains why a measurement of causal strength depends critically on the model being assumed.
- Figure 2: An example demonstrating the differential causal effect measure (DCE) in a linear system.
- Figure 3: An example demonstrating the DCE in a nonlinear system.
- Figure 4: The DCE applied to a time series problem.

## Instructions
To generate all figures (as .png files), you just need to run `main.m`. Also all the subfiles can be run independently.

We wrote the code with Matlab 2022a.
We did not optimize the code for efficiency; However, the code takes roughly 30 seconds to generate all of the figures. 

## Citation
If you use this project, please cite:
```
@article{butler2022differential,
  title={A Differential Measure of the Strength of Causation},
  author={Butler, Kurt and Feng, Guanchao and Djuri{\'c}, Petar M},
  journal={IEEE Signal Processing Letters},
  year={2022},
  publisher={IEEE},
  url={https://ieeexplore.ieee.org/abstract/document/9924562}
}
```
## Code
You can download the code using git:
```
git clone https://github.com/newyorkfishmarket/DCE-paper.git
```
This project is also available as a capsule on [Code Ocean](https://codeocean.com/capsule/5974235/tree/v1):
```
git clone https://git.codeocean.com/capsule-0604411.git
```
