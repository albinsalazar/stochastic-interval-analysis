# stochastic-interval-analysis

The program is used to coarse-grain stochastic chemical reaction networks using overlapping intervals with probabilities.

## Name
Case study: 2-d interval analysis of competition for resources.

## Description
The goal of the analysis is to coarse-grain with intervals a 2-dimensional reaction network and label macro-transitions with probabilities. Each interval represents a range of values for each chemical species, and the product of intervals are referred to as states. The source of a macro-transition is a point when a coordinate enters a new interval, and the target (a goal) is the proceeding interval. Thus, once a sequence of chemical transitions enter a new interval, the probability to reach the next interval is computed using a recurrence relation for the reference coordinate. 


## Visuals
![alt text](https://github.com/albinsalazar/stochastic-interval-analysis/blob/main/2-d-random-walk.pdf)

## Installation
The code can be easily ran locally.


## Roadmap
The next version of the analysis will include coarse-grain approximations for multi-dimensional systems using an ILP forumaltion. 

## Project status
Last updated: 2 May, 2022. The project is under active development.
