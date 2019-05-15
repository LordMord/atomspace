Distributional Values
=====================

An implementation of Distributional Values that might eventually replace all current TruthValues.

Problem
-------

One of the problems this is trying to solve is the combination of fuzzy and propabalistc logic. As both are used within OpenCog it is important to be able to clearly represent and distinguise them.

Representation
--------------

To solve these issues we are using a Dirichlet Distribution where each category coresponds to the bin in a Histogram. As we need to do complex operations in potentially high dimensional spaces with these Histograms we only store the center of a bin instead of a complete Interval. These Mutlidemensional bins are used to store a Joint Probability Distribution over n Concepts in the Histogram.

For storing expclicity conditional Distribution Values we use a second order Histogram where each bin contains not a count but instead a first order Histogram.  These can be used for InheritanceLinks and the like.

To efficently work with these Histograms they utilize a Cover Tree for storage. 
