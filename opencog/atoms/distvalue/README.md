Distributional Values
=====================

An implementation of Distributional Values that might eventually replace all current TruthValues.

Problem
-------

One of the problems this is trying to solve is the combination of fuzzy and propabalistc logic. As both are used within OpenCog it is important to be able to clearly represent and distinguise them.

Representation
--------------

To solve these issues we are using a Dirichlet Distribution where each category is a bin in a Histogram over the range of truthness (0-1).
There is no fixed number of bins, they can have different widths and are allowed to overlapp. So when trying to increase the count of a bin that doesn't exist yet it is created.

Each bin is defined by a Product of Intervals and has a count associated with it. This allows us to represent a Joint Probabilty Distribution over n Concepts.

Another improvment is the addtion of an explicitly conditional DV which can be used for InheritanceLinks and the like. Similar to the basic DV it consists of a Set of bins with the differenc that we associated a full DV to each of these.

For many calculations involving a DV and another DV or an interval we need to calculate how different bins are overlapping.
For example when calculating the count of an Interval multiple overlapping bins might contribute to this count.
When we consider the Intervals to be continuous uniform distributions this calculation is equivalent to caluclatin the conditional probabilties of the bins given an Interval and summing their weighted counts accordingly.

Because of this it might be valuable to consider using other distributions for the Intervals.
