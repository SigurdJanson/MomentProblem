# MomentProblem

Imagine a situation when you want a probability distribution with a specified series of moments. Such a thing can happen e.g. when we investigate statistical properties in Monte Carlo studies. The idea behind this repository is to collect methods to achieve that, i.e. functions that:

* Approximate a distribution from specified moments
  * Identify starting points for approxmation
  * Identify optimal parameters using specified optimization algorithms
* Evaluate the appropriateness of the approximation
* Visualise the results

## Architecture

It's quite simple. One base class provides the full data structure and the basic functionality shared by all approximations. Some of the functionality required to accomplish an approximation is hidden inside additional functions. Technically, they are not part of the class.

![Class Diagram](http://www.plantuml.com/plantuml/proxy?src=https://raw.githubusercontent.com/SigurdJanson/MomentProblem/master/uml/bymomentpdf.puml)

