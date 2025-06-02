# Error Bounds for the Network Scale-Up Method

[Sergio Díaz-Aranda](https://networks.imdea.org/team/imdea-networks-team/people/sergio-diaz-aranda/), [Juan Marcos Ramirez](https://juanmarcosramirez.github.io/), [Mohit Daga](https://www.kth.se/profile/mdaga), [Jaya Prakash Champati](https://www.uvic.ca/ecs/computerscience/people/faculty/profiles/champati-jaya.php), [Jose Aguilar](https://networks.imdea.org/team/imdea-networks-team/people/jose-aguilar/), [Rosa Lillo](https://halweb.uc3m.es/esp/Personal/personas/rlillo/research.html), [Antonio Fernández-Anta].

## Abstract

Epidemiologists and social scientists have used the Network Scale-Up Method (NSUM) for over thirty years to estimate the size of a hidden sub-population within a social network. This method involves querying a subset of network nodes about the number of their neighbors belonging to the hidden sub-population. In general, NSUM assumes that the social network topology and the hidden sub-population distribution are well-behaved; hence, the NSUM estimate is close to the actual value. However, bounds on NSUM estimation errors have not been analytically proven. This paper provides analytical bounds on the error incurred by the two most popular NSUM estimators. These bounds assume that the queried nodes accurately provide their degree and the number of neighbors belonging to the hidden sub-population. Our key findings are twofold. First, we show that when an adversary designs the network and places the hidden sub-population, then the estimate can be a factor of $\Omega(\sqrt{n})$ off from the real value (in a network with $n$ nodes). Second, we also prove error bounds when the underlying network is randomly generated, showing that a small constant factor can be achieved with high probability using samples of logarithmic size $O(\log n)$. We present improved analytical bounds for Erdős–Rényi and Scale-Free networks. Our theoretical analysis is supported by an extensive set of numerical experiments designed to determine the effect of the sample size on the accuracy of the estimates in both synthetic and real networks.

## Real Networks

Real networks can be downloaded from [Standford Large Network Dataset Collection](https://snap.stanford.edu/data/). Specifically we use [Graph Embedding with Self Clustering: Deezer, February 13 2018](https://snap.stanford.edu/data/gemsec-Deezer.html).


