# Error Bounds for the Network Scale-Up Method

[Sergio Díaz-Aranda](https://networks.imdea.org/team/imdea-networks-team/people/sergio-diaz-aranda/), [Juan Marcos Ramirez](https://juanmarcosramirez.github.io/), [Mohit Daga](https://www.kth.se/profile/mdaga), [Jaya Prakash Champati](https://www.uvic.ca/ecs/computerscience/people/faculty/profiles/champati-jaya.php), [Jose Aguilar](https://networks.imdea.org/team/imdea-networks-team/people/jose-aguilar/), [Rosa Lillo](https://halweb.uc3m.es/esp/Personal/personas/rlillo/research.html), [Antonio Fernández-Anta](https://software.imdea.org/es/people/antonio.fernandez/).

## Abstract

Epidemiologists and social scientists have used the Network Scale-Up Method (NSUM) for over thirty years to estimate the size of a hidden sub-population within a social network. This method involves querying a subset of network nodes about the number of their neighbors belonging to the hidden sub-population. In general, NSUM assumes that the social network topology and the hidden sub-population distribution are well-behaved; hence, the NSUM estimate is close to the actual value. However, bounds on NSUM estimation errors have not been analytically proven. This paper provides analytical bounds on the error incurred by the two most popular NSUM estimators. These bounds assume that the queried nodes accurately provide their degree and the number of neighbors belonging to the hidden sub-population. Our key findings are twofold. First, we show that when an adversary designs the network and places the hidden sub-population, then the estimate can be a factor of $\Omega(\sqrt{n})$ off from the real value (in a network with $n$ nodes). Second, we also prove error bounds when the underlying network is randomly generated, showing that a small constant factor can be achieved with high probability using samples of logarithmic size $O(\log n)$. We present improved analytical bounds for Erdős–Rényi and Scale-Free networks. Our theoretical analysis is supported by an extensive set of numerical experiments designed to determine the effect of the sample size on the accuracy of the estimates in both synthetic and real networks.

## Real Networks

Real networks can be downloaded from [Standford Large Network Dataset Collection](https://snap.stanford.edu/data/). Specifically we use [Graph Embedding with Self Clustering: Deezer, February 13 2018](https://snap.stanford.edu/data/gemsec-Deezer.html).

[![DOI](https://zenodo.org/badge/994663840.svg)](https://doi.org/10.5281/zenodo.15575415)


## Bibtex

```
@misc{díazaranda2024errorboundsnetworkscaleup,
      title={Error Bounds for the Network Scale-Up Method}, 
      author={Sergio Díaz-Aranda and Juan Marcos Ramírez and Mohit Daga and Jaya Prakash Champati and José Aguilar and Rosa Elvira Lillo and Antonio Fernández Anta},
      year={2024},
      eprint={2407.10640},
      archivePrefix={arXiv},
      primaryClass={cs.DC},
      url={https://arxiv.org/abs/2407.10640}, 
}
```

## Platform

The code has been executed in a Linux environment, specifically on the Ubuntu 24.04 Operating System.

## Acknowledgments

This paper has been funded by project PID2022-140560OB-I00 (DRONAC) funded by MICIU/AEI /10.13039/501100011033 and ERDF, EU. This research is part of the I+D+i projects PID2022-137243OB-I00 funded by MCIN/AEI/10.13039/501100011033 and European Union NextGenerationEU/PRTR and the project CuidaNSUM of the Instituto de las Mujeres. This initiative has also been partially carried out within the framework of the Recovery, Transformation and Resilience Plan funds, financed by the European Union (Next Generation) through the grant ANTICIPA (INCIBE) and the ENIA 2022 Chairs for the creation of university-industry chairs in AI-AImpulsa: UC3M-Universia. The work of Sergio Díaz-Aranda has been funded by \textit{Comunidad de Madrid} predoctoral grant PIPF-2022/COM-24467.

## License

This code package is licensed under the GNU GENERAL PUBLIC LICENSE (version 3) - see the [LICENSE](LICENSE) file for details


## Author

Juan Marcos Ramírez Rondón. Postdoctoral Researcher. [IMDEA Networks Institute](https://networks.imdea.org/es/). Leganés, 28918, Spain. 


### Contact

[Juan Marcos Ramirez](juan.ramirez@imdea.org)

## Date

June 2nd, 2025
