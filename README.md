## The Pilus Carbohydrates Modeling Project
A Research Project by Xiaotong Zuo, Aug. 2015 - Feb. 2016.

---
###Introduction

Acinetobacter baumannii is an environmental bacterium that has recently become a significant source of hospital-acquired infections. Here we modeled the C-terminal glycosylation structures of the major pilin subunit, PilA, from two Acinetobacter strains, __ACICU__ and __M2__, which shows a substantial reduction in surface area available for antibody binding, suggesting a role for pilin glycosylation in immune escape.

---
###Modeling Algorithm

The glycans attached to PilA ACICU and M2 were modeled using [__PyRosetta__](http://www.pyrosetta.org/) protein structural modeling suite. The initial models were generated in the __Discovery Studio Visualizer__. And then the glycans were modeled using the **_FloppyTail_ Algorithm**. The protocol consists two parts. 

1. In the low-resolution part, a random perturbation of the torsion angles was applied.

2. The structures were then refined in the high-resolution part by applying a more precise perturbation of the torsion angles, the side-chain packing and the minimization. 

With this protocol, 6000 structures of PilA ACICU and 6000 structures of PilA M2 were generated.

(revised on 02/16/2016)
