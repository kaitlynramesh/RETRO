# RETRO
RETRO (REal Time Reference Ordering) infers pseudotime values that are calibrated with experimental sampling time information to recapitulate the kinetics of cell state transitions for a wide variety of trajectories. It advances upon existing pseudotime inference methodology by integrating sampling time and gene expression data to correctly connect cell states that are sequential in a developmental process. RETRO calculates an ensemble of trajectories and scores each one based on its fit along the projection and consistency with sampling time. The optimal trajectory is smoothed using Bézier spline fitting, and the arc lengths along these curves are calibrated with sampling time to estimate cell-specific pseudotime values. We use RETRO pseudotime to infer kinetic parameters for ODE models of gene regulatory interactions in post-myocardial infarction remodeling and in vitro myelopoiesis. 


## Authors
Kaitlyn Ramesh ramesh.ka@northeastern.edu (maintainer),<br>
Mingyang Lu m.lu@northeastern.edu <br><br>
from the Lu Lab @ Northeastern University (https://lusystemsbio.northeastern.edu/)

## Installation
```
install.packages("remotes") 
remotes::install_github("kaitlynramesh/RETRO")
```

## Tutorial
Coming
