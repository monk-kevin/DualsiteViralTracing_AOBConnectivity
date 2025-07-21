# DualsiteViralTracing_AOBConnectivity
These files represent the analysis pipeline I have created to define efferent connectivity patterns between the accessory olfactory bulb (AOB) and downstream limbic areas. Below is a description of the data and the steps in the analysis pipeline. Raw data can be made available upon reasonable request.

## Description of data
The data analyzed are counts of labeled cells using ImageJ/FIJI. Cells were labeled through viral injection of target areas with retrograde viruses carrying one of two fluorophores. Cells were counted as having one, both, or neither fluorophore present. With these data we were able to describe the overall connectivity patterns between the AOB and downstream areas, sex-specific patterns of neural connectivity, and estimates of true co-projection probabilities using ground-truth, probabilistic modeling.

## Analytical pipeline
### Region-specific connectivity patterns
The first set of question relates to describing AOB neurons that project to downstream areas. Namely, we first describe the proportion of neuron that project to specific regions. Moreover, by normalizing spatial location across individual samples, we were able to establish a map of the AOB and determine whether there are distinct anatomical patterns within the AOB as it relates to downstream connectivity. Significant differences were established with non-parametric and parametric statistical tests.

### Sex-specific connectivity patterns
We then determined whether connecitivity patterns differed between male and female mice. To this end we compared distributions between sexes and performed statistical analyses to define differences.

### Co-projection probability analyses
Through the use of bootstrap statistics, we were able to define whether the observed pattern of co-labeling is consistent or inconsistent with a system wherein neurons have structured connectivity patterns with downstream areas.
Moreover, through the use of probabilistic modeling, we created ground-truth models with which we compared our observations to estimate the true co-projection probability of neurons projecting to multiple areas.
