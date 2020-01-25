---
title: "AD Target Gene Report"
author: "Jake Gockley"
date: "11/11/2019"
header-includes:
  - \usepackage{multicol}
  - \newcommand{\btwocol}{\begin{multicols}{2}}
  - \newcommand{\etwocol}{\end{multicols}}
output:
  html_document: default
  pdf_document: default
editor_options: 
  chunk_output_type: console
---




```
## [1] "/home/jgockley/AD_TargetRank"
```

Importing genelist genelists/EmoryTargets.txt as GNAME

```
## Cache found
```


```
## Skipping IGAP based Model and Gene-Ranking: User specified module 1 not included in configs/EmoryTargets.yaml
```


```
## Skipping Differential Expression based Model and Gene-Ranking: User specified module 2 not included in configs/EmoryTargets.yaml
```


```
## Skipping Differential Proteomic Expression based Model and Gene-Ranking: User specified module 3 not included in configs/EmoryTargets.yaml
```


```
## Skipping eQTL based Model and Gene-Ranking: User specified module 4 not included in configs/EmoryTargets.yaml
```


```
## Skipping Coexpression Module based Model and Gene-Clustering: User specified module 5 not included in configs/EmoryTargets.yaml
```


```
## Skipping Coexpression Module based Model and Gene-Clustering: User specified module 6 not included in configs/EmoryTargets.yaml
```


```
## Modeling and ranking genelist genelists/EmoryTargets.txt by Differential Expression evidence
```

```
## Loading required package: foreach
```

```
## 
## Attaching package: 'foreach'
```

```
## The following objects are masked from 'package:purrr':
## 
##     accumulate, when
```

```
## Loading required package: iterators
```

```
## status: getting commit information about: https://api.github.com/repos/jgockley62/AD_TargetRank/git/commits/89c4f0d208ee65c5d21af8486806eaeec75dee8d
## status: getting information about the commit tree
## For each tissue specified in configs/EmoryTargets.yaml the coexpression of genes specified in genelists/EmoryTargets.txt is given. This is measured as a percentage of 1000 pemutations in which the gene was a significant predictor of expression of its comparison gene by a linear model.
```

```
## Cache found
```

```
## Cache found
## Cache found
## Cache found
## Cache found
## Cache found
## Cache found
```

${image?fileName=Module7%2D1%2Epng&align=none&scale=100}
${image?fileName=Module7%2D2%2Epng&align=none&scale=100}
${image?fileName=Module7%2D3%2Epng&align=none&scale=100}
${image?fileName=Module7%2D4%2Epng&align=none&scale=100}
${image?fileName=Module7%2D5%2Epng&align=none&scale=100}
${image?fileName=Module7%2D6%2Epng&align=none&scale=100}

```
## All Listed Genes Expressed in CBE
```

```
## All Listed Genes Expressed in DLPFC
```

```
## All Listed Genes Expressed in FP
```

```
## All Listed Genes Expressed in IFG
```

```
## All Listed Genes Expressed in PHG
```

```
## All Listed Genes Expressed in STG
```

```
## All Listed Genes Expressed in TCX
```
