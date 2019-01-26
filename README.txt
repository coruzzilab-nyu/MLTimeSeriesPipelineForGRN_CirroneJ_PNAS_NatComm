% Copyright (C) 2019 Jacopo Cirrone

Machine learning pipeline for Modeling Root and Shoot time-series data in Arabidopsis to infer, validate and analyze tissue-specific Gene Regulatory Networks, which includes:
- Wrapper Script to parse the data sets to Dynamic Factor Graph(DFG) input format.
- DFG modeling and selecting the combination of hyper-parameters according to the best model accuracy.
- Validate the Gene Regulatory Network of the best model given a gold standard set of edges.
- Computing Area under Precision-Recall Curve.