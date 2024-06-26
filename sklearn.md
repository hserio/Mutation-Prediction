# About the sklearn scaling strategy: Partial Fit API method
This implementation is helpfule in the ML model trainning step when the amount of examples, features, or both propose a timely challenge to process.
This method of "external memory" learning technique fits large input data that cannot fit in the main memory of a computer (RAM)

### Method Overview
1) Stream instances
- yielding in instances (examples) from the data source
2) Feature Extraction
- i.e. transforming arbitrary data into numerical features to be used for ML
- statefull vs stateless vectorizer, if needed
- sklearn.feature_extraction.text.HashingVectorizer for text documents --> we believe this might be especially useful for our sequential data
- DictVectorizer is also a useful representation transformation for training sequence classifiers in Natural Language Processing models that typically work by extracting feature windows around a particular word of interest
3) Incremental algorithm
- although all algorithms cannot learn incrementally (without seeing all data at once), all estimators are valid options (partial_fit)
- learning incrementally from a mini-batch is key to this method which guarantees that at any given moment there is only a stated amount of instances in the main memory
- this step requires some testing for the adequate size of the mini-batch that balances relevancy and memory footprint
- estimator examples: 
  - Classification--> sklearn.linear_model.Perceptron
  - Regression --> sklearn.linear_model.SGDRegressor
  - Clustering --> sklearn.cluster.MiniBatchKMeans
  - Decomposition (feature extraction) --> sklearn.decomposition.MiniBatchDictionaryLearning

Making this pipeline generally more applicable will require a larger dataset that is more represetative of other viruses to train the ML models. 
We acknowledge the challenges that comes with this decision which is why we propose a possible solution to go about this problem to make it more efficient to process.

Ref: https://scikit-learn.org/0.15/modules/scaling_strategies.html 
