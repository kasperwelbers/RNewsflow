CRAN v1.2.1 (Release date: 2019-10-08)
================
  
This update introduces rather big changes, but we have tried to make it as painless as possible. 

* We are working towards dropping support for the tm package style DTMs. RNewsflow started out using the tm package for managing the DTM, but at present the quanteda package is much more suitable. The current version still supports tm, but all functions that use it have been marked as deprecated, and will be removed in the future.
* This major change was a nice opportunity to do some refactoring. The now deprecated functions used a style in which function and argument names used dots to separate words. The new functions instead use underscores. 
* The newsflow_compare function (or newsflow.compare) is now just a thin wrapper for the more versatile compare_documents function. 
* The document comparison functions (and the tcrossprod_sparse backend) now also return margin attributes, such has the number of documents each document has been compared to, and the sum weight. Essentially this is the row/column information that is normally available in an adjacency matrix, but with the special filtering features for comparing documents over time, within groups or using thresholds, this information is lost. (We're currently using this information to experiment with some similarity threshold estimation functions with promising results, so this will probably be added in a next update).
