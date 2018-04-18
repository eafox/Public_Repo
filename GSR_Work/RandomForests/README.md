# Random and Gradient Forests

#### Background
Random and gradient forests are a system of decision trees. This method is useful as an alternative to a glm for parsing the effect of multiple parameters on a response variable.  
This folder includes examples for the following packages:  
* [tree](https://cran.r-project.org/web/packages/tree/index.html)
* [randomForest](https://cran.r-project.org/web/packages/randomForest/index.html)
* [gradientForest](http://gradientforest.r-forge.r-project.org/)
* [LFMM](http://membres-timc.imag.fr/Olivier.Francois/lfmm/software.htm)

## Folders and Contents
* __Tutorials__  
Includes data and code for tutorials (sourced online, copied, and appropriately cited) for each of the packages mentioned above. ALL TUTORIALS IN THIS FOLDER ARE COPIED FROM OTHER SOURCES ONLINE. Have included the scripts I used to run each of them to check they worked.  
Also includes a script with information on installing the various packages.  
* __eDNA_Practice__  
_Data Cleaning and Mining_  
Combine_DataSets.R - Take all 16 primer result ASV files and combine them into one data set with columns representing unique ASVs and each row representing a sample.  
Collapse_byTaxonomy.R - Can take the combined ASV file from the team drive and collapse the data down to a specific taxonomic level (phylum, order, class, etc.).   
Get_Layers.R - Script showing how some of the data layers were downloaded and how values were extracted.   
_Tutorials_  
These tutorials use the combined ASV data from the team drive to show an example of how to use the various techniques with our data.  
eDNA_Tree_Tutorial.R  
eDNA_RandomForest_Tutorial.R  
eDNA_GradientForest_Tutorial.R  




#### Useful Papers
1. [Statistical Modeling: The Two Cultures](https://projecteuclid.org/euclid.ss/1009213726)  
2. [Ecological genomics meets community-level modelling of
biodiversity: mapping the genomic landscape of current and
future environmental adaptation](http://onlinelibrary.wiley.com/doi/10.1111/ele.12376/full)  
3. [Genomic divergence across ecological gradients in the Central
African rainforest songbird (Andropadus virens)](http://onlinelibrary.wiley.com/doi/10.1111/mec.14270/full)
