# dry-earth

Code repository for research into the Australian Millennium Drought. Measured and modeled hydrological variables are gathered and analysed in a variety of ways. The analyses make use of techniques like normalization, auto-correlation and lagged correlations between time series.

## code
* The first set of .sh scripts is used to access and download data from multiple sources: hydrological simulations from the [eartH2Observe project](https://wci.earth2observe.eu/) with download\_e2o\_script.sh, a digital elevation model with download\_dem\_script.sh and remotely sensed MODIS data with download\_modis\_script.sh. The scripts also clip data to the spatio-temporal extent of this study.
* A second set of .R scripts constructs the datasets. dataset\_building.R combines the downloaded data into tidy R datasets. prepare\_gw\_data.R presents the construction of a groundwater dataset (not runnable but presupplied (see data section)). The spatial\_classification.R script produces regionalization tables that are joined to the datasets. 
* The last set of .R scripts contain the analyses and operate on the constructed R datasets. They also export the results in the form of figures and tables. analysis1.R is only observation-based and therefore empirical. analysis\_of\_streamflow is complementary to this. analysis2.R has the goal to analyse the accordance of the models with observations and therefore also model-based.

## data
All datasets are hosted on [dataverseNL](https://dataverse.nl/dataverse/geosciences). This concerns the pre-supplied data like the groundwater database, streamflow data and SILO observations, but also the intermediate results of the download scripts and also the higher level R datasets on which the analyses scripts run. 

## reproducability
* The first possibility is to run everything by invoking the meta\_script.sh, which needs only the scripts and presupplied data listed in it. In this case the download scripts will take care of the additional earth2observe and MODIS data and everything will unfold from scratch.
* The second possibility is to directly download the higher level datasets (All\_gw\_data.rds, Set1\_p05deg.RData, Set2\_p25deg.RData), the presupplied streamflow data and the regionalization tables and run only the analysis script that reproduces your desired results.
