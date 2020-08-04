# arctic-plants-busco

Code to identify BUSCO/other markers in ancient samples

To recreate the analyses follow the instructions below:

> This has been tested on R 3.6.3

Clone the repo with:

```bash
git clone https://github.com/genomewalker/arctic-plants-busco.git
cd arctic-plants-busco
```

Then let's install the packages we used to plot the figures. First start R to get renv installed:

```
R
```

> If you open the project file `protein-coverage.Rproj` in Rstudio it will perform the same steps.

If everything went well, [renv](https://rstudio.github.io/renv/articles/renv.html) will be installed and you will get a message like:

```
* Installing renv 0.11.0 ... Done!
Successfully installed and loaded renv 0.11.0.
* Project '~/Desktop/repos/arctic-plants-busco' loaded. [renv 0.11.0]
```

And restore the environment:

```r
renv::restore()
q()
```
