# HPR

This is the repository for the \texttt{HPR} R package, which accompanies Chase, Taylor, and 
Boonstra's paper on horseshoe process regression (HPR; 2022+). Interested users should
examine that paper, along with the manual and vignette posted here, for examples and explanations
for how to use the package.

Users should first install \texttt{cmdstanr} and \texttt{cmdstan} before attempting to install \texttt{HPR}. To do so,
we first recommend starting a fresh R session. Then, run the commands `install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))`, followed by `library(cmdstanr)` and `check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)`.

If your toolchain is configured correctly, then run: `install_cmdstan(cores = 2)`.

To check that \texttt{cmdstan} was installed correctly, run: `cmdstan_path()` and `cmdstan_version()`.

At this point, you should be all set. If you have trouble with any of the steps above, please
see the \texttt{cmdstanr} documentation for more information:

- https://mc-stan.org/cmdstanr/

- https://mc-stan.org/cmdstanr/articles/cmdstanr.html

After \texttt{cmdstanr} and texttt{cmdstan} have successfully been installed, our package can be installed
by running the command `install.packages("devtools")` followed by 
`devtools::install_github("elizabethchase/HPR")`. At this point, you should be all set.

If you are still having trouble, we recommend checking that you have updated to the latest
versions of both R and RStudio and then attempt reinstallation.
