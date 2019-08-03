#The central limit theorem tells us that no matter what distribution things have, 
#the sampling distribution tends to be normal if the sample is large enough (n > 30).

library(tidyverse)
library(dplyr)
library(magrittr)
library(broom)
library(ggplot2)
library(qqplotr)
library(purrr)
library(gridExtra)
library(GWASTools)
load("/home/brunettt/PERU_GWAS_Rasika_Ayo/PERU_Ayo_Rasika_GENESIS_PC_pheno_covar_data_updated_sex_covars_normalized_addedIgEspecificPhenos_caseControlSplitsAdded.Rdata")
# transformations: square/quadratic, cube root, log10, natural log, squart root, square root of the square root, reciprical, log of log -- all should have Shapiro-Wilks test for normality

pdf()

checkThese <- c('hx2_cleaned', 'ses')
tblFull <- as_tibble(pData(scanAnnot))
tblSub <- tblFull %>% select(checkThese)


plotGraphs <- function(graphObjs){
  grid.arrange(graphObjs[[1]], graphObjs[[3]],
               graphObjs[[4]], graphObjs[[5]],
               graphObjs[[6]], graphObjs[[7]],
               graphObjs[[8]], graphObjs[[9]],
               graphObjs[[10]],
               nrow = 3, ncol = 3)
}

transformSuite <- function(tmpTbl, variable){
  tmpVar <- c('sqrt', 'squareQuadratic', 'cubeRoot', 'log10', 'ln', 'sqrtOfsqrt', 'reciprocal', 'loglog')
  #tmpTbl %>% ggplot(aes_(x=variable)) + geom_density() # original
  tmpTbl %<>% mutate(!!tmpVar[1] := sqrt(!!(sym(variable))))
  tmpTbl %<>% mutate(!!tmpVar[2] := (!!(sym(variable)))^2)
  tmpTbl %<>% mutate(!!tmpVar[3] := (!!(sym(variable)))^(1/3))
  tmpTbl %<>% mutate(!!tmpVar[4] := log10((!!(sym(variable)))))
  tmpTbl %<>% mutate(!!tmpVar[5] := log((!!(sym(variable)))))
  tmpTbl %<>% mutate(!!tmpVar[6] := sqrt(sqrt((!!(sym(variable))))))
  tmpTbl %<>% mutate(!!tmpVar[7] := 1/(!!(sym(variable))))
  tmpTbl %<>% mutate(!!tmpVar[8] := log(log((!!(sym(variable))))))
  tmpTbl %>% names() %>% map(~ ggplot(tmpTbl, aes_string(x=.x)) + geom_density() + 
                               ggtitle(paste("W=", shapiro.test(tmpTbl %>% pull(.x)) %>% glance() %>% pull(statistic) %>% round(digits=3),
                                              "pValue=",  shapiro.test(tmpTbl %>% pull(.x)) %>% glance() %>% pull(p.value) %>% formatC(format = "e", digits=3), 
                                             sep = " ")) + theme(plot.title = element_text(hjust = 0.5)))

}


#x is a tibble of single column, !! means to unquote a variable, therefore not evaluate it, while the sym() while sym() removes the quotation marks around an input string
transformChecks <- function(tmpTbl, variable){
  tmpTbl %<>% drop_na()
  if(tmpTbl %>% filter((!!sym(variable)) < 5) %>% nrow() > 0){ #make sure no zeros otherwise fails ln and log10 transforms, so add pseudocount
    tmpVar = paste(variable, 'pseudoCount', sep = '_') 
    tmpTbl %<>% mutate(!!tmpVar := (!!sym(variable)) + 5) # := for mutate allows for dynamic assignment of new parameter names which is why I can dynamically assign tmpVar as long as I use !!
    transformSuite(tmpTbl, tmpVar)

  }else{
    print("zeroCheck: PASS")
    transformSuite(tmpTbl, variable)
  }
  # TO DO: check for negatives
}

  
# from the tibble subset extract the column names as a list and iterativly apply each column name in the list to the the tranformation function (note that . means one value in output of names)
graphObjs <- tblSub %>% names() %>% map(~ transformChecks(tblSub[.x], .x))
graphs <- graphObjs %>% map( ~ plotGraphs(.x))


