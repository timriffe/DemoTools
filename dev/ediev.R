smooth_age_zigzag_dd <- function(
    Value, #population counts by single year of age
    Age, #age index for Value: should be a continuous range of integer values (we will not check that)
    OAG = TRUE, #is the last age group an open interval?
    ret_ggPlot = TRUE, #if TRUE, we return a list(data=ValueOut, ggplot=myPlot) with a plot for checking results
    #if FALSE, we return the smoothed series ValueOut
    legendPlot = NULL #legend for the ggPlot
)
{
  #Extension of Feeney (2013) "Removing Zig-Zag from Age Data" for single age data
  #in the spirit of Dalkhat M. Ediev (2021) "A model of age heaping with applications to population graduation that retains informative demographic variation"
  #Main differences with Feeney's algorithm is that, apart for working with single age data
  #to suppose that age heaping can start 2 years before (and symmetrically 2 years after)
  #and not only 1 year before and after.
  #For example for age 40, we suppose that heaping occurs symmetrically, in the same proportion, at
  #ages 39 and 41, as Feeney supposed for five-years age groups with adjacent age groups,
  #but in the same spirit we suppose that
  #heaping can occur also at ages 38 and 42 towards attracting age 40.
  #We also suppose that heaping can only occurs for ages multiple of 5 and 10
  if (OAG) {
    #if the last age is an open interval, we drop it
    AgeRangeOut <- Age[1:(length(Age)-1)]
  }
  #determine the ages with the variable parameters
  #For example we start with ages 3 and 4 for heaping at age 5 and ages 6 and 7 use symmetrical values
  #generally speaking heaping can start at ages x-2, x-1 and then symmetrically at ages x+1, x+2 towards attracting age x
  #we suppose that the first and the last age we can treat have at least values for 2 ages before and after
  #that means that if the first age in Age is 9 or 10, we will skip the first
  getAgeRange <- function(AgeRangeOut) {
    ageRangeMin <- max(5, trunc((AgeRangeOut[1] + 2) / 5) * 5)
    if (ageRangeMin < AgeRangeOut[1] + 2) ageRangeMin <- ageRangeMin + 5
    ageRangeMax <- trunc((AgeRangeOut[length(AgeRangeOut)] - 2) / 5) * 5
    if (ageRangeMax > AgeRangeOut[length(AgeRangeOut)] - 2) ageRangeMax <- ageRangeMax - 5
    return (c(ageRangeMin, ageRangeMax))
  }
  ageRange <- getAgeRange (AgeRangeOut)
  ageRangeMin <- ageRange[1]
  ageRangeMax <- ageRange[2]
  #number of ages multiples of 5
  nAges_mult <- trunc ((ageRangeMax - ageRangeMin) / 5) + 1
  if (nAges_mult < 2) stop ("Age range should include at least two consecutive ages multiple of 5")
  #initial values for the parameters
  #if first attracting age is 5, age 3 has initial proportion of 0.1, age 4 of 0.2
  #these two parameters will get 'optimized' by R optim
  #the same parameters are applied, symmetrically, for age 6 (0.2) and age 7 (0.1)
  param <- rep (c(0.1, 0.2), nAges_mult)
  param_min <- rep(0, nAges_mult * 2)
  param_max <- rep(9, nAges_mult * 2)
  
  computeNewPop <- function(data, par, Age, ageRangeMin, nAges_mult, ret_Props=FALSE) {
    if (ret_Props) props <- rep(NA, length(data))
    for (ageIter in (1:nAges_mult)) {
      currFirstAge <- ageRangeMin + (ageIter - 1) * 5 - 2
      indexInData <- currFirstAge - Age[1] + 1
      if (ret_Props) {
        props[indexInData] <- par[(ageIter - 1) * 2 + 1]
        props[indexInData+1] <- par[(ageIter - 1) * 2 + 2]
        props[indexInData+3] <- par[(ageIter - 1) * 2 + 2]
        props[indexInData+4] <- par[(ageIter - 1) * 2 + 1]
      } else {
        data[indexInData + 2] <- data[indexInData + 2] -
          (data[indexInData] + data[indexInData + 4]) * par[(ageIter - 1) * 2 + 1] -
          (data[indexInData + 1] + data[indexInData + 3]) * par[(ageIter - 1) * 2 + 2]
        data[indexInData] <- data[indexInData] * (1 + par[(ageIter - 1) * 2 + 1])
        data[indexInData+1] <- data[indexInData+1] * (1 + par[(ageIter - 1) * 2 + 2])
        data[indexInData+3] <- data[indexInData+3] * (1 + par[(ageIter - 1) * 2 + 2])
        data[indexInData+4] <- data[indexInData+4] * (1 + par[(ageIter - 1) * 2 + 1])
      }
    }
    
    if (ret_Props) {
      return (props)
    } else {
      return (data)
    }
  }
  
  optim_fun <- function(data, par) {
    popIter <- computeNewPop(data, par, Age, ageRangeMin, nAges_mult)
    totPop <- sum (popIter)
    objective_to_min <- 0
    for (ageIter in (1:nAges_mult)) {
      currFirstAge <- ageRangeMin + (ageIter - 1) * 5 - 2
      indexInData <- currFirstAge - Age[1] + 1
      term1 <- ((popIter[indexInData + 2] - (popIter[indexInData + 1] + popIter[indexInData + 3]) / 2)) / totPop
      term2 <- ((popIter[indexInData + 2] - (popIter[indexInData] + popIter[indexInData + 4]) / 2)) / totPop
      objective_to_min <- objective_to_min + term1^2 + term2^2
    }
    return (objective_to_min)
  }
  
  oldOptimValue <- 10^1000
  
  for (rep in (1:10)) {
    #apply optimization up to 10 times, using parameters obtained in the previous step, just in case...
    parResult <- optim(data=Value, par=param, lower=param_min, upper=param_max, fn=optim_fun, method = "L-BFGS-B")
    if (parResult$value >= oldOptimValue) break
    param <- parResult$par
    oldOptimValue <- parResult$value
  }
  valueOut <- computeNewPop(Value, parResult$par, Age, ageRangeMin, nAges_mult)
  
  if (ret_ggPlot) {
    library2 <- function(pack) {
      #pack<-"tcltk2"
      if( !(pack %in% installed.packages()))
      {install.packages(pack)}
      library(pack,character.only = TRUE)
    }
    library2("ggplot2")
    library2("scales")
    
    getScales <- function(Value, ageRangeMin, ageRangeMax) {
      maxCount <- max(Value) * 1.0
      power <- 0
      while (maxCount > 1) {
        maxCount <- maxCount / 10.0
        power <- power + 1
      }
      maxCount <- ceiling(maxCount * 10) * 10^(power-1)
      leftScale <- seq(0, maxCount, maxCount / 4)
      rightScale <- seq(0, 1, 0.25)
      xScale <- seq(ageRangeMin-5, ageRangeMax+5, 5)
      return (list(left=leftScale, right=rightScale, x=xScale))
    }
    
    df <- data.frame(Age=Age, val=Value)
    df$type <- 'Original'
    temp <- data.frame(Age=Age, val=valueOut)
    temp$type <- 'Corrected'
    df <- rbind(df, temp)
    props <- computeNewPop(Value, parResult$par, Age, ageRangeMin, nAges_mult, ret_Prop = TRUE)
    temp <- data.frame(Age=Age, val=props)
    temp$type <- 'Attraction'
    #transform into proportions
    temp$val <- temp$val / (1 + temp$val)
    scale_axes <- getScales (Value, ageRangeMin, ageRangeMax)
    maxY_left <- scale_axes$left[length(scale_axes$left)]
    #adjust values in preparation for second axis (we will correct the labels later on in scale_y_continuous)
    temp$val <- temp$val * maxY_left
    df <- rbind(df, temp)
    df$type <- factor (df$type, levels=c("Original", "Corrected", "Attraction"))
    myLabelNames <- c("Original", "Corrected", "Attraction factors")
    
    myPlot <- ggplot(df, aes(x=Age, y=val, color=type, alpha=type, size=type)) + geom_line()
    myPlot <- myPlot + scale_y_continuous(breaks = scale_axes$left, name="population counts",
                                          sec.axis = sec_axis(as.formula(paste("~./", maxY_left, sep="")), name="attraction proportions",
                                                              breaks=scale_axes$right),
                                          labels = scales::label_number())
    myPlot <- myPlot + scale_x_continuous(breaks = scale_axes$x)
    myPlot <- myPlot + scale_colour_manual(values=c("red", "blue", "grey40"), labels=myLabelNames)
    myPlot <- myPlot + scale_size_manual(values=c(1, 2, 2), labels=myLabelNames)
    myPlot <- myPlot + scale_alpha_manual(values=c(1, 1, 0.2), labels=myLabelNames)
    myPlot <- myPlot + theme_bw() + theme(legend.key = element_blank(), legend.title=element_blank(), legend.position="bottom")
    myPlot <- myPlot + theme_text(1)
    myPlot <- myPlot + labs(color="", alpha="", size="")
    if (!is.null(legendPlot)) myPlot <- myPlot + ggtitle(legendPlot)
    
    return (list(data=valueOut, ggplot=myPlot))
  } else {
    return (valueOut)
  }
}

#data for Colombia 1973
col1973 <- structure(list(
  Age = 0:99,
  Value = c(493538L, 497954L, 580000L, 
            610204L, 609266L, 609670L, 599184L, 628352L, 624302L, 549158L, 
            624112L, 527810L, 611266L, 541084L, 490342L, 466650L, 435310L, 
            423912L, 444468L, 334732L, 345390L, 301686L, 319468L, 312100L, 
            267828L, 293902L, 230690L, 230770L, 244100L, 178590L, 278284L, 
            136552L, 201120L, 218514L, 159054L, 217800L, 156786L, 162270L, 
            203560L, 141258L, 247758L, 94344L, 166104L, 164336L, 116024L, 
            182858L, 104268L, 102916L, 136222L, 90550L, 180598L, 68132L, 
            109788L, 108058L, 81438L, 116446L, 77230L, 62844L, 73906L, 49812L, 
            136486L, 37506L, 61658L, 69702L, 43570L, 72760L, 36930L, 34102L, 
            40052L, 24126L, 68924L, 16004L, 32988L, 40868L, 18406L, 31836L, 
            14222L, 11824L, 16206L, 6972L, 24294L, 5016L, 7392L, 6906L, 6244L, 
            9424L, 3806L, 3286L, 3478L, 1912L, 5260L, 916L, 1540L, 1142L, 
            1266L, 1674L, 602L, 468L, 752L, 2360L)),
  class = "data.frame",
  row.names = c(NA,-100L))

#Alexandropol 1897
Alex1897 <- structure(list(
  Age = 0:99,
  Value = c(376.7579251, 342.1757925, 
          398.3717579, 337.8530259, 314.0778098, 290.3025937, 355.1440922, 
          334, 345, 360, 392, 260.0432277, 370, 307.5936599, 290.3025937, 
          247.074928, 281.6570605, 268.6887608, 298.9481268, 164.9423631, 
          333.5302594, 1401.253602, 1729.783862, 1643.32853, 1669.26513, 
          700.9654179, 324.8847262, 273.0115274, 234.1066282, 104.4236311, 
          471.8587896, 108.7463977, 186.556196, 151.9740634, 113.0691643, 
          389.7262248, 151.9740634, 126.037464, 156.29683, 48.22766571, 
          458.8904899, 48.22766571, 108.7463977, 56.87319885, 65.51873199, 
          277.3342939, 82.80979827, 52.55043228, 91.45533141, 48.22766571, 
          363.7896254, 26.61383285, 52.55043228, 39.58213256, 39.58213256, 
          199.5244957, 56.87319885, 39.58213256, 48.22766571, 5, 285.9798271, 
          17.96829971, 26.61383285, 17.96829971, 22.29106628, 78.4870317, 
          22.29106628, 22.29106628, 26.61383285, 9.322766571, 121.7146974, 
          0.677233429, 9.322766571, 5, 5, 30.93659942, 13.64553314, 9.322766571, 
          9.322766571, 13.64553314, 35.25936599, 9.322766571, 9.322766571, 
          5, 0.677233429, 17.96829971, 0.677233429, 5, 0.677233429, 9.322766571, 
          17.96829971, 5, 5, 5, 5, 0.677233429, 0.677233429, 0.677233429, 
          0.677233429, 0.677233429)),
  class = "data.frame",
  row.names = c(NA, -100L))

#call the function for Colombia
resCol <- smooth_age_zigzag_dd(col1973$Value, col1973$Age, OAG = TRUE, ret_ggPlot = TRUE, legendPlot="Colombia 1973")
#show the plot for Colombia
resCol$ggplot
#call the function for Alexandropol
resAlex <- smooth_age_zigzag_dd(Alex1897$Value, Alex1897$Age, OAG = TRUE, ret_ggPlot = TRUE, legendPlot="Alexandropol 1897")
#show the plot
resAlex$ggplot
