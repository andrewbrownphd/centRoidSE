#' @export centroidMethod
centroidMethod <- function(input, 
                           decimals=NULL, 
                           resolution=250, 
                           type="OR", 
                           distribution="z",
                           df=NULL, 
                           plot=FALSE, 
                           ci=95, 
                           pInfo=NULL)
{
  
  # Check for legal values
  # df can be added to 'input' or put in its own column, 
  # but must be specified for distribution=t
  
  if(length(input) == 4 & !is.null(df)) 
    stop("Four items entered for input and df also specified.")
  
  # set df if input has four items
  
  if(length(input) == 4) df <- as.numeric(input[4])
  
  #make sure df is specified
  
  distribution <- tolower(distribution)
  
  if(is.null(df) & distribution=="t") 
    stop("Degrees of freedom must be specified for t distribution.")
  if(!(distribution %in% c("z","t"))) 
    stop("Distribution must specify 'z' or 't'.")
  
  #check to make sure type is appropriately specified
  
  if(!(type %in% c("OR","log(OR)","beta","mean")))
    stop(paste("Type must be \"OR\", \"log(OR)\", \"beta\", or \"mean\".",
               "Note: \"log(OR)\", \"beta\", and \"mean\" produce the same calculations."))
  
  if(type=="OR" & distribution=="t")
  {
    warning("type=\"OR\" is incompatible with distribution=\"t\"; setting distribution= to \"z\"")
    distribution <- "z"
  }
  
  #Process decimals
  
  if(is.null(decimals))
  {
    for (i in 1:3)
    {
      # Extract decimals if unspecified. Note that numbers entered with trailing zeroes will truncate
      # unless entered as character
      decimals[i] <- nchar(strsplit(as.character(input[i]), 
                                    split=".", fixed=T)[[1]][2]) #
      if(is.na(decimals[i])){
        warning(paste("Extracting decimals from input value",
                      i,
                      "resulted in an NA. This may occur if the input has no decimal places",
                      "(e.g., an OR input as '1' instead of '1.0') or if the value is an integer.",
                      "Setting the decimal for that input value to zero."))
        decimals[i] <- 0
      }
    }
  } else {
    # If only one decimal is specified, then enter it in the vector for each of three elements
    if(length(decimals) == 1) decimals <- c(decimals,decimals,decimals)
  }
  if(decimals[1] != decimals[2] | decimals[2] != decimals[3])
    warning(paste("Unequal decimals extracted for input: ",paste(input,collapse=",")))
  
  #Convert input to numeric in case entered as character
  
  input <- as.numeric(input)
  sumna <- function(a) { return(sum(is.na(a))) }
  if(sumna(input) > 0)
    stop("Unable to convert inputs into numbers; check input for numeric values.")
  
  #Calculate SE from lower and upper CIs
  
  if(distribution=="z")
  {
    crit<-qnorm((100-(100-ci)/2)/100)
  }
  
  if(distribution=="t")
  {
    crit<-qt((100-(100-ci)/2)/100,df=df)
  }
  
  if(type == "OR")
  {
    #gets four most extreme rounding combinations, but may make impossible SEs
    ll <- (log(input[1] + (5*10^(-1-decimals[1]))) #increase OR
           - log(input[2] - (5*10^(-1-decimals[2])))) /crit #decrease LCL
    uu <- (log(input[3] + (5*10^(-1-decimals[3]))) #increase UCL
           - log(input[1] - (5*10^(-1-decimals[1]))))/crit #decrease OR
    lu <- (log(input[1] - (5*10^(-1-decimals[1]))) #decrease OR
           - log(input[2] + (5*10^(-1-decimals[2])))) /crit #increase LCL
    ul <- (log(input[3] - (5*10^(-1-decimals[3]))) #decrease UCL
           - log(input[1] + (5*10^(-1-decimals[1]))))/crit #increase OR
    
    l <- min(ll, uu, lu, ul)
    
    # this bounds SE to be >= 0
    l <- ifelse(l < 0, 0, l)
    u <- max(ll, uu, lu, ul)
    
    SE_l <- (log(input[1])-log(input[2]))/crit #backcalc SE from LCL
    SE_u <- (log(input[3])-log(input[1]))/crit #backcalc SE from UCL
    
    #Sequence of possible log(OR)
    thetas <- seq(input[1]-(5*10^(-1-decimals[1])),
                  input[1]+(5*10^(-1-decimals[1])),
                  length.out=resolution) 
    thetas <- log(thetas) #on log scale, search area
  }
  
  if(type %in% c("log(OR)","mean","beta"))
  {
    
    # Back calc SE from LCL and UCL 
    SE_l <- (input[1]-input[2])/crit
    SE_u <- (input[3]-input[1])/crit
    
    #set l to min and u to max of the two SE estimates
    l <- min(SE_l, SE_u)
    u <- max(SE_l, SE_u)
    
    #Sequence of possible point estimates
    thetas <- seq(input[1]-(5*10^(-1-decimals[1])),
                  input[1]+(5*10^(-1-decimals[1])),
                  length.out=resolution) # search area
  }
  
  #Sequence of possible SEs between low and high estimates from CI
  
  if(round(SE_l,digits=8) == round(SE_u,digits=8))
    warning(paste0("Lower and upper SE estimates are equal. Point estimate: ",
                   sprintf(paste0("%.",decimals[1],"f"),input[1]), "; SE: ",SE_l))
  
  SEs <- seq(l,u,length.out=resolution)  # SE for log OR, search area
  
  #Test matrix; 1 means satisfies all rounding criteria, 0 otherwise
  
  mat <- matrix(0, length(thetas), length(SEs))
  
  if(is.null(pInfo) & distribution=="z")
  {
    for(i in 1:length(thetas))
      for(j in 1:length(SEs))
      {
        pe <- thetas[i] ; se <- SEs[j]
        x1 <- round(exp(pe),decimals[1])
        x2 <- round(exp(pe - crit*se),decimals[2]) #lower CI
        x3 <- round(exp(pe + crit*se),decimals[3]) #upper CI
        
        #check rounding
        if(x1==input[1] & x2==input[2] & x3==input[3]) mat[i,j] <- 1
      }
  }
  
  if(is.null(pInfo) & distribution=="t")
  {
    for(i in 1:length(thetas))
      for(j in 1:length(SEs))
      {
        pe <- thetas[i] ; se <- SEs[j]
        x1 <- round(pe,decimals[1])
        x2 <- round(pe - crit*se,decimals[2]) #lower CI
        x3 <- round(pe + crit*se,decimals[3]) #upper CI
        
        #check rounding
        if(x1==input[1] & x2==input[2] & x3==input[3]) mat[i,j] <- 1
      }
  }
  
  if(!is.null(pInfo) & distribution=="z")
  {
    for(i in 1:length(thetas))
      for(j in 1:length(SEs))
      {
        pe <- thetas[i] ; se <- SEs[j]
        x1 <- round(exp(pe),decimals[1]) 
        x2 <- round(exp(pe - crit*se),decimals[2]) #lower CI
        x3 <- round(exp(pe + crit*se),decimals[3]) #upper CI
        x4 <- 2*(1-pnorm(abs(pe/se))) #p-value if z distribution
        
        #check rounding
        if(x1==input[1] & x2==input[2] & x3==input[3])
        {
          #check p-value against p value operator
          if(pInfo[3] == "==" & get(pInfo[3])(round(x4,as.numeric(pInfo[2])),as.numeric(pInfo[1]))) mat[i,j] <- 1 
          if(pInfo[3] %in% c("<",">") & get(pInfo[3])(x4,as.numeric(pInfo[1]))) mat[i,j] <- 1
        }
      }
  }
  
  if(!is.null(pInfo) & distribution=="t")
  {
    for(i in 1:length(thetas))
      for(j in 1:length(SEs))
      {
        pe <- thetas[i] ; se <- SEs[j]
        x1 <- round(pe,decimals[1]) 
        x2 <- round(pe - crit*se,decimals[2]) #lower CI
        x3 <- round(pe + crit*se,decimals[3]) #upper CI
        x4 <- 2*(1-pt(abs(pe/se),df=df)) #p-value if t distribution
        
        #check rounding
        if(x1==input[1] & x2==input[2] & x3==input[3])
        {
          #check p-value against p value operator
          if(pInfo[3] == "==" & get(pInfo[3])(round(x4,as.numeric(pInfo[2])),as.numeric(pInfo[1]))) mat[i,j] <- 1 
          if(pInfo[3] %in% c("<",">") & get(pInfo[3])(x4,as.numeric(pInfo[1]))) mat[i,j] <- 1
        }
      }
  }
  
  m <- sum(mat)
  onepoints <- matrix(NA, m, 2)
  k <- 0
  for(i in 1:length(thetas))
    for(j in 1:length(SEs))
    {
      if(mat[i,j]==1)
      {
        k <- k + 1
        onepoints[k,] <- c(thetas[i], SEs[j])
      }
    }
  centroid <- apply(onepoints,2,mean) # Might be a way to make this more efficient, but this is simplest
  
  if(plot)
  {
    if(type=="OR")
    {
      plot(0, 0, type="n", xlim=c(min(thetas),max(thetas)),
           ylim=c(min(SEs),max(SEs)), xlab="log(OR)", ylab="SE of log(OR)",
           main="Values within log(OR) and SE bounds")
      for(i in 1:length(thetas))
        for(j in 1:length(SEs))
          if(mat[i,j]==1)
            points(thetas[i], SEs[j], pch=".",col="gray")
      #
      points(centroid[1],centroid[2],pch=16,col=2)
      points(log(input[1]),SE_l,pch="L",col=4)
      points(log(input[1]),SE_u,pch="U",col=4)
      points(log(input[1]),mean(c(l,u)),pch=16)
    }
    if(type %in% c("log(OR)","mean","beta"))
    {
      plot(0, 0, type="n", xlim=c(min(thetas),max(thetas)),
           ylim=c(min(SEs),max(SEs)), xlab="Point Estimate",
           ylab="SE of Point Estimate",
           main="Values within Point Estimate and SE bounds")
      for(i in 1:length(thetas))
        for(j in 1:length(SEs))
          if(mat[i,j]==1)
            points(thetas[i], SEs[j], pch=".",col="gray")
      #
      points(centroid[1],centroid[2],pch=16,col=2)
      points(input[1],l,pch="L",col=4)
      points(input[1],u,pch="U",col=4)
      points(input[1],mean(c(l,u)),pch=16)
    }
  }
  
  
  if(type=="OR") 
  {
    names(centroid) <- c("log(OR)","SE")
  } else {
    names(centroid) <- c(type,"SE")
  }
  
  if(sumna(centroid)==2) # Both NaNs
    message(noquote("There is no pair of values that satisfies the given constraints; check input."))
  return(centroid)
}
#
