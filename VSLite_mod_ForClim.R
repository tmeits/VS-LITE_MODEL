###################################################################################################
## Journal title: Forward modeling of tree-ring width improves simulation of forest growth responses to drought 
## Journal: Agricultural and Forest Metereology
## Authors: Marco Mina, Dario Martin-Benito, Harald Bugmann, Maxime Cailleret
## Affiliation: ETH Zurich, Department of Environmental Sciences, Forest Ecology, Switzerland
##
##
## THIS R-SCRIPT CONTAINS 3 FUNCTIONS: 
##  1) Scaled daylenght function for estimation of proxy for climatological insolation (gE)
##  2) ForClim Water Balance Submodel for calculation of soil moisture (M) and potential evapotranspiration (potEv)  
##  3) Modified VS-LITE MODEL for calculation of growth responses to temperature (gT), moisture (gM), overall growth respose (Gr)
##     and the simulated ring-width index (trw)
##
## The R version of the original VS-Lite model can be found at 
##   https://github.com/suztolwinskiward/VSLiteR.git (accessed on 07/11/2015). 
##
###################################################################################################

#--------SCALED DAYLENGTH FUNCTION
# Computes gE (insolation factor) from latitude according to Gates (1980) Biophysical ecology. Springer, New York
# gE = numeric vector [1:12] with the number of daylight hours depending on the latitude
Compute_gE<-function (phi)
{
  gE = rep(NA,12)
  # Compute normalized daylength (neglecting small difference in calculation for leap-years)
  ndays = c(0,31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  cdays = cumsum(ndays)
  doy=1:365
  phi[phi > 90 | phi < -90] <- NA
  P <- asin(0.39795 * cos(0.2163108 + 2 * atan(0.9671396 * tan(0.0086 * (doy - 186)))))
  a <- (sin(0.8333 * pi/180) + sin(phi * pi/180) * sin(P))/(cos(phi * pi/180) * cos(P))
  a <- pmin(pmax(a, -1), 1)
  DL <- 24 - (24/pi) * acos(a)
  ndl=DL/max(DL)
  for (t in 1:12){gE[t] = mean(DL[(cdays[t]+1):cdays[t+1]])}
  return(gE)
}

#########################################################################

#--------FORCLIM WATER BALANCE SUBMODEL
# Calculate the Soil Moisture (M or SM) and potential evapotranspiration (PET or potEv) based on the soil water balance model described in 
# Bugmann and Cramer (1998), a modified version of the Thornthwaite and Mather water balance model (T&MWB)
#
# References: 
# Bugmann and Cramer (1998). Forest Ecology and Management. 103, 247-263
# Thornthwaite and Mather (1957). Publications in Climatology. 10, 183-311 
# Federer (1982). Water Resour. Res. 18, 355-362
#
# Input data required:
# syear = start (first) year of simulation; eyear = end (last) year of simulation
# Temp = monthly mean temperature values (°C) 
# Precip = monthly precipitation sums (cm)
# kBS = bucket size (i.e., water holding capacity, in cm) ; kLat = latitude 

forclim_WB <- function(syear,eyear,Temp,Precip,kBS,kLat)
{
  
  kCw = 12      # Maximum evaporation rate
  kIcpt = 0.3   # Interception rate of the rainfall
    
  iyear = c(syear:eyear)
  nyrs = length(iyear)
  
  # Storage for growth response output variables (size [12 x Nyears]):
  M <- matrix(NA,nyrs,12)
  potEv <- matrix(NA,nyrs,12)
  
  kA=kB=kLatPtr=gPET=gPi=gS=gD=gE=SM=NULL 
  
  #Parameters for the calculation of the heat index (empirical function; see Bugmann and Cramer 1998; T&MWB)
  kPM = 1.6; k1 = 0.2; k2 = 1.514; k3 = 6.75E-7; k4 = -7.71E-5; k5 = 0.01792; k6 = 0.49239;
  kA[1] = 1.1226; kA[2]=0.9859; kA[3]=1.0454; kA[4]=0.9708; kA[5]=0.9605; kA[6]=0.9185; kA[7]=0.9669; kA[8]=0.9892; kA[9]=0.9900;          
  kA[10]=1.0600; kA[11]=1.0815; kA[12]=1.1444; kB[1]=-7.3094E-3; kB[2]=-3.8701E-3; kB[3]=-4.9231E-4; kB[4]=3.5179E-3; kB[5]=7.1453E-3; 
  kB[6] = 8.4718E-3; kB[7]=7.6410E-3; kB[8]=4.9436E-3; kB[9]=1.2000E-3; kB[10]=-2.6256E-3; kB[11] = -6.3692E-3; kB[12] = -8.6598E-3
  
  # Adjustment of PET based on the slope and aspect (when such data are available at the study site) 
  # Slope and aspect parameter [0=flat; -2=north >30°; +2=south >30°; -1=north <30°; +1=south <30°]
  # For the calculation of the slope and aspect parameter see Schumacher (2004) Ph.D. Thesis, ETH Zьrich 
  #   http://dx.doi.org/10.3929/ethz-a-004818825
  kSlAsp = 0    
  if (kSlAsp > 0) {kPMod = 1 + kSlAsp * 0.125} else {kPMod = 1 + kSlAsp * 0.063}
  
  #Calculation of Soil Moisture (SM) for each month and each year
  for (cyear in c(1:length(iyear)))
  {      
    kHi = 0; for (t in 1:12) {kHi = kHi + (max(c(0, k1 * Temp[cyear,t])))^k2}
    kC = ((k3 * kHi + k4) * kHi + k5) * kHi + k6
    for (t in 1:12) {kLatPtr[t] = kA[t] + kB[t] * kLat} 
    for (t in 1:12) {SM[t] = kBS}   # Initial values for state variables
    gPs = 0                         # Infiltrated precipitation
    
    uAET=NULL; uDr=NULL; uSM=NULL
    for (t in 1:12)
    {
      gPET[t] = kPM * kPMod * kLatPtr[t] * (10 / kHi * max(c(Temp[cyear,t],0)))^kC
      gPi[t] = min(c(kIcpt * Precip[cyear,t]/10, gPET[t])) #intercepted precipitation
      gPs = Precip[cyear,t]/10 - gPi[t] #infiltrated precipitation 
      
      gS[t] = min(c(kBS, kCw)) * SM[t] / kBS; #supply function
      gD[t] = gPET[t] - gPi[t] #demand function
      gE[t] = min(c(gS[t], gD[t])) #evapotranspiration from the soil
      
      ## Error catching for avoiding negative soil moistures
      if (t <= 11) {SM[t + 1] = max(c(min(c(SM[t] + (gPs - gE[t]), kBS)), 0))}   
      
      ## Actual, potential evapotranspiration and soil moisture:
      uAET[t] = gE[t] + gPi[t];
      potEv[cyear,t] = gPET[t]
      M[cyear,t] = SM[t]
      
    }  # end loop for month
  }   # end loop for year
  return(list(M,potEv))
}

#########################################################################


#----------------------------MODIFIED VS-LITE MODEL
#
# Computes growth response to moisture (gM), growth response to temperature (gT), soil moisture (M), overall growth response (Gr),
# simulated ring-width (trw), and potential evapotranspiration (potEv) using temperature response parameters (T1 and T2s),
# moisture response parameters (M1 and M2) for the four seasons, monthly temperature (T) and precipitation (M), and soil bucket size (kBS, in cm)

# Inputs:
#   syear = start year of simulation.
#   eyear = end year of simulation.
#   phi = site latitude (in degrees N)
#   T1 = scalar temperature threshold below which growth response to temperature is zero (in deg. C)
#   T2s = parameter reflecting the shape of the Gompertz curve which influences the growth rate of y as a function of x  
#   M1 = scalar soil moisture threshold below which growth response to soil moisture is zero (in % of the bucket size)
#   M2 = scalar soil moisture threshold above which growth response to soil moisture is one (in % of the bucket size)
#   Note that parameters T1, T2s, seasonal M1, M2 can be estimated simultaneously with an optimization procedure
#	(e.g., by maximizing the correlation between simulated and measured ring-width series; we suggest the use of the library DEoptim for such purpose)
#
#   T   = (12 x Nyrs) matrix of observed mean monthly temperatures (in degrees C)
#   P   = (12 x Nyrs) matrix of observed monthly precipitation sum (in cm)
#   kBs = site-specific water holding capacity (i.e., bucket size, in cm)

# References:
# Mina et al (subm.). Agr.For.Met. 
# Tolwinski-Ward et al (2011). Clim. Dyn. 36:2419-2439 (in MATLAB)

VSlite_v4 <- function(T1,T2s,M1sp,M2sp,M1su,M2su,M1wi,M2wi,M1fa,M2fa,syear,eyear,phi,Temp,Precip,kBS){
  
  iyear = seq(syear,eyear) #index of simulation year
  nyrs = length(syear:eyear) #number of simulated years
  Precip = Precip
  Temp = Temp

  ###### Pre-allocate storage for outputs #####
  Gr = matrix(NA,nyrs,12)
  gT = matrix(NA,nyrs,12)
  gM = matrix(NA,nyrs,12)
  M = matrix(NA,nyrs,12)
  
  gE = Compute_gE(phi) #scaled monthly proxy for insolation
  kLat <- phi #latitude
  
  #Compute potential evapotranspiration and soil moisture:
  potEv <- forclim_WB(syear,eyear,Temp, Precip,kBS,kLat)[[2]]
  M = forclim_WB(syear,eyear,Temp,Precip,kBS,kLat)[[1]]
  
  ###########################################################################
  # Calculate monthly growth response to T & M, and overall growth response G:
  # cyear = year the model is currently working on
  
  for (cyear in c(1:length(iyear))){      # cycle over years
    for (t in 1:12){  # cycle over months within a year
      
      ### Calculate Growth Response functions gT(t)
      # First, temperature growth response with Gompetz s-shaped function:
      x = Temp[cyear,t]
      e <- exp(1)
      As=1  # asymptote (optimal growth)
      gT[cyear,t] <- max (0, As * exp(-exp(T2s * e * (T1 - x)/As + 1))  )
      
      ### Calculate Growth Response functions gM(t) for the different seasons:
      if(t %in% c(1:2,12)) #winter
      {
        x = M[cyear,t]
        if (x < M1wi){gM[cyear,t] = 0}
        if (x >= M1wi & x <= M2wi){gM[cyear,t] = (x - M1wi)/(M2wi - M1wi)}
        if (x >= M2wi){gM[cyear,t] = 1}
      }
      
      if(t %in% c(3:5)) #spring
      {
        x = M[cyear,t]
        if (x < M1sp){gM[cyear,t] = 0}
        if (x >= M1sp & x <= M2sp){gM[cyear,t] = (x - M1sp)/(M2sp - M1sp)}
        if (x >= M2sp){gM[cyear,t] = 1}
      }
      
      if(t %in% c(6:8)) #summer
      {
        x = M[cyear,t]
        if (x < M1su){gM[cyear,t] = 0}
        if (x >= M1su & x <= M2su){gM[cyear,t] = (x - M1su)/(M2su - M1su)}
        if (x >= M2su){gM[cyear,t] = 1}
      }
      if(t %in% c(9:11)) #fall
      {
        x = M[cyear,t]
        if (x < M1fa){gM[cyear,t] = 0}
        if (x >= M1fa & x <= M2fa){gM[cyear,t] = (x - M1fa)/(M2fa - M1fa)}
        if (x >= M2fa){gM[cyear,t] = 1}
      }
      
      # Compute overall monthly growth rate:
      Gr[cyear,t] = gE[t]*min(gT[cyear,t],gM[cyear,t])
      
    } #end monthly cycle    
  } #end yearly cycle
  
  #Integration window parameters
  I_0 = -4    #start month (-4 means September of the previous year)
  I_f = 12    #end month (12 means December of the current year)
  
  
  #Compute simulated series of ring-width index from growth responses
  width = matrix(nrow=length(syear:eyear), ncol=1,NA)
  if (phi>0){ # for sites in the Northern Hemisphere only
    if (I_0<0){ # if we include part of the previous year in each year's modeled growth:
      startmo = 13+I_0
      endmo = I_f
      # use average of growth data across modeled years to estimate first year's growth due to previous year:
      width[1] = sum(Gr[1,1:endmo]) + mean(Gr[,startmo:12])
      for (cyear in 2:length(syear:eyear)){
        width[cyear] = sum(Gr[cyear-1,c(startmo:12)]) + sum(Gr[cyear,1:endmo])
      }
    }else{ # no inclusion of last year's growth conditions in estimates of this year's growth:
      startmo = I_0+1
      endmo = I_f
      for (cyear in 1:length(syear:eyear)){
        width[cyear] = sum(Gr[cyear,startmo:endmo])}
    }
  }

  trw = ((width-mean(width))/apply(width,2,sd)) #proxy series standardized width
  
  return(list(gM, gT, Gr,trw))
}



