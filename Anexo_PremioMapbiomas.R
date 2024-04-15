### Script Consolidado
rm(list = ls())
library(tidyverse)
library(summarytools)
library(plm)
library(car)
library(lmtest)
library(stargazer)
library(fixest)
library(sandwich)
library(writexl)
library(readxl)
library(reshape)


## ETAPAS COMUNS ####

# FUNÇÃO OLEA-PFLUEGUER e ANDERSON RUBIN ####-----

# Olea
olea_f_fixest <- function(endog, instrument, data, controls = NULL, fixed_effects=NULL,  vcov.type = "cluster",...){
  
  
  data = data[complete.cases(data[,c(endog,instrument,controls, fixed_effects)]),]
  if(is.null(controls))
    controls.formula = "" else controls.formula = paste("+",paste(controls, collapse="+"))
    
    if(is.null(fixed_effects))
      fe.formula = "" else fe.formula = paste("|",paste(fixed_effects, collapse="+"))
      
      instrument.formula = paste(instrument, collapse = "+")
      
      fs.formula = as.formula(paste(endog, "~", instrument.formula, controls.formula, fe.formula))
      
      fs.reg = feols(fs.formula, data = data, vcov =vcov.type)
      
      summ= summary(fs.reg,...)  
      
      #print(summ$cov.scaled)
      
      vscaled = summ$cov.scaled[instrument, instrument]
      
      coef.inst = coefficients(fs.reg)[instrument]
      
      nobs = fs.reg$nobs -  fs.reg$nparams
      mat_instr = c()
      
      
      for(inst in c(endog,instrument))
        mat_instr =cbind(mat_instr, resid(feols(as.formula(paste(inst, "~ 1", controls.formula, fe.formula)), data =data),na.rm=F))
      
      vl = (t(mat_instr[,1])%*%mat_instr[,-1])
      
      QZ = (t(mat_instr[,-1])%*%mat_instr[,-1])
      if(!is.matrix(vscaled))
        vscaled = as.matrix(vscaled)
      num  = vl%*%solve(QZ)%*%t(vl)
      
      decomp = eigen(QZ)
      
      if(length(decomp$values)==1)
        sand = matrix(sqrt(QZ)) else sand = decomp$vectors%*%diag(sqrt(decomp$values))%*%t(decomp$vectors)
      
      return(num/(sum(diag(sand%*%vscaled%*%sand))))
      
      
}
## Anderson Rubin Test 
anderson_rubin_test <- function(outcome, endogenous, instruments, vcov, data, beta_0=0, controls = c(), intercept = T, weights = NULL, cluster = NULL, ...)
{
  #Dof for AR test
  dof = length(instruments)
  
  if(length(controls)>0)
    data.kept = data[,c(outcome, endogenous,instruments, controls)] else data.kept = data[,c(outcome, endogenous,instruments)]
    
    keep_ind = complete.cases(data.kept)
    
    data.kept = data.kept[keep_ind,]
    
    Nobs = nrow(data.kept)
    
    if(intercept)
    {
      data.kept = cbind(data.kept, "intercept" = 1)
      controls = c(controls, "intercept")
    }
    
    #We will pool outcome and endogenous now
    data.pooled = rbind(data.kept, data.kept)
    
    data.pooled$pool_variable = c(data.kept[,outcome], data.kept[,endogenous])
    data.pooled$variable_indicator = as.factor(c(rep("reduced_form",    Nobs),rep("first_stage",Nobs)))
    
    
    
    #Constructing the formula for regression
    if(length(controls)>0)
      formula = paste("pool_variable ~ -1 + ", paste(paste("variable_indicator", instruments, sep = ":"),collapse = "+"), "+", paste(paste("variable_indicator", controls, sep = ":"),collapse = "+")) else  formula = paste("pool_variable ~ -1 +", paste(paste("variable_indicator", instruments, sep = ":"),collapse = "+")) 
    
    
    if(is.null(weights))
      pool.model = lm(as.formula(formula), data.pooled) else  pool.model = lm(as.formula(formula), data.pooled, weights = rep(weights[keep_ind],2)) 
    
    coefs = pool.model$coefficients
    
    if(!is.null(cluster))
      vcov_model = vcov(pool.model,cluster = rep(data[keep_ind, cluster],2), ...) else vcov_model = vcov(pool.model, ...)
    
    
    lin_vec = 1*grepl(paste("reduced_form", instruments, sep = ":"), names(coefs)) - 
      beta_0*grepl(paste("first_stage", instruments, sep = ":"), names(coefs))
    
    #constructing test statistic
    val = (coefs%*%lin_vec)
    vcov_lin = t(lin_vec)%*%vcov_model%*%lin_vec
    
    ar = val%*%solve(vcov_lin)%*%val
    
    pvalue =1 - pchisq(ar, dof)
    
    return(list("AR test statistic" = ar, "P-value" = pvalue, "Nobs" = Nobs, "Dof" = dof))
}

#Anderson Rubin Confidence Intervalor (CI)
anderson_rubin_ci <- function(outcome, endogenous, instruments, vcov, data, grid_beta, confidence = 0.95, controls = c(), intercept = T, weights = NULL, cluster = NULL, ...)
{
  #Dof for AR test
  dof = length(instruments)
  
  if(length(controls)>0)
    data.kept = data[,c(outcome, endogenous,instruments, controls)] else data.kept = data[,c(outcome, endogenous,instruments)]
    
    keep_ind = complete.cases(data.kept)
    
    data.kept = data.kept[keep_ind,]
    
    Nobs = nrow(data.kept)
    
    if(intercept)
    {
      data.kept = cbind(data.kept, "intercept" = 1)
      controls = c(controls, "intercept")
    }
    
    #We will pool outcome and endogenous now
    data.pooled = rbind(data.kept, data.kept)
    
    data.pooled$pool_variable = c(data.kept[,outcome], data.kept[,endogenous])
    data.pooled$variable_indicator = as.factor(c(rep("reduced_form",    Nobs),rep("first_stage",Nobs)))
    
    
    
    #Constructing the formula for regression
    if(length(controls)>0)
      formula = paste("pool_variable ~ -1 + ", paste(paste("variable_indicator", instruments, sep = ":"),collapse = "+"), "+", paste(paste("variable_indicator", controls, sep = ":"),collapse="+")) else  formula = paste("pool_variable ~ -1 +", paste(paste("variable_indicator", instruments, sep = ":"),collapse = "+")) 
    
    
    if(is.null(weights))
      pool.model = lm(as.formula(formula), data.pooled) else  pool.model = lm(as.formula(formula), data.pooled, weights = rep(weights[keep_ind],2)) 
    
    coefs = pool.model$coefficients
    
    if(!is.null(cluster))
      vcov_model = vcov(pool.model,cluster = rep(data[keep_ind, cluster],2), ...) else vcov_model = vcov(pool.model, ...)
    
    p1 = grepl(paste("reduced_form", instruments, sep = ":"), names(coefs))
    p2 = grepl(paste("first_stage", instruments, sep = ":"), names(coefs))
    
    acc_vec = c() 
    
    #Looping over grid
    for(beta in grid_beta)
    {
      
      lin_vec = p1 - beta*p2
      #constructing test statistic
      val = (coefs%*%lin_vec)
      vcov_lin = t(lin_vec)%*%vcov_model%*%lin_vec
      
      ar = val%*%solve(vcov_lin)%*%val
      
      pvalue = pchisq(ar, dof)
      
      if(pvalue<= confidence)
        acc_vec = c(acc_vec, T) else acc_vec = c(acc_vec, F)
      
    }
    
    if(sum(acc_vec)==0)
      return("Confidence set is empty!") else {
        
        vec_region_start = c()
        vec_region_end = c()
        if(acc_vec[1]==T)
        {
          warning("Lower boundary point was accepted! Perhaps decrease grid lower bound to see what happens?")
          vec_region_start = grid_beta[1]
        }
        
        if(acc_vec[length(acc_vec)]==T)
        {
          warning("Upper boundary point was accepted! Perhaps increase grid upper bound to see what happens?")
          vec_region_end = grid_beta[length(acc_vec)]
        }
        
        vec_region_start = c(vec_region_start, grid_beta[c(F,diff(acc_vec)==1)]  )
        vec_region_end = c(grid_beta[c(diff(acc_vec)==-1, F)],vec_region_end)
        
        CI.text = paste(paste("[",vec_region_start, ",", vec_region_end, "]"),collapse = " U ")
        
        return(CI.text)
      }
    
    
}

# BASE WIDE ## -----

db_mun_wide = read_excel("db_mun_wide.xlsx")

# BASE LONG ## ----------
db_mun_long = reshape(db_mun_wide,
                      idvar = "CodIBGE",
                      varying = 6:447,
                      sep="",
                      timevar = "Year",
                      times = c(2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,
                                2014,2015,2016,2017,2018,2020),
                      new.row.names = 1:94656,
                      direction = "long")

## Agregando dados de Infrações Ambientais 

Infracoes_LONG <- read_excel("Infracoes_LONG.xlsx")

Infracoes_LONG = dplyr::rename(Infracoes_LONG, Year = "year", Fines = "Real_Value", N_Fines = "Number")


  # Juntando à base LONG 
db_mun_long = left_join(db_mun_long, Infracoes_LONG, by = c("CodIBGE","Year"), keep = FALSE)


db_mun_long$N_Fines[is.na(db_mun_long$N_Fines)] = 0 # se for NA é porque não houve multa
db_mun_long$Fines[is.na(db_mun_long$Fines)] = 0 # se for NA é porque não houve multa


 ## Convertendo importações que vieram como caracter e não número 
db_mun_long = db_mun_long %>% mutate_at(c(
                                          "DesmatPri",
                                         "GestArealIHS","GestAreal","PeCreal","PeCrealIHS","VAagroReal",
                                         "AreaPlant","POP", "ReceitaLimpa", "PIBtot", "ReceitaTrib", "ReceitaTribCotas","ReceitaTribCaPa","ReceitaTotal") 
                                       , as.numeric)

  ## Criando os Logaritmos 
db_mun_long2 = db_mun_long %>% mutate(lAreaPlant = log(AreaPlant+1),
                                      lAreaPast=log(AreaPast+1),
                                      lGarimpo=log(Garimpo+1),
                                      lPIBtot=log(PIBtot),
                                      lPOP=log(POP),
                                      lPeCreal=log(PeCreal+1),
                                      lGestAreal=log(GestAreal+1),
                                      lReceitaTotal=log(ReceitaTotal),
                                      lReceitaLimpa=log(ReceitaLimpa),
                                      lReceitaTrib =log(ReceitaTrib),
                                      lReceitaTribCotas=log(ReceitaTribCotas),
                                      lReceitaTribCaPa=log(ReceitaTribCaPa),
                                      lFines=log(N_Fines+1),
                                      lVAagro=log(VAagroReal),
                                      lDesmatPri=log(DesmatPri+1))

db_mun_long2$lPIBtot[is.nan(db_mun_long2$lPIBtot)] = NA
db_mun_long2$lPIBtot[is.infinite(db_mun_long2$lPIBtot)] = NA

db_mun_long2$lReceitaTotal[is.nan(db_mun_long2$lReceitaTotal)] = NA 
db_mun_long2$lReceitaTotal[is.infinite(db_mun_long2$lReceitaTotal)] = NA

db_mun_long2$lReceitaLimpa[is.nan(db_mun_long2$lReceitaLimpa)] = NA 
db_mun_long2$lReceitaLimpa[is.infinite(db_mun_long2$lReceitaLimpa)] = NA

db_mun_long2$lReceitaTrib[is.nan(db_mun_long2$lReceitaTrib)] = NA
db_mun_long2$lReceitaTrib[is.infinite(db_mun_long2$lReceitaTrib)] = NA

db_mun_long2$lReceitaTrib[is.nan(db_mun_long2$lReceitaTrib)] = NA
db_mun_long2$lReceitaTrib[is.infinite(db_mun_long2$lReceitaTrib)] = NA

db_mun_long2$lReceitaTribCotas[is.nan(db_mun_long2$lReceitaTribCotas)] = NA
db_mun_long2$lReceitaTribCotas[is.infinite(db_mun_long2$lReceitaTribCotas)] = NA

db_mun_long2$lReceitaTribCaPa[is.nan(db_mun_long2$lReceitaTribCaPa)] = NA
db_mun_long2$lReceitaTribCaPa[is.infinite(db_mun_long2$lReceitaTribCaPa)] = NA

db_mun_long2$VAagroReal[is.infinite(db_mun_long2$VAagroReal)] = NA 
db_mun_long2$VAagroReal[is.nan(db_mun_long2$VAagroReal)] = NA  


  ## Criando LAGs
db_mun_long2 = db_mun_long2 %>% subset(Year<2020) %>% group_by(CodIBGE) %>% mutate(lag_lPeCreal = dplyr::lag(lPeCreal), 
                                              lag_PeCrealIHS = dplyr::lag(PeCrealIHS), 
                                              lag_ConselhoLimpo = dplyr::lag(ConselhoLimpo),
                                              lag_Conselho = dplyr::lag(Conselho),
                                              lag_lGestAreal = dplyr::lag(lGestAreal), 
                                              lag_GestArealIHS = dplyr::lag(GestArealIHS), 
                                              lag_lReceitaTotal = dplyr::lag(lReceitaTotal), 
                                              lag_lReceitaLimpa = dplyr::lag(lReceitaLimpa), 
                                              lag_lReceitaTrib = dplyr::lag(lReceitaTrib),
                                              lag_lReceitaTribCotas = dplyr::lag(lReceitaTribCotas),
                                              lag_lReceitaTribCaPa = dplyr::lag(lReceitaTribCaPa),
                                              lag_lFines = dplyr::lag(lFines),
                                              lag_D_AreaProt = dplyr::lag(D_AreaProt),
                                              lag_AreaProt = dplyr::lag(AreaProt),
                                              lag_LeiFloresta = dplyr::lag(LeiFloresta),
                                              lag_MPAmazon = dplyr::lag(MPAmazon),
                                              lag_MPCerrado = dplyr::lag(MPCerrado))


### 1.0 BRASIL ######


{
  # Selecionado variávels que vão para o modelo: harmonizando numero de observações

   db_br = select(db_mun_long2,"CodIBGE","UF","Year","Nota","Remanescente",
               "DesmatPri","lDesmatPri","DesmatPriIHS",
               "lag_lPeCreal","lag_PeCrealIHS",
               "lag_lGestAreal","lag_GestArealIHS", 
               "lag_lReceitaTotal","lag_lReceitaTrib",
               "lag_ConselhoLimpo", "lag_Conselho",
               "Conselho","lAreaPlant","lAreaPast","lPOP","lPIBtot","lGarimpo", "lVAagro",
               "lag_LeiFloresta","lag_lFines","lag_D_AreaProt","lag_AreaProt","lag_AreaProt") %>% na.omit()
               

  n_distinct(db_br$CodIBGE)
  
  
  db_br = db_br %>% group_by(CodIBGE) %>% 
    mutate(DesmatMean = mean(DesmatPri), 
           DesmatStd = sd(DesmatPri),
          DesmatPriNorm = case_when(DesmatMean > 0 ~ (DesmatPri - DesmatMean)/DesmatStd)) %>%
    na.omit()
  
  db_br = panel(db_br, ~CodIBGE+Year) 

  n_distinct(db_br$CodIBGE)
  
  }


  ## Estatisticas descritivas

db_br_descr = db_br %>% as.data.frame() %>% mutate(PeC = exp(lag_lPeCreal)-1, 
                               DesmatPri = exp(lDesmatPri)-1, 
                               AreaPlant =  exp(lAreaPlant)-1,
                               AreaPast = exp(lAreaPast)-1,
                               POP = exp(lPOP),
                               Garimpo = exp(lGarimpo)-1,
                               Multas = exp(lag_lFines)-1
                              ) 

db_br_descr %>% select("DesmatPri","PeC","AreaPlant","AreaPast","POP","Garimpo","Multas",
                                           "lag_lPeCreal","lDesmatPri","lAreaPlant","lAreaPast","lPOP","lGarimpo","lag_lFines") %>% 
                                descr(order = "p",round.digits = 2) %>% View()

db_br_descr %>% select("AreaPast","lAreaPlant") %>% descr(round.digits=2)




# 1.1 IV CONSELHO IRRESTRITO ## ----------------------------------------------------------

## Primeiro Estágio 
etable(feols(lDesmatPri~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_Conselho, cluster=~CodIBGE, db_br), stage = 1, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPri="DesmatPri",
                lAreaPast = "AreaPast",                                         
                "lag_Conselho"="Lag_COMUMA",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()



# Segundos Estágios 
etable(feols(lDesmatPri~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_Conselho, cluster=~CodIBGE, db_br), stage = 2, fitstat = ~n,
       dict = c(lDesmatPri="DesmatPri",
                lAreaPast = "AreaPast",                                         
                "lag_Conselho"="Lag_COMUMA",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()


# Robustez ## --------------


### Conselho Restrito ###

etable(feols(lDesmatPri~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_ConselhoLimpo, cluster=~CodIBGE, db_br), stage = 1, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPri="DesmatPri",
                lAreaPast = "AreaPast",                                         
                "lag_ConselhoLimpo"="Lag_COMUMArestrito",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()

# Segundo Estágio 
etable(feols(lDesmatPri~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_ConselhoLimpo, cluster=~CodIBGE, db_br), stage = 2, fitstat = ~n,
       dict = c(lDesmatPri="DesmatPri",
                lAreaPast = "AreaPast",                                         
                "lag_ConselhoLimpo"="Lag_COMUMArestrito",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()


### Normalizado ###
etable(feols(DesmatPriNorm~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_Conselho, cluster=~CodIBGE, db_br), stage = 1, fitstat = ~ivwald+ivwald.p+n,
       dict = c(
                lAreaPast = "AreaPast",                                         
                "lag_Conselho"="Lag_COMUMA",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()

 # Segundo Estágio 
etable(feols(DesmatPriNorm~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_Conselho, cluster=~CodIBGE, db_br), stage = 2, fitstat = ~n,
       dict = c(
         lAreaPast = "AreaPast",                                         
         "lag_Conselho"="Lag_COMUMA",
         "lag_lPeCreal" = "Lag_PeC",
         lAreaPlant="AreaPlant",                                     
         lPOP="POP",
         lGarimpo="Garimpo",
         lag_lFines = "Lag_Multas",
         lag_D_AreaProt = "Lag_DAreaProt",
         lag_LeiFloresta = "Lag_LeiFloresta",
         CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()


### IHS no Y ### 

etable(feols(DesmatPriIHS~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_Conselho, cluster=~CodIBGE, db_br), stage = 1, fitstat = ~ivwald+ivwald.p+n,
       dict = c(
                lAreaPast = "AreaPast",                                         
                "lag_Conselho"="Lag_COMUMA",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()


etable(feols(DesmatPriIHS~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_Conselho, cluster=~CodIBGE, db_br), stage = 2, fitstat = ~ivwald+ivwald.p+n,
       dict = c(
         lAreaPast = "AreaPast",                                         
         "lag_Conselho"="Lag_COMUMA",
         "lag_lPeCreal" = "Lag_PeC",
         lAreaPlant="AreaPlant",                                     
         lPOP="POP",
         lGarimpo="Garimpo",
         lag_lFines = "Lag_Multas",
         lag_D_AreaProt = "Lag_DAreaProt",
         lag_LeiFloresta = "Lag_LeiFloresta",
         CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()


## Forma reduzida 

etable(feols(lDesmatPri~csw(lag_Conselho,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta)|
               CodIBGE+UF^Year, db_br), fitstat = ~ivwald+ivwald.p+n,
       dict = c(
         lAreaPast = "AreaPast",                                         
         
         "lag_lPeCreal" = "Lag_PeC",
         lAreaPlant="AreaPlant",                                     
         lPOP="POP",
         lGarimpo="Garimpo",
         lag_lFines = "Lag_Multas",
         lag_D_AreaProt = "Lag_DAreaProt",
         lag_LeiFloresta = "Lag_LeiFloresta",
         CodIBGE= "Municipio")) %>% View()


  # OLS
etable(feols(lDesmatPri~csw(lag_lPeCreal,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta)|
               CodIBGE+UF^Year, cluster=~CodIBGE, db_br))




# Exercício Extra: Gestão Ambiental--------------------------------------------

etable(feols(lDesmatPri~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta)|
               CodIBGE+UF^Year|lag_lGestAreal~lag_Conselho, cluster=~CodIBGE, db_br), stage = 1, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPri="DesmatPri",
                lAreaPast = "AreaPast",                                         
                "lag_Conselho"="Lag_COMUMA",
                "lag_lGestAreal" = "Lag_GestA",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()

# Segundo Estágio

etable(feols(lDesmatPri~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta)|
               CodIBGE+UF^Year|lag_lGestAreal~lag_Conselho, cluster=~CodIBGE, db_br), stage = 2, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPri="DesmatPri",
                lAreaPast = "AreaPast",                                         
                "lag_Conselho"="Lag_COMUMA",
                "lag_lGestAreal" = "Lag_GestA",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()

### 2.0 AMAZÔNIA - MAPBIOMAS #####----------------------------------

desmatAmazon = read_excel("DesmatAmazonia.xlsx")
amazon_long = reshape(desmatAmazon, idvar = "CodIBGE",
                       varying = 3:70,
                       sep="",
                       timevar = "Year",
                       times = c(2004:2020),
                       new.row.names = 1:9586,
                       direction = "long")

n_distinct(amazon_long$CodIBGE)

amazon_long= amazon_long %>% mutate_at(c("DesmatINPE","DesmatINPEihs"), 
                                       as.numeric) %>% na.omit() %>% subset(Year<2020)
  # já removendo NAs para deixar base igual a do INPE

n_distinct(amazon_long$CodIBGE)

db_amazon_long = right_join(db_mun_long2, amazon_long, by = c("CodIBGE","Year"), keep = FALSE)

  ## Incluindo Logs Adicionais para o desmatamento exclusivo do Bioma Amazônico no Mapbiomas (sufixo AM) e INPE

db_amazon_long = db_amazon_long %>% mutate(lDesmatPriAM=log(DesmatPriAM+1),
                                          lDesmatINPE=log(DesmatINPE+1))
                                          

## Convertendo em Painel 
{
  db_amazon = select(db_amazon_long,"CodIBGE","UF","Year","Nota","RemanescenteAM","DesmatPriAM","lDesmatPriAM","DesmatPriAMihs","DesmatINPE","lDesmatINPE","DesmatINPEihs",
                   "lag_lPeCreal","lag_PeCrealIHS","lag_lGestAreal","lag_GestArealIHS",
                   "lag_lReceitaTrib","lag_Conselho", "lag_ConselhoLimpo",
                   "Conselho", "lag_Conselho", "lVAagro", "lag_lFines", "lag_D_AreaProt", "lag_lReceitaTotal",
                   "lAreaPlant","lAreaPast","lPOP", "lPIBtot","lGarimpo","lag_LeiFloresta","lag_MPAmazon") %>% na.omit()
  
  n_distinct(db_amazon$CodIBGE)
  
  db_amazon = db_amazon %>% group_by(CodIBGE) %>% 
    mutate(DesmatMean = mean(DesmatPriAM), 
           DesmatStd = sd(DesmatPriAM),
           DesmatPriAMn = case_when(DesmatMean > 0 ~ (DesmatPriAM - DesmatMean)/DesmatStd),
           DesmatMeanINPE = mean(DesmatINPE), 
           DesmatStdINPE = sd(DesmatINPE),
           DesmatINPEn = case_when(DesmatMeanINPE > 0 ~ (DesmatINPE - DesmatMeanINPE)/DesmatStdINPE)) %>%
    na.omit()
  
  n_distinct(db_amazon$CodIBGE)
  
  db_amazon = panel(db_amazon, ~CodIBGE+Year) 
}



## Estatisticas descritivas


db_amazon_descr = db_amazon %>% as.data.frame() %>% mutate(PeC = exp(lag_lPeCreal)-1, 
                                                   DesmatPri = exp(lDesmatPriAM)-1, 
                                                   DesmatINPE=exp(lDesmatINPE)-1,
                                                   AreaPlant =  exp(lAreaPlant)-1,
                                                   AreaPast = exp(lAreaPast)-1,
                                                   POP = exp(lPOP),
                                                   Garimpo = exp(lGarimpo)-1,
                                                   Multas = exp(lag_lFines)-1
                                                   ) 


  select("DesmatPriAM", "DesmatINPE","PeC","AreaPlant","AreaPast","POP","Garimpo","Multas", 
            "lDesmatPriAM","lDesmatINPE","lag_lPeCreal","lAreaPlant","lAreaPast","lPOP","lGarimpo","lag_lFines") %>% 
  descr(order = "p",round.digits = 2)

db_amazon_descr %>% select("AreaPast","lAreaPlant") %>% descr(round.digits=2)


# 2.1 IV - CONSELHO IRRESTRITO ----

## Primeiro Estágio 
etable(feols(lDesmatPriAM~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta,lag_MPAmazon)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_Conselho, cluster=~CodIBGE, db_amazon), stage = 1, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPri="DesmatPriAM",
                lAreaPast = "AreaPast",                                         
                "lag_Conselho"="Lag_COMUMA",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()



{
  y = feols(lDesmatPriAM~1|CodIBGE+UF^Year, db_amazon)
  y_resid= y$residuals
  
  x_lagendog = feols(lag_lPeCreal~1|CodIBGE+UF^Year, db_amazon)
  x_lagendog_resid= x_lagendog$residuals
  
  z_laginstrument = feols(lag_Conselho~1|CodIBGE+UF^Year, db_amazon)
  z_laginstrument_resid = z_laginstrument$residuals
  
  x_areapast = feols(lAreaPast~1|CodIBGE+UF^Year, db_amazon)
  x_areapast_resid =  x_areapast$residuals
  
  x_garimpo = feols(lGarimpo~1|CodIBGE+UF^Year, db_amazon)
  x_garimpo_resid =  x_garimpo$residuals
  
  x_vaagro = feols(lVAagro~1|CodIBGE+UF^Year, db_amazon)
  x_vaagro_resid =  x_vaagro$residuals
  
  x_areaplant = feols(lAreaPlant~1|CodIBGE+UF^Year, db_amazon)
  x_areaplant_resid =  x_areaplant$residuals
  
  x_pop = feols(lPOP~1|CodIBGE+UF^Year, db_amazon)
  x_pop_resid =  x_pop$residuals
  
  x_conselho = feols(lag_Conselho~1|CodIBGE+UF^Year, db_amazon)
  x_conselho_resid =  x_conselho$residuals
  
  x_lei = feols(lag_LeiFloresta~1|CodIBGE+UF^Year, db_amazon)
  x_lei_resid =  x_lei$residuals
  
  x_areaprot = feols(lag_D_AreaProt~1|CodIBGE+UF^Year, db_amazon)
  x_areaprot_resid =  x_areaprot$residuals
  
  x_multa = feols(lag_lFines~1|CodIBGE+UF^Year, db_amazon)
  x_multa_resid =  x_multa$residuals
  
  x_MPAmazon = feols(lag_MPAmazon~1|CodIBGE+UF^Year, db_amazon)
  x_MPAmazon_resid =  x_MPAmazon$residuals
  
  CodIBGE = db_amazon$CodIBGE
  
  ## data frame para inputar no teste
  db_AR_Test = data.frame(y_resid,
                          x_lagendog_resid,
                          z_laginstrument_resid,
                          x_areapast_resid,
                          x_areaplant_resid,
                          x_vaagro_resid,
                          x_pop_resid,
                          x_garimpo_resid,
                          x_conselho_resid,
                          x_lei_resid,
                          x_areaprot_resid,
                          x_multa_resid,
                          x_MPAmazon_resid,
                          CodIBGE)
  
  
  # Nada
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid"
                               
                                             
                                             
                                , vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid"
                            
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  
  
  
  
  # AreaPlant
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  # AreaPlant, AreaPast 
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  
  # AreaPlant, AreasPast, POP
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                               endogenous = "x_lagendog_resid",
                               instruments = "z_laginstrument_resid",
                               controls = c("x_areaplant_resid",
                                            "x_areapast_resid",
                                            "x_pop_resid"
                                            
                                            
                               ), vcovCL,
                               data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  
  
  # # AreaPlant, AreasPast, POP, Garimpo
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  #AreaPlant, AreasPast, POP, Garimpo, Multa 
  
  
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                         
                                             "x_multa_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                        
                                         "x_multa_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  # AreaPlant, AreasPast, POP, Garimpo, Multa, AreaProt
  
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                   
                                             "x_multa_resid",
                                             "x_areaprot_resid"
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                  
                                         "x_multa_resid",
                                         "x_areaprot_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  
  # AreaPlant, AreasPast, POP, Garimpo, Multa, AreaProt, Lei
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                          
                                             "x_multa_resid",
                                             "x_areaprot_resid",
                                             "x_lei_resid"
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                    
                                         "x_multa_resid",
                                         "x_areaprot_resid",
                                         "x_lei_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  
  # AreaPlant, AreasPast, POP, Garimpo, Multa, AreaProt, MPAmazon
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                            
                                             "x_multa_resid",
                                             "x_areaprot_resid",
                                             "x_lei_resid",
                                             "x_MPAmazon_resid"
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                       
                                         "x_multa_resid",
                                         "x_areaprot_resid",
                                         "x_lei_resid",
                                         "x_MPAmazon_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
       }


etable(feols(lDesmatPriAM~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta,lag_MPAmazon)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_Conselho, cluster=~CodIBGE, db_amazon), stage = 2, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPriAM="DesmatPriAM",
                lAreaPast = "AreaPast",                                         
                "lag_Conselho"="Lag_COMUMA",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()





# Robustez------------------------------------------------------------ 


## Conselho Restrito ## 

## Primeiro Estágio 
etable(feols(lDesmatPriAM~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta,lag_MPAmazon)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_ConselhoLimpo, cluster=~CodIBGE, db_amazon), stage = 1, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPri="DesmatPriAM",
                lAreaPast = "AreaPast",                                         
                "lag_ConselhoLimpo"="Lag_COMUMArestrito",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()


{
  y = feols(lDesmatPriAM~1|CodIBGE+UF^Year, db_amazon)
  y_resid= y$residuals
  
  x_lagendog = feols(lag_lPeCreal~1|CodIBGE+UF^Year, db_amazon)
  x_lagendog_resid= x_lagendog$residuals
  
  z_laginstrument = feols(lag_ConselhoLimpo~1|CodIBGE+UF^Year, db_amazon)
  z_laginstrument_resid = z_laginstrument$residuals
  
  x_areapast = feols(lAreaPast~1|CodIBGE+UF^Year, db_amazon)
  x_areapast_resid =  x_areapast$residuals
  
  x_garimpo = feols(lGarimpo~1|CodIBGE+UF^Year, db_amazon)
  x_garimpo_resid =  x_garimpo$residuals
  
  x_vaagro = feols(lVAagro~1|CodIBGE+UF^Year, db_amazon)
  x_vaagro_resid =  x_vaagro$residuals
  
  x_areaplant = feols(lAreaPlant~1|CodIBGE+UF^Year, db_amazon)
  x_areaplant_resid =  x_areaplant$residuals
  
  x_pop = feols(lPOP~1|CodIBGE+UF^Year, db_amazon)
  x_pop_resid =  x_pop$residuals
  
  x_conselho = feols(lag_Conselho~1|CodIBGE+UF^Year, db_amazon)
  x_conselho_resid =  x_conselho$residuals
  
  x_lei = feols(lag_LeiFloresta~1|CodIBGE+UF^Year, db_amazon)
  x_lei_resid =  x_lei$residuals
  
  x_areaprot = feols(lag_D_AreaProt~1|CodIBGE+UF^Year, db_amazon)
  x_areaprot_resid =  x_areaprot$residuals
  
  x_multa = feols(lag_lFines~1|CodIBGE+UF^Year, db_amazon)
  x_multa_resid =  x_multa$residuals
  
  x_MPAmazon = feols(lag_MPAmazon~1|CodIBGE+UF^Year, db_amazon)
  x_MPAmazon_resid =  x_MPAmazon$residuals
  
  CodIBGE = db_amazon$CodIBGE
  
  ## data frame para inputar no teste
  db_AR_Test = data.frame(y_resid,
                          x_lagendog_resid,
                          z_laginstrument_resid,
                          x_areapast_resid,
                          x_areaplant_resid,
                          x_vaagro_resid,
                          x_pop_resid,
                          x_garimpo_resid,
                          x_conselho_resid,
                          x_lei_resid,
                          x_areaprot_resid,
                          x_multa_resid,
                          x_MPAmazon_resid,
                          CodIBGE)
  
  # Nada
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid"
                             
                                             
                                             
                                , vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                           
                            
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  
  
  
  
  
  # AreaPlant
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  # AreaPlant, AreaPast 
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  
  # AreaPlant, AreasPast, POP
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  
  
  # # AreaPlant, AreasPast, POP, Garimpo
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  #AreaPlant, AreasPast, POP, Garimpo, Multa 
  
  
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                             
                                             "x_multa_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                         
                                         "x_multa_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  # AreaPlant, AreasPast, POP, Garimpo, Multa, AreaProt
  
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                             
                                             "x_multa_resid",
                                             "x_areaprot_resid"
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                         
                                         "x_multa_resid",
                                         "x_areaprot_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  
  # AreaPlant, AreasPast, POP, Garimpo, Multa, AreaProt, Lei
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                             
                                             "x_multa_resid",
                                             "x_areaprot_resid",
                                             "x_lei_resid"
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                         
                                         "x_multa_resid",
                                         "x_areaprot_resid",
                                         "x_lei_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  
  # AreaPlant, AreasPast, POP, Garimpo, Multa, AreaProt, MPAmazon
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                             
                                             "x_multa_resid",
                                             "x_areaprot_resid",
                                             "x_lei_resid",
                                             "x_MPAmazon_resid"
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                         
                                         "x_multa_resid",
                                         "x_areaprot_resid",
                                         "x_lei_resid",
                                         "x_MPAmazon_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
       }

etable(feols(lDesmatPriAM~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta,lag_MPAmazon)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_ConselhoLimpo, cluster=~CodIBGE, db_amazon), stage = 2, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPri="DesmatPriAM",
                lAreaPast = "AreaPast",                                         
                "lag_ConselhoLimpo"="Lag_COMUMArestrito",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()



### Normalizado ###

etable(feols(DesmatPriAMn~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta,lag_MPAmazon)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_Conselho, cluster=~CodIBGE, db_amazon), stage = 1, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPri="DesmatPriAM",
                lAreaPast = "AreaPast",                                         
                "lag_Conselho"="Lag_COMUMA",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()


{
  y = feols(DesmatPriAMn~1|CodIBGE+UF^Year, db_amazon)
  y_resid= y$residuals
  
  x_lagendog = feols(lag_lPeCreal~1|CodIBGE+UF^Year, db_amazon)
  x_lagendog_resid= x_lagendog$residuals
  
  z_laginstrument = feols(lag_Conselho~1|CodIBGE+UF^Year, db_amazon)
  z_laginstrument_resid = z_laginstrument$residuals
  
  x_areapast = feols(lAreaPast~1|CodIBGE+UF^Year, db_amazon)
  x_areapast_resid =  x_areapast$residuals
  
  x_garimpo = feols(lGarimpo~1|CodIBGE+UF^Year, db_amazon)
  x_garimpo_resid =  x_garimpo$residuals
  
  x_vaagro = feols(lVAagro~1|CodIBGE+UF^Year, db_amazon)
  x_vaagro_resid =  x_vaagro$residuals
  
  x_areaplant = feols(lAreaPlant~1|CodIBGE+UF^Year, db_amazon)
  x_areaplant_resid =  x_areaplant$residuals
  
  x_pop = feols(lPOP~1|CodIBGE+UF^Year, db_amazon)
  x_pop_resid =  x_pop$residuals
  
  x_conselho = feols(lag_Conselho~1|CodIBGE+UF^Year, db_amazon)
  x_conselho_resid =  x_conselho$residuals
  
  x_lei = feols(lag_LeiFloresta~1|CodIBGE+UF^Year, db_amazon)
  x_lei_resid =  x_lei$residuals
  
  x_areaprot = feols(lag_D_AreaProt~1|CodIBGE+UF^Year, db_amazon)
  x_areaprot_resid =  x_areaprot$residuals
  
  x_multa = feols(lag_lFines~1|CodIBGE+UF^Year, db_amazon)
  x_multa_resid =  x_multa$residuals
  
  x_MPAmazon = feols(lag_MPAmazon~1|CodIBGE+UF^Year, db_amazon)
  x_MPAmazon_resid =  x_MPAmazon$residuals
  
  CodIBGE = db_amazon$CodIBGE
  
  ## data frame para inputar no teste
  db_AR_Test = data.frame(y_resid,
                          x_lagendog_resid,
                          z_laginstrument_resid,
                          x_areapast_resid,
                          x_areaplant_resid,
                          x_vaagro_resid,
                          x_pop_resid,
                          x_garimpo_resid,
                          x_conselho_resid,
                          x_lei_resid,
                          x_areaprot_resid,
                          x_multa_resid,
                          x_MPAmazon_resid,
                          CodIBGE)
  
  
  # AreaPlant, AreasPast, POP, Garimpo
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  # AreaPlant, AreasPast, POP, Garimpo, Multa, AreaProt, MPAmazon
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                             
                                             "x_multa_resid",
                                             "x_areaprot_resid",
                                             "x_lei_resid",
                                             "x_MPAmazon_resid"
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                         
                                         "x_multa_resid",
                                         "x_areaprot_resid",
                                         "x_lei_resid",
                                         "x_MPAmazon_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
       }

etable(feols(DesmatPriAMn~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta,lag_MPAmazon)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_Conselho, cluster=~CodIBGE, db_amazon), stage = 2, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPri="DesmatPriAM",
                lAreaPast = "AreaPast",                                         
                "lag_ConselhoLimpo"="Lag_COMUMArestrito",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()


### IHS no Y ###

etable(feols(DesmatPriAMihs~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta,lag_MPAmazon)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_Conselho, cluster=~CodIBGE, db_amazon), stage = 1, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPri="DesmatPriAM",
                lAreaPast = "AreaPast",                                         
                "lag_Conselho"="Lag_COMUMA",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()


{
  y = feols(DesmatPriAMihs~1|CodIBGE+UF^Year, db_amazon)
  y_resid= y$residuals
  
  x_lagendog = feols(lag_lPeCreal~1|CodIBGE+UF^Year, db_amazon)
  x_lagendog_resid= x_lagendog$residuals
  
  z_laginstrument = feols(lag_Conselho~1|CodIBGE+UF^Year, db_amazon)
  z_laginstrument_resid = z_laginstrument$residuals
  
  x_areapast = feols(lAreaPast~1|CodIBGE+UF^Year, db_amazon)
  x_areapast_resid =  x_areapast$residuals
  
  x_garimpo = feols(lGarimpo~1|CodIBGE+UF^Year, db_amazon)
  x_garimpo_resid =  x_garimpo$residuals
  
  x_vaagro = feols(lVAagro~1|CodIBGE+UF^Year, db_amazon)
  x_vaagro_resid =  x_vaagro$residuals
  
  x_areaplant = feols(lAreaPlant~1|CodIBGE+UF^Year, db_amazon)
  x_areaplant_resid =  x_areaplant$residuals
  
  x_pop = feols(lPOP~1|CodIBGE+UF^Year, db_amazon)
  x_pop_resid =  x_pop$residuals
  
  x_conselho = feols(lag_Conselho~1|CodIBGE+UF^Year, db_amazon)
  x_conselho_resid =  x_conselho$residuals
  
  x_lei = feols(lag_LeiFloresta~1|CodIBGE+UF^Year, db_amazon)
  x_lei_resid =  x_lei$residuals
  
  x_areaprot = feols(lag_D_AreaProt~1|CodIBGE+UF^Year, db_amazon)
  x_areaprot_resid =  x_areaprot$residuals
  
  x_multa = feols(lag_lFines~1|CodIBGE+UF^Year, db_amazon)
  x_multa_resid =  x_multa$residuals
  
  x_MPAmazon = feols(lag_MPAmazon~1|CodIBGE+UF^Year, db_amazon)
  x_MPAmazon_resid =  x_MPAmazon$residuals
  
  CodIBGE = db_amazon$CodIBGE
  
  ## data frame para inputar no teste
  db_AR_Test = data.frame(y_resid,
                          x_lagendog_resid,
                          z_laginstrument_resid,
                          x_areapast_resid,
                          x_areaplant_resid,
                          x_vaagro_resid,
                          x_pop_resid,
                          x_garimpo_resid,
                          x_conselho_resid,
                          x_lei_resid,
                          x_areaprot_resid,
                          x_multa_resid,
                          x_MPAmazon_resid,
                          CodIBGE)
  
  
  # AreaPlant, AreasPast, POP, Garimpo
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  # AreaPlant, AreasPast, POP, Garimpo, Multa, AreaProt, MPAmazon
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                             
                                             "x_multa_resid",
                                             "x_areaprot_resid",
                                             "x_lei_resid",
                                             "x_MPAmazon_resid"
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                         
                                         "x_multa_resid",
                                         "x_areaprot_resid",
                                         "x_lei_resid",
                                         "x_MPAmazon_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
       }

etable(feols(DesmatPriAMihs~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta,lag_MPAmazon)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_Conselho, cluster=~CodIBGE, db_amazon), stage = 2, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPri="DesmatPriAM",
                lAreaPast = "AreaPast",                                         
                "lag_ConselhoLimpo"="Lag_COMUMArestrito",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()


## OLS 

etable(feols(lDesmatPriAM~csw(lag_lPeCreal,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta,lag_MPAmazon)|
              CodIBGE+UF^Year, cluster=~CodIBGE, db_amazon)) %>% View()


## Forma reduzida  

etable(feols(lDesmatPriAM~csw(lag_Conselho,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta,lag_MPAmazon)|
               CodIBGE+UF^Year, cluster=~CodIBGE, db_amazon)) %>% View()






# Exercício Extra: Gestão Ambiental ----------------------------------------------------


etable(feols(lDesmatPriAM~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta,lag_MPAmazon)|
               CodIBGE+UF^Year|lag_lGestAreal~lag_Conselho, cluster=~CodIBGE, db_amazon), stage = 1, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPri="DesmatPriAM",
                lAreaPast = "AreaPast",                                         
                "lag_Conselho"="Lag_COMUMA",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()


{
  y = feols(lDesmatPriAM~1|CodIBGE+UF^Year, db_amazon)
  y_resid= y$residuals
  
  x_lagendog = feols(lag_lGestAreal~1|CodIBGE+UF^Year, db_amazon)
  x_lagendog_resid= x_lagendog$residuals
  
  z_laginstrument = feols(lag_Conselho~1|CodIBGE+UF^Year, db_amazon)
  z_laginstrument_resid = z_laginstrument$residuals
  
  x_areapast = feols(lAreaPast~1|CodIBGE+UF^Year, db_amazon)
  x_areapast_resid =  x_areapast$residuals
  
  x_garimpo = feols(lGarimpo~1|CodIBGE+UF^Year, db_amazon)
  x_garimpo_resid =  x_garimpo$residuals
  
  x_vaagro = feols(lVAagro~1|CodIBGE+UF^Year, db_amazon)
  x_vaagro_resid =  x_vaagro$residuals
  
  x_areaplant = feols(lAreaPlant~1|CodIBGE+UF^Year, db_amazon)
  x_areaplant_resid =  x_areaplant$residuals
  
  x_pop = feols(lPOP~1|CodIBGE+UF^Year, db_amazon)
  x_pop_resid =  x_pop$residuals
  
  x_conselho = feols(lag_Conselho~1|CodIBGE+UF^Year, db_amazon)
  x_conselho_resid =  x_conselho$residuals
  
  x_lei = feols(lag_LeiFloresta~1|CodIBGE+UF^Year, db_amazon)
  x_lei_resid =  x_lei$residuals
  
  x_areaprot = feols(lag_D_AreaProt~1|CodIBGE+UF^Year, db_amazon)
  x_areaprot_resid =  x_areaprot$residuals
  
  x_multa = feols(lag_lFines~1|CodIBGE+UF^Year, db_amazon)
  x_multa_resid =  x_multa$residuals
  
  x_MPAmazon = feols(lag_MPAmazon~1|CodIBGE+UF^Year, db_amazon)
  x_MPAmazon_resid =  x_MPAmazon$residuals
  
  CodIBGE = db_amazon$CodIBGE
  
  ## data frame para inputar no teste
  db_AR_Test = data.frame(y_resid,
                          x_lagendog_resid,
                          z_laginstrument_resid,
                          x_areapast_resid,
                          x_areaplant_resid,
                          x_vaagro_resid,
                          x_pop_resid,
                          x_garimpo_resid,
                          x_conselho_resid,
                          x_lei_resid,
                          x_areaprot_resid,
                          x_multa_resid,
                          x_MPAmazon_resid,
                          CodIBGE)
  
  
  # AreaPlant, AreasPast, POP, Garimpo
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  # AreaPlant, AreasPast, POP, Garimpo, Multa, AreaProt, MPAmazon
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                             
                                             "x_multa_resid",
                                             "x_areaprot_resid",
                                             "x_lei_resid",
                                             "x_MPAmazon_resid"
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                         
                                         "x_multa_resid",
                                         "x_areaprot_resid",
                                         "x_lei_resid",
                                         "x_MPAmazon_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
       }

etable(feols(lDesmatPriAM~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta,lag_MPAmazon)|
               CodIBGE+UF^Year|lag_lGestAreal~lag_Conselho, cluster=~CodIBGE, db_amazon), stage = 2, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPri="DesmatPriAM",
                lAreaPast = "AreaPast",                                         
                "lag_ConselhoLimpo"="Lag_COMUMA",
                "lag_lGestAreal" = "Lag_GestA",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()




### 3.0 AMAZÔNIA - INPE  #####
# 3.1 IV - CONSELHO IRRESTRITO ----

## Primeiro Estágio 
etable(feols(lDesmatINPE~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta,lag_MPAmazon)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_Conselho, cluster=~CodIBGE, db_amazon), stage = 1, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPriINPE="DesmatINPE",
                lAreaPast = "AreaPast",                                         
                "lag_Conselho"="Lag_COMUMA",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()



{
  y = feols(lDesmatINPE~1|CodIBGE+UF^Year, db_amazon)
  y_resid= y$residuals
  
  x_lagendog = feols(lag_lPeCreal~1|CodIBGE+UF^Year, db_amazon)
  x_lagendog_resid= x_lagendog$residuals
  
  z_laginstrument = feols(lag_Conselho~1|CodIBGE+UF^Year, db_amazon)
  z_laginstrument_resid = z_laginstrument$residuals
  
  x_areapast = feols(lAreaPast~1|CodIBGE+UF^Year, db_amazon)
  x_areapast_resid =  x_areapast$residuals
  
  x_garimpo = feols(lGarimpo~1|CodIBGE+UF^Year, db_amazon)
  x_garimpo_resid =  x_garimpo$residuals
  
  x_vaagro = feols(lVAagro~1|CodIBGE+UF^Year, db_amazon)
  x_vaagro_resid =  x_vaagro$residuals
  
  x_areaplant = feols(lAreaPlant~1|CodIBGE+UF^Year, db_amazon)
  x_areaplant_resid =  x_areaplant$residuals
  
  x_pop = feols(lPOP~1|CodIBGE+UF^Year, db_amazon)
  x_pop_resid =  x_pop$residuals
  
  x_conselho = feols(lag_Conselho~1|CodIBGE+UF^Year, db_amazon)
  x_conselho_resid =  x_conselho$residuals
  
  x_lei = feols(lag_LeiFloresta~1|CodIBGE+UF^Year, db_amazon)
  x_lei_resid =  x_lei$residuals
  
  x_areaprot = feols(lag_D_AreaProt~1|CodIBGE+UF^Year, db_amazon)
  x_areaprot_resid =  x_areaprot$residuals
  
  x_multa = feols(lag_lFines~1|CodIBGE+UF^Year, db_amazon)
  x_multa_resid =  x_multa$residuals
  
  x_MPAmazon = feols(lag_MPAmazon~1|CodIBGE+UF^Year, db_amazon)
  x_MPAmazon_resid =  x_MPAmazon$residuals
  
  CodIBGE = db_amazon$CodIBGE
  
  ## data frame para inputar no teste
  db_AR_Test = data.frame(y_resid,
                          x_lagendog_resid,
                          z_laginstrument_resid,
                          x_areapast_resid,
                          x_areaplant_resid,
                          x_vaagro_resid,
                          x_pop_resid,
                          x_garimpo_resid,
                          x_conselho_resid,
                          x_lei_resid,
                          x_areaprot_resid,
                          x_multa_resid,
                          x_MPAmazon_resid,
                          CodIBGE)
  
  
  # NAda
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  
  
  
  # AreaPlant
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  # AreaPlant, AreaPast 
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  
  # AreaPlant, AreasPast, POP
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  
  
  # # AreaPlant, AreasPast, POP, Garimpo
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  #AreaPlant, AreasPast, POP, Garimpo, Multa 
  
  
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                             
                                             "x_multa_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                         
                                         "x_multa_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  # AreaPlant, AreasPast, POP, Garimpo, Multa, AreaProt
  
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                             
                                             "x_multa_resid",
                                             "x_areaprot_resid"
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                         
                                         "x_multa_resid",
                                         "x_areaprot_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  
  # AreaPlant, AreasPast, POP, Garimpo, Multa, AreaProt, Lei
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                             
                                             "x_multa_resid",
                                             "x_areaprot_resid",
                                             "x_lei_resid"
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                         
                                         "x_multa_resid",
                                         "x_areaprot_resid",
                                         "x_lei_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  
  # AreaPlant, AreasPast, POP, Garimpo, Multa, AreaProt, MPAmazon
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                             
                                             "x_multa_resid",
                                             "x_areaprot_resid",
                                             "x_lei_resid",
                                             "x_MPAmazon_resid"
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                         
                                         "x_multa_resid",
                                         "x_areaprot_resid",
                                         "x_lei_resid",
                                         "x_MPAmazon_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
       }


etable(feols(lDesmatINPE~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta,lag_MPAmazon)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_Conselho, cluster=~CodIBGE, db_amazon), stage = 2, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPriINPE="DesmatINPE",
                lAreaPast = "AreaPast",                                         
                "lag_Conselho"="Lag_COMUMA",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()





# Robustez------------------------------------------------------------ 


## Conselho Restrito ## 

## Primeiro Estágio 
etable(feols(lDesmatINPE~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta,lag_MPAmazon)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_ConselhoLimpo, cluster=~CodIBGE, db_amazon), stage = 1, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPriINPE="DesmatINPE",
                lAreaPast = "AreaPast",                                         
                "lag_ConselhoLimpo"="Lag_COMUMArestrito",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()


{
  y = feols(lDesmatINPE~1|CodIBGE+UF^Year, db_amazon)
  y_resid= y$residuals
  
  x_lagendog = feols(lag_lPeCreal~1|CodIBGE+UF^Year, db_amazon)
  x_lagendog_resid= x_lagendog$residuals
  
  z_laginstrument = feols(lag_ConselhoLimpo~1|CodIBGE+UF^Year, db_amazon)
  z_laginstrument_resid = z_laginstrument$residuals
  
  x_areapast = feols(lAreaPast~1|CodIBGE+UF^Year, db_amazon)
  x_areapast_resid =  x_areapast$residuals
  
  x_garimpo = feols(lGarimpo~1|CodIBGE+UF^Year, db_amazon)
  x_garimpo_resid =  x_garimpo$residuals
  
  x_vaagro = feols(lVAagro~1|CodIBGE+UF^Year, db_amazon)
  x_vaagro_resid =  x_vaagro$residuals
  
  x_areaplant = feols(lAreaPlant~1|CodIBGE+UF^Year, db_amazon)
  x_areaplant_resid =  x_areaplant$residuals
  
  x_pop = feols(lPOP~1|CodIBGE+UF^Year, db_amazon)
  x_pop_resid =  x_pop$residuals
  
  x_conselho = feols(lag_Conselho~1|CodIBGE+UF^Year, db_amazon)
  x_conselho_resid =  x_conselho$residuals
  
  x_lei = feols(lag_LeiFloresta~1|CodIBGE+UF^Year, db_amazon)
  x_lei_resid =  x_lei$residuals
  
  x_areaprot = feols(lag_D_AreaProt~1|CodIBGE+UF^Year, db_amazon)
  x_areaprot_resid =  x_areaprot$residuals
  
  x_multa = feols(lag_lFines~1|CodIBGE+UF^Year, db_amazon)
  x_multa_resid =  x_multa$residuals
  
  x_MPAmazon = feols(lag_MPAmazon~1|CodIBGE+UF^Year, db_amazon)
  x_MPAmazon_resid =  x_MPAmazon$residuals
  
  CodIBGE = db_amazon$CodIBGE
  
  ## data frame para inputar no teste
  db_AR_Test = data.frame(y_resid,
                          x_lagendog_resid,
                          z_laginstrument_resid,
                          x_areapast_resid,
                          x_areaplant_resid,
                          x_vaagro_resid,
                          x_pop_resid,
                          x_garimpo_resid,
                          x_conselho_resid,
                          x_lei_resid,
                          x_areaprot_resid,
                          x_multa_resid,
                          x_MPAmazon_resid,
                          CodIBGE)
  
  
  # Nada
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                              
                                             
                                             
                                vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                     
                          
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  
  
  
  # AreaPlant
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  # AreaPlant, AreaPast 
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  
  # AreaPlant, AreasPast, POP
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  
  
  # # AreaPlant, AreasPast, POP, Garimpo
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  #AreaPlant, AreasPast, POP, Garimpo, Multa 
  
  
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                             
                                             "x_multa_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                         
                                         "x_multa_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  # AreaPlant, AreasPast, POP, Garimpo, Multa, AreaProt
  
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                             
                                             "x_multa_resid",
                                             "x_areaprot_resid"
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                         
                                         "x_multa_resid",
                                         "x_areaprot_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  
  # AreaPlant, AreasPast, POP, Garimpo, Multa, AreaProt, Lei
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                             
                                             "x_multa_resid",
                                             "x_areaprot_resid",
                                             "x_lei_resid"
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                         
                                         "x_multa_resid",
                                         "x_areaprot_resid",
                                         "x_lei_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  
  # AreaPlant, AreasPast, POP, Garimpo, Multa, AreaProt, MPAmazon
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                             
                                             "x_multa_resid",
                                             "x_areaprot_resid",
                                             "x_lei_resid",
                                             "x_MPAmazon_resid"
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                         
                                         "x_multa_resid",
                                         "x_areaprot_resid",
                                         "x_lei_resid",
                                         "x_MPAmazon_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
       }

etable(feols(lDesmatINPE~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta,lag_MPAmazon)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_ConselhoLimpo, cluster=~CodIBGE, db_amazon), stage = 2, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPriINPE="DesmatINPE",
                lAreaPast = "AreaPast",                                         
                "lag_ConselhoLimpo"="Lag_COMUMArestrito",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()



### Normalizado ###

etable(feols(DesmatINPEn~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta,lag_MPAmazon)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_Conselho, cluster=~CodIBGE, db_amazon), stage = 1, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPri="DesmatINPE",
                lAreaPast = "AreaPast",                                         
                "lag_Conselho"="Lag_COMUMA",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()


{
  y = feols(DesmatINPEn~1|CodIBGE+UF^Year, db_amazon)
  y_resid= y$residuals
  
  x_lagendog = feols(lag_lPeCreal~1|CodIBGE+UF^Year, db_amazon)
  x_lagendog_resid= x_lagendog$residuals
  
  z_laginstrument = feols(lag_Conselho~1|CodIBGE+UF^Year, db_amazon)
  z_laginstrument_resid = z_laginstrument$residuals
  
  x_areapast = feols(lAreaPast~1|CodIBGE+UF^Year, db_amazon)
  x_areapast_resid =  x_areapast$residuals
  
  x_garimpo = feols(lGarimpo~1|CodIBGE+UF^Year, db_amazon)
  x_garimpo_resid =  x_garimpo$residuals
  
  x_vaagro = feols(lVAagro~1|CodIBGE+UF^Year, db_amazon)
  x_vaagro_resid =  x_vaagro$residuals
  
  x_areaplant = feols(lAreaPlant~1|CodIBGE+UF^Year, db_amazon)
  x_areaplant_resid =  x_areaplant$residuals
  
  x_pop = feols(lPOP~1|CodIBGE+UF^Year, db_amazon)
  x_pop_resid =  x_pop$residuals
  
  x_conselho = feols(lag_Conselho~1|CodIBGE+UF^Year, db_amazon)
  x_conselho_resid =  x_conselho$residuals
  
  x_lei = feols(lag_LeiFloresta~1|CodIBGE+UF^Year, db_amazon)
  x_lei_resid =  x_lei$residuals
  
  x_areaprot = feols(lag_D_AreaProt~1|CodIBGE+UF^Year, db_amazon)
  x_areaprot_resid =  x_areaprot$residuals
  
  x_multa = feols(lag_lFines~1|CodIBGE+UF^Year, db_amazon)
  x_multa_resid =  x_multa$residuals
  
  x_MPAmazon = feols(lag_MPAmazon~1|CodIBGE+UF^Year, db_amazon)
  x_MPAmazon_resid =  x_MPAmazon$residuals
  
  CodIBGE = db_amazon$CodIBGE
  
  ## data frame para inputar no teste
  db_AR_Test = data.frame(y_resid,
                          x_lagendog_resid,
                          z_laginstrument_resid,
                          x_areapast_resid,
                          x_areaplant_resid,
                          x_vaagro_resid,
                          x_pop_resid,
                          x_garimpo_resid,
                          x_conselho_resid,
                          x_lei_resid,
                          x_areaprot_resid,
                          x_multa_resid,
                          x_MPAmazon_resid,
                          CodIBGE)
  
  
  # AreaPlant, AreasPast, POP, Garimpo
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  # AreaPlant, AreasPast, POP, Garimpo, Multa, AreaProt, MPAmazon
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                             
                                             "x_multa_resid",
                                             "x_areaprot_resid",
                                             "x_lei_resid",
                                             "x_MPAmazon_resid"
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                         
                                         "x_multa_resid",
                                         "x_areaprot_resid",
                                         "x_lei_resid",
                                         "x_MPAmazon_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
       }

etable(feols(DesmatINPEn~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta,lag_MPAmazon)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_Conselho, cluster=~CodIBGE, db_amazon), stage = 2, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPri="DesmatINPE",
                lAreaPast = "AreaPast",                                         
                "lag_ConselhoLimpo"="Lag_COMUMArestrito",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()


### IHS no Y ###

etable(feols(DesmatINPEihs~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta,lag_MPAmazon)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_Conselho, cluster=~CodIBGE, db_amazon), stage = 1, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPri="DesmatINPE",
                lAreaPast = "AreaPast",                                         
                "lag_Conselho"="Lag_COMUMA",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()


{
  y = feols(DesmatINPEihs~1|CodIBGE+UF^Year, db_amazon)
  y_resid= y$residuals
  
  x_lagendog = feols(lag_lPeCreal~1|CodIBGE+UF^Year, db_amazon)
  x_lagendog_resid= x_lagendog$residuals
  
  z_laginstrument = feols(lag_Conselho~1|CodIBGE+UF^Year, db_amazon)
  z_laginstrument_resid = z_laginstrument$residuals
  
  x_areapast = feols(lAreaPast~1|CodIBGE+UF^Year, db_amazon)
  x_areapast_resid =  x_areapast$residuals
  
  x_garimpo = feols(lGarimpo~1|CodIBGE+UF^Year, db_amazon)
  x_garimpo_resid =  x_garimpo$residuals
  
  x_vaagro = feols(lVAagro~1|CodIBGE+UF^Year, db_amazon)
  x_vaagro_resid =  x_vaagro$residuals
  
  x_areaplant = feols(lAreaPlant~1|CodIBGE+UF^Year, db_amazon)
  x_areaplant_resid =  x_areaplant$residuals
  
  x_pop = feols(lPOP~1|CodIBGE+UF^Year, db_amazon)
  x_pop_resid =  x_pop$residuals
  
  x_conselho = feols(lag_Conselho~1|CodIBGE+UF^Year, db_amazon)
  x_conselho_resid =  x_conselho$residuals
  
  x_lei = feols(lag_LeiFloresta~1|CodIBGE+UF^Year, db_amazon)
  x_lei_resid =  x_lei$residuals
  
  x_areaprot = feols(lag_D_AreaProt~1|CodIBGE+UF^Year, db_amazon)
  x_areaprot_resid =  x_areaprot$residuals
  
  x_multa = feols(lag_lFines~1|CodIBGE+UF^Year, db_amazon)
  x_multa_resid =  x_multa$residuals
  
  x_MPAmazon = feols(lag_MPAmazon~1|CodIBGE+UF^Year, db_amazon)
  x_MPAmazon_resid =  x_MPAmazon$residuals
  
  CodIBGE = db_amazon$CodIBGE
  
  ## data frame para inputar no teste
  db_AR_Test = data.frame(y_resid,
                          x_lagendog_resid,
                          z_laginstrument_resid,
                          x_areapast_resid,
                          x_areaplant_resid,
                          x_vaagro_resid,
                          x_pop_resid,
                          x_garimpo_resid,
                          x_conselho_resid,
                          x_lei_resid,
                          x_areaprot_resid,
                          x_multa_resid,
                          x_MPAmazon_resid,
                          CodIBGE)
  
  
  # AreaPlant, AreasPast, POP, Garimpo
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  # AreaPlant, AreasPast, POP, Garimpo, Multa, AreaProt, MPAmazon
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                             
                                             "x_multa_resid",
                                             "x_areaprot_resid",
                                             "x_lei_resid",
                                             "x_MPAmazon_resid"
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                         
                                         "x_multa_resid",
                                         "x_areaprot_resid",
                                         "x_lei_resid",
                                         "x_MPAmazon_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
       }

etable(feols(DesmatINPEihs~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta,lag_MPAmazon)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_Conselho, cluster=~CodIBGE, db_amazon), stage = 2, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPri="DesmatINPE",
                lAreaPast = "AreaPast",                                         
                "lag_ConselhoLimpo"="Lag_COMUMArestrito",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()



## OLS 

etable(feols(lDesmatINPE~csw(lag_lPeCreal,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta,lag_MPAmazon)|
               CodIBGE+UF^Year, cluster=~CodIBGE, db_amazon)) %>% View()


## Forma reduzida  

etable(feols(lDesmatINPE~csw(lag_Conselho,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta,lag_MPAmazon)|
               CodIBGE+UF^Year, cluster=~CodIBGE, db_amazon)) %>% View()



# Exercício Extra: Gestão Ambiental ----------------------------------------------------


etable(feols(lDesmatINPE~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta,lag_MPAmazon)|
               CodIBGE+UF^Year|lag_lGestAreal~lag_Conselho, cluster=~CodIBGE, db_amazon), stage = 1, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPriINPE="DesmatINPE",
                lAreaPast = "AreaPast",                                         
                "lag_Conselho"="Lag_COMUMA",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()


{
  y = feols(lDesmatINPE~1|CodIBGE+UF^Year, db_amazon)
  y_resid= y$residuals
  
  x_lagendog = feols(lag_lGestAreal~1|CodIBGE+UF^Year, db_amazon)
  x_lagendog_resid= x_lagendog$residuals
  
  z_laginstrument = feols(lag_Conselho~1|CodIBGE+UF^Year, db_amazon)
  z_laginstrument_resid = z_laginstrument$residuals
  
  x_areapast = feols(lAreaPast~1|CodIBGE+UF^Year, db_amazon)
  x_areapast_resid =  x_areapast$residuals
  
  x_garimpo = feols(lGarimpo~1|CodIBGE+UF^Year, db_amazon)
  x_garimpo_resid =  x_garimpo$residuals
  
  x_vaagro = feols(lVAagro~1|CodIBGE+UF^Year, db_amazon)
  x_vaagro_resid =  x_vaagro$residuals
  
  x_areaplant = feols(lAreaPlant~1|CodIBGE+UF^Year, db_amazon)
  x_areaplant_resid =  x_areaplant$residuals
  
  x_pop = feols(lPOP~1|CodIBGE+UF^Year, db_amazon)
  x_pop_resid =  x_pop$residuals
  
  x_conselho = feols(lag_Conselho~1|CodIBGE+UF^Year, db_amazon)
  x_conselho_resid =  x_conselho$residuals
  
  x_lei = feols(lag_LeiFloresta~1|CodIBGE+UF^Year, db_amazon)
  x_lei_resid =  x_lei$residuals
  
  x_areaprot = feols(lag_D_AreaProt~1|CodIBGE+UF^Year, db_amazon)
  x_areaprot_resid =  x_areaprot$residuals
  
  x_multa = feols(lag_lFines~1|CodIBGE+UF^Year, db_amazon)
  x_multa_resid =  x_multa$residuals
  
  x_MPAmazon = feols(lag_MPAmazon~1|CodIBGE+UF^Year, db_amazon)
  x_MPAmazon_resid =  x_MPAmazon$residuals
  
  CodIBGE = db_amazon$CodIBGE
  
  ## data frame para inputar no teste
  db_AR_Test = data.frame(y_resid,
                          x_lagendog_resid,
                          z_laginstrument_resid,
                          x_areapast_resid,
                          x_areaplant_resid,
                          x_vaagro_resid,
                          x_pop_resid,
                          x_garimpo_resid,
                          x_conselho_resid,
                          x_lei_resid,
                          x_areaprot_resid,
                          x_multa_resid,
                          x_MPAmazon_resid,
                          CodIBGE)
  
  
  # AreaPlant, AreasPast, POP, Garimpo
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  # AreaPlant, AreasPast, POP, Garimpo, Multa, AreaProt, MPAmazon
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                             
                                             "x_multa_resid",
                                             "x_areaprot_resid",
                                             "x_lei_resid",
                                             "x_MPAmazon_resid"
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                         
                                         "x_multa_resid",
                                         "x_areaprot_resid",
                                         "x_lei_resid",
                                         "x_MPAmazon_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
       }

etable(feols(lDesmatINPE~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta,lag_MPAmazon)|
               CodIBGE+UF^Year|lag_lGestAreal~lag_Conselho, cluster=~CodIBGE, db_amazon), stage = 2, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPri="DesmatINPE",
                lAreaPast = "AreaPast",                                         
                "lag_ConselhoLimpo"="Lag_COMUMA",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()



### 4.0 CERRADO #####

desmatCerrado = read_excel("DesmatCerrado.xlsx")

cerrado_long = reshape(desmatCerrado, idvar = "CodIBGE",
                        varying = 3:36,
                        sep="",
                        timevar = "Year",
                        times = c(2004:2020),
                        new.row.names = 1:36762,
                        direction = "long")

cerrado_long= cerrado_long %>% mutate_at(c("DesmatPriCE"), as.numeric) %>% subset(Year<2020)

n_distinct(cerrado_long$CodIBGE)

db_cerrado_long = right_join(db_mun_long2, cerrado_long, by = c("CodIBGE","Year"), keep = FALSE)

n_distinct(db_cerrado_long$CodIBGE)

## Incluindo Logs Adicionais para o desmatamento exclusivo do Bioma Cerrado (sufixo CE)

db_cerrado_long = db_cerrado_long %>% mutate(
                                           lDesmatPriCE=log(DesmatPriCE+1))

## Convertendo em Painel

{

db_cerrado = select(db_cerrado_long,"CodIBGE","UF","Year","Nota","RemanescenteCE","DesmatPriCE","lDesmatPriCE", "DesmatPriCEihs",
                    "lag_lPeCreal","lag_PeCrealIHS","lag_lGestAreal","lag_GestArealIHS", 
                    "lag_lReceitaTrib","lag_lReceitaTotal","lag_D_AreaProt","lag_lFines",
                    "lag_ConselhoLimpo", "lag_Conselho","Conselho","lVAagro",
                    "lAreaPlant","lAreaPast","lPOP", "lPIBtot","lGarimpo","lag_MPCerrado","lag_LeiFloresta") %>% na.omit()

  n_distinct(db_cerrado$CodIBGE)  
  
db_cerrado = db_cerrado %>% group_by(CodIBGE) %>% 
  mutate(DesmatMean = mean(DesmatPriCE), 
         DesmatStd = sd(DesmatPriCE),
         DesmatPriCEn = case_when(DesmatMean > 0 ~ (DesmatPriCE - DesmatMean)/DesmatStd)) %>%  na.omit()


n_distinct(db_cerrado$CodIBGE)

db_cerrado = panel(db_cerrado, ~CodIBGE+Year) 

}


## Estatisticas descritivas

db_cerrado_descr = db_cerrado %>% as.data.frame() %>% mutate(PeC = exp(lag_lPeCreal)-1, 
                                                           DesmatPri = exp(lDesmatPriCE)-1, 
                                                         
                                                           AreaPlant =  exp(lAreaPlant)-1,
                                                           AreaPast = exp(lAreaPast)-1,
                                                           POP = exp(lPOP),
                                                           Garimpo = exp(lGarimpo)-1,
                                                           Multas = exp(lag_lFines)-1) 


  select("DesmatPriCE","PeC","AreaPlant","POP","Garimpo","Multas", 
         "lDesmatPriCE","lag_lPeCreal","lAreaPast","lPOP","lGarimpo","lag_lFines") %>% 
  descr(order = "p",round.digits = 2)


db_cerrado_descr %>% select("lDesmatPriCE","AreaPast","lAreaPlant") %>% descr(round.digits=2)



# 4.1 IV - CONSELHO IRRESTRITO ----

## Primeiro Estágio 
etable(feols(lDesmatPriCE~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta,lag_MPCerrado)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_Conselho, cluster=~CodIBGE, db_cerrado), stage = 1, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPriINPE="DesmatPriCE",
                lAreaPast = "AreaPast",                                         
                "lag_Conselho"="Lag_COMUMA",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()



{
  y = feols(lDesmatPriCE~1|CodIBGE+UF^Year, db_cerrado)
  y_resid= y$residuals
  
  x_lagendog = feols(lag_lPeCreal~1|CodIBGE+UF^Year, db_cerrado)
  x_lagendog_resid= x_lagendog$residuals
  
  z_laginstrument = feols(lag_Conselho~1|CodIBGE+UF^Year, db_cerrado)
  z_laginstrument_resid = z_laginstrument$residuals
  
  x_areapast = feols(lAreaPast~1|CodIBGE+UF^Year, db_cerrado)
  x_areapast_resid =  x_areapast$residuals
  
  x_garimpo = feols(lGarimpo~1|CodIBGE+UF^Year, db_cerrado)
  x_garimpo_resid =  x_garimpo$residuals
  
  x_vaagro = feols(lVAagro~1|CodIBGE+UF^Year, db_cerrado)
  x_vaagro_resid =  x_vaagro$residuals
  
  x_areaplant = feols(lAreaPlant~1|CodIBGE+UF^Year, db_cerrado)
  x_areaplant_resid =  x_areaplant$residuals
  
  x_pop = feols(lPOP~1|CodIBGE+UF^Year, db_cerrado)
  x_pop_resid =  x_pop$residuals
  
  x_conselho = feols(lag_Conselho~1|CodIBGE+UF^Year, db_cerrado)
  x_conselho_resid =  x_conselho$residuals
  
  x_lei = feols(lag_LeiFloresta~1|CodIBGE+UF^Year, db_cerrado)
  x_lei_resid =  x_lei$residuals
  
  x_areaprot = feols(lag_D_AreaProt~1|CodIBGE+UF^Year, db_cerrado)
  x_areaprot_resid =  x_areaprot$residuals
  
  x_multa = feols(lag_lFines~1|CodIBGE+UF^Year, db_cerrado)
  x_multa_resid =  x_multa$residuals
  
  x_MPCerrado = feols(lag_MPCerrado~1|CodIBGE+UF^Year, db_cerrado)
  x_MPCerrado_resid =  x_MPCerrado$residuals
  
  CodIBGE = db_cerrado$CodIBGE
  
  ## data frame para inputar no teste
  db_AR_Test = data.frame(y_resid,
                          x_lagendog_resid,
                          z_laginstrument_resid,
                          x_areapast_resid,
                          x_areaplant_resid,
                          x_vaagro_resid,
                          x_pop_resid,
                          x_garimpo_resid,
                          x_conselho_resid,
                          x_lei_resid,
                          x_areaprot_resid,
                          x_multa_resid,
                          x_MPCerrado_resid,
                          CodIBGE)
  
  # Nada
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                              vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                          
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  
  # AreaPlant
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  # AreaPlant, AreaPast 
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  
  # AreaPlant, AreasPast, POP
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  
  
  # # AreaPlant, AreasPast, POP, Garimpo
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  #AreaPlant, AreasPast, POP, Garimpo, Multa 
  
  
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                             
                                             "x_multa_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                         
                                         "x_multa_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  # AreaPlant, AreasPast, POP, Garimpo, Multa, AreaProt
  
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                             
                                             "x_multa_resid",
                                             "x_areaprot_resid"
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                         
                                         "x_multa_resid",
                                         "x_areaprot_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  
  # AreaPlant, AreasPast, POP, Garimpo, Multa, AreaProt, Lei
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                             
                                             "x_multa_resid",
                                             "x_areaprot_resid",
                                             "x_lei_resid"
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                         
                                         "x_multa_resid",
                                         "x_areaprot_resid",
                                         "x_lei_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  
  # AreaPlant, AreasPast, POP, Garimpo, Multa, AreaProt, MPCerrado
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                             
                                             "x_multa_resid",
                                             "x_areaprot_resid",
                                             "x_lei_resid",
                                             "x_MPCerrado_resid"
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                         
                                         "x_multa_resid",
                                         "x_areaprot_resid",
                                         "x_lei_resid",
                                         "x_MPCerrado_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
       }


etable(feols(lDesmatPriCE~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta,lag_MPCerrado)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_Conselho, cluster=~CodIBGE, db_cerrado), stage = 2, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPriINPE="DesmatPriCE",
                lAreaPast = "AreaPast",                                         
                "lag_Conselho"="Lag_COMUMA",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()





# Robustez------------------------------------------------------------ 


## Conselho Restrito ## 

## Primeiro Estágio 
etable(feols(lDesmatPriCE~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta,lag_MPCerrado)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_ConselhoLimpo, cluster=~CodIBGE, db_cerrado), stage = 1, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPriINPE="DesmatPriCE",
                lAreaPast = "AreaPast",                                         
                "lag_ConselhoLimpo"="Lag_COMUMArestrito",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()


{
  y = feols(lDesmatPriCE~1|CodIBGE+UF^Year, db_cerrado)
  y_resid= y$residuals
  
  x_lagendog = feols(lag_lPeCreal~1|CodIBGE+UF^Year, db_cerrado)
  x_lagendog_resid= x_lagendog$residuals
  
  z_laginstrument = feols(lag_ConselhoLimpo~1|CodIBGE+UF^Year, db_cerrado)
  z_laginstrument_resid = z_laginstrument$residuals
  
  x_areapast = feols(lAreaPast~1|CodIBGE+UF^Year, db_cerrado)
  x_areapast_resid =  x_areapast$residuals
  
  x_garimpo = feols(lGarimpo~1|CodIBGE+UF^Year, db_cerrado)
  x_garimpo_resid =  x_garimpo$residuals
  
  x_vaagro = feols(lVAagro~1|CodIBGE+UF^Year, db_cerrado)
  x_vaagro_resid =  x_vaagro$residuals
  
  x_areaplant = feols(lAreaPlant~1|CodIBGE+UF^Year, db_cerrado)
  x_areaplant_resid =  x_areaplant$residuals
  
  x_pop = feols(lPOP~1|CodIBGE+UF^Year, db_cerrado)
  x_pop_resid =  x_pop$residuals
  
  x_conselho = feols(lag_ConselhoLimpo~1|CodIBGE+UF^Year, db_cerrado)
  x_conselho_resid =  x_conselho$residuals
  
  x_lei = feols(lag_LeiFloresta~1|CodIBGE+UF^Year, db_cerrado)
  x_lei_resid =  x_lei$residuals
  
  x_areaprot = feols(lag_D_AreaProt~1|CodIBGE+UF^Year, db_cerrado)
  x_areaprot_resid =  x_areaprot$residuals
  
  x_multa = feols(lag_lFines~1|CodIBGE+UF^Year, db_cerrado)
  x_multa_resid =  x_multa$residuals
  
  x_MPCerrado = feols(lag_MPCerrado~1|CodIBGE+UF^Year, db_cerrado)
  x_MPCerrado_resid =  x_MPCerrado$residuals
  
  CodIBGE = db_cerrado$CodIBGE
  
  ## data frame para inputar no teste
  db_AR_Test = data.frame(y_resid,
                          x_lagendog_resid,
                          z_laginstrument_resid,
                          x_areapast_resid,
                          x_areaplant_resid,
                          x_vaagro_resid,
                          x_pop_resid,
                          x_garimpo_resid,
                          x_conselho_resid,
                          x_lei_resid,
                          x_areaprot_resid,
                          x_multa_resid,
                          x_MPCerrado_resid,
                          CodIBGE)
  
  # Nada
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid"
                                
                                , vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  
  
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid"
                           
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  
  
  # AreaPlant
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  # AreaPlant, AreaPast 
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  
  # AreaPlant, AreasPast, POP
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  
  
  # # AreaPlant, AreasPast, POP, Garimpo
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  #AreaPlant, AreasPast, POP, Garimpo, Multa 
  
  
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                             
                                             "x_multa_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                         
                                         "x_multa_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  # AreaPlant, AreasPast, POP, Garimpo, Multa, AreaProt
  
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                             
                                             "x_multa_resid",
                                             "x_areaprot_resid"
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                         
                                         "x_multa_resid",
                                         "x_areaprot_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  
  # AreaPlant, AreasPast, POP, Garimpo, Multa, AreaProt, Lei
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                             
                                             "x_multa_resid",
                                             "x_areaprot_resid",
                                             "x_lei_resid"
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                         
                                         "x_multa_resid",
                                         "x_areaprot_resid",
                                         "x_lei_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  
  # AreaPlant, AreasPast, POP, Garimpo, Multa, AreaProt, MPCerrado
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                             
                                             "x_multa_resid",
                                             "x_areaprot_resid",
                                             "x_lei_resid",
                                             "x_MPCerrado_resid"
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                         
                                         "x_multa_resid",
                                         "x_areaprot_resid",
                                         "x_lei_resid",
                                         "x_MPCerrado_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
       }

etable(feols(lDesmatPriCE~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta,lag_MPCerrado)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_ConselhoLimpo, cluster=~CodIBGE, db_cerrado), stage = 2, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPriINPE="DesmatPriCE",
                lAreaPast = "AreaPast",                                         
                "lag_ConselhoLimpo"="Lag_COMUMArestrito",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()



### Normalizado ###

etable(feols(DesmatPriCEn~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta,lag_MPCerrado)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_Conselho, cluster=~CodIBGE, db_cerrado), stage = 1, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPri="DesmatPriCE",
                lAreaPast = "AreaPast",                                         
                "lag_Conselho"="Lag_COMUMA",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()


{
  y = feols(DesmatPriCEn~1|CodIBGE+UF^Year, db_cerrado)
  y_resid= y$residuals
  
  x_lagendog = feols(lag_lPeCreal~1|CodIBGE+UF^Year, db_cerrado)
  x_lagendog_resid= x_lagendog$residuals
  
  z_laginstrument = feols(lag_Conselho~1|CodIBGE+UF^Year, db_cerrado)
  z_laginstrument_resid = z_laginstrument$residuals
  
  x_areapast = feols(lAreaPast~1|CodIBGE+UF^Year, db_cerrado)
  x_areapast_resid =  x_areapast$residuals
  
  x_garimpo = feols(lGarimpo~1|CodIBGE+UF^Year, db_cerrado)
  x_garimpo_resid =  x_garimpo$residuals
  
  x_vaagro = feols(lVAagro~1|CodIBGE+UF^Year, db_cerrado)
  x_vaagro_resid =  x_vaagro$residuals
  
  x_areaplant = feols(lAreaPlant~1|CodIBGE+UF^Year, db_cerrado)
  x_areaplant_resid =  x_areaplant$residuals
  
  x_pop = feols(lPOP~1|CodIBGE+UF^Year, db_cerrado)
  x_pop_resid =  x_pop$residuals
  
  x_conselho = feols(lag_Conselho~1|CodIBGE+UF^Year, db_cerrado)
  x_conselho_resid =  x_conselho$residuals
  
  x_lei = feols(lag_LeiFloresta~1|CodIBGE+UF^Year, db_cerrado)
  x_lei_resid =  x_lei$residuals
  
  x_areaprot = feols(lag_D_AreaProt~1|CodIBGE+UF^Year, db_cerrado)
  x_areaprot_resid =  x_areaprot$residuals
  
  x_multa = feols(lag_lFines~1|CodIBGE+UF^Year, db_cerrado)
  x_multa_resid =  x_multa$residuals
  
  x_MPCerrado = feols(lag_MPCerrado~1|CodIBGE+UF^Year, db_cerrado)
  x_MPCerrado_resid =  x_MPCerrado$residuals
  
  CodIBGE = db_cerrado$CodIBGE
  
  ## data frame para inputar no teste
  db_AR_Test = data.frame(y_resid,
                          x_lagendog_resid,
                          z_laginstrument_resid,
                          x_areapast_resid,
                          x_areaplant_resid,
                          x_vaagro_resid,
                          x_pop_resid,
                          x_garimpo_resid,
                          x_conselho_resid,
                          x_lei_resid,
                          x_areaprot_resid,
                          x_multa_resid,
                          x_MPCerrado_resid,
                          CodIBGE)
  
  
  # AreaPlant, AreasPast, POP, Garimpo
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  # AreaPlant, AreasPast, POP, Garimpo, Multa, AreaProt, MPCerrado
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                             
                                             "x_multa_resid",
                                             "x_areaprot_resid",
                                             "x_lei_resid",
                                             "x_MPCerrado_resid"
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                         
                                         "x_multa_resid",
                                         "x_areaprot_resid",
                                         "x_lei_resid",
                                         "x_MPCerrado_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
       }

etable(feols(DesmatPriCEn~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta,lag_MPCerrado)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_Conselho, cluster=~CodIBGE, db_cerrado), stage = 2, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPri="DesmatPriCE",
                lAreaPast = "AreaPast",                                         
                "lag_ConselhoLimpo"="Lag_COMUMArestrito",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()


### IHS no Y ###

etable(feols(DesmatPriCEihs~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta,lag_MPCerrado)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_Conselho, cluster=~CodIBGE, db_cerrado), stage = 1, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPri="DesmatPriCE",
                lAreaPast = "AreaPast",                                         
                "lag_Conselho"="Lag_COMUMA",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()


{
  y = feols(DesmatPriCEihs~1|CodIBGE+UF^Year, db_cerrado)
  y_resid= y$residuals
  
  x_lagendog = feols(lag_lPeCreal~1|CodIBGE+UF^Year, db_cerrado)
  x_lagendog_resid= x_lagendog$residuals
  
  z_laginstrument = feols(lag_Conselho~1|CodIBGE+UF^Year, db_cerrado)
  z_laginstrument_resid = z_laginstrument$residuals
  
  x_areapast = feols(lAreaPast~1|CodIBGE+UF^Year, db_cerrado)
  x_areapast_resid =  x_areapast$residuals
  
  x_garimpo = feols(lGarimpo~1|CodIBGE+UF^Year, db_cerrado)
  x_garimpo_resid =  x_garimpo$residuals
  
  x_vaagro = feols(lVAagro~1|CodIBGE+UF^Year, db_cerrado)
  x_vaagro_resid =  x_vaagro$residuals
  
  x_areaplant = feols(lAreaPlant~1|CodIBGE+UF^Year, db_cerrado)
  x_areaplant_resid =  x_areaplant$residuals
  
  x_pop = feols(lPOP~1|CodIBGE+UF^Year, db_cerrado)
  x_pop_resid =  x_pop$residuals
  
  x_conselho = feols(lag_Conselho~1|CodIBGE+UF^Year, db_cerrado)
  x_conselho_resid =  x_conselho$residuals
  
  x_lei = feols(lag_LeiFloresta~1|CodIBGE+UF^Year, db_cerrado)
  x_lei_resid =  x_lei$residuals
  
  x_areaprot = feols(lag_D_AreaProt~1|CodIBGE+UF^Year, db_cerrado)
  x_areaprot_resid =  x_areaprot$residuals
  
  x_multa = feols(lag_lFines~1|CodIBGE+UF^Year, db_cerrado)
  x_multa_resid =  x_multa$residuals
  
  x_MPCerrado = feols(lag_MPCerrado~1|CodIBGE+UF^Year, db_cerrado)
  x_MPCerrado_resid =  x_MPCerrado$residuals
  
  CodIBGE = db_cerrado$CodIBGE
  
  ## data frame para inputar no teste
  db_AR_Test = data.frame(y_resid,
                          x_lagendog_resid,
                          z_laginstrument_resid,
                          x_areapast_resid,
                          x_areaplant_resid,
                          x_vaagro_resid,
                          x_pop_resid,
                          x_garimpo_resid,
                          x_conselho_resid,
                          x_lei_resid,
                          x_areaprot_resid,
                          x_multa_resid,
                          x_MPCerrado_resid,
                          CodIBGE)
  
  
  # AreaPlant, AreasPast, POP, Garimpo
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  # AreaPlant, AreasPast, POP, Garimpo, Multa, AreaProt, MPCerrado
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                             
                                             "x_multa_resid",
                                             "x_areaprot_resid",
                                             "x_lei_resid",
                                             "x_MPCerrado_resid"
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                         
                                         "x_multa_resid",
                                         "x_areaprot_resid",
                                         "x_lei_resid",
                                         "x_MPCerrado_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
       }

etable(feols(DesmatPriCEihs~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta,lag_MPCerrado)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_Conselho, cluster=~CodIBGE, db_cerrado), stage = 2, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPri="DesmatPriCE",
                lAreaPast = "AreaPast",                                         
                "lag_ConselhoLimpo"="Lag_COMUMArestrito",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()


## OLS 

etable(feols(lDesmatPriCE~csw(lag_lPeCreal,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta,lag_MPCerrado)|
               CodIBGE+UF^Year, cluster=~CodIBGE, db_cerrado)) %>% View()


## Forma Reduzida

etable(feols(lDesmatPriCE~csw(lag_Conselho,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta,lag_MPCerrado)|
               CodIBGE+UF^Year, cluster=~CodIBGE, db_cerrado)) %>% View()






# Exercício Extra: Gestão Ambiental ----------------------------------------------------


etable(feols(lDesmatPriCE~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta,lag_MPCerrado)|
               CodIBGE+UF^Year|lag_lGestAreal~lag_Conselho, cluster=~CodIBGE, db_cerrado), stage = 1, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPriINPE="DesmatPriCE",
                lAreaPast = "AreaPast",                                         
                "lag_Conselho"="Lag_COMUMA",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()


{
  y = feols(lDesmatPriCE~1|CodIBGE+UF^Year, db_cerrado)
  y_resid= y$residuals
  
  x_lagendog = feols(lag_lGestAreal~1|CodIBGE+UF^Year, db_cerrado)
  x_lagendog_resid= x_lagendog$residuals
  
  z_laginstrument = feols(lag_Conselho~1|CodIBGE+UF^Year, db_cerrado)
  z_laginstrument_resid = z_laginstrument$residuals
  
  x_areapast = feols(lAreaPast~1|CodIBGE+UF^Year, db_cerrado)
  x_areapast_resid =  x_areapast$residuals
  
  x_garimpo = feols(lGarimpo~1|CodIBGE+UF^Year, db_cerrado)
  x_garimpo_resid =  x_garimpo$residuals
  
  x_vaagro = feols(lVAagro~1|CodIBGE+UF^Year, db_cerrado)
  x_vaagro_resid =  x_vaagro$residuals
  
  x_areaplant = feols(lAreaPlant~1|CodIBGE+UF^Year, db_cerrado)
  x_areaplant_resid =  x_areaplant$residuals
  
  x_pop = feols(lPOP~1|CodIBGE+UF^Year, db_cerrado)
  x_pop_resid =  x_pop$residuals
  
  x_conselho = feols(lag_Conselho~1|CodIBGE+UF^Year, db_cerrado)
  x_conselho_resid =  x_conselho$residuals
  
  x_lei = feols(lag_LeiFloresta~1|CodIBGE+UF^Year, db_cerrado)
  x_lei_resid =  x_lei$residuals
  
  x_areaprot = feols(lag_D_AreaProt~1|CodIBGE+UF^Year, db_cerrado)
  x_areaprot_resid =  x_areaprot$residuals
  
  x_multa = feols(lag_lFines~1|CodIBGE+UF^Year, db_cerrado)
  x_multa_resid =  x_multa$residuals
  
  x_MPCerrado = feols(lag_MPCerrado~1|CodIBGE+UF^Year, db_cerrado)
  x_MPCerrado_resid =  x_MPCerrado$residuals
  
  CodIBGE = db_cerrado$CodIBGE
  
  ## data frame para inputar no teste
  db_AR_Test = data.frame(y_resid,
                          x_lagendog_resid,
                          z_laginstrument_resid,
                          x_areapast_resid,
                          x_areaplant_resid,
                          x_vaagro_resid,
                          x_pop_resid,
                          x_garimpo_resid,
                          x_conselho_resid,
                          x_lei_resid,
                          x_areaprot_resid,
                          x_multa_resid,
                          x_MPCerrado_resid,
                          CodIBGE)
  
  
  # AreaPlant, AreasPast, POP, Garimpo
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  # AreaPlant, AreasPast, POP, Garimpo, Multa, AreaProt, MPCerrado
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                             
                                             "x_multa_resid",
                                             "x_areaprot_resid",
                                             "x_lei_resid",
                                             "x_MPCerrado_resid"
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                         
                                         "x_multa_resid",
                                         "x_areaprot_resid",
                                         "x_lei_resid",
                                         "x_MPCerrado_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
       }

etable(feols(lDesmatPriCE~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta,lag_MPCerrado)|
               CodIBGE+UF^Year|lag_lGestAreal~lag_Conselho, cluster=~CodIBGE, db_cerrado), stage = 2, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPri="DesmatPriCE",
                lAreaPast = "AreaPast",                                         
                "lag_ConselhoLimpo"="Lag_COMUMA",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()


### 5.0 MATA ATLÂNTICA ---------------------------------------------------------

## preaprando para filtrar apenas os municipios da Mata Atlântica
desmatMata = read_excel("DesmatMata.xlsx")

mata_long = reshape(desmatMata, idvar = "CodIBGE",
                     varying = 3:36,
                     sep="",
                     timevar = "Year",
                     times = c(2004:2020),
                     new.row.names = 1:52275,
                     direction = "long")

mata_long= mata_long %>% mutate_at(c("DesmatPriMA"), as.numeric)

db_mata_long = right_join(db_mun_long2, mata_long, by = c("CodIBGE","Year"), keep = FALSE)

db_mata_long = db_mata_long %>% mutate(lDesmatPriMA=log(DesmatPriMA+1))


## Convertendo em Painel
{
   
  db_mata = select(db_mata_long,"CodIBGE","UF","Year","Nota", "RemanescenteMA","DesmatPriMA","lDesmatPriMA","DesmatPriMAihs",
                      "lag_lPeCreal","lag_PeCrealIHS","lag_lGestAreal","lag_GestArealIHS", "lag_lFines","lag_D_AreaProt","lag_lReceitaTotal",
                      "lag_lReceitaTrib","lag_ConselhoLimpo", "lag_Conselho","Conselho",
                      "lAreaPlant","lAreaPast","lPOP","lPIBtot","lGarimpo","lag_LeiFloresta","lVAagro") %>% na.omit()
  
  n_distinct(db_mata$CodIBGE)

  db_mata = db_mata %>% group_by(CodIBGE) %>% 
    mutate(DesmatMean = mean(DesmatPriMA), 
           DesmatStd = sd(DesmatPriMA),
           DesmatPriMAn = case_when(DesmatMean > 0 ~ (DesmatPriMA - DesmatMean)/DesmatStd)) %>%
    na.omit()
  
  n_distinct(db_mata$CodIBGE)
  
    db_mata = panel(db_mata, ~CodIBGE+Year)  
}



## Estatisticas descritivas

db_mata_descr = db_mata %>% as.data.frame() %>% mutate(PeC = exp(lag_lPeCreal)-1, 
                                                             DesmatPri = exp(lDesmatPriMA)-1, 
                                                             
                                                             AreaPlant =  exp(lAreaPlant)-1,
                                                             AreaPast = exp(lAreaPast)-1,
                                                             POP = exp(lPOP),
                                                             Garimpo = exp(lGarimpo)-1,
                                                             Multas = exp(lag_lFines)-1)
 
  select("DesmatPriMA","PeC","AreaPlant","POP","Garimpo","Multas", 
         "lDesmatPriMA","lag_lPeCreal","lAreaPast","lPOP","lGarimpo","lag_lFines") %>% 
  descr(order = "p",round.digits = 2)

  
  db_mata_descr %>% select("lDesmatPriMA","AreaPast","lAreaPlant") %>% descr(round.digits=2)


# 4.1 IV - CONSELHO IRRESTRITO ----

## Primeiro Estágio 
etable(feols(lDesmatPriMA~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_Conselho, cluster=~CodIBGE, db_mata), stage = 1, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPriINPE="DesmatPriMA",
                lAreaPast = "AreaPast",                                         
                "lag_Conselho"="Lag_COMUMA",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()



{
  y = feols(lDesmatPriMA~1|CodIBGE+UF^Year, db_mata)
  y_resid= y$residuals
  
  x_lagendog = feols(lag_lPeCreal~1|CodIBGE+UF^Year, db_mata)
  x_lagendog_resid= x_lagendog$residuals
  
  z_laginstrument = feols(lag_Conselho~1|CodIBGE+UF^Year, db_mata)
  z_laginstrument_resid = z_laginstrument$residuals
  
  x_areapast = feols(lAreaPast~1|CodIBGE+UF^Year, db_mata)
  x_areapast_resid =  x_areapast$residuals
  
  x_garimpo = feols(lGarimpo~1|CodIBGE+UF^Year, db_mata)
  x_garimpo_resid =  x_garimpo$residuals
  
  x_vaagro = feols(lVAagro~1|CodIBGE+UF^Year, db_mata)
  x_vaagro_resid =  x_vaagro$residuals
  
  x_areaplant = feols(lAreaPlant~1|CodIBGE+UF^Year, db_mata)
  x_areaplant_resid =  x_areaplant$residuals
  
  x_pop = feols(lPOP~1|CodIBGE+UF^Year, db_mata)
  x_pop_resid =  x_pop$residuals
  
  x_conselho = feols(lag_Conselho~1|CodIBGE+UF^Year, db_mata)
  x_conselho_resid =  x_conselho$residuals
  
  x_lei = feols(lag_LeiFloresta~1|CodIBGE+UF^Year, db_mata)
  x_lei_resid =  x_lei$residuals
  
  x_areaprot = feols(lag_D_AreaProt~1|CodIBGE+UF^Year, db_mata)
  x_areaprot_resid =  x_areaprot$residuals
  
  x_multa = feols(lag_lFines~1|CodIBGE+UF^Year, db_mata)
  x_multa_resid =  x_multa$residuals

  CodIBGE = db_mata$CodIBGE
  
  ## data frame para inputar no teste
  db_AR_Test = data.frame(y_resid,
                          x_lagendog_resid,
                          z_laginstrument_resid,
                          x_areapast_resid,
                          x_areaplant_resid,
                          x_vaagro_resid,
                          x_pop_resid,
                          x_garimpo_resid,
                          x_conselho_resid,
                          x_lei_resid,
                          x_areaprot_resid,
                          x_multa_resid,
                          x_MPCerrado_resid,
                          CodIBGE)
  
  # AreaPlant
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  # AreaPlant, AreaPast 
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  
  # AreaPlant, AreasPast, POP
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  
  
  # # AreaPlant, AreasPast, POP, Garimpo
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  #AreaPlant, AreasPast, POP, Garimpo, Multa 
  
  
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                             
                                             "x_multa_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                         
                                         "x_multa_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  # AreaPlant, AreasPast, POP, Garimpo, Multa, AreaProt
  
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                             
                                             "x_multa_resid",
                                             "x_areaprot_resid"
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                         
                                         "x_multa_resid",
                                         "x_areaprot_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  
  # AreaPlant, AreasPast, POP, Garimpo, Multa, AreaProt, Lei
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                             
                                             "x_multa_resid",
                                             "x_areaprot_resid",
                                             "x_lei_resid"
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                         
                                         "x_multa_resid",
                                         "x_areaprot_resid",
                                         "x_lei_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  
  # AreaPlant, AreasPast, POP, Garimpo, Multa, AreaProt, MPCerrado
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                             
                                             "x_multa_resid",
                                             "x_areaprot_resid",
                                             "x_lei_resid",
                                             "x_MPCerrado_resid"
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                         
                                         "x_multa_resid",
                                         "x_areaprot_resid",
                                         "x_lei_resid",
                                         "x_MPCerrado_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
       }


etable(feols(lDesmatPriMA~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_Conselho, cluster=~CodIBGE, db_mata), stage = 2, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPriINPE="DesmatPriMA",
                lAreaPast = "AreaPast",                                         
                "lag_Conselho"="Lag_COMUMA",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()





# Robustez------------------------------------------------------------ 


## Conselho Restrito ## 

## Primeiro Estágio 
etable(feols(lDesmatPriMA~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_ConselhoLimpo, cluster=~CodIBGE, db_mata), stage = 1, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPriINPE="DesmatPriMA",
                lAreaPast = "AreaPast",                                         
                "lag_ConselhoLimpo"="Lag_COMUMArestrito",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()


{
  y = feols(lDesmatPriMA~1|CodIBGE+UF^Year, db_mata)
  y_resid= y$residuals
  
  x_lagendog = feols(lag_lPeCreal~1|CodIBGE+UF^Year, db_mata)
  x_lagendog_resid= x_lagendog$residuals
  
  z_laginstrument = feols(lag_ConselhoLimpo~1|CodIBGE+UF^Year, db_mata)
  z_laginstrument_resid = z_laginstrument$residuals
  
  x_areapast = feols(lAreaPast~1|CodIBGE+UF^Year, db_mata)
  x_areapast_resid =  x_areapast$residuals
  
  x_garimpo = feols(lGarimpo~1|CodIBGE+UF^Year, db_mata)
  x_garimpo_resid =  x_garimpo$residuals
  
  x_vaagro = feols(lVAagro~1|CodIBGE+UF^Year, db_mata)
  x_vaagro_resid =  x_vaagro$residuals
  
  x_areaplant = feols(lAreaPlant~1|CodIBGE+UF^Year, db_mata)
  x_areaplant_resid =  x_areaplant$residuals
  
  x_pop = feols(lPOP~1|CodIBGE+UF^Year, db_mata)
  x_pop_resid =  x_pop$residuals
  
  x_conselho = feols(lag_Conselho~1|CodIBGE+UF^Year, db_mata)
  x_conselho_resid =  x_conselho$residuals
  
  x_lei = feols(lag_LeiFloresta~1|CodIBGE+UF^Year, db_mata)
  x_lei_resid =  x_lei$residuals
  
  x_areaprot = feols(lag_D_AreaProt~1|CodIBGE+UF^Year, db_mata)
  x_areaprot_resid =  x_areaprot$residuals
  
  x_multa = feols(lag_lFines~1|CodIBGE+UF^Year, db_mata)
  x_multa_resid =  x_multa$residuals
  
  x_MPCerrado = feols(lag_MPCerrado~1|CodIBGE+UF^Year, db_mata)
  x_MPCerrado_resid =  x_MPCerrado$residuals
  
  CodIBGE = db_mata$CodIBGE
  
  ## data frame para inputar no teste
  db_AR_Test = data.frame(y_resid,
                          x_lagendog_resid,
                          z_laginstrument_resid,
                          x_areapast_resid,
                          x_areaplant_resid,
                          x_vaagro_resid,
                          x_pop_resid,
                          x_garimpo_resid,
                          x_conselho_resid,
                          x_lei_resid,
                          x_areaprot_resid,
                          x_multa_resid,
                          x_MPCerrado_resid,
                          CodIBGE)
  
  
  # AreaPlant, AreasPast, POP, Garimpo
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  # AreaPlant, AreasPast, POP, Garimpo, Multa, AreaProt, MPCerrado
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                             
                                             "x_multa_resid",
                                             "x_areaprot_resid",
                                             "x_lei_resid",
                                             "x_MPCerrado_resid"
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                         
                                         "x_multa_resid",
                                         "x_areaprot_resid",
                                         "x_lei_resid",
                                         "x_MPCerrado_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
       }

etable(feols(lDesmatPriMA~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_ConselhoLimpo, cluster=~CodIBGE, db_mata), stage = 2, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPriINPE="DesmatPriMA",
                lAreaPast = "AreaPast",                                         
                "lag_ConselhoLimpo"="Lag_COMUMArestrito",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()



### Normalizado ###

etable(feols(DesmatPriMAn~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_Conselho, cluster=~CodIBGE, db_mata), stage = 1, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPri="DesmatPriMA",
                lAreaPast = "AreaPast",                                         
                "lag_Conselho"="Lag_COMUMA",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()


{
  y = feols(DesmatPriMAn~1|CodIBGE+UF^Year, db_mata)
  y_resid= y$residuals
  
  x_lagendog = feols(lag_lPeCreal~1|CodIBGE+UF^Year, db_mata)
  x_lagendog_resid= x_lagendog$residuals
  
  z_laginstrument = feols(lag_Conselho~1|CodIBGE+UF^Year, db_mata)
  z_laginstrument_resid = z_laginstrument$residuals
  
  x_areapast = feols(lAreaPast~1|CodIBGE+UF^Year, db_mata)
  x_areapast_resid =  x_areapast$residuals
  
  x_garimpo = feols(lGarimpo~1|CodIBGE+UF^Year, db_mata)
  x_garimpo_resid =  x_garimpo$residuals
  
  x_vaagro = feols(lVAagro~1|CodIBGE+UF^Year, db_mata)
  x_vaagro_resid =  x_vaagro$residuals
  
  x_areaplant = feols(lAreaPlant~1|CodIBGE+UF^Year, db_mata)
  x_areaplant_resid =  x_areaplant$residuals
  
  x_pop = feols(lPOP~1|CodIBGE+UF^Year, db_mata)
  x_pop_resid =  x_pop$residuals
  
  x_conselho = feols(lag_Conselho~1|CodIBGE+UF^Year, db_mata)
  x_conselho_resid =  x_conselho$residuals
  
  x_lei = feols(lag_LeiFloresta~1|CodIBGE+UF^Year, db_mata)
  x_lei_resid =  x_lei$residuals
  
  x_areaprot = feols(lag_D_AreaProt~1|CodIBGE+UF^Year, db_mata)
  x_areaprot_resid =  x_areaprot$residuals
  
  x_multa = feols(lag_lFines~1|CodIBGE+UF^Year, db_mata)
  x_multa_resid =  x_multa$residuals
  
  x_MPCerrado = feols(lag_MPCerrado~1|CodIBGE+UF^Year, db_mata)
  x_MPCerrado_resid =  x_MPCerrado$residuals
  
  CodIBGE = db_mata$CodIBGE
  
  ## data frame para inputar no teste
  db_AR_Test = data.frame(y_resid,
                          x_lagendog_resid,
                          z_laginstrument_resid,
                          x_areapast_resid,
                          x_areaplant_resid,
                          x_vaagro_resid,
                          x_pop_resid,
                          x_garimpo_resid,
                          x_conselho_resid,
                          x_lei_resid,
                          x_areaprot_resid,
                          x_multa_resid,
                          x_MPCerrado_resid,
                          CodIBGE)
  
  
  # AreaPlant, AreasPast, POP, Garimpo
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  # AreaPlant, AreasPast, POP, Garimpo, Multa, AreaProt, MPCerrado
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                             
                                             "x_multa_resid",
                                             "x_areaprot_resid",
                                             "x_lei_resid",
                                             "x_MPCerrado_resid"
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                         
                                         "x_multa_resid",
                                         "x_areaprot_resid",
                                         "x_lei_resid",
                                         "x_MPCerrado_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
       }

etable(feols(DesmatPriMAn~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_Conselho, cluster=~CodIBGE, db_mata), stage = 2, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPri="DesmatPriMA",
                lAreaPast = "AreaPast",                                         
                "lag_ConselhoLimpo"="Lag_COMUMArestrito",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()


### IHS no Y ###

etable(feols(DesmatPriMAihs~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_Conselho, cluster=~CodIBGE, db_mata), stage = 1, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPri="DesmatPriMA",
                lAreaPast = "AreaPast",                                         
                "lag_Conselho"="Lag_COMUMA",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()


{
  y = feols(DesmatPriMAihs~1|CodIBGE+UF^Year, db_mata)
  y_resid= y$residuals
  
  x_lagendog = feols(lag_lPeCreal~1|CodIBGE+UF^Year, db_mata)
  x_lagendog_resid= x_lagendog$residuals
  
  z_laginstrument = feols(lag_Conselho~1|CodIBGE+UF^Year, db_mata)
  z_laginstrument_resid = z_laginstrument$residuals
  
  x_areapast = feols(lAreaPast~1|CodIBGE+UF^Year, db_mata)
  x_areapast_resid =  x_areapast$residuals
  
  x_garimpo = feols(lGarimpo~1|CodIBGE+UF^Year, db_mata)
  x_garimpo_resid =  x_garimpo$residuals
  
  x_vaagro = feols(lVAagro~1|CodIBGE+UF^Year, db_mata)
  x_vaagro_resid =  x_vaagro$residuals
  
  x_areaplant = feols(lAreaPlant~1|CodIBGE+UF^Year, db_mata)
  x_areaplant_resid =  x_areaplant$residuals
  
  x_pop = feols(lPOP~1|CodIBGE+UF^Year, db_mata)
  x_pop_resid =  x_pop$residuals
  
  x_conselho = feols(lag_Conselho~1|CodIBGE+UF^Year, db_mata)
  x_conselho_resid =  x_conselho$residuals
  
  x_lei = feols(lag_LeiFloresta~1|CodIBGE+UF^Year, db_mata)
  x_lei_resid =  x_lei$residuals
  
  x_areaprot = feols(lag_D_AreaProt~1|CodIBGE+UF^Year, db_mata)
  x_areaprot_resid =  x_areaprot$residuals
  
  x_multa = feols(lag_lFines~1|CodIBGE+UF^Year, db_mata)
  x_multa_resid =  x_multa$residuals
  
  x_MPCerrado = feols(lag_MPCerrado~1|CodIBGE+UF^Year, db_mata)
  x_MPCerrado_resid =  x_MPCerrado$residuals
  
  CodIBGE = db_mata$CodIBGE
  
  ## data frame para inputar no teste
  db_AR_Test = data.frame(y_resid,
                          x_lagendog_resid,
                          z_laginstrument_resid,
                          x_areapast_resid,
                          x_areaplant_resid,
                          x_vaagro_resid,
                          x_pop_resid,
                          x_garimpo_resid,
                          x_conselho_resid,
                          x_lei_resid,
                          x_areaprot_resid,
                          x_multa_resid,
                          x_MPCerrado_resid,
                          CodIBGE)
  
  
  # AreaPlant, AreasPast, POP, Garimpo
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  # AreaPlant, AreasPast, POP, Garimpo, Multa, AreaProt, MPCerrado
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                             
                                             "x_multa_resid",
                                             "x_areaprot_resid",
                                             "x_lei_resid",
                                             "x_MPCerrado_resid"
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                         
                                         "x_multa_resid",
                                         "x_areaprot_resid",
                                         "x_lei_resid",
                                         "x_MPCerrado_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
       }

etable(feols(DesmatPriMAihs~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta)|
               CodIBGE+UF^Year|lag_lPeCreal~lag_Conselho, cluster=~CodIBGE, db_mata), stage = 2, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPri="DesmatPriMA",
                lAreaPast = "AreaPast",                                         
                "lag_ConselhoLimpo"="Lag_COMUMArestrito",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()

## OLS 

etable(feols(lDesmatPriMA~csw(lag_lPeCreal,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta)|
               CodIBGE+UF^Year, cluster=~CodIBGE, db_mata)) %>% View()

## reduzida 

etable(feols(lDesmatPriMA~csw(lag_Conselho,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta)|
               CodIBGE+UF^Year, cluster=~CodIBGE, db_mata)) %>% View()



# Exercício Extra: Gestão Ambiental ----------------------------------------------------


etable(feols(lDesmatPriMA~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta)|
               CodIBGE+UF^Year|lag_lGestAreal~lag_Conselho, cluster=~CodIBGE, db_mata), stage = 1, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPriINPE="DesmatPriMA",
                lAreaPast = "AreaPast",                                         
                "lag_Conselho"="Lag_COMUMA",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()


{
  y = feols(lDesmatPriMA~1|CodIBGE+UF^Year, db_mata)
  y_resid= y$residuals
  
  x_lagendog = feols(lag_lGestAreal~1|CodIBGE+UF^Year, db_mata)
  x_lagendog_resid= x_lagendog$residuals
  
  z_laginstrument = feols(lag_Conselho~1|CodIBGE+UF^Year, db_mata)
  z_laginstrument_resid = z_laginstrument$residuals
  
  x_areapast = feols(lAreaPast~1|CodIBGE+UF^Year, db_mata)
  x_areapast_resid =  x_areapast$residuals
  
  x_garimpo = feols(lGarimpo~1|CodIBGE+UF^Year, db_mata)
  x_garimpo_resid =  x_garimpo$residuals
  
  x_vaagro = feols(lVAagro~1|CodIBGE+UF^Year, db_mata)
  x_vaagro_resid =  x_vaagro$residuals
  
  x_areaplant = feols(lAreaPlant~1|CodIBGE+UF^Year, db_mata)
  x_areaplant_resid =  x_areaplant$residuals
  
  x_pop = feols(lPOP~1|CodIBGE+UF^Year, db_mata)
  x_pop_resid =  x_pop$residuals
  
  x_conselho = feols(lag_Conselho~1|CodIBGE+UF^Year, db_mata)
  x_conselho_resid =  x_conselho$residuals
  
  x_lei = feols(lag_LeiFloresta~1|CodIBGE+UF^Year, db_mata)
  x_lei_resid =  x_lei$residuals
  
  x_areaprot = feols(lag_D_AreaProt~1|CodIBGE+UF^Year, db_mata)
  x_areaprot_resid =  x_areaprot$residuals
  
  x_multa = feols(lag_lFines~1|CodIBGE+UF^Year, db_mata)
  x_multa_resid =  x_multa$residuals
  
  x_MPCerrado = feols(lag_MPCerrado~1|CodIBGE+UF^Year, db_mata)
  x_MPCerrado_resid =  x_MPCerrado$residuals
  
  CodIBGE = db_mata$CodIBGE
  
  ## data frame para inputar no teste
  db_AR_Test = data.frame(y_resid,
                          x_lagendog_resid,
                          z_laginstrument_resid,
                          x_areapast_resid,
                          x_areaplant_resid,
                          x_vaagro_resid,
                          x_pop_resid,
                          x_garimpo_resid,
                          x_conselho_resid,
                          x_lei_resid,
                          x_areaprot_resid,
                          x_multa_resid,
                          x_MPCerrado_resid,
                          CodIBGE)
  
  
  # AreaPlant, AreasPast, POP, Garimpo
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid"
                                             
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
  # AreaPlant, AreasPast, POP, Garimpo, Multa, AreaProt, MPCerrado
  AR.Test = anderson_rubin_test(outcome = "y_resid",
                                endogenous = "x_lagendog_resid",
                                instruments = "z_laginstrument_resid",
                                controls = c("x_areaplant_resid",
                                             "x_areapast_resid",
                                             "x_pop_resid",
                                             "x_garimpo_resid",
                                             
                                             "x_multa_resid",
                                             "x_areaprot_resid",
                                             "x_lei_resid",
                                             "x_MPCerrado_resid"
                                             
                                ), vcovCL,
                                data =  db_AR_Test, cluster = "CodIBGE") %>% print()
  
  grid_beta = seq(-100, 100, 0.01)
  AR.CI = anderson_rubin_ci(outcome = "y_resid",
                            endogenous = "x_lagendog_resid",
                            instruments = "z_laginstrument_resid",
                            controls = c("x_areaplant_resid",
                                         "x_areapast_resid",
                                         "x_pop_resid",
                                         "x_garimpo_resid",
                                         
                                         "x_multa_resid",
                                         "x_areaprot_resid",
                                         "x_lei_resid",
                                         "x_MPCerrado_resid")
                            ,
                            vcovCL, db_AR_Test,
                            cluster = "CodIBGE",
                            grid_beta = grid_beta,
                            confidence = 0.90) %>% show()
  
       }

etable(feols(lDesmatPriMA~csw(1,lAreaPlant,lAreaPast,lPOP,lGarimpo,lag_lFines,lag_D_AreaProt,lag_LeiFloresta)|
               CodIBGE+UF^Year|lag_lGestAreal~lag_Conselho, cluster=~CodIBGE, db_mata), stage = 2, fitstat = ~ivwald+ivwald.p+n,
       dict = c(lDesmatPri="DesmatPriMA",
                lAreaPast = "AreaPast",                                         
                "lag_ConselhoLimpo"="Lag_COMUMA",
                "lag_lPeCreal" = "Lag_PeC",
                lAreaPlant="AreaPlant",                                     
                lPOP="POP",
                lGarimpo="Garimpo",
                lag_lFines = "Lag_Multas",
                lag_D_AreaProt = "Lag_DAreaProt",
                lag_LeiFloresta = "Lag_LeiFloresta",
                CodIBGE="Municipio"), 
       order = c("Lag_COMUMA","Lag_PeC")) %>% View()


