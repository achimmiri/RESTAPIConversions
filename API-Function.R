library(jsonlite)
library("hmdbQuery")
library("digest")
library("RMassBank")
library(devtools)
library("MetaboAnalystR")
library("janitor")

##### REST API Functions
##########################################################
##########################################################
PuInKtoSM<-function(getINK)
{
  ########## Functions return canoncial smiles	
  url<- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/"
  out<-tryCatch({jsonlite::fromJSON(paste0(url,getINK, "/JSON"))} ,error = function(x) {return(NA)})
  prop.names  <-tryCatch({out$PC_Compounds$props[[1]][[1]]},error = function(x) {return(NA)})
  prop.values <- tryCatch({out$PC_Compounds$props[[1]][[2]]},error = function(x) {return(NA)})
  sm <-tryCatch({grep("smiles", prop.names[,"label"], ignore.case = TRUE)},error = function(x) {return(NA)})
  csmiles<-c()
  if(length(sm) >= 1) {
    can <- tryCatch({grep("canonical", prop.names[,"name"], ignore.case = TRUE)},error= function(x) {return(NA)})
    can1<-tryCatch({prop.values[sm[1],"sval"]},warning= function(x) {return(NA)})
    csmiles<-c(csmiles,can1)
  }else{
    csmiles<-c(csmiles,NA)
  }
  return(csmiles)

}
############################################################
############################################################
PuInKtoSM1<-function(getINK)
{
  ###### Functions return Isomeric smiles
  url<- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/"
  out<-tryCatch({jsonlite::fromJSON(paste0(url,getINK, "/JSON"))}, error = function(x) {return(NA)})
  prop.names  <-tryCatch({out$PC_Compounds$props[[1]][[1]]}, error = function(x) {return(NA)})
  prop.values <- tryCatch({out$PC_Compounds$props[[1]][[2]]}, error = function(x) {return(NA)})
  sm <-tryCatch({grep("smiles", prop.names[,"label"], ignore.case = TRUE)}, error = function(x) {return(NA)})
  csmiles<-c()
  if(length(sm) == 2) {
    can <- tryCatch({grep("canonical", prop.names[,"name"], ignore.case = TRUE)}, error = function(x) {return(NA)})
    can1<-tryCatch({prop.values[sm[2],"sval"]},warning= function(x) {return(NA)})
    csmiles<-c(csmiles,can1)
  }else{
    csmiles<-c(csmiles,NA)
  }
  return(csmiles)

}
##########################################################################
###########################################################################
PuCIDtoEM<-function(getCID)
{
  ##### Function returns Compound ID to exact mass
  #######################
  URL="https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/"
  URL1=paste0(URL,getCID, "/JSON/?response_type=display")
############################
  ########################
  data=tryCatch({jsonlite::fromJSON(URL1)} ,error = function(x) {return(NA)})
  ########################
  data1<-tibble::enframe(unlist(data))
  data2<-as.data.frame(data1)
  ##########################
  inVa<-which(data2$name %in% "Record.Section.Section.Section.Information.Value.StringWithMarkup.String")
  ##########################
  TEST<-sapply(data2, "[", inVa)
  TEST1<-grep("InChI=",TEST)
  TEST2<-TEST1[1]
  InchI<-TEST[TEST2]
  inchikey<-TEST[TEST2+1]
  #########################
  url<- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/"
  #########################
  out<-tryCatch({jsonlite::fromJSON(paste0(url,inchikey, "/JSON"))} ,error = function(x) {return(NA)})
  #######################
  EMV<-out$PC_Compounds$props[[1]][22,]$value$sval
  #######################
  EMV1<-c()
  ######################
  if(!sjmisc::is_empty(EMV))
  {
    EMV1<-c(EMV1,EMV)
  }else{
    EMV1<-c(EMV1,0)
  }
  #####################
  return(EMV1)
  ####################
}
###########################################################################
###########################################################################
PuCAStoOI<-function(getCAS)
{
  ###CAS: 328-50-7
  ### Function returns CAS to other formats like Inchi,Inchikey,Smiles,CompoundID	
  getCAS1<-stringr::str_replace(getCAS,pattern='CAS:',replacement ="")
  getCAS2<-stringr::str_trim(getCAS1)
  url<- "http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
  out<-tryCatch({jsonlite::fromJSON(paste0(url,getCAS2, "/property/InChIKey"))}, warning = function(x) {return(NA)})
  out1<-tryCatch({jsonlite::fromJSON(paste0(url,getCAS2, "/property/CanonicalSMILES"))}, warning = function(x) {return(NA)})
  out2<-tryCatch({jsonlite::fromJSON(paste0(url,getCAS2, "/property/InChI"))}, warning = function(x) {return(NA)})
  ###########################
  OIN<-tryCatch({out2$PropertyTable$Properties$InChI},error=function(cond){message("List value is empty")})
  OIK<-tryCatch({out$PropertyTable$Properties$InChIKey},error=function(cond){message("List value is empty")})
  OSM<-tryCatch({out1$PropertyTable$Properties$CanonicalSMILES},error=function(cond){message("List value is empty")})
  OCID<-tryCatch({out1$PropertyTable$Properties$CID},error=function(cond){message("List value is empty")})
  ############################
  ##returning Inchi,InchiKey,Smiles,CompoundID...in order
  ###############################
  return(c(OIN[1],OIK[1],OSM[1],OCID[1]))
  ###############################
}
#########################################################################
#########################################################################
ConvPCIDtoOCN<-function(getPCID)
{
  #### Function returns CompoundID to other formats like inchikey,smiles,pubchem_cid,exactmass, formula	
  url<- "https://www.metabolomicsworkbench.org/rest/compound/pubchem_cid/"
  out<-tryCatch({jsonlite::fromJSON(paste0(url,getPCID, "/all"))}, warning = function(x) {return(NA)})
  #########################
  OIK<-tryCatch({out$inchi_key},error=function(cond){message("Inchikey value is empty")})
  OSM<-tryCatch({out$smiles},error=function(cond){message("smiles value is empty")})
  OCID<-tryCatch({out$pubchem_cid},error=function(cond){message("Pubchem CID value is empty")})
  OEM<-tryCatch({out$exactmass},error=function(cond){message("Exact mass value is empty")})
  OFOR<-tryCatch({out$formula},error=function(cond){message("FORMULA value is empty")})
  ################################
  return(c(OIK[1],OSM[1],OCID[1],OEM[1],OFOR[1]))
  ################################

}
##########################################################################
###########################################################################
ConvINKtoOID<-function(getINK)
{
  #### Function returns inchi_key to other formats like smiles,pubchem_cid,exactmass,formula	
  url<- "https://www.metabolomicsworkbench.org/rest/compound/inchi_key/"
  out<-tryCatch({jsonlite::fromJSON(paste0(url,getINK, "/all"))}, warning = function(x) {return(NA)})
  #########################
  OIK<-tryCatch({out$inchi_key},error=function(cond){message("Inchikey value is empty")})
  OSM<-tryCatch({out$smiles},error=function(cond){message("smiles value is empty")})
  OCID<-tryCatch({out$pubchem_cid},error=function(cond){message("Pubchem CID value is empty")})
  OEM<-tryCatch({out$exactmass},error=function(cond){message("Exact mass value is empty")})
  OFOR<-tryCatch({out$formula},error=function(cond){message("FORMULA value is empty")})
  ###########################
  return(c(OIK[1],OSM[1],OCID[1],OEM[1],OFOR[1]))
  ###########################

}
##########################################################################
##########################################################################
ClassSmilesToOntolgy<-function(getSMILE)
{
  #####################
  url="https://gnps-structure.ucsd.edu/classyfire?smiles="
  url1=paste0(url,getSMILE)
  out=jsonlite::fromJSON(url1)
  res=do.call(paste, c(as.list(tryCatch({rev(out$ancestors)},warning=function(cond){message("Classifier could not fecth the information")})), sep = ","))
  ###########
  return(res)
  ############
}
########################################################################
########################################################################
ConvHMDBtoOCN<-function(getHMDB)
{
  url<- "https://www.metabolomicsworkbench.org/rest/compound/hmdb_id/"
  out<-tryCatch({jsonlite::fromJSON(paste0(url,getHMDB, "/all"))}, warning = function(x) {return(NA)})
  #########################
  OIK<-tryCatch({out$inchi_key},error=function(cond){message("Inchikey value is empty")})
  OSM<-tryCatch({out$smiles},error=function(cond){message("smiles value is empty")})
  OCID<-tryCatch({out$pubchem_cid},error=function(cond){message("Pubchem CID value is empty")})
  OEM<-tryCatch({out$exactmass},error=function(cond){message("Exact mass value is empty")})
  OFOR<-tryCatch({out$formula},error=function(cond){message("FORMULA value is empty")})
  return(c(OIK[1],OSM[1],OCID[1],OEM[1],OFOR[1]))
  ############################

  ###########################
}
#########################################################################
#########################################################################
