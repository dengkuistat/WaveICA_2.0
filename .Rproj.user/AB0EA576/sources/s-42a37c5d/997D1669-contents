#' @title WaveICA2.0
#' @description Removing batch effects for metabolomics data without
#' using batch information
#' @author Kui Deng
#' \email{dengkui_stat@163.com}
#' @param data Sample-by-Peak matrix untargeted metabolomics data.
#' @param wf  wavelet functions. The default is "haar".
#' @param Injection_Order  Injection Order of samples.
#' @param alpha The trade-off value between the independence of samples and those
#' of variables in ICA, and should be between 0 and 1.The default is 0.
#' @param Cutoff The threshold of the variation explained by the injection order for
#' independent components. It should be between 0 and 1.
#' @param K The maximal component that ICA decomposes. The default is 10.

#' @return A list that contains the clean data.
#' @export
#' @examples
#' \dontrun{
#' ## load the demo data
#'data(Amide_data, package = "WaveICA")

#'### Data management for original data
#' Amide_data_order <- Amide_data[order(Amide_data$Injection_order),]
#' data_Amide_order <- Amide_data_order[,-c(1:3)]
#' ### Generating group and batch
#' Amide_data_order$group[Amide_data_order$group=="QC"] <- 2
#' group_zong_Amide<-as.numeric(Amide_data_order$group)
#' Injection_order<-Amide_data_order$Injection_order

#'############### run WaveICA2.0
#'data_wave2.0_reconstruct_Amide <- WaveICA_2.0(data=data_Amide_order,Injection_Order=Injection_order,alpha=0,Cutoff=0.1,K=10)
#'data_Amide_zong_wave2.0 <- data_wave2.0_reconstruct_Amide$data_wave


#'data_Amide_sample_wave2.0 <- data_Amide_zong_wave2.0[group_zong_Amide!=2,]
#'data_Amide_qc_wave2.0<-data_Amide_zong_wave2.0[group_zong_Amide==2,]

#' }


WaveICA_2.0<-function(data,wf="haar",Injection_Order,alpha,Cutoff,K){
  ### Wavelet Decomposition
  library(waveslim)
  library(parallel)
  library(ica)
  library(mgcv)
  level<-floor(log(nrow(data),2))
  if (is.null(colnames(data))){
    stop("data must have colnames")
  }
  coef<-list()
  for (k in 1:(level+1)){
    coef[[k]] <-matrix(NA,nrow(data),ncol(data))
  }
  for (j in 1:ncol(data)){
    cat(paste("######Decomposition",j,"########\n"))
    data_temp<-data[,j]
    x_modwt<-modwt(data_temp,wf=wf,n.levels =level)
    for (k in 1:(level+1)){
      coef[[k]][,j]<-x_modwt[[k]]
    }
  }

  ##### ICA
  index<-level+1
  data_wave_ICA<-list()
  for (i in (1:index)){
    cat(paste("######### ICA",i,"#############\n"))
    data_coef<-coef[[i]]
    data_coef_ICA<-unbiased_stICA(X=t(data_coef),k=K,alpha)
    B <- data_coef_ICA$B
    A <- data_coef_ICA$A

    B <- as.data.frame(B)

    ## Gam
    corr <- mclapply(B,function(x){
      corr <- gam(x~s(Injection_Order))
      corr_summary <- summary(corr)
      corr_r <- corr_summary$r.sq
      return(corr_r)
    })
    corr <- unlist(corr)
    label <- which(corr>=Cutoff)
    B_new <- B[,label,drop=F]
    A_new <- A[,label,drop=F]
    Xn = data_coef-t(A_new %*% t(B_new))

    data_wave_ICA[[i]]<-Xn
  }

  ### Wavelet Reconstruction
  index<-ncol(data)
  index1<-length(data_wave_ICA)
  data_coef<-matrix(NA,nrow(data_wave_ICA[[1]]),index1)
  data_wave<-matrix(NA,nrow(data_wave_ICA[[1]]),ncol(data_wave_ICA[[1]]))
  for (i in 1:index){
    cat(paste("######Reconstruction",i,"########\n"))
    for (j in 1:index1){
      data_coef[,j]<-data_wave_ICA[[j]][,i]
    }
    data_temp<-data[,i]
    data_coef<-as.data.frame(data_coef)
    colnames(data_coef)<-c(paste("d",1:(index1-1),sep=""),paste("s",(index1-1),sep=""))
    y<-as.list(data_coef)
    attributes(y)$class<-"modwt"
    attributes(y)$wavelet<-wf
    attributes(y)$boundary<-"periodic"
    data_wave[,i]<-imodwt(y)+mean(data_temp)
  }
  rownames(data_wave)<-rownames(data)
  colnames(data_wave)<-colnames(data)
  return(list(data_wave=data_wave))
}





