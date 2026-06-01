#' @title Format data table into the input of tempted
#' @description This function applies a variety of transformations to the read counts and
#' format the sample by feature table and meta data into a data list
#' that can be used as the input of \code{\link{tempted}} and \code{\link{svd_centralize}}.
#' For data that are not read counts, or data that are not microbiome data,
#' the user can apply their desired transformation to the data before formatting into list.
#' @param featuretable A sample by feature matrix.
#' @param timepoint The time stamp of each sample, matched with the rows of \code{featuretable}.
#' @param subjectID The subject ID of each sample, matched with the rows of \code{featuretable}.
#' @param threshold A threshold for feature filtering for microbiome data.
#' Features with zero value percentage > threshold will be excluded. Default is 0.95.
#' @param pseudo A small number to add to all the counts before
#' normalizing into proportions and log transformation.
#' Default is 1/2 of the smallest non-zero value that is specific for each sample.
#' This pseudo count is added for \code{transform=c("logcomp", "clr", "logit")}.
#' @param transform The transformation applied to the data.
#' \code{"logcomp"} for log of compositions.
#' \code{"comp"} for compositions.
#' \code{"ast"} for arcsine squared transformation.
#' \code{"clr"} for central log ratio transformation.
#' \code{"lfb"} for log 2 fold change over baseline (first time point) transformation.
#' \code{"logit"} for logit transformation.
#' \code{"none"} for no transformation.
#' Default \code{transform="clr"} is recommended for microbiome data.
#' For data that are already transformed, use \code{transform="none"}.
#' @return A length n list of matrices as the input of \code{\link{tempted}} and \code{\link{svd_centralize}}.  Each matrix represents a subject, with columns representing samples from this subject, the first row representing the sampling time points, and the following rows representing the feature values.
#' @seealso Examples can be found in \code{\link{tempted}}.
#' @export
#' @md
format_tempted <- function(featuretable, timepoint, subjectID,
                           threshold=0.95, pseudo=NULL, transform="clr"){
  ntm <- which(table(subjectID)==1)
  if(length(ntm)>0)
    stop(paste('Please remove these subjects with only one time point:',
               paste(names(ntm), collapse=', ')))
  if (length(subjectID)!=nrow(featuretable))
    stop('length of subjectID does not match featuretable!')
  if (length(timepoint)!=nrow(featuretable))
    stop('length of timepoint does not match featuretable!')
  # get pseudo count
  if (is.null(pseudo) & (transform %in% c("clr", "logcomp", "logit"))){
    pseudo <- apply(featuretable, 1, function(x){
      min(x[x!=0])/2
    })
  }
  # keep taxon that has non-zeros in >1-threshold samples
  featuretable <- featuretable[,colMeans(featuretable==0)<=threshold]
  if(transform=='logcomp' | transform=="lfb"){
    featuretable <- featuretable+pseudo
    featuretable <- t(log(featuretable/rowSums(featuretable)))
  }else if(transform=='comp'){
    featuretable <- featuretable
    featuretable <- t(featuretable/rowSums(featuretable))
  }else if(transform=='ast'){
    featuretable <- featuretable
    featuretable <- t(asin(sqrt(featuretable/rowSums(featuretable))))
  }else if(transform=='clr'){
    featuretable <- featuretable+pseudo
    featuretable <- log(featuretable/rowSums(featuretable))
    featuretable <- t(featuretable-rowMeans(featuretable))
  }else if(transform=="lfb"){
    featuretable <- featuretable+pseudo
    featuretable <- t(log(featuretable/rowSums(featuretable), 2))
  }else if(transform=='logit'){
    featuretable <- featuretable+pseudo
    featuretable <- t(featuretable/rowSums(featuretable))
    featuretable <- log(featuretable/(1-featuretable))
  }else if(transform=='none'){
    featuretable <- t(featuretable)
  }else{
    message('Input transformation method is wrong! logcomp is applied instead')
    featuretable <- featuretable+pseudo
    featuretable <- t(log(featuretable/rowSums(featuretable)))
  }
  featuretable <- rbind(timepoint, featuretable)
  rownames(featuretable)[1] <- 'timepoint'
  subID <- unique(subjectID)
  nsub <- length(subID)

  # construct list of data matrices, each element representing one subject
  datlist <- vector("list", length = nsub)
  names(datlist) <- subID

  # Each slice represents an individual (unequal sized matrix).
  for (i in 1:nsub){
    datlist[[i]] <- featuretable[, subjectID==subID[i]]
    datlist[[i]] <- datlist[[i]][,order(datlist[[i]][1,])]
    datlist[[i]] <- datlist[[i]][,!duplicated(datlist[[i]][1,])]
    if(transform=="lfb"){
      datlist[[i]] <- datlist[[i]][,-1, drop=FALSE] - datlist[[i]][,1]
    }
  }
  return(datlist)
}
