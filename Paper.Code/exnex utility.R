ix_sel = function(sel,vnames,matrix2list=FALSE) {
  ns <- nchar(sel)
  # check for exact name matching
  match1 <- sapply(vnames, function(e)
    (substring(paste(e, sel, sep = ""), 1, ns) == sel)  &
      (nchar(e) == ns)
  )
  # check for matching matrix of vector, i.e., with "[" following sel
  match2 <- sapply(vnames, function(e) 
    (substring(paste(e, sel, sep = ""), 1, ns) == sel)  &
      (substring(e, ns + 1, ns + 1) == "[") 
    # &
    #   (!is.element(",",strsplit(e,split="")[[1]]))
  )
  if (matrix2list)
    match3 <- sapply(vnames, function(e) 
      (substring(paste(e, sel, sep = ""), 1, ns) == sel)  &
        (substring(e, ns + 1, ns + 1) == "[") &
        (is.element(",",strsplit(e,split="")[[1]]))
    )
  else
    match3 <- sapply(vnames,function(e) FALSE)
  
  out = NULL
  
  if (any(match1)) {
    out = (1:length(vnames))[match1]
    out = list(out)
    names(out) = sel
  }
  
  if (any(match2)) {
    out = (1:length(vnames))[match2]
    out = list(out)
    names(out) = sel
  }
  
  if (any(match3)) {
    sss = sapply(vnames,function(e) strsplit(e,split=",")[[1]][1])
    sss=sss[match3]
    sss = unique(sss)
    ssslab=lapply(sss,function(e) strsplit(e,split="")[[1]])
    ssslab=lapply(ssslab,function(e) paste(e[e != "["],collapse=""))
    
    out = lapply(sss,function(ee) {
      match3 <- sapply(vnames,
                       function(e)
                         substring(e,1,nchar(ee)) == ee
      )
      return(  (1:length(vnames))[match3])
    }
    )
    names(out) = ssslab
  }
  
  if (is.list(out)) {
    
    outvec = as.vector(sapply(out,function(e) vnames[e]))
    g1 = sapply(outvec,function(e) match(",",strsplit(e,"")[[1]]))
    g2 = sapply(1:length(outvec),function(e) substring(outvec[e],g1[e]+1,g1[e]+1))
    
    if (length(unique(g2))==1)
      out = as.vector(sapply(out,function(e) e))
  }
  
  return(out)
}

sim_jagswb <-  function(pars,
                        out,
                        matrix2list = TRUE) {
  allout = list()
  for (j in 1:length(pars)) {
    
    parsj = pars[j]
    out = cbind(out)
    vnames <- colnames(out)
    out.ix = ix_sel(parsj, vnames, matrix2list = matrix2list)
    
    if (is.list(out.ix))
      out.simsum = lapply(out.ix, function(e)
        out[, e])
    else
      out.simsum = out[, out.ix]
    
    allout[[j]] = out.simsum
  }
  
  names(allout) = pars
  
  return(allout)
}


sum_jagswb <-  function(pars,
                        out,
                        matrix2list = TRUE,
                        sum.cols =  c("mean", "sd", "2.5%", "50%", "97.5%")) {
    allout = list()
  
  for (j in 1:length(pars)) {
    parsj = pars[j]
    out = rbind(out)
    vnames <- rownames(out)
    out.ix = ix_sel(parsj, vnames, matrix2list = matrix2list)
    if (is.list(out.ix))
      out.simsum = lapply(out.ix, function(e)
        out[e, sum.cols])
    else
      out.simsum = rbind(out)[out.ix, sum.cols]
    
    allout[[j]] = out.simsum
    
  }
  names(allout) = pars
  
  allout = lapply(allout, function(e) {
    if (is.list(e) & length(e) == 1)
      return(e[[1]])
    else
      return(e)
  })
  
  return(allout)
}


