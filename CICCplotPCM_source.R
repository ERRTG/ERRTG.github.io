
gamma_r_rec_pcm <- function(pars, r, par.grp){
  if (r == 0) return( 1 )
  if (r > length(pars)| r<0) return( 0 )
  if (r != 0 | r <= length(pars)) return( sum(exp(pars[par.grp==1]) * sapply(1:length(pars[par.grp==1]), FUN = function(x){gamma_r_rec_pcm( pars = pars[par.grp!=1], r = r-x, par.grp = par.grp[par.grp!=1]-1)})) +  gamma_r_rec_pcm(pars = pars[par.grp!=1], r = r, par.grp = par.grp[par.grp!=1]-1))
}

# PCM version af CICCplot funktionen
CICCplot <- function(model, which.item = 1, lower.groups = NULL, all.items = F, grid.items = T, error.bar = F,  color = NULL){
  
  data <- model$X
  betas <- model$betapar
  #print(betas)
  
  n.itemcat <- apply(data, 2, FUN = function(x){max(x, na.rm = T) - min(x, na.rm = T)})
  par.itemgrp <- rep(1:ncol(data), times = n.itemcat)
  
  phi <- NULL
  for(x in 1:ncol(model$X)){phi <- c(phi,-cumsum(betas[par.itemgrp==x]))}
  
  if (is.double(lower.groups)|is.integer(lower.groups)) {if(any(lower.groups>length(betas))){ message("\n lower.group index greater that maximum possible score \n"); return(NA) }}
  if (any(which.item>dim(data)[2])){ message("\n which.item index greater that number of items in the model \n"); return(NA) }
  
  
  if( length(which.item)==1 & is.double(which.item)){
    
    Tot.val <- 0:length(phi)
    
    exp.val <- sapply(Tot.val, FUN = function(R){ 
      l <- par.itemgrp[par.itemgrp!=which.item]
      par.itemgrp_noitem <- ifelse( l > which.item, l-1, l)
      g1 <- gamma_r_rec_pcm(-phi, R, par.itemgrp)
      return( sum( sapply( 1:sum(par.itemgrp==which.item), FUN = function(X){ 
        g2 <- gamma_r_rec_pcm(-phi[par.itemgrp!=which.item], R-X, par.itemgrp_noitem)
        return( X*exp(-phi[par.itemgrp==which.item][X])*g2/g1)})))
    })
    
    data_exp <- data.frame(Tot.val, exp.val)
    
    if (!is.double(lower.groups) & !is.integer(lower.groups)){
      if (lower.groups == "all"){
        
        Tot.val_grp <- 0:length(phi)
        obs.val_grp <- sapply(Tot.val_grp, FUN = function(x){ mean( data[which(rowSums(data) == x), which.item] )})
        var.val_grp <- sapply(Tot.val_grp, FUN = function(x){ var( data[which(rowSums(data) == x), which.item] )})
        n.val_grp <- sapply(Tot.val_grp, FUN = function(x){ length( data[which(rowSums(data) == x), which.item] )})
      }}
    
    if (is.double(lower.groups)|is.integer(lower.groups)){
      
      breaks <- sort(x = unique(c(floor(lower.groups), min(Tot.val))))
      n.groups <- length(breaks)
      
      Tot.val_grp <- rep(NA, times = n.groups)
      obs.val_grp <- rep(NA, times = n.groups)
      var.val_grp <- rep(NA, times = n.groups)
      n.val_grp <- rep(NA, times = n.groups)
      
      for (i in seq_along(breaks)){
        
        if(i != n.groups){
          
          obs.val_grp[i] <- mean( data[which(rowSums(data) %in% breaks[i]:(breaks[i+1]-1)), which.item])
          var.val_grp[i] <- var( data[which(rowSums(data) %in% breaks[i]:(breaks[i+1]-1)), which.item])
          n.val_grp[i] <- length( data[which(rowSums(data) %in% breaks[i]:(breaks[i+1]-1)), which.item])
          Tot.val_grp[i] <- mean( rowSums(data)[which(rowSums(data) %in% breaks[i]:(breaks[i+1]-1))])
          
        } else{
          
          obs.val_grp[i] <- mean( data[which(rowSums(data) %in% breaks[i]:max(Tot.val)), which.item])
          var.val_grp[i] <- var( data[which(rowSums(data) %in% breaks[i]:max(Tot.val)), which.item])
          n.val_grp[i] <- length( data[which(rowSums(data) %in% breaks[i]:max(Tot.val)), which.item])
          Tot.val_grp[i] <- mean( rowSums(data)[which(rowSums(data) %in% breaks[i]:max(Tot.val))])
        }
      }}
    
    data_obs <- data.frame(Tot.val_grp, obs.val_grp, var.val_grp, n.val_grp, CI.bound = NA)
    
    if (error.bar){ data_obs$CI.bound <- 1.96*sqrt(data_obs[,3]/data_obs[,4]) } 
    
    col <- c("Expected" = "black", "Observed" = "red")
    if (!is.null(color)) col <- c("Expected" = color$expected[1], "Observed" = color$observed[1])
    
    p <- ggplot(data = data_exp, aes(x = Tot.val, y= exp.val, color = "Expected")) + 
      geom_line(linetype = "dashed") + geom_point() +
      geom_point(data = data_obs, aes(x = Tot.val_grp, y = obs.val_grp, color = "Observed"), size = 2) +
      theme(legend.title = element_blank(), plot.title = element_text(size = 15,hjust = 0.5)) + scale_colour_manual(values = col) +
      ggtitle(paste0("Item ", which.item)) + xlab("Total Score") + ylab("Conditional Item-Score")  + 
      scale_x_continuous(breaks=Tot.val) + geom_errorbar(data = data_obs, 
                                                         aes(x = Tot.val_grp, y = obs.val_grp, ymin = obs.val_grp - CI.bound, ymax = obs.val_grp + CI.bound, color = "Observed"), 
                                                         width = 0.2, size = 1)
    
    P <- p
  }
  
  if (all.items | length(which.item)>1 ){
    
    which.item.arg <- which.item
    
    if(all.items){
      n.items <- dim(data)[2]
      pp <- list(rep(NA, n.items))
      ii <- 1:n.items}
    
    if((!all.items) & (length(which.item)>1) ){ 
      pp <- list(rep(NA, length(which.item)))
      ii <- which.item}
    
    j <- 1
    
    Tot.val <- 0:length(phi)
    
    # g_r <- sapply(Tot.val, FUN = function(x) gamma_r_rec(betas,x))
    
    for (k in ii) {
      
      which.item <- k 
      l <- par.itemgrp[par.itemgrp!=which.item]
      par.itemgrp_noitem <- ifelse( l > which.item, l-1, l)
      
      exp.val <- sapply(Tot.val, FUN = function(R){
        g1 <- gamma_r_rec_pcm(-phi, R, par.itemgrp)
        return( sum( sapply( 1:sum(par.itemgrp==which.item), FUN = function(X){ 
          g2 <- gamma_r_rec_pcm(-phi[par.itemgrp!=which.item], R-X, par.itemgrp_noitem)
          return( X*exp(-phi[par.itemgrp==which.item][X])*g2/g1)})))
      })
      data_exp <- data.frame(Tot.val, exp.val)
      
      if (!is.double(lower.groups) & !is.integer(lower.groups)){
        
        if (lower.groups == "all"){
          
          Tot.val_grp <- 0:length(phi)
          obs.val_grp <- sapply(Tot.val_grp, FUN = function(x){ mean( data[which(rowSums(data) == x), which.item] )})
          var.val_grp <- sapply(Tot.val_grp, FUN = function(x){ var( data[which(rowSums(data) == x), which.item] )})
          n.val_grp <- sapply(Tot.val_grp, FUN = function(x){ length( data[which(rowSums(data) == x), which.item] )})
        }}
      
      if (is.double(lower.groups)|is.integer(lower.groups)){
        
        breaks <- sort(x=unique(c(floor(lower.groups),min(Tot.val))))
        n.groups <- length(breaks)
        
        Tot.val_grp <- rep(NA, times = n.groups)
        obs.val_grp <- rep(NA, times = n.groups)
        var.val_grp <- rep(NA, times = n.groups)
        n.val_grp <- rep(NA, times = n.groups)
        
        for (i in seq_along(breaks)){
          
          if(i != n.groups){
            
            obs.val_grp[i] <- mean( data[which(rowSums(data) %in% breaks[i]:(breaks[i+1]-1)), which.item])
            var.val_grp[i] <- var( data[which(rowSums(data) %in% breaks[i]:(breaks[i+1]-1)), which.item])
            n.val_grp[i] <- length( data[which(rowSums(data) %in% breaks[i]:(breaks[i+1]-1)), which.item])
            Tot.val_grp[i] <- mean( rowSums(data)[which(rowSums(data) %in% breaks[i]:(breaks[i+1]-1))])
            
          } else{
            
            obs.val_grp[i] <- mean( data[which(rowSums(data) %in% breaks[i]:max(Tot.val)), which.item])
            var.val_grp[i] <- var( data[which(rowSums(data) %in% breaks[i]:max(Tot.val)), which.item])
            n.val_grp[i] <- length( data[which(rowSums(data) %in% breaks[i]:max(Tot.val)), which.item])
            Tot.val_grp[i] <- mean( rowSums(data)[which(rowSums(data) %in% breaks[i]:max(Tot.val))])
          }
        }}
      
      data_obs <- data.frame(Tot.val_grp, obs.val_grp, var.val_grp, n.val_grp, CI.bound = NA)
      
      if (error.bar){ data_obs$CI.bound <- 1.96*sqrt(data_obs[,3]/data_obs[,4]) } 
      
      col <- c("Expected" = "black", "Observed" = "red")
      if (!is.null(color)) col <- c("Expected" = color$expected[1], "Observed" = color$observed[1])
      
      p <- ggplot(data = data_exp, aes(x = Tot.val, y= exp.val, color = "Expected")) + 
        geom_line(linetype = "dashed") + geom_point() +
        geom_point(data = data_obs, aes(x = Tot.val_grp, y = obs.val_grp, color = "Observed"), size = 2) +
        theme(legend.title = element_blank(), plot.title = element_text(size = 15,hjust = 0.5)) + scale_colour_manual(values = col) +
        ggtitle(paste0("Item ", which.item)) + xlab("Total Score") + ylab("Conditional Item-Score")  + 
        scale_x_continuous(breaks=Tot.val) + 
        geom_errorbar(data = data_obs, 
                      aes(x = Tot.val_grp, y = obs.val_grp, ymin = obs.val_grp - CI.bound, ymax = obs.val_grp + CI.bound, color = "Observed"), 
                      width = 0.2, size = 1)
      
      pp[[j]] <- p 
      j <- j+1
      
    }
    
    if (grid.items){ 
      
      if (all.items){ P <- ggarrange(plotlist= pp, common.legend = T, legend = "bottom", ncol = min(2, n.items), nrow = min(2 ,ceiling(n.items/2)), align = "hv")}
      if (!all.items){P <- ggarrange(plotlist= pp, common.legend = T, legend = "bottom", ncol = 2, nrow = min(2 ,ceiling(length(which.item.arg)/2)), align = "hv")}
    }
    if (!grid.items){ P <- pp}
    
  }
  P
}
