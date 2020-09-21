setwd('C:/Users/Curtiss Chapman/Desktop/MPI/fMRI project - taxonomic-thematic/Behavioral_fmri_results/')

x <- read.csv('fMRI taxonomic thematic behavioral participant data (long).csv', header=T, row.names=1)

library(ggplot2)
library(gridExtra)
      library(psycho)

# check subject counterbalancing

      cbind.data.frame(table(x$task_order)/384) # divided by number of total trials to get # of subjects per condition

# get summary data by condition for each participant

      Ss <- as.character(unique(x$S_ID))
      
      S.sum.list <- vector('list', length(Ss)*8)
      ii=1
      for(i in 1:length(Ss)){
            print(Ss[i])
            a <- x[x$S_ID==Ss[i], ]
            Ts <- unique(a$Task)
            for(j in 1:length(Ts)){
                  b <- a[a$Task==Ts[j], ]
                  Cs <- unique(b$Condition)
                  
                  # d prime calculation
                  #     NOTE: addition of 0.5 to each cell (which also adds 1 to the denominator of the proportion)
                  #     is the "log-linear" correction for extreme hit/false alarm rates of 0 and 1 (Hautus, 1995)
                        if(b$Task[1]=="Scrambled"){
                          hit <- nrow(b[b$Accuracy==1 & !is.na(b$Accuracy) & b$Condition=="Mirror",])
                          false.alarm <- nrow(b[b$Accuracy==0 & !is.na(b$Accuracy) & b$Condition!="Mirror",])
                          miss <- nrow(b[b$Accuracy==0 & !is.na(b$Accuracy) & b$Condition=="Mirror",])
                          corr.rej <- nrow(b[b$Accuracy==1 & !is.na(b$Accuracy) & b$Condition!="Mirror",])
                        } else {
                          hit <- nrow(b[b$Accuracy==1 & !is.na(b$Accuracy) & b$Condition=="Congruent",])
                          false.alarm <- nrow(b[b$Accuracy==0 & !is.na(b$Accuracy) & b$Condition!="Congruent",])
                          miss <- nrow(b[b$Accuracy==0 & !is.na(b$Accuracy) & b$Condition=="Congruent",])
                          corr.rej <- nrow(b[b$Accuracy==1 & !is.na(b$Accuracy) & b$Condition!="Congruent",])
                        }
                  
                        d.df <- data.frame(hit, false.alarm, miss, corr.rej)
                        ds <- dprime(d.df$hit, d.df$false.alarm, d.df$miss, d.df$corr.rej, adjusted=T)
                        
                        hit.p <- hit/(hit+miss)
                        fa.p <- false.alarm/(false.alarm+corr.rej)

                        d.prime <- qnorm(hit.p) - qnorm(fa.p)
                  
                  for(k in 1:length(Cs)){
                        c <- b[b$Condition==Cs[k], ]
                    
                        # c$Response <- ifelse(c$Response==1, 2, 1)
                        # c$Accuracy <- ifelse(c$Expected_response==c$Response, 1, 0)
                        # 
                        # count condition errors
                              err <- 1-mean(c$Accuracy, na.rm=T)
                              err.n <- nrow(c[!is.na(c$Accuracy),])
                        # get mean RT and mean RT w/outliers
                              c[c$Accuracy==0, 'RT'] <- NA
                              c$outlier.3SD <- 0
                              c[c$Accuracy==0, 'outlier.3SD'] <- NA
                              
                              RT.m <- mean(c$RT, na.rm=T)
                              RT.sd <- sd(c$RT, na.rm=T)
                              RT.n <- nrow(c[c$Accuracy==1 & !is.na(c$Accuracy),])
                              
                              if(!is.na(RT.m) & !is.na(RT.sd)){
                                tmp <-  c[(c$RT<(RT.m-(3*RT.sd)) | c$RT>(RT.m+(3*RT.sd))) & !is.na(c$RT), ]
                              } else {tmp <- NULL}
                              
                              if(!is.null(nrow(tmp))){
                                    c[(c$RT<(RT.m-(3*RT.sd)) | c$RT>(RT.m+(3*RT.sd))) & !is.na(c$RT), 'outlier.3SD'] <- 1
                              }
                              
                              RT.m.out <- mean(c[c$outlier.3SD==0 & !is.na(c$outlier.3SD), 'RT'], na.rm=T)
                              RT.n.out <- nrow(c[c$outlier.3SD==0 & !is.na(c$outlier.3SD),])
                              out <- nrow(c[c$outlier.3SD==1 & !is.na(c$outlier.3SD),])
                              
                        # put condition info together
                              
                              S.sum.list[[ii]] <- cbind.data.frame(S_ID=Ss[i], task=Ts[j], cond=Cs[k],
                                               d.prime=ds[1], beta=ds[2],
                                               task.err=1-mean(b$Accuracy, na.rm=T), task.RT.m=mean(b[b$Accuracy==1 & !is.na(b$Accuracy), 'RT'], na.rm=T),
                                               cond.err=err, err.n=err.n,
                                               RT.m=RT.m, RT.n=RT.n,
                                               RT.m.out=RT.m.out, RT.n.out=RT.n.out, out=out)
                              ii = ii+1
                  }
                  
            }  
      }
      S.sum.df <- do.call(rbind, S.sum.list)
      
      S.sum.df <- S.sum.df[order(S.sum.df$S_ID, S.sum.df$task, S.sum.df$cond),]
      
# subjects with extreme error amounts or low d'
      # I think exclusion criteria should be:
      #     S has overall error proportion > 40%
      #     S has single condition error proportion > 50% DEFINITELY
      #     S has multiple condition error proportions > 30%
      #     S has d' in task of < .5
      
      S.sum.df[S.sum.df$task.err>.4,] # high overall error
          # YA001
          # YA002
      S.sum.df[S.sum.df$cond.err>.5,c(1:6,9)] # high error in any condition
          # YA001: Them IncongRel
          # YA002: Tax IncongRel
      high.cond.errs <- S.sum.df[S.sum.df$cond.err>.3,c(1:7,10)] # high error in any condition
      high.cond.errs[high.cond.errs$S_ID==unique(high.cond.errs$S_ID)[table(high.cond.errs$S_ID)>1],]
          # YA002: Scrambled Mirror; Tax Cong; Tax IncongRel; Tax IncongUnrel
      S.sum.df[S.sum.df$dprime<1,] # low d' in any task
          # Ss < .5
          #       YA002: d'=.21
      
      # SUBJECTS TO EXCLUDE:
      #     YA001
      #     YA002
      #     YA016
      
      # save df with all participants
            S.sum.df.0 <- S.sum.df 
      
      # exclude subjects
            S.sum.df[S.sum.df$S_ID=="TTY001" | S.sum.df$S_ID=="TTY002" | S.sum.df$S_ID=="TTY016", 
                     c('d.prime','task.err','task.RT.m','cond.err','err.n',
                      'RT.m','RT.n','RT.m.out','RT.n.out','out')] <- NA
            x[x$S_ID=='TTY001' | x$S_ID=="TTY002" | x$S_ID=="TTY016" & !is.na(x$S_ID), 16:18] <- NA
      
# check counterbalancing after exclusion
      cbind.data.frame(table(x$task_order)/384) # divided by number of total trials to get # of subjects per condition
      
      
# summarize across participants

      # create data for within-subject error bars
      
            WS.list <- vector('list', nrow(S.sum.df)/8)
            for(i in 1:(nrow(S.sum.df)/8)){
                  a <- S.sum.df[((i*8)-7):(i*8),]
                  sc.0 <- a[1,]
                  ta.0 <- a[3,]
                  th.0 <- a[6,]
                  sc.0[,c('cond.err', 'RT.m', 'RT.m.out')] <- colMeans(a[1:2,c('cond.err', 'RT.m', 'RT.m.out')], na.rm=T)
                  ta.0[,c('cond.err', 'RT.m', 'RT.m.out')] <- colMeans(a[3:5,c('cond.err', 'RT.m', 'RT.m.out')], na.rm=T)
                  th.0[,c('cond.err', 'RT.m', 'RT.m.out')] <- colMeans(a[6:8,c('cond.err', 'RT.m', 'RT.m.out')], na.rm=T)
                  WS.list[[i]] <- rbind(sc.0, ta.0, th.0)
            }
            WS <- do.call(rbind, WS.list)
            
            T.ms <- lapply(c(1:3), function(f) {
                  a <- colMeans(WS[WS$task==Ts[f],c('cond.err', 'RT.m', 'RT.m.out')], na.rm=T)
                  return(a)
            })
            WS$cond.err.m <- NA
            WS$RT.m.m <- NA
            WS$RT.m.out.m <- NA
            WS[WS$task==Ts[1], c('cond.err.m')] <- T.ms[[1]][1]
            WS[WS$task==Ts[1], c('RT.m.m')] <- T.ms[[1]][2]
            WS[WS$task==Ts[1], c('RT.m.out.m')] <- T.ms[[1]][3]
            WS[WS$task==Ts[2], c('cond.err.m')] <- T.ms[[2]][1]
            WS[WS$task==Ts[2], c('RT.m.m')] <- T.ms[[2]][2]
            WS[WS$task==Ts[2], c('RT.m.out.m')] <- T.ms[[2]][3]
            WS[WS$task==Ts[3], c('cond.err.m')] <- T.ms[[3]][1]
            WS[WS$task==Ts[3], c('RT.m.m')] <- T.ms[[3]][2]
            WS[WS$task==Ts[3], c('RT.m.out.m')] <- T.ms[[3]][3]
            
            WS.se.list <- vector('list', nrow(S.sum.df)/8)
            for(i in 1:(nrow(S.sum.df)/8)){
              a <- S.sum.df[((i*8)-7):(i*8),]
              b <- WS[WS$S_ID==a$S_ID[1],]
              sc.0 <- a[1:2,]
              ta.0 <- a[3:5,]
              th.0 <- a[6:8,]
              sc.0[1,c('cond.err', 'RT.m', 'RT.m.out')] <- sc.0[1,c('cond.err', 'RT.m', 'RT.m.out')] - b[1,c('cond.err', 'RT.m', 'RT.m.out')]
              sc.0[2,c('cond.err', 'RT.m', 'RT.m.out')] <- sc.0[2,c('cond.err', 'RT.m', 'RT.m.out')] - b[1,c('cond.err', 'RT.m', 'RT.m.out')]
              ta.0[1,c('cond.err', 'RT.m', 'RT.m.out')] <- ta.0[1,c('cond.err', 'RT.m', 'RT.m.out')] - b[2,c('cond.err', 'RT.m', 'RT.m.out')]
              ta.0[2,c('cond.err', 'RT.m', 'RT.m.out')] <- ta.0[2,c('cond.err', 'RT.m', 'RT.m.out')] - b[2,c('cond.err', 'RT.m', 'RT.m.out')]
              ta.0[3,c('cond.err', 'RT.m', 'RT.m.out')] <- ta.0[3,c('cond.err', 'RT.m', 'RT.m.out')] - b[2,c('cond.err', 'RT.m', 'RT.m.out')]
              th.0[1,c('cond.err', 'RT.m', 'RT.m.out')] <- th.0[1,c('cond.err', 'RT.m', 'RT.m.out')] - b[3,c('cond.err', 'RT.m', 'RT.m.out')]
              th.0[2,c('cond.err', 'RT.m', 'RT.m.out')] <- th.0[2,c('cond.err', 'RT.m', 'RT.m.out')] - b[3,c('cond.err', 'RT.m', 'RT.m.out')]
              th.0[3,c('cond.err', 'RT.m', 'RT.m.out')] <- th.0[3,c('cond.err', 'RT.m', 'RT.m.out')] - b[3,c('cond.err', 'RT.m', 'RT.m.out')]
              WS.se.list[[i]] <- rbind(sc.0, ta.0, th.0)
            }
            WS.se <- do.call(rbind, WS.se.list)
            
      
      se <- function(x) sqrt(var(x)/length(x))
      
      d.sum <- aggregate(dprime~task*cond, S.sum.df, mean)
      beta.sum <- aggregate(beta~task*cond, S.sum.df, mean)
      err.sum <- aggregate(cond.err~task*cond, S.sum.df, mean)
      err.se <- aggregate(cond.err~task*cond, WS.se, se)
      RT.sum <- aggregate(RT.m~task*cond, S.sum.df, mean)
      RT.se <- aggregate(RT.m~task*cond, WS.se, se)
      RT.out.sum <- aggregate(RT.m.out~task*cond, S.sum.df, mean)
      RT.out.se <- aggregate(RT.m.out~task*cond, WS.se, se)
      RT.n.out.sum <- aggregate(RT.n.out~task*cond, S.sum.df, FUN=function(f) sum(f))
      out.sum <- aggregate(out~task*cond, S.sum.df, FUN=function(f) sum(f))
      err.n.sum <- aggregate(err.n~task*cond, S.sum.df, FUN=function(f) sum(f))
      RT.n.sum <- aggregate(RT.n~task*cond, S.sum.df, FUN=function(f) sum(f))
      
      sum <- cbind.data.frame(err.sum[,1:2], n=length(unique(S.sum.df[!is.na(S.sum.df$RT.m), 'S_ID'])), d.sum[3], beta.sum[3],
                              err.obs=err.n.sum[3][,1], err.sum[3], err.se=err.se[3][,1], 
                              RT.obs=RT.n.sum[3][,1], RT.sum[3], RT.se=RT.se[3][,1], 
                              RT.out.obs=RT.n.out.sum[3][,1], RT.out.sum[3], RT.out.se=RT.out.se[3][,1],
                              out.sum[3], out.p=round(out.sum[3][,1]/(out.sum[3][,1]+RT.n.out.sum[3]), 3)[,1])
            sum <- sum[order(sum$task, sum$cond),]
      
            sum.scram <- sum[sum$task=="Scrambled",]
            sum.tax <- sum[sum$task=="Taxonomic",]
            sum.them <- sum[sum$task=="Thematic",]
            
# visualize subject results
            # the outlier.shape variable in geom_boxplot 
            # keeps the outliersfrom being duplicated by the jitter command
            
            
      # subject boxplots
            
            # errors
                  S.box.err <- ggplot(S.sum.df, aes(x=cond, y=cond.err)) + 
                    geom_boxplot(alpha=.7, outlier.shape = NA) + 
                    geom_jitter(width=.2, alpha=.4) +
                    facet_grid(. ~ task)
                  
                  # including outlier subjects
                        S.box.err.0 <- ggplot(S.sum.df.0, aes(x=cond, y=cond.err)) +
                          geom_boxplot(alpha=.7, outlier.shape = NA) +
                          geom_jitter(width=.2, alpha=.4) +
                          facet_grid(. ~ task)
      
            # RTs (outliers not excluded)
                  S.box.RT <- ggplot(S.sum.df, aes(x=cond, y=RT.m)) + 
                    geom_boxplot(alpha=.7, outlier.shape = NA) + 
                    geom_jitter(width=.2, alpha=.4) +
                    facet_grid(. ~ task)
                  
                  # including outlier subjects
                        S.box.RT.0 <- ggplot(S.sum.df.0, aes(x=cond, y=RT.m)) +
                          geom_boxplot(alpha=.7, outlier.shape = NA) +
                          geom_jitter(width=.2, alpha=.4) +
                          facet_grid(. ~ task)
                        
            # RTs (outliers excluded)
                  S.box.RT.out <- ggplot(S.sum.df, aes(x=cond, y=RT.m.out)) + 
                    geom_boxplot(alpha=.7, outlier.shape = NA) + 
                    geom_jitter(width=.2, alpha=.4) +
                    facet_grid(. ~ task)
                  
                  # including outlier subjects
                        S.box.RT.out.0 <- ggplot(S.sum.df.0, aes(x=cond, y=RT.m.out)) +
                          geom_boxplot(alpha=.7, outlier.shape = NA) +
                          geom_jitter(width=.2, alpha=.4) +
                          facet_grid(. ~ task)
                        
      # subject means + SE
                  
            sum.taxthem <- sum[sum$task!="Scrambled",]
                        
            min.err <- floor(min(sum.taxthem$cond.err - sum.taxthem$err.se)*100)/100
            max.err <- ceiling(max(sum.taxthem$cond.err + sum.taxthem$err.se)*100)/100
            lim.err <- c(min.err, max.err)
            min.RT <- floor(min(sum.taxthem$RT.m - sum.taxthem$RT.se)/100)*100
            max.RT <- ceiling(max(sum.taxthem$RT.m + sum.taxthem$RT.se)/100)*100
            lim.RT <- c(min.RT, max.RT)
            min.RT.out <- floor(min(sum.taxthem$RT.m.out - sum.taxthem$RT.out.se)/100)*100
            max.RT.out <- ceiling(max(sum.taxthem$RT.m.out + sum.taxthem$RT.out.se)/100)*100
            lim.RT.out <- c(min.RT.out, max.RT.out)
            
            
            # errors
                        
                  scram.p1 <- ggplot(sum.scram, aes(x=cond, y=cond.err)) +
                    geom_point(position=position_dodge(.9)) +
                    geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin=cond.err-err.se, ymax=cond.err+err.se)) +
                    xlab("Condition") + ylab("errors") +
                    ggtitle("Scrambled")
                  tax.p1 <- ggplot(sum.tax, aes(x=cond, y=cond.err)) +
                    geom_point(position=position_dodge(.9)) +
                    geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin=cond.err-err.se, ymax=cond.err+err.se)) +
                    ylim(c(lim.err)) +
                    xlab("Condition") + ylab("errors") +
                    ggtitle("Taxonomic")
                  them.p1 <- ggplot(sum.them, aes(x=cond, y=cond.err)) +
                    geom_point(position=position_dodge(.9)) +
                    geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin=cond.err-err.se, ymax=cond.err+err.se)) +
                    ylim(c(lim.err)) +
                    xlab("Condition") + ylab("errors") +
                    ggtitle("Thematic")

            # RTs (raw, no errors)
                  scram.p2 <- ggplot(sum.scram, aes(x=cond, y=RT.m)) +
                    geom_point(position=position_dodge(.9)) +
                    geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin=RT.m-RT.se, ymax=RT.m+RT.se)) +
                    xlab("Condition") + ylab("RT") +
                    ggtitle("Scrambled")
                    tax.p2 <- ggplot(sum.tax, aes(x=cond, y=RT.m)) +
                    geom_point(position=position_dodge(.9)) +
                    geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin=RT.m-RT.se, ymax=RT.m+RT.se)) +
                      ylim(c(lim.RT)) +
                      xlab("Condition") + ylab("RT") +
                    ggtitle("Taxonomic")
                  them.p2 <- ggplot(sum.them, aes(x=cond, y=RT.m)) +
                    geom_point(position=position_dodge(.9)) +
                    geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin=RT.m-RT.se, ymax=RT.m+RT.se)) +
                    ylim(c(lim.RT)) +
                    xlab("Condition") + ylab("RT") +
                    ggtitle("Thematic")

            # RTs (outliers removed)
                  
                  scram.p3 <- ggplot(sum.scram, aes(x=cond, y=RT.m.out)) +
                    geom_point(position=position_dodge(.9)) +
                    geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin=RT.m.out-RT.out.se, ymax=RT.m.out+RT.out.se)) +
                    xlab("Condition") + ylab("RT") +
                    ggtitle("Scrambled")
                  tax.p3 <- ggplot(sum.tax, aes(x=cond, y=RT.m.out)) +
                    geom_point(position=position_dodge(.9)) +
                    geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin=RT.m.out-RT.out.se, ymax=RT.m.out+RT.out.se)) +
                    ylim(c(lim.RT.out)) +
                    xlab("Condition") + ylab("RT") +
                    ggtitle("Taxonomic")
                  them.p3 <- ggplot(sum.them, aes(x=cond, y=RT.m.out)) +
                    geom_point(position=position_dodge(.9)) +
                    geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin=RT.m.out-RT.out.se, ymax=RT.m.out+RT.out.se)) +
                    ylim(c(lim.RT.out)) +
                    xlab("Condition") + ylab("RT") +
                    ggtitle("Thematic")
         
# get summary data by condition for each item
        
        x$item12 <- paste0(x$Pic1, x$Pic2)
        Is <- as.character(unique(x$item12))
        Ts <- unique(x$Task)
        
        I.sum.list <- vector('list', length(Is)*2)
        ii=1
        for(i in 1:length(Is)){
          if(i %% 25 == 0) print(i/length(Is))
          for(j in 1:length(Ts)){
            a <- x[x$item12==Is[i] & x$Task==Ts[j], ]
            if(nrow(a)==0) next
            
            # errors
                  err <- 1-mean(a$Accuracy, na.rm=T)
                  err.n <- nrow(a[!is.na(a$Accuracy),])
            # RTs
                  a[a$Accuracy==0 & !is.na(a$Accuracy), 'RT'] <- NA
                  a$outlier.3SD <- ifelse(!is.na(a$RT), 0, NA)
                  a[a$Accuracy==0 & !is.na(a$Accuracy), 'outlier.3SD'] <- NA
                  RT.m <- mean(a$RT, na.rm=T)
                  RT.sd <- sd(a$RT, na.rm=T)
                  RT.n <- nrow(a[a$Accuracy==1 & !is.na(a$Accuracy),])
                  a[(a$RT<(RT.m-(3*RT.sd)) | a$RT>(RT.m+(3*RT.sd))) & !is.na(a$RT), 'outlier.3SD'] <- 1
                  RT.m.out <- mean(a[a$outlier.3SD==0 & !is.na(a$outlier.3SD), 'RT'], na.rm=T)
                  RT.n.out <- nrow(a[a$outlier.3SD==0 & !is.na(a$outlier.3SD),])
                  out <- nrow(a[a$outlier.3SD==1 & !is.na(a$outlier.3SD),])
                  
            # put the info together
                  I.sum.list[[ii]] <- cbind.data.frame(item1=a$Pic1[1], item2=a$Pic2[1], item12=Is[i], task=Ts[j], cond=a$Condition[1],
                                                       err=err, err.n=err.n,
                                                       RT.m=RT.m, RT.n=RT.n,
                                                       RT.m.out=RT.m.out, RT.n.out=RT.n.out, out=out)
                  ii = ii+1
          }
        }
        I.sum.df <- do.call(rbind, I.sum.list)
        
        I.sum.df <- I.sum.df[order(I.sum.df$err, I.sum.df$RT.m, decreasing = T),]
        
# look for items that were problematic
        
        I.sum.df[I.sum.df$err>.49,]
        I.sum.df[I.sum.df$out>(.25*I.sum.df$err.n),]
        I.sum.df[I.sum.df$out>0,]

# exclude problematic items
        
        
# summarize effects from tasks
        
      scram.S <- S.sum.df.0[S.sum.df.0$task=="Scrambled",]
      tax.S <- S.sum.df.0[S.sum.df.0$task=="Taxonomic",]
      them.S <- S.sum.df.0[S.sum.df.0$task=="Thematic",]
      
      S.fx.list <- vector('list', nrow(scram.S)/2)
      ii=1
      for(i in 1:(nrow(scram.S)/2)){
          sc <- scram.S[((i*2)-1):((i*2)), ]
          ta <- tax.S[((i*3)-2):((i*3)), ]
          th <- them.S[((i*3)-2):((i*3)), ]
          
          sc.0 <- sc[1,]
          ta.0 <- ta[1,]
          th.0 <- th[1,]
          
          sc.0[,c('cond.err', 'RT.m', 'RT.m.out')] <- sc[1,c('cond.err', 'RT.m', 'RT.m.out')] - 
                                                                          sc[2,c('cond.err', 'RT.m', 'RT.m.out')]
          sc.0$cond <- 'Mirr-Rot'
          sc.0$RT.n <- paste0(sc[1,'RT.n'], "_", sc[2,'RT.n'])
          sc.0$RT.n.out <- paste0(sc[1,'RT.n.out'], "_", sc[2,'RT.n.out'])
          sc.0$out <- paste0(sc[1,'out'], "_", sc[2,'out'])
          
          ta.0[,c('cond.err', 'RT.m', 'RT.m.out')] <- ta[3,c('cond.err', 'RT.m', 'RT.m.out')] - 
                                                                          ta[2,c('cond.err', 'RT.m', 'RT.m.out')]
          ta.0$cond <- 'Rel-Unrel'
          ta.0$RT.n <- paste0(ta[3,'RT.n'], "_", ta[2,'RT.n'])
          ta.0$RT.n.out <- paste0(ta[3,'RT.n.out'], "_", ta[2,'RT.n.out'])
          ta.0$out <- paste0(ta[3,'out'], "_", ta[2,'out'])
          
          th.0[,c('cond.err', 'RT.m', 'RT.m.out')] <- th[3,c('cond.err', 'RT.m', 'RT.m.out')] - 
                                                                          th[2,c('cond.err', 'RT.m', 'RT.m.out')]
          th.0$cond <- 'Rel-Unrel'
          th.0$RT.n <- paste0(th[3,'RT.n'], "_", th[2,'RT.n'])
          th.0$RT.n.out <- paste0(th[3,'RT.n.out'], "_", th[2,'RT.n.out'])
          th.0$out <- paste0(th[3,'out'], "_", th[2,'out'])
          
          S.fx.list[[i]] <- rbind(sc.0, ta.0, th.0)
      }
      S.fx <- do.call(rbind, S.fx.list)
        
      
# match up behavioral data & survey data to participants
      
      
# write dfs 
      
      write.csv(S.sum.df.0, "C:/Users/Curtiss Chapman/Desktop/MPI/fMRI project - taxonomic-thematic/Behavioral_fmri_results/TaxThem Behavioral fMRI subject-level results (all subjects).csv")
      write.csv(S.sum.df, "C:/Users/Curtiss Chapman/Desktop/MPI/fMRI project - taxonomic-thematic/Behavioral_fmri_results/TaxThem Behavioral fMRI subject-level results (outliers removed).csv")
      write.csv(sum, "C:/Users/Curtiss Chapman/Desktop/MPI/fMRI project - taxonomic-thematic/Behavioral_fmri_results/TaxThem Behavioral fMRI results across subjects (outliers removed).csv")
      save(S.box.err, S.box.RT, S.box.RT.out,
           scram.p1, scram.p2, scram.p3,
           tax.p1, tax.p2, tax.p3,
           them.p1, them.p2, them.p3, file="TaxThem Behavioral fMRI result plots.Rdata")

      # means       
      # err.plot
      grid.arrange(scram.p1, tax.p1, them.p1, ncol=3)
      # RT.plot
      grid.arrange(scram.p2, tax.p2, them.p2, ncol=3)
      # RT.out.plot
      grid.arrange(scram.p3, tax.p3, them.p3, ncol=3)
      
      
#############################################
###  CODE TAKEAWAYS  ########################
#############################################        
### output dfs
    S.sum.df.0 # task x condition results for all subjects
    S.sum.df # task x condition results excluding subjects with high error rates
    S.fx
        S.fx[S.fx$task=="Scrambled",]
        S.fx[S.fx$task=="Taxonomic",]
        S.fx[S.fx$task=="Thematic",]
    sum     # task x condition results averaged across subjects
    
### output plots
    # boxplots
          S.box.err
          S.box.err.0
          S.box.RT
          S.box.RT.0
          S.box.RT.out
          S.box.RT.out.0
    
    # means       
          # err.plot
              grid.arrange(scram.p1, tax.p1, them.p1, ncol=3)
          # RT.plot
              grid.arrange(scram.p2, tax.p2, them.p2, ncol=3)
          # RT.out.plot
              grid.arrange(scram.p3, tax.p3, them.p3, ncol=3)
    
    