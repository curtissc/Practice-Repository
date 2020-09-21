setwd('C:/Users/Curtiss Chapman/Desktop/MPI/fMRI project - taxonomic-thematic/Behavioral_fMRI_data/')

# grab the appropriate files

      f.0 <- list.files(getwd(), pattern='.txt$')
      f <- f.0[-grep("XXX", f.0)] 
      f <- f[grep("Tax|Them|Scram", f)] 
      
# loop through the files and put them together
      
      x.list <- vector('list', length(f))
      for(i in 1:length(f)){
            if(i %% 10 == 0) print(i/length(f))
            inf <- readLines(paste0(getwd(), '/', f[i]), n=3)
                  inf[1] <- gsub("Participant", "", inf[1])
                  inf[1] <- gsub(":| ", "", inf[1])
                  subj <- unlist(strsplit(inf[1], "_"))[1]
                  date.time <- unlist(strsplit(inf[2], " "))
                  inf[3] <- gsub("Version: ", "", inf[3])
                        CB <- unlist(strsplit(inf[3], "_"))[2] #counterbalance

            a <- read.table(paste0(getwd(), '/', f[i]), skip=3, header=T)
            b <- cbind.data.frame(file=f[i], S_ID=subj, date=date.time[1], time=date.time[2],
                                  CB=CB, a)
            x.list[[i]] <- b
      }
      x <- do.call(rbind, x.list)
            x$Task <- as.character(x$Task)
            x[x$Task=="Tax/Them", 'Task'] <- "Thematic"
      
# do all subjects have all trials?
      table(x$S_ID) # should all be 384
      
# have a loop go through each subject to re-assert their task order 
# and to check that their "CB" counterbalance variable is always the same
      
      # add task order for each subject based on time
          Ss <- unique(x$S_ID)
          
          x2.list <- vector('list', length(Ss))
          for(i in 1:length(Ss)){
               print(Ss[i])
               a <- x[x$S_ID==Ss[i],]
                     a$time <- as.POSIXct(a$time, format="%H:%M:%S")
               b <- a[!duplicated(a$time),]
                     b <- b[order(b$time),]
               a$task_order <- NA
               a$task_order <- ifelse(b[1,'Task']=='Taxonomic', "Tx_S_Tm", ifelse(b[1,'Task']=='Thematic', "Tm_S_Tx", "SOMETHING_WRONG"))
               a$task_order <- paste0(a$task_order, '_', a$CB)
               
               x2.list[[i]] <- a
          }
          x2 <- do.call(rbind, x2.list)
          
          x2$Condition <- as.character(x2$Condition)
          x2[x2$Condition=="IncongRel",'Condition'] <- 'Related'
          x2[x2$Condition=="IncongUnrel",'Condition'] <- 'Unrelated'
          x2$Condition <- factor(x2$Condition, levels=c('Mirror', 'Rotated', 'Congruent', 'Unrelated', 'Related'))
          
      # make sure everyone has a task 1, 2, and 3
          
          table(x2$task_order)/384
          
          write.csv(x2, "C:/Users/Curtiss Chapman/Desktop/MPI/fMRI project - taxonomic-thematic/Behavioral_fmri_results/fMRI taxonomic thematic behavioral participant data (long).csv")
          