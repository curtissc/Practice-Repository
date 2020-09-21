setwd('C:/Users/Curtiss Chapman/Desktop/MPI/fMRI project - taxonomic-thematic/Behavioral_fmri_results/')

x <- read.csv('fMRI taxonomic thematic behavioral participant data (long).csv', header=T, row.names=1)
S.sum.df <- read.csv("C:/Users/Curtiss Chapman/Desktop/MPI/fMRI project - taxonomic-thematic/Behavioral_fmri_results/TaxThem Behavioral fMRI subject-level results (outliers removed).csv", header=T, row.names=1)
sum <- read.csv("C:/Users/Curtiss Chapman/Desktop/MPI/fMRI project - taxonomic-thematic/Behavioral_fmri_results/TaxThem Behavioral fMRI results across subjects (outliers removed).csv", header=T, row.names=1)
load("TaxThem Behavioral fMRI result plots.Rdata")

library(ggplot2)
library(gridExtra)
library(psycho)

# plot data

      # err.plot
          grid.arrange(scram.p1, tax.p1, them.p1, ncol=3)
      # RT.plot
          grid.arrange(scram.p2, tax.p2, them.p2, ncol=3)
      # RT.out.plot
          grid.arrange(scram.p3, tax.p3, them.p3, ncol=3)
      
                  
# statistical analysis
          
      # make sure data are structured correctly
          
          scram <- S.sum.df[S.sum.df$task=="Scrambled",]
          tax <- S.sum.df[S.sum.df$task=="Taxonomic",]
          them <- S.sum.df[S.sum.df$task=="Thematic",]
          
          
      # ANOVA analysis
          
          library(car)
          
          err.lm <- lm(cbind(tax[tax$cond=="Congruent" & !is.na(tax$cond.err), 'cond.err'], 
                             tax[tax$cond=="Related" & !is.na(tax$cond.err), 'cond.err'], 
                             tax[tax$cond=="Unrelated" & !is.na(tax$cond.err), 'cond.err'], 
                             them[them$cond=="Congruent" & !is.na(them$cond.err), 'cond.err'], 
                             them[them$cond=="Related" & !is.na(them$cond.err), 'cond.err'], 
                             them[them$cond=="Unrelated" & !is.na(them$cond.err), 'cond.err'])~1)
          RT.lm <- lm(cbind(tax[tax$cond=="Congruent" & !is.na(tax$RT.m.out), 'RT.m.out'], 
                             tax[tax$cond=="Related" & !is.na(tax$RT.m.out), 'RT.m.out'], 
                             tax[tax$cond=="Unrelated" & !is.na(tax$RT.m.out), 'RT.m.out'], 
                             them[them$cond=="Congruent" & !is.na(them$RT.m.out), 'RT.m.out'], 
                             them[them$cond=="Related" & !is.na(them$RT.m.out), 'RT.m.out'], 
                             them[them$cond=="Unrelated" & !is.na(them$RT.m.out), 'RT.m.out'])~1)

          Task <- factor(c(rep('Tax', 3), rep('Them', 3)))
          Cond <- factor(rep(c('Cong', 'Rel', 'Unrel'), 2))
          
          err.anv <- Anova(err.lm, idata=data.frame(Task, Cond), idesign=~Task*Cond)
          summary(err.anv, multivariate=F)
                # Task        F(1,16)=5.96, p=.03
                # Cond        F(2,32)=8.99, p=.002 (GG corrected)
                # Task x Cond F(2,32)=1.84, p=.18 (GG corrected)
          
          RT.anv <- Anova(RT.lm, idata=data.frame(Task, Cond), idesign=~Task*Cond)
          summary(RT.anv, multivariate=F)
                # Task        F(1,16)=5.45, p=.03
                # Cond        F(2,32)=1.44, p=.25 (GG corrected)
                # Task x Cond F(2,32)=6.57, p=.005 (GG corrected)
          
          # err.plot
                grid.arrange(scram.p1, tax.p1, them.p1, ncol=3)
          # RT.out.plot
                grid.arrange(scram.p2, tax.p2, them.p2, ncol=3)
  
                
          # error contrasts
                # c('TaxCong', 'TaxRel', 'TaxUnrel', 'ThemCong', 'ThemRel', 'ThemUnrel')
                
                t.test((    -1*tax[tax$cond=="Congruent" & !is.na(tax$cond.err), 'cond.err']) +
                         (  -1*tax[tax$cond=="Related" & !is.na(tax$cond.err), 'cond.err']) + 
                         (  -1*tax[tax$cond=="Unrelated" & !is.na(tax$cond.err), 'cond.err']) + 
                         (  1*them[them$cond=="Congruent" & !is.na(them$cond.err), 'cond.err']) + 
                         (  1*them[them$cond=="Related" & !is.na(them$cond.err), 'cond.err']) +
                         (  1*them[them$cond=="Unrelated" & !is.na(them$cond.err), 'cond.err']))
                # tax has fewer errors than them: t(16)=2.44, p=.03
                
                t.test((  1*tax[tax$cond=="Congruent" & !is.na(tax$cond.err), 'cond.err']) +
                      (   -1*tax[tax$cond=="Related" & !is.na(tax$cond.err), 'cond.err']) + 
                      (   0*tax[tax$cond=="Unrelated" & !is.na(tax$cond.err), 'cond.err']) + 
                      (   1*them[them$cond=="Congruent" & !is.na(them$cond.err), 'cond.err']) + 
                      (   -1*them[them$cond=="Related" & !is.na(them$cond.err), 'cond.err']) +
                        ( 0*them[them$cond=="Unrelated" & !is.na(them$cond.err), 'cond.err']))
                # congruent vs. related across tasks: t(16)=2.66, p=.02
                
                t.test((  0*tax[tax$cond=="Congruent" & !is.na(tax$cond.err), 'cond.err']) +
                         (   1*tax[tax$cond=="Related" & !is.na(tax$cond.err), 'cond.err']) + 
                         (   -1*tax[tax$cond=="Unrelated" & !is.na(tax$cond.err), 'cond.err']) + 
                         (   0*them[them$cond=="Congruent" & !is.na(them$cond.err), 'cond.err']) + 
                         (   1*them[them$cond=="Related" & !is.na(them$cond.err), 'cond.err']) +
                         ( -1*them[them$cond=="Unrelated" & !is.na(them$cond.err), 'cond.err']))
                # related vs. unrelated across tasks: t(16)=1.91, p=.07
                
          # RT contrasts
                # c('TaxCong', 'TaxRel', 'TaxUnrel', 'ThemCong', 'ThemRel', 'ThemUnrel')
                
                t.test((    -1*tax[tax$cond=="Congruent" & !is.na(tax$RT.m.out), 'RT.m.out']) +
                         (  -1*tax[tax$cond=="Related" & !is.na(tax$RT.m.out), 'RT.m.out']) + 
                         (  -1*tax[tax$cond=="Unrelated" & !is.na(tax$RT.m.out), 'RT.m.out']) + 
                         (  1*them[them$cond=="Congruent" & !is.na(them$RT.m.out), 'RT.m.out']) + 
                         (  1*them[them$cond=="Related" & !is.na(them$RT.m.out), 'RT.m.out']) +
                         (  1*them[them$cond=="Unrelated" & !is.na(them$RT.m.out), 'RT.m.out']))
                # tax is performed more quickly than them: t(16)=2.33, p=.03
                
                t.test((    0*tax[tax$cond=="Congruent" & !is.na(tax$RT.m.out), 'RT.m.out']) +
                         (  1*tax[tax$cond=="Related" & !is.na(tax$RT.m.out), 'RT.m.out']) + 
                         (  -1*tax[tax$cond=="Unrelated" & !is.na(tax$RT.m.out), 'RT.m.out']) + 
                         (  0*them[them$cond=="Congruent" & !is.na(them$RT.m.out), 'RT.m.out']) + 
                         (  -1*them[them$cond=="Related" & !is.na(them$RT.m.out), 'RT.m.out']) +
                         (  1*them[them$cond=="Unrelated" & !is.na(them$RT.m.out), 'RT.m.out']))
                # unrelated vs. related interaction is different depending on task: t(16)=3.22, p=.005
                    # in taxonomic task, unrelated is faster than related
                    # in thematic task, related is faster than unrelated
