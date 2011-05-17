# validated R version
validatedVersion <- "2.12.2"

# identity of RS and QC sera for the experiment
ref <- c("AVR801", "AVR801F2")
qc <- c("NR717")

# text used to indicate unused wells on a plate
unused <- c("blank", "empty")

# numerical plate acceptance criteria (prefix cr.)
cr.top <- c(0.861, 2.35)
cr.bot <- c(0.209, 0.510)
cr.fm.ref.deviance <- 0.117
cr.rED50 <- c(461, 1010)
cr.fm.qc.deviance <- 0.272
cr.RNA.qc <- c(0.0627, 0.195)

# sample acceptance criteria
cr.smp.toxicity <- 1.4
cr.smp.similarity <- c(1, 1) # placeholder - not currently implemented
cr.smp.deviance <- 1 # placeholder - not currently implemented
cr.smp.uncertainty <- 0.132 # 99% quantile for samples in the validation expt. with
                            # RNA inside [-1, 1] n.b. closed interval
