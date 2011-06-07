# implementation of RNA data reduction for the TNA
# pre-alpha version written 2010 by
# K. Omland, Mergus Analytics, LLC (www.mergusanalytics.com)
#
# additional pre-alpha modifications by R Carnell, Battelle
# April 20, 2011
#
# and by K. Omland
# May 10, 2011

# RNAWrapper performs data reduction for a single plate
# accepting inputObject and producing outputObject according to
#    RNA Algorithm Input_Output Specification rev1 20110421.pdf

RNAWrapper <- function(inputObject, ref, qc) {
    outputObject <- list(plate.info = NA,
                         fm.ref = NA,
                         fm.qc = NA,
                         sample.results = NA)

    outputObject$plate.info <- data.frame(plate = inputObject$plate,
                                          operator = inputObject$operator,
                                          run = inputObject$run,
                                          lab = inputObject$lab)

    # dat.pl is the raw data from the plate
    # unused wells on the plate are removed in the subsetting clause
    ind <- which(!is.na(inputObject$dat$smp) & !tolower(inputObject$dat$smp) %in% unused)
    dat.pl <- inputObject$dat[ind,]

    # assign dummy variables
    # ref and qc may be vectors of length 1 or more (e.g., AVR801 and AVR801F2)
    # anything else treated as tst
    dat.pl$ref <- ifelse(!is.na(dat.pl$smp) & dat.pl$smp %in% ref,
                         1, 0)
    dat.pl$qc  <- ifelse(!is.na(dat.pl$smp) & dat.pl$smp %in% qc,
                         1, 0)
    dat.pl$tst <- ifelse(!is.na(dat.pl$smp) & !(dat.pl$smp %in% c(ref, qc)) &
                         !substr(dat.pl$type, 1, 2) %in% c("NC", "SC"),
                         1, 0)

    # independent and response variables
    dat.pl$log.dln <- log(dat.pl$dln)
    for(i in 1:length(dat.pl$OD)){
        if(dat.pl$OD[i] > 0) {
            dat.pl$log.OD[i] <-  log(dat.pl$OD[i])
        } else {
            warning("all ODs must be positive")
        }
    }

    # check plate acceptance criteria on raw ODs
    outputObject$plate.info$od.accept <- check.ods(dat.pl)

    if (outputObject$plate.info$od.accept) {
        # check plate acceptance criteria on models fit to ref and qc
        ref.qc.info <- check.ref.qc(dat.pl)

        outputObject$fm.ref                       <- ref.qc.info$fm.ref
        outputObject$plate.info$fm.ref.accept.fit <- ref.qc.info$fm.ref.accept.fit
        outputObject$plate.info$fm.ref.accept.est <- ref.qc.info$fm.ref.accept.est
        outputObject$fm.qc                        <- ref.qc.info$fm.qc
        outputObject$plate.info$fm.qc.accept.fit  <- ref.qc.info$fm.qc.accept.fit
        outputObject$plate.info$fm.qc.accept.est  <- ref.qc.info$fm.qc.accept.est

        outputObject$plate.info$accept.plt <- ref.qc.info$fm.ref.accept.fit &
                                              ref.qc.info$fm.ref.accept.est &
                                              (is.na(ref.qc.info$fm.qc.accept.fit) |
                                               ref.qc.info$fm.qc.accept.fit) &
                                              (is.na(ref.qc.info$fm.qc.accept.est) |
                                               ref.qc.info$fm.qc.accept.est)

        # append reactivity threshold to plate.info
        if (!(length(ref.qc.info$fm.ref) == 1 && is.na(ref.qc.info$fm.ref))) {
            coef.fm.ref <- coef(ref.qc.info$fm.ref)
            outputObject$plate.info$react.threshold <-
                unname(coef.fm.ref["b"] +
                       coef.fm.ref["d"] / (3 + sqrt(3)))
        } else {
            outputObject$plate.info$react.threshold <- NA
        }

    } else {
        outputObject$plate.info$accept.plt <- FALSE
    } # end if/else od.accept

    # for samples on passing plates, estimate RNA
    if (outputObject$plate.info$accept.plt == TRUE) {
        smps <- as.character(sort(unique(dat.pl$smp[dat.pl$tst == 1])))
        smps <- smps[!tolower(smps) %in% c("blank", "empty")]

        for (tstsmp in smps) {
            # determine sample toxicity
            toxic <- smp.tox(dat.pl,
                             tstsmp)

            # determine reactivity
            reactive <- smp.react(dat.pl,
                                  tstsmp,
                                  outputObject$plate.info$react.threshold)

            # determine sufficient dilution (for very active samples)
            diluted.enough <- dil.nuf(dat.pl,
                                      tstsmp,
                                      outputObject$plate.info$react.threshold)

            # not currently implemented:
            # fit reference curve constraint model
            # and estimate similarity test statistics

            # fit parallel curves constraint model
            # assumes row 8 used for control treatments
            ind <- which(as.character(dat.pl$smp) %in% c(ref, tstsmp) &
                         dat.pl$type %in% paste("D", 1:7, sep=""))
            fm.prl <- fit.4pl.prl(dat.pl[ind,], outputObject$fm.ref)

            # collect output for all samples on plate
            if (tstsmp == smps[1]) {
                sample.results <- summarize.sample(dat.pl,
                                                  toxic, reactive, diluted.enough,
                                                  fm.prl, tstsmp)
            } else {
                sample.results <-
                    rbind(sample.results,
                          summarize.sample(dat.pl,
                                           toxic, reactive, diluted.enough,
                                           fm.prl, tstsmp))
            }
        } # end tstsmp in smps

        outputObject$sample.results <- sample.results
    } # end if/else plate.basic$accept.plt == TRUE

    return(outputObject)
} # end RNAWrapper()

# implements plate acceptance criteria on ODs for ref
check.ods <- function(dat.pl)
{
    ref.dln <- sort(unique(dat.pl$dln[dat.pl$ref == 1]))
    top <- range(dat.pl$OD[dat.pl$ref == 1 &
                           dat.pl$dln %in% ref.dln[1:2]])
    bot <- range(c(dat.pl$OD[dat.pl$ref == 1 &
                             dat.pl$dln == ref.dln[length(ref.dln)]],
                   dat.pl$OD[dat.pl$type == "NC"]))

    return(top[1] >= cr.top[1] &
           top[2] <= cr.top[2] &
           bot[1] >= cr.bot[1] &
           bot[2] <= cr.bot[2])
} # end check.ods()

# implements plate acceptance criteria on models fit to ref and qc
check.ref.qc <- function(dat.pl)
{
    ref.qc.info <- list(fm.ref = NA,
                    fm.ref.accept.fit = NA,
                    fm.ref.accept.est = NA,
                    fm.qc = NA,
                    fm.qc.accept.fit = NA,
                    fm.qc.accept.est = NA)

    ref.qc.info$fm.ref <- fit.4pl.ref(dat.pl[dat.pl$ref == 1,])

    if (!(length(ref.qc.info$fm.ref) == 1 && is.na(ref.qc.info$fm.ref)))
    {
        ref.qc.info$fm.ref.accept.fit <- (deviance(ref.qc.info$fm.ref) <= cr.fm.ref.deviance)
        ref.qc.info$fm.ref.accept.est <-
            (exp(unname(coef(ref.qc.info$fm.ref)["log.rED50"])) >= cr.rED50[1] &
             exp(unname(coef(ref.qc.info$fm.ref)["log.rED50"])) <= cr.rED50[2])

      if (sum(dat.pl$qc) > 0)
      {
          qc.info <- check.qc(dat.pl, ref.qc.info$fm.ref)
          ref.qc.info$fm.qc <- qc.info$fm.qc
          ref.qc.info$fm.qc.accept.fit <- qc.info$fm.qc.accept.fit
          ref.qc.info$fm.qc.accept.est <- qc.info$fm.qc.accept.est
      }
    } else {
        ref.qc.info$fm.ref.accept.fit <- FALSE
        ref.qc.info$fm.ref.accept.est <- FALSE
    }
    return(ref.qc.info)
} # end check.ref.qc

check.qc <- function(dat.pl, fm.ref)
{
    qc.info <- list(fm.qc = NA,
                   fm.qc.accept.fit = NA,
                   fm.qc.accept.est = NA)
    # isolate just ref + qc
    ind <- which(dat.pl$ref == 1 | (dat.pl$qc == 1 & !is.na(dat.pl$dln)))
    dat.foqc <- dat.pl[ind,]
    # treat qc as tst
    dat.foqc$tst <- 1 - dat.foqc$ref
    qc.info$fm.qc <- fit.4pl.prl(dat.foqc, fm.ref)

    if (!(length(qc.info) == 1 && is.na(qc.info$fm.qc)))
    {
        qc.info$fm.qc.accept.fit <- deviance(qc.info$fm.qc) <= cr.fm.qc.deviance
        qc.info$fm.qc.accept.est <-
            unname(coef(qc.info$fm.qc)["log.RNA"]) >= log(cr.RNA.qc[1]) &
        unname(coef(qc.info$fm.qc)["log.RNA"]) <= log(cr.RNA.qc[2])
    } else {
        qc.info$fm.qc.accept.fit <- FALSE
        qc.info$fm.qc.accept.est <- FALSE
    }
    return(qc.info)
} # end check.qc()

# 4PL parameterization:
#
# independent variable, x, is log(dilution)
#
# shape:
# b for baseline (lower asymptote),
#     ie, log.OD from unprotected wells / fully diluted serum
# d for curve depth, ie, upper asymptote minus lower asymptote
#     (implicit parameter a for asymptote (upper),
#      ie, log.OD from fully protected wells / undiluted)
# s for slope
#
# location:
# log.rED50 = log(ED50[ref])
# log.RNA = log(ED50[tst]/ED50[ref])
# tst: indicator variable, 0 for ref or 1 for tst (including qc)
#
# other notes:
# c for combination: c = -s * log(ED50)
#   slope constrained positive (so numerical routine cannot flip 4PL)

# four parameter logistic function as it appears in documentation
# n.b. retain single character symbols including c
fourpl.prl <- function(b, d, s, log.rED50, log.RNA, tst, x)
{
    s.pos <- abs(s)
    c <- -s.pos * (log.rED50 + tst * log.RNA)
    b + d * 1 / (1 + exp(c + s.pos * x))  # 1 / (1 + exp(c + s.pos * x))
                                          # recognizable as logistic function
} # end fourpl.prl()

# main function to fit fourpl.prl to data,
# called by fit.4pl.ref and fit.4pl.prl
fit.4pl <- function(y, tst, x,
                    b.init, d.init, s.init, log.rED50.init, log.RNA.init)
{
    if (log.RNA.init == 0) ### better to use an explicit flag for fitting fm.ref
    {
        form <- y ~ fourpl.prl(b, d, s, log.rED50, 0, tst, x)
        startlist <- list(b = b.init,
                          d = d.init,
                          s = s.init,
                          log.rED50 = log.rED50.init)
    } else {
        form <- y ~ fourpl.prl(b, d, s, log.rED50, log.RNA, tst, x)
        startlist <- list(b = b.init,
                          d = d.init,
                          s = s.init,
                          log.rED50 = log.rED50.init,
                          log.RNA = log.RNA.init)
    }

    nlsfit <-
        tryCatch(nls(form, start = startlist),
                 error = function(e)
             {
                 return(list(convInfo = list(isConv = FALSE)))
             })

    # if failed convergence, grid search for initial values
    # grid search begins near inits selected above and widens gradually
    # widens log.rED50.init first, then s.init
    if (!nlsfit$convInfo$isConv)
    {
        # n.b. first element is 0 to hold s.init at value selected above
        # while widening log.rED50.init
        for (i in 0:10)
        {
            # widening s.init up ...
            startlist$s <- s.init * exp(i / 10)

            # here first element is 1 to immediately widen log.rED50.init
            for (j1 in 1:10)
            {
                # widening log.rED50.init up ...
                startlist$log.rED50 <- log.rED50.init * exp(j1 / 10)

                nlsfit <- tryCatch(nls(form, start = startlist),
                                   error = function(e){
                                       return(list(convInfo = list(isConv = FALSE)))
                                   })

                if (nlsfit$convInfo$isConv) return(nlsfit)

                # ... widening log.rED50.init down
                startlist$log.rED50 <- log.rED50.init * exp(-j1 / 10)

                nlsfit <- tryCatch(nls(form, start = startlist),
                                   error = function(e){
                                       return(list(convInfo = list(isConv = FALSE)))
                                   })

            if (nlsfit$convInfo$isConv) return(nlsfit)
            } # end j1

            # ... widening s.init down
            startlist$s <- s.init * exp(-i / 10)

            for (j2 in 1:10)
            {
                # widening log.rED50.init up ...
                startlist$log.rED50 <- log.rED50.init * exp(j2 / 10)

                nlsfit <- tryCatch(nls(form, start = startlist),
                                   error = function(e){
                                       return(list(convInfo = list(isConv = FALSE)))
                                   })

                if (nlsfit$convInfo$isConv) return(nlsfit)

                # ... widening log.rED50.init down
                startlist$log.rED50 <- log.rED50.init * exp(-j2 / 10)

                nlsfit <- tryCatch(nls(form, start = startlist),
                                   error = function(e){
                                       return(list(convInfo = list(isConv = FALSE)))
                                   })

                if (nlsfit$convInfo$isConv) return(nlsfit)
            } # end j2

        } # end i

        # convergence not reached
        return(NA)

    } else {
        return(nlsfit)
    }
} # end function fit.4pl()

# a function to fit 4PL to ref alone
fit.4pl.ref <- function(dat.ref)
{
    # initial values
    b.init <- min(dat.ref$log.OD)
    d.init <- max(dat.ref$log.OD) - b.init
    s.init <- 3
    # interpolate log.rED50.init using
    # least dilute spike with log.OD < infl.pt and the next more dilute spike
    # or spikes at one or other end of range if ODs don't cross inflection point
    infl.pt <- b.init + d.init * 0.5
    dat.ref$dist <- dat.ref$log.OD - infl.pt
    # assumes 1:2 serial dilution
    trgt.dln <- max(dat.ref$dln[dat.ref$dist >= 0]) * c(1, 2)
    dat.interp <- dat.ref[dat.ref$dln %in% trgt.dln,]
    coef.interp <- coef(lm(dist ~ log.dln, data = dat.interp))
    log.rED50.init <- unname(-coef.interp["(Intercept)"]/coef.interp["log.dln"])

    return(fit.4pl(y = dat.ref$log.OD,
                   tst = dat.ref$tst,
                   x = dat.ref$log.dln,
                   b.init = b.init,
                   d.init = d.init,
                   s.init = s.init,
                   log.rED50.init = log.rED50.init,
                   log.RNA.init = 0))
} # end fit.4pl.ref()

# a function to fit the parallel curves constraint model to ref and a single tst
fit.4pl.prl <- function(dat.focal, fm.ref)
{
    # log.RNA.init:
    # if all ODs are less than (initial) inflection point, use first two
    # if they're all greater, use last two
    # if they're on both sides, use most dilute still above and next
    infl.pt <- unname(coef(fm.ref)["b"] + coef(fm.ref)["d"] / 2)
    dat.tst.interp <- dat.focal[dat.focal$ref == 0,]
    dat.tst.interp$dist <- dat.tst.interp$log.OD - infl.pt
    if (max(dat.tst.interp$dist) < 0)
    {
        trgt.dln <- sort(dat.tst.interp$dln)[1:2]
    } else {
        if (min(dat.tst.interp$dist) > 0)
        {
            r <- dim(dat.tst.interp)[1]
            trgt.dln <- sort(dat.tst.interp$dln)[(r-1):r]
        } else {
            # assumes 1:2 serial dilution
            trgt.dln <- max(dat.tst.interp$dln[dat.tst.interp$dist >= 0]) * c(1, 2)
        }
    }
    dat.tst.interp <- dat.tst.interp[dat.tst.interp$dln %in% trgt.dln,]
    coef.tst.interp <- coef(lm(dist ~ log.dln, data = dat.tst.interp))
    log.RNA.init <- unname(-coef.tst.interp["(Intercept)"] /
                           coef.tst.interp["log.dln"] -
                           coef(fm.ref)["log.rED50"])

    return(fit.4pl(y = dat.focal$log.OD,
                   tst = dat.focal$tst,
                   x = dat.focal$log.dln,
                   b.init = unname(coef(fm.ref)["b"]),
                   d.init = unname(coef(fm.ref)["d"]),
                   s.init = unname(coef(fm.ref)["s"]),
                   log.rED50.init = unname(coef(fm.ref)["log.rED50"]),
                   log.RNA.init = log.RNA.init)
           )
} # end function fit.4pl.prl

smp.tox <- function(dat.tox, tst)
{
    ind <- which(!is.na(dat.tox$smp) &
                 as.character(dat.tox$smp) == tst &
                 substr(dat.tox$type, 1, 2) == "SC")
    sc.od <- dat.tox$OD[ind]

    ref.dln <- sort(unique(dat.tox$dln[dat.tox$ref == 1]))
    ind <- which(dat.tox$ref == 1 & dat.tox$dln %in% ref.dln[1:2])
    mean.ref.top <- mean(dat.tox$OD[ind])

    smp.toxicity <- data.frame(test.stat = sc.od / mean.ref.top - 1)
    smp.toxicity$toxic <- smp.toxicity$test.stat >= 1/cr.smp.toxicity &
                          smp.toxicity$test.stat <= cr.smp.toxicity
    return(smp.toxicity)
} # end smp.tox()

smp.react <- function(dat.react, tst, thresh)
{
    ind <- which(dat.react$smp == tst & dat.react$type ==  "D1")
    log.OD.ts1 <- dat.react$log.OD[ind]
    return(data.frame(log.OD.ts1 = log.OD.ts1,
                      reactive = (log.OD.ts1 >= thresh)))
} # end smp.react()

dil.nuf <- function(dat.dil.nuf, tst, thresh)
{
    ind <- which(dat.dil.nuf$smp == tst & dat.dil.nuf$type == "D7")
    log.OD.ts7 <- dat.dil.nuf$log.OD[ind]
    return(data.frame(log.OD.ts7 = log.OD.ts7,
                      diluted.enough = (log.OD.ts7 <= thresh)))
} # end dil.nuf

# not currently implemented
#tst.ref.similar <- function(dat.smlr)
#{
#
#    return(data.frame(slope.resids = ,
#                      mean.resids = ))
#} # end tst.ref.similar

summarize.sample <- function(dat.pl, toxic, reactive, diluted.enough, fm.prl, tstsmp)
{
    if (!(length(fm.prl) == 1 && is.na(fm.prl)))
    {
        deviance <- deviance(fm.prl)
        log.RNA.se <- unname(sqrt(diag(vcov(fm.prl))))[5]
        accept.smp <- !toxic$toxic &
                      reactive$reactive &
                      #diluted.enough$diluted.enough &
                      #tst.sim.ref$slope.resids <= cr.smp.similarity[1] &
                      #tst.sim.ref$mean.resids <= cr.smp.similarity[2] &
                      #deviance <= cr.smp.deviance &
                      log.RNA.se <= cr.smp.uncertainty
        log.RNA <- unname(coef(fm.prl)["log.RNA"])
        log.rED50.prl <- unname(coef(fm.prl)["log.rED50"])
        b.prl <- unname(coef(fm.prl)["b"])
        d.prl <- unname(coef(fm.prl)["d"])
        s.prl <- unname(coef(fm.prl)["s"])
    } else {
        deviance <- NA
        log.RNA.se <- NA
        accept.smp <- FALSE
        log.RNA <- NA
        log.rED50.prl <- NA
        log.tED50.prl <- NA
        b.prl <- NA
        d.prl <- NA
        s.prl <- NA
    }
    return(data.frame(smp = tstsmp,
                      col = unique(na.exclude(dat.pl$col[dat.pl$smp == tstsmp])),
                      start.dln = min(na.exclude(dat.pl$dln[dat.pl$smp == tstsmp])),
                      toxic = toxic$toxic,
                      reactive = reactive$reactive,
                      diluted.enough = diluted.enough$diluted.enough,
                      slope.resids = NA, #tst.sim.ref$slope.resids,
                      mean.resids = NA, #tst.sim.ref$mean.resids,
                      deviance = deviance,
                      log.RNA.se = log.RNA.se,
                      accept.smp = accept.smp,
                      log.RNA = log.RNA,
                      log.rED50.prl = log.rED50.prl,
                      log.tED50.prl = log.rED50.prl + log.RNA,
                      b.prl = b.prl,
                      d.prl = d.prl,
                      s.prl = s.prl))
} # end summarize.sample()
