# making a change in the z.hallyarboroughL function to take additional parameters
# pres.pc, temp.pc, pres.pr, temp.pr
z.hallyarboroughL <- function(pres.a, temp.f, gas.sg, pres.pc, temp.pc,
                              pres.pr, temp.pr,
                              n2.frac = 0, co2.frac = 0, h2s.frac = 0, ...)
{
    funcY <- function(y) {
        # implicit equation
        # in some literature A = A1, B = A2, C = A3, D = A4

        - A * pres.pr + (y + y^2 + y^3 - y^4) / (1 - y)^3  - B * y^2 + C * y^D
    }

    if (missing(gas.sg)) {
        print("No gas.sg supplied.")
        if (missing(pres.pc) || missing(temp.pc)) {
            if (missing(pres.pr) || missing(temp.pr)) stop()
            # providing pres.pr and temp.pr
            print("using pres.pr and temp.pr")
            temp.r <- 1 / temp.pr
            # calculate pseudo-reduced
            pres.pc <- pres.a / pres.pr
            temp.pc <- (temp.f + 460) / temp.pr  # worksheet has bug in the Farenheit add
        } else {
            if (missing(pres.pc) || missing(temp.pc)) stop()
            print("Using Ppc and Tpc")
            crit <- calcGasPseudoReduced(pres.a, pres.pc, temp.f, temp.pc)
            pres.pr <- crit$pres.pr
            temp.pr <- crit$temp.pr
            temp.r  <- crit$temp.r
        }
    } else {
        print("gas.sg has been provided")
        crit <- calcCriticals(pres.a, temp.f, gas.sg,
                              co2.frac = 0, h2s.frac = 0, n2.frac = 0)
        pres.pr <- crit$pres.pr
        temp.pr <- crit$temp.pr
        temp.r  <- crit$temp.r
        pres.pc <- crit$pres.pc
        temp.pc <- crit$temp.pc
    }
    # return(res)
    # pres.a = absolute pressure, psia
    # temp.f  = temperature, deg F

    # calculate pseudo-critical pressure and temperature
    # get pseudo-reduced
    # crit <- calcCriticals(pres.a, temp.f, gas.sg, co2.frac, h2s.frac, n2.frac, ...)


    t <- temp.r   # make it easier to read in the equation below

    A <- 0.06125 * t * exp(-1.2 * (1 - t)^2)
    B <- t * (14.76 - 9.76 * t + 4.58 * t^2)
    C <- t * (90.7 - 242.2 * t + 42.4 * t^2)
    D <- 2.18 + 2.82 * t

    All <- rootSolve::uniroot.all(funcY, c(-5.01, 5.99)) # find the root of the equation
    Y <- min(All)                         # minimum value
    z <- A * pres.pr / Y                  # calculate z
    zfactors <- named.list(z, Y, A, B, C, D,
                           pres.pr, temp.pr, pres.pc, temp.pc, temp.r)
    return(zfactors)
}

calcCriticals <- function(pres.a, temp.f, gas.sg,
                          co2.frac = 0, h2s.frac = 0, n2.frac = 0,
                          correction = "")
{
    # pseudocritical temperatures in Ranine

    if (correction == "Brown") {
        brown <- calcGasCriticals.Brown(pres.a, temp.f, gas.sg,
                                        co2.frac = 0, h2s.frac = 0, n2.frac = 0)
        pres.pc <- brown$pres.pc
        temp.pc <- brown$temp.pc
    } else {
        # if (h2s.frac < 0.03 & n2.frac < 0.05) {
        if (h2s.frac > 0.3 & n2.frac > 0.05) {
            # calculate pseudo-criticals.
            # Guo pg 2/22
            # cat("\nsmall H2S and N2\n")
            # Valid for H2S < 3%; N2 < 5% and inorganic < 7%
            pres.pc <- 709.604 - 58.718 * gas.sg      # Eq 2.22
            temp.pc <- 170.491 + 307.344 * gas.sg     # Eq 2.23
        } else {
            # print("Guo 2/23")
            # Pseudo-criticals for impurity corrections
            # Ahmed, Guo pg 2/23
            pres.pc <- 678 - 50 * (gas.sg - 0.5) - 206.7 * n2.frac +
                440 * co2.frac + 606.7 * h2s.frac
            temp.pc <- 326 + 315.7 * (gas.sg - 0.5) - 240 * n2.frac -
                83.3 * co2.frac + 133.3 * h2s.frac
        }
    }

    pr <- calcGasPseudoReduced(pres.a, pres.pc, temp.f, temp.pc)
    pres.pr <- pr$pres.pr
    temp.pr <- pr$temp.pr
    temp.r  <- pr$temp.r

    criticals <- list(pres.pr = pres.pr,
                      temp.pr = temp.pr,
                      temp.r = temp.r,
                      pres.pc = pres.pc,
                      temp.pc = temp.pc)
    return(criticals)
}


calcGasPseudoReduced <- function(pres.a, pres.pc, temp.f, temp.pc) {
    # calculate pseudo-reduced
    pres.pr <- pres.a / pres.pc
    temp.pr <- (temp.f + 460) / temp.pc  # worksheet has bug in the Farenheit add
    temp.r <- 1 / temp.pr                # wrong division in worksheet cell c15
    return(named.list(pres.pr, temp.pr, temp.r))
}

named.list <- function(...) {
    nl <- setNames( list(...) , as.character( match.call()[-1]) )
    # nl <- setNames( list(...) , as.character( match.call()[-1]) )
    nl
}