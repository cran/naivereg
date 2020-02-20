\name{TradeAndGrowthData}
\alias{TradeAndGrowthData}
\title{Trade and growth data}
\description{This is a revisit of the seminal work of Frankel and Romer (1999, AER) ''Does trade cause growth?" using new 2017 data.  NAIVE method is detailed in Fan and Zhong (2018, JBES), "Nonparametric additive instrumental variable estimator: a group shrinkage estimation perspective". We found that trade still has significant effect on growth, but compared to the original studies, the effect of trade seems to be smaller in magnitude and NAIVE is able to select the important IVs which yield unbiased estimation. See details of the study including the invalid IV discussion in Fan and Wu (2020).
  \itemize{
    \item N: log(economically active population in millions).
    \item A: log(land area).
    \item water: water area.
    \item coast: coastline.
    \item arable: arable percentage.
    \item border: land border.
    \item forest: forest percentage.
    \item lang: number of official languages.
    \item pm25: PM2.5.
    \item ww: waterways.
    \item rw: railway.
    \item in_water=T_hat*water.
    \item in_coast=T_hat*coast.
    \item in_arable=T_hat*arable.
    \item in_border=T_hat*border.
    \item in_forest=T_hat*forest.
    \item in_lang=T_hat*lang.
    \item ...
    }
}
\usage{data(TradeAndGrowthData)}
\keyword{datasets}

