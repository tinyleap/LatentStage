#' Snack Set
#'
#' In this dataset, there are 12 kinds of snacks for people to choose from.
#' Participants rated how much they liked each of the 12 snacks using an
#' 11-point Likert scale (1 = "dislike very much" and 11 = "like very much"),
#' which gives us the covariate "Rate". Following that, participants ranked the
#' 12 snacks, from their most favorite (rank = 1) to their least favorite (rank
#' = 12), which accounts for the covariate "Rank".
#'
#' @docType data
"SnackSet"

#' Candy Set
#'
#' In this dataset, there are 6 kinds of candies for people to choose from.
#' Participants rated how much they liked each of the 6 candies using an
#' 11-point Likert scale (1 = "dislike very much" and 11 = "like very much"),
#' which gives us the covariate "Rate". Following that, participants ranked the
#' 6 candies, from their most favorite (rank = 1) to their least favorite (rank
#' = 12), which accounts for the covariate "Rank".
#'
#' @docType data
"CandySet"

#' Consideration Set Simulation 1
#'
#' NumOpts=5;
#' NumPeople=5000;
#' NumCovs=3;
#' BetaConsid=[2 0 -1];
#' BetaSelect=[0 -.5 1];
#' ASCsConsid=[1 2 3 4 5]/5;
#' ASCsSelect=-[4 3 2 1 0]/4; % Last to zero for identification
#'
#' @docType data
"CFS5_1"

#' Consideration Set Simulation 2
#'
#' NumOpts=6;
#' NumPeople=5000;
#' NumCovs=3;
#' BetaConsid=[2 0 -1];
#' BetaSelect=[0 -.5 1];
#' ASCsConsid=[1 2 3 4 5 6]/5;
#' ASCsSelect=-[5 4 3 2 1 0]/4; % Last to zero for identification
#'
#' @docType data
"CFS6_1"

#' Consideration Set Simulation 3
#'
#' Here we generate some data where (1) there are "multiplicities" (an outcome
#' choice vector can be like [2 0 1 0 0]); and (2) some items are "not
#' available". We use the same values as before, but with 10% of the items not
#' available, and an average number of chosen items = 1.8 (instead of 1).
#'
#' @docType data
"CFS5_2"

#' Consideration Set Simulation 4
#'
#' Here we generate some data where (1) there are "multiplicities" (an outcome
#' choice vector can be like [2 0 1 0 0 1]); and (2) some items are "not
#' available". We use the same values as before, but with 10% of the items not
#' available, and an average number of chosen items = 1.8 (instead of 1).
#'
#' @docType data
"CFS6_2"

#' Three Stage
#'
#' Structure of the data is single stage simulation, collapsing multistage model
#' into standard Discrete Heterogeneity Model.
#'
#' @docType data
"threestage"

#' Dating
#'
#' Response variables are choices that apply to different stages. The first
#' stage corresponds to whether or not a user browsed a profile. The second
#' stage corresponds to whether or not a user wrote to a profile. For each
#' stage, we have only one identical covariate "agedif".
#'
#' @docType data
"dating"
