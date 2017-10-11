#' Design for the data attentionshifting_signal
#'
#' @description Design of the an experiment measuring the brain response of 15 participants who have been shown images of neutral and angry faces. Those faces were shown at a different visibility 16ms and 166ms and were displayed either to the left or to the right of a screen. The laterality, sex, age, and 2 measures of anxiety of each subjects are also available. The value of the responses are located in the dataset \link{attentionshifting_signal}.
#'
#' \itemize{
#' \item id : identifier of the subject.
#' \item visibility : time of exposure to the image (16ms: subliminal or 166ms:supraliminal).
#' \item emotion : type of emotion of the image (angry or neutral).
#' \item direction : position of image one the screen (left or right).
#' \item laterality_id : measure of laterality of the subject.
#' \item age : age of the subject.
#' \item sex : sex of the subject.
#' \item STAIS_state : measure of the state of anxiety of the subject.
#' \item STAIS_trait : measure of the personality trait of anxiety of the subject.
#' }
#'
#' @name attentionshifting_design
#' @docType data
#' @usage data(attentionshifting_design)
#' @format A data frame with 120 rows and 10 variables.
NULL



#' Dataset containing th signal for a given electrod (O1) in an experiment using EEG
#'
#' The ERP (averaged EEG signal) for a given electrod (O1) in an experiment in attention shifting sampled at 1024 Hz. The design of the experiment is given in the dataset \link{attentionshifting_design}.
#'
#' \itemize{
#' \item ERP (in muV) of the electrod O1 measured from -200 to 600 timeframes before and after the onset of the stimulus.
#' }
#'
#' @name attentionshifting_signal
#' @docType data
#' @usage data(attentionshifting_signal)
#' @format A data frame with 120 rows and 819 variables.
NULL



#' Control study in psychology.
#'
#' A subset of a control experiment measuring the impulsive approach tendencies towrad physical activity or sedentary behaviors.
#'
#' \itemize{
#' \item id identifier of the subject.
#' \item bmi body mass index.
#' \item age.
#' \item sex.
#' \item condition the experimental condition where the task was to approach physical activity and avoid sedentary behavior (ApSB_AvPA), approach sedentarity behavior and avoid physical activity (ApPA_AvSB), and a control condition (control).
#' \item time pre, post.
#' \item iapa measure of impulsive approach tendencies toward physical activity (dependant variable).
#' \item iasb measure of impulsive approach tendencies toward sedentary behavior (dependant variable).
#' }
#'
#' @name jpah2016
#' @docType data
#' @usage data(jpah2016)
#' @format A data frame with 38 rows and 8 variables.
#' @references Cheval, B., Sarrazin, P., Pelletier, L., & Friese, M. (2016). Effect of retraining approach-avoidance tendencies on an exercise task: A randomized controlled trial. Journal of Physical Activity and Health, 13(12), 1396-1403.
NULL




#' Emergency patients data
#'
#' Observational data from 176 emergency patients with variables :
#'
#' \itemize{
#' \item \code{sex}.
#' \item \code{age}.
#' \item \code{insurance} : the type of insurance, private or semi private (\code{semi_private}) or public (\code{public}).
#' \item \code{LOS} : the length of the stay in days.
#' \item \code{cost} : the cost in CHF.
#' }
#'
#' @name emergencycost
#' @docType data
#' @usage data(emergencycost)
#' @format A data frame with 176 rows and 5 variables.
NULL
