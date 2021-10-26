# Multifactor-rate-modelling

This is a project from my master's curriculum, which I have worked on with Xinglong TIAN, Hamza BERNOUSSI and Cedric NOGUE NOGHA. Our work consists in the analysis of different rates models to choose which of these will be the more adapted to price structured hybrid equity-rates that require the observation of many rates parameters (Libor, CMS, CMS spread).

We have significantly worked on the general framework of multifactor rates modelling where we explain the fundamental ideas behind multifactor rates model and all the theoretical requirements. The core elements of this approach are explained in the HJM model. To see the limitations of a one-factor model ant therefore the motivation of multifactor rates models, we used the Hull White one-factor model to develop a Monte-Carlo pricing of zero-coupons and the G1++ to price fictitious European call options on zero coupons, then we calibrate those prices with a set of given prices. That being done, we enlighten the advantages of the G2++ (two factors) to catch short term and long term structure, then we price swaptions using Monte-Carlo algorithm. The main steps of our work are presented in details in the attached report, and the language used for any numerical result is C++.
