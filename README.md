# Bayesian-randomized-basket-trial-design
Scripts for 'Bayesian randomized basket trial design: a case study from the ultra-rare invasive mould infections'

# Abstract
Invasive mould infections (IMIs) are rare but life-threatening. Regulatory approval for new antifungal drugs requires a well-powered, randomized non-inferiority trial, which is nearly infeasible due to the rarity of IMIs. Additionally, heterogeneity among mould types complicates study design and treatment effect interpretation when a study includes patients infected by different types of pathogens. Despite the success of single-arm oncology basket trials in evaluating treatment effect in multiple disease types, statistical methods for randomized basket trials in non-oncologic settings remain underdeveloped. We propose a robust borrowing strategy to enhance the efficiency of randomized basket trials for IMIs by (i) borrowing treatment effects across mould types while accounting for heterogeneity and (ii) augmenting control arms using external data. The proposed approach increases the efficiency and precision of the treatment effect estimates for various moulds. It also increases the ethical appeal by reducing the number of patients required for the control arm. Using simulation and real-life examples, we demonstrated the proposed approach can significantly increase statistical power and precision while maintaining the family-wise type I error rate at an acceptable level. Our approach offers a substantial improvement over the current practice of pooling different moulds together for inference and is applicable to rare disease trials facing similar accrual and ethical challenges.

# Subfolders
--Paper.Code: scripts to reproduce case study analysis.  

# Code Overview

### Paper.code:
--reg_bma.R/myfun.R/mem functions.R: scripts for implmenting BMA and BMA-2D methods  

--exnex.R/exnex utility.R: scripts for implementing EXNEX and EXNEX-2D methods  

--Re_Analysis of IMI Trial.Rmd: scripts for re_analysis of IMI trial used in Case Study





